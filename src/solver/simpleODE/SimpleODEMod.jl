module SimpleODEMod

using PDESolver
using ArrayViews
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using ForwardDiff
using NonlinearSolvers
using MPI
using Utils
using Input
import ODLCommonTools: get_uninitialized_SolutionData, sview
export SimpleODEData, SimpleODEData_ #getMass, assembleSolution, disassembleSolution
export evalResidual, init, run_simpleode # exported from simpleODE_funcs.jl
export ICDict              # exported from ic.jl
export ode_pre_func, ode_post_func    # exported from simpleODE_func.jl

type ParamType{Tsol, Tres, Tdim} <: AbstractParamType

  f::BufferedIO{IOStream}     # TODO: needed?

  t::Float64
  time::Timings

  function ParamType(mesh, sbp, opts)
    myrank = mesh.myrank
    if DB_LEVEL >= 1
      _f = open("log_$myrank.dat", "w")
      f = BufferedIO(_f)
    else
      f = BufferedIO()  # create a dummy IOStream
    end

    t = 0.0

    time = Timings()

    return new(f, t, time)
  end

end   # end of ParamType type def

typealias ParamType2{Tsol, Tres} ParamType{Tsol, Tres, 2}
typealias ParamType3{Tsol, Tres} ParamType{Tsol, Tres, 3}
typealias ParamTypes Union{ParamType2, ParamType3}
abstract AbstractSimpleODEData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
abstract SimpleODEData{Tsol, Tres, Tdim} <: AbstractSimpleODEData{Tsol, Tres}

@doc """
### SimpleODEMod.SimpleODEData_

  This type is an implementation of the abstract SimpleODEData.  It is
  paramterized by the residual type Tres and the mesh type Tmsh
  because it stores some arrays of those types.  Tres is the 'maximum' type of
  Tsol and Tmsh, where Tsol is the type of the conservative variables.

"""->
type SimpleODEData_{Tsol, Tres, Tdim, Tmsh} <: SimpleODEData{Tsol, Tres, Tdim}

  # params::ParamType{Tdim}
  params::ParamType{Tsol, Tres, Tdim}

  comm::MPI.Comm
  commsize::Int
  myrank::Int

  t::Float64
  res_type::DataType  # type of res
  q::Array{Tsol, 3}
  q_face::Array{Tsol, 4}  # store solution values interpolated to faces
  res::Array{Tres, 3}      # result of computation
  res_vec::Array{Tres, 1}  # result of computation in vector form
  res_edge::Array{Tres, 4} # edge based residual storage
  q_vec::Array{Tres,1}     # initial condition in vector form
  q_bndry::Array{Tsol, 3}  # store solution variables interpolated to 
                          # the boundaries with boundary conditions
  q_face_send::Array{Array{Tsol, 3}, 1}  # send buffers for sending q values
                                         # to other processes
  q_face_recv::Array{Array{Tsol, 3}, 1}  # recieve buffers for q values
  M::Array{Float64, 1}       # mass matrix
  Minv::Array{Float64, 1}    # inverse mass matrix
  Minv3D::Array{Float64, 3}    # inverse mass matrix for application to res, not res_vec
  disassembleSolution::Function # function u_vec -> eqn.q
  assembleSolution::Function    # function : eqn.res -> res_vec
  multiplyA0inv::Function       # multiply an array by inv(A0), where A0
                                # is the coefficient matrix of the time 
                                # derivative
  majorIterationCallback::Function # called before every major (Newton/RK) itr

  function SimpleODEData_(eqn::SimpleODEData_)

    return new()

  end

  function SimpleODEData_(mesh::AbstractMesh, sbp::AbstractSBP, opts)
    println("\nConstruction SimpleODEData object")
    println("  Tsol = ", Tsol)
    println("  Tres = ", Tres)
    println("  Tdim = ", Tdim)
    println("  Tmsh = ", Tmsh)
    numfacenodes = mesh.numNodesPerFace

    eqn = new()  # incomplete initialization
    eqn.comm = mesh.comm
    eqn.commsize = mesh.commsize
    eqn.myrank = mesh.myrank

    eqn.params = ParamType{Tsol, Tres, Tdim}(mesh, sbp, opts)
    eqn.t = 0.0
    eqn.res_type = Tres
    eqn.disassembleSolution = disassembleSolution
    eqn.assembleSolution = assembleSolution
    eqn.majorIterationCallback = majorIterationCallback
    eqn.M = calcMassMatrix(mesh, sbp, eqn)
    eqn.Minv = calcMassMatrixInverse(mesh, sbp, eqn)
    eqn.Minv3D = calcMassMatrixInverse3D(mesh, sbp, eqn)
    eqn.q = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.res = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.res_edge = Array(Tres, 0, 0, 0, 0)
    if mesh.isDG
      eqn.q_vec = reshape(eqn.q, mesh.numDof)
      eqn.res_vec = reshape(eqn.res, mesh.numDof)
    else
      eqn.q_vec = zeros(Tres, mesh.numDof)
      eqn.res_vec = zeros(Tres, mesh.numDof)
    end
    eqn.multiplyA0inv = matVecA0inv

    if mesh.isDG
      eqn.q_face = zeros(Tsol, 1, 2, numfacenodes, mesh.numInterfaces)
      eqn.q_bndry = zeros(Tsol, 1, numfacenodes, mesh.numBoundaryFaces)

      for i=1:mesh.npeers
        eqn.flux_sharedface[i] = zeros(Tres, 1, numfacenodes, 
                                       mesh.peer_face_counts[i])
      end
    else
      eqn.q_face = Array(Tres, 0, 0, 0, 0)
      eqn.q_bndry = Array(Tsol, 0, 0, 0)
    end

    # send and receive buffers
    #TODO: rename buffers to not include face
    eqn.q_face_send = Array(Array{Tsol, 3}, mesh.npeers)
    eqn.q_face_recv = Array(Array{Tsol, 3}, mesh.npeers)
    if mesh.isDG
      if opts["parallel_type"] == 1
        dim2 = numfacenodes
        dim3_send = mesh.peer_face_counts
        dim3_recv = mesh.peer_face_counts
      elseif opts["parallel_type"] == 2
        dim2 = mesh.numNodesPerElement
        dim3_send = mesh.local_element_counts
        dim3_recv = mesh.remote_element_counts
      else
        ptype = opts["parallel_type"]
        throw(ErrorException("Unsupported parallel type requested: $ptype"))
      end
    end
        
    for i=1:mesh.npeers
      eqn.q_face_send[i] = Array(Tsol, mesh.numDofPerNode, dim2, 
                                       dim3_send[i])
      eqn.q_face_recv[i] = Array(Tsol,mesh.numDofPerNode, dim2,
                                      dim3_recv[i])
    end

    return eqn
  end # ends the constructer SimpleODEData_

end # End type SimpleODEData_

# TODO: which of these
include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("simpleODE_funcs.jl")
include("common_funcs.jl")
include("ic.jl")
include("startup_func.jl")

global const PhysicsName = "SimpleODE"
register_physics(PhysicsName, SimpleODEMod, run_simpleode)

@doc """

# TODO
"""->
function get_uninitialized_SolutionData(eqn::SimpleODEData_)

  return SimpleODEData_(eqn)

end

@doc """
### SimpleODEMod.calcMassMatrix

  This function calculate the mass matrix and returns it.
  Beause w are using SBP operators, the mass matrix is diagonal, so it is
  stored in a vector.

  Arguments:
    mesh: AbstractMesh
    sbp: SBP operator
    eqn: an implementation of SimpleODEData. Does not have to be fully initialized.

  Outputs:
    M: vector containing mass matrix

"""->
function calcMassMatrix{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                        sbp::AbstractSBP, eqn::SimpleODEData{Tsol, Tres, Tdim})
# calculate the (diagonal) mass matrix as a vector
# return the vector M

  M = zeros(Tmsh, mesh.numDof)
  for i=1:mesh.numEl
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
        dofnum_k = mesh.dofs[k,j,i]
        # multiplication is faster than division, so do the division here
        # and then multiply solution vector times M
        M[dofnum_k] += (sbp.w[j]/mesh.jac[j,i])
      end
    end
  end

  return M

end     # end of calcMassMatrix function

# functions needed to make it compatible with the NonLinearSolvers module
function matVecA0inv{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                     sbp::AbstractSBP, eqn::SimpleODEData{Tsol, Tres, Tdim},
                     opts, res_arr::AbstractArray{Tsol, 3})

  return nothing
end

function matVecA0{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                  sbp::AbstractSBP, eqn::SimpleODEData{Tsol, Tres, Tdim}, opts,
                  res_arr::AbstractArray{Tsol, 3})

  return nothing
end

function majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSimpleODEData, opts)

  return nothing
end

end # end module
