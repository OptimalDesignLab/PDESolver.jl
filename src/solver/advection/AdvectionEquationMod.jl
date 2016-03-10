module AdvectionEquationMod

push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using ArrayViews
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using ForwardDiff
using Utils
export AdvectionData, AdvectionData_ #getMass, assembleSolution, disassembleSolution
export evalAdvection, init # exported from advectionFunctions.jl
export ICDict              # exported from ic.jl

# include("advectionFunctions.jl")
# include("getMass.jl")


abstract AbstractAdvectionData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
abstract AdvectionData{Tsol, Tres, Tdim} <: AbstractAdvectionData{Tsol, Tres}

@doc """
### AdvectionEquationMod.AdvectionData_

  This type is an implimentation of the abstract AdvectionData.  It is
  paramterized by the residual type Tres and the mesh type Tmsh
  because it stores some arrays of those types.  Tres is the 'maximum' type of
  Tsol and Tmsh, where Tsol is the type of the conservative variables.

"""->

type AdvectionData_{Tsol, Tres, Tdim, Tmsh} <: AdvectionData{Tsol, Tres, Tdim}

  # params::ParamType{Tdim}
  t::Float64
  res_type::DataType  # type of res
  alpha_x::Array{Tsol, 3}
  alpha_y::Array{Tsol, 3}
  q::Array{Tsol, 3}
  qface::Array{Tsol, 4}  # store solution values interpolated to faces
  aux_vars::Array{Tres, 3}  # storage for auxiliary variables 
  flux_parametric::Array{Tsol,4}  # flux in xi direction
  flux_face::Array{Tres, 3}  # flux for each interface, scaled by jacobian
  res::Array{Tres, 3}      # result of computation
  res_vec::Array{Tres, 1}  # result of computation in vector form
  res_edge::Array{Tres, 4} # edge based residual storage
  q_vec::Array{Tres,1}     # initial condition in vector form
  bndryflux::Array{Tsol, 3}  # boundary flux
  M::Array{Float64, 1}       # mass matrix
  Minv::Array{Float64, 1}    # inverse mass matrix
  disassembleSolution::Function # function u_vec -> eqn.q
  assembleSolution::Function    # function : eqn.res -> res_vec
  multiplyA0inv::Function       # multiply an array by inv(A0), where A0
                                # is the coefficient matrix of the time 
                                # derivative
  src_func::SRCType  # functor for source term
  majorIterationCallback::Function # called before every major (Newton/RK) itr

  function AdvectionData_(mesh::PumiMesh2, sbp::AbstractSBP, opts)
    println("\nConstruction AdvectionData object")
    println("  Tsol = ", Tsol)
    println("  Tres = ", Tres)
    println("  Tdim = ", Tdim)
    println("  Tmsh = ", Tmsh)
    eqn = new()  # incomplete initilization
    eqn.t = 0.0
    eqn.res_type = Tres
    eqn.alpha_x = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.alpha_y = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.disassembleSolution = disassembleSolution
    eqn.assembleSolution = assembleSolution
    eqn.majorIterationCallback = majorIterationCallback
    eqn.M = calcMassMatrix(mesh, sbp, eqn)
    eqn.Minv = calcMassMatrixInverse(mesh, sbp, eqn)
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
    eqn.bndryflux = zeros(Tsol, 1, sbp.numfacenodes, mesh.numBoundaryEdges)
    eqn.multiplyA0inv = matVecA0inv

    if mesh.isDG
      eqn.qface = zeros(Tsol, 1, 2, sbp.numfacenodes, mesh.numInterfaces)
      mesh.flux_face = zeros(Tres, 1, sbp.numfacenodes, mesh.numInterfaces)
    else
      eqn.qface = Array(Tres, 0, 0, 0, 0)
      mesh.flux_face = Array(Tres, 0, 0, 0)
    end

    return eqn
  end # ends the constructer AdvectionData_

end # End type AdvectionData_

include("advectionFunctions.jl")
include("common_funcs.jl")
include("boundaryconditions.jl")
include("bc_solvers.jl")
include("ic.jl")
include("GLS.jl")
 include("GLS2.jl")
include("../euler/complexify.jl")
include("source.jl")

@doc """
### AdvectionEquationMod.calcMassMatrix

  This function calculate the mass matrix and returns it.
  Beause w are using SBP operators, the mass matrix is diagonal, so it is
  stored in a vector.

  Arguments:
    mesh: AbstractMesh
    sbp: SBP operator
    eqn: an implementation of AdvectionData. Does not have to be fully initialized.

  Outputs:
    M: vector containing mass matrix

"""->
function calcMassMatrix{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                        sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim})
# calculate the (diagonal) mass matrix as a vector
# return the vector M

  M = zeros(Tmsh, mesh.numDof)

  for i=1:mesh.numEl
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
        dofnum_k = mesh.dofs[k,j,i]
        # multiplication is faster than division, so do the divions here
        # and then multiply solution vector times M
        M[dofnum_k] += (sbp.w[j]/mesh.jac[j,i])
      end
    end
  end

  return M

end     # end of calcMassMatrix function

@doc """
### AdvectionEquationMod.calcMassMatrixInverse

  This function calculates the inverse mass matrix and returns it.
  Because we use SBP operators, the mass matrix is diagonal, so it is stored
  in a vector.  mesh.dofs is used to put the components of the inverse
  mass matrix in the same place as the corresponding values in eqn.res_vec

  Arguments:
    mesh: AbstractMesh
    sbp: SBP operator
    eqn: an implementation of AdvectionData. Does not have to be fully initialized.

  Outputs:
    Minv: vector containing inverse mass matrix

"""->
# used by AdvectionData Constructor
# mid level functions
function calcMassMatrixInverse{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                                  sbp::AbstractSBP, 
                                                  eqn::AdvectionData{Tsol, Tres, Tdim})
# calculate the inverse mass matrix so it can be applied to the entire solution vector
# mass matrix is diagonal, stores in vector eqn.Minv

  Minv = zeros(Tmsh, mesh.numDof)

  for i=1:mesh.numEl
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
        dofnum_k = mesh.dofs[k,j,i]
        # multiplication is faster than division, so do the divisions here
        # and then multiply solution vector times Minv
        Minv[dofnum_k] += (sbp.w[j]/mesh.jac[j,i])
      end
    end
  end

  for i=1:mesh.numDof
    Minv[i] = 1/Minv[i]
  end

  return Minv

end     # end of calcMassMatrixInverse function

# functions needed to make it compatible with the NonLinearSolvers module
function matVecA0inv{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                     sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim},
                     opts, res_arr::AbstractArray{Tsol, 3})

  return nothing
end

function matVecA0{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                  sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim}, opts,
                  res_arr::AbstractArray{Tsol, 3})

  return nothing
end

end # end module
