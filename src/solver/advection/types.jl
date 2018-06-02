# definitions of concrete subtypes of AbstractParamType and AbstractSolutionData

"""
  Subtype of [`AbstractParamType`](@ref).

  **Static Parameters**:

   * Tsol
   * Tres
   *  Tdim

  This is a container passed to all low level function, useful for storing
  miscellaneous parameters or constants
"""
mutable struct ParamType{Tsol, Tres, Tdim} <: AbstractParamType{Tdim}
  LFalpha::Float64  # alpha for the Lax-Friedrich flux
  alpha_x::Float64
  alpha_y::Float64
  alpha_z::Float64
  sin_amplitude::Complex128
  omega::Complex128

  qL_s::Array{Tsol, 1}  # solution vector for a solution grid element
  qR_s::Array{Tsol, 1}  # solution vector for a solution grid element
  qL_f::Array{Tsol, 1}  # solution vector for flux grid element
  qR_f::Array{Tsol, 1}  # solution vector for flux grid element
  resL_s::Array{Tres, 1}  # residual for solution grid element
  resR_s::Array{Tres, 1}  # residual for solution grid element
  resL_f::Array{Tsol, 1}  # residual for a flux grid element
  resR_f::Array{Tsol, 1}  
  f::BufferedIO
  x_design::Array{Tres, 1}  # design variables
  time::Timings
  #=
  # timings
  t_volume::Float64  # time for volume integrals
  t_face::Float64 # time for surface integrals (interior)
  t_source::Float64  # time spent doing source term
  t_sharedface::Float64  # time for shared face integrals
  t_bndry::Float64  # time spent doing boundary integrals
  t_send::Float64  # time spent sending data
  t_wait::Float64  # time spent in MPI_Wait
  t_allreduce::Float64 # time spent in allreduce
  t_jacobian::Float64  # time spent computing jacobian
  t_solve::Float64  # linear solve time
  t_barrier::Float64  # time spent in MPI_Barrier
  t_barrier2::Float64
  t_barrier3::Float64
  t_barriers::Array{Float64, 1}
  =#
  function ParamType(mesh, sbp, opts)
    LFalpha = opts["LFalpha"]
    myrank = mesh.myrank
    if DB_LEVEL >= 1
      _f = open("log_$myrank.dat", "w")
      f = BufferedIO(_f)
    else
      f = BufferedIO()  # create a dummy IOStream
    end

    advection_velocity = opts["advection_velocity"]
    alpha_x = advection_velocity[1]
    alpha_y = advection_velocity[2]
    if Tdim == 3
      alpha_z = advection_velocity[3]
    else
      alpha_z = 1.0
    end

    numNodesPerElement_s = mesh.numNodesPerElement
    if opts["use_staggered_grid"]
      numNodesPerElement_f = mesh.mesh2.numNodesPerElement
    else
      numNodesPerElement_f = numNodesPerElement_s
    end

    qL_s = zeros(Tsol, numNodesPerElement_s)
    qR_s = zeros(Tsol, numNodesPerElement_s)
    qL_f = zeros(Tsol, numNodesPerElement_f)
    qR_f = zeros(Tsol, numNodesPerElement_f)

    resL_s = zeros(Tsol, numNodesPerElement_s)
    resR_s = zeros(Tsol, numNodesPerElement_s)
    resL_f = zeros(Tsol, numNodesPerElement_f)
    resR_f = zeros(Tsol, numNodesPerElement_f)

    # needed for the runtype=660 (CN uadj) objective
    sin_amplitude = 2.0
    omega = 1.0

    x_design = zeros(Tres, 0)  # TODO: get proper size information

    t = Timings()

    return new(LFalpha, alpha_x, alpha_y, alpha_z,
               sin_amplitude, omega,
               qL_s, qR_s, qL_f, qR_f, resL_s, resR_s, resL_f, resR_f,
               f, x_design, t)
  end
end

"""
  Convenient alias for all 2D ParamTypes
"""
ParamType2{Tsol, Tres} =  ParamType{Tsol, Tres, 2}

"""
  Convenient alias for all 3D ParamTypes
"""
ParamType3{Tsol, Tres} =  ParamType{Tsol, Tres, 3}

"""
  All ParamTypes, without static parameters specified
"""
const ParamTypes = Union{ParamType2, ParamType3}

"""
### AdvectionEquationMod.AdvectionData_

  This type is an implementation of the abstract AdvectionData.  It is
  parameterized by Tsol, the datatype of the solution variables and Tmsh,
  the datatype of the mesh variables.
  Tres is the 'maximum' type of Tsol and Tmsh.
  Tdim is the dimensionality of the equation being solve (2 or 3).

  This type is (ultimately) a subtype of [`AbstractSolutionData`](@ref) 
  and contains all the required fields.
"""
mutable struct AdvectionData_{Tsol, Tres, Tdim, Tmsh} <: AdvectionData{Tsol, Tres, Tdim}

  # params::ParamType{Tdim}
  params::ParamType{Tsol, Tres, Tdim}

  comm::MPI.Comm
  commsize::Int
  myrank::Int

  t::Float64
  res_type::DataType  # type of res
  q::Array{Tsol, 3}
  q_face::Array{Tsol, 4}  # store solution values interpolated to faces
  aux_vars::Array{Tres, 3}  # storage for auxiliary variables
  flux_parametric::Array{Tsol,4}  # flux in xi direction
  flux_face::Array{Tres, 3}  # flux for each interface, scaled by jacobian
  res::Array{Tres, 3}      # result of computation
  res_vec::Array{Tres, 1}  # result of computation in vector form
  res_edge::Array{Tres, 4} # edge based residual storage
  q_vec::Array{Tres,1}     # initial condition in vector form
  q_bndry::Array{Tsol, 3}  # store solution variables interpolated to
                          # the boundaries with boundary conditions
  shared_data::Array{SharedFaceData{Tsol}, 1}  # MPI data, including send and receive
                                         # buffers
#  q_face_send::Array{Array{Tsol, 3}, 1}  # send buffers for sending q values
                                         # to other processes
#  q_face_recv::Array{Array{Tsol, 3}, 1}  # recieve buffers for q values
  flux_sharedface::Array{Array{Tres, 3}, 1}  # hold shared face flux
  bndryflux::Array{Tsol, 3}  # boundary flux
  M::Array{Float64, 1}       # mass matrix
  Minv::Array{Float64, 1}    # inverse mass matrix
  Minv3D::Array{Float64, 3}    # inverse mass matrix for application to res, not res_vec
  multiplyA0inv::Function       # multiply an array by inv(A0), where A0
                                # is the coefficient matrix of the time
                                # derivative
  src_func::SRCType  # functor for source term
  flux_func::FluxType  # functor for the face flux
  majorIterationCallback::Function # called before every major (Newton/RK) itr

  function AdvectionData_(mesh::AbstractMesh, sbp::AbstractSBP, opts)
    println("\nConstruction AdvectionData object")
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
    eqn.majorIterationCallback = majorIterationCallback
    eqn.M = calcMassMatrix(mesh, sbp, eqn)
    eqn.Minv = calcMassMatrixInverse(mesh, sbp, eqn)
    eqn.Minv3D = calcMassMatrixInverse3D(mesh, sbp, eqn)
    eqn.q = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.aux_vars = zeros(Tsol, 0, 0, 0)

    if opts["precompute_volume_flux"]
      eqn.flux_parametric = zeros(Tsol, 1, mesh.numNodesPerElement, mesh.numEl,
                                  Tdim)
    else
      eqn.flux_parametric = zeros(Tsol, 0, 0, 0, 0)
    end

    eqn.res = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.res_edge = zeros(Tres, 0, 0, 0, 0)
    if mesh.isDG
      eqn.q_vec = reshape(eqn.q, mesh.numDof)
      eqn.res_vec = reshape(eqn.res, mesh.numDof)
    else
      eqn.q_vec = zeros(Tres, mesh.numDof)
      eqn.res_vec = zeros(Tres, mesh.numDof)
    end

    if opts["precompute_boundary_flux"]
      eqn.bndryflux = zeros(Tsol, 1, numfacenodes, mesh.numBoundaryFaces)
    else
      eqn.bndryflux = zeros(Tsol, 0, 0, 0)
    end

    eqn.multiplyA0inv = matVecA0inv

    if opts["precompute_q_face"]
      eqn.q_face = zeros(Tsol, 1, 2, numfacenodes, mesh.numInterfaces)
    else
      eqn.q_face = zeros(Tsol, 0, 0, 0, 0)
    end

    if opts["precompute_q_bndry"]
      eqn.q_bndry = zeros(Tsol, 1, numfacenodes, mesh.numBoundaryFaces)
    else
      eqn.q_bndry = zeros(Tsol, 0, 0, 0)
    end

    if opts["precompute_face_flux"]
      eqn.flux_face = zeros(Tres, 1, numfacenodes, mesh.numInterfaces)

      if mesh.isDG
        eqn.flux_sharedface = Array(Array{Tres, 3}, mesh.npeers)
        for i=1:mesh.npeers
          eqn.flux_sharedface[i] = zeros(Tres, 1, numfacenodes,
                                         mesh.peer_face_counts[i])
        end
      else
        eqn.flux_sharedface = Array(Array{Tres, 3}, 0)
      end  # end if isDG

    else
      eqn.flux_face = zeros(Tres, 0, 0, 0)
      eqn.flux_sharedface = Array(Array{Tres, 3}, 0)
    end  # end if precompute_face_flux

    if mesh.isDG
      eqn.shared_data = getSharedFaceData(Tsol, mesh, sbp, opts)
    else
      eqn.shared_data = Array(SharedFaceData, 0)
    end

    return eqn
  end # ends the constructor AdvectionData_

end # End type AdvectionData_


import ODLCommonTools.getAllTypeParams

@doc """
### AdvectionEquationMod.getAllTypeParameters

Gets the type parameters for mesh and equation objects.

**Input**

* `mesh` : Object of abstract meshing type.
* `eqn`  : Euler Equation object.
* `opts` : Options dictionary

**Output**

* `tuple` : Tuple of type parameters. Ordering is same as that of the concrete eqn object within this physics module.

"""->
function getAllTypeParams(mesh::AbstractMesh{Tmsh}, eqn::AdvectionData_{Tsol, Tres, Tdim, Tmsh}, opts) where {Tmsh, Tsol, Tres, Tdim}

  tuple = (Tsol, Tres, Tdim, Tmsh)

  return tuple
end
