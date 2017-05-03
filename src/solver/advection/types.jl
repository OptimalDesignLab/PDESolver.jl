# definitions of concrete subtypes of AbstractParamType and AbstractSolutionData

type ParamType{Tsol, Tres, Tdim} <: AbstractParamType{Tdim}
  LFalpha::Float64  # alpha for the Lax-Friedrich flux
  alpha_x::Float64
  alpha_y::Float64
  alpha_z::Float64

  f::BufferedIO
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
    alpha_x = 1.0
#     alpha_x = 0.0
    # Note: alpha_y = 0.0 might be useful for testing out new methods,
    #    but the CI tests will fail unless set to 1.0
    alpha_y = 1.0
#     alpha_y = 0.0
    alpha_z = 1.0


    t = Timings()
    return new(LFalpha, alpha_x, alpha_y, alpha_z, f, t)
  end
end

typealias ParamType2{Tsol, Tres} ParamType{Tsol, Tres, 2}
typealias ParamType3{Tsol, Tres} ParamType{Tsol, Tres, 3}
typealias ParamTypes Union{ParamType2, ParamType3}

@doc """
### AdvectionEquationMod.AdvectionData_

  This type is an implementation of the abstract AdvectionData.  It is
  paramterized by the residual type Tres and the mesh type Tmsh
  because it stores some arrays of those types.  Tres is the 'maximum' type of
  Tsol and Tmsh, where Tsol is the type of the conservative variables.

"""->
type AdvectionData_{Tsol, Tres, Tdim, Tmsh} <: AdvectionData{Tsol, Tres, Tdim}

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
  q_face_send::Array{Array{Tsol, 3}, 1}  # send buffers for sending q values
                                         # to other processes
  q_face_recv::Array{Array{Tsol, 3}, 1}  # recieve buffers for q values
  flux_sharedface::Array{Array{Tres, 3}, 1}  # hold shared face flux
  bndryflux::Array{Tsol, 3}  # boundary flux
  M::Array{Float64, 1}       # mass matrix
  Minv::Array{Float64, 1}    # inverse mass matrix
  Minv3D::Array{Float64, 3}    # inverse mass matrix for application to res, not res_vec
  disassembleSolution::Function # function u_vec -> eqn.q
  assembleSolution::Function    # function : eqn.res -> res_vec
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

    eqn = new()  # incomplete initilization
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
    eqn.flux_parametric = zeros(Tsol, 1, mesh.numNodesPerElement, mesh.numEl, Tdim)
    eqn.res = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.res_edge = Array(Tres, 0, 0, 0, 0)
    if mesh.isDG
      eqn.q_vec = reshape(eqn.q, mesh.numDof)
      eqn.res_vec = reshape(eqn.res, mesh.numDof)
    else
      eqn.q_vec = zeros(Tres, mesh.numDof)
      eqn.res_vec = zeros(Tres, mesh.numDof)
    end
    eqn.bndryflux = zeros(Tsol, 1, numfacenodes, mesh.numBoundaryFaces)
    eqn.multiplyA0inv = matVecA0inv

    if mesh.isDG
      eqn.q_face = zeros(Tsol, 1, 2, numfacenodes, mesh.numInterfaces)
      eqn.flux_face = zeros(Tres, 1, numfacenodes, mesh.numInterfaces)
      eqn.q_bndry = zeros(Tsol, 1, numfacenodes, mesh.numBoundaryFaces)
      eqn.flux_sharedface = Array(Array{Tres, 3}, mesh.npeers)

      for i=1:mesh.npeers
        eqn.flux_sharedface[i] = zeros(Tres, 1, numfacenodes,
                                       mesh.peer_face_counts[i])
      end
    else
      eqn.q_face = Array(Tres, 0, 0, 0, 0)
      eqn.flux_face = Array(Tres, 0, 0, 0)
      eqn.q_bndry = Array(Tsol, 0, 0, 0)
    end

    # send and receive buffers
    #TODO: rename buffers to not include face
    eqn.q_face_send = Array(Array{Tsol, 3}, mesh.npeers)
    eqn.q_face_recv = Array(Array{Tsol, 3}, mesh.npeers)
    if mesh.isDG
      if opts["parallel_data"] == "face"
        dim2 = numfacenodes
        dim3_send = mesh.peer_face_counts
        dim3_recv = mesh.peer_face_counts
      elseif opts["parallel_data"] == "element"
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
  end # ends the constructor AdvectionData_

end # End type AdvectionData_

@doc """
###AdvectionEquationMod.QfluxData

Data type for storing relevant information pertaining to an a functional or an
objective function.

**Members**

*  `is_objective_fn` : Bool whether the functional object is an objective
                       function or not.
*  `geom_faces_functional` : Geometric faces on which the functional is to be
                             computed.
*  `val` : Computed value of the functional
*  `target_qFlux` : Target value for the functional qFlux

The object is constructed using an inner constructor with the following 

**Arguments**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Advection equation object
*  `geom_faces_functional` : Geometric faces on which the functional is to be 
                             computed

"""->

type QfluxData{Topt} <: AbstractOptimizationData
  is_objective_fn::Bool
  geom_faces_functional::AbstractArray{Int,1}
  val::Topt
  target_qflux::Topt
  function QfluxData(mesh, sbp, eqn, opts, geom_faces_functional)

    functional = new()
    functional.is_objective_fn = false
    functional.geom_faces_functional = geom_faces_functional
    functional.val = zero(Topt)
    functional.target_qflux = zero(Topt)
    return functional
  end
end # End type OfluxData
