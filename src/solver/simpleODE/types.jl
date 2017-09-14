# definitions of concrete subtypes of AbstractParamType and AbstractSolutionData

type ParamType{Tsol, Tres, Tdim} <: AbstractParamType

  f::BufferedIO     # TODO: needed?

  t::Float64
  time::Timings

  function ParamType(mesh, sbp, opts)
    myrank = mesh.myrank
    if DB_LEVEL >= 1
      f = BufferedIO("log_$myrank.dat", "w")
    else
      f = BufferedIO(DevNull)  # create a dummy IOStream
    end

    t = 0.0

    time = Timings()

    return new(f, t, time)
  end

end   # end of ParamType type def

typealias ParamType2{Tsol, Tres} ParamType{Tsol, Tres, 2}
typealias ParamType3{Tsol, Tres} ParamType{Tsol, Tres, 3}
typealias ParamTypes Union{ParamType2, ParamType3}

@doc """
### SimpleODEMod.SimpleODEData_

  This type is an implementation of the abstract SimpleODEData.  It is
  parameterized by the residual type Tres and the mesh type Tmsh
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
  shared_data::Array{SharedFaceData{Tsol}, 1}
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

    @assert mesh.commsize == 1

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
    # TODO: update this for parallel
    eqn.shared_data = Array(SharedFaceData{Tsol}, 0)

    return eqn
  end # ends the constructer SimpleODEData_

end # End type SimpleODEData_


import ODLCommonTools.getAllTypeParams

@doc """
### SimpleODEMod.getAllTypeParameters

Gets the type parameters for mesh and equation objects.

**Input**

* `mesh` : Object of abstract meshing type.
* `eqn`  : Euler Equation object.
* `opts` : Options dictionary

**Output**

* `tuple` : Tuple of type parameters. Ordering is same as that of the concrete eqn object within this physics module.

"""->
function getAllTypeParams{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, eqn::SimpleODEData_{Tsol, Tres, Tdim, Tmsh}, opts)

  tuple = (Tsol, Tres, Tdim, Tmsh)

  return tuple
end
