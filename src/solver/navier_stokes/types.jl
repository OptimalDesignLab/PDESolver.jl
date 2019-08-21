# define main types used in this module

"""
  [`AbstractParamTpe`](@ref) for Navier-Stokes physics.
"""
mutable struct ParamType{Tdim, Tsol, Tres, Tmsh} <: AbstractParamType{Tdim}
  # need this in case we ever want to call node level functions in the Euler
  # module
  euler_params::EulerEquationMod.ParamType{Tdim, :conservative, Tsol, Tres, Tmsh}
  # some of these fields alias those in euler_params so the outside world
  # see a consistent view
  f::BufferedIO
  t::Float64  # current time value (for unsteady simulations)
  order::Int  # accuracy of elements (p=1,2,3...)
  time::Timings
  isViscous::Bool  # can this ever be false?
  penalty_relaxation::Float64
  const_tii::Float64
  Re::Float64  # free stream Reynolds number

  function ParamType{Tdim, Tsol, Tres, Tmsh}(euler_eqn::EulerEquationMod.EulerData, mesh, sbp, opts, order::Integer) where {Tdim, Tsol, Tres, Tmsh} 
 

    # incomplete initialization
    params = new{Tdim, Tsol, Tres, Tmsh}()

    params.euler_params = euler_eqn.params
    params.f = params.euler_params.f  # alias
    params.t = 0.0
    params.order = order
    params.time = params.euler_params.time  # alias

    params.penalty_relaxation = opts["Cip"]
    params.isViscous = opts["isViscous"]

    if params.isViscous
      params.const_tii = calcTraceInverseInequalityConst(sbp, mesh.sbpface)
    end

    params.Re = opts["Re"]
    return params
  end

end  # end ParamType definition

"""
  Useful alias for 2D ParamType
"""
const ParamType2 = ParamType{2}

"""
  Useful alias for 3D ParamType
"""
const ParamType3 = ParamType{3}



"""
  Datatype for holding all the data required to solve the Naver-Stokes equations

  In particular, it contains the field `euler_eqn`, which is an
  [`EulerData`](@ref).  This object is used to evaluate the inviscid terms
  and generally allows the `NavierStokesMod` to build on top of
  `EulerEquationMod`.

  Note that many of the corresponding fields between `eqn_euler` and
  the `NSData` are aliased, including:

    * `q`
    * `q_vec`
    * `res`
    * `res_vec`
    * `res_edge`
    * `Minv`
    * `M`
"""
mutable struct NSData_{Tsol, Tres, Tdim, Tmsh} <: NSData{Tsol, Tres, Tdim}

  params::ParamType{Tdim, Tsol, Tres, Tmsh}
  comm::MPI.Comm
  commsize::Int
  myrank::Int

  # the nested EulerData object
  euler_eqn::EulerData_{Tsol, Tres, Tdim, Tmsh, :conservative}

  # required arrays
  q::Array{Tsol,3}        # holds conservative variables for all nodes
  q_vec::Array{Tres,1}    # initial condition in vector form
  res::Array{Tres, 3}     # element-based residual
  res_vec::Array{Tres, 1} # vector form of residual (aliases res in DG case)
  res_edge::Array{Tres, 4}       # edge based residual (not used)

  Minv3D::Array{Float64, 3}       # inverse mass matrix for application to res, not res_vec
  Minv::Array{Float64, 1}         # inverse mass matrix
  M::Array{Float64, 1}            # mass matrix

  # required functions
  multiplyA0inv::Function         # multiply an array by inv(A0), where A0
                                  # is the coefficient matrix of the time derivative
  majorIterationCallback::Function # called before every major (Newton/RK) itr

  # physics-specific stuff
  src_func::SRCType  # functor for the source term
  viscous_flux_func::FluxType  # functor for the viscous flux numerical flux function

  q_face::Array{Tsol, 4}     # q interpolated to face
  q_bndry::Array{Tsol, 3}    # q interpolated to boundary
  bndryflux::Array{Tres, 3}  # inviscid-type boundary flux
  flux_face::Array{Tres, 3}  # inviscid-type face flux
  area_sum::Array{Tmsh, 1}		  # the wet area of each element
  # vecflux_face::Array{Tres, 4}    # stores (u+ - u-)nx*, (numDofs, numNodes, numFaces)
  vecflux_faceL::Array{Tres, 4}     # stores (u+ - u-)nx*, (numDofs, numNodes, numFaces)
  vecflux_faceR::Array{Tres, 4}     # stores (u+ - u-)nx*, (numDofs, numNodes, numFaces)
  vecflux_bndry::Array{Tres, 4}     # stores (u+ - u-)nx*, (numDofs, numNodes, numFaces)

 
  file_dict::Dict{String, IO}  # dictionary of all files used for logging

  function NSData_{Tsol, Tres, Tdim, Tmsh}(mesh::AbstractMesh, sbp::AbstractOperator, opts) where {Tsol, Tres, Tdim, Tmsh}

    eqn = new{Tsol, Tres, Tdim, Tmsh}()  # incomplete initialization

    EulerEquationMod.checkOptions(opts)  # get any default values Euler needs
    eqn.euler_eqn = EulerEquationMod.EulerData_{Tsol, Tres, Tdim, Tmsh, 
                                                :conservative}(mesh, sbp, opts)

    eqn.params = ParamType{Tdim, Tsol, Tres, Tmsh}(eqn.euler_eqn, mesh, sbp, opts, mesh.order)
    eqn.comm = eqn.euler_eqn.comm
    eqn.commsize = eqn.euler_eqn.commsize
    eqn.myrank = eqn.euler_eqn.myrank

    # alias the arrays to avoid having to update them manually
    eqn.q = eqn.euler_eqn.q
    eqn.q_vec = eqn.euler_eqn.q_vec
    eqn.res = eqn.euler_eqn.res
    eqn.res_vec = eqn.euler_eqn.res_vec
    eqn.res_edge = eqn.euler_eqn.res_edge
    eqn.q_face = eqn.euler_eqn.q_face
    eqn.q_bndry = eqn.euler_eqn.q_bndry
    eqn.bndryflux = eqn.euler_eqn.bndryflux
    eqn.flux_face = eqn.euler_eqn.flux_face
    eqn.Minv3D = eqn.euler_eqn.Minv3D
    eqn.Minv = eqn.euler_eqn.Minv
    eqn.M = eqn.euler_eqn.M
    eqn.file_dict = eqn.euler_eqn.file_dict

    # because NSData is solving in the same variables as Euler, we can reuse
    # their A0inv function
    eqn.multiplyA0inv = eqn.euler_eqn.multiplyA0inv

    eqn.majorIterationCallback = majorIterationCallback

   if eqn.params.isViscous
     numfacenodes = mesh.numNodesPerFace
     numfaces = mesh.numInterfaces
     numBndFaces = mesh.numBoundaryFaces
     numvars  = mesh.numDofPerNode
     # eqn.vecflux_face = zeros(Tsol, Tdim, numvars, numfacenodes, numfaces)
     eqn.vecflux_faceL = zeros(Tsol, Tdim, numvars, numfacenodes, numfaces)
     eqn.vecflux_faceR = zeros(Tsol, Tdim, numvars, numfacenodes, numfaces)
     eqn.vecflux_bndry = zeros(Tsol, Tdim, numvars, numfacenodes, numBndFaces)
     eqn.area_sum = zeros(Tmsh, mesh.numEl)
     calcElemSurfaceArea(mesh, sbp, eqn)
   else
     # eqn.vecflux_face  = Array{Tsol}(0, 0, 0, 0)
     eqn.vecflux_faceL = Array{Tsol}(0, 0, 0, 0)
     eqn.vecflux_faceR = Array{Tsol}(0, 0, 0, 0)
     eqn.vecflux_bndry = Array{Tsol}(0, 0, 0, 0)
     eqn.area_sum = Array{Tsol}(0)
   end

   return eqn
 end  # end constructor

end  # end NSData_ type definition


import PDESolver: copyParameters

function copyParameters(eqn_old::NSData, eqn_new::NSData)

  EulerEquationMod.copyParameters(eqn_old.euler_params, eqn_new.euler_params)
  eqn_new.params.order = eqn_old.params.order
  eqn_new.params.isViscous = eqn_old.params.isViscous
  eqn_new.params.penalty_relaxation = eqn_old.params.penalty_relaxation
  eqn_new.params.const_tii = eqn_old.params.const_tii
  eqn_new.params.Re = eqn_old.params.Re

  return nothing
end

