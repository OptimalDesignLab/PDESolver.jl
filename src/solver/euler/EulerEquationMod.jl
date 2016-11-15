# Description:
#   Declare the Euler data object
#   Includes all the files for the Euler module

module EulerEquationMod

include("complexify.jl")
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

using ArrayViews
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using ForwardDiff
using Utils
import ODLCommonTools.sview
using MPI
#using Debugging
# the AbstractEquation type is declared in ODLCommonTools
# every equation will have to declare a new type that is a subtype of AbstractEquation


export AbstractEulerData, EulerData, EulerData_

@doc """
### EulerEquationMod.ParamType

  This type holds the values of any constants or paramters needed during the
  computation.  These paramters can be specified in the opts dictionary or
  have default values set here.  If there is no reasonable default, values
  are initialized to -1

  gamma and R are the independent themodynamic variables

  Whether this type should be immutable or not is an open question

  This type is paramaterized on the dimension of the equation for purposes
  of multiple dispatch

  Static Parameters:
  Tdim : dimensionality of the equation, integer, (used for dispatch)
  var_type : type of variables used used in the weak form, symbol, (used for
             dispatch), currently supported values: :conservative, :entropy
  Tsol : datatype of solution variables q
  Tres : datatype of residual
  Tmsh : datatype of mesh related quantities (mapping jacobian etc.)

  **Fields**
    # these fields have defaults:
    * cv  : specific heat constant
    * R : specific gas constant (J/(Kg*K))
    * gamma : ratio of specific heats
    * gamma_1 : gamma - 1
    # these fields do not have defaults:
    * Ma  : free stream Mach number
    * Re  : free stream Reynolds number
    * aoa : angle of attack (radians)
 
"""->
type ParamType{Tdim, var_type, Tsol, Tres, Tmsh} <: AbstractParamType{Tdim}
  f::IOStream
  t::Float64  # current time value
  order::Int  # accuracy of elements (p=1,2,3...)

  q_vals::Array{Tsol, 1}  # resuable temporary storage for q variables at a node
  q_vals2::Array{Tsol, 1}
  q_vals3::Array{Tsol, 1}
  qg::Array{Tsol, 1}  # reusable temporary storage for boundary condition
  v_vals::Array{Tsol, 1}  # reusable storage for convert back to entropy vars.
  v_vals2::Array{Tsol, 1}

  # numDofPerNode x stencilsize arrays for entropy variables
  w_vals_stencil::Array{Tsol, 2}
  w_vals2_stencil::Array{Tsol, 2}

  res_vals1::Array{Tres, 1}  # reusable residual type storage
  res_vals2::Array{Tres, 1}  # reusable residual type storage

  flux_vals1::Array{Tres, 1}  # reusable storage for flux values
  flux_vals2::Array{Tres, 1}  # reusable storage for flux values

  sat_vals::Array{Tres, 1}  # reusable storage for SAT term

  A0::Array{Tsol, 2}  # reusable storage for the A0 matrix
  A0inv::Array{Tsol, 2}  # reusable storage for inv(A0)
  A1::Array{Tsol, 2}  # reusable storage for a flux jacobian
  A2::Array{Tsol, 2}  # reusable storage for a flux jacobian

  A_mats::Array{Tsol, 3}  # reusable storage for flux jacobians

  Rmat1::Array{Tres, 2}  # reusable storage for a matrix of type Tres
  Rmat2::Array{Tres, 2}

  nrm::Array{Tmsh, 1}  # a normal vector
  nrm2::Array{Tmsh, 1}
  nrm3::Array{Tmsh, 1}

  cv::Float64  # specific heat constant
  R::Float64  # specific gas constant used in ideal gas law (J/(Kg * K))
  gamma::Float64 # ratio of specific heats
  gamma_1::Float64 # = gamma - 1

  Ma::Float64  # free stream Mach number
  Re::Float64  # free stream Reynolds number
  aoa::Float64  # angle of attack
  rho_free::Float64  # free stream density
  E_free::Float64 # free stream energy (4th conservative variable)

  edgestab_gamma::Float64  # edge stabilization parameter
  # debugging options
  writeflux::Bool  # write Euler flux
  writeboundary::Bool  # write boundary data
  writeq::Bool # write solution variables
  use_edgestab::Bool  # use edge stabilization
  use_filter::Bool  # use filtering
  use_res_filter::Bool # use residual filtering

  filter_mat::Array{Float64, 2}  # matrix that performs filtering operation
                                 # includes transformations to/from modal representation

  use_dissipation::Bool  # use artificial dissipation
  dissipation_const::Float64  # constant used for dissipation filter matrix

  tau_type::Int  # type of tau to use for GLS stabilization

  vortex_x0::Float64  # vortex center x coordinate at t=0
  vortex_strength::Float64  # strength of the vortex

  krylov_itr::Int  # Krylov iteration number for iterative solve
  krylov_type::Int # 1 = explicit jacobian, 2 = jac-vec prod

  Rprime::Array{Float64, 2}  # numfaceNodes x numNodesPerElement interpolation matrix
                             # this should live in sbpface instead
  # temporary storage for calcESFaceIntegrals
  A::Array{Tres, 2}
  B::Array{Tres, 3}
  iperm::Array{Int, 1}
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
  t_barrier::Float64  # time spent in MPI_Barrier
  t_jacobian::Float64 # time spend computing Jacobian
  t_solve::Float64 # linear solve time
  =#
  time::Timings

  function ParamType(mesh, sbp, opts, order::Integer)
  # create values, apply defaults
    
    t = 0.0
    myrank = mesh.myrank
    #TODO: don't open a file in non-debug mode
    f = open("log_$myrank.dat", "w")
    q_vals = Array(Tsol, Tdim + 2)
    q_vals2 = Array(Tsol, Tdim + 2)
    q_vals3 = Array(Tsol, Tdim + 2)
    qg = Array(Tsol, Tdim + 2)
    v_vals = Array(Tsol, Tdim + 2)
    v_vals2 = Array(Tsol, Tdim + 2)
 
    w_vals_stencil = Array(Tsol, Tdim + 2, mesh.sbpface.stencilsize)
    w_vals2_stencil = Array(Tsol, Tdim + 2, mesh.sbpface.stencilsize)

    res_vals1 = Array(Tres, Tdim + 2)
    res_vals2 = Array(Tres, Tdim + 2)

    flux_vals1 = Array(Tres, Tdim + 2)
    flux_vals2 = Array(Tres, Tdim + 2)

    sat_vals = Array(Tres, Tdim + 2)

    A0 = zeros(Tsol, Tdim + 2, Tdim + 2)
    A0inv = zeros(Tsol, Tdim + 2, Tdim + 2)
    A1 = zeros(Tsol, Tdim + 2, Tdim + 2)
    A2 = zeros(Tsol, Tdim + 2, Tdim + 2)
    A_mats = zeros(Tsol, Tdim + 2, Tdim + 2, Tdim)

    Rmat1 = zeros(Tres, Tdim + 2, Tdim + 2)
    Rmat2 = zeros(Tres, Tdim + 2, Tdim + 2)

    nrm = zeros(Tmsh, Tdim)
    nrm2 = zeros(nrm)
    nrm3 = zeros(nrm)

    gamma = opts[ "gamma"]
    gamma_1 = gamma - 1
    R = opts[ "R"]
    cv = R/gamma_1

    Ma = opts[ "Ma"]
    Re = opts[ "Re"]
    aoa = opts[ "aoa"]
    E_free = 1/(gamma*gamma_1) + 0.5*Ma*Ma
    rho_free = 1.0

    edgestab_gamma = opts["edgestab_gamma"]

    # debugging options
    writeflux = opts[ "writeflux"]
    writeboundary = opts[ "writeboundary"]
    writeq = opts["writeq"]
    use_edgestab = opts["use_edgestab"]
    if use_edgestab println("edge stabilization enabled") end

    use_filter = opts["use_filter"]
    if use_filter println("solution variables filter enabled") end


    use_res_filter = opts["use_res_filter"]
    if use_res_filter println("residual filter enabled") end

    if use_filter || use_res_filter || opts["use_filter_prec"]
      filter_fname = opts["filter_name"]
      filter_mat = calcFilter(sbp, filter_fname, opts)
    else
      filter_mat = Array(Float64, 0,0)
    end

    use_dissipation = opts["use_dissipation"]
    if use_dissipation println("artificial dissipation enabled") end

    dissipation_const = opts["dissipation_const"]

    tau_type = opts["tau_type"]

    vortex_x0 = opts["vortex_x0"]
    vortex_strength = opts["vortex_strength"]

    krylov_itr = 0
    krylov_type = 1 # 1 = explicit jacobian, 2 = jac-vec prod

    sbpface = mesh.sbpface

    numNodesPerElement = mesh.numNodesPerElement
    Rprime = zeros(size(sbpface.interp, 2), numNodesPerElement)
    # expand into right size (used in SBP Gamma case)
    for i=1:size(sbpface.interp, 1)
      for j=1:size(sbpface.interp, 2)
        Rprime[j, i] = sbpface.interp[i, j]
      end
    end

    A = zeros(Tres, size(Rprime))
    B = zeros(Tres, numNodesPerElement, numNodesPerElement, 2)
    iperm = zeros(Int, size(sbpface.perm, 1))



    time = Timings()
    return new(f, t, order, q_vals, q_vals2, q_vals3,  qg, v_vals, v_vals2, 
               w_vals_stencil, w_vals2_stencil, res_vals1, 
               res_vals2, sat_vals, flux_vals1, 
               flux_vals2, A0, A0inv, A1, A2, A_mats, Rmat1, Rmat2, nrm, 
               nrm2, nrm3,cv, R, 
               gamma, gamma_1, Ma, Re, aoa, 
               rho_free, E_free,
               edgestab_gamma, writeflux, writeboundary, 
               writeq, use_edgestab, use_filter, use_res_filter, filter_mat, 
               use_dissipation, dissipation_const, tau_type, vortex_x0, 
               vortex_strength, 
               krylov_itr, krylov_type,
               Rprime, A, B, iperm,
               time)

    end   # end of ParamType function

end  # end type declaration

# add a layer of abstraction
@doc """
### EulerEquationMod.AbstractEulerData{Tsol, Tres}

  This abstract type should be the supertype of *all* solution data objects
  that are related to the Euler equations.

  It should be used for specify the type of a function argument only when
  the function does no operations on the solution data object itself, it just
  passes it onto other functions that do the work (thus AbstractEulerData 
  should be used for only the highest level functions).

  Another way of saying say it is that this type should only be used when 
  the function only needs to ensure that it is solving the Euler equations,
  but does not care even a little bit about how.
"""->
abstract AbstractEulerData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
#=
Use this type to leverage multiple dispatch.
This type holds any data, including large arrays of solution values, 
  that are specific to the equation for the Euler equations. 
  This includes the solutoin variables q, the fluxes in the xi and eta 
  direction, and the result of the calculation the inverse mass matrix is 
  stored here as well (not sure if that fits better here or in the mesh object)
  things like the coordinate field, the jacobian etc. are stored 
  in the mesh object

The aux_vars array holds all auxiliary variables that are stored over 
  the entire mesh.
Although it can be accessed directly, the preferred method is to use the macros 
  defined in the euler_macros.jl file.
These macros return either a scalar or an ArrayView of the specified indices, 
  depending on if the quantity requested is scalar or vector.
Every time a new variable is added to the array, the size must be updated 
  and a new macro must be created.

The uses of aux_vars should mirror that of eqn.q, in that entire columns 
  should be passed to low level functions and the low level functions 
  use the macros to access individual variables.
The advantages of macros vs functions for access to variables remains unclear
  if aux_vars is a fixed size.
If it is variable sized then macros give the advantage of doing location lookup
  at compile time
=#

#=
  Initial Condition:
    All initial condition functions use conservative variables, and are 
    converted to the variables used to solve the equation.
=#


#=
  Boundary conditions:
    All boundary condition calculations are done in conservative variables.
    If using some other variables, they are converted before being passed
    to the boundary condition functions.

    The functors that compute boundary conditions are gathered from BCDict 
    during init() and stored in mesh.bndry_funcs.  The functors compute the 
    flux from the boundary condition at a node.  The boundary fluxes for all 
    the nodes on the boundary are stored in mesh.bndryflux, to be integrated 
    later using boundaryintegrate!.
=#

#=
  Variable Conversion:
    If solving the equations is some set of variables other than conservative,
    (lets call them the v variables, and the conservative variables the q
     variables), the following functions must be defined:

    # convert to conservative variables at a node
    convertToConservative(params::ParamType, qe::AbstractArray{Tsol, 1},
                          qc::AbstractArray{Tsol, 1}

    # convert from conservative to the other variables at a node
    convertToEntropy(params::ParamType, qe::AbstractArray{Tsol, 1}, 
                     qc::AbstractArray{Tsol, 1})

    # convert from conservative variables to the variables the equation
    # is being solved in
    convertFromConservativeToWorkingVars(params::ParamType{Tdim, :var_type, 
                        qe::AbstractArray{Tsol, 1}, qc::AbstractArray{Tsol, 1})

    # calculate the coefficient matrix of the time derivate (dq/dv) at a node
    calcA0(params::ParamType, q::AbstractArray{Tsol, 1}, 
           A0::AbstractArray{Tsol, 2})
    # calculate inv(dq/dv)
    calcA0inv(params::ParamType, q::AbstractArray{Tsol, 1}, 
              A0::AbstractArray{Tsol, 2})

    # multiply a 3D array by inv(A0)
    matVecA0Inv(mesh, sbp, eqn, opts, arr::AbstractArray{Tsol, 3})

    # calculate the Euler flux at a node
    calcEulerFlux(params::ParamType, q::AbstractArray{Tsol, 1}, 
                  aux_vars::AbstractArray{Tsol, 1}, dir::AbstractArray{Tmsh}, 
                  F::AbstractArray{Tsol, 1})

=#

#=
  Calculating Other Quantities:
    To calculate physical quantities from the solution variables (for which
    the calculation is different for the different sets of variables) at a node,
    a separate function should be used, usnig the params object to dispatch 
    to the right function for the variables being used. 
    
    Some examples:

    calcPressure(params::ParamType, q::AbstractArray{Tsol, 1})
    calcEntropy(params::ParamType, q::AbstractArray{Tsol, 1})
    calcSpeedOfSound(params::ParamType, q::AbstractArray{Tsol, 1})
=#

@doc """
### EulerEquationMod.EulerData

  This type, although abstract, is the type functions should use for their
  input arguments if they do any operations on the solution data object.
  It stores all data used in evaluting the Euler Equations.
  
  It is paramaterized on the types Tsol, the type of the
  conservative variables q, and Tdim, the dimension of the equation

  It should have the following fields:
   * res_type : datatype of residual (depreciated)
    * q  : 3D array holding conservative variables
    * q_vec  : vector to assemble q into
    * aux_vars : 3D array holding auxiliary variables
    * flux_parametric : 4D array [ndof per node, nnodes per element, nelements, Tdim]
             holding the Euler flux in the xi and eta directions
    * res  : 3D array holding residual
    * res_vec   : vector form of res
    * edgestab_alpha : paramater used for edge stabilization, 4d array
    * bndryflux : 3D array holding boundary flux data
    * stabscale : 2D array holding edge stabilization scale factor
    * M : vector holding the mass matrix
    * Minv :  vector holding inverse mass matrix
"""->
abstract EulerData{Tsol, Tres, Tdim, var_type} <: AbstractEulerData{Tsol, Tres}

# high level functions should take in an AbstractEulerData, remaining
# agnostic to the dimensionality of the equation
# Mid level function should take in an EulerData{Tsol, Tdim}, so they
# know the dimensionality of the equation, but have a single method that
# can handle all dimensions

# low level functions should take in EulerData{Tsol, 2} or EulerData{Tsol, 3}
# this allows them to have different methods for different dimension equations.

# now that EulerData is declared, include other files that use it
include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("euler_macros.jl")
include("output.jl")
include("common_funcs.jl")
include("euler_funcs.jl")
include("conversion.jl")
include("euler.jl")
include("ic.jl")
include("bc.jl")
include("stabilization.jl")
include("flux.jl")
# include("artificialViscosity.jl")
# include("constant_diff.jl")
include("GLS2.jl")
include("boundary_functional.jl")
include("adjoint.jl")
include("source.jl")
include("entropy_flux.jl")

@doc """
### EulerEquationMod.EulerData_

  This type is an implimentation of the abstract EulerData.  It is
  paramterized by the residual datatype Tres and the mesh datatype Tmsh
  because it stores some arrays of those types.  Tres is the 'maximum' type of
  Tsol and Tmsh, where Tsol is the type of the conservative variables.
  It is also paremterized by var_type, which should be a symbol describing
  the set of variables stored in eqn.q.  Currently supported values are 
  :conservative and :entropy, which indicate the conservative variables and
  the entropy variables described in 'A New Finite Element Formulation for
  Computational Fluid Dynamics: Part I' by Hughes et al.

  Eventually there will be additional implimentations of EulerData,
  specifically a 3D one.

  Static Parameters:
    Tsol : datatype of variables solution variables, ie. the 
           q vector and array
    Tres : datatype of residual. ie. eltype(res_vec)
    Tdim : dimensionality of equation, integer, (2 or 3, currently only 2 is
           supported).
    Tmsh : datatype of mesh related quantities
    var_type : symbol describing variables used in weak form, (:conservative 
               or :entropy)


"""->
type EulerData_{Tsol, Tres, Tdim, Tmsh, var_type} <: EulerData{Tsol, Tres, Tdim, var_type}  
# hold any constants needed for euler equation, as well as solution and data 
#   needed to calculate it
# Formats of all arrays are documented in SBP.
# Only the constants are initilized here, the arrays are not.

  # this is the ParamType object that uses the same variables as
  # the EulerData_ object
  params::ParamType{Tdim, var_type, Tsol, Tres, Tmsh}
  comm::MPI.Comm
  commsize::Int
  myrank::Int

  # we include a ParamType object of all variable types, because occasionally
  # we need to do a calculation in  variables other than var_type
  # params (above) typically points to the same object as one of these
  params_conservative::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh}
  params_entropy::ParamType{Tdim, :entropy, Tsol, Tres, Tmsh}

  # the following arrays hold data for all nodes
  q::Array{Tsol,3}  # holds conservative variables for all nodes
  q_face::Array{Tsol, 4}  # store solution values interpolated to faces
  q_bndry::Array{Tsol, 3}  # store solution variables interpolated to 
  q_vec::Array{Tres,1}            # initial condition in vector form
  # hold fluxes in all directions
  # [ndof per node by nnodes per element by num element by num dimensions]
  aux_vars::Array{Tres, 3}        # storage for auxiliary variables 
  aux_vars_face::Array{Tres,3}    # storage for aux variables interpolated
                                  # to interior faces
  aux_vars_sharedface::Array{Array{Tres, 3}, 1}  # storage for aux varables interpolate
                                       # to shared faces
  aux_vars_bndry::Array{Tres,3}   # storage for aux variables interpolated 
                                  # to the boundaries
  flux_parametric::Array{Tsol,4}  # flux in xi and eta direction
  q_face_send::Array{Array{Tsol, 3}, 1}  # send buffers for sending q values
                                         # to other processes
  q_face_recv::Array{Array{Tsol, 3}, 1}  # recieve buffers for q values

  flux_face::Array{Tres, 3}  # flux for each interface, scaled by jacobian
  flux_sharedface::Array{Array{Tres, 3}, 1}  # hold shared face flux
  res::Array{Tres, 3}             # result of computation
  res_vec::Array{Tres, 1}         # result of computation in vector form
  Axi::Array{Tsol,4}               # Flux Jacobian in the xi-direction
  Aeta::Array{Tsol,4}               # Flux Jacobian in the eta-direction
  res_edge::Array{Tres, 4}       # edge based residual used for stabilization
                                  # numdof per node x nnodes per element x
				  # numEl x num edges per element

  edgestab_alpha::Array{Tmsh, 4}  # alpha needed by edgestabilization
                                  # Tdim x Tdim x nnodesPerElement x numEl
  bndryflux::Array{Tsol, 3}       # boundary flux
  stabscale::Array{Tsol, 2}       # stabilization scale factor

  # artificial dissipation operator:
  #   a square numnodes x numnodes matrix for every element
  dissipation_mat::Array{Tmsh, 3}  

  Minv::Array{Float64, 1}         # inverse mass matrix
  M::Array{Float64, 1}            # mass matrix

  # TODO: consider overloading getField instead of having function as
  #       fields
  disassembleSolution::Function   # function: q_vec -> eqn.q
  assembleSolution::Function      # function : eqn.res -> res_vec
  multiplyA0inv::Function         # multiply an array by inv(A0), where A0
                                  # is the coefficient matrix of the time 
				  # derivative
  majorIterationCallback::Function # called before every major (Newton/RK) itr

  src_func::SRCType  # functor for the source term
  flux_func::FluxType  # functor for the face flux
  volume_flux_func::FluxType  # functor for the volume flux numerical flux
                              # function
# minorIterationCallback::Function # called before every residual evaluation

  # inner constructor
  function EulerData_(mesh::AbstractMesh, sbp::AbstractSBP, opts)

    println("\nConstruction EulerData object")
    println("  Tsol = ", Tsol)
    println("  Tres = ", Tres)
    println("  Tdim = ", Tdim)
    println("  Tmsh = ", Tmsh)
    eqn = new()  # incomplete initialization

    eqn.comm = mesh.comm
    eqn.commsize = mesh.commsize
    eqn.myrank = mesh.myrank

    numfacenodes = mesh.numNodesPerFace

    vars_orig = opts["variable_type"]
    opts["variable_type"] = :conservative
    eqn.params_conservative = ParamType{Tdim, :conservative, Tsol, Tres, Tmsh}(
                                       mesh, sbp, opts, mesh.order)
    opts["variable_type"] = :entropy
    eqn.params_entropy = ParamType{Tdim, :entropy, Tsol, Tres, Tmsh}(
                                       mesh, sbp, opts, mesh.order)

    opts["variable_type"] = vars_orig
    if vars_orig == :conservative
      eqn.params = eqn.params_conservative
    elseif vars_orig == :entropy
      eqn.params = eqn.params_entropy
    else
      println(STDERR, "Warning: variable_type not recognized")
    end
    eqn.disassembleSolution = disassembleSolution
    eqn.assembleSolution = assembleSolution
    eqn.multiplyA0inv = matVecA0inv
    eqn.majorIterationCallback = majorIterationCallback

    eqn.Minv = calcMassMatrixInverse(mesh, sbp, eqn)
    eqn.M = calcMassMatrix(mesh, sbp, eqn)


    jac_type = opts["jac_type"]::Int
    if opts["use_dissipation"] || opts["use_dissipation_prec"]
      dissipation_name = opts["dissipation_name"]
      eqn.dissipation_mat = calcDissipationOperator(mesh, sbp, eqn, opts,
                                                    dissipation_name)
    else
      eqn.dissipation_mat = Array(Tmsh, 0, 0, 0)
    end
    
    # Must initialize them because some datatypes (BigFloat) 
    #   don't automatically initialize them
    # Taking a sview(A,...) of undefined values is illegal
    # I think its a bug that Array(Float64, ...) initializes values
    eqn.q = zeros(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)

    #TODO: don't store these, recalculate as needed
    eqn.Axi = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, sbp.numnodes,
                    mesh.numEl)
    eqn.Aeta = zeros(eqn.Axi)
    eqn.aux_vars = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.flux_parametric = zeros(Tsol, mesh.numDofPerNode, sbp.numnodes, 
                                mesh.numEl, Tdim)
    eqn.res = zeros(Tres, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)

    if opts["use_edge_res"]
      eqn.res_edge = zeros(Tres, mesh.numDofPerNode, sbp.numnodes, mesh.numEl, 
                           mesh.numTypePerElement[2])
    else
      eqn.res_edge = zeros(Tres, 0, 0, 0, 0)
    end

    if mesh.isDG
      eqn.q_vec = reshape(eqn.q, mesh.numDof)
      eqn.res_vec = reshape(eqn.res, mesh.numDof)
    else
      eqn.q_vec = zeros(Tres, mesh.numDof)
      eqn.res_vec = zeros(Tres, mesh.numDof)
    end

    eqn.edgestab_alpha = zeros(Tmsh,Tdim,Tdim,sbp.numnodes, mesh.numEl)
    if mesh.isDG
      eqn.q_face = zeros(Tsol, mesh.numDofPerNode, 2, numfacenodes, mesh.numInterfaces)
      eqn.flux_face = zeros(Tres, mesh.numDofPerNode, numfacenodes, mesh.numInterfaces)
      eqn.q_bndry = zeros(Tsol, mesh.numDofPerNode, numfacenodes, mesh.numBoundaryFaces)
      eqn.aux_vars_face = zeros(Tres, 1, numfacenodes, mesh.numInterfaces)
      eqn.aux_vars_bndry = zeros(Tres, 1, numfacenodes, mesh.numBoundaryFaces)
    else
      eqn.q_face = Array(Tres, 0, 0, 0, 0)
      eqn.flux_face = Array(Tres, 0, 0, 0)
      eqn.q_bndry = Array(Tsol, 0, 0, 0)
      eqn.aux_vars_face = zeros(Tres, 0, 0, 0)
      eqn.aux_vars_bndry = zeros(Tres, 0, 0, 0)
    end
    eqn.bndryflux = zeros(Tsol, mesh.numDofPerNode, numfacenodes, 
                          mesh.numBoundaryFaces)

    # send and receive buffers
    #TODO: rename buffers to not include face
    eqn.q_face_send = Array(Array{Tsol, 3}, mesh.npeers)
    eqn.q_face_recv = Array(Array{Tsol, 3}, mesh.npeers)
    eqn.flux_sharedface = Array(Array{Tres, 3}, mesh.npeers)
    eqn.aux_vars_sharedface = Array(Array{Tres, 3}, mesh.npeers)
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
      for i=1:mesh.npeers
        eqn.q_face_send[i] = Array(Tsol, mesh.numDofPerNode, dim2, 
                                         dim3_send[i])
        eqn.q_face_recv[i] = Array(Tsol,mesh.numDofPerNode, dim2,
                                        dim3_recv[i])
        eqn.flux_sharedface[i] = Array(Tres, mesh.numDofPerNode, numfacenodes, 
                                       mesh.peer_face_counts[i])
        eqn.aux_vars_sharedface[i] = Array(Tres, mesh.numDofPerNode, 
                                        numfacenodes, mesh.peer_face_counts[i])
      end
    end
   
    if eqn.params.use_edgestab
      eqn.stabscale = zeros(Tres, sbp.numnodes, mesh.numInterfaces) 
      calcEdgeStabAlpha(mesh, sbp, eqn)
    else
      eqn.stabscale = Array(Tres, 0, 0)
      eqn.edgestab_alpha = Array(Tmsh, 0, 0, 0, 0)
    end

    println("Tres = ", Tres)

    return eqn

  end  # end of constructor

end  # end of type declaration

@doc """
### EulerEquationMod.calcMassMatrixInverse

  This function calculates the inverse mass matrix and returns it.
  Because we use SBP operators, the mass matrix is diagonal, so it is stored
  in a vector.  mesh.dofs is used to put the components of the inverse
  mass matrix in the same place as the corresponding values in eqn.res_vec

  Arguments:
    mesh: AbstractMesh
    sbp: SBP operator
    eqn: an implementation of EulerData. Does not have to be fully initialized.

  Outputs:
    Minv: vector containing inverse mass matrix

"""->
# used by EulerData Constructor
# mid level functions
function calcMassMatrixInverse{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                                  sbp::AbstractSBP, 
                                                  eqn::EulerData{Tsol, Tres, Tdim})
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

@doc """
### EulerEquationMod.calcMassMatrix

  This function calculate the mass matrix and returns it.
  Beause w are using SBP operators, the mass matrix is diagonal, so it is
  stored in a vector.

  Arguments:
    mesh: AbstractMesh
    sbp: SBP operator
    eqn: an implementation of EulerData. Does not have to be fully initialized.

  Outputs:
    M: vector containing mass matrix

"""->
function calcMassMatrix{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                           sbp::AbstractSBP, 
                                           eqn::EulerData{Tsol, Tres, Tdim})
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
### EulerEquationMod.applyMassMatrixInverse

  This function multiplies eqn.res_vec (the residual in vector form  by eqn.Minv,
  the diagonal mass matrix.  This is a very fast function because all values
  are precomputed and stored linearly in memory.

  This is a mid level function, and does the correct thing regardless of the
  dimension of the equation.

  Aliasing restrictions: none
"""->
# mid level function (although it doesn't really need to Tdim)
function applyMassMatrixInverse{Tsol, Tres, Tdim}(eqn::EulerData{Tsol, Tres, Tdim}, 
                                            res_vec::AbstractVector{Tsol})
  # apply the inverse mass matrix stored eqn to res_vec

  ndof = length(res_vec)
  for i=1:ndof
    res_vec[i] *= eqn.Minv[i]
  end

  return nothing
end



end # end module
