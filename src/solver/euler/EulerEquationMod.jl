# Description:
#   Declare the Euler data object
#   Includes all the files for the Euler module

module EulerEquationMod

include("complexify.jl")

using ArrayViews
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using ForwardDiff

# the AbstractEquation type is declared in ODLCommonTools
# every equation will have to declare a new type that is a subtype of AbstractEquation

export AbstractEulerData, EulerData, EulerData_

#=
Use this type to leverage multiple dispatch.
This type holds any data, including large arrays of solution values, 
  that are specific to the equation for the Euler equations. 
  This includes the conservative variables q, the fluxes in the xi and eta 
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
type ParamType{Tdim, var_type, Tsol, Tres, Tmsh} <: AbstractParamType
  t::Float64  # current time value
  order::Int  # accuracy of elements (p=1,2,3...)

  q_vals::Array{Tsol, 1}  # resuable temporary storage for q variables at a node

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


  vortex_x0::Float64  # vortex center x coordinate at t=0
  vortex_strength::Float64  # strength of the vortex

  krylov_itr::Int  # Krylov iteration number for iterative solve
  krylov_type::Int # 1 = explicit jacobian, 2 = jac-vec prod
  function ParamType(sbp, opts, order::Integer)
  # create values, apply defaults
    
    t = 0.0
    # get() = get(dictionary, key, default)

    q_vals = Array(Tsol, 4)

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

    vortex_x0 = opts["vortex_x0"]
    vortex_strength = opts["vortex_strength"]

    krylov_itr = 0
    krylov_type = 1 # 1 = explicit jacobian, 2 = jac-vec prod

    return new(t, order, q_vals, cv, R, gamma, gamma_1, Ma, Re, aoa, rho_free, E_free, 
               edgestab_gamma, writeflux, writeboundary, writeq, use_edgestab, 
               use_filter, use_res_filter, filter_mat, use_dissipation,  
               dissipation_const, vortex_x0, vortex_strength, krylov_itr, krylov_type)

    end   # end of ParamType function

end  # end type declaration

# add a layer of abstraction - although this migh be unnecessary
abstract AbstractEulerData{Tsol} <: AbstractSolutionData{Tsol}

@doc """
### EulerEquationMod.EulerData

  This type, although abstract, is the type functions should use for their
  input arguments.  It stores all data used in evaluting the Euler Equations.
  
  It is paramaterized on the types Tsol, the type of the
  conservative variables q, and Tdim, the dimension of the equation

  **Fields**
   * res_type : datatype of residual (depreciated)
    * q  : 3D array holding conservative variables
    * aux_vars : 3D array holding auxiliary variables
    * flux_parametric : 4D array [ndof per node, nnodes per element, nelements, Tdim]
             holding the Euler flux in the xi and eta directions
    * res  : 3D array holding residual
    * res_vec   : vector form of res
    * q_vec  : initial condition vector
    * edgestab_alpha : paramater used for edge stabilization, 4d array
    * bndryflux : 3D array holding boundary flux data
    * stabscale : 2D array holding edge stabilization scale factor
    * Minv :  vector holding inverse mass matrix
"""->
abstract EulerData {Tsol, Tdim, Tres, var_type} <: AbstractEulerData{Tsol}
#TODO: static parameter order is inconsistent with EulerData_

# high level functions should take in an AbstractEulerData, remaining
# agnostic to the dimensionality of the equation
# Mid level function should take in an EulerData{Tsol, Tdim}, so they
# know the dimensionality of the equation, but have a single method that
# can handle all dimensions

# low level functions should take in EulerData{Tsol, 2} or EulerData{Tsol, 3}
# this allows them to have different methods for different dimension equations.

# now that EulerData is declared, include other files that use it

include("euler_macros.jl")
include("output.jl")
include("common_funcs.jl")
include("euler.jl")
include("euler_funcs.jl")
include("ic.jl")
include("bc.jl")
include("stabilization.jl")
include("SUPG.jl")
# include("artificialViscosity.jl")
# include("constant_diff.jl")


@doc """
### EulerEquationMod.EulerData_

  This type is an implimentation of the abstract EulerData.  It is
  paramterized by the residual type Tres and the mesh type Tmsh
  because it stores some arrays of those types.  Tres is the 'maximum' type of
  Tsol and Tmsh, where Tsol is the type of the conservative variables.

  Eventually there will be additional implimentation of EulerData,
  specifically a 3D one.

  Static Parameters:
    Tsol : datatype of variables solution variables, ie. the 
           q vector and array
    Tres : datatype of residual. ie. eltype(res_vec)
    Tdim : dimensionality of equation, integer, (2 or 3)
    Tmsh : datatype of mesh related quantities
    var_type : type of variables used in weak form, symbol, (:conservative 
               or :entropy)


"""->
type EulerData_{Tsol, Tres, Tdim, Tmsh, var_type} <: EulerData{Tsol, Tdim, Tres, var_type}  
# hold any constants needed for euler equation, as well as solution and data 
#   needed to calculate it
# Formats of all arrays are documented in SBP.
# Only the constants are initilized here, the arrays are not.

  params::ParamType{Tdim, var_type, Tsol, Tres, Tmsh}

  # the following arrays hold data for all nodes
  q::Array{Tsol,3}  # holds conservative variables for all nodes
  q_vec::Array{Tres,1}            # initial condition in vector form
  # hold fluxes in all directions
  # [ndof per node by nnodes per element by num element by num dimensions]
  aux_vars::Array{Tres, 3}        # storage for auxiliary variables 
  flux_parametric::Array{Tsol,4}  # flux in xi and eta direction
  res::Array{Tres, 3}             # result of computation
  res_vec::Array{Tres, 1}         # result of computation in vector form
  Axi::Array{Tsol,4}               # Flux Jacobian in the xi-direction
  Aeta::Array{Tsol,4}               # Flux Jacobian in the eta-direction
  res_edge::Array{Tres, 4}       # edge based residual used for stabilization
                                  # numdof per node x nnodes per element x
				  # numEl x num edges per element

  edgestab_alpha::Array{Tmsh, 4}  # alpha needed by edgestabilization
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
  convertToEntropyVars::Function  # function eqn.q -> eqn.q entropy variables
  convertToConsVars::Function     # convert to eqn.q -> eqn.q in conservative
  multiplyA0inv::Function         # multiply an array by inv(A0), where A0
                                  # is the coefficient matrix of the time 
				  # derivative
  majorIterationCallback::Function # called before every major (Newton/RK) itr
# minorIterationCallback::Function # called before every residual evaluation

  # inner constructor
  function EulerData_(mesh::PumiMesh2, sbp::SBPOperator, opts)

    println("\nConstruction EulerData object")
    println("  Tsol = ", Tsol)
    println("  Tres = ", Tres)
    println("  Tdim = ", Tdim)
    println("  Tmsh = ", Tmsh)
    eqn = new()  # incomplete initialization

    eqn.params = ParamType{Tdim, var_type, Tsol, Tres, Tmsh}(sbp, opts, mesh.order)
    eqn.disassembleSolution = disassembleSolution
    eqn.assembleSolution = assembleSolution
    eqn.convertToEntropyVars = convertToEntropy
    eqn.convertToConsVars = convertToConservative
    eqn.multiplyA0inv = matVecA0inv
    eqn.majorIterationCallback = majorIterationCallback

    eqn.Minv = calcMassMatrixInverse(mesh, sbp, eqn)
    eqn.M = calcMassMatrix(mesh, sbp, eqn)

    calcEdgeStabAlpha(mesh, sbp, eqn)

    jac_type = opts["jac_type"]::Int
    if opts["use_dissipation"] || opts["use_dissipation_prec"]
      dissipation_name = opts["dissipation_name"]
      eqn.dissipation_mat = calcDissipationOperator(mesh, sbp, eqn, 
                                                    dissipation_name, opts)
    else
      eqn.dissipation_mat = Array(Tmsh, 0, 0, 0)
    end
    
    # Must initialize them because some datatypes (BigFloat) 
    #   don't automatically initialize them
    # Taking a view(A,...) of undefined values is illegal
    # I think its a bug that Array(Float64, ...) initializes values
    eqn.q = zeros(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
    eqn.Axi = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, sbp.numnodes,
                    mesh.numEl)
    eqn.Aeta = zeros(eqn.Axi)
    eqn.aux_vars = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.flux_parametric = zeros(Tsol, mesh.numDofPerNode, sbp.numnodes, 
                                mesh.numEl, Tdim)
    eqn.res = zeros(Tres, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
    eqn.res_vec = zeros(Tres, mesh.numDof)

    if opts["use_edge_res"]
      eqn.res_edge = zeros(Tres, mesh.numDofPerNode, sbp.numnodes, mesh.numEl, mesh.numTypePerElement[2])
    else
      eqn.res_edge = zeros(Tres, 0, 0, 0, 0)
    end

    eqn.q_vec = zeros(Tres, mesh.numDof)
    eqn.bndryflux = zeros(Tsol, mesh.numDofPerNode, sbp.numfacenodes, 
                          mesh.numBoundaryEdges)
    eqn.stabscale = zeros(Tres, sbp.numnodes, mesh.numInterfaces)
    

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
function calcMassMatrixInverse{Tmsh,  Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                                  sbp::SBPOperator, 
                                                  eqn::EulerData{Tsol, Tdim})
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
function calcMassMatrix{Tmsh,  Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                           sbp::SBPOperator, 
                                           eqn::EulerData{Tsol, Tdim})
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

end # end module
