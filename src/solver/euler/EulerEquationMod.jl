module EulerEquationMod

include("complexify.jl")

using ArrayViews
using PDESolverCommon
using SummationByParts
using PdePumiInterface
using ForwardDiff

include(joinpath(Pkg.dir("PDESolver"), "src/tools/misc.jl"))

# the AbstractEquation type is declared in CommonTypes
# every equation will have to declare a new type that is a subtype of AbstractEquation

export AbstractEulerData, EulerData, EulerData_




# use this type to leverage multiple disbatch
# this type holds any data, including large arrays of solution values, that are specific to the equation
# for the Euler equations, this includes the conservative variables q, the fluxes in the xi and eta direction, and the result of the calculation
# the inverse mass matrix is stored here as well (not sure if that fits better here or in the mesh object)
# things like the coordinate field, the jacobian etc. are stored in the mesh objec

# the aux_vars array holds all auxiliary variables that are stored over the entire mesh
# although it can be accessed directly, the preferred method is to use the macros
# defined in the euler_macros.jl file
# these macros return either a scalar or an Array View of the specified indices, depending on if the quantity requested is scalar or vector
# every time a new variable is added to the array, the size must be updated
# and a new macro must be created.
# The uses of aux_vars should mirror that of eqn.q, in that entire columns should be
# passed to low level functions and the low level functions use the macros to access
# individual variables
# the advantages of macros vs functions for access to variables remains unclear
# if aux_vars is a fixed size
# if it is variable sized then macros give the advantage of doing location lookup
# at compile time

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
immutable ParamType{Tdim} 
  cv::Float64  # specific heat constant
  R::Float64  # specific gas constant used in ideal gas law (J/(Kg * K))
  gamma::Float64 # ratio of specific heats
  gamma_1::Float64 # = gamma - 1

  Ma::Float64  # free stream Mach number
  Re::Float64  # free stream Reynolds number
  aoa::Float64  # angle of attack
  rho_free::Float64  # free stream density
  E_free::Float64 # free stream energy (4th conservative variable)

  # debugging options
  writeflux::Bool  # write Euler flux
  writeboundary::Bool  # write boundary data
  writeq::Bool # write solution variables
  function ParamType(opts)
  # create values, apply defaults

    # get() = get(dictionary, key, default)
    gamma = opts[ "gamma"]
    gamma_1 = gamma - 1
    R = opts[ "R"]
    cv = R/gamma_1

    Ma = opts[ "Ma"]
    Re = opts[ "Re"]
    aoa = opts[ "aoa"]
    rho_free = opts[ "rho_free"]
    E_free = opts[ "E_free"]


    # debugging options
    writeflux = opts[ "writeflux"]
    writeboundary = opts[ "writeboundary"]
    writeq = opts["writeq"]

    return new(cv, R, gamma, gamma_1, Ma, Re, aoa, rho_free, E_free, writeflux, writeboundary, writeq)

  end

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
    * F_xi : 4D array [ndof per node, nnodes per element, nelements, Tdim]
             holding the Euler flux in the xi and eta directions
    * res  : 3D array holding residual
    * SL   : vector form of res
    * SL0  : initial condition vector
    * edgestab_alpha : paramater used for edge stabilization, 4d array
    * bndryflux : 3D array holding boundary flux data
    * stabscale : 2D array holding edge stabilization scale factor
    * Minv :  vector holding inverse mass matrix
"""->
abstract EulerData {Tsol, Tdim} <: AbstractEulerData{Tsol}

# high level functions should take in an AbstractEulerData, remaining
# agnostic to the dimensionality of the equation
# Mid level function should take in an EulerData{Tsol, Tdim}, so they
# know the dimensionality of the equation, but have a single method that
# can handle all dimensions

# low level functions should take in EulerData{Tsol, 2} or EulerData{Tsol, 3}
# this allows them to have different methods for different dimension equations.

@doc """
### EulerEquationMod.EulerData_

  This type is an implimentation of the abstract EulerData.  It is
  paramterized by the residual type Tres and the mesh type Tmsh
  because it stores some arrays of those types.  Tres is the 'maximum' type of
  Tsol and Tmsh, where Tsol is the type of the conservative variables.

  Eventually there will be additional implimentation of EulerData,
  specifically a 3D one.

"""->
type EulerData_{Tsol, Tres, Tdim, Tmsh} <: EulerData{Tsol, Tdim}  # hold any constants needed for euler equation, as well as solution and data needed to calculate it
# formats of all arrays are documented in SBP
# only the constants are initilized here, the arrays are not

  params::ParamType{Tdim}
  res_type::DataType  # type of res

  # the following arrays hold data for all nodes
  q::Array{Tsol,3}  # holds conservative variables for all nodes
  # hold fluxes in all directions
  # [ndof per node by nnodes per element by num element by num dimensions]
  aux_vars::Array{Tsol, 3}  # storage for auxiliary variables 
  F_xi::Array{Tsol,4}  # flux in xi direction
#  F_eta::Array{Tsol,3} # flux in eta direction
  res::Array{Tres, 3}  # result of computation
  SL::Array{Tres, 1}  # result of computation in vector form
  SL0::Array{Tres,1}  # initial condition in vector form


  edgestab_alpha::Array{Tmsh, 4} # alpha needed by edgestabilization
  bndryflux::Array{Tsol, 3}  # boundary flux
  stabscale::Array{Tsol, 2}  # stabilization scale factor

  Minv::Array{Float64, 1}  # invese mass matrix

  # inner constructor
#  function EulerData(mesh::PumiMesh2, sbp::SBPOperator, T2::DataType)
  function EulerData_(mesh::PumiMesh2, sbp::SBPOperator, opts)

    eqn = new()  # incomplete initilization

    eqn.params = ParamType{Tdim}(opts)
#=
    eqn.gamma = 1.4
    eqn.gamma_1 = eqn.gamma - 1
    eqn.R = 287.058  # specific gas constant (unit J/(kg * K)
    eqn.cv = eqn.R/(eqn.gamma - 1)
=#
    eqn.res_type = Tres

    calcMassMatrixInverse(mesh, sbp, eqn)
    calcEdgeStabAlpha(mesh, sbp, eqn)
    
    # must initialize them because some datatypes (BigFloat) 
    # don't automatically initializes them
    # taking a view(A,...) of undefined values is illegal
    # I think its a bug that Array(Float64, ...) initailizes values
    eqn.q = zeros(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
    eqn.aux_vars = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.F_xi = zeros(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl, Tdim)
#    eqn.F_eta = Array(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
  #  eqn.res = Array(T2, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
    eqn.res = zeros(Tres, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
    eqn.SL = zeros(Tres, mesh.numDof)
    eqn.SL0 = zeros(Tres, mesh.numDof)

    eqn.bndryflux = zeros(Tsol, mesh.numDofPerNode, sbp.numfacenodes, mesh.numBoundaryEdges)
    eqn.stabscale = zeros(Tres, sbp.numnodes, mesh.numInterfaces)

    #println("typeof(operator.Q[1]) = ", typeof(operator.Q[1]))
    #type_of_sbp = typeof(operator.Q[1])  # a little hackish
    #return EulerData(cv, R, gamma, bigQT_xi, bigQT_eta, Array(type_of_sbp,0,0,0), Array(type_of_sbp, 0,0,0), Array(type_of_sbp, 0,0,0), Array(type_of_sbp,0,0,0))
    #return EulerData(cv, R, gamma, bigQT_xi, bigQT_eta)

    println("eqn.res_type = ", eqn.res_type)
    return eqn
  end  # end of constructor

end  # end of type declaration




# now that EulerData is defined, include other files that use it

include("euler_macros.jl")
include("output.jl")
include("common_funcs.jl")
#include("sbp_interface.jl")
include("euler.jl")
include("ic.jl")
include("bc.jl")
include("stabilization.jl")
#include("constant_diff.jl")




@doc """
### EulerEquationMod.calcMassMatrixInverse

  This function calculates the inverse mass matrix and stores it in eqn.Minv.
  Because we use SBP operators, the mass matrix is diagonal, so it is stored
  in a vector.  mesh.dofs is used to put the components of the inverse
  mass matrix in the same place as the corresponding values in eqn.SL


"""->
# used by EulerData Constructor
# mid level functions
function calcMassMatrixInverse{Tmsh, Tsbp, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerData{Tsol, Tdim} )
# calculate the inverse mass matrix so it can be applied to the entire solution vector
# mass matrix is diagonal, stores in vector eqn.Minv

  eqn.Minv = zeros(Tmsh, mesh.numDof)

  for i=1:mesh.numEl
#    dofnums_i =  getGlobalNodeNumbers(mesh, i)
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
#	dofnum_k = dofnums_i[k,j]
        dofnum_k = mesh.dofs[k,j,i]
	# multiplication is faster than division, so do the divions here
	# and then multiply solution vector times Minv
	eqn.Minv[dofnum_k] += 1/(sbp.w[j]*mesh.jac[j,i])

#	eqn.Minv[dofnum_k] *= 1/(sbp.w[j])
      end
    end
  end

  return nothing

end


end # end module
