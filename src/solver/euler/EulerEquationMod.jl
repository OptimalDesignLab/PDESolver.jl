module EulerEquationMod

using ArrayViews
using PDESolverCommon
using SummationByParts
using PdePumiInterface
using ForwardDiff
# the AbstractEquation type is declared in CommonTypes
# every equation will have to declare a new type that is a subtype of AbstractEquation

export AbstractEulerEquation, EulerEquation, ConcreteEulerEquation




# use this type to leverage multiple disbatch
# this type holds any data, including large arrays of solution values, that are specific to the equation
# for the Euler equations, this includes the conservative variables q, the fluxes in the xi and eta direction, and the result of the calculation
# the inverse mass matrix is stored here as well (not sure if that fits better here or in the mesh object)
# things like the coordinate field, the jacobian etc. are stored in the mesh objec

abstract AbstractEulerEquation{Tsol} <: AbstractEquation{Tsol}
abstract EulerEquation{Tsol, Tdim} <: AbstractEulerEquation{Tsol}

# use AbstractEulerEquation for high level functions
# use EulerEquation for mid level functions
# don't use ConcreteEulerEquation (unless you really need to dispatch based on Tres)


type ConcreteEulerEquation{Tsol, Tres, Tdim} <: EulerEquation{Tsol, Tdim}  # hold any constants needed for euler equation, as well as solution and data needed to calculate it
# formats of all arrays are documented in SBP
# only the constants are initilized here, the arrays are not
# Tsol is solution conservative variable data type, Tres is solution data type
# Tdim is dimensionality of the equation
  cv::Float64  # specific heat constant
  R::Float64  # gas constant used in ideal gas law
  gamma::Float64 # ratio of specific heats
  res_type::DataType  # type of res

  # the following arrays hold data for all nodes
  q::Array{Tsol,3}  # holds conservative variables for all nodes
  # flux in all directions
  # [ndof per node by nnodes per element by numelements by num dimensions]
  F_xi::Array{Tsol,4}
#  F_eta::Array{Tsol,3} # flux in eta direction
  res::Array{Tres, 3}  # result of computation
  SL::Array{Tres, 1}  # result of computation in vector form
  SL0::Array{Tres,1}  # initial condition in vector form


  edgestab_alpha::Array{Float64, 4} # alpha needed by edgestabilization
  bndryflux::Array{Tsol, 3}  # boundary flux
  stabscale::Array{Tsol, 2}  # stabilization scale factor

  Minv::Array{Float64, 1}  # invese mass matrix

  # inner constructor
#  function EulerEquation(mesh::PumiMesh2, sbp::SBPOperator, T2::DataType)
  function ConcreteEulerEquation(mesh::PumiMesh2, sbp::SBPOperator)

    eqn = new()  # incomplete initilization

    eqn.gamma = 1.4
    eqn.R = 287.058  # specific gas constant (unit J/(kg * K)
    eqn.cv = eqn.R/(eqn.gamma - 1)
    eqn.res_type = Tres

    calcMassMatrixInverse(mesh, sbp, eqn)
    calcEdgeStabAlpha(mesh, sbp, eqn)
    # these variables get overwritten every iteration, so its safe to 
    # leave them without values
    eqn.q = Array(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
    eqn.F_xi = Array(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl, Tdim)
#    eqn.F_eta = Array(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
  #  eqn.res = Array(T2, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
    eqn.res = Array(Tres, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
    eqn.SL = Array(Tres, mesh.numDof)
    eqn.SL0 = Array(Tres, mesh.numDof)

    eqn.bndryflux = Array(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numBoundaryEdges)
    eqn.stabscale = Array(Tsol, sbp.numnodes, mesh.numInterfaces)

    #println("typeof(operator.Q[1]) = ", typeof(operator.Q[1]))
    #type_of_sbp = typeof(operator.Q[1])  # a little hackish
    #return EulerEquation(cv, R, gamma, bigQT_xi, bigQT_eta, Array(type_of_sbp,0,0,0), Array(type_of_sbp, 0,0,0), Array(type_of_sbp, 0,0,0), Array(type_of_sbp,0,0,0))
    #return EulerEquation(cv, R, gamma, bigQT_xi, bigQT_eta)

    println("eq.gamma = ", eqn.gamma, " eqn.R = ", eqn.R, " eqn.cv = ", eqn.cv)
    println("eqn.res_type = ", eqn.res_type)
    return eqn
  end  # end of constructor

end  # end of type declaration


# now that EulerEquation is defined, include other files that use it
include("common_funcs.jl")
include("sbp_interface.jl")
include("euler.jl")
include("ic.jl")
include("bc.jl")
include("stabilization.jl")




# used by EulerEquation Constructor
function calcMassMatrixInverse(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation )
# calculate the inverse mass matrix so it can be applied to the entire solution vector
# mass matrix is diagonal, stores in vector eqn.Minv

  eqn.Minv = ones(Float64, mesh.numDof)

  for i=1:mesh.numEl
    dofnums_i =  getGlobalNodeNumbers(mesh, i)
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
	dofnum_k = dofnums_i[k,j]
	# multiplication is faster than division, so do the divions here
	# and then multiply solution vector times Minv
	eqn.Minv[dofnum_k] *= 1/(sbp.w[j]*mesh.jac[j,i])
      end
    end
  end

  return nothing

end

# used by EulerEQuation Constructor
function calcEdgeStabAlpha(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation)
# calculate alpha, needed by edge stabilization

  numEl = getNumEl(mesh)

  eqn.edgestab_alpha = Array(Float64,2,2,sbp.numnodes,numEl)
  dxidx = mesh.dxidx
  jac = mesh.jac

  # calculating alpha, required by edgestabilize!
  # this canbe compuated apriori
  for k = 1:numEl
    for i = 1:sbp.numnodes
      for di1 = 1:2
        for di2 = 1:2
          eqn.edgestab_alpha[di1,di2,i,k] = (dxidx[di1,1,i,k].*dxidx[di2,1,i,k] + dxidx[di1,2,i,k].*dxidx[di2,2,i,k])*jac[i,k]
        end
      end
    end
  end

  return nothing
end




end # end module
