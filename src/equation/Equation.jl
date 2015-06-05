module Equation

using SummationByParts
using PdePumiInterface
# make this a module?

export EulerEquation

abstract AbstractEquation

# declare a concrete type that is a subtype of Equation for every type of equation
# the concrete type can hold constants needed for the equation
# do not put any functions in this file.  All functions that take a concrete 
# type defined in this file should be in a separate file associated with the
# solver for that equation

# use this type to leverage multiple disbatch
type EulerEquation <: AbstractEquation  # hold any constants needed for euler equation, as well as solution and data needed to calculate it
# formats of all arrays are documented in SBP
# only the constants are initilized here, the arrays are not
  cv::Float64  # specific heat constant
  R::Float64  # gas constant used in ideal gas law
  gamma::Float64 # ratio of specific heats

  # the following arrays hold data for all nodes
  q::AbstractArray{Float64,3}  # holds conservative variables for all nodes
  F_xi::AbstractArray{Float64,3}  # flux in xi direction
  F_eta::AbstractArray{Float64,3} # flux in eta direction
  res::AbstractArray{Float64,3}  # result of computation

  edgestab_alpha::AbstractArray{Float64, 4} # alpha needed by edgestabilization

  Minv::AbstractArray{Float64, 1}  # invese mass matrix
  function EulerEquation(mesh::PumiMesh2, sbp::SBPOperator)
  # construct bigQ_zi and bigA_eta
  # this only works for first order

  eqn = new()  # incomplete initilization
  #=
  bigQT_xi = ones(4*operator.numnodes, 4*operator.numnodes)
  bigQT_eta = ones(4*operator.numnodes,4*operator.numnodes)

  for i=1:3
    i_i = 4*(i-1) + 1
  #   println("i_i = ", i_i)
    for j=1:3
      j_j = 4(j-1) + 1
  #     println("j_j = ", j_j)
      bigQT_xi[i_i:(i_i+3), j_j:(j_j+3)] *= operator.Q[j,i,1]
      bigQT_eta[i_i:(i_i+3), j_j:(j_j+3)] *= operator.Q[j,i,2]
    end
  end
  =#

#bigQT_xi = Array(Float64, 0, 0)
#bigQT_eta = Array(Float64, 0, 0)
eqn.gamma = 1.4
eqn.R = 287.058  # specific gas constant (unit J/(kg * K)
eqn.cv = eqn.R/(eqn.gamma - 1)

calcMassMatrixInverse(mesh, sbp, eqn)
calcEdgeStabAlpha(mesh, sbp, eqn)
# these variables get overwritten every iteration, so its safe to 
# leave them without values
eqn.q = Array(Float64, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
eqn.F_xi = Array(Float64, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
eqn.F_eta = Array(Float64, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
eqn.res = Array(Float64, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)

#println("typeof(operator.Q[1]) = ", typeof(operator.Q[1]))
#type_of_sbp = typeof(operator.Q[1])  # a little hackish
#return EulerEquation(cv, R, gamma, bigQT_xi, bigQT_eta, Array(type_of_sbp,0,0,0), Array(type_of_sbp, 0,0,0), Array(type_of_sbp, 0,0,0), Array(type_of_sbp,0,0,0))
#return EulerEquation(cv, R, gamma, bigQT_xi, bigQT_eta)

println("eq.gamma = ", eqn.gamma, " eqn.R = ", eqn.R, " eqn.cv = ", eqn.cv)
return eqn
end  # end of constructor



end  # end of type declaration


#=
function EulerEquation(operator::SBPOperator)
# construct bigQ_zi and bigA_eta
# this only works for first order

#=
bigQT_xi = ones(4*operator.numnodes, 4*operator.numnodes)
bigQT_eta = ones(4*operator.numnodes,4*operator.numnodes)

for i=1:3
  i_i = 4*(i-1) + 1
#   println("i_i = ", i_i)
  for j=1:3
    j_j = 4(j-1) + 1
#     println("j_j = ", j_j)
    bigQT_xi[i_i:(i_i+3), j_j:(j_j+3)] *= operator.Q[j,i,1]
    bigQT_eta[i_i:(i_i+3), j_j:(j_j+3)] *= operator.Q[j,i,2]
  end
end
=#

bigQT_xi = Array(Float64, 0, 0)
bigQT_eta = Array(Float64, 0, 0)
gamma = 1.4
R = 287.058  # specific gas constant (unit J/(kg * K)
cv = R/(gamma - 1)

#println("typeof(operator.Q[1]) = ", typeof(operator.Q[1]))
type_of_sbp = typeof(operator.Q[1])  # a little hackish
return EulerEquation(cv, R, gamma, bigQT_xi, bigQT_eta, Array(type_of_sbp,0,0,0), Array(type_of_sbp, 0,0,0), Array(type_of_sbp, 0,0,0), Array(type_of_sbp,0,0,0))
#return EulerEquation(cv, R, gamma, bigQT_xi, bigQT_eta)
end
=#



function calcMassMatrixInverse(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation )
# calculate the inverse mass matrix so it can be applied to the entire solution vector

  eqn.Minv = ones(Float64, mesh.numDof)

  for i=1:mesh.numEl
    dofnums_i =  getGlobalNodeNumbers(mesh, i)
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
	dofnum_k = dofnums_i[k,j]
	# multiplication is faster than division, so do the divions here
	# and then multiply solution vector times Minv
	eqn.Minv[dofnum_k] *= 1/sbp.w[j]
      end
    end
  end

  return nothing

end


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
