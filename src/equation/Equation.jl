# make this a module?


abstract AbstractEquation

# declare a concrete type that is a subtype of Equation for every type of equation
# the concrete type can hold constants needed for the equation
# do not put any functions in this file.  All functions that take a concrete 
# type defined in this file should be in a separate file associated with the
# solver for that equation

# use this type to leverage multiple disbatch
type EulerEquation <: AbstractEquation  # hold any constants needed for euler equation
  cv::Float64  # specific heat constant
  R::Float64  # gas constant used in ideal gas law
  gamma::Float64 # ratio of specific heats
  bigQT_xi::Array{Float64, 2}  # big parent element stiffness matrix transposed (Q_zi double under bar transposed
  bigQT_eta::Array{Float64,2}  # big parent element stiffness matrix transposed (Q_eta double under  bar transposed

  # use default inner constructor
end

function EulerEquation(operator::SBPOperator)
# construct bigQ_zi and bigA_eta
# this only works for first order

bigQT_xi = ones(4*operator.numnodes, 4*operator.numnodes)
bigQT_eta = ones(4*operator.numnodes,4*operator.numnodes)

for i=1:3
  i_i = 4*(i-1) + 1
  println("i_i = ", i_i)
  for j=1:3
    j_j = 4(j-1) + 1
    println("j_j = ", j_j)
    bigQT_xi[i_i:(i_i+3), j_j:(j_j+3)] *= operator.Q[j,i,1]
    bigQT_eta[i_i:(i_i+3), j_j:(j_j+3)] *= operator.Q[j,i,2]
  end
end

gamma = 1.4
R = 287.058  # specific gas constant (unit J/(kg * K)
cv = R/(gamma - 1)

return EulerEquation(cv, R, gamma, bigQT_xi, bigQT_eta)

end
