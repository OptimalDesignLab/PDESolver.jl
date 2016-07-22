@doc """
### EulerEquationMod.applySourceTerm

  This function updates eqn.res with the source term.  

  Inputs: 
    mesh
    sbp
    eqn
    opts
    src_func:  the functor that returns the value of the source term at a node
               This functor must have the signature:
               src_func(q, coords, params, t)
               where coords is a vector of length 2 or 3 containing the x and y 
               coordinates of the node, params is the ParamType, t is the 
               current time, and q is the vector to be populated with the 
               source term values.

  Outputs: none

  Aliasing restrictions: params.q_vals cannot be in use

"""->
function applySourceTerm(mesh,sbp, eqn, opts, src_func::SRCType)
#TODO: check that the k loop vectorizes
  weights = sbp.w
  q_vals = eqn.params.q_vals
  t = eqn.params.t

  for i=1:mesh.numEl
    jac_i = sview(mesh.jac, :, i)
    res_i = sview(eqn.res, :, :, i)
    for j=1:mesh.numNodesPerElement
      coords_j = sview(mesh.coords, :, j, i)
      src_func(q_vals, coords_j, eqn.params, t)
      fac = weights[j]/jac_i[j]
      for k=1:mesh.numDofPerNode
        res_i[k, j] += fac*q_vals[k]
      end
    end
  end

  return nothing
end





type SRC0 <: SRCType  # dummy source functor, it should nevery actually be called
end

@doc """
### AdvectionEquationMod.SRC0

  This is the zero source term.  This is the default of source term
  is specified
"""->
type SRCExp <: SRCType
end

function call(obj::SRCExp, q::AbstractVector, coords::AbstractVector, params::ParamType{2}, t)
  x = coords[1]
  y = coords[2]
  gamma_1 = params.gamma_1
  # a and b are parameters determine the scale and offset of the solution
  # the must be the same here and in the boundary condition
  af = 1/5  # a = 1/af
  b = 0.01

  q[1] = 2*y*af*exp(2*x*y*af + b) + 3*x*af*exp(3*x*y*af + b)
  q[2] = 3*y*af*exp(3*x*y*af + b) + 5*y*af*exp(5*x*y*af + b) + 4*x*af*exp(4*x*y*af+b)
  q[3] = 4*y*af*exp(4*x*y*af + b) + 10*x*af*exp(5*x*y*af + b)
  q[4] = 6*y*af*(1/gamma_1 + 1.5)*exp(6*x*y*af + b) + 2*y*af*exp(4*x*y*af + b) + 7*x*af*(1/gamma_1 + 1.5)*exp(7*x*y*af + b) + 2.5*x*af*exp(5*x*y*af + b)

  # make this work for other variables
  convertFromNaturalToWorkingVars(params, q, q)
end

function call(obj::SRCExp, q::AbstractVector, coords::AbstractVector, params::ParamType{3}, t)

  x = coords[1]
  y = coords[2]
  z = coords[3]
  gamma_1 = params.gamma_1
  a = 1/500
  b = 0.01
  c = 20
  d = 0.25
  d3 = d*d*d

#  c1 = c + 1
#  c2 = c + 2
#  c3 = c + 3

  t2 = a*x*y*z*2.0;
  t3 = exp(t2);
  t4 = exp(b);
  t5 = d*d;
  t6 = a*x*y*z*5.0;
  t7 = exp(t6);
  t8 = a*x*y*z*4.0;
  t9 = exp(t8);
  t10 = a*c*x*y*z;
  t11 = exp(t10);
  t12 = a*x*y*z*6.0;
  t13 = exp(t12);
  t14 = a*x*y*z;
  t15 = exp(t14);
  t16 = a*x*y*z*7.0;
  t17 = exp(t16);
  t18 = a*x*y*z*3.0;
  t19 = exp(t18);
  q[1] = a*d*t3*t4*(y*z*2.0+t3*x*y*4.0+t15*x*z*3.0)
  q[2] = a*t4*(c*t11*y*z+t5*t7*x*y*5.0+t5*t9*x*z*4.0+t5*t19*y*z*3.0)
  q[3] = a*t4*(c*t11*x*z+t5*t13*x*y*6.0+t5*t7*x*z*5.0+t5*t9*y*z*4.0)
  q[4] = a*t4*(c*t11*x*y+t5*t17*x*y*7.0+t5*t13*x*z*6.0+t5*t13*y*z*6.0)
  q[5] = a*d*t4*t15*( (t11*y*z*2.0+c*t11*y*z*2.0+t5*t7*x*y*6.0+t3*t11*x*y*6.0) + (t5*t17*x*y*8.0+t5*t9*x*z*5.0+t5*t13*x*z*7.0+t11*t15*x*z*4.0) + (t5*t7*y*z*6.0+t5*t17*y*z*8.0+t5*t19*y*z*4.0+t5*x*y*exp(a*x*y*z*9.0)*1.0E1) + t5*x*z*exp(a*x*y*z*8.0)*9.0+c*t3*t11*x*y*2.0+c*t11*t15*x*z*2.0)*(1.0/2.0)+(a*d*t4*t11*t15*(y*z+c*y*z+t3*x*y*3.0+t15*x*z*2.0+c*t3*x*y+c*t15*x*z))/gamma_1

#=
  q[1] = d*2*af*y*z*exp(2*af*x*y*z + b) + d*3*af*x*z*exp(3*af*x*y*z + b) + d*4*af*x*y*exp(4*af*x*y*z + b)

  q[2] = d*d*3*af*y*z*exp(3*af*x*y*z + b) + c*af*y*z*exp(c*af*x*y*z + b) + d*d*4*af*x*z*exp(4*af*x*y*z + b) + d*d*5*af*x*y*exp(5*af*x*y*z + b)

  q[3] = d*d*4*af*y*z*exp(4*af*x*y*z + b) + d*d*5*af*x*z*exp(5*af*x*y*z + b) + c*af*x*z*exp(c*af*x*y*z + b) + d*d*6*af*x*y*exp(6*af*x*y*z + b)

  q[4] = d*d*5*af*y*z*exp(5*af*x*y*z + b) + d*d*6*af*x*z*exp(6*af*x*y*z + b) + d*d*7*af*x*y*exp(7*af*x*y*z + b) + c*af*x*y*exp(c*af*x*y*z + b)
  
  q[5] = ( d*c1*af*y*z*(1/gamma_1 + 1)*exp(c1*af*x*y*z + b) + d3*2*af*y*z*exp(4*af*x*y*z + b) + d3*3*af*y*z*exp(6*af*x*y*z + b) + d3*4*af*y*z*exp(8*af*x*y*z + b) ) +
         ( d*c2*af*x*z*(1/gamma_1 + 1)*exp(c2*af*x*y*z + b) + d3*2.5*af*x*z*exp(5*af*x*y*z + b) + d3*3.5*af*x*z*exp(7*af*x*y*z + b) + d3*4.5*af*x*z*exp(9*af*x*y*z + b) ) +
         ( d*c3*af*x*y*(1/gamma_1 + 1)*exp(c3*af*x*y*z + b) + d3*3*af*x*y*z + b) + d3*3*af*x*y*exp(6*af*x*y*z + b) + d3*4*af*x*y*exp(8*af*x*y*z + b) + d3*5*af*x*y*exp(10*af*x*y*z + b)
=#
  return nothing
end

@doc """
### EulerEquationMod.SRCDict

  Stores the functors that evaluate source terms at a node.  Every new 
  functor should be added to this dictonary

  All functors must have the signature:

  src_func(q, coords, params::ParamType, t)

  where coords is the vector of length 2 containing the x and y coordinates
  of the node, t is the current time, and q is the vector to be populated with
  the source term values.
"""->
global const SRCDict = Dict{ASCIIString, SRCType}(
"SRCExp" => SRCExp(),
"SRC0" => SRC0(),
)


@doc """
### EulerEquationMod.getSRCFunctors

  This function gets the functor specified by opts["SRCname"] and stores
  it to the equation object.  Currently one 1 source functor is allowed.

"""->
function getSRCFunctors(mesh::AbstractMesh, sbp::AbstractSBP, 
                        eqn::EulerData, opts)

  # currently we only allow 1 source functor
  eqn.src_func = SRCDict[opts["SRCname"]]
  println("using source term functor ", eqn.src_func)
  return nothing
end


