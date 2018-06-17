#
#TODO: As far as we know `q`, `q_x`, `q_xx`, we can move
# the rest into a function, which is less error prone, and 
# will save a lot of space
#
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
    jac_i = ro_sview(mesh.jac, :, i)
    res_i = sview(eqn.res, :, :, i)
    for j=1:mesh.numNodesPerElement
      coords_j = ro_sview(mesh.coords, :, j, i)
      src_func(q_vals, coords_j, eqn.params, t)
      fac = weights[j]/jac_i[j]
      for k=1:mesh.numDofPerNode
        res_i[k, j] += fac*q_vals[k]
      end
    end
  end

  return nothing
end





mutable struct SRC0 <: SRCType  # dummy source functor, it should nevery actually be called
end

@doc """
### AdvectionEquationMod.SRC0

  This is the zero source term.  This is the default of source term
  is specified
"""->
mutable struct SRCExp <: SRCType
end

function (obj::SRCExp)(q::AbstractVector, coords::AbstractVector, params::ParamType{2}, t)
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

# declare some global variables
# this is somewhat namespace polluting
global const MMSExp_a = 1/500
global const MMSExp_b = 0.01
global const MMSExp_c1 = 1
global const MMSExp_c2 = 2
global const MMSExp_c3 = 3
global const MMSExp_c4 = 4
global const MMSExp_c5 = 20
global const MMSExp_d1 = 1
global const MMSExp_d2 = 0.05
global const MMSExp_d3 = 0.15
global const MMSExp_d4 = 0.25
global const MMSExp_d5 = 1

function (obj::SRCExp)(q::AbstractVector, coords::AbstractVector, params::ParamType{3}, t)

  x = coords[1]
  y = coords[2]
  z = coords[3]

  # constant parameters
  gamma_1 = params.gamma_1
  a = MMSExp_a
  b = MMSExp_b
  c1 = MMSExp_c1 
  c2 = MMSExp_c2
  c3 = MMSExp_c3
  c4 = MMSExp_c4
  c5 = MMSExp_c5
  d1 = MMSExp_d1
  d2 = MMSExp_d2
  d3 = MMSExp_d3
  d4 = MMSExp_d4
  d5 = MMSExp_d5

  t2 = exp(b);
  t3 = a*c2*x*y*z;
  t4 = exp(t3);
  t5 = a*c4*x*y*z;
  t6 = exp(t5);
  t7 = c4*d4*t6*x*y;
  t8 = a*c3*x*y*z;
  t9 = exp(t8);
  t10 = c3*d3*t9*x*z;
  t11 = a*c5*x*y*z;
  t12 = exp(t11);
  t13 = 1.0/d1;
  t16 = a*c1*x*y*z;
  t14 = exp(-t16);
  t15 = c2*d2*t4*y*z;
  t17 = b-t16;
  t18 = exp(t17);
  t19 = d2*d2;
  t20 = a*c2*x*y*z*2.0;
  t21 = exp(t20);
  t22 = d3*d3;
  t23 = a*c3*x*y*z*2.0;
  t24 = exp(t23);
  t25 = d4*d4;
  t26 = a*c4*x*y*z*2.0;
  t27 = exp(t26);
  t28 = b+t11;
  t29 = exp(t28);
  t30 = c2*t19*t21;
  t31 = c3*t22*t24;
  t32 = c4*t25*t27;
  t33 = t30+t31+t32;
  t34 = t19*t21;
  t35 = t22*t24;
  t36 = t25*t27;
  t37 = t34+t35+t36;
  t38 = 1.0/gamma_1;
  t39 = c1-c4;
  t40 = exp(-a*t39*x*y*z);
  t41 = c1-c3;
  t42 = exp(-a*t41*x*y*z);
  t43 = d5*t29;
  t44 = d5*t29*t38;
  t45 = t13*t18*t37*(1.0/2.0);
  t46 = t43+t44+t45;
  t47 = c1-c2;
  t48 = exp(-a*t47*x*y*z);
  q[1] = a*t2*(t7+t10+t15);

  q[2] = a*c5*d5*t2*t12*y*z + a*d2*t2*t4*t13*t14*(t7+t10 - c1*d4*t6*x*y + c2*d4*t6*x*y - c1*d3*t9*x*z + c2*d3*t9*x*z - c1*d2*t4*y*z + c2*d2*t4*y*z*2.0);

  q[3] = a*c5*d5*t2*t12*x*z + a*d3*t2*t9*t13*t14*( (t7+t15 - c1*d4*t6*x*y + c3*d4*t6*x*y - c1*d3*t9*x*z) + (c3*d3*t9*x*z*2.0 - c1*d2*t4*y*z + c3*d2*t4*y*z) );

  q[4] = a*c5*d5*t2*t12*x*y + a*d4*t2*t6*t13*t14*( (t10+t15 - c1*d4*t6*x*y + c4*d4*t6*x*y*2.0 - c1*d3*t9*x*z) + (c4*d3*t9*x*z - c1*d2*t4*y*z + c4*d2*t4*y*z) );

  q[5] = d4*t13*t40*(a*t13*t18*t33*x*y + a*c5*d5*t29*x*y + a*c5*d5*t29*t38*x*y - a*c1*t13*t18*t37*x*y*(1.0/2.0)) + d3*t13*t42*(a*t13*t18*t33*x*z + a*c5*d5*t29*x*z + a*c5*d5*t29*t38*x*z - a*c1*t13*t18*t37*x*z*(1.0/2.0)) + d2*t13*t48*(a*t13*t18*t33*y*z + a*c5*d5*t29*y*z + a*c5*d5*t29*t38*y*z - a*c1*t13*t18*t37*y*z*(1.0/2.0)) - a*d4*t13*t39*t40*t46*x*y - a*d3*t13*t41*t42*t46*x*z - a*d2*t13*t46*t47*t48*y*z;

  return nothing
end

"""
  Functor for source term corresponding to ICPeriodicMMS
"""
mutable struct SRCPeriodicMMS <: SRCType
end

function (obj::SRCPeriodicMMS)(q::AbstractVector, coords::AbstractVector, 
              params::ParamType{2}, t)

  x = coords[1]
  y = coords[2]
  gamma_1 = params.gamma_1

  t4 = t*2.0;
  t2 = -t4+x+y;
  t3 = 3.141592653589793*t2;
  t5 = cos(t3);
  t6 = sin(t3);
  t7 = t6+1.5E1;
  t8 = 3.141592653589793*gamma_1*t5*t7*(1.0/5.0E1);
  q[1] = 0
  q[2] = t8;
  q[3] = t8;
  q[4] = 3.141592653589793*gamma_1*t5*t7*(1.0/2.5E1);

  return nothing
end

function (obj::SRCPeriodicMMS)(q::AbstractVector, coords::AbstractVector, 
              params::ParamType{3}, t)

  x = coords[1]
  y = coords[2]
  z = coords[3]
  gamma_1 = params.gamma_1
  gamma = params.gamma
#=
  t2 = x+y+z;
  t3 = 3.141592653589793*t2;
  t4 = cos(t3);
  t5 = sin(t3);
  t6 = t5*(1.0/1.0E1);
  t7 = t6+2.0;
  t8 = 3.141592653589793*gamma_1*t4*t7*(1.0/5.0);
=#
#  tmp = pi*gamma_1*cos(pi*(x + y + z))*(sin(pi*(x + y + z))/10 + 2)/5
#=
  q[1] = 0;
  q[2] = 0
  q[3] = 0
  q[4] = 0
  q[5] = 0;
=#

  t5 = t*2.0;
  t2 = -t5+x+y+z;
  t3 = 3.141592653589793*t2;
  t4 = cos(t3);
  t6 = gamma_1*2.5E1;
  t7 = sin(t3);
  t8 = gamma_1*t7*2.0;
  t9 = t6+t8+1.0E1;
  t10 = 3.141592653589793*t4*t9*(1.0/1.0E2);
  q[1] = 3.141592653589793*t4*(1.0/1.0E1);
  q[2] = t10;
  q[3] = t10;
  q[4] = t10;
  q[5] = 3.141592653589793*t4*(gamma_1*7.5E1+t7*2.0+gamma_1*t7*6.0+4.0E1)*(1.0/1.0E2);

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
global const SRCDict = Dict{String, SRCType}(
"SRCExp" => SRCExp(),
"SRCPeriodicMMS" => SRCPeriodicMMS(),
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


