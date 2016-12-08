# eigenvectors for flux jacobians and scaling vectors that make
# dq/dw = Y*Y.', where Y are the eigenvectors of the conservative variable
# flux jacobian.

# Most of the function bodies were generated from the Matlab scripts 
# Euler_eigenscaling.m and Euler_eigenscaling2d.m in the doc directory
# The scripts are substantially more readable than the generated code

#------------------------------------------------------------------------------
# Flux jacobians

"""
  This function calculate the Euler flux jacobian in the x direction from 
  the conservative variables.  Methods are available for 2 and 3 dimensions

  Note that this code calculates pressure internally (ie. it does not call
  calcPressure).

  Inputs:
    params: a ParamType
    q: a vector of conservative variables

  Inputs/Outputs:
    A: matrix to be populated with the flux jacobian (overwritten)

  Aliasing restrictions: none
"""
function calcA(params::ParamType{2, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q5 = q[4]
  gami = params.gamma_1
  gamma = params.gamma


  # auto-generated code
  t2 = q2*q2;
  t3 = 1.0/q1;
  t4 = 1.0/(q1*q1);
  t5 = gami*t2;
  t6 = q3*q3;
  t7 = gami*t6;

  R[1,1] = 0;
  R[2,1] = t4*(t2*-2.0+t5+t7)*0.5;
  R[3,1] = -q2*q3*t4;
  R[4,1] = t3*t4*q2*(t5+t7-q1*q5-gami*q1*q5);

  R[1,2] = 1.0;
  R[2,2] = -q2*t3*(gami-2.0);
  R[3,2] = q3*t3;
  R[4,2] = t4*(t7+gami*t2*3.0-q1*q5*2.0-gami*q1*q5*2.0)*-0.5;

  R[1,3] = 0;
  R[2,3] = -gami*q3*t3;
  R[3,3] = q2*t3;
  R[4,3] = -gami*q2*q3*t4;

  R[1,4] = 0;
  R[2,4] = gami;
  R[3,4] = 0;
  R[4,4] = q2*t3*(gamma);


  return nothing
end


function calcA(params::ParamType{3, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = q2*q2;
  t3 = 1.0/q1;
  t4 = t3*t3
  t5 = q2*t3;
  t6 = q3*q3;
  t7 = q4*q4;
  t8 = gami*t6;
  t9 = gami*t7;
  R[1,1] = 0;
  R[1,2] = 1.0;
  R[1,3] = 0;
  R[1,4] = 0;
  R[1,5] = 0;
  R[2,1] = t4*(t2*-2.0+t8+t9+gami*t2)*0.5;
  R[2,2] = t3*(q2*2.0-gami*q2);
  R[2,3] = -gami*q3*t3;
  R[2,4] = -gami*q4*t3;
  R[2,5] = gami;
  R[3,1] = -q2*q3*t4;
  R[3,2] = q3*t3;
  R[3,3] = t5;
  R[3,4] = 0;
  R[3,5] = 0;
  R[4,1] = -q2*q4*t4;
  R[4,2] = q4*t3;
  R[4,3] = 0;
  R[4,4] = t5;
  R[4,5] = 0;
  R[5,1] = 1.0/(q1*q1*q1)*q2*(gami*t2*2.0+gami*t6*2.0+gami*t7*2.0-q1*q5*2.0-gami*q1*q5*2.0)*0.5;
  R[5,2] = t4*(t8+t9+gami*t2*3.0-q1*q5*2.0-gami*q1*q5*2.0)*-0.5;
  R[5,3] = -gami*q2*q3*t4;
  R[5,4] = -gami*q2*q4*t4;
  R[5,5] = t3*(q2+gami*q2);

  return nothing
end
"""
  Like calcA, but in the y directions.  See that function for details
"""
function calcB(params::ParamType{2, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q5 = q[4]
  gami = params.gamma_1
  gamma = params.gamma

  t2 = 1.0/q1;
  t3 = t2*t2;
  t4 = q3*q3;
  t5 = q2*q2;
  t6 = gami*t5;
  t7 = gami*t4;

  R[1,1] = 0;
  R[2,1] = -q2*q3*t3;
  R[3,1] = t3*(t4*-2.0+t6+t7)*0.5;
  R[4,1] = t2*t3*q3*(t6+t7-q1*q5-gami*q1*q5);

  R[1,2] = 0;
  R[2,2] = q3*t2;
  R[3,2] = -gami*q2*t2;
  R[4,2] = -gami*q2*q3*t3;

  R[1,3] = 1.0;
  R[2,3] = q2*t2;
  R[3,3] = t2*(q3*2.0-gami*q3);
  R[4,3] = t3*(t6+gami*t4*3.0-q1*q5*2.0-gami*q1*q5*2.0)*-0.5;

  R[1,4] = 0;
  R[2,4] = 0;
  R[3,4] = gami;
  R[4,4] = q3*t2*(gamma);

  return nothing
end


function calcB(params::ParamType{3, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = 1.0/q1;
  t3 = 1.0/(q1*q1);
  t4 = q3*q3;
  t5 = q3*t2;
  t6 = q2*q2;
  t7 = gami*t6;
  t8 = gami*t4;
  t9 = q4*q4;
  t10 = gami*t9;
  R[1,1] = 0;
  R[1,2] = 0;
  R[1,3] = 1.0;
  R[1,4] = 0;
  R[1,5] = 0;
  R[2,1] = -q2*q3*t3;
  R[2,2] = t5;
  R[2,3] = q2*t2;
  R[2,4] = 0;
  R[2,5] = 0;
  R[3,1] = t3*(t4*-2.0+t7+t8+t10)*0.5;
  R[3,2] = -gami*q2*t2;
  R[3,3] = -q3*t2*(gami-2.0);
  R[3,4] = -gami*q4*t2;
  R[3,5] = gami;
  R[4,1] = -q3*q4*t3;
  R[4,2] = 0;
  R[4,3] = q4*t2;
  R[4,4] = t5;
  R[4,5] = 0;
  R[5,1] = 1.0/(q1*q1*q1)*q3*(t7+t8+t10-q1*q5-gami*q1*q5);
  R[5,2] = -gami*q2*q3*t3;
  R[5,3] = t3*(t7+t10+gami*t4*3.0-q1*q5*2.0-gami*q1*q5*2.0)*-0.5;
  R[5,4] = -gami*q3*q4*t3;
  R[5,5] = q3*t2*(gamma);

  return nothing
end

"""
  Like calcA, but in the z direction.  See that function for details
"""
function calcC(params::ParamType{3, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  t2 = 1.0/q1;
  t3 = 1.0/(q1*q1);
  t4 = q4*t2;
  t5 = q4*q4;
  t6 = q2*q2;
  t7 = gami*t6;
  t8 = q3*q3;
  t9 = gami*t8;
  t10 = gami*t5;
  R[1,1] = 0;
  R[1,2] = 0;
  R[1,3] = 0;
  R[1,4] = 1.0;
  R[1,5] = 0;
  R[2,1] = -q2*q4*t3;
  R[2,2] = t4;
  R[2,3] = 0;
  R[2,4] = q2*t2;
  R[2,5] = 0;
  R[3,1] = -q3*q4*t3;
  R[3,2] = 0;
  R[3,3] = t4;
  R[3,4] = q3*t2;
  R[3,5] = 0;
  R[4,1] = t3*(t5*-2.0+t7+t9+t10)*0.5;
  R[4,2] = -gami*q2*t2;
  R[4,3] = -gami*q3*t2;
  R[4,4] = -q4*t2*(gami-2.0);
  R[4,5] = gami;
  R[5,1] = 1.0/(q1*q1*q1)*q4*(t7+t9+t10-q1*q5-gami*q1*q5);
  R[5,2] = -gami*q2*q4*t3;
  R[5,3] = -gami*q3*q4*t3;
  R[5,4] = t3*(t7+t9+gami*t5*3.0-q1*q5*2.0-gami*q1*q5*2.0)*-0.5;
  R[5,5] = q4*t2*(gamma);

  return nothing
end


#------------------------------------------------------------------------------
# Eigenvalues
"""
  This function calculates the eigenvalues for the Euler flux jacobian in the
  x direction corresponding to the eigenvectors from calcEvecsx.

  Inputs:
    params: a ParamType
    q: a vector of conservative variables

  Inputs/Outputs:
    Lambda: vector to be populated with the eigenvalues (overwritten)

  Aliasing restrictions: none
"""
function calcEvalsx(params::ParamType{2, :conservative}, q::AbstractVector, Lambda::AbstractVector)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q5 = q[4]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code

  t2 = 1.0/q1;
  t3 = q2*t2;
  t4 = q2*q2;
  t5 = t4*0.5;
  t6 = q3*q3;
  t7 = t6*0.5;
  t8 = t5+t7;
  t9 = q5-t2*t8;
  t10 = gamma;
  t11 = gami*t2*t9*t10;
  t12 = sqrt(t11);

  Lambda[1] = t3;
  Lambda[2] = t3;
  Lambda[3] = t3+t12;
  Lambda[4] = t3-t12;

  return nothing
end

function calcEvalsx(params::ParamType{3, :conservative}, q::AbstractVector, Lambda::AbstractVector)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = 1.0/q1;
  t3 = q2*t2;
  t4 = gamma;
  t5 = q2*q2;
  t6 = t5*0.5;
  t7 = q3*q3;
  t8 = t7*0.5;
  t9 = q4*q4;
  t10 = t9*0.5;
  t11 = t6+t8+t10;
  t12 = q5-t2*t11;
  t13 = gami*t2*t4*t12;
  t14 = sqrt(t13);
  Lambda[1] = t3;
  Lambda[2] = t3;
  Lambda[3] = t3;
  Lambda[4] = t3+t14;
  Lambda[5] = t3-t14;

  return nothing
end
"""
  Like calcEvalsx, but in the y direction.  See that function for details.
"""
function calcEvalsy(params::ParamType{2, :conservative}, q::AbstractVector, Lambda::AbstractVector)

  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q5 = q[4]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = 1.0/q1;
  t3 = q3*t2;
  t4 = q2*q2;
  t5 = t4*0.5;
  t6 = q3*q3;
  t7 = t6*0.5;
  t8 = t5+t7;
  t9 = q5-t2*t8;
  t10 = gamma;
  t11 = gami*t2*t9*t10;
  t12 = sqrt(t11);

  Lambda[1] = t3;
  Lambda[2] = t3;
  Lambda[3] = t3+t12;
  Lambda[4] = t3-t12; 
  
  return nothing
end


function calcEvalsy(params::ParamType{3, :conservative}, q::AbstractVector, Lambda::AbstractVector)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  t2 = 1.0/q1;
  t3 = q3*t2;
  t4 = sqrt(2.0);
  t5 = t2*t2;
  t6 = gamma;
  t7 = q2*q2;
  t8 = q3*q3;
  t9 = q4*q4;
  t10 = t7+t8+t9-q1*q5*2.0;
  t11 = sqrt(-gami*t5*t6*t10);
  t12 = q1*t4*t11*0.5;
  Lambda[1] = t3;
  Lambda[2] = t3;
  Lambda[3] = t3;
  Lambda[4] = t2*(q3+t12);
  Lambda[5] = t2*(q3-t12);

  return nothing
end

"""
  Like calcEvalsx, but in the z direction.  See that function for details
"""
function calcEvalsz(params::ParamType{3, :conservative}, q::AbstractVector, Lambda::AbstractVector)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = 1.0/q1;
  t3 = q4*t2;
  t4 = sqrt(2.0);
  t5 = 1.0/(q1*q1);
  t6 = gamma;
  t7 = q2*q2;
  t8 = q3*q3;
  t9 = q4*q4;
  t10 = t7+t8+t9-q1*q5*2.0;
  t11 = sqrt(-gami*t5*t6*t10);
  t12 = q1*t4*t11*0.5;
  Lambda[1] = t3;
  Lambda[2] = t3;
  Lambda[3] = t3;
  Lambda[4] = t2*(q4+t12);
  Lambda[5] = t2*(q4-t12);

  return nothing
end

#------------------------------------------------------------------------------
# Eigenvectors

"""
  This function calculates the (right) eigenvectors of the Euler flux jacobian
  in the x direction at a given state q.  Methods are available in 2 and 3 dimensions

  Inputs:
    params: a ParamType
    q: a vector of conservative variables

  Inputs/Outputs:
    R: Matrix whose columns will be populated with the eigenvectors (overwritten)

  Aliasing restrictions: none
"""
function calcEvecsx(params::ParamType{2, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q5 = q[4]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = 1.0/q1;
  t3 = sqrt(2.0)*0.5;
  t4 = q2*q2;
  t5 = t4*0.5;
  t6 = q3*q3;
  t7 = t6*0.5;
  t8 = t5+t7;
  t14 = t2*t8;
  t9 = q5-t14;
  t10 = gamma;
  t11 = gami*t2*t9*t10;
  t16 = sqrt(t11);
  t12 = 1.0/t16;
  t13 = q1*t3*t12;
  t15 = q2*t2;
  t17 = q3*t3*t12;
  t18 = t2*t2;
  t19 = 1.0/gami;
  t20 = t4*t18;
  t21 = t6*t18;
  t22 = t20+t21;
  t23 = gami*t22*0.5;
  t24 = t11+t23;
  t25 = t19*t24;
  t26 = q2*t2*t16;
  t27 = q1*t3*t12

  R[1,1] = 1.0;
  R[2,1] = t15;
  R[3,1] = q3*t2;
  R[4,1] = t4*t18*0.5+t6*t18*0.5;

  R[1,2] = 0;
  R[2,2] = 0;
  R[3,2] = -q1;
  R[4,2] = -q3;

  R[1,3] = t13;
  R[2,3] = t27*(t15+t16);
  R[3,3] = t17;
  R[4,3] = t27*(t25+t26);

  R[1,4] = t13;
  R[2,4] = t27*(t15-t16);
  R[3,4] = t17;
  R[4,4] = t27*(t25-t26);

  return nothing
end

function calcEvecsx(params::ParamType{3, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = 1.0/q1;
  t3 = sqrt(2.0);
  t4 = gamma;
  t5 = q2*q2;
  t6 = t5*0.5;
  t7 = q3*q3;
  t8 = t7*0.5;
  t9 = q4*q4;
  t10 = t9*0.5;
  t11 = t6+t8+t10;
  t17 = t2*t11;
  t12 = q5-t17;
  t13 = gami*t2*t4*t12;
  t14 = 1.0/sqrt(t13);
  t15 = q1*t3*t14*0.5;
  t16 = q2*t2;
  t18 = sqrt(t13);
  t19 = q3*t3*t14*0.5;
  t20 = q4*t3*t14*0.5;
  t21 = 1.0/(q1*q1);
  t22 = 1.0/gami;
  t23 = t5*t21;
  t24 = t7*t21;
  t25 = t9*t21;
  t26 = t23+t24+t25;
  t27 = gami*t26*0.5;
  t28 = t13+t27;
  t29 = t22*t28;
  t30 = q2*t2*t18;
  R[1,1] = 1.0;
  R[1,2] = 0;
  R[1,3] = 0
  R[1,4] = t15;
  R[1,5] = t15;
  R[2,1] = t16;
  R[2,2] = 0;
  R[2,3] = 0;
  R[2,4] = q1*t3*t14*(t16+t18)*0.5;
  R[2,5] = q1*t3*t14*(t16-t18)*0.5;
  R[3,1] = q3*t2;
  R[3,2] = 0;
  R[3,3] = -q1;
  R[3,4] = t19;
  R[3,5] = t19;
  R[4,1] = q4*t2;
  R[4,2] = q1;
  R[4,3] = 0;
  R[4,4] = t20;
  R[4,5] = t20;
  R[5,1] = t5*t21*0.5+t7*t21*0.5+t9*t21*0.5;
  R[5,2] = q4;
  R[5,3] = -q3;
  R[5,4] = q1*t3*t14*(t29+t30)*0.5;
  R[5,5] = q1*t3*t14*(t29-t30)*0.5;

  return nothing
end

"""
  Like calcEvecsx, but for the y direction.  See that function for details
"""
function calcEvecsy(params::ParamType{2, :conservative}, q::AbstractVector, R::AbstractMatrix)

  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q5 = q[4]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = 1.0/q1;
  t3 = sqrt(2.0);
  t4 = q2*q2;
  t5 = t4*0.5;
  t6 = q3*q3;
  t7 = t6*0.5;
  t8 = t5+t7;
  t14 = t2*t8;
  t9 = q5-t14;
  t10 = gamma;
  t11 = gami*t2*t9*t10;
  t17 = sqrt(t11);
  t12 = 1.0/t17;
  t12a = t3*t12*0.5;
  t13 = q1*t12a
  t15 = q2*t12a;
  t16 = q3*t2;
  t18 = t2*t2;
  t19 = 1.0/gami;
  t20 = t4*t18;
  t21 = t6*t18;
  t22 = t20+t21;
  t23 = gami*t22*0.5;
  t24 = t11+t23;
  t25 = t19*t24;
  t26 = q3*t2*t17;
  t27 = q1*t3*t12*0.5

  R[1,1] = 1.0;
  R[2,1] = q2*t2;
  R[3,1] = t16;
  R[4,1] = t4*t18*0.5+t6*t18*0.5;

  R[1,2] = 0;
  R[2,2] = q1;
  R[3,2] = 0
  R[4,2] = q2;

  R[1,3] = t13;
  R[2,3] = t15;
  R[3,3] = t27*(t16+t17);
  R[4,3] = t27*(t25+t26);

  R[1,4] = t13;
  R[2,4] = t15;
  R[3,4] = t27*(t16-t17);
  R[4,4] = t27*(t25-t26);

  return nothing
end


function calcEvecsy(params::ParamType{3, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  # auto-generated code
  t2 = 1.0/q1;
  t3 = sqrt(2.0);
  t4 = gamma;
  t5 = q2*q2;
  t6 = t5*0.5;
  t7 = q3*q3;
  t8 = t7*0.5;
  t9 = q4*q4;
  t10 = t9*0.5;
  t11 = t6+t8+t10;
  t16 = t2*t11;
  t12 = q5-t16;
  t13 = gami*t2*t4*t12;
  t14 = 1.0/sqrt(t13);
  t15 = q1*t3*t14*0.5;
  t17 = q2*t3*t14*0.5;
  t18 = q3*t2;
  t19 = sqrt(t13);
  t20 = q4*t3*t14*0.5;
  t21 = 1.0/(q1*q1);
  t22 = 1.0/gami;
  t23 = t5*t21;
  t24 = t7*t21;
  t25 = t9*t21;
  t26 = t23+t24+t25;
  t27 = gami*t26*0.5;
  t28 = t13+t27;
  t29 = t22*t28;
  t30 = q3*t2*t19;
  R[1,1] = 0
  R[1,2] = 1.0;
  R[1,3] = 0;
  R[1,4] = t15;
  R[1,5] = t15;
  R[2,1] = 0;
  R[2,2] = q2*t2;
  R[2,3] = q1;
  R[2,4] = t17;
  R[2,5] = t17;
  R[3,1] = 0;
  R[3,2] = t18;
  R[3,3] = 0;
  R[3,4] = q1*t3*t14*(t18+t19)*0.5;
  R[3,5] = q1*t3*t14*(t18-t19)*0.5;
  R[4,1] = -q1;
  R[4,2] = q4*t2;
  R[4,3] = 0;
  R[4,4] = t20;
  R[4,5] = t20;
  R[5,1] = -q4;
  R[5,2] = t5*t21*0.5+t7*t21*0.5+t9*t21*0.5;
  R[5,3] = q2;
  R[5,4] = q1*t3*t14*(t29+t30)*0.5;
  R[5,5] = q1*t3*t14*(t29-t30)*0.5;

  return nothing
end

"""
  Like calcEvecx, but in the z direction.  See that function for details.
"""
function calcEvecsz(params::ParamType{3, :conservative}, q::AbstractVector, R::AbstractMatrix)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  # auto generated code
  t2 = 1.0/q1;
  t3 = sqrt(2.0);
  t4 = gamma;
  t5 = q2*q2;
  t6 = t5*0.5;
  t7 = q3*q3;
  t8 = t7*0.5;
  t9 = q4*q4;
  t10 = t9*0.5;
  t11 = t6+t8+t10;
  t16 = t2*t11;
  t12 = q5-t16;
  t13 = gami*t2*t4*t12;
  t14 = 1.0/sqrt(t13);
  t15 = q1*t3*t14*0.5;
  t17 = q2*t3*t14*0.5;
  t18 = q3*t3*t14*0.5;
  t19 = q4*t2;
  t20 = sqrt(t13);
  t21 = 1.0/(q1*q1);
  t22 = 1.0/gami;
  t23 = t5*t21;
  t24 = t7*t21;
  t25 = t9*t21;
  t26 = t23+t24+t25;
  t27 = gami*t26*0.5;
  t28 = t13+t27;
  t29 = t22*t28;
  t30 = q4*t2*t20;
  R[1,1] = 0;
  R[1,2] = 0;
  R[1,3] = 1.0;
  R[1,4] = t15;
  R[1,5] = t15;
  R[2,1] = 0;
  R[2,2] = -q1;
  R[2,3] = q2*t2;
  R[2,4] = t17;
  R[2,5] = t17;
  R[3,1] = q1;
  R[3,2] = 0;
  R[3,3] = q3*t2;
  R[3,4] = t18;
  R[3,5] = t18;
  R[4,1] = 0;
  R[4,2] = 0;
  R[4,3] = t19;
  R[4,4] = q1*t3*t14*(t19+t20)*0.5;
  R[4,5] = q1*t3*t14*(t19-t20)*0.5;
  R[5,1] = q3;
  R[5,2] = -q2;
  R[5,3] = t5*t21*0.5+t7*t21*0.5+t9*t21*0.5;
  R[5,4] = q1*t3*t14*(t29+t30)*0.5;
  R[5,5] = q1*t3*t14*(t29-t30)*0.5;

  return nothing
end

#------------------------------------------------------------------------------
# Eigenvector scaling

"""
  This function calculates the diagonal scaling matrix S^2 such that:
    df/dq = (Y*S)*Lambda*inv(Y*S)  and
    dq/dw = (Y*S)*(Y*S).'

  where df/dx is the Euler flux jacobian in the x direction, Y and Lambda are
  its eigenvectors and eigenvalues, respective (which can be calculated via 
  calcEvecsx and calcEvalsx),
  and dq/dw is the jacobian of the conservative variables wrt the IR entropy
  variables.  Recall that eigenvectors are defined up to a multiplicitive 
  constant, so Y*S is also a matrix of eigenvectors.

  This scaling was originally reported by Merriam: An Entropy-Based Approach
  to Nonlinear Stability (NASA TR).

  Methods are available for 2 and 3 dimensions.

  Aliasing restrictions: none
"""
function calcEScalingx(params::ParamType{2, :conservative}, q::AbstractVector, S::AbstractVector)

  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q5 = q[4]
  gami = params.gamma_1
  gamma = params.gamma


  # auto generated code
  t2 = 1.0/(q1*q1*q1);
  t3 = q2*q2;
  t4 = q3*q3;
  t6 = q1*q5*2.0;
  t5 = t3+t4-t6;
  t6 = -gami*t2*t5*0.5;


  S[1] = (gami*q1)/(gamma);
  S[2] = t6
  S[3] = t6
  S[4] = t6;

  return nothing
end

function calcEScalingx(params::ParamType{3, :conservative}, q::AbstractVector, S::AbstractVector)

  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  t2 = 1.0/(q1*q1*q1);
  t3 = q2*q2;
  t4 = q3*q3;
  t5 = q4*q4;
  t7 = q1*q5*2.0;
  t6 = t3+t4+t5-t7;
  S[1] = (gami*q1)/(gamma);
  S[2] = gami*t2*t6*-0.5;
  S[3] = gami*t2*t6*-0.5;
  S[4] = gami*t2*t6*-0.5;
  S[5] = gami*t2*t6*-0.5;

  return nothing
end

"""
  Like calcEScalingx, but in the y direction.  See that function for details.
"""
function calcEScalingy(params::ParamType{2, :conservative}, q::AbstractVector, S::AbstractVector)
  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q5 = q[4]
  gami = params.gamma_1
  gamma = params.gamma


  t2 = 1.0/(q1*q1*q1);
  t3 = q2*q2;
  t4 = q3*q3;
  t6 = q1*q5*2.0;
  t5 = t3+t4-t6;
  t6 = -gami*t2*t5*0.5;

  S[1] = (gami*q1)/(gamma);
  S[2] = t6
  S[3] = t6;
  S[4] = t6;

  return nothing
end

function calcEScalingy(params::ParamType{3, :conservative}, q::AbstractVector, S::AbstractVector)

  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  t2 = 1.0/(q1*q1*q1);
  t3 = q2*q2;
  t4 = q3*q3;
  t5 = q4*q4;
  t7 = q1*q5*2.0;
  t6 = t3+t4+t5-t7;
  S[1] = gami*t2*t6*-0.5;
  S[2] = (gami*q1)/(gamma);
  S[3] = gami*t2*t6*-0.5;
  S[4] = gami*t2*t6*-0.5;
  S[5] = gami*t2*t6*-0.5;

  return nothing
end

"""
  Like calcEScalingx, but in the y direction.  See that function for details.
"""
function calcEScalingz(params::ParamType{3, :conservative}, q::AbstractVector, S::AbstractVector)

  # unpack variables
  q1 = q[1]
  q2 = q[2]
  q3 = q[3]
  q4 = q[4]
  q5 = q[5]
  gami = params.gamma_1
  gamma = params.gamma

  # auto-generated code
  t2 = 1.0/(q1*q1*q1);
  t3 = q2*q2;
  t4 = q3*q3;
  t5 = q4*q4;
  t7 = q1*q5*2.0;
  t6 = t3+t4+t5-t7;
  S[1] = gami*t2*t6*-0.5;
  S[2] = gami*t2*t6*-0.5;
  S[3] = (gami*q1)/(gamma);
  S[4] = gami*t2*t6*-0.5;
  S[5] = gami*t2*t6*-0.5;

  return nothing
end
