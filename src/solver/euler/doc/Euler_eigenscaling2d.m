% Calculating the eigenvalues, eigenvectors, and scaling matrix
% such that (Y*S)*Lambda*inv(Y*s) = A (the euler flux jacobian
% and (Y*S)*(Y*S)^T = dq/dw
% The scaling idea was shown in 1D by Merriam in "An Entropy-Based Approach
% to Nonlinear Stability"
% The scaling exists in the x, y, and z directions, but not for an
% arbitrary direction (ie if A is the jacobian of the flux in some 
% direction specify by a normal vector [nx ny nz]
% The eigenvectors are from Pulliams VKI notes and are used to 
% compute the scaling and generate the code (for consistency with the 2D
% case)
% Script written by Jared Crean, December, 2016

clear
clc
close all

syms q1 q2 q3 q4 q5  % conservative variables
syms u v w  % primative variables
syms U p a % functions, but treat them as variables initially
syms nx ny nz  % direction vector
syms gamma gami % constants
syms k1 k2 k3  % normalized components of nx, ny, nz
syms rad2 H
syms phi theta a1 alpha  % Pulliams variables

q = [q1; q2; q3; q4; q5]

% flux jacobian
A = [0 k1 k2 0;
     -u*theta + k1*phi*phi, 0 + theta - (gamma-2)*k1*u, k2*u - (gami)*k1*v, gami*k1;
     -v*theta + k2*phi*phi, k1*v - gami*k2*u, 0 + theta - (gamma-2)*k2*v, gami*k2;
     theta*(phi*phi - a1), k1*a1 - gami*u*theta, k2*a1 - gami*v*theta, gamma*theta + 0]
 
% k1exp = nx/sqrt(nx*nx + ny*ny + nz*nz)
% k2exp = ny/sqrt(nx*nx + ny*ny + nz*nz)
% k2exp = nz/sqrt(nx*nx + ny*ny + nz*nz)
thetaexpr = k1*u + k2*v;
a1expr = gamma*q5/q1 - phi*phi;
phiexpr = sqrt(0.5*gami*(u*u + v*v));
alphaexpr = q1/(rad2*a)
k1exp = nx;
k2exp = ny;
k3exp = nz;
Hexpr = q5/q1 + p/q1;
Uexpr = (q2*nx + q3*ny + q5*nz)/q1;
aexpr = sqrt(gamma*p/q1);
pexpr = gami*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1);
uexpr = q2/q1;
vexpr = q3/q1;


% arrays used for substitutions
sym_vars = [theta, a1, phi, alpha, u, v, H, U, a, p, gamma, rad2, nx, ny, nz, q4]
sym_exprs = [thetaexpr, a1expr, phiexpr, alphaexpr, uexpr, vexpr, Hexpr Uexpr, aexpr, pexpr, gami+1, sqrt(2), k1, k2, 0, 0]


 for i=1:length(sym_vars)
     sym_i = sym_vars(i);
     sym_expr = sym_exprs(i);
     A = subs(A, sym_i, sym_expr);
 end

 A = simplify(A)

% Note: k1 and k2 must be components of a *unit* normal vector
R = [1 0 alpha alpha;
     u, k2*q1, alpha*(u + k1*a), alpha*(u - k1*a);
     v, -k1*q1, alpha*(v + k2*a), alpha*(v - k2*a);
     phi*phi/gami, q1*(k2*u - k1*v), alpha*( (phi*phi + a*a)/gami + a*theta), alpha*( (phi*phi + a*a)/gami - a*theta)]
     
Lambda = [U, U, U + a, U - a]

% convert everything to conservative variables
 for i=1:length(sym_vars)
     sym_i = sym_vars(i);
     sym_expr = sym_exprs(i);
     R = subs(R, sym_i, sym_expr);
     Lambda = subs(Lambda, sym_i, sym_expr);
 end
R = simplify(R);

disp('before subs')
A2 = R*diag(Lambda)*inv(R);
A2 = simplify(A2);
A2 = subs(A2, k1*k1 + k2*k2, 1);
A2 = simplify(A2);
A;
diff = simplify(A - A2);
diff = subs(diff, k1*k1 + k2*k2, 1);
diff = simplify(diff)

 
 A0 = [q1, q2, q3, q4, q5;
       q2,  q2*q2/q1  + p, q2*q3/q1, q2*q4/q1, H*q2;
       q3, q2*q3/q1,  q3*q3/q1 + p, q4*q3/q1, H*q3;
       q4, q2*q4/q1, q3*q4/q1, q4*q4/q1 + p, H*q4;
       q5, H*q2, H*q3, H*q4, q1*H*H - a*a*p/gami];
   
 idx = [1 2 3 5];
 A0 = A0(idx, idx)
 for i=1:length(sym_vars)
     sym_i = sym_vars(i);
     sym_expr = sym_exprs(i);
     A0 = subs(A0, sym_i, sym_expr);
 end
 
 A0 = simplify(A0);
 
  % x direction
  sym_vars4 = [k1 k2 k3];
  sym_vals4 = [1 0 0];
  Rx = R;
  Lambdax = Lambda;
  Ax = A;
  for i=1:length(sym_vars4)
     sym_i = sym_vars4(i);
     sym_expr = sym_vals4(i);
     Rx = subs(Rx, sym_i, sym_expr);
     Lambdax = subs(Lambdax, sym_i, sym_expr);
     Ax = subs(Ax, sym_i, sym_expr);
  end
  Lambdax = subs(Lambdax, gamma, gami+1);
  Rx = simplify(Rx);
  
  % verify flux jacobian
  Ax2 = simplify(Rx*diag(Lambdax)*inv(Rx));
  diffx = simplify(Ax2 - Ax)

  S2x = inv(Rx)*A0*inv(Rx).';
  S2x = simplify(S2x)
  
  ccode(Rx, 'file', 'Rx2d.c')
  ccode(diag(S2x), 'file', 'S2x2d.c')
  ccode(Lambdax, 'file', 'Lx2d.c')
  ccode(Ax2, 'file', 'A2d.c')


  % pick a direction - y
  sym_vals4 = [0 1 0];
  Ry = R;
  Lambday = Lambda;
  Ay = A;
  for i=1:length(sym_vars4)
     sym_i = sym_vars4(i);
     sym_expr = sym_vals4(i);
     Ry = subs(Ry, sym_i, sym_expr);
     Lambday = subs(Lambday, sym_i, sym_expr);
     Ay = subs(Ay, sym_i, sym_expr);
  end
  Lambday = subs(Lambday, gamma, gami+1);
  Rx = simplify(Rx);
  S2y = inv(Ry)*A0*inv(Ry).';
  S2y = simplify(S2y)

  Ay2 = simplify(Ry*diag(Lambday)*inv(Ry));
  diffy = simplify(Ay2 - Ay)

  ccode(Ry, 'file', 'Ry2d.c')
  ccode(diag(S2y), 'file', 'S2y2d.c') 
  ccode(Lambday, 'file', 'Ly2d.c')
  ccode(Ay, 'file', 'B2d.c')