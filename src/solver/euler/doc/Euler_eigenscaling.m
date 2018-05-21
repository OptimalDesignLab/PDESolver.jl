% Calculating the eigenvalues, eigenvectors, and scaling matrix
% such that (Y*S)*Lambda*inv(Y*s) = A (the euler flux jacobian
% and (Y*S)*(Y*S)^T = dq/dw
% The scaling idea was shown in 1D by Merriam in "An Entropy-Based Approach
% to Nonlinear Stability"
% The scaling exists in the x, y, and z directions, but not for an
% arbitrary direction (ie if A is the jacobian of the flux in some 
% direction specify by a normal vector [nx ny nz]
% The  R1 and R2 eigenvectors are from: 
% Warming, Beam, Hyett: Diagonalization and Simultaneous Symmetrization of 
% the Gas-Dynamic Matrices
% The R3 eigenvectors are from Pulliams VKI notes and are used to 
% compute the scaling and generate the code (for consistency with the 2D
% case)
% Script written by Jared Crean, December, 2016

clear
clc
close all

syms q1 q2 q3 q4 q5  % conservative variables
%syms rho u v w  % primative variables
syms U p a u v w % functions, but treat them as variables initially
syms nx ny nz  % direction vector
syms gamma gami % constants
syms k1 k2 k3  % normalized components of nx, ny, nz
syms rad2 H
syms phi theta alpha beta  % Pulliams variables

q = [q1; q2; q3; q4; q5]


% k1exp = nx/sqrt(nx*nx + ny*ny + nz*nz)
% k2exp = ny/sqrt(nx*nx + ny*ny + nz*nz)
% k2exp = nz/sqrt(nx*nx + ny*ny + nz*nz)
thetaexpr = k1*u + k2*v + k3*w;
%a1expr = gamma*q5/q1 - phi*phi;
phiexpr = sqrt(0.5*gami*(u*u + v*v + w*w));
alphaexpr = q1/(rad2*a)
betaexpr = 1/(rad2*q1*a)
k1exp = nx
k2exp = ny
k3exp = nz
Hexpr = q5/q1 + p/q1
Uexpr = (q2*nx + q3*ny + q4*nz)/q1
aexpr = sqrt(gamma*p/q1)
pexpr = gami*(q5 - 0.5*(q2*q2 + q3*q3 + q4*q4)/q1)
uexpr = q2/q1;
vexpr = q3/q1;
wexpr = q4/q1;

% arrays used for substitutions
sym_vars = [theta, phi, alpha, beta, k1, k2, k3, H, U, a, p, gamma, rad2, u, v, w]
sym_exprs = [thetaexpr, phiexpr, alphaexpr, betaexpr, k1exp, k2exp, k3exp, Hexpr, Uexpr, aexpr, pexpr, gami+1, sqrt(2), uexpr, vexpr, wexpr]

% From Warming and Beam paper
T = [k1,  k2,   k3, q1/(rad2*a),   q1/(rad2*a);
      0, -k3 ,  k2, k1/rad2    ,   -k1/rad2;
     k3,   0 , -k1, k2/rad2,       -k2/rad2;
     -k2, k1,   0,  k3/rad2,       -k3/rad2;
     0,    0,   0,  q1*a/rad2,     q1*a/rad2]

 M = [1 0 0 0 0; 
     q2/q1 q1 0 0 0;
     q3/q1 0 q1 0 0;
     q4/q1 0 0 q1 0;
     0.5*(q2*q2 + q3*q3 + q4*q4)/(q1*q1) q2 q3 q4 1/gami]
 
 % conservative variable flux jacobian right eigenvectors
 R = M*T;
 
 % eigenvalues
 Lambda = [U, U, U, U + a, U - a];
 
 % diff is zero, showing it is ok to remove the factors of rad2
 A = simplify(R*diag(Lambda)*inv(R));
 T2 = T;
 T2 = subs(T2, rad2, 1);
 R2 = M*T2;
 A2 = R2*diag(Lambda)*inv(R2);
 A2 = simplify(A2);
 diff = simplify(A - A2)
 
 
 % eigenvectors from Pulliams VKI Notes
 R3 = [k1, k2, k3, alpha, alpha;
       k1*u, k2*u - k3*q1, k3*u + k2*q1,  alpha*(u + k1*a), alpha*(u - k1*a);
       k1*v + k3*q1, k2*v, k3*v - k1*q1, alpha*(v + k2*a), alpha*(v - k2*a);
       k1*w - k2*q1,  k2*w + k1*q1, k3*w, alpha*(w + k3*a), alpha*(w - k3*a);
       k1*phi*phi/gami + q1*(k3*v - k2*w),  k2*phi*phi/gami + q1*(k1*w - k3*u), k3*phi*phi/gami + q1*(k2*u - k1*v), alpha*( (phi*phi + a*a)/gami + theta*a), alpha*( (phi*phi + a*a)/gami - theta*a)];
  A3 = R3*diag(Lambda)*inv(R3)
      
 R3inv = [k1*(1 - phi*phi/(a*a)) - (k3*v - k2*w)/q1, k1*gami*u/(a*a), k1*gami*v/(a*a) + k3/q1, k1*gami*w/(a*a) - k2/q1, -k1*gami/(a*a);
          k2*(1 - phi*phi/(a*a)) - (k1*w - k3*u)/q1, k2*gami*u/(a*a) - k3/q1, k2*gami*v/(a*a), k2*gami*w/(a*a) + k1/q1, -k2*gami/(a*a);
          k3*(1 - phi*phi/(a*a)) - (k2*u - k1*v)/q1, k3*gami*u/(a*a) + k2/q1, k3*gami*v/(a*a) - k1/q1, k3*gami*w/(a*a), -k3*gami/(a*a);
          beta*(phi*phi - theta*a), -beta*gami*u - k1*a, -beta*(gami*v - k2*a),  -beta*(gami*w - k3*a), beta*gami;
          beta*(phi*phi + theta*a), -beta*gami*u + k1*a, -beta*(gami*v + k2*a), -beta*(gami*w + k3*a), beta*gami];
 % use R3 instead of R
 R = R3
 
 A0 = [q1, q2, q3, q4, q5;
       q2,  q2*q2/q1  + p, q2*q3/q1, q2*q4/q1, H*q2;
       q3, q2*q3/q1,  q3*q3/q1 + p, q4*q3/q1, H*q3;
       q4, q2*q4/q1, q3*q4/q1, q4*q4/q1 + p, H*q4;
       q5, H*q2, H*q3, H*q4, q1*H*H - a*a*p/gami]
   
 % convert everything to conservative variables
 for i=1:length(sym_vars)
     sym_i = sym_vars(i);
     sym_expr = sym_exprs(i);
     R = subs(R, sym_i, sym_expr);
     R3inv = subs(R3inv, sym_i, sym_expr);

     Lambda = subs(Lambda, sym_i, sym_expr);
     A0 = subs(A0, sym_i, sym_expr);
     A = subs(A, sym_i, sym_expr);
     A2 = subs(A2, sym_i, sym_expr);
     A3 = subs(A3, sym_i, sym_expr);
 end
 
 % verify the different eigenvectors all produce the same flux jacobian
 A = simplify(A);
 A2 = simplify(A2);
 A3 = simplify(A3);
 diff3 = simplify(A - A3)
 
 
 % pick a direction - x
 sym_vars4 = [k1 k2 k3 nx ny nz];
 sym_vals4 = [nx ny nz 1 0 0];
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
 Rx = simplify(Rx);
 Lambdax = subs(Lambdax, gamma, gami+1);
 S2x = inv(Rx)*A0*inv(Rx).';
 S2x = simplify(S2x)
 
 ccode(Rx, 'file', 'Rx.c')
 ccode(diag(S2x), 'file', 'S2x.c')
 ccode(Lambdax, 'file', 'Lx.c')
 Ax2 = Ax;  % matlab complains if the name of the sym passed to ccode is Ax
 ccode(Ax2, 'file', 'A.c')
 
 
 % pick a direction - y
 sym_vals4 = [nx ny nz 0 1 0];
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
 Ry = simplify(Ry);
 Lambday = simplify(Lambday);
 S2y = inv(Ry)*A0*inv(Ry).';
 S2y = simplify(S2y)
 Ay = simplify(Ay);
 
 ccode(Ry, 'file', 'Ry.c')
 ccode(diag(S2y), 'file', 'S2y.c')
 ccode(Lambday, 'file', 'Ly.c')
 ccode(Ay, 'file', 'B.c')


 
 % pick a direction - z
 sym_vals4 = [nx ny nz 0 0 1]
 Rz = R;
 Lambdaz = Lambda;
 Az = A;
 for i=1:length(sym_vars4)
     sym_i = sym_vars4(i);
     sym_expr = sym_vals4(i);
     Rz = subs(Rz, sym_i, sym_expr);
     Lambdaz = subs(Lambdaz, sym_i, sym_expr);
     Az = subs(Az, sym_i, sym_expr);
 end
 Rz = simplify(Rz);
 Lambdaz = simplify(Lambdaz);
 S2z = inv(Rz)*A0*inv(Rz).';
 S2z = simplify(S2z)
 Az = simplify(Az);

 
 ccode(Rz, 'file', 'Rz.c')
 ccode(diag(S2z), 'file', 'S2z.c') 
 ccode(Lambdaz, 'file', 'Lz.c')
 ccode(Az, 'file', 'C.c')

 
 