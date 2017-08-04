# Crank-Nicolson (CN)

Crank-Nicolson is an implicit time marching scheme.
While it adds a measure of complexity, it is unconditionally stable. 
Note that this does not mean that accuracy is unaffected by time step size, 
degree of the operator (p = 1 or p = 2, etc.), or element size.

## Selecting CN & Problem Setup

CN is selected with `"run_type" => 20` in your problem's input options file.
As this is an unsteady problem, you will also need to specify your 
  time step size (`"delta_t"` option) and your final solution time (`"t_max"` option)


## Forward-in-time

\begin{aligned}
\underbrace{\mathbf{I} - \frac{\Delta t}{2} \underbrace{\frac{\partial R\left(u^{i+1}\right)}{\partial u^{i+1}}}_{\text{jac from physicsJac}} \
}_{\text{cnJac}} \
= \\
\underbrace{- \left( \
  \underbrace{u^{i+1}}_{\text{eqn\_nextstep.q\_vec}} - \ 
    \frac{\Delta t}{2} \underbrace{R\left(u^{i+1}\right)}_{\text{eqn\_nextstep.res\_vec}} - \ 
  \underbrace{u^i}_{\text{eqn.q\_vec}} - \ 
    \frac{\Delta t}{2} \underbrace{R\left(u^i\right)}_{\text{eqn.res\_vec}} \
\right)}_{\text{cnRhs}}
\end{aligned}

This equation is solved with PDESolver's Newton's method.

## Backward-in-time (unsteady adjoint)

The unstead adjoint derivation starts with the generic Lagrangian equation:

\begin{equation}
\mathcal{L}(u, \psi) = \psi^T R(u) + J(u)
\end{equation}

In the discrete context of CN, all of these variables are global-in-time.
That is, the adjoint vector contains the adjoint at time step 1 concatenated with 
  the adjoint at time step 2, and so on, until time step $n$.
Therefore, in this document we will rewrite the Lagrangian using bolded symbols to indicate 
  that a vector or matrix is global-in-time, as there will also be corresponding variables
  specific to a particular time step:

\begin{equation}
\boldsymbol{\mathcal{L}}(\boldsymbol{u}, \boldsymbol{\psi}) = \boldsymbol{\psi}^T \boldsymbol{R}(\boldsymbol{u}) + \boldsymbol{J}(\boldsymbol{u})
\end{equation}

The global-in-time residual discretized according to the Crank-Nicolson method is:

$\boldsymbol{R(\boldsymbol{u})} = \begin{bmatrix} u_1 - u_0 - \frac{\Delta t}{2} R(u_1) - \frac{\Delta t}{2} R(u_0) \\ u_2 - u_1 - \frac{\Delta t}{2} R(u_2) - \frac{\Delta t}{2} R(u_1) \\ \vdots \\ u_i - u_{i-1} - \frac{\Delta t}{2} R(u_i) - \frac{\Delta t}{2} R(u_{i-1}) \\ u_{i+1} - u_{i} - \frac{\Delta t}{2} R(u_{i+1}) - \frac{\Delta t}{2} R(u_{i}) \\ \vdots \\ u_n - u_{n-1} - \frac{\Delta t}{2} R(u_n) - \frac{\Delta t}{2} R(u_{n-1}) \end{bmatrix}$

The global-in-time adjoint vector is:

$\boldsymbol{\psi}^T = [\psi_1, \psi_2, \dots, \psi_i, \psi_{i+1}, \dots, \psi_n]$

Taking the derivative of the Lagrangian with respect to the state at step $i$ yields:


### Initial Condition


### Direct Solve

### Checkpointing







