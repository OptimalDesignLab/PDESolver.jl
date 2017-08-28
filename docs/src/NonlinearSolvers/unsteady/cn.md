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

The unsteady adjoint derivation starts with the generic Lagrangian equation:

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

And finally, the global-in-time objective function vector is:

$\boldsymbol{J}^T = [J_1, J_2, \dots, J_i, J_{i+1}, \dots, J_n]$

Therefore, the full discrete Lagrangian is:

$\boldsymbol{\mathcal{L}}(\boldsymbol{u}, \boldsymbol{\psi}) = \boldsymbol{\psi}^T \boldsymbol{R(\boldsymbol{u})} + \boldsymbol{J}(\boldsymbol{u}) = \begin{bmatrix} \psi_1 \left( u_1 - u_0 - \frac{\Delta t}{2} R(u_1) - \frac{\Delta t}{2} R(u_0) \right) \\ \psi_2 \left( u_2 - u_1 - \frac{\Delta t}{2} R(u_2) - \frac{\Delta t}{2} R(u_1) \right) \\ \vdots \\ \psi_i \left( u_i - u_{i-1} - \frac{\Delta t}{2} R(u_i) - \frac{\Delta t}{2} R(u_{i-1}) \right) \\ \psi_{i+1} \left( u_{i+1} - u_{i} - \frac{\Delta t}{2} R(u_{i+1}) - \frac{\Delta t}{2} R(u_{i}) \right) \\ \vdots \\ \psi_n \left( u_n - u_{n-1} - \frac{\Delta t}{2} R(u_n) - \frac{\Delta t}{2} R(u_{n-1}) \right) \end{bmatrix} + \begin{bmatrix} J(u_1) \\ J(u_2) \\ \vdots \\ J(u_i) \\ J(u_{i+1}) \\ \vdots \\ J(u_n) \end{bmatrix}$

Taking the derivative of the Lagrangian with respect to the state at step $i$ yields:

$\frac{\partial \boldsymbol{L}}{\partial u_i} = \underbrace{\psi_i^T - \psi_i^T \frac{\Delta t}{2} \frac{\partial R(u_i)}{\partial u_i}}_{\text{contribution from }\boldsymbol{R}(u_i)} - \underbrace{\psi_{i+1}^T - \psi_{i+1}^T \frac{\Delta t}{2} \frac{\partial R(u_i)}{\partial u}}_{\text{contribution from }\boldsymbol{R}(u_{i+1})} + \frac{\partial J(u_i)}{\partial u_i}= 0^T$

Or, rearranging:

$\frac{\partial \boldsymbol{L}}{\partial u_i} = (\psi_i - \psi_{i+1}) - \frac{\Delta t}{2} \left( \frac{\partial R(u_i)}{\partial u_i} \right)^T (\psi_i + \psi_{i+1}) + \frac{\partial J(u_i)}{\partial u_i} = 0$

### Initial Condition

The derivative of the Lagrangian with respect to the state at the final step $i = n$ is:

$\frac{\partial \boldsymbol{L}}{\partial u_n} = \psi_n - \frac{\Delta t}{2} \left( \frac{\partial R(u_n)}{\partial u_n} \right)^T \psi_n + \frac{\partial J(u_n)}{\partial u_n} = 0$

Therefore, the value of the adjoint at time step n, which is the initial condition for the reverse sweep, is:

$\psi_n = \left( \left(I - \frac{\Delta t}{2} \frac{\partial R(u_n)}{\partial u_n} \right)^T \right)^{-1} \left( - \frac{\partial J(u_n)}{\partial u_n} \right)^T$


### Direct Solve

The method of performing a direct solve to advance the CN reverse sweep (as opposed to using Newton's method to converge each time step) starts with the restatement of the derivative of the Lagrangian at time step $i$:

$\frac{\partial \boldsymbol{L}}{\partial u_i} = \underbrace{\psi_i^T - \psi_i^T \frac{\Delta t}{2} \frac{\partial R(u_i)}{\partial u_i}}_{\text{contribution from }\boldsymbol{R}(u_i)} - \underbrace{\psi_{i+1}^T - \psi_{i+1}^T \frac{\Delta t}{2} \frac{\partial R(u_i)}{\partial u}}_{\text{contribution from }\boldsymbol{R}(u_{i+1})} + \frac{\partial J(u_i)}{\partial u_i}= 0^T$

Rearranging:

$\left[ \psi_i - \frac{\Delta t}{2} \left( \frac{\partial R(u_i)}{\partial u_i} \right)^T \psi_i \right] - \left[ \psi_{i+1} + \frac{\Delta t}{2} \left( \frac{\partial R(u_i)}{\partial u_i} \right)^T \psi_{i+1} \right] + \frac{\partial J(u_i)}{\partial u_i} = 0$

Grouping terms to isolate $\psi_i$:

$\left[ I - \frac{\Delta t}{2} \left( \frac{\partial R(u_i)}{\partial u_i} \right)^T \right] \psi_i = \left[ \psi_{i+1} + \frac{\Delta t}{2} \left( \frac{\partial R(u_i)}{\partial u_i} \right)^T \psi_{i+1} \right] - \frac{\partial J(u_i)}{\partial u_i}$

Solving for $\psi_i$:

$\psi_i = \left[ I - \frac{\Delta t}{2} \left( \frac{\partial R(u_i)}{\partial u_i} \right)^T \right]^{-1} \left( \left[ \psi_{i+1} + \frac{\Delta t}{2} \left( \frac{\partial R(u_i)}{\partial u_i} \right)^T \psi_{i+1} \right] - \frac{\partial J(u_i)}{\partial u_i} \right)$


### Checkpointing

Currently, all time steps are checkpointed. 
Eventually, Revolve will be implemented, for which a separate Julia package has been developed. 
See [here](http://dl.acm.org/citation.cfm?id=347846) for the publication discussing the Revolve algorithm.


### Global-in-time Jacobian

For reference, the structure of the global-in-time Jacobian is shown here.
It should never be formed except in the course of debugging very simple use cases, 
  but it can be helpful for visualizing the matrix form of CN for all space and time.




