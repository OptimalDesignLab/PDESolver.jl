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
\underbrace{\mathbf{I} - \frac{\Delta t}{2} \underbrace{\frac{\partial R\left(u^{i+1}\right)}{\partial u^{i+1}}}_{\text{jac from \texttt{physicsJac}}} \
}_{\text{\texttt{cnJac}}} \
= \\
\underbrace{- \left( \
  \underbrace{u^{i+1}}_{\text{\texttt{eqn\_nextstep.q\_vec}}} - \ 
    \frac{\Delta t}{2} \underbrace{R\left(u^{i+1}\right)}_{\text{\texttt{eqn\_nextstep.res\_vec}}} - \ 
  \underbrace{u^i}_{\text{\texttt{eqn.q\_vec}}} - \ 
    \frac{\Delta t}{2} \underbrace{R\left(u^i\right)}_{\text{\texttt{eqn.res\_vec}}} \
\right)}_{\text{\texttt{cnRhs}}}
\end{aligned}

This equation is solved with PDESolver's Newton's method.

## Backward-in-time (unsteady adjoint)
