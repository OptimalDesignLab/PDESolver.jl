# Unsteady NonlinearSolver Documentation

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
