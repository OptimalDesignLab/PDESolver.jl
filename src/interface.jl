# functions that each physics module must implement
"""
  This function evalutes dq/dt = R(q).  For steady problems it evalutes R(q)
  at some state q.  The state is stored in eqn.q, and eqn.res is populated with
  R(q).  Note that these are 3 dimensional arrays.  The physics modules only
  interact with the 3 dimensional arrays, never the vectors eqn.q_vec and
  eqn.res_vec.  Each physics module must implement this function for its
  own subtype of AbstractSolutionData (ie. with a more specific type for
  the eqn argument and equallty specific types for the other arguments).
  This is important because evalResidual is common to all physics modules,
  so a user or some other part of the code can call evalResidual(mesh, sbp
  eqn, opts), and Julia's multiple dispatch will figure out the right method
  to call based on the type of the eqn argument.

  The evaluation of the residual R(q) should depend only on the data stored in
  mesh, sbp, eqn, and opts, and any data that depend on q should be recalculated
  every time the function is called.  This function is used as a building block
  by other parts of the solver, particularly the NonlinearSolvers.  See
  interfaces.md for details

  **Inputs**

   * mesh: an AbstractMesh describing the mesh on which to solve the physics
   * sbp: an SBP operator
   * eqn: a subtype of AbstractSolution data, used to store all of the data used by the physics module
   * opts: the options dictionary
   * t: the current time value, defaults to 0.0

  #TODO: list required options keys
"""
function evalResidual(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict, t = 0.0)

  throw(ErrorException("Generic fallback evalResidual reached: did you forget to extend evalResidual() with a new method for your AbstractSolutionData?"))

  return nothing
end

function evalHomotopy(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict, res::Abstract3DArray, t = 0.0)

  throw(ErrorException("Generic fallback evalHomotopy reached: did you forget to extend evalHomotopy() with a new method for your AbstractSolutionData?"))

  return nothing
end

"""
  Similar function to [`evalResidual`](@ref), but instead of computing the
  residual, it computes and assembles the Jacobian of the residual with
  respect to `eqn.q` into the system matrix.

  The functions [`assembleElement`](@ref) and [`assembleInterface`](@ref)
  in `NonlinearSolvers` should be used to assembles the Jacobians of
  individual elements and interfaces into the system matrix.

  Physics modules that can compute element Jacobians directly
  should extend this function with a new method, otherwise the coloring
  algorithm with be used with [`evalResidual`](@ref) to compute the
  system matrix.  For physics modules that support multiple formuations,
  this function should throw an exception if called with an unsupported
  formulation.

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * eqn: AbstractSolutionData (physics modules should specialize this
          argument)
   * opts: options dictionary
   * assembler: object that must be passed to `assembleElement` and 
                `assembleInterface`

  **Keyword Arguments**

   * start_comm: whether or not to start parallel communication,
                 default false because the functions Nonlinear solvers
                 have already done parallel communication when this
                 function gets called.

  **Options Keys**

  TODO: describe the key that controls whether this function gets called
        or not.
"""
function evalJacobian(mesh::AbstractMesh, sbp::AbstractSBP,
                      eqn::AbstractSolutionData, opts::Dict, 
                      assembler::AssembleElementData, t=0.0;
                      start_comm=false)


  throw(ErrorException("Generic fallback evalJacobian reached: did you forget to extend evalJacobian() with a new method for your AbstractSolutionData?"))
  return nothing
end

"""
  Evaluates the Jacobian of the [`evalHomotopy`](@ref) multiplied by the
  homotopy parameter lambda.  Conceptually similar to [`evalJacobian`](@ref).

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an AbstractSBP
   * eqn: an AbstractSolutionData (physics modules should specialize this
          argument)
   * opts: options dictionary
   * assembler: an object used to assemble contributions into the Jacobian
   * lambda: the homotopy parameter
"""
function evalHomotopyJacobian(mesh::AbstractMesh, sbp::AbstractSBP,
                              eqn::AbstractSolutionData, opts::Dict, 
                              assembler::AssembleElementData, lambda::Number)


  throw(ErrorException("Generic fallback evalHomotopyJacobian reached: did you forget to extend evalHomotopyJacobian() with a new method for your AbstractSolutionData?"))

  return nothing
end

#TODO: debugging
function evalHomotopy(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict, t = 0.0)

  evalHomotopy(mesh, sbp, eqn, opts, eqn.res, t)

  return nothing
end

"""
  This function evaluates the mesh.numDofPerNode x mesh.numDofPerNode block
  diagonal of the Jacobian of the strong form of the sptail residual.  Note that
  the residual is written

  du/dt = -(Q * f) + SAT

  Note the negative sign.

  Currently this function neglects the SAT terms (both interface and boundary
  conditions)

  This function is most useful with [`AssembleDiagJacData`](@ref), but
  can be used with any [`AssembleElementData`](@ref).

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * eqn: AbstractSolutionData (physics modules should specialize this
          argument)
   * opts: options dictionary
   * assembler: object that must be passed to `assembleElement` and 
                `assembleInterface`
"""
function evalJacobianStrongDiag(mesh::AbstractMesh, sbp::AbstractSBP,
                      eqn::AbstractSolutionData, opts::Dict, 
                      assembler::AssembleElementData, t=0.0;
                      start_comm=false)


  throw(ErrorException("Generic fallback evalJacobianStrongDiag reached: did you forget to extend evalJacobian() with a new method for your AbstractSolutionData?"))

end


"""
  High level function that evaluates the given functional
  This function is agnostic to the type of the functional being
  computed and calls a mid level functional-type specific function for the 
  actual evaluation.

  The functional is evaluated at the state in eqn.q_vec.

  **Arguments**

  *  `mesh` :  Abstract mesh object
  *  `sbp`  : Summation-By-Parts operator
  *  `eqn`  : AbstractSolutionData object
  *  `opts` : Options dictionary
  *  `functionalData` : Object of type AbstractFunctional. 
                        This type determines the functional being computed and
                        holds all the relevant data.
"""
function evalFunctional{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},
                        sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts,
                        functionalData::AbstractFunctional)

  error("Generic fallback for evalFunctional() reached: did you forget to extend evalFunctional() with a new method for you AbstractSolutionData")

end


"""
  Computes a 3D array of hte derivative of a functional wrt eqn.q.

  The derivative is evaluated at the state in eqn.q_vec.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functionalData: AbstractIntegralFunctional to evaluate

  **Inputs/Outputs**

   * func_deriv_arr: array to hold derivative of function wrt eqn.q, same
                     size as equation.q

  **Options Keys**

  This funciton is not compatible with `precompute_q_bndry` = false
"""
function evalFunctionalDeriv{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractSBP,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           functionalData::AbstractIntegralFunctional,
                           func_deriv_arr::Abstract3DArray)

  error("Generic fallback for evalFunctionalDeriv() reached: did you forget to extend evalFunctionalDeriv() with a new method for you AbstractSolutionData")

end
