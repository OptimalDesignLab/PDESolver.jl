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


function evalHomotopy(mesh::AbstractMesh, sbp::AbstractSBP,
                      eqn::AbstractSolutionData, opts::Dict,
                      res::Abstract3DArray, t = 0.0)

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
  This function evaluates the Jacobian of the strong form of the spatial
  residual.  Note that the residual is written

  du/dt = -(Q * f) + SAT

  Note the negative sign.

  Currently this function neglects the SAT terms (both interface and boundary
  conditions)

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * eqn: AbstractSolutionData (physics modules should specialize this
          argument)
   * opts: options dictionary
   * assembler: object that must be passed to `assembleElement` and 
                `assembleInterface`
"""
function evalJacobianStrong(mesh::AbstractMesh, sbp::AbstractSBP,
                      eqn::AbstractSolutionData, opts::Dict, 
                      assembler::AssembleElementData, t=0.0;
                      start_comm=false)


  throw(ErrorException("Generic fallback evalJacobianStrongDiag reached: did you forget to extend evalJacobian() with a new method for your AbstractSolutionData?"))

end


"""
  Creates a functional object.

  Each physics modules should extend this function with a new method,
  specializing the `eqn` argument.

  **Inputs**

  * `mesh` : Abstract PUMI mesh
  * `sbp`  : Summation-by-parts operator
  * `eqn`  : AbstractSolutionData object
  * `opts` : Options dictionary
  * `functional_name`: the name of the functional
  * `functional_bcs`: the boundary condition numbers the functional is
                      computed on.

  **Outputs**

   * functional: an [`AbstractFunctional`](@ref).  Usually an
                 [`AbstractIntegralFunctional`](@ref), although this is not
                 required.
"""
function createFunctional(mesh::AbstractMesh, sbp::AbstractSBP,
                        eqn::AbstractSolutionData, opts,
                        functional_name::AbstractString,
                        functional_bcs::Vector{I}) where I<:Integer



  error("generic fallback for createFunctional() reached.  Did you forget to extend createFunctional() with a new method for your AbstractSolutionData?")

end


"""
  High level function that evaluates the given functional
  This function is agnostic to the type of the functional being
  computed and calls a mid level functional-type specific function for the 
  actual evaluation.

  This is the method physics modules must implemented.  Users should call
  the other method, which takes care of starting parallel communication.
  On entry, `eqn.q` will contain the solution.  If the
  functional requires parallel communication of the solution, it will have
  been started before this function is called, and this function is responsible
  for finishing parallel communication by calling [`finishExchangeData`](@ref).

  This function should perform any global reduction required (ex. calling
  `MPI.Allreduce`).


  **Arguments**

  *  `mesh` :  Abstract mesh object
  *  `sbp`  : Summation-By-Parts operator
  *  `eqn`  : AbstractSolutionData object
  *  `opts` : Options dictionary
  *  `functionalData` : Object of type AbstractFunctional. 
                        This type determines the functional being computed and
                        holds all the relevant data.
"""
function _evalFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts,
            functionalData::AbstractFunctional) where {Tmsh, Tsol}

  error("Generic fallback for _evalFunctional() reached: did you forget to extend evalFunctional() with a new method for you AbstractSolutionData")

end


"""
  Performs reverse-mode differentiation of a functional with respect to the
  metrics.  When this function is called, the solution will be in `eqn.q` and
  parallel communication will already have started if required by the functional.
  This function should update the `_bar` fields of the mesh.  These fields
  should never be overwritten.
  It is the callers responsiblity to zero out these fields beforehand, if
  required.  Mesh implementation should provide a function `zeroBarArrays`
  to do this.  See the documentation of [`_evalFunctional`](@ref) for
  information on parallel communication.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functionalData: the [`AbstractIntegralFunctional`](@ref)
"""
function _evalFunctionalDeriv_m(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractSBP,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           functionalData::AbstractFunctional
                           ) where {Tmsh, Tsol}


  error("Generic fallback for evalFunctionalDeriv_m() reached: did you forget to extend evalFunctionalDeriv_m() with a new method for you AbstractSolutionData")

end


"""
  Computes a 3D array of the derivative of a functional wrt eqn.q.

  The derivative is evaluated at the state in eqn.q.  See the documentation
  for [`_evalFunctional`](@ref) for information on parallel communication.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functionalData: AbstractIntegralFunctional to evaluate

  **Inputs/Outputs**

   * func_deriv_arr: array to hold derivative of function wrt eqn.q, same
                     size as equation.q.  To be overwritten.

  **Options Keys**

  This function is not compatible with `precompute_q_bndry` = false in most
  cases
"""
function _evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractSBP,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           functionalData::AbstractIntegralFunctional,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}

  error("Generic fallback for evalFunctionalDeriv_q() reached: did you forget to extend evalFunctionalDeriv_q() with a new method for you AbstractSolutionData")

end


"""
  This function recalculates all quantities that depend on the mesh metrics.
  This function handles changes in the metric *values* only (not changed in
  mesh topology).

  This function should be called after the mesh is warped (ie. using mesh
  movement).

  Every physics module should extend this function with a new method
  specializing the `eqn` argument type.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function updateMetricDependents(mesh::AbstractMesh, sbp::AbstractSBP,
                                 eqn::AbstractSolutionData, opts)


  error("Generic fallback for updateMetricDependents() reached: did you forget to extend updateMetricDependents() with a new method for you AbstractSolutionData")

  return nothing
end


"""
  This function takes fully initialzed objects and solves the partial
  differential equation.

  Every physics should extend this function with a new method, specialzing the
  `eqn` argument type.

  Objects can be created with the [`createObjects`](@ref) function.

  The functions in the `SolverCommon` module are helpful in writing
  this function.

  **Inputs**
  
   * mesh: an AbstractMesh
   * sbp: an SBP Operator
   * eqn: an AbstractSolutionData
   * opts: the options dictionary
   * pmesh: mesh for preconditioning (optional)

  **Outputs**

   * mesh
   * sbp
   * eqn: on exit, eqn.q_vec should have the converged solution in it.
   * opts

  **Options Keys**

   * Relfunc_name: also writes vtk files called "solution_relfunc"
                  if key not present, ignored
                   TODO: fix that
   * IC_name
   * calc_error: also write vtk files called "solution_error"
   * calc_trunc_error
   * perturb_ic
   * calc_dt

    For options like calc_dt and Relfunc_name, it is very important that
    the computed quantity be saved to the options dictionary for use later
    in the code (ie. and not passed directly to another function).  The
    code won't restart correctly if this happens.


"""
function solvePDE(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                  opts::Dict, pmesh::AbstractMesh=mesh)

  error("Generic fallback for solvePDE() reaches: did you forget to extend solvePDE with a new method for your AbstractSolutionData?")

end


"""
  This function performs the reverse mode calculation of the residual with
  respect to the metrics.  The derivative should be evaluated at the solution
  in `eqn.q`.  The following fields of the mesh
  are updated with the result:

   * dxidx_bar
   * jac_bar
   * nrm_bndry_bar
   * nrm_face_bar
   * coords_bndry_bar

  It is the caller's responsibility to zero out these arrays before calling
  this function, if required.  Mesh implementations should provide a
  `zeroBarArrays` function to do this.

  Note: the mesh and equation objects must have the `need_adjoint` option set
  to true in the options dictionary used to create them.

  This function must finish parallel communication for `eqn.shared_data` and
  `eqn.shared_data_bar` by calling either [`finishExchangeData`](@ref) or
  [`finishExchangeData_rev`](@ref).  The parallel data setting of
  `eqn.shared_data` and `eqn.shared_data_bar` are always the same, and
  (currently) always equal to PARALLEL_DATA_ELEMENT.

  **Inputs**

   * mesh:
   * sbp
   * eqn: when this function is called, `eqn.res_bar` must already have the 3D
          array form of the seed vector in it
   * opts
   * t: optional argument, gives current time value

  See the other method of this function if the seed vector is in 1D form.
"""
function evalResidual_revm(mesh::AbstractMesh, sbp::AbstractSBP,
                           eqn::AbstractSolutionData, opts::Dict, t=0.0)

  error("Generic fallback for evalResidual_revm: did you forget to extend evalResidual_revm with a new method for your AbstractSolutionData")

end


"""
  This function performs reverse mode calculation of the residual with
  respect to the solution.  The `eqn.q_bar` array is updated with the result.
  It is the callers responsibility to zero out this array before calling this
  function, if required.

  Note: the mesh and equation objects must have the `need_adjoint` option set
  to true in the options dictionary used to create them.

  This function must finish parallel communication for `eqn.shared_data` and
  `eqn.shared_data_bar` by calling [`finishExchangeData_rev2`](@ref).
  The parallel data setting of`eqn.shared_data` and `eqn.shared_data_bar` are
  always the same.

  In this function `eqn.shared_data_bar` contains the part of `res_bar`
  that is sent in parallel, not `q_bar`.  This allows each process to compute
  the *local* contribution to `eqn.q_bar`, and makes hte parallel communication
  more efficient.  This function should call [`finishExchangeData_rev2`](@ref)
  to finish the parallel comunication.


  **Inputs**

   * mesh
   * sbp
   * eqn: `eqn.res_bar` must have the 3D array form of the seed vector in it
           when this function is called
   * opts
   * t: optional argument giving current time value

"""
function evalResidual_revq(mesh::AbstractMesh, sbp::AbstractSBP,
                           eqn::AbstractSolutionData, opts::Dict, t=0.0)


  error("Generic fallback for evalResidual_revq: did you forget to extend evalResidual_revq with a new method for your AbstractSolutionData")

end
