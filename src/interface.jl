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
  doc/interfaces.md for details

  Inputs:
    mesh: an AbstractMesh describing the mesh on which to solve the physics
    sbp: an SBP operator
    eqn: a subtype of AbstractSolution data, used to store all of the data used by the physics module
    opts: the options dictionary
    t: the current time value, defaults to 0.0
"""
function evalResidual(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict, t = 0.0)

  throw(ErrorException("Generic fallback evalResidual reached: did you forget to extend evalResidual with a new method for your AbstractSolutionData?"))

  return nothing
end

function init(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict, pmesh=mesh)

  throw(ErrorException("Generic fallback init() reached: did you forget to extend init() with a new method for your AbstractSolutionData?"))

  return nothing

end

function eqn_deepcopy(eqn::AbstractSolutionData, mesh::AbstractMesh, sbp::AbstractSBP, opts::Dict)

  throw(ErrorException("Generic fallback eqn_deepcopy() reached: eqn_deepcopy() needs to be extended with a new method for each physics module."))

  return nothing

end
