# calculation of steady adjoint

using PETSc2

"""
  Use an existing [`LinearSolver`](@ref) to compute the adjoint.
  A preconditioner and linear operator can be obtained from 
  [`NonlinearSolver.getNewtonPCandLO`](@ref).

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * ls: a LinearSolver object
   * functionalData: an [`AbstractFunctional`](@ref) to compute the adjoint of
   
  **Inputs/Outputs**

   * adjoint_vec: Array{Tsol, 1} of length mesh.numDof

  **Keyword Arguments**

   * recalc_jac: recalc the linear operator inside `ls`, default false
   * recalc_pc: recalc the preconditioner inside `ls`, default false
   * start_comm: do parallel communication before recomputing the pc and/or jac,
                 default false
"""
function calcAdjoint{Tmsh, Tsol, Tres}(mesh::AbstractDGMesh{Tmsh},
                  sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol, Tres}, opts,
                  ls::LinearSolver, functionalData::AbstractFunctional,
                  adjoint_vec::Array{Tsol,1}; recalc_jac=false,
                  recalc_pc=false, start_comm=false)
 

  # recalc operators if requested
  ctx_residual = (evalResidual,)
  if recalc_jac && recalc_pc
    calcPCandLO(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0, start_comm=start_comm)
  elseif recalc_jac
    calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0, start_comm=start_comm)
  elseif recalc_pc
    calcPC(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0, start_comm=start_comm)
  end

  #TODO: store this somewhere?
  func_deriv = zeros(Tsol, mesh.numDof)

  # 3D array into which func_deriv_arr gets interpolated
  func_deriv_arr = zeros(eqn.q)

  # Calculate df/dq_bndry on edges where the functional is calculated and put
  # it back in func_deriv_arr
  evalFunctionalDeriv(mesh, sbp, eqn, opts, functionalData, func_deriv_arr)

  # Assemble func_deriv
  array3DTo1D(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)
  scale!(func_deriv, -1.0)

  # do transpose solve
  _adjoint_vec = zeros(real(Tsol), length(adjoint_vec))
  linearSolveTranspose(ls, real(func_deriv), _adjoint_vec)
  copy!(adjoint_vec, _adjoint_vec)

  # Output/Visualization options for Adjoint
  if opts["write_adjoint"]
    outname = string("adjoint_vec_", mesh.myrank,".dat")
    f = open(outname, "w")
    for i = 1:length(adjoint_vec)
      println(f, real(adjoint_vec[i]))
    end
    close(f)
  end
  
  if opts["write_adjoint_vis"]
    saveSolutionToMesh(mesh, adjoint_vec)
    fname = "adjoint_field"
    writeVisFiles(mesh, fname)
  end

  return nothing
end
