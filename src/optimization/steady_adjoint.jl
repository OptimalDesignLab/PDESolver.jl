# calculation of steady adjoint

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

   * recompute_jac: recompute the linear operator inside `ls`, default false
   * recompute_pc: recompute the preconditioner inside `ls`, default false
"""
function calcAdjoint{Tmsh, Tsol, Tres}(mesh::AbstractDGMesh{Tmsh},
                  sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol, Tres}, opts,
                  ls::LinearSolver, functionalData::AbstractFunctional,
                  adjoint_vec::Array{Tsol,1}; recompute_jac=false,
                  recompute_pc=false)
 

  # recompute operators if requested
  if recompute_jac && recompute_pc
    calcPCandLO(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0, start_comm=true)
  elseif recompute_jac
    calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0, start_comm=true)
  elseif recompute_pc
    calcPC(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0, start_comm=true)
  end

  #TODO: store this somewhere?
  func_deriv = zeros(Tsol, mesh.numDof)

  # 3D array into which func_deriv_arr gets interpolated
  func_deriv_arr = zeros(eqn.q)

  # Calculate df/dq_bndry on edges where the functional is calculated and put
  # it back in func_deriv_arr
  evalFunctionalDeriv(mesh, sbp, eqn, opts, functionalData, func_deriv_arr)

  # Assemble func_deriv
  assembleSolution(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)
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
