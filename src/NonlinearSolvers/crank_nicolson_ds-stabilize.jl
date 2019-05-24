
function stabilizeCNDSLO(lo, mesh, sbp, eqn, opts, ctx_residual, t)

  println(BSTDOUT, "        stabilizeCNDSLO called")

  # TODO: need to put stab_assembler in ctx_residual

  evalJacobianStrong(mesh, sbp, eqn, opts, stab_assembler, t)

  # TODO: real(tmp_imag)?
  # TODO: clipJacData?
  # TODO: stab_A?

  filterDiagJac(mesh, opts, real(tmp_imag), clipJacData, stab_A, eigs_to_remove="neg")

  # TODO: here is where we then add the stabilizer to the lo

  return numberOfEigchgs

end
