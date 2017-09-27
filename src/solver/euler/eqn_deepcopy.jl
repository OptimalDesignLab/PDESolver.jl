import ODLCommonTools.eqn_deepcopy
"""
  EulerEquationMod.eqn_deepcopy

  This function performs a proper deepcopy (unlike julia's builtin deepcopy) 
    on an Euler equation object.
  It preserves reference topology (i.e. q & q_vec pointing to same array in DG schemes).

    Inputs:
      eqn
      mesh
      sbp
      opts

    Outputs:
      eqn_copy

    One reason for doing this is this case:
      a = rand(2,2)
      b = a
      a[3] = 8
      b[3] == 8
      this is because 'a[3] =' is actually setindex!

"""
function eqn_deepcopy{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData_{Tsol, Tres}, opts::Dict)

  # The eqn object has over 100 fields, so it is necessary to write a better approach for 
  #   copying than explicitly copying every named field
  # This is the second write of eqn_deepcopy; the first version explicitly copied each field. 
  #   Was done for SimpleODE and Advection before Euler made it obvious that that was intractable.

  # 1: call constructor on eqn_copy

  param_tuple = getAllTypeParams(mesh, eqn, opts)

  eqn_copy = EulerData_{param_tuple...}(mesh, sbp, opts, open_files=false)

  eqn_copy = eqn_deepcopy_fields(mesh, sbp, eqn, eqn_copy, opts)

  return eqn_copy

end


