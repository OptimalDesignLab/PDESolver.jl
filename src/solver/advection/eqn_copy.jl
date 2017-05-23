# export eqn_deepcopy
import PDESolver.eqn_deepcopy

# this is because of this case:
#   a = rand(2,2)
#   b = a
#   a[3] = 8
#   b[3] == 8
#   this is because 'a[3] =' is actually setindex!
"""
  AdvectionEquationMod.eqn_deepcopy

  This function performs a proper deepcopy (unlike julia's builtin deepcopy) on an equation object.
  It preserves reference topology (i.e. q & q_vec pointing to same array in DG schemes).

    Inputs:
      eqn
      mesh
      sbp
      opts

    Outputs:
      eqn_copy

"""
# TODO: check q/q_vec, res/res_vec
# TODO: tests
#  1. ensure copy matches
#  2. ensure changes to eqn don't affect eqn_copy
#  3. ensure changes to eqn_copy.q change eqn_copy.q_vec, same for res

function eqn_deepcopy{Tmsh, Tsol, Tres, Tdim}(eqn::AdvectionData_{Tsol, Tres, Tdim}, mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, opts::Dict)

  # 1: call constructor on eqn_copy

  # eqn_copy = AdvectionData_{Tsol, Tres, Tdim, Tmsh, var_type}(mesh, sbp, opts)
  eqn_copy = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)

  # should be zero'ed out

  # 2: copy over fields

  # immutable fields
  eqn_copy.commsize = eqn.commsize
  eqn_copy.myrank = eqn.myrank
  eqn_copy.t = eqn.t

  # mutable fields that can't use copy
  #   (confirmed in REPL that changes to eqn do not affect eqn_copy after assignment)
  eqn_copy.res_type = eqn.res_type
  eqn_copy.disassembleSolution = eqn.disassembleSolution
  eqn_copy.assembleSolution = eqn.assembleSolution
  eqn_copy.multiplyA0inv = eqn.multiplyA0inv

  # mutable array fields
  eqn_copy.q_face = copy(eqn.q_face)
  eqn_copy.aux_vars = copy(eqn.aux_vars)
  eqn_copy.flux_parametric = copy(eqn.flux_parametric)
  eqn_copy.flux_face = copy(eqn.flux_face)
  eqn_copy.res_edge = copy(eqn.res_edge)
  eqn_copy.q_bndry = copy(eqn.q_bndry)
  eqn_copy.bndryflux = copy(eqn.bndryflux)
  eqn_copy.M = copy(eqn.M)
  eqn_copy.Minv = copy(eqn.Minv)
  eqn_copy.Minv3D = copy(eqn.Minv3D)

  # special case: array of arrays
  eqn_copy.q_face_send = copy(eqn.q_face_send)
  for i = 1:length(eqn.q_face_send)
    eqn_copy.q_face_send[i] = copy(eqn.q_face_send[i])
  end
  eqn_copy.q_face_recv = copy(eqn.q_face_recv)
  for i = 1:length(eqn.q_face_recv)
    eqn_copy.q_face_recv[i] = copy(eqn.q_face_send[i])
  end
  if mesh.isDG
    eqn_copy.flux_sharedface = copy(eqn.flux_sharedface)
    for i = 1:length(eqn.flux_sharedface)
      eqn_copy.flux_sharedface[i] = copy(eqn.flux_sharedface[i])
    end
  end

  # Notes about DG/CG:
  #   only defined for DG:
  #     flux_sharedface
  #     flux_func?
  
  # special mutable fields
  eqn_copy.q = copy(eqn.q)
  eqn_copy.res = copy(eqn.res)
  # need to preserve topology -> don't copy over q_vec & res_vec
  if mesh.isDG
    eqn_copy.q_vec = reshape(eqn_copy.q, mesh.numDof)
    eqn_copy.res_vec = reshape(eqn_copy.res, mesh.numDof)
  else
    eqn_copy.q_vec = copy(eqn.q_vec)
    eqn_copy.res_vec = copy(eqn.res_vec)
  end

  # params type: immutable
  eqn_copy.params.LFalpha = eqn.params.LFalpha
  eqn_copy.params.alpha_x = eqn.params.alpha_x
  eqn_copy.params.alpha_y = eqn.params.alpha_y
  eqn_copy.params.alpha_z = eqn.params.alpha_z
  eqn_copy.params.sin_amplitude = eqn.params.sin_amplitude
  eqn_copy.params.f = eqn.params.f
  eqn_copy.params.time = eqn.params.time

  # other fields, mutable but ok to use assignment operator
  eqn_copy.src_func = eqn.src_func
  if mesh.isDG
    eqn_copy.flux_func = eqn.flux_func
  end
  eqn_copy.majorIterationCallback = eqn.majorIterationCallback

  return eqn_copy

end

