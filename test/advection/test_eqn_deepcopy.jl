# call order:
# TICON/test/runtests.jl
#  -> TICON/test/advection/runtests.jl

# This test suite is intended to test the advection physics module's 
#   implementation of eqn_deepcopy

# 5 tests:
#   1. ensure eqn_copy matches eqn
#   2. ensure eqn_copy pointers do not equal eqn's
#   3. ensure that changes to eqn do not affect eqn_copy
#   4. ensure that changes to eqn_copy.q also change q_vec
#   5. ensure that changes to eqn_copy.res also change res_vec

global const test_eqn_copy_inputfile = "input_vals_channel_DG_foreqncopy.jl"

function test_eqn_copy(mesh, sbp, eqn, opts)

  # should we step through RK4 or CN a few steps before trying eqn_deepcopy?

  eqn_copy = eqn_deepcopy(mesh, sbp, eqn, opts)

  facts("--- Testing eqn_deepcopy, DG q/q_vec & res/res_vec equality ---") do

    @fact mesh.isDG --> true
    @fact (pointer(eqn.q) == pointer(eqn.q_vec)) --> true
    @fact (pointer(eqn.res) == pointer(eqn.res_vec)) --> true
    @fact (pointer(eqn_copy.q) == pointer(eqn_copy.q_vec)) --> true
    @fact (pointer(eqn_copy.res) == pointer(eqn_copy.res_vec)) --> true

  end

  facts("--- Testing eqn_deepcopy, copy phase, values ---") do

    @fact (eqn_copy.commsize == eqn.commsize) --> true
    @fact (eqn_copy.myrank == eqn.myrank) --> true
    @fact (eqn_copy.t == eqn.t) --> true
    @fact (eqn_copy.res_type == eqn.res_type) --> true
    @fact (eqn_copy.multiplyA0inv == eqn.multiplyA0inv) --> true
    @fact (eqn_copy.params.LFalpha == eqn.params.LFalpha) --> true
    @fact (eqn_copy.params.alpha_x == eqn.params.alpha_x) --> true
    @fact (eqn_copy.params.alpha_y == eqn.params.alpha_y) --> true
    @fact (eqn_copy.params.alpha_z == eqn.params.alpha_z) --> true
    # commented out since this is an artifact of some uadj adv BC code used to test.
    # @fact (eqn_copy.params.sin_amplitude == eqn.params.sin_amplitude) --> true
    @fact (eqn_copy.params.f == eqn.params.f) --> true
    @fact (eqn_copy.params.time == eqn.params.time) --> true
    @fact (eqn_copy.src_func == eqn.src_func) --> true
    @fact (eqn_copy.flux_func == eqn.flux_func) --> true
    @fact (eqn_copy.majorIterationCallback == eqn.majorIterationCallback) --> true

    # arrays
    @fact (eqn_copy.q_face[1] == eqn.q_face[1]) --> true
    @fact (eqn_copy.flux_parametric[1] == eqn.flux_parametric[1]) --> true
    @fact (eqn_copy.flux_face[1] == eqn.flux_face[1]) --> true
    @fact (eqn_copy.q_bndry[1] == eqn.q_bndry[1]) --> true
    @fact (eqn_copy.bndryflux[1] == eqn.bndryflux[1]) --> true
    @fact (eqn_copy.M[1] == eqn.M[1]) --> true
    @fact (eqn_copy.Minv[1] == eqn.Minv[1]) --> true
    @fact (eqn_copy.Minv3D[1] == eqn.Minv3D[1]) --> true
    @fact (eqn_copy.q[1] == eqn.q[1]) --> true
    @fact (eqn_copy.q_vec[1] == eqn.q_vec[1]) --> true
    @fact (eqn_copy.res[1] == eqn.res[1]) --> true
    @fact (eqn_copy.res_vec[1] == eqn.res_vec[1]) --> true

    # Zero length arrays?
    if length(eqn.aux_vars) > 0
      @fact (eqn_copy.aux_vars[1] == eqn.aux_vars[1]) --> true
    end
    if length(eqn.res_edge) > 0
      @fact (eqn_copy.res_edge[1] == eqn.res_edge[1]) --> true
    end
    # note: q_face_send and q_face_recv are now in ParallelData, see jcrean for details. 
    #       Their purpose is now served by eqn.shared_data
    if length(eqn.shared_data) > 0
      @fact (eqn_copy.shared_data[1] == eqn.shared_data[1]) --> true
    end
    if length(eqn.flux_sharedface) > 0
      @fact (eqn_copy.flux_sharedface[1] == eqn.flux_sharedface[1]) --> true
    end
    
  end

  facts("--- Testing eqn_deepcopy, copy phase, pointers of arrays ---") do

    @fact (pointer(eqn_copy.q_face) == pointer(eqn.q_face)) --> false
    @fact (pointer(eqn_copy.aux_vars) == pointer(eqn.aux_vars)) --> false
    @fact (pointer(eqn_copy.flux_parametric) == pointer(eqn.flux_parametric)) --> false
    @fact (pointer(eqn_copy.flux_face) == pointer(eqn.flux_face)) --> false
    @fact (pointer(eqn_copy.res_edge) == pointer(eqn.res_edge)) --> false
    @fact (pointer(eqn_copy.q_bndry) == pointer(eqn.q_bndry)) --> false
    @fact (pointer(eqn_copy.bndryflux) == pointer(eqn.bndryflux)) --> false
    @fact (pointer(eqn_copy.M) == pointer(eqn.M)) --> false
    @fact (pointer(eqn_copy.Minv) == pointer(eqn.Minv)) --> false
    @fact (pointer(eqn_copy.Minv3D) == pointer(eqn.Minv3D)) --> false
    # note: q_face_send and q_face_recv are now in ParallelData, see jcrean for details. 
    #       Their purpose is now served by eqn.shared_data
    @fact (pointer(eqn_copy.shared_data) == pointer(eqn.shared_data)) --> false
    @fact (pointer(eqn_copy.flux_sharedface) == pointer(eqn.flux_sharedface)) --> false
    @fact (pointer(eqn_copy.q) == pointer(eqn.q)) --> false
    @fact (pointer(eqn_copy.q_vec) == pointer(eqn.q_vec)) --> false
    @fact (pointer(eqn_copy.res) == pointer(eqn.res)) --> false
    @fact (pointer(eqn_copy.res_vec) == pointer(eqn.res_vec)) --> false
      
  end

  rand!(eqn.q_face)
  rand!(eqn.flux_parametric)
  rand!(eqn.flux_face)
  rand!(eqn.q_bndry)
  rand!(eqn.bndryflux)
  rand!(eqn.M)
  rand!(eqn.Minv)
  rand!(eqn.Minv3D)
  rand!(eqn.q)
  rand!(eqn.res)

  facts("--- Testing eqn_deepcopy, eqn change phase ---") do

    @fact (eqn_copy.q_face[1] == eqn.q_face[1]) --> false
    @fact (eqn_copy.flux_parametric[1] == eqn.flux_parametric[1]) --> false
    @fact (eqn_copy.flux_face[1] == eqn.flux_face[1]) --> false
    @fact (eqn_copy.q_bndry[1] == eqn.q_bndry[1]) --> false
    @fact (eqn_copy.bndryflux[1] == eqn.bndryflux[1]) --> false
    @fact (eqn_copy.M[1] == eqn.M[1]) --> false
    @fact (eqn_copy.Minv[1] == eqn.Minv[1]) --> false
    @fact (eqn_copy.Minv3D[1] == eqn.Minv3D[1]) --> false
    @fact (eqn_copy.q[1] == eqn.q[1]) --> false
    @fact (eqn_copy.q_vec[1] == eqn.q_vec[1]) --> false
    @fact (eqn_copy.res[1] == eqn.res[1]) --> false
    @fact (eqn_copy.res_vec[1] == eqn.res_vec[1]) --> false

  end

  rand!(eqn_copy.q)
  rand!(eqn_copy.res)

  facts("--- Testing eqn_deepcopy, eqn_copy change phase ---") do

    @fact (eqn_copy.q[1] == eqn_copy.q_vec[1]) --> true
    @fact (eqn_copy.res[1] == eqn_copy.res_vec[1]) --> true

  end

end   # end of function test_eqn_copy

add_func2!(AdvectionTests, test_eqn_copy, test_eqn_copy_inputfile, [TAG_SHORTTEST])

#=
Fields of advection eqn object:

  commsize
  myrank
  t
  res_type
  multiplyA0inv
  params.LFalpha
  params.alpha_x
  params.alpha_y
  params.alpha_z
  params.sin_amplitude
  params.f
  params.time
  src_func
  flux_func
  majorIterationCallback

  # arrays
  q_face
  aux_vars
  flux_parametric
  flux_face
  res_edge
  q_bndry
  bndryflux
  M
  Minv
  Minv3D
  q_face_send
  q_face_recv
  flux_sharedface
  q
  q_vec
  res
  res_vec
  
=#
