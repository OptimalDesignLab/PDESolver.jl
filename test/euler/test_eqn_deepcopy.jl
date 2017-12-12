# call order:
# TICON/test/runtests.jl
#   -> TICON/test/euler/runtests.jl
#
# This test suite is intended to test the euler physics module's implementation of eqn_deepcopy
#
# 5 types of tests:
#   1. ensure eqn_copy matches eqn
#   2. ensure eqn_copy pointers do not equal eqn's (for array fields)
#   3. ensure that changes to eqn do not affect eqn copy
#   4. ensure that changes to eqn_copy.q also change q_vec
#   5. ensure that changes to eqn_copy.res also change res_vec

# TODO? CG tests

# Note about params/params_conservative/params_entropy:
#   Tests are done below for params and params_entropy, and not params_conservative,
#   since there shouldn't be a case where params_conservative works and params_entropy doesn't

global const test_eqn_copy_inputfile = "input_vals_channel_dg.jl"

function test_eqn_copy()

  mesh, sbp, eqn, opts = run_euler(test_eqn_copy_inputfile)

  eqn_copy = eqn_deepcopy(mesh, sbp, eqn, opts)

  facts("--- Testing eqn_deepcopy, DG q/q_vec & res/res_vec equality ---") do

    @fact mesh.isDG --> true
    @fact (pointer(eqn.q) == pointer(eqn.q_vec)) --> true
    @fact (pointer(eqn.res) == pointer(eqn.res_vec)) --> true
    @fact (pointer(eqn_copy.q) == pointer(eqn_copy.q_vec)) --> true
    @fact (pointer(eqn_copy.res) == pointer(eqn_copy.res_vec)) --> true

  end

  #----------------------
  # Test type #1: eqn_copy field values equal eqn's
  facts("--- Testing eqn_deepcopy, copy phase, values ---") do

    # testing 1st level

    # note: these three commented out because the double equals operator reverts back to triple equals 
    #   (test if they're same object). So the double equals will return false
    # @fact (eqn_copy.params == eqn.params) --> true
    # @fact (eqn_copy.params_conservative == eqn.params_conservative) --> true
    # @fact (eqn_copy.params_entropy == eqn.params_entropy) --> true
    @fact (eqn_copy.comm == eqn.comm) --> true
    @fact (eqn_copy.commsize == eqn.commsize) --> true
    @fact (eqn_copy.myrank == eqn.myrank) --> true
    @fact (eqn_copy.q == eqn.q) --> true
    @fact (eqn_copy.q_face == eqn.q_face) --> true
    @fact (eqn_copy.q_bndry == eqn.q_bndry) --> true
    @fact (eqn_copy.q_vec == eqn.q_vec) --> true
    @fact (eqn_copy.aux_vars == eqn.aux_vars) --> true
    @fact (eqn_copy.aux_vars_face == eqn.aux_vars_face) --> true
    @fact (eqn_copy.aux_vars_sharedface == eqn.aux_vars_sharedface) --> true
    @fact (eqn_copy.aux_vars_bndry == eqn.aux_vars_bndry) --> true
    @fact (eqn_copy.flux_parametric == eqn.flux_parametric) --> true
    # note: q_face_send and q_face_recv purposes are now served by eqn.shared_data
    @fact (eqn_copy.shared_data == eqn.shared_data) --> true
    @fact (eqn_copy.flux_face == eqn.flux_face) --> true
    @fact (eqn_copy.flux_sharedface == eqn.flux_sharedface) --> true
    @fact (eqn_copy.res == eqn.res) --> true
    @fact (eqn_copy.res_vec == eqn.res_vec) --> true
    @fact (eqn_copy.Axi == eqn.Axi) --> true
    @fact (eqn_copy.Aeta == eqn.Aeta) --> true
    @fact (eqn_copy.res_edge == eqn.res_edge) --> true
    @fact (eqn_copy.edgestab_alpha == eqn.edgestab_alpha) --> true
    @fact (eqn_copy.bndryflux == eqn.bndryflux) --> true
    @fact (eqn_copy.stabscale == eqn.stabscale) --> true
    @fact (eqn_copy.Minv3D == eqn.Minv3D) --> true
    @fact (eqn_copy.Minv == eqn.Minv) --> true
    @fact (eqn_copy.M == eqn.M) --> true
    @fact (eqn_copy.multiplyA0inv == eqn.multiplyA0inv) --> true
    @fact (eqn_copy.majorIterationCallback == eqn.majorIterationCallback) --> true
    @fact (eqn_copy.src_func == eqn.src_func) --> true
    @fact (eqn_copy.flux_func == eqn.flux_func) --> true
    @fact (eqn_copy.volume_flux_func == eqn.volume_flux_func) --> true
    @fact (eqn_copy.face_element_integral_func == eqn.face_element_integral_func) --> true

    # testing 2nd level, params
    @fact (eqn_copy.params.f == eqn.params.f) --> true
    @fact (eqn_copy.params.t == eqn.params.t) --> true
    @fact (eqn_copy.params.order == eqn.params.order) --> true
    @fact (eqn_copy.params.q_vals == eqn.params.q_vals) --> true
    @fact (eqn_copy.params.q_vals2 == eqn.params.q_vals2) --> true
    @fact (eqn_copy.params.q_vals3 == eqn.params.q_vals3) --> true
    @fact (eqn_copy.params.qg == eqn.params.qg) --> true
    @fact (eqn_copy.params.v_vals == eqn.params.v_vals) --> true
    @fact (eqn_copy.params.v_vals2 == eqn.params.v_vals2) --> true
    @fact (eqn_copy.params.Lambda == eqn.params.Lambda) --> true
    @fact (eqn_copy.params.w_vals_stencil == eqn.params.w_vals_stencil) --> true
    @fact (eqn_copy.params.w_vals2_stencil == eqn.params.w_vals2_stencil) --> true
    @fact (eqn_copy.params.res_vals1 == eqn.params.res_vals1) --> true
    @fact (eqn_copy.params.res_vals2 == eqn.params.res_vals2) --> true
    @fact (eqn_copy.params.res_vals3 == eqn.params.res_vals3) --> true
    @fact (eqn_copy.params.flux_vals1 == eqn.params.flux_vals1) --> true
    @fact (eqn_copy.params.flux_vals2 == eqn.params.flux_vals2) --> true
    @fact (eqn_copy.params.sat_vals == eqn.params.sat_vals) --> true
    @fact (eqn_copy.params.A0 == eqn.params.A0) --> true
    @fact (eqn_copy.params.A0inv == eqn.params.A0inv) --> true
    @fact (eqn_copy.params.A1 == eqn.params.A1) --> true
    @fact (eqn_copy.params.A2 == eqn.params.A2) --> true
    @fact (eqn_copy.params.S2 == eqn.params.S2) --> true
    @fact (eqn_copy.params.A_mats == eqn.params.A_mats) --> true
    @fact (eqn_copy.params.Rmat1 == eqn.params.Rmat1) --> true
    @fact (eqn_copy.params.Rmat2 == eqn.params.Rmat2) --> true
    @fact (eqn_copy.params.P == eqn.params.P) --> true
    @fact (eqn_copy.params.nrm == eqn.params.nrm) --> true
    @fact (eqn_copy.params.nrm2 == eqn.params.nrm2) --> true
    @fact (eqn_copy.params.nrm3 == eqn.params.nrm3) --> true
    @fact (eqn_copy.params.h == eqn.params.h) --> true
    @fact (eqn_copy.params.cv == eqn.params.cv) --> true
    @fact (eqn_copy.params.R == eqn.params.R) --> true
    @fact (eqn_copy.params.gamma == eqn.params.gamma) --> true
    @fact (eqn_copy.params.gamma_1 == eqn.params.gamma_1) --> true
    @fact (eqn_copy.params.Ma == eqn.params.Ma) --> true
    @fact (eqn_copy.params.Re == eqn.params.Re) --> true
    @fact (eqn_copy.params.aoa == eqn.params.aoa) --> true
    @fact (eqn_copy.params.rho_free == eqn.params.rho_free) --> true
    @fact (eqn_copy.params.E_free == eqn.params.E_free) --> true
    @fact (eqn_copy.params.edgestab_gamma == eqn.params.edgestab_gamma) --> true
    @fact (eqn_copy.params.writeflux == eqn.params.writeflux) --> true
    @fact (eqn_copy.params.writeboundary == eqn.params.writeboundary) --> true
    @fact (eqn_copy.params.writeq == eqn.params.writeq) --> true
    @fact (eqn_copy.params.use_edgestab == eqn.params.use_edgestab) --> true
    @fact (eqn_copy.params.use_filter == eqn.params.use_filter) --> true
    @fact (eqn_copy.params.use_res_filter == eqn.params.use_res_filter) --> true
    @fact (eqn_copy.params.filter_mat == eqn.params.filter_mat) --> true
    @fact (eqn_copy.params.use_dissipation == eqn.params.use_dissipation) --> true
    @fact (eqn_copy.params.dissipation_const == eqn.params.dissipation_const) --> true
    @fact (eqn_copy.params.tau_type == eqn.params.tau_type) --> true
    @fact (eqn_copy.params.vortex_x0 == eqn.params.vortex_x0) --> true
    @fact (eqn_copy.params.vortex_strength == eqn.params.vortex_strength) --> true
    @fact (eqn_copy.params.krylov_itr == eqn.params.krylov_itr) --> true
    @fact (eqn_copy.params.krylov_type == eqn.params.krylov_type) --> true
    @fact (eqn_copy.params.Rprime == eqn.params.Rprime) --> true
    @fact (eqn_copy.params.A == eqn.params.A) --> true
    @fact (eqn_copy.params.B == eqn.params.B) --> true
    @fact (eqn_copy.params.iperm == eqn.params.iperm) --> true
    @fact (eqn_copy.params.time == eqn.params.time) --> true

    # testing 2nd level, params_entropy
    @fact (eqn_copy.params_entropy.f == eqn.params_entropy.f) --> true
    @fact (eqn_copy.params_entropy.t == eqn.params_entropy.t) --> true
    @fact (eqn_copy.params_entropy.order == eqn.params_entropy.order) --> true
    @fact (eqn_copy.params_entropy.q_vals == eqn.params_entropy.q_vals) --> true
    @fact (eqn_copy.params_entropy.q_vals2 == eqn.params_entropy.q_vals2) --> true
    @fact (eqn_copy.params_entropy.q_vals3 == eqn.params_entropy.q_vals3) --> true
    @fact (eqn_copy.params_entropy.qg == eqn.params_entropy.qg) --> true
    @fact (eqn_copy.params_entropy.v_vals == eqn.params_entropy.v_vals) --> true
    @fact (eqn_copy.params_entropy.v_vals2 == eqn.params_entropy.v_vals2) --> true
    @fact (eqn_copy.params_entropy.Lambda == eqn.params_entropy.Lambda) --> true
    @fact (eqn_copy.params_entropy.w_vals_stencil == eqn.params_entropy.w_vals_stencil) --> true
    @fact (eqn_copy.params_entropy.w_vals2_stencil == eqn.params_entropy.w_vals2_stencil) --> true
    @fact (eqn_copy.params_entropy.res_vals1 == eqn.params_entropy.res_vals1) --> true
    @fact (eqn_copy.params_entropy.res_vals2 == eqn.params_entropy.res_vals2) --> true
    @fact (eqn_copy.params_entropy.res_vals3 == eqn.params_entropy.res_vals3) --> true
    @fact (eqn_copy.params_entropy.flux_vals1 == eqn.params_entropy.flux_vals1) --> true
    @fact (eqn_copy.params_entropy.flux_vals2 == eqn.params_entropy.flux_vals2) --> true
    @fact (eqn_copy.params_entropy.sat_vals == eqn.params_entropy.sat_vals) --> true
    @fact (eqn_copy.params_entropy.A0 == eqn.params_entropy.A0) --> true
    @fact (eqn_copy.params_entropy.A0inv == eqn.params_entropy.A0inv) --> true
    @fact (eqn_copy.params_entropy.A1 == eqn.params_entropy.A1) --> true
    @fact (eqn_copy.params_entropy.A2 == eqn.params_entropy.A2) --> true
    @fact (eqn_copy.params_entropy.S2 == eqn.params_entropy.S2) --> true
    @fact (eqn_copy.params_entropy.A_mats == eqn.params_entropy.A_mats) --> true
    @fact (eqn_copy.params_entropy.Rmat1 == eqn.params_entropy.Rmat1) --> true
    @fact (eqn_copy.params_entropy.Rmat2 == eqn.params_entropy.Rmat2) --> true
    @fact (eqn_copy.params_entropy.P == eqn.params_entropy.P) --> true
    @fact (eqn_copy.params_entropy.nrm == eqn.params_entropy.nrm) --> true
    @fact (eqn_copy.params_entropy.nrm2 == eqn.params_entropy.nrm2) --> true
    @fact (eqn_copy.params_entropy.nrm3 == eqn.params_entropy.nrm3) --> true
    @fact (eqn_copy.params_entropy.h == eqn.params_entropy.h) --> true
    @fact (eqn_copy.params_entropy.cv == eqn.params_entropy.cv) --> true
    @fact (eqn_copy.params_entropy.R == eqn.params_entropy.R) --> true
    @fact (eqn_copy.params_entropy.gamma == eqn.params_entropy.gamma) --> true
    @fact (eqn_copy.params_entropy.gamma_1 == eqn.params_entropy.gamma_1) --> true
    @fact (eqn_copy.params_entropy.Ma == eqn.params_entropy.Ma) --> true
    @fact (eqn_copy.params_entropy.Re == eqn.params_entropy.Re) --> true
    @fact (eqn_copy.params_entropy.aoa == eqn.params_entropy.aoa) --> true
    @fact (eqn_copy.params_entropy.rho_free == eqn.params_entropy.rho_free) --> true
    @fact (eqn_copy.params_entropy.E_free == eqn.params_entropy.E_free) --> true
    @fact (eqn_copy.params_entropy.edgestab_gamma == eqn.params_entropy.edgestab_gamma) --> true
    @fact (eqn_copy.params_entropy.writeflux == eqn.params_entropy.writeflux) --> true
    @fact (eqn_copy.params_entropy.writeboundary == eqn.params_entropy.writeboundary) --> true
    @fact (eqn_copy.params_entropy.writeq == eqn.params_entropy.writeq) --> true
    @fact (eqn_copy.params_entropy.use_edgestab == eqn.params_entropy.use_edgestab) --> true
    @fact (eqn_copy.params_entropy.use_filter == eqn.params_entropy.use_filter) --> true
    @fact (eqn_copy.params_entropy.use_res_filter == eqn.params_entropy.use_res_filter) --> true
    @fact (eqn_copy.params_entropy.filter_mat == eqn.params_entropy.filter_mat) --> true
    @fact (eqn_copy.params_entropy.use_dissipation == eqn.params_entropy.use_dissipation) --> true
    @fact (eqn_copy.params_entropy.dissipation_const == eqn.params_entropy.dissipation_const) --> true
    @fact (eqn_copy.params_entropy.tau_type == eqn.params_entropy.tau_type) --> true
    @fact (eqn_copy.params_entropy.vortex_x0 == eqn.params_entropy.vortex_x0) --> true
    @fact (eqn_copy.params_entropy.vortex_strength == eqn.params_entropy.vortex_strength) --> true
    @fact (eqn_copy.params_entropy.krylov_itr == eqn.params_entropy.krylov_itr) --> true
    @fact (eqn_copy.params_entropy.krylov_type == eqn.params_entropy.krylov_type) --> true
    @fact (eqn_copy.params_entropy.Rprime == eqn.params_entropy.Rprime) --> true
    @fact (eqn_copy.params_entropy.A == eqn.params_entropy.A) --> true
    @fact (eqn_copy.params_entropy.B == eqn.params_entropy.B) --> true
    @fact (eqn_copy.params_entropy.iperm == eqn.params_entropy.iperm) --> true
    @fact (eqn_copy.params_entropy.time == eqn.params_entropy.time) --> true

  end   # end of facts check on values of eqn_copy matching eqn

  facts("--- Testing eqn_deepcopy, copy phase, pointers of arrays ---") do

    # 1st level fields of type Array
    @fact (pointer(eqn_copy.q) == pointer(eqn.q)) --> false
    @fact (pointer(eqn_copy.q_face) == pointer(eqn.q_face)) --> false
    @fact (pointer(eqn_copy.q_bndry) == pointer(eqn.q_bndry)) --> false
    @fact (pointer(eqn_copy.q_vec) == pointer(eqn.q_vec)) --> false
    @fact (pointer(eqn_copy.aux_vars) == pointer(eqn.aux_vars)) --> false
    @fact (pointer(eqn_copy.aux_vars_face) == pointer(eqn.aux_vars_face)) --> false
    @fact (pointer(eqn_copy.aux_vars_sharedface) == pointer(eqn.aux_vars_sharedface)) --> false
    @fact (pointer(eqn_copy.aux_vars_bndry) == pointer(eqn.aux_vars_bndry)) --> false
    @fact (pointer(eqn_copy.flux_parametric) == pointer(eqn.flux_parametric)) --> false
    # note: q_face_send and q_face_recv purposes are now served by eqn.shared_data
    @fact (pointer(eqn_copy.shared_data) == pointer(eqn.shared_data)) --> false
    @fact (pointer(eqn_copy.flux_face) == pointer(eqn.flux_face)) --> false
    @fact (pointer(eqn_copy.flux_sharedface) == pointer(eqn.flux_sharedface)) --> false
    @fact (pointer(eqn_copy.res) == pointer(eqn.res)) --> false
    @fact (pointer(eqn_copy.res_vec) == pointer(eqn.res_vec)) --> false
    @fact (pointer(eqn_copy.Axi) == pointer(eqn.Axi)) --> false
    @fact (pointer(eqn_copy.Aeta) == pointer(eqn.Aeta)) --> false
    @fact (pointer(eqn_copy.res_edge) == pointer(eqn.res_edge)) --> false
    @fact (pointer(eqn_copy.edgestab_alpha) == pointer(eqn.edgestab_alpha)) --> false
    @fact (pointer(eqn_copy.bndryflux) == pointer(eqn.bndryflux)) --> false
    @fact (pointer(eqn_copy.stabscale) == pointer(eqn.stabscale)) --> false
    @fact (pointer(eqn_copy.Minv3D) == pointer(eqn.Minv3D)) --> false
    @fact (pointer(eqn_copy.Minv) == pointer(eqn.Minv)) --> false
    @fact (pointer(eqn_copy.M) == pointer(eqn.M)) --> false

    # 2nd level fields of type Array (just going to check params)
    @fact (pointer(eqn_copy.params.q_vals) == pointer(eqn.params.q_vals)) --> false
    @fact (pointer(eqn_copy.params.q_vals2) == pointer(eqn.params.q_vals2)) --> false
    @fact (pointer(eqn_copy.params.q_vals3) == pointer(eqn.params.q_vals3)) --> false
    @fact (pointer(eqn_copy.params.qg) == pointer(eqn.params.qg)) --> false
    @fact (pointer(eqn_copy.params.v_vals) == pointer(eqn.params.v_vals)) --> false
    @fact (pointer(eqn_copy.params.v_vals2) == pointer(eqn.params.v_vals2)) --> false
    @fact (pointer(eqn_copy.params.Lambda) == pointer(eqn.params.Lambda)) --> false
    @fact (pointer(eqn_copy.params.w_vals_stencil) == pointer(eqn.params.w_vals_stencil)) --> false
    @fact (pointer(eqn_copy.params.w_vals2_stencil) == pointer(eqn.params.w_vals2_stencil)) --> false
    @fact (pointer(eqn_copy.params.res_vals1) == pointer(eqn.params.res_vals1)) --> false
    @fact (pointer(eqn_copy.params.res_vals2) == pointer(eqn.params.res_vals2)) --> false
    @fact (pointer(eqn_copy.params.res_vals3) == pointer(eqn.params.res_vals3)) --> false
    @fact (pointer(eqn_copy.params.flux_vals1) == pointer(eqn.params.flux_vals1)) --> false
    @fact (pointer(eqn_copy.params.flux_vals2) == pointer(eqn.params.flux_vals2)) --> false
    @fact (pointer(eqn_copy.params.sat_vals) == pointer(eqn.params.sat_vals)) --> false
    @fact (pointer(eqn_copy.params.A0) == pointer(eqn.params.A0)) --> false
    @fact (pointer(eqn_copy.params.A0inv) == pointer(eqn.params.A0inv)) --> false
    @fact (pointer(eqn_copy.params.A1) == pointer(eqn.params.A1)) --> false
    @fact (pointer(eqn_copy.params.A2) == pointer(eqn.params.A2)) --> false
    @fact (pointer(eqn_copy.params.S2) == pointer(eqn.params.S2)) --> false
    @fact (pointer(eqn_copy.params.A_mats) == pointer(eqn.params.A_mats)) --> false
    @fact (pointer(eqn_copy.params.Rmat1) == pointer(eqn.params.Rmat1)) --> false
    @fact (pointer(eqn_copy.params.Rmat2) == pointer(eqn.params.Rmat2)) --> false
    @fact (pointer(eqn_copy.params.P) == pointer(eqn.params.P)) --> false
    @fact (pointer(eqn_copy.params.nrm) == pointer(eqn.params.nrm)) --> false
    @fact (pointer(eqn_copy.params.nrm2) == pointer(eqn.params.nrm2)) --> false
    @fact (pointer(eqn_copy.params.nrm3) == pointer(eqn.params.nrm3)) --> false
    @fact (pointer(eqn_copy.params.filter_mat) == pointer(eqn.params.filter_mat)) --> false
    @fact (pointer(eqn_copy.params.Rprime) == pointer(eqn.params.Rprime)) --> false
    @fact (pointer(eqn_copy.params.A) == pointer(eqn.params.A)) --> false
    @fact (pointer(eqn_copy.params.B) == pointer(eqn.params.B)) --> false
    @fact (pointer(eqn_copy.params.iperm) == pointer(eqn.params.iperm)) --> false

    @fact (pointer(eqn_copy.params_entropy.q_vals) == pointer(eqn.params_entropy.q_vals)) --> false
    @fact (pointer(eqn_copy.params_entropy.q_vals2) == pointer(eqn.params_entropy.q_vals2)) --> false
    @fact (pointer(eqn_copy.params_entropy.q_vals3) == pointer(eqn.params_entropy.q_vals3)) --> false
    @fact (pointer(eqn_copy.params_entropy.qg) == pointer(eqn.params_entropy.qg)) --> false
    @fact (pointer(eqn_copy.params_entropy.v_vals) == pointer(eqn.params_entropy.v_vals)) --> false
    @fact (pointer(eqn_copy.params_entropy.v_vals2) == pointer(eqn.params_entropy.v_vals2)) --> false
    @fact (pointer(eqn_copy.params_entropy.Lambda) == pointer(eqn.params_entropy.Lambda)) --> false
    @fact (pointer(eqn_copy.params_entropy.w_vals_stencil) == pointer(eqn.params_entropy.w_vals_stencil)) --> false
    @fact (pointer(eqn_copy.params_entropy.w_vals2_stencil) == pointer(eqn.params_entropy.w_vals2_stencil)) --> false
    @fact (pointer(eqn_copy.params_entropy.res_vals1) == pointer(eqn.params_entropy.res_vals1)) --> false
    @fact (pointer(eqn_copy.params_entropy.res_vals2) == pointer(eqn.params_entropy.res_vals2)) --> false
    @fact (pointer(eqn_copy.params_entropy.res_vals3) == pointer(eqn.params_entropy.res_vals3)) --> false
    @fact (pointer(eqn_copy.params_entropy.flux_vals1) == pointer(eqn.params_entropy.flux_vals1)) --> false
    @fact (pointer(eqn_copy.params_entropy.flux_vals2) == pointer(eqn.params_entropy.flux_vals2)) --> false
    @fact (pointer(eqn_copy.params_entropy.sat_vals) == pointer(eqn.params_entropy.sat_vals)) --> false
    @fact (pointer(eqn_copy.params_entropy.A0) == pointer(eqn.params_entropy.A0)) --> false
    @fact (pointer(eqn_copy.params_entropy.A0inv) == pointer(eqn.params_entropy.A0inv)) --> false
    @fact (pointer(eqn_copy.params_entropy.A1) == pointer(eqn.params_entropy.A1)) --> false
    @fact (pointer(eqn_copy.params_entropy.A2) == pointer(eqn.params_entropy.A2)) --> false
    @fact (pointer(eqn_copy.params_entropy.S2) == pointer(eqn.params_entropy.S2)) --> false
    @fact (pointer(eqn_copy.params_entropy.A_mats) == pointer(eqn.params_entropy.A_mats)) --> false
    @fact (pointer(eqn_copy.params_entropy.Rmat1) == pointer(eqn.params_entropy.Rmat1)) --> false
    @fact (pointer(eqn_copy.params_entropy.Rmat2) == pointer(eqn.params_entropy.Rmat2)) --> false
    @fact (pointer(eqn_copy.params_entropy.P) == pointer(eqn.params_entropy.P)) --> false
    @fact (pointer(eqn_copy.params_entropy.nrm) == pointer(eqn.params_entropy.nrm)) --> false
    @fact (pointer(eqn_copy.params_entropy.nrm2) == pointer(eqn.params_entropy.nrm2)) --> false
    @fact (pointer(eqn_copy.params_entropy.nrm3) == pointer(eqn.params_entropy.nrm3)) --> false
    @fact (pointer(eqn_copy.params_entropy.filter_mat) == pointer(eqn.params_entropy.filter_mat)) --> false
    @fact (pointer(eqn_copy.params_entropy.Rprime) == pointer(eqn.params_entropy.Rprime)) --> false
    @fact (pointer(eqn_copy.params_entropy.A) == pointer(eqn.params_entropy.A)) --> false
    @fact (pointer(eqn_copy.params_entropy.B) == pointer(eqn.params_entropy.B)) --> false
    @fact (pointer(eqn_copy.params_entropy.iperm) == pointer(eqn.params_entropy.iperm)) --> false

  end   # end of facts checking pointers of arrays

  rand!(eqn.q_face)
  rand!(eqn.flux_parametric)
  rand!(eqn.flux_face)
  rand!(eqn.q)
  rand!(eqn.res)
  rand!(eqn.M)
  rand!(eqn.Minv)
  rand!(eqn.Minv3D)
  rand!(eqn.bndryflux)
  rand!(eqn.q_bndry)

  facts("--- Testing eqn_deepcopy, eqn change phase ---") do

    @fact (eqn_copy.q_face[1] == eqn.q_face[1]) --> false
    @fact (eqn_copy.flux_parametric[1] == eqn.flux_parametric[1]) --> false
    @fact (eqn_copy.flux_face[1] == eqn.flux_face[1]) --> false
    @fact (eqn_copy.q[1] == eqn.q[1]) --> false
    @fact (eqn_copy.res[1] == eqn.res[1]) --> false
    @fact (eqn_copy.M[1] == eqn.M[1]) --> false
    @fact (eqn_copy.Minv[1] == eqn.Minv[1]) --> false
    @fact (eqn_copy.Minv3D[1] == eqn.Minv3D[1]) --> false
    @fact (eqn_copy.bndryflux[1] == eqn.bndryflux[1]) --> false
    @fact (eqn_copy.q_bndry[1] == eqn.q_bndry[1]) --> false

  end

  rand!(eqn_copy.q)
  rand!(eqn_copy.res)

  facts("--- Testing eqn_deepcopy, eqn_copy change phase ---") do

    @fact (eqn_copy.q[1] == eqn_copy.q_vec[1]) --> true
    @fact (eqn_copy.res[1] == eqn_copy.res_vec[1]) --> true

  end


end   # end of function test_eqn_copy

add_func1!(EulerTests, test_eqn_copy, [TAG_SHORTTEST])

#=
All fields of Euler's eqn obj, as of 20170522

params  EulerEquationMod.ParamType{2,:conservative,Float64,Float64,Float64}  false
comm  MPI.Comm  false
commsize  Int64  true
myrank  Int64  true
params_conservative  EulerEquationMod.ParamType{2,:conservative,Float64,Float64,Float64}  false
params_entropy  EulerEquationMod.ParamType{2,:entropy,Float64,Float64,Float64}  false
q  Array{Float64,3}  false
q_face  Array{Float64,4}  false
q_bndry  Array{Float64,3}  false
q_vec  Array{Float64,1}  false
aux_vars  Array{Float64,3}  false
aux_vars_face  Array{Float64,3}  false
aux_vars_sharedface  Array{Array{Float64,3},1}  false
aux_vars_bndry  Array{Float64,3}  false
flux_parametric  Array{Float64,4}  false
q_face_send  Array{Array{Float64,3},1}  false
q_face_recv  Array{Array{Float64,3},1}  false
flux_face  Array{Float64,3}  false
flux_sharedface  Array{Array{Float64,3},1}  false
res  Array{Float64,3}  false
res_vec  Array{Float64,1}  false
Axi  Array{Float64,4}  false
Aeta  Array{Float64,4}  false
res_edge  Array{Float64,4}  false
edgestab_alpha  Array{Float64,4}  false
bndryflux  Array{Float64,3}  false
stabscale  Array{Float64,2}  false
dissipation_mat  Array{Float64,3}  false
Minv3D  Array{Float64,3}  false
Minv  Array{Float64,1}  false
M  Array{Float64,1}  false
multiplyA0inv  Function  false
majorIterationCallback  Function  false
src_func  EulerEquationMod.SRC0  false
flux_func  EulerEquationMod.RoeFlux  false
volume_flux_func  EulerEquationMod.StandardFlux  false
face_element_integral_func  EulerEquationMod.ESLFFaceIntegral  false
== params ==
f  IOStream  false
t  Float64  true
order  Int64  true
q_vals  Array{Float64,1}  false
q_vals2  Array{Float64,1}  false
q_vals3  Array{Float64,1}  false
qg  Array{Float64,1}  false
v_vals  Array{Float64,1}  false
v_vals2  Array{Float64,1}  false
Lambda  Array{Float64,1}  false
w_vals_stencil  Array{Float64,2}  false
w_vals2_stencil  Array{Float64,2}  false
res_vals1  Array{Float64,1}  false
res_vals2  Array{Float64,1}  false
res_vals3  Array{Float64,1}  false
flux_vals1  Array{Float64,1}  false
flux_vals2  Array{Float64,1}  false
sat_vals  Array{Float64,1}  false
A0  Array{Float64,2}  false
A0inv  Array{Float64,2}  false
A1  Array{Float64,2}  false
A2  Array{Float64,2}  false
S2  Array{Float64,1}  false
A_mats  Array{Float64,3}  false
Rmat1  Array{Float64,2}  false
Rmat2  Array{Float64,2}  false
P  Array{Float64,2}  false
nrm  Array{Float64,1}  false
nrm2  Array{Float64,1}  false
nrm3  Array{Float64,1}  false
h  Float64  true
cv  Float64  true
R  Float64  true
gamma  Float64  true
gamma_1  Float64  true
Ma  Float64  true
Re  Float64  true
aoa  Float64  true
rho_free  Float64  true
E_free  Float64  true
edgestab_gamma  Float64  true
writeflux  Bool  true
writeboundary  Bool  true
writeq  Bool  true
use_edgestab  Bool  true
use_filter  Bool  true
use_res_filter  Bool  true
filter_mat  Array{Float64,2}  false
use_dissipation  Bool  true
dissipation_const  Float64  true
tau_type  Int64  true
vortex_x0  Float64  true
vortex_strength  Float64  true
krylov_itr  Int64  true
krylov_type  Int64  true
Rprime  Array{Float64,2}  false
A  Array{Float64,2}  false
B  Array{Float64,3}  false
iperm  Array{Int64,1}  false
time  Utils.Timings  false
== params_conservative ==
f  IOStream  false
t  Float64  true
order  Int64  true
q_vals  Array{Float64,1}  false
q_vals2  Array{Float64,1}  false
q_vals3  Array{Float64,1}  false
qg  Array{Float64,1}  false
v_vals  Array{Float64,1}  false
v_vals2  Array{Float64,1}  false
Lambda  Array{Float64,1}  false
w_vals_stencil  Array{Float64,2}  false
w_vals2_stencil  Array{Float64,2}  false
res_vals1  Array{Float64,1}  false
res_vals2  Array{Float64,1}  false
res_vals3  Array{Float64,1}  false
flux_vals1  Array{Float64,1}  false
flux_vals2  Array{Float64,1}  false
sat_vals  Array{Float64,1}  false
A0  Array{Float64,2}  false
A0inv  Array{Float64,2}  false
A1  Array{Float64,2}  false
A2  Array{Float64,2}  false
S2  Array{Float64,1}  false
A_mats  Array{Float64,3}  false
Rmat1  Array{Float64,2}  false
Rmat2  Array{Float64,2}  false
P  Array{Float64,2}  false
nrm  Array{Float64,1}  false
nrm2  Array{Float64,1}  false
nrm3  Array{Float64,1}  false
h  Float64  true
cv  Float64  true
R  Float64  true
gamma  Float64  true
gamma_1  Float64  true
Ma  Float64  true
Re  Float64  true
aoa  Float64  true
rho_free  Float64  true
E_free  Float64  true
edgestab_gamma  Float64  true
writeflux  Bool  true
writeboundary  Bool  true
writeq  Bool  true
use_edgestab  Bool  true
use_filter  Bool  true
use_res_filter  Bool  true
filter_mat  Array{Float64,2}  false
use_dissipation  Bool  true
dissipation_const  Float64  true
tau_type  Int64  true
vortex_x0  Float64  true
vortex_strength  Float64  true
krylov_itr  Int64  true
krylov_type  Int64  true
Rprime  Array{Float64,2}  false
A  Array{Float64,2}  false
B  Array{Float64,3}  false
iperm  Array{Int64,1}  false
time  Utils.Timings  false
== params_entropy ==
f  IOStream  false
t  Float64  true
order  Int64  true
q_vals  Array{Float64,1}  false
q_vals2  Array{Float64,1}  false
q_vals3  Array{Float64,1}  false
qg  Array{Float64,1}  false
v_vals  Array{Float64,1}  false
v_vals2  Array{Float64,1}  false
Lambda  Array{Float64,1}  false
w_vals_stencil  Array{Float64,2}  false
w_vals2_stencil  Array{Float64,2}  false
res_vals1  Array{Float64,1}  false
res_vals2  Array{Float64,1}  false
res_vals3  Array{Float64,1}  false
flux_vals1  Array{Float64,1}  false
flux_vals2  Array{Float64,1}  false
sat_vals  Array{Float64,1}  false
A0  Array{Float64,2}  false
A0inv  Array{Float64,2}  false
A1  Array{Float64,2}  false
A2  Array{Float64,2}  false
S2  Array{Float64,1}  false
A_mats  Array{Float64,3}  false
Rmat1  Array{Float64,2}  false
Rmat2  Array{Float64,2}  false
P  Array{Float64,2}  false
nrm  Array{Float64,1}  false
nrm2  Array{Float64,1}  false
nrm3  Array{Float64,1}  false
h  Float64  true
cv  Float64  true
R  Float64  true
gamma  Float64  true
gamma_1  Float64  true
Ma  Float64  true
Re  Float64  true
aoa  Float64  true
rho_free  Float64  true
E_free  Float64  true
edgestab_gamma  Float64  true
writeflux  Bool  true
writeboundary  Bool  true
writeq  Bool  true
use_edgestab  Bool  true
use_filter  Bool  true
use_res_filter  Bool  true
filter_mat  Array{Float64,2}  false
use_dissipation  Bool  true
dissipation_const  Float64  true
tau_type  Int64  true
vortex_x0  Float64  true
vortex_strength  Float64  true
krylov_itr  Int64  true
krylov_type  Int64  true
Rprime  Array{Float64,2}  false
A  Array{Float64,2}  false
B  Array{Float64,3}  false
iperm  Array{Int64,1}  false
time  Utils.Timings  false

=#
