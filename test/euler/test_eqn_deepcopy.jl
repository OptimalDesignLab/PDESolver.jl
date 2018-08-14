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

  mesh, sbp, eqn, opts = solvePDE(test_eqn_copy_inputfile)

  eqn_copy = eqn_deepcopy(mesh, sbp, eqn, opts)

  @testset "--- Testing eqn_deepcopy, DG q/q_vec & res/res_vec equality ---" begin

    @test ( mesh.isDG )== true
    @test ( (pointer(eqn.q) == pointer(eqn.q_vec)) )== true
    @test ( (pointer(eqn.res) == pointer(eqn.res_vec)) )== true
    @test ( (pointer(eqn_copy.q) == pointer(eqn_copy.q_vec)) )== true
    @test ( (pointer(eqn_copy.res) == pointer(eqn_copy.res_vec)) )== true

  end

  #----------------------
  # Test type #1: eqn_copy field values equal eqn's
  @testset "--- Testing eqn_deepcopy, copy phase, values ---" begin

    # testing 1st level

    # note: these three commented out because the double equals operator reverts back to triple equals 
    #   (test if they're same object). So the double equals will return false
    # @test ( (eqn_copy.params == eqn.params) )==
    # @test eqn_copy.params == eqn.params
    # @test ( (eqn_copy.params_conservative == eqn.params_conservative) )==
    # @test eqn_copy.params_conservative == eqn.params_conservative
    # @test ( (eqn_copy.params_entropy == eqn.params_entropy) )==
    # @test eqn_copy.params_entropy == eqn.params_entropy
    @test ( (eqn_copy.comm == eqn.comm) )== true
    @test ( (eqn_copy.commsize == eqn.commsize) )== true
    @test ( (eqn_copy.myrank == eqn.myrank) )== true
    @test ( (eqn_copy.q == eqn.q) )== true
    @test ( (eqn_copy.q_face == eqn.q_face) )== true
    @test ( (eqn_copy.q_bndry == eqn.q_bndry) )== true
    @test ( (eqn_copy.q_vec == eqn.q_vec) )== true
    @test ( (eqn_copy.aux_vars == eqn.aux_vars) )== true
    @test ( (eqn_copy.aux_vars_face == eqn.aux_vars_face) )== true
    @test ( (eqn_copy.aux_vars_sharedface == eqn.aux_vars_sharedface) )== true
    @test ( (eqn_copy.aux_vars_bndry == eqn.aux_vars_bndry) )== true
    @test ( (eqn_copy.flux_parametric == eqn.flux_parametric) )== true
    # note: q_face_send and q_face_recv purposes are now served by eqn.shared_data
    @test ( (eqn_copy.shared_data == eqn.shared_data) )== true
    @test ( (eqn_copy.flux_face == eqn.flux_face) )== true
    @test ( (eqn_copy.flux_sharedface == eqn.flux_sharedface) )== true
    @test ( (eqn_copy.res == eqn.res) )== true
    @test ( (eqn_copy.res_vec == eqn.res_vec) )== true
    @test ( (eqn_copy.Axi == eqn.Axi) )== true
    @test ( (eqn_copy.Aeta == eqn.Aeta) )== true
    @test ( (eqn_copy.res_edge == eqn.res_edge) )== true
    @test ( (eqn_copy.edgestab_alpha == eqn.edgestab_alpha) )== true
    @test ( (eqn_copy.bndryflux == eqn.bndryflux) )== true
    @test ( (eqn_copy.stabscale == eqn.stabscale) )== true
    @test ( (eqn_copy.Minv3D == eqn.Minv3D) )== true
    @test ( (eqn_copy.Minv == eqn.Minv) )== true
    @test ( (eqn_copy.M == eqn.M) )== true
    @test ( (eqn_copy.multiplyA0inv == eqn.multiplyA0inv) )== true
    @test ( (eqn_copy.majorIterationCallback == eqn.majorIterationCallback) )== true
    @test ( (eqn_copy.src_func == eqn.src_func) )== true
    @test ( (eqn_copy.flux_func == eqn.flux_func) )== true
    @test ( (eqn_copy.volume_flux_func == eqn.volume_flux_func) )== true
    @test ( (eqn_copy.face_element_integral_func == eqn.face_element_integral_func) )== true

    # testing 2nd level, params
    @test ( (eqn_copy.params.f == eqn.params.f) )== true
    @test ( (eqn_copy.params.t == eqn.params.t) )== true
    @test ( (eqn_copy.params.order == eqn.params.order) )== true
    @test ( (eqn_copy.params.q_vals == eqn.params.q_vals) )== true
    @test ( (eqn_copy.params.q_vals2 == eqn.params.q_vals2) )== true
    @test ( (eqn_copy.params.q_vals3 == eqn.params.q_vals3) )== true
    @test ( (eqn_copy.params.qg == eqn.params.qg) )== true
    @test ( (eqn_copy.params.v_vals == eqn.params.v_vals) )== true
    @test ( (eqn_copy.params.v_vals2 == eqn.params.v_vals2) )== true
    @test ( (eqn_copy.params.Lambda == eqn.params.Lambda) )== true
    @test ( (eqn_copy.params.w_vals_stencil == eqn.params.w_vals_stencil) )== true
    @test ( (eqn_copy.params.w_vals2_stencil == eqn.params.w_vals2_stencil) )== true
    @test ( (eqn_copy.params.res_vals1 == eqn.params.res_vals1) )== true
    @test ( (eqn_copy.params.res_vals2 == eqn.params.res_vals2) )== true
    @test ( (eqn_copy.params.res_vals3 == eqn.params.res_vals3) )== true
    @test ( (eqn_copy.params.flux_vals1 == eqn.params.flux_vals1) )== true
    @test ( (eqn_copy.params.flux_vals2 == eqn.params.flux_vals2) )== true
    @test ( (eqn_copy.params.sat_vals == eqn.params.sat_vals) )== true
    @test ( (eqn_copy.params.A0 == eqn.params.A0) )== true
    @test ( (eqn_copy.params.A0inv == eqn.params.A0inv) )== true
    @test ( (eqn_copy.params.A1 == eqn.params.A1) )== true
    @test ( (eqn_copy.params.A2 == eqn.params.A2) )== true
    @test ( (eqn_copy.params.S2 == eqn.params.S2) )== true
    @test ( (eqn_copy.params.A_mats == eqn.params.A_mats) )== true
    @test ( (eqn_copy.params.Rmat1 == eqn.params.Rmat1) )== true
    @test ( (eqn_copy.params.Rmat2 == eqn.params.Rmat2) )== true
    @test ( (eqn_copy.params.P == eqn.params.P) )== true
    @test ( (eqn_copy.params.nrm == eqn.params.nrm) )== true
    @test ( (eqn_copy.params.nrm2 == eqn.params.nrm2) )== true
    @test ( (eqn_copy.params.nrm3 == eqn.params.nrm3) )== true
    @test ( (eqn_copy.params.h == eqn.params.h) )== true
    @test ( (eqn_copy.params.cv == eqn.params.cv) )== true
    @test ( (eqn_copy.params.R == eqn.params.R) )== true
    @test ( (eqn_copy.params.gamma == eqn.params.gamma) )== true
    @test ( (eqn_copy.params.gamma_1 == eqn.params.gamma_1) )== true
    @test ( (eqn_copy.params.Ma == eqn.params.Ma) )== true
    @test ( (eqn_copy.params.aoa == eqn.params.aoa) )== true
    @test ( (eqn_copy.params.rho_free == eqn.params.rho_free) )== true
    @test ( (eqn_copy.params.E_free == eqn.params.E_free) )== true
    @test ( (eqn_copy.params.edgestab_gamma == eqn.params.edgestab_gamma) )== true
    @test ( (eqn_copy.params.writeflux == eqn.params.writeflux) )== true
    @test ( (eqn_copy.params.writeboundary == eqn.params.writeboundary) )== true
    @test ( (eqn_copy.params.writeq == eqn.params.writeq) )== true
    @test ( (eqn_copy.params.use_edgestab == eqn.params.use_edgestab) )== true
    @test ( (eqn_copy.params.use_filter == eqn.params.use_filter) )== true
    @test ( (eqn_copy.params.use_res_filter == eqn.params.use_res_filter) )== true
    @test ( (eqn_copy.params.filter_mat == eqn.params.filter_mat) )== true
    @test ( (eqn_copy.params.use_dissipation == eqn.params.use_dissipation) )== true
    @test ( (eqn_copy.params.dissipation_const == eqn.params.dissipation_const) )== true
    @test ( (eqn_copy.params.tau_type == eqn.params.tau_type) )== true
    @test ( (eqn_copy.params.vortex_x0 == eqn.params.vortex_x0) )== true
    @test ( (eqn_copy.params.vortex_strength == eqn.params.vortex_strength) )== true
    @test ( (eqn_copy.params.krylov_itr == eqn.params.krylov_itr) )== true
    @test ( (eqn_copy.params.krylov_type == eqn.params.krylov_type) )== true
    @test ( (eqn_copy.params.Rprime == eqn.params.Rprime) )== true
    @test ( (eqn_copy.params.A == eqn.params.A) )== true
    @test ( (eqn_copy.params.B == eqn.params.B) )== true
    @test ( (eqn_copy.params.iperm == eqn.params.iperm) )== true
    @test ( (eqn_copy.params.time == eqn.params.time) )== true

    # testing 2nd level, params_entropy
    @test ( (eqn_copy.params_entropy.f == eqn.params_entropy.f) )== true
    @test ( (eqn_copy.params_entropy.t == eqn.params_entropy.t) )== true
    @test ( (eqn_copy.params_entropy.order == eqn.params_entropy.order) )== true
    @test ( (eqn_copy.params_entropy.q_vals == eqn.params_entropy.q_vals) )== true
    @test ( (eqn_copy.params_entropy.q_vals2 == eqn.params_entropy.q_vals2) )== true
    @test ( (eqn_copy.params_entropy.q_vals3 == eqn.params_entropy.q_vals3) )== true
    @test ( (eqn_copy.params_entropy.qg == eqn.params_entropy.qg) )== true
    @test ( (eqn_copy.params_entropy.v_vals == eqn.params_entropy.v_vals) )== true
    @test ( (eqn_copy.params_entropy.v_vals2 == eqn.params_entropy.v_vals2) )== true
    @test ( (eqn_copy.params_entropy.Lambda == eqn.params_entropy.Lambda) )== true
    @test ( (eqn_copy.params_entropy.w_vals_stencil == eqn.params_entropy.w_vals_stencil) )== true
    @test ( (eqn_copy.params_entropy.w_vals2_stencil == eqn.params_entropy.w_vals2_stencil) )== true
    @test ( (eqn_copy.params_entropy.res_vals1 == eqn.params_entropy.res_vals1) )== true
    @test ( (eqn_copy.params_entropy.res_vals2 == eqn.params_entropy.res_vals2) )== true
    @test ( (eqn_copy.params_entropy.res_vals3 == eqn.params_entropy.res_vals3) )== true
    @test ( (eqn_copy.params_entropy.flux_vals1 == eqn.params_entropy.flux_vals1) )== true
    @test ( (eqn_copy.params_entropy.flux_vals2 == eqn.params_entropy.flux_vals2) )== true
    @test ( (eqn_copy.params_entropy.sat_vals == eqn.params_entropy.sat_vals) )== true
    @test ( (eqn_copy.params_entropy.A0 == eqn.params_entropy.A0) )== true
    @test ( (eqn_copy.params_entropy.A0inv == eqn.params_entropy.A0inv) )== true
    @test ( (eqn_copy.params_entropy.A1 == eqn.params_entropy.A1) )== true
    @test ( (eqn_copy.params_entropy.A2 == eqn.params_entropy.A2) )== true
    @test ( (eqn_copy.params_entropy.S2 == eqn.params_entropy.S2) )== true
    @test ( (eqn_copy.params_entropy.A_mats == eqn.params_entropy.A_mats) )== true
    @test ( (eqn_copy.params_entropy.Rmat1 == eqn.params_entropy.Rmat1) )== true
    @test ( (eqn_copy.params_entropy.Rmat2 == eqn.params_entropy.Rmat2) )== true
    @test ( (eqn_copy.params_entropy.P == eqn.params_entropy.P) )== true
    @test ( (eqn_copy.params_entropy.nrm == eqn.params_entropy.nrm) )== true
    @test ( (eqn_copy.params_entropy.nrm2 == eqn.params_entropy.nrm2) )== true
    @test ( (eqn_copy.params_entropy.nrm3 == eqn.params_entropy.nrm3) )== true
    @test ( (eqn_copy.params_entropy.h == eqn.params_entropy.h) )== true
    @test ( (eqn_copy.params_entropy.cv == eqn.params_entropy.cv) )== true
    @test ( (eqn_copy.params_entropy.R == eqn.params_entropy.R) )== true
    @test ( (eqn_copy.params_entropy.gamma == eqn.params_entropy.gamma) )== true
    @test ( (eqn_copy.params_entropy.gamma_1 == eqn.params_entropy.gamma_1) )== true
    @test ( (eqn_copy.params_entropy.Ma == eqn.params_entropy.Ma) )== true
    @test ( (eqn_copy.params_entropy.aoa == eqn.params_entropy.aoa) )== true
    @test ( (eqn_copy.params_entropy.rho_free == eqn.params_entropy.rho_free) )== true
    @test ( (eqn_copy.params_entropy.E_free == eqn.params_entropy.E_free) )== true
    @test ( (eqn_copy.params_entropy.edgestab_gamma == eqn.params_entropy.edgestab_gamma) )== true
    @test ( (eqn_copy.params_entropy.writeflux == eqn.params_entropy.writeflux) )== true
    @test ( (eqn_copy.params_entropy.writeboundary == eqn.params_entropy.writeboundary) )== true
    @test ( (eqn_copy.params_entropy.writeq == eqn.params_entropy.writeq) )== true
    @test ( (eqn_copy.params_entropy.use_edgestab == eqn.params_entropy.use_edgestab) )== true
    @test ( (eqn_copy.params_entropy.use_filter == eqn.params_entropy.use_filter) )== true
    @test ( (eqn_copy.params_entropy.use_res_filter == eqn.params_entropy.use_res_filter) )== true
    @test ( (eqn_copy.params_entropy.filter_mat == eqn.params_entropy.filter_mat) )== true
    @test ( (eqn_copy.params_entropy.use_dissipation == eqn.params_entropy.use_dissipation) )== true
    @test ( (eqn_copy.params_entropy.dissipation_const == eqn.params_entropy.dissipation_const) )== true
    @test ( (eqn_copy.params_entropy.tau_type == eqn.params_entropy.tau_type) )== true
    @test ( (eqn_copy.params_entropy.vortex_x0 == eqn.params_entropy.vortex_x0) )== true
    @test ( (eqn_copy.params_entropy.vortex_strength == eqn.params_entropy.vortex_strength) )== true
    @test ( (eqn_copy.params_entropy.krylov_itr == eqn.params_entropy.krylov_itr) )== true
    @test ( (eqn_copy.params_entropy.krylov_type == eqn.params_entropy.krylov_type) )== true
    @test ( (eqn_copy.params_entropy.Rprime == eqn.params_entropy.Rprime) )== true
    @test ( (eqn_copy.params_entropy.A == eqn.params_entropy.A) )== true
    @test ( (eqn_copy.params_entropy.B == eqn.params_entropy.B) )== true
    @test ( (eqn_copy.params_entropy.iperm == eqn.params_entropy.iperm) )== true
    @test ( (eqn_copy.params_entropy.time == eqn.params_entropy.time) )== true

  end   # end of facts check on values of eqn_copy matching eqn

  @testset "--- Testing eqn_deepcopy, copy phase, pointers of arrays ---" begin

    # 1st level fields of type Array
    @test ( (pointer(eqn_copy.q) == pointer(eqn.q)) )== false
    @test ( (pointer(eqn_copy.q_face) == pointer(eqn.q_face)) )== false
    @test ( (pointer(eqn_copy.q_bndry) == pointer(eqn.q_bndry)) )== false
    @test ( (pointer(eqn_copy.q_vec) == pointer(eqn.q_vec)) )== false
    @test ( (pointer(eqn_copy.aux_vars) == pointer(eqn.aux_vars)) )== false
    @test ( (pointer(eqn_copy.aux_vars_face) == pointer(eqn.aux_vars_face)) )== false
    @test ( (pointer(eqn_copy.aux_vars_sharedface) == pointer(eqn.aux_vars_sharedface)) )== false
    @test ( (pointer(eqn_copy.aux_vars_bndry) == pointer(eqn.aux_vars_bndry)) )== false
    @test ( (pointer(eqn_copy.flux_parametric) == pointer(eqn.flux_parametric)) )== false
    # note: q_face_send and q_face_recv purposes are now served by eqn.shared_data
    @test ( (pointer(eqn_copy.shared_data) == pointer(eqn.shared_data)) )== false
    @test ( (pointer(eqn_copy.flux_face) == pointer(eqn.flux_face)) )== false
    @test ( (pointer(eqn_copy.flux_sharedface) == pointer(eqn.flux_sharedface)) )== false
    @test ( (pointer(eqn_copy.res) == pointer(eqn.res)) )== false
    @test ( (pointer(eqn_copy.res_vec) == pointer(eqn.res_vec)) )== false
    @test ( (pointer(eqn_copy.Axi) == pointer(eqn.Axi)) )== false
    @test ( (pointer(eqn_copy.Aeta) == pointer(eqn.Aeta)) )== false
    @test ( (pointer(eqn_copy.res_edge) == pointer(eqn.res_edge)) )== false
    @test ( (pointer(eqn_copy.edgestab_alpha) == pointer(eqn.edgestab_alpha)) )== false
    @test ( (pointer(eqn_copy.bndryflux) == pointer(eqn.bndryflux)) )== false
    @test ( (pointer(eqn_copy.stabscale) == pointer(eqn.stabscale)) )== false
    @test ( (pointer(eqn_copy.Minv3D) == pointer(eqn.Minv3D)) )== false
    @test ( (pointer(eqn_copy.Minv) == pointer(eqn.Minv)) )== false
    @test ( (pointer(eqn_copy.M) == pointer(eqn.M)) )== false

    # 2nd level fields of type Array (just going to check params)
    @test ( (pointer(eqn_copy.params.q_vals) == pointer(eqn.params.q_vals)) )== false
    @test ( (pointer(eqn_copy.params.q_vals2) == pointer(eqn.params.q_vals2)) )== false
    @test ( (pointer(eqn_copy.params.q_vals3) == pointer(eqn.params.q_vals3)) )== false
    @test ( (pointer(eqn_copy.params.qg) == pointer(eqn.params.qg)) )== false
    @test ( (pointer(eqn_copy.params.v_vals) == pointer(eqn.params.v_vals)) )== false
    @test ( (pointer(eqn_copy.params.v_vals2) == pointer(eqn.params.v_vals2)) )== false
    @test ( (pointer(eqn_copy.params.Lambda) == pointer(eqn.params.Lambda)) )== false
    @test ( (pointer(eqn_copy.params.w_vals_stencil) == pointer(eqn.params.w_vals_stencil)) )== false
    @test ( (pointer(eqn_copy.params.w_vals2_stencil) == pointer(eqn.params.w_vals2_stencil)) )== false
    @test ( (pointer(eqn_copy.params.res_vals1) == pointer(eqn.params.res_vals1)) )== false
    @test ( (pointer(eqn_copy.params.res_vals2) == pointer(eqn.params.res_vals2)) )== false
    @test ( (pointer(eqn_copy.params.res_vals3) == pointer(eqn.params.res_vals3)) )== false
    @test ( (pointer(eqn_copy.params.flux_vals1) == pointer(eqn.params.flux_vals1)) )== false
    @test ( (pointer(eqn_copy.params.flux_vals2) == pointer(eqn.params.flux_vals2)) )== false
    @test ( (pointer(eqn_copy.params.sat_vals) == pointer(eqn.params.sat_vals)) )== false
    @test ( (pointer(eqn_copy.params.A0) == pointer(eqn.params.A0)) )== false
    @test ( (pointer(eqn_copy.params.A0inv) == pointer(eqn.params.A0inv)) )== false
    @test ( (pointer(eqn_copy.params.A1) == pointer(eqn.params.A1)) )== false
    @test ( (pointer(eqn_copy.params.A2) == pointer(eqn.params.A2)) )== false
    @test ( (pointer(eqn_copy.params.S2) == pointer(eqn.params.S2)) )== false
    @test ( (pointer(eqn_copy.params.A_mats) == pointer(eqn.params.A_mats)) )== false
    @test ( (pointer(eqn_copy.params.Rmat1) == pointer(eqn.params.Rmat1)) )== false
    @test ( (pointer(eqn_copy.params.Rmat2) == pointer(eqn.params.Rmat2)) )== false
    @test ( (pointer(eqn_copy.params.P) == pointer(eqn.params.P)) )== false
    @test ( (pointer(eqn_copy.params.nrm) == pointer(eqn.params.nrm)) )== false
    @test ( (pointer(eqn_copy.params.nrm2) == pointer(eqn.params.nrm2)) )== false
    @test ( (pointer(eqn_copy.params.nrm3) == pointer(eqn.params.nrm3)) )== false
    @test ( (pointer(eqn_copy.params.filter_mat) == pointer(eqn.params.filter_mat)) )== false
    @test ( (pointer(eqn_copy.params.Rprime) == pointer(eqn.params.Rprime)) )== false
    @test ( (pointer(eqn_copy.params.A) == pointer(eqn.params.A)) )== false
    @test ( (pointer(eqn_copy.params.B) == pointer(eqn.params.B)) )== false
    @test ( (pointer(eqn_copy.params.iperm) == pointer(eqn.params.iperm)) )== false

    @test ( (pointer(eqn_copy.params_entropy.q_vals) == pointer(eqn.params_entropy.q_vals)) )== false
    @test ( (pointer(eqn_copy.params_entropy.q_vals2) == pointer(eqn.params_entropy.q_vals2)) )== false
    @test ( (pointer(eqn_copy.params_entropy.q_vals3) == pointer(eqn.params_entropy.q_vals3)) )== false
    @test ( (pointer(eqn_copy.params_entropy.qg) == pointer(eqn.params_entropy.qg)) )== false
    @test ( (pointer(eqn_copy.params_entropy.v_vals) == pointer(eqn.params_entropy.v_vals)) )== false
    @test ( (pointer(eqn_copy.params_entropy.v_vals2) == pointer(eqn.params_entropy.v_vals2)) )== false
    @test ( (pointer(eqn_copy.params_entropy.Lambda) == pointer(eqn.params_entropy.Lambda)) )== false
    @test ( (pointer(eqn_copy.params_entropy.w_vals_stencil) == pointer(eqn.params_entropy.w_vals_stencil)) )== false
    @test ( (pointer(eqn_copy.params_entropy.w_vals2_stencil) == pointer(eqn.params_entropy.w_vals2_stencil)) )== false
    @test ( (pointer(eqn_copy.params_entropy.res_vals1) == pointer(eqn.params_entropy.res_vals1)) )== false
    @test ( (pointer(eqn_copy.params_entropy.res_vals2) == pointer(eqn.params_entropy.res_vals2)) )== false
    @test ( (pointer(eqn_copy.params_entropy.res_vals3) == pointer(eqn.params_entropy.res_vals3)) )== false
    @test ( (pointer(eqn_copy.params_entropy.flux_vals1) == pointer(eqn.params_entropy.flux_vals1)) )== false
    @test ( (pointer(eqn_copy.params_entropy.flux_vals2) == pointer(eqn.params_entropy.flux_vals2)) )== false
    @test ( (pointer(eqn_copy.params_entropy.sat_vals) == pointer(eqn.params_entropy.sat_vals)) )== false
    @test ( (pointer(eqn_copy.params_entropy.A0) == pointer(eqn.params_entropy.A0)) )== false
    @test ( (pointer(eqn_copy.params_entropy.A0inv) == pointer(eqn.params_entropy.A0inv)) )== false
    @test ( (pointer(eqn_copy.params_entropy.A1) == pointer(eqn.params_entropy.A1)) )== false
    @test ( (pointer(eqn_copy.params_entropy.A2) == pointer(eqn.params_entropy.A2)) )== false
    @test ( (pointer(eqn_copy.params_entropy.S2) == pointer(eqn.params_entropy.S2)) )== false
    @test ( (pointer(eqn_copy.params_entropy.A_mats) == pointer(eqn.params_entropy.A_mats)) )== false
    @test ( (pointer(eqn_copy.params_entropy.Rmat1) == pointer(eqn.params_entropy.Rmat1)) )== false
    @test ( (pointer(eqn_copy.params_entropy.Rmat2) == pointer(eqn.params_entropy.Rmat2)) )== false
    @test ( (pointer(eqn_copy.params_entropy.P) == pointer(eqn.params_entropy.P)) )== false
    @test ( (pointer(eqn_copy.params_entropy.nrm) == pointer(eqn.params_entropy.nrm)) )== false
    @test ( (pointer(eqn_copy.params_entropy.nrm2) == pointer(eqn.params_entropy.nrm2)) )== false
    @test ( (pointer(eqn_copy.params_entropy.nrm3) == pointer(eqn.params_entropy.nrm3)) )== false
    @test ( (pointer(eqn_copy.params_entropy.filter_mat) == pointer(eqn.params_entropy.filter_mat)) )== false
    @test ( (pointer(eqn_copy.params_entropy.Rprime) == pointer(eqn.params_entropy.Rprime)) )== false
    @test ( (pointer(eqn_copy.params_entropy.A) == pointer(eqn.params_entropy.A)) )== false
    @test ( (pointer(eqn_copy.params_entropy.B) == pointer(eqn.params_entropy.B)) )== false
    @test ( (pointer(eqn_copy.params_entropy.iperm) == pointer(eqn.params_entropy.iperm)) )== false

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

  @testset "--- Testing eqn_deepcopy, eqn change phase ---" begin

    @test ( (eqn_copy.q_face[1] == eqn.q_face[1]) )== false
    @test ( (eqn_copy.flux_parametric[1] == eqn.flux_parametric[1]) )== false
    @test ( (eqn_copy.flux_face[1] == eqn.flux_face[1]) )== false
    @test ( (eqn_copy.q[1] == eqn.q[1]) )== false
    @test ( (eqn_copy.res[1] == eqn.res[1]) )== false
    @test ( (eqn_copy.M[1] == eqn.M[1]) )== false
    @test ( (eqn_copy.Minv[1] == eqn.Minv[1]) )== false
    @test ( (eqn_copy.Minv3D[1] == eqn.Minv3D[1]) )== false
    @test ( (eqn_copy.bndryflux[1] == eqn.bndryflux[1]) )== false
    @test ( (eqn_copy.q_bndry[1] == eqn.q_bndry[1]) )== false

  end

  rand!(eqn_copy.q)
  rand!(eqn_copy.res)

  @testset "--- Testing eqn_deepcopy, eqn_copy change phase ---" begin

    @test ( (eqn_copy.q[1] == eqn_copy.q_vec[1]) )== true
    @test ( (eqn_copy.res[1] == eqn_copy.res_vec[1]) )== true

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
