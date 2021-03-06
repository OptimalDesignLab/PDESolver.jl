# test jacobian calculation functions

"""
  Test the jacobian of individual terms.  This only tests the SBP Omega
  operators due to test time limits.
"""
function test_jac_terms()

  fname = "input_vals_jac2d.jl"
  fname3 = "input_vals_jac3d.jl"
  mesh, sbp, eqn, opts = run_solver(fname)
  mesh3, sbp3, eqn3, opts3 = run_solver(fname3)

#=
  # SBPOmega, Petsc Mat
  fname4 = "input_vals_jac_tmp.jl"
  opts_tmp = read_input_file(fname3)
  opts_tmp["jac_type"] = 3
  opts_tmp["operator_type"] = "SBPOmega"
  make_input(opts_tmp, fname4)
  mesh4, sbp4, eqn4, opts4 = run_solver(fname4)
=#

  # list of boundary conditions to test revm method in 2D
  bclist_revm_2d = ["noPenetrationBC", "FreeStreamBC", "ExpBC",
                    "isentropicVortexBC", "noPenetrationESBC", "ZeroFluxBC",
                    "wedge20BC"]
  bclist_revm_3d = [                   "FreeStreamBC"]
  bclist_revq_2d = ["noPenetrationBC", "FreeStreamBC", "ExpBC",
                    "isentropicVortexBC", "noPenetrationESBC", "ZeroFluxBC",
                    "wedge20BC"]
  bclist_revq_3d = [                      "FreeStreamBC",]
  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)
  Tmsh = eltype(mesh.jac)


  @testset "----- Testing jacobian -----" begin

    func = EulerEquationMod.FluxDict["HLLFlux"]
    func_diff = EulerEquationMod.FluxDict_diff["HLLFlux"]
    func_revm = EulerEquationMod.FluxDict_revm["HLLFlux"]
    func_revq = EulerEquationMod.FluxDict_revq["HLLFlux"]

    nrm = [0.45, 0.55]
    nrm2 = -nrm
#=
    #q = Complex128[2.0, 3.0, 4.0, 7.0]
    #qg = q + 0.1
    q = Complex128[1.1, 0.64405, 0.73405, 2.2] 
    qg = q + 0.1
    q = Complex128[2.0, 3.0, 4.0, 7.0]
    qg = q + 0.1


    test_2RWave(eqn.params, q, qg, nrm)
    test_2RWave_revq(eqn.params, q, qg, nrm)
    test_2RWave_revm(eqn.params, q, qg, nrm)
    #test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
    #test_2flux_revm(eqn.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn.params, q, qg, nrm, func, func_revq)
=#


    
#    q = Complex128[2.0, 3.0, 4.0, 7.0]
#    qg = q + 0.1
#    q = Complex128[1.1, -0.72405, -0.82405, 2.2] 
#    qg = q + 0.1
    q = Complex128[1.1, -0.681, 0.47, 2.2]
    qg = q + 0.1

    test_ad_inner(eqn.params, q, qg, nrm2, func, func_diff)

#    error("stop here")

    test_pressure(eqn.params)
    test_pressure(eqn3.params)

    test_eulerflux(eqn.params)
    test_eulerflux_revm(eqn.params)
    test_eulerflux_revq(eqn.params)
    test_eulerflux(eqn3.params)
    test_eulerflux_revm(eqn3.params)
    test_eulerflux_revq(eqn3.params)

    test_logavg()

    test_IRA0(eqn.params)
    test_IRA0(eqn3.params)

    test_SCurvilinear(eqn.params, sbp)

    nrm = [0.45, 0.55]
    nrm2 = -nrm
    coords = Complex128[1.1, 1.2]


    println("testing all positive eigenvalues")

    # 2 point flux functions
    func = EulerEquationMod.FluxDict["RoeFlux"]
    func_diff = EulerEquationMod.FluxDict_diff["RoeFlux"]
    func_revm = EulerEquationMod.FluxDict_revm["RoeFlux"]
    func_revq = EulerEquationMod.FluxDict_revq["RoeFlux"]
   
    func2 = EulerEquationMod.FluxDict["LFFlux"]
    func2_diff = EulerEquationMod.FluxDict_diff["LFFlux"]
    func2_revm = EulerEquationMod.FluxDict_revm["LFFlux"]
    func2_revq = EulerEquationMod.FluxDict_revq["LFFlux"]
 
    func3 = EulerEquationMod.FluxDict["IRFlux"]
    func3_diff = EulerEquationMod.FluxDict_diff["IRFlux"]
    func3_revm = EulerEquationMod.FluxDict_revm["IRFlux"]
    func3_revq = EulerEquationMod.FluxDict_revq["IRFlux"]
    
    func4 = EulerEquationMod.FluxDict["IRSLFFlux"]
    func4_diff = EulerEquationMod.FluxDict_diff["IRSLFFlux"]
    func4_revm = EulerEquationMod.FluxDict_revm["IRSLFFlux"]
    func4_revq = EulerEquationMod.FluxDict_revq["IRSLFFlux"]
 
    func5 = EulerEquationMod.FluxDict["HLLFlux"]
    func5_diff = EulerEquationMod.FluxDict_diff["HLLFlux"]
    func5_revm = EulerEquationMod.FluxDict_revm["HLLFlux"]
    func5_revq = EulerEquationMod.FluxDict_revq["HLLFlux"]


    # Abstract Entropy Kernels
    lf_kernel = EulerEquationMod.LFKernel{Tsol, Tres, Tmsh}(mesh.numDofPerNode, 2*mesh.numDofPerNode)
    lf_kernel3 = EulerEquationMod.LFKernel{Tsol, Tres, Tmsh}(mesh3.numDofPerNode, 2*mesh3.numDofPerNode)
    
    q = Complex128[2.0, 3.0, 4.0, 7.0]
    qg = q + 0.1

    test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
    test_2flux_revm(eqn.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn.params, q, qg, nrm, func, func_revq)

    test_lambda(eqn.params, q, nrm)
    test_lambdasimple(eqn.params, q, qg, nrm)
    test_ad_inner(eqn.params, q, qg, nrm, func2, func2_diff)
    test_2flux_revm(eqn.params, q, qg, nrm, func2, func2_revm)
    test_2flux_revq(eqn.params, q, qg, nrm, func2, func2_revq)


    # make sure arrays are zerod out
    test_ad_inner(eqn.params, q, qg, nrm, func2, func2_diff)
    test_ad_inner(eqn.params, q, qg, nrm, func3, func3_diff)
    test_2flux_revq(eqn.params, q, qg, nrm, func3, func3_revq, test_multid=true)
    test_2flux_revm(eqn.params, q, qg, nrm, func3, func3_revm, test_multid=true)
   
    println("testing applyEntropyKernel_diagE")
    test_ad_inner(eqn.params, q, qg, nrm, func4, func4_diff)
    test_2flux_revq(eqn.params, q, qg, nrm, func4, func4_revq)
    test_2flux_revm(eqn.params, q, qg, nrm, func4, func4_revm)
 
    test_EntropyKernel(eqn.params, lf_kernel)
    test_EntropyKernel_revq(eqn.params, lf_kernel)
    test_EntropyKernel_revm(eqn.params, lf_kernel)

    test_ad_inner(eqn.params, q, qg, nrm, func5, func5_diff)
    test_2flux_revq(eqn.params, q, qg, nrm, func5, func5_revq)
    test_2flux_revm(eqn.params, q, qg, nrm, func5, func5_revm)
 

    # test boundary conditions
    for bcname in bclist_revm_2d
      bc = EulerEquationMod.BCDict[bcname](mesh, eqn)
      bc_revm = EulerEquationMod.BCDict_revm[bcname](mesh, eqn)
      test_bc_revm(eqn.params, bc, bc_revm)
    end

    for bcname in bclist_revq_2d
      bc = EulerEquationMod.BCDict[bcname](mesh, eqn)
      bc_revq = EulerEquationMod.BCDict_revq[bcname](mesh, eqn)
      test_bc_revq(eqn.params, bc, bc_revq)
    end

    # test boundary conditions
    for bcname in bclist_revm_3d
      bc = EulerEquationMod.BCDict[bcname](mesh3, eqn3)
      bc_revm = EulerEquationMod.BCDict_revm[bcname](mesh3, eqn3)
      test_bc_revm(eqn3.params, bc, bc_revm)
    end

    for bcname in bclist_revq_3d
      bc = EulerEquationMod.BCDict[bcname](mesh3, eqn3)
      bc_revq = EulerEquationMod.BCDict_revq[bcname](mesh3, eqn3)
      test_bc_revq(eqn3.params, bc, bc_revq)
    end


    # test common_funcs
    test_common_func_rev(eqn.params, coords, EulerEquationMod.calcExp,
                                             EulerEquationMod.calcExp_rev)
 
    test_common_func_rev(eqn.params, coords, EulerEquationMod.calcIsentropicVortex,
                                              EulerEquationMod.calcIsentropicVortex_rev)
  
    println("testing all negative eigenvalues")
    q = Complex128[2.0, 3.0, 4.0, 7.0]
    qg = q + 0.1
    test_ad_inner(eqn.params, q, qg, nrm2, func, func_diff)
    test_2flux_revm(eqn.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn.params, q, qg, nrm, func, func_revq)


    test_lambda(eqn.params, q, nrm2)
    test_lambdasimple(eqn.params, q, qg, nrm2)
    
    test_2RWave(eqn.params, q, qg, nrm)
    test_2RWave_revq(eqn.params, q, qg, nrm)
    test_2RWave_revm(eqn.params, q, qg, nrm)


    test_ad_inner(eqn.params, q, qg, nrm2, func5, func5_diff)
    test_2flux_revq(eqn.params, q, qg, nrm, func5, func5_revq)
    test_2flux_revm(eqn.params, q, qg, nrm, func5, func5_revm)
 


    println("testing lambda1 entropy fix")
    q = Complex128[1.1, -0.72405, -0.82405, 2.2] 
    qg = q + 0.1
    test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
    test_2flux_revm(eqn.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn.params, q, qg, nrm, func, func_revq)
   
    test_2RWave(eqn.params, q, qg, nrm)
    test_2RWave_revq(eqn.params, q, qg, nrm)
    test_2RWave_revm(eqn.params, q, qg, nrm)

    test_ad_inner(eqn.params, q, qg, nrm, func5, func5_diff)
    test_2flux_revq(eqn.params, q, qg, nrm, func5, func5_revq)
    test_2flux_revm(eqn.params, q, qg, nrm, func5, func5_revm)
 
    println("testing lambda2 entropy fix")
    q = Complex128[1.1, 0.64405, 0.73405, 2.2] 
    qg = q + 0.1
    test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn.params, q, qg, nrm2, func, func_diff)
    test_2flux_revm(eqn.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn.params, q, qg, nrm, func, func_revq)

    test_ad_inner(eqn.params, q, qg, nrm, func5, func5_diff)
    test_ad_inner(eqn.params, q, qg, nrm2, func5, func5_diff)
    test_2flux_revq(eqn.params, q, qg, nrm, func5, func5_revq)
    test_2flux_revm(eqn.params, q, qg, nrm, func5, func5_revm)
 

    println("testing lambda3 entropy fix")
    q = Complex128[1.1, -0.681, 0.47, 2.2]
    qg = q + 0.1
    test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn.params, q, qg, nrm2, func, func_diff)
    test_2flux_revm(eqn.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn.params, q, qg, nrm, func, func_revq)

    test_ad_inner(eqn.params, q, qg, nrm, func5, func5_diff)
    test_ad_inner(eqn.params, q, qg, nrm2, func5, func5_diff)
    test_2flux_revq(eqn.params, q, qg, nrm, func5, func5_revq)
    test_2flux_revm(eqn.params, q, qg, nrm, func5, func5_revm)
 

    test_faceElementIntegral(eqn.params, mesh.sbpface, func3, func3_diff)
    
    test_entropyPenalty(eqn.params, mesh.sbpface, lf_kernel)

    #--------------------------------------------------------------------------
    # 3D
    nrm = [0.45, 0.55, 0.65]
    nrm2 = -nrm

    println("testing all positive eigenvalues")
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm, func, func_diff)
    test_2flux_revm(eqn3.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn3.params, q, qg, nrm, func, func_revq)
    test_lambda(eqn3.params, q, nrm)
    test_lambdasimple(eqn3.params, q, qg, nrm)

    test_EntropyKernel(eqn3.params, lf_kernel3)
    test_EntropyKernel_revq(eqn3.params, lf_kernel3)
    test_EntropyKernel_revm(eqn3.params, lf_kernel3)

    # test with nd > required nd
    lf_kernel3 = EulerEquationMod.LFKernel{Tsol, Tmsh, Tres}(mesh3.numDofPerNode, 2*mesh.numDofPerNode + 2)
    test_EntropyKernel(eqn3.params, lf_kernel3)

    println("testing all negative eigenvalues")
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm2, func, func_diff)
    test_2flux_revm(eqn3.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn3.params, q, qg, nrm, func, func_revq)
    test_lambda(eqn3.params, q, nrm2)
    test_lambdasimple(eqn3.params, q, qg, nrm2)


    println("testing lambda1 entropy fix")
    # lambda1 entropy fix active
    q = Complex128[1.05, 0.9, 1.2, -1.3, 2.5]
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn3.params, q, qg, nrm2, func, func_diff)
    test_2flux_revm(eqn3.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn3.params, q, qg, nrm, func, func_revq)

    println("testing lambda2 entropy fix")
    # lambda1 entropy fix active
    q = Complex128[1.05, -1.1, -1.2, -1.3, 2.5] 
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn3.params, q, qg, nrm2, func, func_diff)
    test_2flux_revm(eqn3.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn3.params, q, qg, nrm, func, func_revq)

    println("testing lambda3 entropy fix")
    # lambda3 entropy fix active
    q = Complex128[1.05, -0.52, -0.47, -0.36, 8.5]
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn3.params, q, qg, nrm2, func, func_diff)
    test_2flux_revm(eqn3.params, q, qg, nrm, func, func_revm)
    test_2flux_revq(eqn3.params, q, qg, nrm, func, func_revq)

    test_ad_inner(eqn3.params, q, qg, nrm, func3, func3_diff)
    test_2flux_revq(eqn3.params, q, qg, nrm, func3, func3_revq, test_multid=true)
    test_2flux_revm(eqn3.params, q, qg, nrm, func3, func3_revm, test_multid=true)


#    test_faceflux_revm(mesh, sbp, eqn, opts)


    println("\ntesting jac assembly 2d")
    test_jac_assembly(mesh, sbp, eqn, opts)

    opts_tmp = copy(opts)
    test_jac_homotopy(mesh, sbp, eqn, opts_tmp)


    println("\ntesting jac assembly 3d")
    test_jac_assembly(mesh3, sbp3, eqn3, opts3)


   opts_sc = read_input_file(fname)
   opts_sc["addShockCapturing"] = true
   opts_sc["homotopy_shock_capturing"] = true
   opts_sc["homotopy_shock_sensor"] = "SensorBO"
   opts_sc["shock_sensor_name"] = "SensorHHO"
   mesh, sbp, eqn, opts = run_solver(opts_sc)

   test_jac_homotopy(mesh, sbp, eqn, opts_sc)

  end

  return nothing
end


add_func1!(EulerTests, test_jac_terms, [TAG_SHORTTEST, TAG_JAC, TAG_TMP])


"""
  Tests assembling the jacobian of all the different operators
"""
function test_jac_terms_long()

  @testset "----- Testing additional Jacobian calculation -----" begin

    # starting point for different configurations
    fname = "input_vals_jac2d.jl"
    fname3 = "input_vals_jac3d.jl"

    mesh3, sbp3, eqn3, opts3 = solvePDE("input_vals_jaccurve3d.jl")
    test_sbp_cartesian_revm(mesh3, sbp3, eqn3, opts3)
#=
    #TESTING
    # SBPDiagonalE, SparseMatrixCSC, SBPParabolicReducedSC
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 2  # was 2
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["order"] = 1
    opts_tmp["use_lps"] = true
    opts_tmp["addShockCapturing"] = true
    opts_tmp["shock_capturing_name"] = "SBPParabolicReduced"
    opts_tmp["shock_sensor_name"] = "SensorOddBO"
    opts_tmp["need_adjoint"] = true

    #opts_tmp["use_lps"] = true
    make_input(opts_tmp, fname4)
    mesh9, sbp9, eqn9, opts9 = run_solver(fname4)

    testQx(mesh9, sbp9, eqn9, opts9)
=#


    # SBPGamma, Petsc Mat
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    make_input(opts_tmp, fname4)
    mesh4, sbp4, eqn4, opts4 = run_solver(fname4)


    # SBPDiagonalE, Petsc Mat
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["order"] = 2
#    opts_tmp["write_dofs"] = true
    make_input(opts_tmp, fname4)
    mesh5, sbp5, eqn5, opts5 = run_solver(fname4)


    # SBPDiagonalE, SparseMatrixCSC
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 2
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["use_lps"] = true
    make_input(opts_tmp, fname4)
    mesh6, sbp6, eqn6, opts6 = run_solver(fname4)

    # SBPDiagonalE, Petsc Mat, use_Minv
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["order"] = 2
    opts_tmp["use_Minv"] = true
    make_input(opts_tmp, fname4)
    mesh7, sbp7, eqn7, opts7 = run_solver(fname4)

    # SBPOmega, Petsc Mat
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPOmega"
    make_input(opts_tmp, fname4)
    mesh8, sbp8, eqn8, opts8 = run_solver(fname4)

    # SBPOmega, SparseMatrixCSC
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 2
    opts_tmp["operator_type"] = "SBPOmega"
    make_input(opts_tmp, fname4)
    mesh9, sbp9, eqn9, opts9 = run_solver(fname4)

    # SBPOmega, SparseMatrixCSC, ES
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 2
    opts_tmp["operator_type"] = "SBPOmega"
    opts_tmp["volume_integral_type"] = 2
    opts_tmp["Volume_flux_name"] = "IRFlux"
    opts_tmp["Flux_name"] = "IRFlux"
    opts_tmp["FaceElementIntegral_name"] = "ESLFFaceIntegral"
    make_input(opts_tmp, fname4)
    mesh10, sbp10, eqn10, opts10 = run_solver(fname4)

    # SBPOmega, Petsc Mat, ES
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPOmega"
    opts_tmp["volume_integral_type"] = 2
    opts_tmp["Volume_flux_name"] = "IRFlux"
    opts_tmp["face_integral_type"] = 2
    opts_tmp["Flux_name"] = "IRFlux"
    opts_tmp["FaceElementIntegral_name"] = "ESLFFaceIntegral"
    make_input(opts_tmp, fname4)
    mesh11, sbp11, eqn11, opts11 = run_solver(fname4)

    # SBPDiagonalE, Petsc Mat, ES
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["volume_integral_type"] = 2
    opts_tmp["Volume_flux_name"] = "IRFlux"
    opts_tmp["face_integral_type"] = 1
    opts_tmp["Flux_name"] = "IRSLFFlux"
    make_input(opts_tmp, fname4)
    mesh12, sbp12, eqn12, opts12 = run_solver(fname4)


    # SBPDiagonalE, SparseMatrixCSC, SBPParabolicReducedSC
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 2  # was 2
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["order"] = 1
    opts_tmp["use_lps"] = true
    opts_tmp["addShockCapturing"] = true
    #opts_tmp["shock_capturing_name"] = "SBPParabolicReduced"
    #opts_tmp["shock_sensor_name"] = "SensorOddBO"
    opts_tmp["shock_capturing_name"] = "SBPParabolicReduced"
    opts_tmp["shock_sensor_name"] = "SensorBO"

    opts_tmp["need_adjoint"] = true

    #opts_tmp["use_lps"] = true
    make_input(opts_tmp, fname4)
    mesh13, sbp13, eqn13, opts13 = run_solver(fname4)
    
    test_jac_general(mesh13, sbp13, eqn13, opts13)



    # test various matrix and operator combinations
    println("testing mode 4")
    test_jac_general(mesh4, sbp4, eqn4, opts4)
    test_BDiagPC(mesh4, sbp4, eqn4, opts4)
    test_BJacobiPC(mesh4, sbp4, eqn4, opts4)

    opts4["preallocate_jacobian_coloring"] = true
    test_jac_general(mesh4, sbp4, eqn4, opts4, is_prealloc_exact=false, set_prealloc=false)
    
    println("testing mode 5")
    test_jac_general(mesh5, sbp5, eqn5, opts5)
    # run the test twice to make sure the arrays are zeroed out properly
    println("testing mode 5 twice")
    test_jac_general(mesh5, sbp5, eqn5, opts5)

    println("testing mode 6")
    test_jac_general(mesh6, sbp6, eqn6, opts6)
    test_diagjac(mesh6, sbp6, eqn6, opts6)
 
    println("testing mode 7")
    test_jac_general(mesh7, sbp7, eqn7, opts7)
 
    println("testing mode 8")
    test_jac_general(mesh8, sbp8, eqn8, opts8)
  
    opts4["preallocate_jacobian_coloring"] = true
    test_jac_general(mesh8, sbp8, eqn8, opts8, is_prealloc_exact=true, set_prealloc=false)

    println("testing mode 9")
    test_jac_general(mesh9, sbp9, eqn9, opts9)
    test_diagjac(mesh9, sbp9, eqn9, opts9)
    test_strongdiagjac(mesh9, sbp9, eqn9, opts9)

    println("testing mode 10")
    test_jac_general(mesh10, sbp10, eqn10, opts10)
 
    println("testing mode 11")
    test_jac_general(mesh11, sbp11, eqn11, opts11)
 
    println("testing mode 12")
    test_jac_general(mesh12, sbp12, eqn12, opts12)


    println("testing mode 13")
    test_revm_product2(mesh13, sbp13, eqn13, opts13)
    test_revq_product(mesh13, sbp13, eqn13, opts13)
    test_jac_general(mesh13, sbp13, eqn13, opts13)


    # test shock capturing schemes
    println("\ntesting shock capturing")
    shock_capturing_schemes = ["SBPParabolic"]
    for sc in shock_capturing_schemes

      Tsol = eltype(eqn4.q); Tres = eltype(eqn4.res)

      println("case4")
      sensor = EulerEquationMod.ShockSensorVelocity{Tsol, Tres}(mesh4, sbp4, opts4)
      capture = EulerEquationMod.ShockCapturingDict[sc]{Tsol, Tres}(mesh4,
                                                    sbp4, eqn4, opts4, sensor)
      test_sbp_cartesian(mesh4, sbp4, eqn4, opts4)
      test_shock_capturing_jac(mesh4, sbp4, eqn4, opts4, capture)
      test_shock_capturing_jac(mesh4, sbp4, eqn4, opts4, capture, freeze=true)


      println("case5")
      sensor = EulerEquationMod.ShockSensorVelocity{Tsol, Tres}(mesh5, sbp5, opts5)
      capture = EulerEquationMod.ShockCapturingDict[sc]{Tsol, Tres}(mesh5,
                                                    sbp5, eqn5, opts5, sensor) 
      test_sbp_cartesian(mesh5, sbp5, eqn5, opts5)
      test_shock_capturing_jac(mesh5, sbp5, eqn5, opts5, capture)

      println("case6")
      sensor = EulerEquationMod.ShockSensorVelocity{Tsol, Tres}(mesh6, sbp6, opts6)
      capture = EulerEquationMod.ShockCapturingDict[sc]{Tsol, Tres}(mesh6,
                                                    sbp6, eqn6, opts6, sensor) 
      test_sbp_cartesian(mesh6, sbp6, eqn6, opts6)
      test_shock_capturing_jac(mesh6, sbp6, eqn6, opts6, capture)

      println("case8")
      sensor = EulerEquationMod.ShockSensorVelocity{Tsol, Tres}(mesh8, sbp8, opts8)
      capture = EulerEquationMod.ShockCapturingDict[sc]{Tsol, Tres}(mesh8,
                                                      sbp8, eqn8, opts8, sensor) 
      test_sbp_cartesian(mesh8, sbp8, eqn8, opts8)
      test_shock_capturing_jac(mesh8, sbp8, eqn8, opts8, capture)

      println("case9")
      sensor = EulerEquationMod.ShockSensorVelocity{Tsol, Tres}(mesh9, sbp9, opts9)
      capture = EulerEquationMod.ShockCapturingDict[sc]{Tsol, Tres}(mesh9,
                                                      sbp9, eqn9, opts9, sensor) 
      test_sbp_cartesian(mesh9, sbp9, eqn9, opts9)
      test_shock_capturing_jac(mesh9, sbp9, eqn9, opts9, capture)

    end


    # only need to do a single element type because all element volume matrices
    # have the same structure
    println("\nTesting volume shock capturing schemes")

    shock_capturing_schemes_volume = ["ElementProjection", "Volume"]
    for sc in shock_capturing_schemes_volume
      println("sc = ", sc)

      Tsol = eltype(eqn4.q); Tres = eltype(eqn4.res)

      println("case4")
      sensor = EulerEquationMod.ShockSensorVelocity{Tsol, Tres}(mesh4, sbp4, opts4)
      capture = EulerEquationMod.ShockCapturingDict[sc]{Tsol, Tres}(mesh4,
                                                    sbp4, eqn4, opts4, sensor)
      test_sbp_cartesian(mesh4, sbp4, eqn4, opts4)
      test_sbp_cartesian_revm(mesh4, sbp4, eqn4, opts4)
      test_shock_capturing_jac(mesh4, sbp4, eqn4, opts4, capture)

      if sc == "Volume"
        test_shock_capturing_jac(mesh4, sbp4, eqn4, opts4, capture, freeze=true)
      end
    end


    # test revm products

    # regular Roe scheme
    println("\n\nTesting Roe scheme")
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["IC_name"] = "ICExp"
    opts_tmp["BC1_name"] = "FreeStreamBC"  # BC with reverse mode functor
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["order"] = 2
    opts_tmp["need_adjoint"] = true
    make_input(opts_tmp, fname4)
    mesh_r1, sbp_r1, eqn_r1, opts_r1 = run_solver(fname4)

    test_revm_product(mesh_r1, sbp_r1, eqn_r1, opts_r1)
    test_revq_product(mesh_r1, sbp_r1, eqn_r1, opts_r1)

    # test 2d BC that depends on mesh.coords_bndry
    println("\n\nTesting 2D boundary condition")
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname)
    opts_tmp["IC_name"] = "ICExp"
    opts_tmp["BC1_name"] = "isentropicVortexBC"  # BC with reverse mode functor
    opts_tmp["operator_type"] = "SBPOmega"
    opts_tmp["order"] = 2
    opts_tmp["need_adjoint"] = true
    opts_tmp["addShockCapturing"] = true
    opts_tmp["shock_capturing_name"] = "ElementProjection"
    make_input(opts_tmp, fname4)
    mesh_r2, sbp_r2, eqn_r2, opts_r2 = run_solver(fname4)
    opts_r2["addShockCapturing"] = false

    test_revm_product(mesh_r2, sbp_r2, eqn_r2, opts_r2)
    test_revq_product(mesh_r2, sbp_r2, eqn_r2, opts_r2)

    # the ExpIC looks enough like a shock that the indicator says it is one
    opts_r2["addShockCapturing"] = true
    test_jac_general(mesh_r2, sbp_r2, eqn_r2, opts_r2)
    opts_r2["addShockCapturing"] = false



    # SBPOmega ES scheme
    println("\n\ntesting SBP Omega ES scheme")
    println("testing ES scheme")
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["IC_name"] = "ICExp"
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPOmega"
    opts_tmp["BC1_name"] = "FreeStreamBC"  # BC with reverse mode functor
    opts_tmp["volume_integral_type"] = 2
    opts_tmp["Volume_flux_name"] = "IRFlux"
    opts_tmp["face_integral_type"] = 2
    opts_tmp["Flux_name"] = "IRFlux"
    opts_tmp["FaceElementIntegral_name"] = "ESLFFaceIntegral"
    opts_tmp["need_adjoint"] = true
    make_input(opts_tmp, fname4)
    mesh_r3, sbp_r3, eqn_r3, opts_r3 = run_solver(fname4)

    test_revm_product(mesh_r3, sbp_r3, eqn_r3, opts_r3)
    test_revq_product(mesh_r3, sbp_r3, eqn_r3, opts_r3)

    # SBPDiagonalE ES scheme
    println("\n\ntesting diagonalE ES scheme")
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["IC_name"] = "ICExp"
    opts_tmp["BC1_name"] = "FreeStreamBC"  # BC with reverse mode functor
    opts_tmp["volume_integral_type"] = 2
    opts_tmp["Volume_flux_name"] = "IRFlux"
    opts_tmp["Flux_name"] = "IRSLFFlux"
    opts_tmp["order"] = 2
    opts_tmp["need_adjoint"] = true
    make_input(opts_tmp, fname4)
    mesh_r4, sbp_r4, eqn_r4, opts_r4 = run_solver(fname4)

    test_revm_product(mesh_r4, sbp_r4, eqn_r4, opts_r4)
    test_revq_product(mesh_r4, sbp_r4, eqn_r4, opts_r4)

    # shock capturing
    println("\n\ntesting diagonalE ES scheme with shock capturing")
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["IC_name"] = "ICExp"
    opts_tmp["BC1_name"] = "FreeStreamBC"  # BC with reverse mode functor
    opts_tmp["volume_integral_type"] = 2
    opts_tmp["Volume_flux_name"] = "IRFlux"
    opts_tmp["Flux_name"] = "IRSLFFlux"
    opts_tmp["order"] = 1
    opts_tmp["need_adjoint"] = true
    opts_tmp["addShockCapturing"] = true
    opts_tmp["shock_capturing_name"] = "Volume"
    opts_tmp["shock_sensor_name"] = "SensorHHO"
    opts_tmp["addVolumeIntegrals"] = false
    opts_tmp["addFaceIntegrals"] = false
    opts_tmp["addBoundaryIntegrals"] = false
    opts_tmp["use_lps"] = true
    make_input(opts_tmp, fname4)

    mesh_r5, sbp_r5, eqn_r5, opts_r5 = run_solver(fname4)

   test_revm_product(mesh_r5, sbp_r5, eqn_r5, opts_r5)
   test_revm_product(mesh_r5, sbp_r5, eqn_r5, opts_r5; freeze=true)
   test_revq_product(mesh_r5, sbp_r5, eqn_r5, opts_r5)
   println("\nTesting revq frozen")
   test_revq_product(mesh_r5, sbp_r5, eqn_r5, opts_r5; freeze=true)

  end

  return nothing
end

add_func1!(EulerTests, test_jac_terms_long, [TAG_LONGTEST, TAG_JAC])


#------------------------------------------------------------------------------
# functions that test individual functionality

"""
  Returns an array of the specified size with random values for the real part
  and zeros for the imaginary part
"""
function rand_realpart(dims...)

  a = rand(Complex128, dims...)
  for i=1:length(a)
    a[i] = real(a[i])
  end

  return a
end

function rand_realpart!(a::AbstractArray)

  rand!(a)
  for i=1:length(a)
    a[i] = real(a[i])
  end
end


function test_pressure(params::AbstractParamType{Tdim}) where Tdim


  if Tdim == 2
    q = Complex128[2.0, 3.0, 4.0, 7.0]
  else
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
  end

  numDofPerNode = length(q)

  p_dot = zeros(q)

  h = 1e-20
  pert = Complex128(0, h)
  for j=1:numDofPerNode
    q[j] += pert
    p = EulerEquationMod.calcPressure(params, q)
    p_dot[j] = imag(p)/h
    q[j] -= pert
  end

  p_dot2 = zeros(q)
  EulerEquationMod.calcPressure_diff(params, q, p_dot2)

  @test isapprox( maximum(abs.(p_dot - p_dot2)), 0.0) atol=1e-13

  # test reverse mode
  q_bar = zeros(q)

  EulerEquationMod.calcPressure_revq(params, q, q_bar, 1.0)

  @test maximum(abs.(q_bar - p_dot)) < 1e-13


  # test hessian
  p_dot = zeros(Complex128, numDofPerNode)
  p_hess = zeros(Complex128, numDofPerNode, numDofPerNode)
  p_hessc = zeros(Complex128, numDofPerNode, numDofPerNode)

  h = 1e-20 
  pert = Complex128(0, h)

  for i=1:numDofPerNode
    q[i] += pert
    EulerEquationMod.calcPressure_diff(params, q, p_dot)
    p_hessc[:, i] = imag(p_dot)./h
    q[i] -= pert
  end

  EulerEquationMod.calcPressure_hess(params, q, p_hess)

  @test maximum(abs.(p_hess - p_hessc)) < 1e-13


  # test calcPressure_diff_revq
  # use hessian to test reverse mode of p_dot, use regular calcPressure_revq
  # to test reverse mode of pressure
  p_dot_bar = rand_realpart(size(q))

  q_bar = p_hess.'*p_dot_bar

  q_bar2 = zeros(Complex128, length(q))
  EulerEquationMod.calcPressure_diff_revq(params, q, q_bar2, p_dot_bar, 0)

  @test maximum(abs.(q_bar - q_bar2)) < 1e-13

  fill!(q_bar, 0); fill!(q_bar2, 0); fill!(p_dot_bar, 0)

  EulerEquationMod.calcPressure_revq(params, q, q_bar, 2.0)
  EulerEquationMod.calcPressure_diff_revq(params, q, q_bar2, p_dot_bar, 2.0)

  @test maximum(abs.(q_bar - q_bar2)) < 1e-13


  return nothing
end



"""
  Tests the differentiated version of a single point flux
  function (eg. the Euler flux)
"""
function test_eulerflux(params::AbstractParamType{Tdim}) where Tdim

  if Tdim == 2
    nrm = [0.45, 0.55]
    q = Complex128[2.0, 3.0, 4.0, 7.0]
  else
    nrm = [0.45, 0.55, 0.65]
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
  end

  numDofPerNode = length(q)

  aux_vars = Complex128[0.0]
  F = zeros(Complex128, numDofPerNode)

  res = zeros(Complex128, numDofPerNode, numDofPerNode)

  h =1e-20
  pert = Complex128(0, h)
  for j=1:numDofPerNode
    q[j] += pert
    EulerEquationMod.calcEulerFlux(params, q, aux_vars, nrm, F)

    for k=1:numDofPerNode
      res[k, j] = imag(F[k])/h
    end

    q[j] -= pert
  end


  res2 = zeros(res)
  EulerEquationMod.calcEulerFlux_diff(params, q, aux_vars, nrm, res2)

  @test isapprox( maximum(abs.(res - res2)), 0.0) atol=1e-13

  # test that res2 is summed into
  res2_orig = copy(res2)
  EulerEquationMod.calcEulerFlux_diff(params, q, aux_vars, nrm, res2)

  @test maximum(abs.(res2 - 2*res2_orig)) < 1e-13
end

  
function test_eulerflux_revm(params::AbstractParamType{Tdim}) where Tdim

  if Tdim == 2
    nrm = Complex128[0.45, 0.55]
    q = Complex128[2.0, 3.0, 4.0, 7.0]
  else
    nrm = Complex128[0.45, 0.55, 0.65]
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
  end

  numDofPerNode = length(q)

  aux_vars = Complex128[0.0]
  
  F_bar = rand_realpart(numDofPerNode)
  F_dot = zeros(Complex128, numDofPerNode)
  nrm_bar = zeros(Complex128, Tdim)
  nrm_dot = rand_realpart(Tdim)


  h =1e-20
  pert = Complex128(0, h)


  nrm .+= pert*nrm_dot
  EulerEquationMod.calcEulerFlux(params, q, aux_vars, nrm, F_dot)
  val = sum(imag(F_dot)/h .* F_bar)
  nrm .-= pert*nrm_dot

  EulerEquationMod.calcEulerFlux_revm(params, q, aux_vars,
                                      nrm, nrm_bar, F_bar)

  val2 = sum(nrm_bar.*nrm_dot)
  @test abs(val - val2) < 1e-13

  # test accumulation behavior
  nrm_bar_orig = copy(nrm_bar)
  EulerEquationMod.calcEulerFlux_revm(params, q, aux_vars,
                                      nrm, nrm_bar, F_bar)
  @test maximum(abs.(nrm_bar - 2*nrm_bar_orig)) < 1e-13


  return nothing
end


function test_eulerflux_revq(params::AbstractParamType{Tdim}) where Tdim

  if Tdim == 2
    nrm = Complex128[0.45, 0.55]
    q = Complex128[2.0, 3.0, 4.0, 7.0]
  else
    nrm = Complex128[0.45, 0.55, 0.65]
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
  end

  numDofPerNode = length(q)

  aux_vars = Complex128[0.0]
  
  F_bar = rand_realpart(numDofPerNode)
  F_dot = zeros(Complex128, numDofPerNode)
  q_dot = rand_realpart(numDofPerNode)
  q_bar = zeros(Complex128, numDofPerNode)


  h =1e-20
  pert = Complex128(0, h)
  q .+= pert*q_dot
  EulerEquationMod.calcEulerFlux(params, q, aux_vars, nrm, F_dot)
  val = sum(imag(F_dot)/h .* F_bar)
  q .-= pert*q_dot


  EulerEquationMod.calcEulerFlux_revq(params, q, q_bar, aux_vars, nrm, F_bar)
  val2 = sum(q_bar .* q_dot)

  @test abs(val - val2) < 1e-13


  q_bar_orig = copy(q_bar)
  EulerEquationMod.calcEulerFlux_revq(params, q, q_bar, aux_vars, nrm, F_bar)
  @test maximum(abs.(2*q_bar_orig - q_bar)) < 1e-13

  return nothing
end


"""
  Test getLambdaMax differentiated versions
"""
function test_lambda(params::AbstractParamType{Tdim}, qL::AbstractVector,
                     nrm::AbstractVector) where Tdim


  lambda_dot = zeros(qL)
  lambda_dot3 = zeros(qL)
  EulerEquationMod.getLambdaMax_diff(params, qL, nrm, lambda_dot)
  EulerEquationMod.getLambdaMax_diff(params, qL, lambda_dot3)

  lambda_dot2 = zeros(qL)
  lambda_dot4 = zeros(qL)
  h=1e-20
  pert = Complex128(0, h)
  for i=1:length(lambda_dot)
    println("i = ", i)
    qL[i] += pert
    lambda_dot2[i] = imag(EulerEquationMod.getLambdaMax(params, qL, nrm))/h
    lambda_dot4[i] = imag(EulerEquationMod.getLambdaMax(params, qL))/h
    qL[i] -= pert
  end

  @test isapprox( norm(lambda_dot - lambda_dot2), 0.0) atol=1e-13
  @test isapprox( norm(lambda_dot3 - lambda_dot4), 0.0) atol=1e-13

  # test accumulation behavior
  lambda_dot_orig = copy(lambda_dot)
  lambda_dot3_orig = copy(lambda_dot3)
  EulerEquationMod.getLambdaMax_diff(params, qL, nrm, lambda_dot)
  EulerEquationMod.getLambdaMax_diff(params, qL, lambda_dot3)
  @test isapprox( norm(lambda_dot - lambda_dot_orig), 0.0) atol=1e-13
  @test isapprox( norm(lambda_dot3 - lambda_dot3_orig), 0.0) atol=1e-13

  # revq
  qL_bar = zeros(qL)
  qL2_bar = zeros(qL)
  qL_barc = zeros(qL)
  qL2_barc = zeros(qL)

  for i=1:length(qL)
    qL[i] += pert
    qL_barc[i] = imag(EulerEquationMod.getLambdaMax(params, qL, nrm))/h
    qL2_barc[i] = imag(EulerEquationMod.getLambdaMax(params, qL))/h
    qL[i] -= pert
  end

  EulerEquationMod.getLambdaMax_revq(params, qL, qL_bar, nrm, 2)
  EulerEquationMod.getLambdaMax_revq(params, qL, qL2_bar, 2)

  @test maximum(abs.(qL_bar - 2*qL_barc)) < 1e-14
  @test maximum(abs.(qL2_bar - 2*qL2_barc)) < 1e-14

  # test accumulation behavior
  qL_bar_orig = copy(qL_bar)
  qL2_bar_orig = copy(qL2_bar)
  EulerEquationMod.getLambdaMax_revq(params, qL, qL_bar, nrm, 2)
  EulerEquationMod.getLambdaMax_revq(params, qL, qL2_bar, 2)
  @test maximum(abs.(qL_bar - 2*qL_bar_orig)) < 1e-14
  @test maximum(abs.(qL2_bar - 2*qL2_bar_orig)) < 1e-14


  # revm
  nrmc = zeros(Complex128, length(nrm)); copy!(nrmc, nrm)
  nrm_bar = zeros(nrmc)
  nrm_barc = zeros(nrmc)

  for i=1:length(nrm)
    nrmc[i] += pert
    nrm_barc[i] = imag(EulerEquationMod.getLambdaMax(params, qL, nrmc))/h
    nrmc[i] -= pert
  end

  EulerEquationMod.getLambdaMax_revm(params, qL, nrmc, nrm_bar, 1)

  @test maximum(abs.(nrm_bar - nrm_barc)) < 1e-14

  nrm_bar_orig = copy(nrm_bar)
  EulerEquationMod.getLambdaMax_revm(params, qL, nrmc, nrm_bar, 1)
  @test maximum(abs.(nrm_bar - 2*nrm_bar_orig)) < 1e-14

  return nothing
end

function test_lambdasimple(params::AbstractParamType{Tdim}, qL::AbstractVector,
                           qR::AbstractVector,
                           nrm::AbstractVector) where Tdim


  lambda_dotL = zeros(qL)
  lambda_dotR = zeros(qL)
  EulerEquationMod.getLambdaMaxSimple_diff(params, qL, qR, nrm, lambda_dotL, lambda_dotR)

  lambda_dotL2 = zeros(qL)
  lambda_dotR2 = zeros(qL)
  h = 1e-20
  pert = Complex128(0, h)
  for i=1:length(lambda_dotL)
    qL[i] += pert
    lambda_dotL2[i] = imag(EulerEquationMod.getLambdaMaxSimple(params, qL, qR, nrm))/h
    qL[i] -= pert
  end

  for i=1:length(lambda_dotL)
    qR[i] += pert
    lambda_dotR2[i] = imag(EulerEquationMod.getLambdaMaxSimple(params, qL, qR, nrm))/h
    qR[i] -= pert
  end


  @test isapprox( norm(lambda_dotL - lambda_dotL2), 0.0) atol=1e-13
  @test isapprox( norm(lambda_dotR - lambda_dotR2), 0.0) atol=1e-13

  # revm
  nrmc = zeros(Complex128, length(nrm)); copy!(nrmc, nrm)
  nrm_bar = zeros(nrmc)
  nrm_barc = zeros(nrmc)

  for i=1:length(nrm)
    nrmc[i] += pert
    nrm_barc[i] = imag(EulerEquationMod.getLambdaMaxSimple(params, qL, qR, nrmc))/h
    nrmc[i] -= pert
  end

  EulerEquationMod.getLambdaMaxSimple_revm(params, qL, qR, nrmc, nrm_bar, 1)

  @test maximum(abs.(nrm_bar - nrm_barc)) < 1e-14

  nrm_bar_orig = copy(nrm_bar)
  EulerEquationMod.getLambdaMaxSimple_revm(params, qL, qR, nrmc, nrm_bar, 1)
  @test maximum(abs.(nrm_bar - 2*nrm_bar_orig)) < 1e-14


  # revq
  qL_bar = zeros(qL)
  qR_bar = zeros(qR)
  qL_barc = zeros(qL)
  qR_barc = zeros(qR)

  for i=1:length(qL)
    qL[i] += pert
    qL_barc[i] = imag(EulerEquationMod.getLambdaMaxSimple(params, qL, qR, nrm))/h
    qL[i] -= pert

    qR[i] += pert
    qR_barc[i] = imag(EulerEquationMod.getLambdaMaxSimple(params, qL, qR, nrm))/h
    qR[i] -= pert
  end

  EulerEquationMod.getLambdaMaxSimple_revq(params, qL, qL_bar, qR, qR_bar,
                                           nrm, 2)

  @test maximum(abs.(qL_bar - 2*qL_barc)) < 1e-14
  @test maximum(abs.(qR_bar - 2*qR_barc)) < 1e-14

  #TODO: test accumulation

  qL_bar_orig = copy(qL_bar); qR_bar_orig = copy(qR_bar)
  EulerEquationMod.getLambdaMaxSimple_revq(params, qL, qL_bar, qR, qR_bar,
                                           nrm, 2)
  @test maximum(abs.(qL_bar - 2*qL_bar_orig)) < 1e-14
  @test maximum(abs.(qR_bar - 2*qR_bar_orig)) < 1e-14



  return nothing
end


"""
  Test the differentiated versions of IRA0
"""
function test_IRA0(params::AbstractParamType{Tdim}) where {Tdim}

  h = 1e-20
  pert = Complex128(0, h)

  if Tdim == 2
    q = Complex128[1.0, 0.3, 0.4, 7.0]
  else
    q = Complex128[1.0, 0.3, 0.4, 0.5, 13.0]
  end

  numDofPerNode = length(q)
  nd = numDofPerNode*2
  pert_vec = rand_realpart((numDofPerNode, nd))


  for i=1:2  # run test twice to make sure internal arrays are zeroed out
    A0c = zeros(Complex128, numDofPerNode, numDofPerNode)
    A0_dotc = zeros(Float64, numDofPerNode, numDofPerNode, nd)

    A0 = zeros(Complex128, numDofPerNode, numDofPerNode)
    A0_dot = zeros(Complex128, numDofPerNode, numDofPerNode, nd)
    q_dot = zeros(Complex128, numDofPerNode, nd)


    for i=1:nd
      q .+= pert*pert_vec[:, i]
      EulerEquationMod.getIRA0(params, q, A0c)
      q .-= pert*pert_vec[:, i]

      A0_dotc[:, :, i] = imag(A0c)/h
    end

    EulerEquationMod.getIRA0_diff(params, q, pert_vec, A0, A0_dot)

    #println("maxdiff = ", maximum(abs.(A0_dot - A0_dotc)) )
    @test maximum(abs.(A0_dot - A0_dotc)) < 1e-12

    # revq
    # The jacobian is 3d, so the test is val = A0_ijk * v_i * w_j * c_k
    # For complex step, use c_k as perturbation seed, sum v_i w_j afterward
    # For reverse mode, use v_i w_j as A0_bar, sum with c_k afterwards
    vi = rand(numDofPerNode)
    wj = rand(numDofPerNode)
    ck = rand(numDofPerNode)
    A0 = zeros(Complex128, numDofPerNode, numDofPerNode)
    A02 = zeros(Complex128, numDofPerNode, numDofPerNode)
    A0_bar = zeros(Complex128, numDofPerNode, numDofPerNode)
    q_bar = zeros(Complex128, numDofPerNode)

    # complex step
    for i=1:numDofPerNode
      q[i] += pert*ck[i]
    end

    EulerEquationMod.getIRA0(params, q, A0)
    A0_dot = imag(A0)/h

    for i=1:numDofPerNode
      q[i] -= pert*ck[i]
    end

    val = 0.0
    for i=1:numDofPerNode
      for j=1:numDofPerNode
        val += A0_dot[i, j]*vi[i]*wj[j]
      end
    end

    # reverse mode
    for i=1:numDofPerNode
      for j=1:numDofPerNode
        A0_bar[i, j] = vi[i]*wj[j]
      end
    end
    EulerEquationMod.getIRA0_revq(params, q, q_bar, A02, A0_bar)

    val2 = sum(q_bar.*ck)

    @test abs(val - val2) < 1e-12
    @test maximum(abs.(A0 - A02)) < 1e-13

    # test accumulation
    q_bar_orig = copy(q_bar)
    EulerEquationMod.getIRA0_revq(params, q, q_bar, A02, A0_bar)
    @test maximum(abs.(q_bar - 2*q_bar_orig)) < 1e-13

  end

  return nothing
end


"""
  Test the differentiated version of an `AbstractEntropyKernel`
"""
function test_EntropyKernel(params::AbstractParamType{Tdim},
                  kernel::EulerEquationMod.AbstractEntropyKernel) where {Tdim}


  h = 1e-20
  pert = Complex128(0, h)

  if Tdim == 2
    q = Complex128[1.0, 0.3, 0.4, 7.0]
    delta_w = Complex128[0.1, 0.2, 0.3, 0.4]
    nrm = [1.0, 2.0]
  else
    q = Complex128[1.0, 0.3, 0.4, 0.5, 13.0]
    delta_w = Complex128[0.1, 0.2, 0.3, 0.4, 0.5]
    nrm = [1.0, 2.0, 3.0]
  end

  numDofPerNode = length(q)
  nd = 2*numDofPerNode
  pert_vec_q = rand_realpart((numDofPerNode, nd))
  pert_vec_w = rand_realpart((numDofPerNode, nd))
  flux = zeros(Complex128, numDofPerNode)
  fluxc = zeros(Complex128, numDofPerNode)
  flux_dotc = zeros(Complex128, numDofPerNode, nd)
  flux_dot = zeros(Complex128, numDofPerNode, nd)

  for i=1:2  # run test twice to make sure intermediate arrays are zeroed out
    fill!(flux_dot, 0.0)
    fill!(flux_dotc, 0.0)
    for i=1:nd
      q .+= pert*pert_vec_q[:, i]
      delta_w .+= pert*pert_vec_w[:, i]
      EulerEquationMod.applyEntropyKernel(kernel, params, q, delta_w, nrm, fluxc)
      q .-= pert*pert_vec_q[:, i]
      delta_w .-= pert*pert_vec_w[:, i]

      flux_dotc[:, i] = imag(fluxc)/h
    end

    EulerEquationMod.applyEntropyKernel_diff(kernel, params, q, pert_vec_q, delta_w, pert_vec_w, nrm, flux, flux_dot)

    @test maximum(abs.(flux - real(fluxc))) < 1e-14
    @test maximum(abs.(flux_dot - flux_dotc)) < 3e-12



  end

  return nothing
end


function test_EntropyKernel_revq(params::AbstractParamType{Tdim},
                  kernel::EulerEquationMod.AbstractEntropyKernel) where {Tdim}

  h = 1e-20
  pert = Complex128(0, h)

  if Tdim == 2
    q = Complex128[1.0, 0.3, 0.4, 7.0]
    delta_w = Complex128[0.1, 0.2, 0.3, 0.4]
    nrm = [1.0, 2.0]
  else
    q = Complex128[1.0, 0.3, 0.4, 0.5, 13.0]
    delta_w = Complex128[0.1, 0.2, 0.3, 0.4, 0.5]
    nrm = [1.0, 2.0, 3.0]
  end

  numDofPerNode = length(q)

  for i=1:2  # make sure intermeidate arrays are zeroed out
    q_dot = rand_realpart(numDofPerNode)
    delta_w_dot = rand_realpart(numDofPerNode)
    q_bar = zeros(Complex128, numDofPerNode)
    delta_w_bar = zeros(Complex128, numDofPerNode)
    F_bar = rand_realpart(numDofPerNode)
    flux = zeros(Complex128, numDofPerNode)
    fluxc = zeros(Complex128, numDofPerNode)

    q .+= pert*q_dot
    delta_w .+= pert*delta_w_dot
    EulerEquationMod.applyEntropyKernel(kernel, params, q, delta_w, nrm, fluxc)
    q .-= pert*q_dot
    delta_w .-= pert*delta_w_dot
    val = sum(F_bar.*imag(fluxc)/h)

    EulerEquationMod.applyEntropyKernel_revq(kernel, params, q, q_bar, delta_w, delta_w_bar, nrm, flux, F_bar)

    val2 = sum(delta_w_bar.*delta_w_dot +  q_bar.*q_dot)

    @test abs(val - val2) < 3e-12
    @test maximum(abs.(flux - fluxc)) < 1e-13

    # test accumulation
    q_bar_orig = copy(q_bar)
    delta_w_bar_orig = copy(delta_w_bar)

    EulerEquationMod.applyEntropyKernel_revq(kernel, params, q, q_bar, delta_w, delta_w_bar, nrm, flux, F_bar)

    @test maximum(abs.(flux - fluxc)) < 1e-12
    @test maximum(abs.(q_bar - 2*q_bar_orig)) < 1e-12
    @test maximum(abs.(delta_w_bar - 2*delta_w_bar_orig)) < 1e-12
  end

  return nothing
end


function test_EntropyKernel_revm(params::AbstractParamType{Tdim},
                  kernel::EulerEquationMod.AbstractEntropyKernel) where {Tdim}

  h = 1e-20
  pert = Complex128(0, h)

  if Tdim == 2
    q = Complex128[1.0, 0.3, 0.4, 7.0]
    delta_w = Complex128[0.1, 0.2, 0.3, 0.4]
    nrm = [1.0, 2.0]
  else
    q = Complex128[1.0, 0.3, 0.4, 0.5, 13.0]
    delta_w = Complex128[0.1, 0.2, 0.3, 0.4, 0.5]
    nrm = [1.0, 2.0, 3.0]
  end

  numDofPerNode = length(q)

  for i=1:2

    nrmc = zeros(Complex128, length(nrm)); copy!(nrmc, nrm)
    F_bar = rand_realpart(numDofPerNode)
    flux = zeros(Complex128, numDofPerNode)
    nrm_dot = rand_realpart(length(nrm))
    nrm_bar = zeros(Complex128, length(nrm))

    nrmc .+= pert*nrm_dot
    EulerEquationMod.applyEntropyKernel(kernel, params, q, delta_w, nrmc, flux)
    nrmc .-= pert*nrm_dot
    val = sum(F_bar.*imag(flux)/h)

    EulerEquationMod.applyEntropyKernel_revm(kernel, params, q, delta_w, nrmc, nrm_bar, flux, F_bar)
    val2 = sum(nrm_bar .* nrm_dot)

    @test abs(val - val2) < 1e-12

    # test accumulation
    nrm_bar_orig = copy(nrm_bar)
    flux_orig = copy(flux)
    EulerEquationMod.applyEntropyKernel_revm(kernel, params, q, delta_w, nrmc, nrm_bar, flux, F_bar)

    @test maximum(abs.(flux - flux_orig)) < 1e-13
    @test maximum(abs.(nrm_bar - 2*nrm_bar_orig)) < 1e-12

  end

  return nothing
end



"""
  Test a differentiated numerical flux function via complex step of the
  original function

  **Inputs**

   * params: a Params object
   * qL: left state
   * qR: right state
   * nrm: normal vector
   * func: the original function
   * func_diff: the differentiated version
"""
function test_ad_inner(params::AbstractParamType{Tdim}, qL, qR, nrm,
                       func, func_diff, output=false) where Tdim

  # compute jacobian with complex step and AD, compare results
  numDofPerNode = length(qL)

  aux_vars = Complex128[0.0]
  F = zeros(Complex128, numDofPerNode)

  resL = zeros(Complex128, numDofPerNode, numDofPerNode)
  resR = zeros(Complex128, numDofPerNode, numDofPerNode)
  h = 1e-20
  pert = Complex128(0, h)

  for j=1:numDofPerNode
    qL[j] += pert
    func(params, qL, qR, aux_vars, nrm, F)

    for k=1:numDofPerNode
      resL[k, j] = imag(F[k])/h
    end

    qL[j] -= pert
  end


  for j=1:numDofPerNode
    qR[j] += pert
    func(params, qL, qR, aux_vars, nrm, F)

    for k=1:numDofPerNode
      resR[k, j] = imag(F[k])/h
    end

    qR[j] -= pert
  end

  # AD version
  resL2 = zeros(resL)
  resR2 = zeros(resR)

  func_diff(params, qL, qR, aux_vars, nrm, resL2, resR2)

  if output
    println("\n----- Comparing left results -----")
    println("resL = \n", real(resL))
    println("resL2 = \n", real(resL2))
     println("diff = \n", real(resL - resL2))
    println("max diff = \n", maximum(abs.(resL - resL2)))

    println("\n----- Comparing right results -----")
    println("resR = \n", real(resR))
    println("resR2 = \n", real(resR2))
    println("diff = \n", real(resR - resR2))
    println("max diff = \n", maximum(abs.(resR - resR2)))
  end


  @test isapprox( maximum(abs.(resL - resL2)), 0.0) atol=1e-12
  @test isapprox( maximum(abs.(resR - resR2)), 0.0) atol=1e-12

  # test resL,R are summed into
  resL2_orig = copy(resL2)
  resR2_orig = copy(resR2)
  func_diff(params, qL, qR, aux_vars, nrm, resL2, resR2)

  @test maximum(abs.(resL2 - 2*resL2_orig)) < 1e-12
  @test maximum(abs.(resR2 - 2*resR2_orig)) < 1e-12

  return nothing
end

"""
  Test the reverse mode with respect to q of a 2 point flux function

  **Inputs**
  
   * params: a ParamType
   * qL: left state
   * qR: right state
   * nrm: normal vector
   * func: flux function
   * func_revq: reverse mode (wrt q) flux function
"""
function test_2flux_revq(params::AbstractParamType{Tdim}, qL, qR, nrm, func,
                         func_revq; test_multid=false) where {Tdim}


  h = 1e-20
  pert = Complex128(0, h)

  numDofPerNode = length(qL)
  # test the single direction version
  qL_dot = rand_realpart(size(qL))
#  qL_dot = zeros(numDofPerNode); qL_dot[1] = 1
  qR_dot = rand_realpart(size(qR))
#  qR_dot = zeros(numDofPerNode)
  qL_bar = zeros(qL)
  qR_bar = zeros(qR)
  aux_vars = Complex128[]
  flux = zeros(Complex128, numDofPerNode)
  F_bar = rand_realpart(numDofPerNode)
#   F_bar = zeros(Complex128, numDofPerNode); F_bar[1] = 1


  for i=1:1 #2  # run test twice to make sure all intermediate arrays are zeroed out
    # compute F_bar.'* dF/dq * q_dot using forward and reverse mode
    qL += pert*qL_dot
    func(params, qL, qR, aux_vars, nrm, flux)
    qL -= pert*qL_dot
    val_c = sum(F_bar.*imag(flux)/h)

    qR += pert*qR_dot
    func(params, qL, qR, aux_vars, nrm, flux)
    qR -= pert*qR_dot
    val_c += sum(F_bar.*imag(flux)/h)

    func_revq(params, qL, qL_bar, qR, qR_bar, aux_vars, nrm, F_bar)
    val = sum(qL_bar.*qL_dot) + sum(qR_bar.*qR_dot)

    println("val_c = ", val_c)
    println("val = ", val)
    @test abs(val - val_c) < 1e-12

    # test qL_bar is summed into
    qL_bar_orig = copy(qL_bar)
    qR_bar_orig = copy(qR_bar)
    func_revq(params, qL, qL_bar, qR, qR_bar, aux_vars, nrm, F_bar)

    @test maximum(abs.(qL_bar - 2*qL_bar_orig)) < 1e-12
    @test maximum(abs.(qR_bar - 2*qR_bar_orig)) < 1e-12

    fill!(qL_bar, 0.0); fill!(qR_bar, 0.0)
  end

  if test_multid
    flux = zeros(Complex128, numDofPerNode, Tdim)
    F_bar = rand_realpart(numDofPerNode, Tdim)
    nrm2 = zeros(Complex128, Tdim, Tdim)
    nrm2[:, 1] = nrm
    for i=2:Tdim
      nrm2[:, i] = nrm2[:, i-1] + 1
    end


    for i=1:2
      qL += pert*qL_dot
      func(params, qL, qR, aux_vars, nrm2, flux)
      qL -= pert*qL_dot
      val_c = 0.0
      for i=1:Tdim
        val_c += sum(F_bar[:, i].*imag(flux[:, i])/h)
      end

      qR += pert*qR_dot
      func(params, qL, qR, aux_vars, nrm2, flux)
      qR -= pert*qR_dot
      for i=1:Tdim
        val_c += sum(F_bar[:, i].*imag(flux[:, i])/h)
      end

      func_revq(params, qL, qL_bar, qR, qR_bar, aux_vars, nrm2, F_bar)
      val = sum(qL_bar.*qL_dot) + sum(qR_bar.*qR_dot)

      @test abs(val_c - val) < 1e-12

      # test qL_bar is summed into
      qL_bar_orig = copy(qL_bar)
      qR_bar_orig = copy(qR_bar)
      func_revq(params, qL, qL_bar, qR, qR_bar, aux_vars, nrm2, F_bar)

      @test maximum(abs.(qL_bar - 2*qL_bar_orig)) < 1e-12
      @test maximum(abs.(qR_bar - 2*qR_bar_orig)) < 1e-12

      fill!(qL_bar, 0.0); fill!(qR_bar, 0.0)
    end
  end

  return nothing
end


"""
  Test the reverse mode with respect to nrm of a 2 point flux function

  **Inputs**
  
   * params: a ParamType
   * qL: left state
   * qR: right state
   * nrm: normal vector
   * func: flux function
   * func_revm: reverse mode (wrt nrm) flux function
"""
function test_2flux_revm(params::AbstractParamType{Tdim}, qL, qR, nrm,
                         func::FluxType, func_revm::FluxType_revm;
                         test_multid=false) where {Tdim}

  h = 1e-20
  pert = Complex128(0, h)

  numDofPerNode = length(qL)
  flux = zeros(Complex128, numDofPerNode)
  F_bar = rand_realpart(numDofPerNode)
  nrm_dot = rand_realpart(Tdim)
  nrm_bar = zeros(Complex128, Tdim)
  aux_vars = Complex128[]


  # compute nrm_bar * df/dnrm * nrm_dot using forward and reverse mode
  for i=1:2  # run test twice to make sure intermeidate arrays are zeroed out
    nrm += pert*nrm_dot
    func(params, qL, qR, aux_vars, nrm, flux)
    nrm -= pert*nrm_dot
    val_c = sum(F_bar.*imag(flux)/h)

    func_revm(params, qL, qR, aux_vars, nrm, nrm_bar, F_bar)

    val = sum(nrm_bar.*nrm_dot)

    @test abs(val - val_c) < 1e-13

    # test nrm_bar is summed into
    nrm_bar_orig = copy(nrm_bar)
    func_revm(params, qL, qR, aux_vars, nrm, nrm_bar, F_bar)

    @test maximum(abs.(nrm_bar - 2*nrm_bar_orig)) < 1e-13

    fill!(nrm_bar, 0.0)
  end

  if test_multid
    flux = zeros(Complex128, numDofPerNode, Tdim)
    nrm2 = zeros(Complex128, Tdim, Tdim)
    nrm2[:, 1] = nrm
    for i=2:Tdim
      nrm2[:, i] = nrm2[:, i-1] + 1
    end
    nrm2_bar = zeros(nrm2)

    flux = zeros(Complex128, numDofPerNode, Tdim)
    F_bar = rand_realpart(numDofPerNode, Tdim)
    nrm2_dot = rand_realpart(Tdim, Tdim)


    for i=1:2
      nrm2 += pert*nrm2_dot
      func(params, qL, qR, aux_vars, nrm2, flux)
      nrm2 -= pert*nrm2_dot
      val_c = sum(F_bar.*imag(flux)/h)

      func_revm(params, qL, qR, aux_vars, nrm2, nrm2_bar, F_bar)

      val = sum(nrm2_bar.*nrm2_dot)

      @test abs(val - val_c) < 1e-13

      # test nrm_bar is summed into
      nrm_bar_orig = copy(nrm2_bar)
      func_revm(params, qL, qR, aux_vars, nrm2, nrm2_bar, F_bar)

      @test maximum(abs.(nrm2_bar - 2*nrm_bar_orig)) < 1e-13


      fill!(nrm2_bar, 0.0)
    end
  end

  return nothing
end


"""
  Test the reverse mode of a function in common_funcs.jl
"""
function test_common_func_rev(params::AbstractParamType{Tdim},
                              coords::AbstractArray{Tmsh, 1}, func, func_rev
                             ) where {Tdim, Tmsh}

  println("testing common function ", func)

  h = 1e-20
  pert = Complex128(0, h)


  coords_dot = rand_realpart(Tdim)
  coords_bar = zeros(Tmsh, Tdim)
  q_bar = rand_realpart(params.numDofPerNode)
  q = zeros(Complex128, params.numDofPerNode)

  coords .+= pert*coords_dot
  func(params, coords, q)
  coords .-= pert*coords_dot
  val1 = sum(imag(q)/h .* q_bar)

  func_rev(params, coords, coords_bar, q_bar)
  val2 = sum(coords_bar .* coords_dot)

  @test abs(val1 - val2) < 1e-13

  # test accumulation
  coords_bar_orig = copy(coords_bar)
  func_rev(params, coords, coords_bar, q_bar)

  @test maximum(abs.(coords_bar - 2*coords_bar_orig)) < 1e-13


  return nothing
end


function test_2RWave(params::AbstractParamType{Tdim}, qL, qR, nrm) where Tdim

  # compute jacobian with complex step and AD, compare results
  numDofPerNode = length(qL)
  h = 1e-20
  pert = Complex128(0, h)

  sL_dot = zeros(Complex128, numDofPerNode, 2)
  sR_dot = zeros(Complex128, numDofPerNode, 2)

  sL_dot2 = zeros(Complex128, numDofPerNode, 2)
  sR_dot2 = zeros(Complex128, numDofPerNode, 2)


  for i=1:numDofPerNode
    
    qL[i] += pert
    sL, sR = EulerEquationMod.calc2RWaveSpeeds(params, qL, qR, nrm)
    sL_dot[i, 1] = imag(sL)/h
    sR_dot[i, 1] = imag(sR)/h
    qL[i] -= pert
    

    
    qR[i] += pert
    sL, sR = EulerEquationMod.calc2RWaveSpeeds(params, qL, qR, nrm)
    sL_dot[i, 2] = imag(sL)/h
    sR_dot[i, 2] = imag(sR)/h
    qR[i] -= pert
    
  end

  # call this twice to make sure intermediate arrays are zeroed out if needed
  EulerEquationMod.calc2RWaveSpeeds_diff(params, qL, qR, nrm, sL_dot2, sR_dot2)
  EulerEquationMod.calc2RWaveSpeeds_diff(params, qL, qR, nrm, sL_dot2, sR_dot2)

  maxdiffL = maximum(abs.(sL_dot - sL_dot2))
  maxdiffR = maximum(abs.(sR_dot - sR_dot2))

  @test maxdiffL < 1e-12
  @test maxdiffR < 1e-12


  return nothing
end


function test_2RWave_revq(params::AbstractParamType{Tdim}, qL, qR, nrm) where {Tdim}

  numDofPerNode = length(qL)
  h = 1e-20
  pert = Complex128(0, h)
  qL_dot = rand_realpart(length(qL))
  #qL_dot = zeros(numDofPerNode); qL_dot[1] = 1
  qR_dot = rand_realpart(length(qR))
  #qR_dot = zeros(numDofPerNode)
  qL_bar = zeros(Complex128, numDofPerNode)
  qR_bar = zeros(Complex128, numDofPerNode)

  sL_bar = rand()
  sR_bar = rand()
  #sR_bar = 0.0

  # complex step
  qL .+= pert*qL_dot
  qR .+= pert*qR_dot
  sL, sR = EulerEquationMod.calc2RWaveSpeeds(params, qL, qR, nrm)
  qL .-= pert*qL_dot
  qR .-= pert*qR_dot
  sL_dot = imag(sL)./h
  sR_dot = imag(sR)./h

  val1 = sL_dot*sL_bar + sR_dot*sR_bar
  # reverse mode
  EulerEquationMod.calc2RWaveSpeeds_revq(params, qL, qL_bar, qR, qR_bar, nrm,
                                         sL_bar, sR_bar)

  println("qL_bar = ", real(qL_bar))
  val2 = sum(qL_bar.*qL_dot) + sum(qR_bar .* qR_dot)

  println("val1 = ", val1)
  println("val2 = ", val2)
  @test abs(val1 - val2) < 1e-13

  return nothing

end


function test_2RWave_revm(params::AbstractParamType{Tdim}, qL, qR, _nrm) where {Tdim}

  nrm = zeros(Complex128, length(_nrm)); copy!(nrm, _nrm)
  numDofPerNode = length(qL)
  h = 1e-20
  pert = Complex128(0, h)
  nrm_dot = rand_realpart(length(nrm))
  nrm_bar = zeros(Complex128, length(nrm))

  sL_bar = rand()
  sR_bar = rand()

  # complex step
  nrm .+= pert*nrm_dot
  sL, sR = EulerEquationMod.calc2RWaveSpeeds(params, qL, qR, nrm)
  nrm .-= pert*nrm_dot

  sL_dot = imag(sL)./h
  sR_dot = imag(sR)./h

  val1 = sL_dot*sL_bar + sR_dot*sR_bar

  # reverse mode
  EulerEquationMod.calc2RWaveSpeeds_revm(params, qL, qR, nrm, nrm_bar,
                                         sL_bar, sR_bar)

  val2 = sum(nrm_bar.*nrm_dot)

  println("val1 = ", val1)
  println("val2 = ", val2)
  @test abs(val1 - val2) < 1e-13

  return nothing

end






"""
  Test reverse mode of BC functor
"""
function test_bc_revm(params::AbstractParamType{Tdim},
                      func::EulerEquationMod.BCType,
                      func_revm::EulerEquationMod.BCType_revm) where {Tdim}

  println("testing BC function ", func, " in ", Tdim, " dimensions")
  h = 1e-20
  pert = Complex128(0, h)

  if Tdim == 2
    q = Complex128[1.0, 0.3, 0.4, 7.0]
    nrm = Complex128[1.0, 2.0]
    coords = Complex128[1.0, 1.1]
  else
    q = Complex128[1.0, 0.3, 0.4, 0.5, 13.0]
    nrm = Complex128[1.0, 2.0, 3.0]
    coords = Complex128[1.0, 1.1, 1.2]
  end

  numDofPerNode = length(q)
  aux_vars = Complex128[]
  nrm_bar = zeros(Complex128, Tdim)
  nrm_dot = rand_realpart(Tdim)
  flux_bar = rand_realpart(numDofPerNode)
  flux_dot = zeros(Complex128, numDofPerNode)
  coords_bar = zeros(Complex128, Tdim)
  coords_dot = rand_realpart(Tdim)

  nrm .+= pert*nrm_dot
  coords .+= pert*coords_dot
  func(params, q, aux_vars, coords, nrm, flux_dot)
  nrm .-= pert*nrm_dot
  coords .-= pert*coords_dot
  val = sum(flux_bar .* imag(flux_dot)/h)

  func_revm(params, q, aux_vars, coords, coords_bar, nrm, nrm_bar, flux_bar)
  val2 = sum(nrm_bar .* nrm_dot + coords_bar .* coords_dot)

  @test abs(val - val2) < 1e-13

  # test accumulation
  nrm_bar_orig = copy(nrm_bar)
  func_revm(params, q, aux_vars, coords, coords_bar, nrm, nrm_bar, flux_bar)
 
  @test maximum(abs.(nrm_bar - 2*nrm_bar_orig)) < 1e-13

  return nothing
end


function test_bc_revq(params::AbstractParamType{Tdim},
                      func::EulerEquationMod.BCType,
                      func_revq::EulerEquationMod.BCType_revq) where {Tdim}

  println("testing BC function ", func, " in ", Tdim, " dimensions")
  h = 1e-20
  pert = Complex128(0, h)

  if Tdim == 2
    q = Complex128[1.0, 0.3, 0.4, 7.0]
    nrm = Complex128[1.0, 2.0]
    coords = Complex128[1.0, 1.1]
  else
    q = Complex128[1.0, 0.3, 0.4, 0.5, 13.0]
    nrm = Complex128[1.0, 2.0, 3.0]
    coords = Complex128[1.0, 1.1, 1.2]
  end

  numDofPerNode = length(q)
  aux_vars = Complex128[]
  q_dot = rand_realpart(numDofPerNode)
  q_bar = zeros(Complex128, numDofPerNode)
  flux_bar = rand_realpart(numDofPerNode)
  flux = zeros(Complex128, numDofPerNode)

  q .+= pert*q_dot
  func(params, q, aux_vars, coords, nrm, flux)
  q .-= pert*q_dot
  val1 = sum(flux_bar .* imag(flux)/h)

  func_revq(params, q, q_bar, aux_vars, coords, nrm, flux_bar)
  val2 = sum(q_bar .* q_dot)

  @test abs(val1 - val2) < 1e-13

  return nothing
end


function test_faceflux_revm(mesh, sbp, eqn, opts)

    h = 1e-20
    pert = Complex128(0, h)
    EulerEquationMod.init_revm(mesh, sbp, eqn, opts)

    fill!(mesh.nrm_face_bar, 0)
    nrm_face_dot = zeros(Complex128, size(mesh.nrm_face_bar))

    rand!(eqn.flux_face_bar)
    flux_face_dot = zeros(Complex128, size(eqn.nrm_face))

    mesh.nrm_face .+= pert*nrm_face_dot
    EulerEquationMod.calcFaceFlux(mesh, sbp, eqn, opts, eqn.flux_func, mesh.interfaces, eqn.flux_face)
    mesh.nrm_face .-= pert*nrm_face_dot
    val = sum(imag(eqn.flux_face)/h .* eqn.flux_face_bar)


    @assert opts["face_integral_type"] == 1
    EulerEquationMod.calcFaceFlux_revm(mesh, sbp, eqn, eqn.flux_func_bar,
                                       mesh.interfaces, eqn.flux_face_bar)

    val2 = sum(nrm_face_dot .* mesh.nrm_face_bar)

    @test abs(val - val2) < 1e-12


    fill!(mesh.nrm_face_bar, 0)
    fill!(eqn.flux_face_bar, 0)
    return nothing
  end


function test_logavg()

  h = 1e-20
  pert = Complex128(0, h)

  nd = 4
  data = EulerEquationMod.LogAvgData{Float64, Float64}(nd)
  aL_dot = zeros(nd); aR_dot = zeros(nd)
  a_avg_dotL = zeros(nd); a_avg_dotR = zeros(nd)

  aL = 1.0
  aR = 2.0


  a_dotL_cs = imag(EulerEquationMod.logavg(aL + pert, aR))/h
  a_dotR_cs = imag(EulerEquationMod.logavg(aL, aR + pert))/h
  fill!(aL_dot, 1)
  fill!(aR_dot, 1)
  EulerEquationMod.logavg_diff(data, aL, aL_dot, aR, aR_dot, a_avg_dotL, a_avg_dotR)
  for i=1:nd
    @test isapprox(a_dotL_cs, a_avg_dotL[i]) atol=1e-13
    @test isapprox(a_dotR_cs, a_avg_dotR[i]) atol=1e-13
  end

  # reverse mode
  aL_bar, aR_bar = EulerEquationMod.logavg_rev(aL, aR, 1.0)
  @test abs(a_dotL_cs - aL_bar) < 1e-13
  @test abs(a_dotR_cs - aR_bar) < 1e-13

  # test other branch of if statement
  aL = 1 + 10.0^-4
  aR = 1

  a_dotL_cs = imag(EulerEquationMod.logavg(aL + pert, aR))/h
  a_dotR_cs = imag(EulerEquationMod.logavg(aL, aR + pert))/h
  fill!(aL_dot, 1)
  fill!(aR_dot, 1)
  EulerEquationMod.logavg_diff(data, aL, aL_dot, aR, aR_dot, a_avg_dotL, a_avg_dotR)
  for i=1:nd
    @test isapprox(a_dotL_cs, a_avg_dotL[i]) atol=1e-13
    @test isapprox(a_dotR_cs, a_avg_dotR[i]) atol=1e-13
  end


  # reverse mode
  aL_bar, aR_bar = EulerEquationMod.logavg_rev(aL, aR, 1.0)
  @test abs(a_dotL_cs - aL_bar) < 1e-13
  @test abs(a_dotR_cs - aR_bar) < 1e-13


  return nothing
end


function test_faceElementIntegral(params::AbstractParamType{Tdim},
                   sbpface::AbstractFace, func::FluxType,
                   func_diff::FluxType_diff) where {Tdim}


  h = 1e-20
  pert = Complex128(0, h)


  if Tdim == 2
    q = Complex128[1.0, 0.3, 0.4, 7.0]
    nrm = [1.0, 2.0]
  else
    q = Complex128[1.0, 0.3, 0.4, 0.5, 13.0]
    nrm = [1.0, 2.0, 3.0]
  end

  qL = zeros(Complex128, params.numDofPerNode, params.numNodesPerElement)
  qR = zeros(Complex128, params.numDofPerNode, params.numNodesPerElement)
  nrm_face = zeros(Tdim, params.numNodesPerFace)
  for i=1:params.numNodesPerElement
    qL[:, i] = q + (i-1)*0.1
    qR[:, i] = q + i*0.1
  end

  for i=1:params.numNodesPerFace
    nrm_face[:, i] = nrm + (i-1)*0.1
  end
  aux_vars = zeros(Complex128, 0, 0)
  resL = zeros(qL)
  resR = zeros(qR)
  jacLL = zeros(Complex128, params.numDofPerNode, params.numDofPerNode,
                params.numNodesPerElement, params.numNodesPerElement)
  jacLR = copy(jacLL)
  jacRL = copy(jacLL)
  jacRR = copy(jacLL)

  #TODO: test all configurations
  iface = Interface(1, 2, 1, 1, 1)

  qL_dot = rand_realpart(size(qL))
  qR_dot = rand_realpart(size(qR))
#  qL_dot = zeros(Complex128, params.numDofPerNode, params.numNodesPerElement)
#  qR_dot = zeros(Complex128, params.numDofPerNode, params.numNodesPerElement)

  qL_dot[1, 1] = 1

  # complex step
  qL .+= pert*qL_dot
  qR .+= pert*qR_dot

  EulerEquationMod.calcECFaceIntegral(params, sbpface, iface, qL, qR, aux_vars,
                                      nrm_face, func, resL, resR)
  valLc = imag(resL)/h
  valRc = imag(resR)/h

  EulerEquationMod.calcECFaceIntegral_diff(params, sbpface, iface, qL, qR,
                       aux_vars, nrm_face, func_diff, jacLL, jacLR, jacRL, jacRR)


  valL = zeros(qL)
  valR = zeros(qR)

  for q=1:params.numNodesPerElement
    for p=1:params.numNodesPerElement
      for i=1:params.numDofPerNode
        for j=1:params.numDofPerNode
          valL[i, p] += jacLL[i, j, p, q]*qL_dot[j, q]
          valL[i, p] += jacLR[i, j, p, q]*qR_dot[j, q]
          valR[i, p] += jacRL[i, j, p, q]*qL_dot[j, q]
          valR[i, p] += jacRR[i, j, p, q]*qR_dot[j, q]
        end
      end
    end
  end

  @test maximum(abs.(valL - valLc)) < 1e-13
  @test maximum(abs.(valR - valRc)) < 1e-13

  return nothing
end


function test_entropyPenalty(params::AbstractParamType{Tdim},
                   sbpface::AbstractFace,
                   kernel::EulerEquationMod.AbstractEntropyKernel) where {Tdim}


  h = 1e-20
  pert = Complex128(0, h)


  if Tdim == 2
    q = Complex128[1.0, 0.3, 0.4, 7.0]
    nrm = [1.0, 2.0]
  else
    q = Complex128[1.0, 0.3, 0.4, 0.5, 13.0]
    nrm = [1.0, 2.0, 3.0]
  end

  qL = zeros(Complex128, params.numDofPerNode, params.numNodesPerElement)
  qR = zeros(Complex128, params.numDofPerNode, params.numNodesPerElement)
  nrm_face = zeros(Tdim, params.numNodesPerFace)
  for i=1:params.numNodesPerElement
    qL[:, i] = q + (i-1)*0.1
    qR[:, i] = q + i*0.1
  end

  for i=1:params.numNodesPerFace
    nrm_face[:, i] = nrm + (i-1)*0.1
  end
  aux_vars = zeros(Complex128, 0, 0)

  for i=1:2  # run tests twice to make sure intermediate arrays are zeroed out
    resL = zeros(qL)
    resR = zeros(qR)
    jacLL = zeros(Complex128, params.numDofPerNode, params.numDofPerNode,
                  params.numNodesPerElement, params.numNodesPerElement)
    jacLR = copy(jacLL)
    jacRL = copy(jacLL)
    jacRR = copy(jacLL)

    #TODO: test all configurations
    iface = Interface(1, 2, 1, 1, 1)

    qL_dot = rand_realpart(size(qL))
    qR_dot = rand_realpart(size(qR))
  #  qL_dot = zeros(Complex128, params.numDofPerNode, params.numNodesPerElement)
  #  qR_dot = zeros(Complex128, params.numDofPerNode, params.numNodesPerElement)

  #  qL_dot[1, 1] = 1

    # complex step
    qL .+= pert*qL_dot
    qR .+= pert*qR_dot

    EulerEquationMod.calcEntropyPenaltyIntegral(params, sbpface, iface, kernel,
                                        qL, qR, aux_vars, nrm_face, resL, resR)
    valLc = imag(resL)/h
    valRc = imag(resR)/h

    qL .-= pert*qL_dot
    qR .-= pert*qR_dot


    EulerEquationMod.calcEntropyPenaltyIntegral_diff(params, sbpface, iface,
                  kernel, qL, qR, aux_vars, nrm_face, jacLL, jacLR, jacRL, jacRR)

    valL = zeros(qL)
    valR = zeros(qR)

    for q=1:params.numNodesPerElement
      for p=1:params.numNodesPerElement
        for i=1:params.numDofPerNode
          for j=1:params.numDofPerNode
            valL[i, p] += jacLL[i, j, p, q]*qL_dot[j, q]
            valL[i, p] += jacLR[i, j, p, q]*qR_dot[j, q]
            valR[i, p] += jacRL[i, j, p, q]*qL_dot[j, q]
            valR[i, p] += jacRR[i, j, p, q]*qR_dot[j, q]
          end
        end
      end
    end


    @test maximum(abs.(valL - valLc)) < 1e-12
    @test maximum(abs.(valR - valRc)) < 1e-12

    # test accumulation
    jacLL_orig = copy(jacLL)
    jacLR_orig = copy(jacLR)
    jacRL_orig = copy(jacRL)
    jacRR_orig = copy(jacRR)
    EulerEquationMod.calcEntropyPenaltyIntegral_diff(params, sbpface, iface,
                  kernel, qL, qR, aux_vars, nrm_face, jacLL, jacLR, jacRL, jacRR)


    @test maximum(abs.(jacLL - 2*jacLL_orig)) < 1e-13
    @test maximum(abs.(jacLR - 2*jacLR_orig)) < 1e-13
    @test maximum(abs.(jacRL - 2*jacRL_orig)) < 1e-13
    @test maximum(abs.(jacRR - 2*jacRR_orig)) < 1e-13
  end

  return nothing
end


"""
  Test reverse mode of calcSCurvilinear
"""
function test_SCurvilinear(params::AbstractParamType{Tdim}, sbp) where {Tdim}

  h = 1e-20
  pert = Complex128(0, h)

  numNodesPerElement = params.numNodesPerElement
  dxidx = rand_realpart(Tdim, Tdim, numNodesPerElement)
  Sx = zeros(Complex128, numNodesPerElement, numNodesPerElement, Tdim)

  dxidx_dot = rand_realpart(Tdim, Tdim, numNodesPerElement)
  dxidx_bar = zeros(Complex128, Tdim, Tdim, numNodesPerElement)
  Sx_bar = rand_realpart(numNodesPerElement, numNodesPerElement, Tdim)


  dxidx .+= pert*dxidx_dot
  Utils.calcSCurvilinear(sbp, dxidx, Sx)
  dxidx .-= pert*dxidx_dot
  val = sum( imag(Sx)/h .* Sx_bar)

  Utils.calcSCurvilinear_rev(sbp, dxidx, dxidx_bar, Sx, Sx_bar)
  val2 = sum(dxidx_bar .* dxidx_dot)

  @test abs(val - val2) < 1e-13

  # test accumulation
  dxidx_bar_orig = copy(dxidx_bar)
  Utils.calcSCurvilinear_rev(sbp, dxidx, dxidx_bar, Sx, Sx_bar)

  @test maximum(abs.(dxidx_bar - 2.*dxidx_bar_orig)) < 1e-13


  return nothing
end


"""
  Test the volume, face, and boundary terms of the explicitly computed Jacobian
  individually.  Small meshes only!
"""
function test_jac_assembly(mesh, sbp, eqn, _opts)

  opts = copy(_opts)  # don't modify the original

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]

  # test volume terms
  opts["addBoundaryIntegrals"] = false
  opts["addFaceIntegrals"] = false
  
  jac1 = SparseMatrixCSC(mesh, Float64, COLORING, LinearSolvers.getFaceType(mesh.sbpface))
  jac2 = zeros(mesh.numDof, mesh.numDof)
  assembler = Jacobian._AssembleElementData(jac2, mesh, sbp, eqn, opts)

  # compute jacobian via coloring
  opts["calc_jac_explicit"] = false
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  evalJacobian(mesh, sbp, eqn, opts, assembler)

  jac1d = full(jac1)

  @test isapprox( maximum(abs.(jac1d - jac2)), 0.0) atol=1e-13


  # test face integrals
  opts["addVolumeIntegrals"] = false
  opts["addFaceIntegrals"] = true
  fill!(jac1, 0.0)
  fill!(jac2, 0.0)
  # compute jacobian via coloring
  opts["calc_jac_explicit"] = false
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  evalJacobian(mesh, sbp, eqn, opts, assembler)

  jac1d = full(jac1)

  @test isapprox( maximum(abs.(jac1d - jac2)), 0.0) atol=1e-13

  # test boundary integral
  # test face integrals
  opts["addVolumeIntegrals"] = false
  opts["addFaceIntegrals"] = false
  opts["addBoundaryIntegrals"] = true
  fill!(jac1, 0.0)
  fill!(jac2, 0.0)
  # compute jacobian via coloring
  opts["calc_jac_explicit"] = false
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  evalJacobian(mesh, sbp, eqn, opts, assembler)

  jac1d = full(jac1)

  @test isapprox( maximum(abs.(jac1d - jac2)), 0.0) atol=1e-13


  return nothing
end


"""
  Test the entire jacobian assembly, for any type of jacobian matrix

  is_prealloc_exact: test that the jacobian preallocation is exact (Petsc only)
  set_prealloc: if true, set the preallocation of the jacobian to be tight
                for the explicitly computed jacobian
                if false, use the value currently in the dictionary
"""
function test_jac_general(mesh, sbp, eqn, opts; is_prealloc_exact=true, set_prealloc=true)

  #opts["addVolumeIntegrals"] = false
  #opts["addFaceIntegrals"] = false
  #opts["addBoundaryIntegrals"] = false
  #opts["addStabilization"] = false


  # use a spatially varying solution
  srand(1234)
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q .+= 0.01*rand(size(eqn.q))

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  println("getting flux function named ", opts["Flux_name"])
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]
  eqn.volume_flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Volume_flux_name"]]



  startSolutionExchange(mesh, sbp, eqn, opts)

  # get coloring linear operator
  opts["calc_jac_explicit"] = false
  pc1, lo1 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)

  # get explicitly computed linear operator
  opts["calc_jac_explicit"] = true
  val_orig = opts["preallocate_jacobian_coloring"]
  if set_prealloc
    opts["preallocate_jacobian_coloring"] = false
  end
  pc2, lo2 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
  opts["preallocate_jacobian_coloring"] = val_orig

  jac1 = getBaseLO(lo1).A
  jac2 = getBaseLO(lo2).A

  assembler = Jacobian._AssembleElementData(getBaseLO(lo2).A, mesh, sbp, eqn, opts)

  # computing coloring jacobian
  opts["calc_jac_explicit"] = false
  println("calculating regular jacobian"); flush(STDOUT)
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true

  evalJacobian(mesh, sbp, eqn, opts, assembler)

  assembly_begin(jac1, MAT_FINAL_ASSEMBLY)
  assembly_begin(jac2, MAT_FINAL_ASSEMBLY)

  assembly_end(jac1, MAT_FINAL_ASSEMBLY)
  assembly_end(jac2, MAT_FINAL_ASSEMBLY)


  # multiply against a random vector to make sure the jacobian is
  # the same
  for i=1:10
    x = rand(PetscScalar, mesh.numDof)
    b1 = zeros(PetscScalar, mesh.numDof)
    b2 = zeros(PetscScalar, mesh.numDof)
    b3 = zeros(PetscScalar, mesh.numDof)

    t = 0.0
    applyLinearOperator(lo1, mesh, sbp, eqn, opts, ctx_residual, t, x, b1)
    applyLinearOperator(lo2, mesh, sbp, eqn, opts, ctx_residual, t, x, b2)
    evaldRdqProduct(mesh, sbp, eqn, opts, x, b3)

    println("b1 - b2 = ", maximum(abs.((b1 - b2))))
    println("b2 - b3 = ", maximum(abs.((b2 - b3))))
    println("b1 - b3 = ", maximum(abs.((b1 - b3))))
#    println("b2[1] = ", b2[1])
#    println("b3[1] = ", b3[1])
    @test isapprox( maximum(abs.((b1 - b2))), 0.0) atol=1e-11
    @test isapprox( maximum(abs.((b2 - b3))), 0.0) atol=1e-11
  end

  A = getBaseLO(lo2).A
  if typeof(A) <: PetscMat
    matinfo = MatGetInfo(A, PETSc2.MAT_LOCAL)
    if is_prealloc_exact
      @test ( matinfo.nz_unneeded )== 0
    else
      @test  matinfo.nz_unneeded  > 0
    end
  end

  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)

  return nothing
end

function test_jac_homotopy(mesh, sbp, eqn, opts)

  println("\nTesting homotopy jacobian")

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff["RoeFlux"]
  opts["homotopy_addBoundaryIntegrals"] = true
#=
  res1 = zeros(eqn.res)
  res2 = zeros(eqn.res)
  println("\ncomputing regular homotopy dissipation")
#  EulerEquationMod.calcHomotopyDiss(mesh, sbp, eqn, opts, res1)
  println("\ncomputing new homotopy dissipation")
  h = 1e-20
  pert = Complex128(0, h)
  eqn.q[1] += pert
  EulerEquationMod.calcHomotopyDiss(mesh, sbp, eqn, opts, res2)
  eqn.q[1] -= pert
=#
#=
  println("diffnorm = ", vecnorm(res1 - res2))
  println("res1 = \n", res1)
  println("res2 = \n", res2)
  println("diff = \n", res1 - res2)
  @assert vecnorm(res1 - res2) < 1e-13
=#

  function _evalHomotopy(mesh, sbp, eqn, opts, t)
    NonlinearSolvers._homotopyPhysics(hdata, mesh, sbp, eqn, opts, t)
  end


  startSolutionExchange(mesh, sbp, eqn, opts)

  Tsol = eltype(eqn.q)
  Tjac = real(eltype(eqn.res))
  hdata = NonlinearSolvers.HomotopyData{Tsol, Tjac}(mesh, sbp, eqn, opts,
                                            evalResidual, evalHomotopy)
  hdata.lambda = 0
  println("constructing first operator")
  opts["calc_jac_explicit"] = false
  pc1, lo1 = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)
  lo1.hdata = hdata

  println("constructing second operator")
  opts["calc_jac_explicit"] = true
  pc2, lo2 = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)
  lo2.hdata = hdata

  jac1 = getBaseLO(lo1).A
  jac2 = getBaseLO(lo2).A

  ctx_residual = (_evalHomotopy,)
  println("\nevaluating jacobians")

  opts["calc_jac_explicit"] = false
  println("calculating regular jacobian"); flush(STDOUT)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
 
  # was lambda = lo2.hdata.lambda
  calcLinearOperator(lo2, mesh, sbp, eqn, opts, ctx_residual, 0.0)

  assembly_begin(jac1, MAT_FINAL_ASSEMBLY)
  assembly_begin(jac2, MAT_FINAL_ASSEMBLY)

  # multiply against a random vector to make sure the jacobian is
  # the same
  for i=1:10
    x = rand(PetscScalar, mesh.numDof)
    b1 = zeros(PetscScalar, mesh.numDof)
    b2 = zeros(PetscScalar, mesh.numDof)

    t = 0.0
    applyLinearOperator(lo1, mesh, sbp, eqn, opts, ctx_residual, t, x, b1)
    applyLinearOperator(lo2, mesh, sbp, eqn, opts, ctx_residual, t, x, b2)

    @test isapprox( norm(b1 - b2), 0.0) atol=1e-12
  end

  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)


  # test enabling/disabling terms
  opts_orig = copy(opts)

  NonlinearSolvers.disableTerms(hdata, opts)
  NonlinearSolvers.enableTerms(hdata, opts)

  @test length(keys(opts)) == length(keys(opts_orig))
  for (key, val) in opts_orig
    @test opts[key] == val
  end

  return nothing
end

function test_diagjac(mesh, sbp, eqn, opts)
# compare the diagonal Jacobian with the diagonal of an explicitly computed
# jacobian

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]


  startSolutionExchange(mesh, sbp, eqn, opts)

  #----------------------------------------------------------------------------
  # construct the matrices/assemblers
  nblocks = mesh.numEl
  blocksize = mesh.numDofPerNode*mesh.numNodesPerElement


  opts["calc_jac_explicit"] = true

  jac1 = NonlinearSolvers.DiagJac(Complex128, blocksize, mesh.numEl)
  assem1 = NonlinearSolvers.AssembleDiagJacData(mesh, sbp, eqn, opts, jac1)


  opts["preallocate_jacobian_coloring"] = false
  pc2, lo2 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
  jac2 = getBaseLO(lo2).A

  assem2 = Jacobian._AssembleElementData(getBaseLO(lo2).A, mesh, sbp, eqn, opts)

 
  evalJacobian(mesh, sbp, eqn, opts, assem1)
  evalJacobian(mesh, sbp, eqn, opts, assem2)

  # compare the matrices themselves
  # zero out off block-diagonal part of jac2

  jac2_full = full(jac2)
#  println("jac2_full = \n", jac2_full)
#  println("jac1 = \n", jac1)

  # because of the format of DiagJac, it isn't generally true that the blocks
  # are in the same order as in the real Jacobian, but it works for the simple
  # order where each node of each element is numbered 1:n
  for block=1:nblocks
    idx = ((block - 1)*blocksize + 1):(block*blocksize)
    @test isapprox( norm(jac2_full[idx, idx] - jac1.A[:, :, block]), 0.0) atol=1e-13
  end

  # zero out the off diagonal parts of jac2_full
  tmp = zeros(blocksize, blocksize, nblocks)
  for block=1:nblocks
    idx = ((block - 1)*blocksize + 1):(block*blocksize)
    tmp[:, :, block] = jac2_full[idx, idx]
  end

  fill!(jac2_full, 0.0)

  for block=1:nblocks
    idx = ((block - 1)*blocksize + 1):(block*blocksize)
    jac2_full[idx, idx] = tmp[:, :, block]
  end


  # test multiplication
  for i=1:10
    x = rand(mesh.numDof)
    b = zeros(Complex128, mesh.numDof)
    
    NonlinearSolvers.diagMatVec(jac1, mesh, x, b)
    b2 = jac2_full * x

    @test isapprox( norm(b2 - b), 0.0) atol=1e-12
  end


  return nothing
end

function test_strongdiagjac(mesh, sbp, eqn, _opts)
# for a uniform flow (and BCs), the strong form and weak form should be
# equivalent

  opts = copy(_opts)  # dont modify the original

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICRho1E2U3"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]

  startSolutionExchange(mesh, sbp, eqn, opts)

  #----------------------------------------------------------------------------
  # construct the matrices/assemblers
  nblocks = mesh.numEl
  blocksize = mesh.numDofPerNode*mesh.numNodesPerElement

  opts["addBoundaryIntegrals"] = false
  opts["addStabilization"] = false
  opts["addFaceIntegrals"] = false
  opts["Q_transpose"] = false

  # complex step jacobian (coloring)
  opts["calc_jac_explicit"] = false
  opts["preallocate_jacobian_coloring"] = true
  pc1, lo1 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
  jac1 = getBaseLO(lo1).A
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # diagonal jacobian
  opts["calc_jac_explicit"] = true
  jac2 = NonlinearSolvers.DiagJac(Complex128, blocksize, nblocks)
  assem2 = NonlinearSolvers.AssembleDiagJacData(mesh, sbp, eqn, opts, jac2) 
  evalJacobianStrong(mesh, sbp, eqn, opts, assem2)

  # test multiplication
  for i=1:10
    x = rand(mesh.numDof)
    b2 = zeros(Complex128, mesh.numDof)
    
    b = jac1 * x
    NonlinearSolvers.diagMatVec(jac2, mesh, x, b2)

    @test isapprox( norm(b2 - b), 0.0) atol=1e-12
  end


  return nothing
end

"""
  Test psi^T dRdm product
"""
function test_revm_product(mesh, sbp, eqn, opts; freeze=false)

  h = 1e-20
  pert = Complex128(0, h)

  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q .+= 0.01.*rand(size(eqn.q))  # add a little noise, to make jump across
                                     # interfaces non-zero
  
  # to test freezing the viscosity, have two different qs and evaluate the
  # Jacobian using the frozen viscosity while a different q in eqn.q
  if freeze
    q_frozen = copy(eqn.q)
    EulerEquationMod.freezeViscosity(mesh, sbp, eqn, opts)
    eqn.q .+= 0.01*rand(size(eqn.q))
  end


  # fields: dxidx, jac, nrm_bndry, nrm_face, coords_bndry

  res_bar = rand_realpart(mesh.numDof)

  dxidx_dot       = rand_realpart(size(mesh.dxidx))
  jac_dot         = rand_realpart(size(mesh.jac))
  nrm_bndry_dot   = rand_realpart(size(mesh.nrm_bndry))
  nrm_face_dot    = rand_realpart(size(mesh.nrm_face_bar))
  coords_bndry_dot = rand_realpart(size(mesh.coords_bndry))


  zeroBarArrays(mesh)

  mesh.dxidx        .+= pert*dxidx_dot
  mesh.jac          .+= pert*jac_dot
  mesh.nrm_bndry    .+= pert*nrm_bndry_dot
  mesh.nrm_face     .+= pert*nrm_face_dot
  mesh.coords_bndry .+= pert*coords_bndry_dot
  updateMetricDependents(mesh, sbp, eqn, opts)

  fill!(eqn.res, 0)
  evalResidual(mesh, sbp, eqn, opts)
  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  val = sum(imag(eqn.res_vec)/h .* res_bar)

  mesh.dxidx        .-= pert*dxidx_dot
  mesh.jac          .-= pert*jac_dot
  mesh.nrm_bndry    .-= pert*nrm_bndry_dot
  mesh.nrm_face     .-= pert*nrm_face_dot
  mesh.coords_bndry .-= pert*coords_bndry_dot
  updateMetricDependents(mesh, sbp, eqn, opts)

  evalResidual_revm(mesh, sbp, eqn, opts, res_bar)
  val2 = sum(mesh.dxidx_bar .* dxidx_dot)              +
         sum(mesh.jac_bar .* jac_dot)                  +
         sum(mesh.nrm_bndry_bar .* nrm_bndry_dot)      +
         sum(mesh.nrm_face_bar .* nrm_face_dot)        +
         sum(mesh.coords_bndry_bar .* coords_bndry_dot)

  println("val = ", real(val))
  println("val2 = ", real(val2))
  println("max dxidx_bar = ", maximum(abs.(mesh.dxidx_bar)))
  println("max jac_bar = ", maximum(abs.(mesh.jac_bar)))
  println("max nrm_bndry_bar = ", maximum(abs.(mesh.nrm_bndry_bar)))
  println("max nrm_face_bar = ", maximum(abs.(mesh.nrm_face_bar)))
  println("max coords_bndry_bar = ", maximum(abs.(mesh.coords_bndry_bar)))
  @test abs(val - val2) < 1e-12

  # test accumulation behavior


  dxidx_bar_orig        = copy(mesh.dxidx_bar)
  jac_bar_orig          = copy(mesh.jac_bar)
  nrm_bndry_bar_orig    = copy(mesh.nrm_bndry_bar)
  nrm_face_bar_orig     = copy(mesh.nrm_face_bar)
  coords_bndry_bar_orig = copy(mesh.coords_bndry_bar)

  evalResidual_revm(mesh, sbp, eqn, opts, res_bar)

  @test maximum(abs.(mesh.dxidx_bar - 2*dxidx_bar_orig)) < 1e-13
  @test maximum(abs.(mesh.jac_bar - 2*jac_bar_orig)) < 1e-13
  @test maximum(abs.(mesh.nrm_bndry_bar - 2*nrm_bndry_bar_orig)) < 1e-13
  @test maximum(abs.(mesh.nrm_face_bar - 2*nrm_face_bar_orig)) < 1e-13
  @test maximum(abs.(mesh.coords_bndry_bar - 2*coords_bndry_bar_orig)) < 1e-13

#=
  if opts["face_integral_type"] == 1
    println("testing type 1 face integrals")
    zeroBarArrays(mesh)
    fill!(eqn.res, 0)

    mesh.nrm_face .+= pert*nrm_face_dot
    EulerEquationMod.calcFaceIntegral_nopre(mesh, sbp, eqn, opts, eqn.flux_func, mesh.interfaces)
    mesh.nrm_face .-= pert*nrm_face_dot
    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    val1 = sum(imag(eqn.res_vec)/h .* res_bar)

    EulerEquationMod.calcFaceIntegral_nopre_revm(mesh, sbp, eqn, opts, eqn.flux_func_bar, mesh.interfaces)
    val2 = sum(mesh.nrm_face_bar .* nrm_face_dot)

    println("val = ", real(val))
    println("val2 = ", real(val2))
    @test abs(val - val2) < 1e-13
  end
=#


  if opts["volume_integral_type"] == 2
    # test SplitFormLinear and SplitFormCurvilinear

    # linear
    zeroBarArrays(mesh)
    fill!(eqn.res, 0)

    mesh.dxidx .+= pert*dxidx_dot
    EulerEquationMod.calcVolumeIntegralsSplitFormLinear(mesh, sbp, eqn,
                                                   opts, eqn.volume_flux_func)

    mesh.dxidx .-= pert*dxidx_dot
    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    val = sum(imag(eqn.res_vec)/h .* res_bar)

    EulerEquationMod.calcVolumeIntegralsSplitFormLinear_revm(
        mesh, sbp, eqn, opts, eqn.volume_flux_func, eqn.volume_flux_func_revm)
    val2 = sum(mesh.dxidx_bar .* dxidx_dot)

    @test abs(val - val2) < 1e-13

    # test accumulation
    dxidx_bar_orig = copy(mesh.dxidx_bar)
    EulerEquationMod.calcVolumeIntegralsSplitFormLinear_revm(
        mesh, sbp, eqn, opts, eqn.volume_flux_func, eqn.volume_flux_func_revm)

    @test maximum(abs.(mesh.dxidx_bar - 2*dxidx_bar_orig)) < 1e-13


    # curvilinear
    zeroBarArrays(mesh)
    fill!(eqn.res, 0)

    mesh.dxidx .+= pert*dxidx_dot
    EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear(mesh, sbp, eqn,
                                                   opts, eqn.volume_flux_func)
    mesh.dxidx .-= pert*dxidx_dot
    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    val = sum(imag(eqn.res_vec)/h .* res_bar)

    EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear_revm(
        mesh, sbp, eqn, opts, eqn.volume_flux_func, eqn.volume_flux_func_revm)

    val2 = sum(mesh.dxidx_bar .* dxidx_dot)

    @test abs(val - val2) < 1e-13

    # test accumulation
    dxidx_bar_orig = copy(mesh.dxidx_bar)
    EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear_revm(
        mesh, sbp, eqn, opts, eqn.volume_flux_func, eqn.volume_flux_func_revm)

    @test maximum(abs.(mesh.dxidx_bar - 2*dxidx_bar_orig)) < 1e-13

  end

  if freeze
    EulerEquationMod.thawViscosity(mesh, sbp, eqn, opts)
  end




  return nothing
end


"""
  Tests the reverse mode product by perturbing mesh.vert_coords and
  recalculating the metrics from that
"""
function test_revm_product2(mesh, sbp, eqn, opts)

  h = 1e-20
  pert = Complex128(0, h)

  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q .+= 0.01.*rand(size(eqn.q))  # add a little noise, to make jump across
                                     # interfaces non-zero
  
  vert_coords_dot = rand_realpart(size(mesh.vert_coords))

  res_bar = rand_realpart(mesh.numDof)

  # complex step
  mesh.vert_coords .+= pert*vert_coords_dot
  recalcCoordinatesAndMetrics(mesh, sbp, opts)
  evalResidual(mesh, sbp, eqn, opts)
  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  val1 = sum(imag(eqn.res_vec)/h .* res_bar)

  # undo perturbation
  for i=1:length(mesh.vert_coords)
    mesh.vert_coords[i] = real(mesh.vert_coords[i])
  end
  recalcCoordinatesAndMetrics(mesh, sbp, opts)

  # reverse mode

  zeroBarArrays(mesh)
  evalResidual_revm(mesh, sbp, eqn, opts, res_bar)
  getAllCoordinatesAndMetrics_rev(mesh, sbp, opts)

  val2 = sum(mesh.vert_coords_bar .* vert_coords_dot)

  println("val1 = ", val1)
  println("val2 = ", val2)

  @test abs(val1 - val2) < max( min(abs(val1), abs(val2))*1e-12, 1e-12)

  return nothing
end


"""
  Test the psi^T dR/dq product
"""
function test_revq_product(mesh, sbp, eqn, opts; freeze=false)

  srand(1234)
  h = 1e-20
  pert = Complex128(0, h)

  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
#  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q_vec .+= 0.01.*rand(size(eqn.q_vec))  # add a little noise, to make jump across
                                     # interfaces non-zero

  # to test freezing the viscosity, have two different qs and evaluate the
  # Jacobian using the frozen viscosity while a different q in eqn.q
  if freeze
    q_frozen = copy(eqn.q)
    EulerEquationMod.freezeViscosity(mesh, sbp, eqn, opts)
    eqn.q .+= 0.01*rand(size(eqn.q))
  end


  # fields: dxidx, jac, nrm_bndry, nrm_face, coords_bndry

  res_vec_bar = rand_realpart(mesh.numDof)
  q_vec_dot = rand_realpart(mesh.numDof)
#  q_vec_dot = zeros(mesh.numDof)
#  q_vec_dot[1] = 1
  q_vec_bar = zeros(Complex128, mesh.numDof)

  fill!(eqn.res, 0)
  eqn.q_vec .+= pert*q_vec_dot
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  evalResidual(mesh, sbp, eqn, opts)
  eqn.q_vec .-= pert*q_vec_dot
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  val = sum(res_vec_bar .* imag(eqn.res_vec)/h)

  evalResidual_revq(mesh, sbp, eqn, opts, res_vec_bar, q_vec_bar)
  val2 = sum(q_vec_bar .* q_vec_dot)

  println("val = ", val)
  println("val2 = ", val2)
  @test abs(val - val2) < 5e-12

  #TODO: test accumulation behavior
#=
  if opts["volume_integral_type"] == 2
    # test SplitFormLinear and SplitFormCurvilinear

    # linear
    fill!(eqn.res, 0)
    fill!(eqn.q_bar, 0)

    q_dot = rand_realpart(size(eqn.q))
    q_bar = zeros(Complex128, size(eqn.q))
    res_bar = rand_realpart(size(eqn.res))

    eqn.q .+= pert*q_dot
    EulerEquationMod.calcVolumeIntegralsSplitFormLinear(mesh, sbp, eqn,
                                                   opts, eqn.volume_flux_func)
    eqn.q .-= pert*q_dot
    val = sum(res_bar .* imag(eqn.res)/h)

    copy!(eqn.res_bar, res_bar)
    EulerEquationMod.calcVolumeIntegralsSplitFormLinear_revq(
        mesh, sbp, eqn, opts, eqn.volume_flux_func, eqn.volume_flux_func_revq)

    val2 = sum( q_dot .* eqn.q_bar)

    println("val = ", val)
    println("val2 = ", val2)
    @test abs(val - val2) < 1e-13

    # test accumulation
    q_bar_orig = copy(eqn.q_bar)
    EulerEquationMod.calcVolumeIntegralsSplitFormLinear_revq(
        mesh, sbp, eqn, opts, eqn.volume_flux_func, eqn.volume_flux_func_revq)
    
    @test maximum(abs.(eqn.q_bar - 2.*q_bar_orig)) < 1e-13

    # curvilinear

    fill!(eqn.res, 0)
    fill!(eqn.q_bar, 0)

    eqn.q .+= pert*q_dot
    EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear(mesh, sbp, eqn,
                                                   opts, eqn.volume_flux_func)
    eqn.q .-= pert*q_dot
    val = sum(res_bar .* imag(eqn.res)/h)

    copy!(eqn.res_bar, res_bar)
    EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear_revq(
        mesh, sbp, eqn, opts, eqn.volume_flux_func, eqn.volume_flux_func_revq)

    val2 = sum(q_dot .* eqn.q_bar)

    println("val = ", val)
    println("val2 = ", val2)
    @test abs(val - val2) < 1e-13

    # test accumulation
    q_bar_orig = copy(eqn.q_bar)
    EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear_revq(
        mesh, sbp, eqn, opts, eqn.volume_flux_func, eqn.volume_flux_func_revq)
    
    @test maximum(abs.(eqn.q_bar - 2.*q_bar_orig)) < 1e-13


  end
=#

  if freeze
    EulerEquationMod.thawViscosity(mesh, sbp, eqn, opts)
  end


  return nothing
end


"""
  Test NewtonBDiagPC
"""
function test_BDiagPC(mesh, sbp, eqn, opts)

  opts["calc_jac_explicit"] = true
  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q .+= 0.01*rand(size(eqn.q))

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  println("getting flux function named ", opts["Flux_name"])
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]
  eqn.volume_flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Volume_flux_name"]]


  bs = mesh.numDofPerNode*mesh.numNodesPerElement
  diag_jac = NonlinearSolvers.DiagJac(Float64, bs, mesh.numEl)
  assem = NonlinearSolvers.AssembleDiagJacData(mesh, sbp, eqn, opts, diag_jac)

  evalJacobian(mesh, sbp, eqn, opts, assem)

  pc = NonlinearSolvers.NewtonBDiagPC(mesh, sbp, eqn, opts)
  calcPC(pc, mesh, sbp, eqn, opts, (evalResidual,), 0.0)

  @testset "NewtonBDiagPC" begin
    # test inverse: use applyPC to multiply by A^-1, then use diagMatVec to
    # undo it
    x = rand(mesh.numDof)
    x_orig = copy(x)
    b = zeros(mesh.numDof)
    applyPC(pc, mesh, sbp, eqn, opts, 0.0, x, b)
    fill!(x, 0)
    NonlinearSolvers.diagMatVec(diag_jac, mesh, b, x)

    @test maximum(abs.(x - x_orig)) < 1e-12

    # test factored multiplication
    @test pc.is_factored
    b = zeros(mesh.numDof)
    b2 = zeros(mesh.numDof)
    NonlinearSolvers.diagMatVec(diag_jac, mesh, x, b)
    NonlinearSolvers.applyBDiagPCInv(pc, mesh, sbp, eqn, opts, x, b2)
  
    @test maximum(abs.(b2 - b)) < 1e-12

    # test transpose inverse
    x = rand(mesh.numDof)
    x_orig = copy(x)
    b = zeros(mesh.numDof)
    applyPCTranspose(pc, mesh, sbp, eqn, opts, 0.0, x, b)
    fill!(x, 0)
    NonlinearSolvers.diagMatVec(diag_jac, mesh, b, x, trans=true)

    @test maximum(abs.(x - x_orig)) < 1e-12

    # test factored multiplication
    @test pc.is_factored
    b = zeros(mesh.numDof)
    b2 = zeros(mesh.numDof)
    NonlinearSolvers.diagMatVec(diag_jac, mesh, x, b, trans=true)
    NonlinearSolvers.applyBDiagPCInv(pc, mesh, sbp, eqn, opts, x, b2, trans=true)
 
  end


  return nothing
end


function test_BJacobiPC(mesh, sbp, eqn, opts)

  opts["calc_jac_explicit"] = true

  const fac = 5000  # factor for diagonal
  # functions for test
  function eval_jac(mesh, sbp, eqn, opts, assem)
    evalJacobian(mesh, sbp, eqn, opts, assem)

    bs = mesh.numDofPerNode*mesh.numNodesPerElement
    res_jac = zeros(mesh.numDofPerNode, mesh.numDofPerNode,
                    mesh.numNodesPerElement, mesh.numNodesPerElement)
    for i=1:mesh.numNodesPerElement
      for j=1:mesh.numDofPerNode
        res_jac[j, j, i, i] = fac
      end
    end

    for i=1:mesh.numEl
      assembleElement(assem, mesh, i, res_jac)
    end

    return nothing
  end

  function eval_jacvec(mesh, sbp, eqn, opts, x, b)

    evaldRdqProduct(mesh, sbp, eqn, opts, x, b)
    for i=1:length(b)
      b[i] += fac*x[i]
    end

    return nothing
  end

  function eval_jacTvec(mesh, sbp, eqn, opts, x, b)

    evalResidual_revq(mesh, sbp, eqn, opts, x, b)
    for i=1:length(b)
      b[i] += fac*x[i]
    end

    return nothing
  end
 
  
  #TODO: pull this out into a separate function
  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q .+= 0.01*rand(size(eqn.q))

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  println("getting flux function named ", opts["Flux_name"])
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]
  eqn.volume_flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Volume_flux_name"]]

  
  pc = NonlinearSolvers.NewtonBJacobiPC(mesh, sbp, eqn, opts, itermax=10,
                                        res_tol=1e-8, verbose=true)
  NonlinearSolvers.setEvalJacobian(pc, eval_jac, eval_jacvec, eval_jacTvec)

  ctx_residual = (evalResidual,)
  calcPC(pc, mesh, sbp, eqn, opts, ctx_residual, 0.0)

  # test that it solves diagonally dominant problem
  b = ones(Float64, mesh.numDof)
  x = zeros(Float64, mesh.numDof)
  applyPC(pc, mesh, sbp, eqn, opts, 0.0, b, x)

  b2 = zeros(Float64, mesh.numDof)
  eval_jacvec(mesh, sbp, eqn, opts, x, b2)

  println("maxdiff = ", maximum(abs.(b - b2)))
  @test maximum(abs.(b - b2)) < 2e-8


  # test transpose solve
  fill!(x, 0)
  applyPCTranspose(pc, mesh, sbp, eqn, opts, 0.0, b, x)

  b2 = zeros(Float64, mesh.numDof)
  eval_jacTvec(mesh, sbp, eqn, opts, x, b2)

  println("maxdiff = ", maximum(abs.(b - b2)))
  @test maximum(abs.(b - b2)) < 2e-8



  return nothing
end

"""
  Check if two floating point values are approximately equal, using a
  either an absolule or relative tolerance
  
  **Inputs**

   * val1
   * val2
   * tol: the tolerance
"""
function checkvals(val1, val2, tol)

  pass = abs(val1 - val2) < tol ||
         abs(val1 - val2) < tol*min(abs(val1), abs(val2))

  if !pass
    println("val1 = ", val1)
    println("val2 = ", val2)
    println("diff = ", val1 - val2)
  end

  return pass
end


function test_sbp_cartesian(mesh::AbstractMesh, sbp, eqn, opts)

  h = 1e-20
  pert = Complex128(0, h)

  # use a spatially varying solution
  #icfunc = EulerEquationMod.ICDict["ICExp"]
  #icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q .+= 0.01*rand(size(eqn.q))

  #q_dot = rand_realpart(size)

  Dx = zeros(Complex128, mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  wxi = zeros(Complex128, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  res = zeros(Complex128, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)

  t2_dot = zeros(Complex128, mesh.numDofPerNode, mesh.numDofPerNode, mesh.dim,
                 mesh.numNodesPerElement, mesh.numNodesPerElement)
  t3_dot = zeros(Complex128, mesh.numDofPerNode, mesh.numDofPerNode,
                 mesh.numNodesPerElement,
                         mesh.numNodesPerElement)
  t4_dot = zeros(Complex128, mesh.numDofPerNode, mesh.numDofPerNode, mesh.dim,
                 mesh.numNodesPerElement, mesh.numNodesPerElement)

  q2 = zeros(Complex128, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  res2 = zeros(Complex128, mesh.numDofPerNode, mesh.numNodesPerElement)

  for i=1:mesh.numEl

    q_dot = rand_realpart(mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerElement,
                           mesh.numNodesPerElement)
    q_i = eqn.q[:, :, i]
    dxidx_i = sview(mesh.dxidx, :, :, :, i)
    jac_i = sview(mesh.jac, :, i)

    # use differentiated version (4D -> 5D)
    EulerEquationMod.calcDx(sbp, dxidx_i, jac_i, Dx)
    EulerEquationMod.applyOperatorJac(Dx, q_dot, t2_dot)

    for q=1:mesh.numNodesPerElement
#      println("q = ", q)
      for j=1:mesh.numDofPerNode
#        println("  j = ", j)
        q_i .+= pert*q_dot[:, j, :, q]
        fill!(res, 0)
#        println("q = ", q_i)
#        println("dxidx = ", dxidx_i)
#        println("jac = ", jac_i)
        EulerEquationMod.applyDx(sbp, q_i, dxidx_i, jac_i, wxi, res)
        q_i .-= pert*q_dot[:, j, :, q]
        for d=1:mesh.dim
#          println("    d = ", d)
#          println("     t2_dot = \n", real(t2_dot[:, j, d, :, q]))
#          println("     res_dot = \n", res[:, :, d])
          @test maximum(abs.(t2_dot[:, j, d, :, q] - imag(res[:, :, d])./h)) < 5e-12
        end
      end
    end

    
    # test the second method (5D -> 4D)
    q_dot2 = rand_realpart(mesh.numDofPerNode, mesh.numDofPerNode, mesh.dim,
                           mesh.numNodesPerElement, mesh.numNodesPerElement)
    EulerEquationMod.applyOperatorJac(Dx, q_dot2, t3_dot)

    for d=1:mesh.dim
      q2_d = sview(q2, :, :, d)
      copy!(q2, q_i)
    end

    for q=1:mesh.numNodesPerElement
      for j=1:mesh.numDofPerNode
        for d=1:mesh.dim
          q2[:, :, d] .+= pert*q_dot2[:, j, d, :, q]
        end

        fill!(res2, 0)
        EulerEquationMod.applyDx(sbp, q2, dxidx_i, jac_i, wxi, res2)

        for d=1:mesh.dim
          q2[:, :, d] .-= pert*q_dot2[:, j, d, :, q]
        end

        @test maximum(abs.(t3_dot[:, j, :, q] - imag(res2)/h)) < 5e-12
      end
    end

    # test the 3rd methods: (3D -> 5D)
    # test consistency with 4D method
    q_dot3 = rand_realpart(mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerElement)
    q_dot3b = zeros(mesh.numDofPerNode, mesh.numDofPerNode,
                    mesh.numNodesPerElement, mesh.numNodesPerElement)

    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        for p=1:mesh.numDofPerNode
          q_dot3b[p, k, j, j] = q_dot3[p, k, j]
        end
      end
    end
    EulerEquationMod.applyOperatorJac(Dx, q_dot3, t4_dot)
    EulerEquationMod.applyOperatorJac(Dx, q_dot3b, t2_dot)

    @test maximum(abs.(t4_dot - t2_dot)) < 5e-12
    
  end

  return nothing
end

function test_sbp_cartesian_revm(mesh, sbp, eqn, opts)

  Tres = Complex128
  w_dot = rand_realpart(mesh.numDofPerNode, mesh.numNodesPerElement)
  wx_bar = rand_realpart(mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  dxidx_dot = rand_realpart(mesh.dim, mesh.dim, mesh.numNodesPerElement)
  jac_dot = rand_realpart(mesh.numNodesPerElement)
  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  wx1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  wx2 = copy(wx1); wx3 = copy(wx1); wx4 = copy(wx1)

  w_bar1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  w_bar2 = copy(w_bar1); w_bar3 = copy(w_bar1); w_bar4 = copy(w_bar1)

  dxidx_bar1 = zeros(Tres, mesh.dim, mesh.dim, mesh.numNodesPerElement)
  jac_bar1 = zeros(mesh.numNodesPerElement)
  dxidx_bar2 = copy(dxidx_bar1); dxidx_bar3 = copy(dxidx_bar1); dxidx_bar4 = copy(dxidx_bar1)
  jac_bar2 = copy(jac_bar1); jac_bar3 = copy(jac_bar1)

  w_i = rand_realpart((mesh.numDofPerNode, mesh.numNodesPerElement))
  op = SummationByParts.Subtract()

  h = 1e-20
  pert = Complex128(0, h)
  for i=1: mesh.numEl

    #w_i = sview(eqn.q, :, :, i)
    dxidx = sview(mesh.dxidx, :, :, :, i)
    jac = sview(mesh.jac, :, i)

    # perturb all the functions at once
    w_i .+= pert*w_dot
    dxidx .+= pert*dxidx_dot
    jac .+= pert*jac_dot
    EulerEquationMod.applyQxTransposed(sbp, w_i, dxidx, work, wx1, op)
    EulerEquationMod.applyDx(sbp, w_i, dxidx, jac, work, wx2, op)
    EulerEquationMod.applyDxTransposed(sbp, w_i, dxidx, jac, work, wx3, op)
    EulerEquationMod.applyQx(sbp, w_i, dxidx, work, wx4, op)
    w_i .-= pert*w_dot
    dxidx .-= pert*dxidx_dot
    jac .-= pert*jac_dot

    # test applyQxTransposed, first method
    EulerEquationMod.applyQxTransposed_revm(sbp, w_i, w_bar1, dxidx, dxidx_bar1, work, wx_bar, op)
    wx_dot = imag(wx1)./h
    val1 = sum(wx_dot .* wx_bar)
    val2 = sum(dxidx_bar1 .* dxidx_dot) + sum(w_bar1 .* w_dot)

    @test checkvals(val1, val2, 1e-12)

    # applyDx, first method
    EulerEquationMod.applyDx_revm(sbp, w_i, w_bar2, dxidx, dxidx_bar2, jac,
                                  jac_bar2, work, wx_bar, op)
    wx_dot = imag(wx2)./h
    val1 = sum(wx_dot .* wx_bar)
    val2 = sum(dxidx_bar2 .* dxidx_dot) + sum(jac_bar2 .* jac_dot) + sum(w_bar2 .* w_dot)

    @test checkvals(val1, val2, 1e-12)

    # applyDxTransposed, first method
    EulerEquationMod.applyDxTransposed_revm(sbp, w_i, w_bar3, dxidx, dxidx_bar3,
                                            jac, jac_bar3, work, wx_bar, op)
    wx_dot = imag(wx3)./h
    val1 = sum(wx_dot .* wx_bar)
    val2 = sum(dxidx_bar3 .* dxidx_dot) + sum(jac_bar3 .* jac_dot) + sum(w_bar3 .* w_dot)

    @test checkvals(val1, val2, 1e-12)

    
    # applyQx, first method
    EulerEquationMod.applyQx_revm(sbp, w_i, w_bar4, dxidx, dxidx_bar4, work,
                                  wx_bar, op)
    wx_dot = imag(wx4)./h
    val1 = sum(wx_dot .* wx_bar)
    val2 = sum(dxidx_bar4 .* dxidx_dot) + sum(w_bar4 .* w_dot)

    @test checkvals(val1, val2, 1e-12)
    
  end

  test_sbp_cartesian2_revm(mesh, sbp, eqn, opts)

  return nothing
end


function test_sbp_cartesian2_revm(mesh, sbp, eqn, opts)
  # test the second method (3D -> 2D) of all the sbp cartesian methods

  println("\ntesting second methods")
  Tres = Complex128
  w_dot = rand_realpart(mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  wx_bar = rand_realpart(mesh.numDofPerNode, mesh.numNodesPerElement)
  dxidx_dot = rand_realpart(mesh.dim, mesh.dim, mesh.numNodesPerElement)
  jac_dot = rand_realpart(mesh.numNodesPerElement)
  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  wx1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  wx2 = copy(wx1); wx3 = copy(wx1); wx4 = copy(wx1)

  w_bar1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  w_bar2 = copy(w_bar1); w_bar3 = copy(w_bar1); w_bar4 = copy(w_bar1)

  dxidx_bar1 = zeros(Tres, mesh.dim, mesh.dim, mesh.numNodesPerElement)
  jac_bar1 = zeros(mesh.numNodesPerElement)
  dxidx_bar2 = copy(dxidx_bar1); dxidx_bar3 = copy(dxidx_bar1); dxidx_bar4 = copy(dxidx_bar1)
  jac_bar2 = copy(jac_bar1); jac_bar3 = copy(jac_bar1)

  w_i = rand_realpart((mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim))
  op = SummationByParts.Subtract()

  h = 1e-20
  pert = Complex128(0, h)
  for i=1:mesh.numEl

    #w_i = sview(eqn.q, :, :, i)
    dxidx = sview(mesh.dxidx, :, :, :, i)
    jac = sview(mesh.jac, :, i)

    # perturb all the functions at once
    w_i .+= pert*w_dot
    dxidx .+= pert*dxidx_dot
    jac .+= pert*jac_dot
    EulerEquationMod.applyQxTransposed(sbp, w_i, dxidx, work, wx1, op)
    EulerEquationMod.applyDx(sbp, w_i, dxidx, jac, work, wx2, op)
    EulerEquationMod.applyDxTransposed(sbp, w_i, dxidx, jac, work, wx3, op)
    EulerEquationMod.applyQx(sbp, w_i, dxidx, work, wx4, op)
    w_i .-= pert*w_dot
    dxidx .-= pert*dxidx_dot
    jac .-= pert*jac_dot

 
    # test applyQxTransposed, second method
    EulerEquationMod.applyQxTransposed_revm(sbp, w_i, w_bar1, dxidx, dxidx_bar1, work, wx_bar, op)
    wx_dot = imag(wx1)./h
    val1 = sum(wx_dot .* wx_bar)
    val2 = sum(dxidx_bar1 .* dxidx_dot) + sum(w_bar1 .* w_dot)

    @test checkvals(val1, val2, 1e-12)

    # applyDx, second method
    EulerEquationMod.applyDx_revm(sbp, w_i, w_bar2, dxidx, dxidx_bar2, jac,
                                  jac_bar2, work, wx_bar, op)
    wx_dot = imag(wx2)./h
    val1 = sum(wx_dot .* wx_bar)
    val2 = sum(dxidx_bar2 .* dxidx_dot) + sum(jac_bar2 .* jac_dot) + sum(w_bar2 .* w_dot)

    @test checkvals(val1, val2, 1e-12)

    # applyDxTransposed, second method
    EulerEquationMod.applyDxTransposed_revm(sbp, w_i, w_bar3, dxidx, dxidx_bar3,
                                            jac, jac_bar3, work, wx_bar, op)
    wx_dot = imag(wx3)./h
    val1 = sum(wx_dot .* wx_bar)
    val2 = sum(dxidx_bar3 .* dxidx_dot) + sum(jac_bar3 .* jac_dot) + sum(w_bar3 .* w_dot)

    @test checkvals(val1, val2, 1e-12)

     # applyQx, second method
    EulerEquationMod.applyQx_revm(sbp, w_i, w_bar4, dxidx, dxidx_bar4, work, wx_bar, op)
    wx_dot = imag(wx4)./h
    val1 = sum(wx_dot .* wx_bar)
    val2 = sum(dxidx_bar4 .* dxidx_dot) + sum(w_bar4 .* w_dot)

    @test checkvals(val1, val2, 1e-12)


  end

  return nothing
end



"""
  Tests jacobian of `AbstractFaceShockCapturing` against complex step.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * capture: `AbstractFaceShockCapturing`

  **Keyword Arguments**

   * use_sensor: if true use `ShockSensorPP`, otherwise use `ShockSensorEverywher`
"""
function test_shock_capturing_jac(mesh, sbp, eqn::EulerData{Tsol, Tres}, _opts,
      capture::EulerEquationMod.AbstractShockCapturing; freeze=false) where {Tsol, Tres}

  opts = copy(_opts)  # don't modify the original
  opts["addVolumeIntegrals"] = false
  opts["addFaceIntegrals"] = false
  opts["addBoundaryIntegrals"] = false
  opts["addStabilization"] = false

  setShockSensorAlpha(eqn, 2)
  h = 1e-20
  pert = Complex128(0, h)

  # use a spatially varying solution
  srand(1234)
  icfunc = EulerEquationMod.ICDict["ICRho1E2U3"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
#  if use_sensor
    # this only works for the 2 element mesh
#    eqn.q[1] += 0.02
#    eqn.q[25] += 0.01
#  else
#    sensor = EulerEquationMod.ShockSensorEverywhere{Tsol, Tres}(mesh, sbp, opts)
#  end
  
  eqn.q .+= 0.01*rand(size(eqn.q))
  eqn.shock_capturing = capture

  # to test freezing the viscosity, have two different qs and evaluate the
  # Jacobian using the frozen viscosity while a different q in eqn.q
  if freeze
    q_frozen = copy(eqn.q)
    EulerEquationMod.freezeViscosity(mesh, sbp, eqn, opts)
    eqn.q .+= 0.01*rand(size(eqn.q))
  end


  # set a boundary condition where getDirichletState depends on qface
  opts["BC1_name"] = "noPenetrationBC"
  EulerEquationMod.getBCFunctors(mesh, sbp, eqn, opts)

  # get explicitly computed linear operator
  opts["calc_jac_explicit"] = true
  opts["addShockCapturing"] = true
#  val_orig = opts["preallocate_jacobian_coloring"]
#  opts["preallocate_jacobian_coloring"] = true
  pc, lo = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
#  opts["preallocate_jacobian_coloring"] = val_orig

  jac = getBaseLO(lo).A
  assembler = Jacobian._AssembleElementData(jac, mesh, sbp, eqn, opts)

#  println("mesh.dofs[:, :, 1] = \n", mesh.dofs[:, :, 1])
#  println("mesh.dofs[:, :, 5] = \n", mesh.dofs[:, :, 5])

  # complex step

  q_dot = rand_realpart(size(eqn.q_vec))
#  q_dot = zeros(length(eqn.q_vec)); q_dot[1] = 1
  res_dot = zeros(eqn.res_vec)
  evaldRdqProduct(mesh, sbp, eqn, opts, q_dot, res_dot)

  # explictly compute jacobian
#  EulerEquationMod.applyShockCapturing_diff(mesh, sbp, eqn, opts, sensor, capture, assembler)
  evalJacobian(mesh, sbp, eqn, opts, assembler)

  assembly_begin(jac, MAT_FINAL_ASSEMBLY)
  assembly_end(jac, MAT_FINAL_ASSEMBLY)


  b = zeros(mesh.numDof)
  t = 0.0
  ctx_residual = (evalResidual,)
  applyLinearOperator(lo, mesh, sbp, eqn, opts, ctx_residual, t, q_dot, b)

  if freeze
    EulerEquationMod.thawViscosity(mesh, sbp, eqn, opts)
  end


#  dofs1 = vec(mesh.dofs[:, :, 1])
#  println("res_dot = ", real(res_dot))
#  println("b = ", real(b))
#  dofs5 = vec(mesh.dofs[:, :, 5])
#  b = jac*q_dot
#  println("diff = ", real(b - res_dot))
#  println("\nb = ", real(b))
#  println("\nres_dot = ", real(res_dot))
#  println("diff1 = \n", real(b[dofs1] - res_dot[dofs1]))
#  println("diff5 = \n", real(b[dofs5] - res_dot[dofs5]))
#  println("max(b) = ", maximum(abs.(b)))
#  println("max(res_dot) = ", maximum(abs.(res_dot)))
#  println("maxdiff = ",  maximum(abs.(b - res_dot)))
  @test maximum(abs.(b - res_dot)) < 1e-11

  #println("b = \n", real(b))
  #println("res_dot = \n", real(res_dot))
  #println("q_dot = ", real(q_dot))
  #println("full(A) = ", full(jac))
  

  if typeof(jac) <: PetscMat
    matinfo = MatGetInfo(jac, PETSc2.MAT_LOCAL)
#    if is_prealloc_exact
      @test ( matinfo.nz_unneeded )== 0
#    else
#      @test  matinfo.nz_unneeded  > 0
#    end
  end

  setShockSensorAlpha(eqn, 1)

  return nothing 
end
