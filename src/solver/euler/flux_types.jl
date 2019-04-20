# types needed by flux functions
# these need to be declared in a separate file from the flux functions
# themselves because they need to be defind *before* ParamType

"""
  Data for [`RoeSolver`](@ref) and related methods
"""
struct RoeFluxData{Tsol, Tres, Tmsh}
  # primal flux
  dq::Vector{Tsol}
  sat::Vector{Tres}
  roe_vars::Vector{Tres}
  euler_flux::Vector{Tres}
  v_vals::Vector{Tsol}
  nrm2::Vector{Tmsh}

  # revm
  euler_flux_bar::Vector{Tres}
  sat_bar::Vector{Tres}

  # revq
  v_vals_bar::Vector{Tres}
  dq_bar::Vector{Tres}
  roe_vars_bar::Vector{Tres}

  # diff
  roe_vars_dot::Vector{Tres}

  function RoeFluxData{Tsol, Tres, Tmsh}(numDofPerNode::Integer, dim::Integer) where {Tsol, Tres, Tmsh}

    # primal flux
    dq = zeros(Tsol, numDofPerNode)
    sat = zeros(Tres, numDofPerNode)
    roe_vars = zeros(Tres, numDofPerNode)
    euler_flux = zeros(Tres, numDofPerNode)
    v_vals = zeros(Tsol, numDofPerNode)
    nrm2 = zeros(Tmsh, dim)

    # revm
    euler_flux_bar = zeros(Tres, numDofPerNode)
    sat_bar = zeros(Tres, numDofPerNode)

    # revq
    v_vals_bar = zeros(Tres, numDofPerNode)
    dq_bar = zeros(Tres, numDofPerNode)
    roe_vars_bar = zeros(Tres, 4)

    # diff
    roe_vars_dot = zeros(Tres, 22)  # 22 needed in 3D

    obj = new(dq, sat, roe_vars, euler_flux, v_vals, nrm2,
              euler_flux_bar, sat_bar,
              v_vals_bar, dq_bar, roe_vars_bar,
              roe_vars_dot)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end

end

"""
  Data for [`calcSat`](@ref) and related methods 
"""
struct CalcSatData{Tres}
  E1dq::Vector{Tres}
  E2dq::Vector{Tres}
  E3dq::Vector{Tres}
  E4dq::Vector{Tres}

  # no arrays are used for the diff version

  # revq
  E1dq_bar::Vector{Tres}
  E2dq_bar::Vector{Tres}
  E3dq_bar::Vector{Tres}
  E4dq_bar::Vector{Tres}

  function CalcSatData{Tres}(numDofPerNode::Integer) where {Tres}

    E1dq = zeros(Tres, numDofPerNode)
    E2dq = zeros(Tres, numDofPerNode)
    E3dq = zeros(Tres, numDofPerNode)
    E4dq = zeros(Tres, numDofPerNode)

    # revq
    E1dq_bar = zeros(Tres, numDofPerNode)
    E2dq_bar = zeros(Tres, numDofPerNode)
    E3dq_bar = zeros(Tres, numDofPerNode)
    E4dq_bar = zeros(Tres, numDofPerNode)



    obj = new(E1dq, E2dq, E3dq, E4dq,
              E1dq_bar, E2dq_bar, E3dq_bar, E4dq_bar)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end



"""
  Data needed by [`calcLFFlux`](@ref) and related method

  **Static Parameters**

   * Tres
"""
struct LFFluxData{Tres}
  fluxL::Vector{Tres}
  fluxR::Vector{Tres}
  F_dotL::Array{Tres, 2}
  F_dotR::Array{Tres, 2}
  lambda_dotL::Vector{Tres}
  lambda_dotR::Vector{Tres}

  function LFFluxData{Tres}(numDofPerNode::Integer, nd::Integer) where {Tres}

    @assert numDofPerNode == nd  # not sure what happens if this is not true
    fluxL = zeros(Tres, numDofPerNode)
    fluxR = zeros(Tres, numDofPerNode)
    F_dotL = zeros(Tres, numDofPerNode, nd)
    F_dotR = zeros(Tres, numDofPerNode, nd)
    lambda_dotL = zeros(Tres, numDofPerNode)
    lambda_dotR = zeros(Tres, numDofPerNode)

    obj = new(fluxL, fluxR, F_dotL, F_dotR, lambda_dotL, lambda_dotR)

    assertArraysUnique(obj); assertFieldsConcrete(obj)
    
    return obj
  end
end


"""
  Data needed by [`logavg_diff`](@ref)

  **Static Parameters**

   * Tl: datatype of left state
   * Tr: datatype of right state
"""
struct LogAvgData{Tl, Tr}
  xi_dotL::Vector{Tl}
  xi_dotR::Vector{Tr}
  f_dotL::Vector{Tl}
  f_dotR::Vector{Tr}
  u_dotL::Vector{Tl}
  u_dotR::Vector{Tr}
  F_dotL::Vector{Tl}
  F_dotR::Vector{Tr}

  function LogAvgData{Tl, Tr}(nd::Integer) where {Tl, Tr}
    xi_dotL = zeros(Tl, nd)
    xi_dotR = zeros(Tr, nd)
    f_dotL = zeros(Tl, nd)
    f_dotR = zeros(Tr, nd)
    u_dotL = zeros(Tl, nd)
    u_dotR = zeros(Tr, nd)
    F_dotL = zeros(Tl, nd)
    F_dotR = zeros(Tr, nd)

    return new(xi_dotL, xi_dotR, f_dotL, f_dotR, u_dotL, u_dotR, F_dotL, F_dotR)
  end
end

"""
  This struct holds all the temporary arrays needed to calculate the IR flux
"""
struct IRFluxData{Tsol}
  pL_dot::Vector{Tsol}
  pR_dot::Vector{Tsol}
  z1L_dot::Vector{Tsol}
  z2L_dot::Vector{Tsol}
  z3L_dot::Vector{Tsol}
  z4L_dot::Vector{Tsol}
  z5L_dot::Vector{Tsol}

  z1R_dot::Vector{Tsol}
  z2R_dot::Vector{Tsol}
  z3R_dot::Vector{Tsol}
  z4R_dot::Vector{Tsol}
  z5R_dot::Vector{Tsol}

  z5avg_dotL::Vector{Tsol}
  z5avg_dotR::Vector{Tsol}
  z4avg_dotL::Vector{Tsol}
  z4avg_dotR::Vector{Tsol}
  z1avg_dotL::Vector{Tsol}
  z1avg_dotR::Vector{Tsol}

  rho_hat_dotL::Vector{Tsol}
  rho_hat_dotR::Vector{Tsol}
  u_hat_dotL::Vector{Tsol}
  u_hat_dotR::Vector{Tsol}
  v_hat_dotL::Vector{Tsol}
  v_hat_dotR::Vector{Tsol}
  w_hat_dotL::Vector{Tsol}
  w_hat_dotR::Vector{Tsol}
  p1_hat_dotL::Vector{Tsol}
  p1_hat_dotR::Vector{Tsol}
  h_hat_dotL::Vector{Tsol}
  h_hat_dotR::Vector{Tsol}
  avgdata::LogAvgData{Tsol, Tsol}

  """
    Constructor for IRFluxData.

    **Inputs**

     * nd: number of derivative components
  """
  function IRFluxData{Tsol}(nd::Integer) where {Tsol}

    #TODO: consider making these views of an array to get spatial locality
    pL_dot = zeros(Tsol, nd)
    pR_dot = zeros(Tsol, nd)

    z1L_dot = zeros(Tsol, nd)
    z2L_dot = zeros(Tsol, nd)
    z3L_dot = zeros(Tsol, nd)
    z4L_dot = zeros(Tsol, nd)
    z5L_dot = zeros(Tsol, nd)

    z1R_dot = zeros(Tsol, nd)
    z2R_dot = zeros(Tsol, nd)
    z3R_dot = zeros(Tsol, nd)
    z4R_dot = zeros(Tsol, nd)
    z5R_dot = zeros(Tsol, nd)

    z5avg_dotL = zeros(Tsol, nd)
    z5avg_dotR = zeros(Tsol, nd)
    z4avg_dotL = zeros(Tsol, nd)
    z4avg_dotR = zeros(Tsol, nd)
    z1avg_dotL = zeros(Tsol, nd)
    z1avg_dotR = zeros(Tsol, nd)

    rho_hat_dotL = zeros(Tsol, nd)
    rho_hat_dotR = zeros(Tsol, nd)
    u_hat_dotL = zeros(Tsol, nd)
    u_hat_dotR = zeros(Tsol, nd)
    v_hat_dotL = zeros(Tsol, nd)
    v_hat_dotR = zeros(Tsol, nd)
    w_hat_dotL = zeros(Tsol, nd)
    w_hat_dotR = zeros(Tsol, nd)
    p1_hat_dotL = zeros(Tsol, nd)
    p1_hat_dotR = zeros(Tsol, nd)
    h_hat_dotL = zeros(Tsol, nd)
    h_hat_dotR = zeros(Tsol, nd)

    logdata = LogAvgData{Tsol, Tsol}(nd)

    obj =  new(pL_dot, pR_dot, z1L_dot, z2L_dot, z3L_dot, z4L_dot, z5L_dot,
               z1R_dot, z2R_dot, z3R_dot, z4R_dot, z5R_dot,
               z5avg_dotL, z5avg_dotR, z4avg_dotL, z4avg_dotR, z1avg_dotL,
               z1avg_dotR,
               rho_hat_dotL, rho_hat_dotR, u_hat_dotL, u_hat_dotR,
               v_hat_dotL, v_hat_dotR, w_hat_dotL, w_hat_dotR,
               p1_hat_dotL, p1_hat_dotR, h_hat_dotL,h_hat_dotR,
               logdata)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


"""
  Temporary arrays for [`applyEntropyKernel_diagE`](@ref)
"""
struct ApplyEntropyKernel_diagEData{Tsol, Tres}
  q_avg::Vector{Tsol}
  
  # diff method
  q_avg_dot::Matrix{Tres}
  F::Vector{Tres}
  F_dot::Matrix{Tres}

  # for applyEntropyKernel_inner
  vL::Vector{Tsol}
  vR::Vector{Tsol}
  F_tmp::Vector{Tres}

  # diff method
  delta_w_dot::Matrix{Tsol}

  # rev method
  q_avg_bar::Vector{Tres}
  delta_w::Vector{Tsol}
  delta_w_bar::Vector{Tsol}
  A0inv::Matrix{Tsol}

  function ApplyEntropyKernel_diagEData{Tsol, Tres}(numDofPerNode::Integer, nd::Integer) where {Tsol, Tres}

    q_avg = zeros(Tsol, numDofPerNode)
    vL = zeros(Tsol, numDofPerNode)
    vR = zeros(Tsol, numDofPerNode)
    F_tmp = zeros(Tres, numDofPerNode)

    # diff  method
    q_avg_dot = zeros(Tsol, numDofPerNode, nd)
    F = zeros(Tres, numDofPerNode)
    F_dot = zeros(Tres, numDofPerNode, nd)

    delta_w_dot = zeros(Tsol, numDofPerNode, nd)

    # rev
    q_avg_bar = zeros(Tres, numDofPerNode)
    delta_w = zeros(Tsol, numDofPerNode)
    delta_w_bar = zeros(Tsol, numDofPerNode)
    A0inv = zeros(Tsol, numDofPerNode, numDofPerNode)



    obj = new(q_avg, q_avg_dot, F, F_dot, vL, vR, F_tmp, delta_w_dot,
              q_avg_bar, delta_w, delta_w_bar, A0inv)

    assertArraysUnique(obj); assertFieldsConcrete(obj)
    return obj
  end
end

"""
  Temporary storage for [`getLambdaMaxSimple`](@ref) and differentiated versions
"""
struct GetLambdaMaxSimpleData{Tsol}
  q_avg::Vector{Tsol}

  function GetLambdaMaxSimpleData{Tsol}(numDofPerNode::Integer) where {Tsol}

    q_avg = zeros(Tsol, numDofPerNode)

    obj = new(q_avg)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end

"""
  Temporary storage for [`getLambdaMax`](@ref) and differentiated versions.
"""
struct GetLambdaMaxData{Tsol}
  p_dot::Vector{Tsol}

  function GetLambdaMaxData{Tsol}(numDofPerNode::Integer) where {Tsol}

    p_dot = zeros(Tsol, numDofPerNode)

    return new(p_dot)
  end
end


struct CalcVorticityData{Tsol, Tres, Tmsh}
  dxidx_unscaled::Array{Tmsh, 3}
  velocities::Matrix{Tsol}
  velocity_deriv::Array{Tmsh, 3}
  velocity_deriv_xy::Array{Tmsh, 3}

  function CalcVorticityData{Tsol, Tres, Tmsh}(numNodesPerElement::Integer, dim::Integer) where {Tsol, Tres, Tmsh}

    dxidx_element = zeros(Tmsh, dim, dim, numNodesPerElement)
    velocities = zeros(Tsol, dim, numNodesPerElement)
    velocity_deriv = zeros(Tsol, dim, numNodesPerElement, dim)
    velocity_deriv_xy = zeros(Tres, dim, dim, numNodesPerElement)

    obj = new(dxidx_element, velocities, velocity_deriv, velocity_deriv_xy)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


struct ContractResEntropyVarsData{Tsol}
  w_vals::Vector{Tsol}

  function ContractResEntropyVarsData{Tsol}(numDofPerNode::Integer) where {Tsol}

    w_vals = zeros(Tsol, numDofPerNode)

    return new(w_vals)
  end
end

"""
  Temporary storage for[`getTau`](@ref)
"""
struct GetTauData{Tsol, Tres}
  AjAk::Matrix{Tsol}
  flux_term::Matrix{Tsol}
  A0inv::Matrix{Tsol}
  tmp_mat::Matrix{Tres}

  # second method
  B_d::Matrix{Tres}
  B_p::Matrix{Tres}
  A_mat_hat::Array{Tsol, 3}
  A0::Matrix{Tsol}
  tmp_mat2::Matrix{Tsol}

  function GetTauData{Tsol, Tres}(numDofPerNode::Integer, dim::Integer) where {Tsol, Tres}

    AjAk = zeros(Tsol, numDofPerNode, numDofPerNode)
    flux_term = zeros(Tsol, numDofPerNode, numDofPerNode)
    A0inv = zeros(Tsol, numDofPerNode, numDofPerNode)
    tmp_mat = zeros(Tres, numDofPerNode, numDofPerNode)

    B_d = zeros(Tres, numDofPerNode, numDofPerNode)
    B_p = zeros(Tres, numDofPerNode, numDofPerNode)
    A_mat_hat = zeros(Tsol, numDofPerNode, numDofPerNode, dim)
    A0 = zeros(Tsol, numDofPerNode, numDofPerNode)
    tmp_mat2 = zeros(Tsol, numDofPerNode, numDofPerNode)


    obj = new(AjAk, flux_term, A0inv, tmp_mat,
              B_d, B_p, A_mat_hat, A0, tmp_mat2)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end





#------------------------------------------------------------------------------
# structs for face element integrals

"""
  Temporary storage for [`calcECFaceIntegral`](@ref)
"""
struct CalcECFaceIntegralData{Tres, Tmsh}
  fluxD::Matrix{Tres}
  nrmD::Matrix{Tmsh}

  # for SparseFace method
  flux_tmp::Vector{Tres}

  # diff method
  fL_dotD::Array{Tres, 3}
  fR_dotD::Array{Tres, 3}
  fL_dot::Array{Tres, 2}
  fR_dot::Array{Tres, 2}

  # rev methods
  fluxD_bar::Matrix{Tres}

  function CalcECFaceIntegralData{Tres, Tmsh}(numDofPerNode::Integer, dim::Integer) where {Tres, Tmsh}
    fluxD = zeros(Tres, numDofPerNode, dim)
    nrmD = zeros(Tmsh, dim, dim)
    flux_tmp = zeros(Tres, numDofPerNode)

    fL_dotD = zeros(Tres, numDofPerNode, numDofPerNode, dim)
    fR_dotD = zeros(Tres, numDofPerNode, numDofPerNode, dim)

    fL_dot = zeros(Tres, numDofPerNode, numDofPerNode)
    fR_dot = zeros(Tres, numDofPerNode, numDofPerNode)

    fluxD_bar = zeros(Tres, numDofPerNode, dim)

    obj = new(fluxD, nrmD, flux_tmp, fL_dotD, fR_dotD, fL_dot, fR_dot,
              fluxD_bar)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


struct CalcEntropyPenaltyIntegralData{Tsol, Tres}
  wL::Matrix{Tsol}
  wR::Matrix{Tsol}
  wL_i::Vector{Tsol}
  wR_i::Vector{Tsol}
  qL_i::Vector{Tsol}
  qR_i::Vector{Tsol}
  flux::Vector{Tres}
  A0::Matrix{Tsol}

  # sparseface method
  delta_w::Vector{Tsol}
  q_avg::Vector{Tsol}
  res_vals::Vector{Tres}

  #------------------
  # diff method
  q_avg_dot::Matrix{Tsol}
  delta_w_dot::Matrix{Tsol}
  flux_dot_i::Matrix{Tres}  # used to compute the flux

  # arrays in the format required by SBP
  flux_dotL::Array{Tres, 3}
  flux_dotR::Array{Tres, 3}

  jacLL_tmp::Array{Tres, 4}
  jacLR_tmp::Array{Tres, 4}
  jacRL_tmp::Array{Tres, 4}
  jacRR_tmp::Array{Tres, 4}

  A0invL::Matrix{Tsol}
  A0invR::Matrix{Tsol}

  #---------------
  # rev methods
  flux_bar::Vector{Tres}
  qL_bar_i::Vector{Tres}
  qR_bar_i::Vector{Tres}
  wL_bar_i::Vector{Tres}
  wR_bar_i::Vector{Tres}
  delta_w_bar::Vector{Tres}
  q_avg_bar::Vector{Tres}
  wL_bar::Matrix{Tres} 
  wR_bar::Matrix{Tres}




  function CalcEntropyPenaltyIntegralData{Tsol, Tres}(numDofPerNode::Integer,
                              numNodesPerFace::Integer,
                              stencilsize::Integer, numNodesPerElement::Integer,
                              nd::Integer) where {Tsol, Tres}

    wL = zeros(Tsol, numDofPerNode, stencilsize)
    wR = zeros(Tsol, numDofPerNode, stencilsize)
    wL_i = zeros(Tsol, numDofPerNode)
    wR_i = zeros(Tsol, numDofPerNode)
    qL_i = zeros(Tsol, numDofPerNode)
    qR_i = zeros(Tsol, numDofPerNode)
    flux = zeros(Tres, numDofPerNode)
    A0 = zeros(Tsol, numDofPerNode, numDofPerNode)

    delta_w = zeros(Tsol, numDofPerNode)
    q_avg = zeros(Tsol, numDofPerNode)
    res_vals = zeros(Tres, numDofPerNode)


    # diff method
    q_avg_dot = zeros(Tsol, numDofPerNode, nd)
    delta_w_dot = zeros(Tsol, numDofPerNode, nd)
    flux_dot_i = zeros(Tres, numDofPerNode, nd)

    # arrays needed by SBP function
    flux_dotL = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerFace)
    flux_dotR = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerFace)

    jacLL_tmp = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement, numNodesPerElement)
    jacLR_tmp = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement, numNodesPerElement)
    jacRL_tmp = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement, numNodesPerElement)
    jacRR_tmp = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement, numNodesPerElement)

    A0invL = zeros(Tsol, numDofPerNode, numDofPerNode)
    A0invR = zeros(Tsol, numDofPerNode, numDofPerNode)

    # rev methods
    flux_bar = zeros(Tres, numDofPerNode)

    qL_bar_i = zeros(Tres, numDofPerNode)
    qR_bar_i = zeros(Tres, numDofPerNode)
    wL_bar_i = zeros(Tres, numDofPerNode)
    wR_bar_i = zeros(Tres, numDofPerNode)
    delta_w_bar = zeros(Tres, numDofPerNode)
    q_avg_bar = zeros(Tres, numDofPerNode)
    wL_bar = zeros(Tres, numDofPerNode, stencilsize)
    wR_bar = zeros(Tres, numDofPerNode, stencilsize)



    obj = new(wL, wR, wL_i, wR_i, qL_i, qR_i, flux, A0, delta_w, q_avg,
              res_vals,
             q_avg_dot, delta_w_dot, flux_dot_i, flux_dotL, flux_dotR,
             jacLL_tmp, jacLR_tmp, jacRL_tmp, jacRR_tmp, A0invL, A0invR,
             flux_bar, qL_bar_i, qR_bar_i, wL_bar_i, wR_bar_i, delta_w_bar,
             q_avg_bar, wL_bar, wR_bar)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


"""
  Temporary storage for [`interpolateElementStaggered`](@ref)
"""
struct InterpolateElementStaggeredData{Tsol}
  wvars_s::Matrix{Tsol}
  wvars_f::Matrix{Tsol}

  function InterpolateElementStaggeredData{Tsol}(numDofPerNode::Integer,
              numNodesPerElement_s::Integer, numNodesPerElement_f) where {Tsol}

    wvars_s = zeros(Tsol, numDofPerNode, numNodesPerElement_s)
    wvars_f = zeros(Tsol, numDofPerNode, numNodesPerElement_f)

    obj = new(wvars_s, wvars_f)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end

"""
  Temporary storage for [`calcEulerFlux`](@ref) and differentiated versions
"""
struct CalcEulerFluxData{Tsol}
  p_dot::Vector{Tsol}

  function CalcEulerFluxData{Tsol}(numDofPerNode::Integer) where {Tsol}

    p_dot = zeros(Tsol, numDofPerNode)

    obj = new(p_dot)

    return obj
  end
end


"""
  Temporary storage for all BC functors.  This one is a bit odd, because
  we can't store these arrays in the BC functors themselves, because when
  doing boundary conditions for explicitly computed jacobians, the
  Tsol of the eqn object and the type of the solution the boundary condition
  functor operates on are different.

  This design needs to be rethought.
"""
struct BCData{Tsol, Tres}

  qg::Vector{Tsol}
  v_vals::Vector{Tsol}

  sat::Vector{Tres}
  dq::Vector{Tsol}
  roe_vars::Vector{Tsol}
  euler_flux::Vector{Tres}

  # reverse mode
  q_bar::Vector{Tres}
  qg_bar::Vector{Tres}

  function BCData{Tsol, Tres}(numDofPerNode::Integer) where {Tsol, Tres}

    qg = zeros(Tsol, numDofPerNode)
    v_vals = zeros(Tsol, numDofPerNode)
    sat = zeros(Tsol, numDofPerNode)
    dq = zeros(Tsol, numDofPerNode)
    roe_vars = zeros(Tsol, numDofPerNode)
    euler_flux = zeros(Tres, numDofPerNode)

    # reverse mode
    q_bar = zeros(Tres, numDofPerNode)
    qg_bar = zeros(Tres, numDofPerNode)

    obj = new(qg, v_vals, sat, dq, roe_vars, euler_flux,
              q_bar, qg_bar)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


#------------------------------------------------------------------------------
# entropy kernel functions for entropy-stable scheme

struct LW2Kernel{Tsol, Tres, Tmsh} <: AbstractEntropyKernel
  nrm::Array{Tmsh, 1}
  P::Array{Tmsh, 2}
  Y::Array{Tsol, 2}  # eigenvectors
  Lambda::Array{Tsol, 1}  # eigenvalues
  S2::Array{Tsol, 1}  # scaling for the eigensystem
  q_tmp::Array{Tsol, 1}
  tmp1::Array{Tres, 1}
  tmp2::Array{Tres, 1}

  function LW2Kernel{Tsol, Tres, Tmsh}(numDofPerNode::Integer, dim::Integer) where {Tsol, Tres, Tmsh}
    nrm = zeros(Tmsh, dim)
    P = zeros(Tmsh, numDofPerNode, numDofPerNode)
    Y = zeros(Tsol, numDofPerNode, numDofPerNode)
    Lambda = zeros(Tsol, numDofPerNode)
    S2 = zeros(Tsol, numDofPerNode)
    q_tmp = zeros(Tsol, numDofPerNode)
    tmp1 = zeros(Tres, numDofPerNode)
    tmp2 = zeros(Tres, numDofPerNode)

    obj =  new(nrm, P, Y, Lambda, S2, q_tmp, tmp1, tmp2)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end

"""
  Data for Lax-Friedrich entropy dissipation
"""
mutable struct LFKernel{Tsol, Tres, Tmsh} <: AbstractEntropyKernel
  A0::Array{Tsol, 2}

  # diff method
  A0_dot::Array{Tsol, 3}
  t1::Array{Tsol, 1}
  t1_dot::Array{Tsol, 2}
  lambda_max_dot::Array{Tres, 1}

  # revm method
  t1_bar::Array{Tres, 1}
  A0_bar::Array{Tres, 2}

  function LFKernel{Tsol, Tres, Tmsh}(numDofPerNode::Integer, nd::Integer) where {Tsol, Tres, Tmsh}

    A0 = zeros(Tsol, numDofPerNode, numDofPerNode)

    # diff method
    A0_dot = zeros(Tsol, numDofPerNode, numDofPerNode, nd)
    t1 = zeros(Tsol, numDofPerNode)
    t1_dot = zeros(Tsol, numDofPerNode, nd)
    lambda_max_dot = zeros(Tres, numDofPerNode)

    # revq method
    t1_bar = zeros(Tres, numDofPerNode)
    A0_bar = zeros(Tres, numDofPerNode, numDofPerNode)
   

    obj = new(A0, A0_dot, t1, t1_dot, lambda_max_dot, t1_bar, A0_bar)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


mutable struct IdentityKernel{Tsol, Tres, Tmsh} <: AbstractEntropyKernel

  function IdentityKernel{Tsol, Tres, Tmsh}() where {Tsol, Tres, Tmsh}

    return new()
  end
end

struct GetIRA0Data{Tsol}
  p_dot::Array{Tsol, 1}

  function GetIRA0Data{Tsol}(numDofPerNode::Integer) where {Tsol}

    p_dot = zeros(Tsol, numDofPerNode)

    return new(p_dot)
  end
end


#------------------------------------------------------------------------------
# structs for functions that loop over faces

"""
  Temporary storage for [`getFaceElementIntegral`](@ref) and related functions
"""
struct FaceElementIntegralData{Tsol, Tres}
  # staggered grid things
  resL_s::Matrix{Tres}
  resR_s::Matrix{Tres}
  resL_f::Matrix{Tres}
  resR_f::Matrix{Tres}
  qvars_f::Matrix{Tsol}

  function FaceElementIntegralData{Tsol, Tres}(numDofPerNode::Integer,
                                         numNodesPerElement_s::Integer,
                                         numNodesPerElement_f::Integer) where {Tsol, Tres}

    resL_s = zeros(Tres, numDofPerNode, numNodesPerElement_s)
    resR_s = zeros(Tres, numDofPerNode, numNodesPerElement_s)

    resL_f = zeros(Tres, numDofPerNode, numNodesPerElement_f)
    resR_f = zeros(Tres, numDofPerNode, numNodesPerElement_f)
    qvars_f = zeros(Tsol, numDofPerNode, numNodesPerElement_f)

    obj = new(resL_s, resR_s, resL_f, resR_f, qvars_f)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


"""
  Temporary storage for [`calcFaceIntegrals`](@ref) and related functions
  (including shared face integral functions
"""
struct CalcFaceIntegralsData{Tsol, Tres}
  q_faceL::Matrix{Tsol}
  q_faceR::Matrix{Tsol}

  # diff functions
  flux_dotL::Array{Tres, 3}
  flux_dotR::Array{Tsol, 3}
  res_jacLL::Array{Tres, 4}
  res_jacLR::Array{Tres, 4}
  res_jacRL::Array{Tres, 4}
  res_jacRR::Array{Tres, 4}

  function CalcFaceIntegralsData{Tsol, Tres}(numDofPerNode::Integer,
                    numNodesPerFace::Integer, numNodesPerElement::Integer) where {Tsol, Tres}

    q_faceL = zeros(Tsol, numDofPerNode, numNodesPerFace)
    q_faceR = zeros(Tsol, numDofPerNode, numNodesPerFace)

    flux_dotL = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerFace)
    flux_dotR = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerFace)

    res_jacLL = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement,
                            numNodesPerElement)
    res_jacLR = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement,
                            numNodesPerElement)
    res_jacRL = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement,
                            numNodesPerElement)
    res_jacRR = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement,
                            numNodesPerElement)


    obj = new(q_faceL, q_faceR, flux_dotL, flux_dotR, res_jacLL, res_jacLR, res_jacRL, res_jacRR)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


"""
  Temporary storage for [`calcVolumeIntegrals`](@ref) and
  [`calcVolumeIntegralsSplitForm`](@ref).
"""
struct CalcVolumeIntegralsData{Tres, Tmsh}

  # calcVolumeIntegrals
  flux_jac::Array{Tres, 4}
  res_jac::Array{Tres, 4}
  nrm::Vector{Tmsh}

  # calcVolumeIntegralsSplitFormLinear
  nrmD::Matrix{Tmsh}
  F_d::Matrix{Tres}
  nrmD_bar::Matrix{Tmsh}
  F_d_bar::Matrix{Tres}

  S::Array{Float64, 3}  # SBP S matrix (should elementtype really be Float64?)

  # calcVolumeIntegralsSplitFormCurvilinear
  Sx::Array{Tmsh, 3}
  Sx_bar::Array{Tmsh, 3}

  # staggered grid methods
  res_s::Matrix{Tres}
  res_f::Matrix{Tres}

  function CalcVolumeIntegralsData{Tres, Tmsh}(numDofPerNode::Integer,
                  dim::Integer, numNodesPerElement_s::Integer,
                  numNodesPerElement_f::Integer, sbp::AbstractOperator) where {Tres, Tmsh}

    # calcVolumeIntegrals
    flux_jac = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement_s, dim)
    res_jac = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement_s, numNodesPerElement_s)
    nrm = zeros(Tmsh, dim)

    # split form stuff
    nrmD = zeros(Tmsh, dim, dim)
    F_d = zeros(Tres, numDofPerNode, dim)
    nrmD_bar = zeros(Tmsh, dim, dim)
    F_d_bar = zeros(Tres, numDofPerNode, dim)


    # the staggered grid calculation only uses the curvilinear method, so
    # make S numNodesPerElement_s, to match the sbp operator on the solution
    # grid (S is not used if they don't match>
    S = zeros(Float64, numNodesPerElement_s, numNodesPerElement_s, dim)
    Sx = zeros(Tmsh, numNodesPerElement_f, numNodesPerElement_f, dim)
    Sx_bar = zeros(Tmsh, numNodesPerElement_f, numNodesPerElement_f, dim)

    for i=1:dim
      S[:, :, i] = 0.5*(sbp.Q[:, :, i] - sbp.Q[:, :, i].')
    end

    res_s = zeros(Tres, numDofPerNode, numNodesPerElement_s)
    res_f = zeros(Tres, numDofPerNode, numNodesPerElement_f)

    obj = new(flux_jac, res_jac, nrm,
              nrmD, F_d, nrmD_bar, F_d_bar,  S, Sx, Sx_bar, res_s, res_f)

    assertArraysUnique(obj); assertFieldsConcrete(obj)

    return obj
  end
end


#------------------------------------------------------------------------------
# Stabilization structs

"""
  Data for Local Projection Stabiization (LPS).

  **Fields**

   * P: the projection matrix that captures the high frequency modes
        (I - L L^T H) in the paper
   * alpha: coefficient in front of the stabilization term
   * entropy_vars: specifies what variables to apply the stabiization to.
                   This is an abstractly-typed field, so it should be access
                   through a function barrier.
  
"""
struct LPSData{Tsol, Tres}
  P::Matrix{Float64}
  alpha::Float64
  entropy_vars::AbstractVariables

  t1::Matrix{Tsol}
  t2::Matrix{Tres}
  A0::Matrix{Tsol}
  dir::Vector{Tres}

  t1_jac::Array{Tsol, 3}
  t2_jac::Array{Tsol, 3}
  t3_jac::Array{Tres, 3}
  A0_diff::Array{Tsol, 3}
  q_dot::Matrix{Float64}
  lambdas::Vector{Tres}
  lambdas_jac::Matrix{Tres}
  lambda_dot_p::Vector{Tres}

  t1_bar::Matrix{Tres}
  t2_bar::Matrix{Tres}
  tmp::Vector{Tres}
  tmp_bar::Vector{Tres}
  A0_bar::Matrix{Tres}

  dir_bar::Vector{Tres}
  function LPSData{Tsol, Tres}(mesh, sbp, opts) where {Tsol, Tres}

    @unpack mesh numDofPerNode numNodesPerElement dim
    P = getLPSMatrix(sbp)
    alpha = 1.0
    #entropy_vars = ConservativeVariables()
    entropy_vars = IRVariables()

    t1 = zeros(Tsol, numDofPerNode, numNodesPerElement)
    t2 = zeros(Tres, numDofPerNode, numNodesPerElement)
    A0 = zeros(Tsol, numDofPerNode, numDofPerNode)
    dir = zeros(Tres, dim)

    
    t1_jac = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement)
    t2_jac = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement)
    t3_jac = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement)
    A0_diff = zeros(Tsol, numDofPerNode, numDofPerNode, numDofPerNode)
    q_dot = zeros(numDofPerNode, numDofPerNode)
    for i=1:numDofPerNode
      q_dot[i, i] = 1
    end

    lambdas = zeros(Tres, numNodesPerElement)
    lambdas_jac = zeros(Tres, numDofPerNode, numNodesPerElement)
    lambda_dot_p = zeros(Tres, numDofPerNode)


    t1_bar = zeros(Tres, numDofPerNode, numNodesPerElement)
    t2_bar = zeros(Tres, numDofPerNode, numNodesPerElement)
    tmp = zeros(Tres, numDofPerNode)
    tmp_bar = zeros(Tres, numDofPerNode)
    A0_bar = zeros(Tres, numDofPerNode, numDofPerNode)

  
    dir_bar = zeros(Tres, dim)

    return new(P, alpha, entropy_vars,
               t1, t2, A0, dir,
               t1_jac, t2_jac, t3_jac, A0_diff, q_dot, lambdas, lambdas_jac,
               lambda_dot_p,
               t1_bar, t2_bar, tmp, tmp_bar, A0_bar,
               dir_bar)
  end
end


