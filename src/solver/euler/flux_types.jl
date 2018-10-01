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

    # diff
    roe_vars_dot = zeros(Tres, 22)  # 22 needed in 3D

    obj = new(dq, sat, roe_vars, euler_flux, v_vals, nrm2,
              euler_flux_bar, sat_bar,
              roe_vars_dot)

    assertArraysUnique(obj)

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

  function CalcSatData{Tres}(numDofPerNode::Integer) where {Tres}

    E1dq = zeros(Tres, numDofPerNode)
    E2dq = zeros(Tres, numDofPerNode)
    E3dq = zeros(Tres, numDofPerNode)
    E4dq = zeros(Tres, numDofPerNode)


    obj = new(E1dq, E2dq, E3dq, E4dq)

    assertArraysUnique(obj)

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

  function LFFluxData{Tres}(numDofPerNode::Integer, nd::Integer) where {Tres}

    @assert numDofPerNode == nd  # not sure what happens if this is not true
    fluxL = zeros(Tres, numDofPerNode)
    fluxR = zeros(Tres, numDofPerNode)
    F_dotL = zeros(Tres, numDofPerNode, nd)
    F_dotR = zeros(Tres, numDofPerNode, nd)

    obj = new(fluxL, fluxR, F_dotL, F_dotR)

    assertArraysUnique(obj)
    
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

    assertArraysUnique(obj)

    return obj
  end
end


"""
  Temporary arrays for [`getEntropyLFStab`](@ref)
"""
struct GetEntropyLFStabData{Tsol}
  q_avg::Vector{Tsol}
  
  # for getEntropyLFStab_inner
  vL::Vector{Tsol}
  vR::Vector{Tsol}
  A0::Matrix{Tsol}

  function GetEntropyLFStabData{Tsol}(numDofPerNode::Integer) where {Tsol}

    q_avg = zeros(Tsol, numDofPerNode)
    vL = zeros(Tsol, numDofPerNode)
    vR = zeros(Tsol, numDofPerNode)
    A0 = zeros(Tsol, numDofPerNode, numDofPerNode)

    obj = new(q_avg, vL, vR, A0)

    assertArraysUnique(obj)
    return obj
  end
end


struct GetLambdaMaxSimpleData{Tsol}
  q_avg::Vector{Tsol}

  function GetLambdaMaxSimpleData{Tsol}(numDofPerNode::Integer) where {Tsol}

    q_avg = zeros(Tsol, numDofPerNode)

    obj = new(q_avg)

    assertArraysUnique(obj)

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

  function CalcECFaceIntegralData{Tres, Tmsh}(numDofPerNode::Integer, dim::Integer) where {Tres, Tmsh}
    fluxD = zeros(Tres, numDofPerNode, dim)
    nrmD = zeros(Tmsh, dim, dim)
    flux_tmp = zeros(Tres, numDofPerNode)

    obj = new(fluxD, nrmD, flux_tmp)

    assertArraysUnique(obj)

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

  function CalcEntropyPenaltyIntegralData{Tsol, Tres}(numDofPerNode::Integer,
                                    stencilsize::Integer) where {Tsol, Tres}

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

    obj = new(wL, wR, wL_i, wR_i, qL_i, qR_i, flux, A0, delta_w, q_avg, res_vals)

    assertArraysUnique(obj)

    return obj
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

    assertArraysUnique(obj)

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

    assertArraysUnique(obj)

    return obj
  end
end
