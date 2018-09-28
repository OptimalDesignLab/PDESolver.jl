# types needed by flux functions
# these need to be declared in a separate file from the flux functions
# themselves because they need to be defind *before* ParamType

"""
  Data needed by [`calcLFFlux_diff`](@ref)

  **Static Parameters**

   * Tres
"""
struct LFFluxData{Tres}
  F_dotL::Array{Tres, 2}
  F_dotR::Array{Tres, 2}

  function LFFluxData{Tres}(numDofPerNode::Integer, nd::Integer) where {Tres}

    @assert numDofPerNode == nd  # not sure what happens if this is not true
    F_dotL = zeros(Tres, numDofPerNode, nd)
    F_dotR = zeros(Tres, numDofPerNode, nd)

    obj = new(F_dotL, F_dotR)

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


