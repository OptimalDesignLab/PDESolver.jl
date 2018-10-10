# differentiated version of IR_stab.jl

function getIRA0_diff(params::ParamType{2}, 
                      q::AbstractArray{Tsol,1},
                      q_dot::AbstractArray{Tsol, 2},
                      A0::AbstractArray{Tsol, 2},
                      A0_dot::AbstractArray{Tsol, 3}) where Tsol

  nd = size(q_dot, 2)  # number of dot vectors
  numDofPerNode = size(q, 1)
  
  @debug1 begin
    @assert size(q, 1) == size(q_dot, 1)
    @assert size(A0_dot, 3) >= nd
    @assert size(A0, 1) == size(A0_dot, 1)
    @assert size(A0, 2) == size(A0_dot, 2)
  end


  p_dot = params.get_ira0data.p_dot
  gamma = params.gamma
  gamma_1 = params.gamma_1
  
  # p_dot is needed for i=1:nd, so compute the derivative with respect to
  # q[1:4] and then multiply with q_dot[:, i] later
  p = calcPressure_diff(params, q, p_dot)

  rho = q[1]
  rhou = q[2]
  rhov = q[3]
  rhoe = q[4]

  rhoinv = 1/rho

  h = (rhoe + p)*rhoinv  # this isn't really h, but including the factor of
                         # 1/rho is convenient
  a2 = gamma*p*rhoinv  # speed of sound

  A0[1,1] = rho
  A0[2,1] = rhou
  A0[3,1] = rhov
  A0[4,1] = rhoe

  A0[1,2] = rhou
  A0[2,2] = rhou*rhou*rhoinv + p
  A0[3,2] = rhou*rhov*rhoinv
  A0[4,2] = rhou*h

  A0[1,3] = rhov
  A0[2,3] = rhou*rhov*rhoinv
  A0[3,3] = rhov*rhov*rhoinv + p
  A0[4,3] = rhov*h

  A0[1,4] = rhoe
  A0[2,4] = h*rhou
  A0[3,4] = h*rhov
  A0[4,4] = rho*h*h - a2*p/gamma_1

  for i=1:nd

    # compute p_dot
    p_dot_i = zero(Tsol)
    for j=1:numDofPerNode
      p_dot_i += p_dot[j]*q_dot[j, i]
    end

    rho_dot  = q_dot[1, i]
    rhou_dot = q_dot[2, i]
    rhov_dot = q_dot[3, i]
    rhoe_dot = q_dot[4, i]

    rhoinv = 1/rho
    rhoinv_dot = -rhoinv*rhoinv*rho_dot

    h = (rhoe + p)*rhoinv
    h_dot = (rhoe + p)*rhoinv_dot + (rhoe_dot + p_dot_i)*rhoinv

    a2 = gamma*p*rhoinv
    a2_dot = gamma*(p_dot_i*rhoinv + p*rhoinv_dot)

    A0_dot[1,1,i] = rho_dot
    A0_dot[2,1,i] = rhou_dot
    A0_dot[3,1,i] = rhov_dot
    A0_dot[4,1,i] = rhoe_dot

    A0_dot[1,2,i] = rhou_dot
    A0_dot[2,2,i] = 2*rhou*rhoinv*rhou_dot + rhou*rhou*rhoinv_dot + p_dot_i
    A0_dot[3,2,i] = rhov*rhoinv*rhou_dot   + rhou*rhoinv*rhov_dot + rhou*rhov*rhoinv_dot
    A0_dot[4,2,i] = h*rhou_dot + rhou*h_dot

    A0_dot[1,3,i] = rhov_dot
    A0_dot[2,3,i] = rhov*rhoinv*rhou_dot   + rhou*rhoinv*rhov_dot + rhou*rhov*rhoinv_dot
    A0_dot[3,3,i] = 2*rhov*rhoinv*rhov_dot + rhov*rhov*rhoinv_dot + p_dot_i
    A0_dot[4,3,i] = h*rhov_dot + rhov*h_dot

    A0_dot[1,4,i] = rhoe_dot
    A0_dot[2,4,i] = rhou*h_dot  + h*rhou_dot
    A0_dot[3,4,i] = rhov*h_dot  + h*rhov_dot
    A0_dot[4,4,i] = h*h*rho_dot + rho*2*h*h_dot - (p*a2_dot + p_dot_i*a2)/gamma_1
  end

  return nothing
end

function getIRA0_diff(params::ParamType{3},
                      q::AbstractArray{Tsol,1},
                      q_dot::AbstractArray{Tsol, 2},
                      A0::AbstractArray{Tsol, 2},
                      A0_dot::AbstractArray{Tsol, 3}) where Tsol

  nd = size(q_dot, 2)  # number of dot vectors
  numDofPerNode = size(q, 1)

  @debug1 begin
    @assert size(q, 1) == size(q_dot, 1)
    @assert size(A0_dot, 3) == nd
    @assert size(A0, 1) == size(A0_dot, 1)
    @assert size(A0, 2) == size(A0_dot, 2)
  end

  p_dot = zeros(Tsol, numDofPerNode)
  gamma = params.gamma
  gamma_1 = params.gamma_1
  
  # p_dot is needed for i=1:nd, so compute the derivative with respect to
  # q[1:4] and then multiply with q_dot[:, i] later
  p = calcPressure_diff(params, q, p_dot)

  rho = q[1]
  rhou = q[2]
  rhov = q[3]
  rhow = q[4]
  rhoe = q[5]

  rhoinv = 1/rho

  h = (rhoe + p)*rhoinv
  a2 = gamma*p*rhoinv  # speed of sound

  A0[1,1] = rho
  A0[2,1] = rhou
  A0[3,1] = rhov
  A0[4,1] = rhow
  A0[5,1] = rhoe

  A0[1,2] = rhou
  A0[2,2] = rhou*rhou*rhoinv + p
  A0[3,2] = rhou*rhov*rhoinv
  A0[4,2] = rhou*rhow*rhoinv
  A0[5,2] = rhou*h

  A0[1,3] = rhov
  A0[2,3] = rhou*rhov/rho
  A0[3,3] = rhov*rhov*rhoinv + p
  A0[4,3] = rhov*rhow*rhoinv
  A0[5,3] = rhov*h


  A0[1,4] = rhow
  A0[2,4] = rhow*rhou*rhoinv
  A0[3,4] = rhow*rhov*rhoinv
  A0[4,4] = rhow*rhow*rhoinv + p
  A0[5,4] = rhow*h

  A0[1,5] = rhoe
  A0[2,5] = h*rhou
  A0[3,5] = h*rhov
  A0[4,5] = h*rhow
  A0[5,5] = rho*h*h - a2*p/gamma_1



  for i=1:nd

    # compute p_dot
    p_dot_i = zero(Tsol)
    for j=1:numDofPerNode
      p_dot_i += p_dot[j]*q_dot[j, i]
    end

    rho_dot  = q_dot[1, i]
    rhou_dot = q_dot[2, i]
    rhov_dot = q_dot[3, i]
    rhow_dot = q_dot[4, i]
    rhoe_dot = q_dot[5, i]

    rhoinv = 1/rho
    rhoinv_dot = -1*rhoinv*rhoinv*rho_dot

    h = (rhoe + p)*rhoinv
    h_dot = (rhoe + p)*rhoinv_dot + (rhoe_dot + p_dot_i)*rhoinv

    a2 = gamma*p*rhoinv
    a2_dot = gamma*(p_dot_i*rhoinv + p*rhoinv_dot)

    A0_dot[1,1,i] = rho_dot
    A0_dot[2,1,i] = rhou_dot
    A0_dot[3,1,i] = rhov_dot
    A0_dot[4,1,i] = rhow_dot
    A0_dot[5,1,i] = rhoe_dot

    A0_dot[1,2,i] = rhou_dot
    A0_dot[2,2,i] = 2*rhou*rhoinv*rhou_dot + rhou*rhou*rhoinv_dot + p_dot_i
    A0_dot[3,2,i] = rhov*rhoinv*rhou_dot   + rhou*rhoinv*rhov_dot + rhou*rhov*rhoinv_dot
    A0_dot[4,2,i] = rhow*rhoinv*rhou_dot   + rhou*rhoinv*rhow_dot + rhou*rhow*rhoinv_dot
    A0_dot[5,2,i] = h*rhou_dot + rhou*h_dot

    A0_dot[1,3,i] = rhov_dot
    A0_dot[2,3,i] = rhov*rhoinv*rhou_dot   + rhou*rhoinv*rhov_dot + rhou*rhov*rhoinv_dot
    A0_dot[3,3,i] = 2*rhov*rhoinv*rhov_dot + rhov*rhov*rhoinv_dot + p_dot_i
    A0_dot[4,3,i] = rhov*rhoinv*rhow_dot   + rhow*rhoinv*rhov_dot + rhov*rhow*rhoinv_dot
    A0_dot[5,3,i] = h*rhov_dot + rhov*h_dot

    A0_dot[1,4,i] = rhow_dot
    A0_dot[2,4,i] = rhow*rhoinv*rhou_dot   + rhou*rhoinv*rhow_dot + rhou*rhow*rhoinv_dot
    A0_dot[3,4,i] = rhow*rhoinv*rhov_dot   + rhov*rhoinv*rhow_dot + rhow*rhov*rhoinv_dot
    A0_dot[4,4,i] = 2*rhow*rhoinv*rhow_dot + rhow*rhow*rhoinv_dot + p_dot_i
    A0_dot[5,4,i] = h*rhow_dot + rhow*h_dot


    A0_dot[1,5,i] = rhoe_dot
    A0_dot[2,5,i] = rhou*h_dot  + h*rhou_dot
    A0_dot[3,5,i] = rhov*h_dot  + h*rhov_dot
    A0_dot[4,5,i] = rhow*h_dot  + h*rhow_dot
    A0_dot[5,5,i] = h*h*rho_dot + rho*2*h*h_dot - (p*a2_dot + p_dot_i*a2)/gamma_1
  end

  return nothing
end



