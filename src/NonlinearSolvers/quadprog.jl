
"""
  returns matrix `Jacpert` such that `u.'*sym(Jac + Jacpert)*u` is strictly
  positive

  **Inputs**

   * `u`: an given vector whose dimensions are consistent with `Jac`
   * `Jac`: matrix that needs to be perturbed
   * `A`: a work vector needed by this function (overwritten).  The
          element type should be the "maximum" type of the element types
          of `u` and `Jac`


  **Inputs/Outputs**

   * `Jacpert`: matrix perturbation
"""
function findStablePerturbation!(Jac::AbstractMatrix,
                                 u::AbstractVector,
                                 A::AbstractVector{T},
                                 eigs_to_remove::String) where T

  @assert( size(Jac,1) == size(Jac,2) == length(u) )

  n = size(Jac,1)
  # compute baseline product, 0.5*u.'*(Jac^T + Jac)*u
  # This is b in the derivation (without the negative sign)
  # Note that 0.5*(Jac^T + Jac) is J_sym.
  prod = zero(T)
  for i = 1:n
    for j = 1:n
      prod += 0.5*(Jac[i,j] + Jac[j,i])*u[i]*u[j]
    end
  end

  if prod < 0
    # nothing to do
    # println("prod > 0 check hit, not stabilizing")
    stab_term_vecnorm = 0.0
    return stab_term_vecnorm
  end

  # println("prod <= 0, now stabilizing")

  # Form the KKT constraint Jacobian in A
  # size: A = zeros(div(n*(n+1),2))
  # A is a row matrix. It looks like this:
  #   [u1^2, 2*u2*u1, u2^2, 2*u3*u1, 2*u3*u2, u3^2, ...]
  # This comes from the fact that we are doing a u^T*Jac*u product,
  #   resulting in a scalar. Since Jac is symmetric, the off-diagonal
  #   entries are multiplied by 2.
  for i = 1:n
    A[div(i*(i-1),2)+i] = u[i]*u[i]
    for j = 1:(i-1)
      A[div(i*(i-1),2)+j] = 2.0*u[i]*u[j]
    end
  end

  # The stabilization matrix is symmetric. So we only need to consider
  #   one triangle of it (here, lower triangle w/ diagonal).
  # This lower triangle, in vector form, is given by
  #   s = -((u^T*Jac*u)/(A*A^T)) * A^T.
  # The term inside the parentheses is a scalar. In other words, one entry is
  #   s_1 = -((u^T*Jac*u)/(A*A^T)) * A^T_1.
  # Here, we are scaling A by that scalar factor so that it becomes s.
  # Note: dot(A, A) == A*A^T

  # old: A *= -prod/dot(A,A)
  scale!(A, -prod/dot(A, A))        # divide by zero! root of NaN.

  # Form the stabilization matrix from its entries on the lower triangle,
  #   now stored in A. This stabilization matrix is what needs to be returned.
  # We are reusing the Jac matrix that was passed in (element strong volume Jacobian),
  #   because that is what the assembly routine in stabilizeCNDSLO assembles into the
  #   full Jacobian.
  # It is to be ADDED to the full Jacobian, unlike the clipJac method, which is subtracted
  #   by scaling the stabilization contribution.
  # Historical note: We used to have '+='s below, not '='s. This was because this function
  #   performed the assembly into the Lorenz Jacobian. 
  #   Now assembly is handled in the calling function.
  for i = 1:n
    # diagonal elements of the stab matrix
    Jac[i,i] = A[div(i*(i-1),2)+i]
    for j = 1:(i-1)
      # off diagonal elements of the stab matrix
      Jac[i,j] = A[div(i*(i-1),2)+j]
      Jac[j,i] = A[div(i*(i-1),2)+j]
    end
  end

  stab_term_vecnorm = vecnorm(Jac)

  return stab_term_vecnorm

end     # end function findStablePerturbation!

"""
  modifies matrix `Jac` such that `u.'*sym(Jac)*u` is strictly positive

  **Inputs**

   * `u`: an given vector whose dimensions are consistent with `Jac`
   * `A`: a work vector needed by this function (overwritten).  The
          element type should be the "maximum" type of the element types
          of `u` and `Jac`

  **Inputs/Outputs**

   * `Jac`: matrix being modified
"""
function removeUnstableModes!(Jac::AbstractMatrix,
                              u::AbstractVector,
                              A::AbstractVector{T}) where T
  @assert( size(Jac,1) == size(Jac,2) == length(u) )
  n = size(Jac,1)
  # compute baseline product, 0.5*u.'*(Jac^T + Jac)*u
  prod = zero(T)
  for i = 1:n
    for j = 1:n
      prod += 0.5*(Jac[i,j] + Jac[j,i])*u[i]*u[j]
    end
  end
  if prod > 0
    # nothing to do
    return
  end
  # array A stores the entries in the contraint Jacobian
#  A = zeros(div(n*(n+1),2))
  for i = 1:n
    A[div(i*(i-1),2)+i] = u[i]*u[i]
    for j = 1:(i-1)
      A[div(i*(i-1),2)+j] = 2.0*u[i]*u[j]
    end
  end
  scale!(A, -prod/dot(A,A))         # flattening lower triangular into vector
#  A *= -prod/dot(A,A)
  for i = 1:n
    Jac[i,i] += A[div(i*(i-1),2)+i]       # diagonal entries
    for j = 1:(i-1)
      Jac[i,j] += A[div(i*(i-1),2)+j]     # off-diagonal entries
      Jac[j,i] += A[div(i*(i-1),2)+j]
    end
  end

  return nothing
end

"""
  returns matrix `Jacpert` such that `u.'*sym(Jac + Jacpert)*u` is strictly
  positive

  **Inputs**

   * `u`: an given vector whose dimensions are consistent with `Jac`
   * `Jac`: matrix that needs to be perturbed
   * `A`: a work vector needed by this function (overwritten).  The
          element type should be the "maximum" type of the element types
          of `u` and `Jac`


  **Inputs/Outputs**

   * `Jacpert`: matrix perturbation
"""
function findStablePerturbation_explicit!(Jac::AbstractMatrix,
                                 u::AbstractVector,
                                 A::AbstractVector{T},
                                 eigs_to_remove::String) where T

# function findStablePerturbation!(Jac::AbstractMatrix,
                                 # u::AbstractVector,
                                 # A::AbstractVector{T},
                                 # Jacpert::AbstractMatrix) where T
  @assert( size(Jac,1) == size(Jac,2) == length(u) )
  # @assert( size(Jac,1) == size(Jac,2) == size(Jacpert,1) == size(Jacpert,2)
           # == length(u) )
  # println(":::::::::::: entering findStablePert")

  scale_u = 1e100
  # scale_u = 1e60
  # scale!(u, scale_u)
  
  u_check = dot(u,u)

  if u_check < 1e-12     # use norm instead of dot?
    # println("u sufficiently small, after scale. u_check: $u_check. u: ", u')
    global STAB_ctr_usmallandnotstabing
    STAB_ctr_usmallandnotstabing = STAB_ctr_usmallandnotstabing+1

    scale!(u, 1/scale_u)
    return
  end


  n = size(Jac,1)
  # compute baseline product, 0.5*u.'*(Jac^T + Jac)*u
  prod = zero(T)
  for i = 1:n
    for j = 1:n
      prod += 0.5*(Jac[i,j] + Jac[j,i])*u[i]*u[j]
      # println(" calcing prod. i = $i, j = $j, prod = $prod")
    end
  end

  #TODO TODO: prod < 0 for eigs_to_remove == "neg"???
  if prod > 0
    # nothing to do
    # println("prod > 0 check hit, not stabilizing")
    global STAB_ctr_prodgt0andnotstabing
    STAB_ctr_prodgt0andnotstabing = STAB_ctr_prodgt0andnotstabing+1

    scale!(u, 1/scale_u)
    return
  end

  # println("prod <= 0, now stabilizing")
  global STAB_ctr_prodle0andstabing
  STAB_ctr_prodle0andstabing = STAB_ctr_prodle0andstabing+1
  # array A stores the entries in the contraint Jacobian
#  A = zeros(div(n*(n+1),2))

  for i = 1:n
    A[div(i*(i-1),2)+i] = u[i]*u[i]
    for j = 1:(i-1)
      A[div(i*(i-1),2)+j] = 2.0*u[i]*u[j]
    end
  end

#  A *= -prod/dot(A,A)
  scale!(A, -prod/dot(A, A))        # divide by zero! root of NaN.

  # fill!(Jacpert, 0.0)

  #=
  for i = 1:n
    Jacpert[i,i] += A[div(i*(i-1),2)+i]
    for j = 1:(i-1)
      Jacpert[i,j] += A[div(i*(i-1),2)+j]
      Jacpert[j,i] += A[div(i*(i-1),2)+j]
    end
  end
  =#

  for i = 1:n
    Jac[i,i] += A[div(i*(i-1),2)+i]
    for j = 1:(i-1)
      Jac[i,j] += A[div(i*(i-1),2)+j]
      Jac[j,i] += A[div(i*(i-1),2)+j]
    end
  end


  # JEH meeting 20181004: Don't need to scale Jac here. u is already scaled??? Think about it.

  # scale!(Jac, (1/scale_u)^4)      # need to unscale Jacpert by the amount used to scale u, to the 4th power
  # scale!(Jacpert, (1/scale_u)^4)      # need to unscale Jacpert by the amount used to scale u, to the 4th power
                                      # this comes from:
                                      #   A = f(u^2)
                                      #   prod = f(u^2), and A scaled by prod

  scale!(u, 1/scale_u)

end
