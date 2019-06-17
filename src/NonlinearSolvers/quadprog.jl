
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

  if eigs_to_remove == "neg"
    scale!(Jac, -1.0)
  elseif eigs_to_remove == "pos"
    # do nothing
  elseif eigs_to_remove == "none"
    return
  else
    error("eigs_to_remove specified incorrectly.")
  end

  
  n = size(Jac,1)
  # compute baseline product, 0.5*u.'*(Jac^T + Jac)*u
  prod = zero(T)
  for i = 1:n
    for j = 1:n
      prod += 0.5*(Jac[i,j] + Jac[j,i])*u[i]*u[j]
    end
  end

  #TODO TODO: prod < 0 for eigs_to_remove == "neg"???
  if prod > 0
  # if prod < 0
    # nothing to do
    # println("prod > 0 check hit, not stabilizing")
    return
  end

  # println("prod <= 0, now stabilizing")

  # array A stores the entries in the contraint Jacobian
  # A = zeros(div(n*(n+1),2))
  for i = 1:n
    A[div(i*(i-1),2)+i] = u[i]*u[i]
    for j = 1:(i-1)
      A[div(i*(i-1),2)+j] = 2.0*u[i]*u[j]
    end
  end

  # A *= -prod/dot(A,A)
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

  if eigs_to_remove == "neg"    # TODO ???
    scale!(Jac, -1.0)
  end

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
