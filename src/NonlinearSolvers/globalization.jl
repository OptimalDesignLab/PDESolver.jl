# globalization.jl
# this file contains the methods for applying differnent globalizations
# to Newtons method


#------------------------------------------------------------------------------
# Inexact-Newton-Krylov
#------------------------------------------------------------------------------
@doc """
### NonlinearSolvers.updateKrylov

  This function up dates the relative toleranced for the linear solve when doing
  inexact Newton-Krylov.

  Inputs/Outputs:
    newton_data: the NewtonData object

"""->
function updateKrylov(newton_data::NewtonData)

  norm_i = newton_data.res_norm_i
  norm_i_1 = newton_data.res_norm_i_1
  gamma = newton_data.krylov_gamma
  reltol = newton_data.ls.reltol*(norm_i/norm_i_1)^gamma
  setTolerances(newton_data.ls, reltol, -1, -1, -1)
#  println("updating krylov reltol to ", reltol)

  return nothing
end


#------------------------------------------------------------------------------
# Psuedo-Transient Continuation (aka. Implicit Euler)
#------------------------------------------------------------------------------

# res_norm_i_1 and res_norm_i, updated inside newton's method
# This is a bad way to pass parameters to updateEuler(), but I can't find a
# better way that retains the composability of the linear operator abstraction.
global const EulerConstants = Float64[0, 0]

"""
  Reinitailize the underlying data.  This function should be called every
  time newtonInner() is entered
"""
function clearEulerConstants()
  fill!(EulerConstants, 0.0)
end

"""
  Record the most recent residual norm

  **Inputs**

   * res_norm_i: must be a real number
"""
function recordEulerResidual(res_norm_i)
  EulerConstants[1] = EulerConstants[2]
  EulerConstants[2] = res_norm_i
end

"""
  Get the two most recent residual norms

  **Outputs**

   * res_norm_i_1: second most recent residual norm
   * res_norm_i: most recent residual norm
"""
function getEulerConstants()
  return EulerConstants[1], EulerConstants[2]
end

"""
  Check if Euler globalization is initialized
"""
function isEulerInitialized()
  if EulerConstants[1] == EulerConstants[2] == 0
    return false
  else
    return true
  end
end


@doc """
### NonlinearSolvers.initEuler

  This function initializes the data needed to do Psudo-Transient Continuation 
  globalization (aka. Implicit Euler) of Newton's method, using a spatially 
  varying pseudo-timestep.

  Updates the jacobian with a diagonal term, as though the jac was the 
  jacobian of this function:
  (u - u_i_1)/delta_t + f(u)
  where f is the original residual and u_i_1 is the previous step solution

  The timestep varies according to tau/(1 + sqrt(det(jac))).

  This globalization is activated using the option `newton_globalize_euler`.
  The initial value of the scaling factor tau is specified by the option 
  `euler_tau`.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Outputs:
    tau_l: the timestep factor
    tau_vec: the vector (of length numDof) that is added to the diagonal 
             of the jacobian.

"""->
function initEuler(mesh, sbp, eqn, opts)

  tau_l = opts["euler_tau"]  # initailize tau to something
  tau_vec = zeros(mesh.numDof)
  calcTauVec(mesh, sbp, eqn, opts, tau_l, tau_vec)

  return tau_l, tau_vec
end

@doc """
### NonlinearSolver.calcTauVec

  This function calculates the spatially varying vector for the timestep.

  Inputs:
    mesh
    sbp
    eqn
    opts
    tau:  timestep factor

  Inputs/Outputs: 
    tau_vec: vector (of length numDof) populated with tau*spatially varying
             factor

  Aliasing restrictions:  none
"""->
function calcTauVec(mesh, sbp, eqn, opts, tau, tau_vec)
# calculate the spatially varying pseudo-timestep
  #TODO: make tau_vec = 1/tau_vec, so we don't have to do fp division when
  #      applying it
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        dof = mesh.dofs[k, j, i]
        tau_vec[dof] = tau/(1 + sqrt(real(mesh.jac[j, i])))
#        tau_vec[dof] = tau
      end
    end
  end

  println(BSTDOUT, "average tau value = ", mean(tau_vec))

  return nothing

end


@doc """
### NonlinearSolver.updateEuler

  This function updates the tau paraemter for the implicit Euler
  globalization.  tau_vec is also updated

  Inputs/Outputs:
    lo: any kind of Newton PC or LO object

"""->
function updateEuler(lo::NewtonLinearObject)
  # updates the tau parameter for the Implicit Euler globalization
  # norm_i is the residual step norm, norm_i_1 is the previous residual norm

  tau_l_old = lo.tau_l
  res_norm_i_1, res_norm_i = getEulerConstants()

  # on the first Jacobian calculate, don't update
  if res_norm_i_1 == 0
    return nothing
  end

  # update tau
  lo.tau_l = lo.tau_l * res_norm_i_1/res_norm_i
  
  tau_update = lo.tau_l/tau_l_old
  println(BSTDOUT, "tau_update factor = ", tau_update)
  for i=1:length(lo.tau_vec)
    lo.tau_vec[i] *= tau_update
  end

  return nothing
end

@doc """
### NonlinearSolvers.applyEuler

  This function updates the jacobian matrix with the term from the implicit 
  Euler globalization.  The term is eqn.M/tau_vec.  Methods are available for
  dense, sparse, Petsc jacobians, as well as jacobian-vector products.

 ** Inputs**
    mesh
    sbp
    eqn
    opts

  **Inputs/Outputs**
    lo:  a [`NewtonHasMat`](@ref) t contaiing tau_vec and the Jacobian matrix

"""->
function applyEuler(mesh, sbp, eqn, opts, lo::NewtonHasMat)
# this allocations memory every time
# should there be a reusable array for this?
# maybe something in lo?
# for explicitly stored jacobian only

  if !isEulerInitialized()
    return nothing
  end

  println(BSTDOUT, "applying Implicit Euler globalization")
  println(BSTDOUT, "average tau value = ", mean(lo.tau_vec))

  lo2 = getBaseLO(lo)
#  println("euler globalization tau = ", lo.tau_l)
  # create the indices

  val = [1/lo.tau_l]
  idx = PetscInt[0]
  idy = PetscInt[0]
  #TODO: replace this with MatDiagonalSet
  for i=1:mesh.numDof
    idx[1] = i + mesh.dof_offset
    idy[1] = i + mesh.dof_offset
    val[1] = -eqn.M[i]/lo.tau_vec[i]
    set_values1!(lo2.A, idx, idy, val, ADD_VALUES)
#    PetscMatSetValues(jac, idx, idy, val, ADD_VALUES)
  end


  return nothing
end

"""
  Updates a jacobian-vector product with the effects of the globalization

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * vec: vector that the jacobian is being multiplied by
   * lo: a linaer operator or preconditioner, usually a matrix-free one

  **Inputs/Outputs**

   * b: result vector, updated
"""
function applyEuler(mesh, sbp, eqn, opts, vec::AbstractArray, 
                    lo, b::AbstractArray)
# apply the diagonal update term to the jacobian vector product

  if !isEulerInitialized()
    return nothing
  end


  for i=1:mesh.numDof
    b[i] -= eqn.M[i]*(1/lo.tau_vec[i])*vec[i]
  end

  return nothing
end


# Globalize by adding a diagonal term to the Jacobian.  Rather crude.
#------------------------------------------------------------------------------
  
function addDiagonal(mesh, sbp, eqn, jac)
# add the mass matrix to the jacobian

  for i=1:mesh.numDof
    idx = PetscInt[i-1]
    idy = PetscInt[i-1]
    vals = [100*eqn.M[i]]

#     println("adding ", vals, " to jacobian entry ", i, ",", i)
    PetscMatSetValues(jac, idx, idy, vals, ADD_VALUES)
  end

  return nothing

end


#------------------------------------------------------------------------------
