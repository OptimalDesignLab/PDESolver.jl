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

  # this function provides the opportunity to update the krylov tolerances.
  # Currently we don't change the tolerances, because leaving the relative
  # tolerance in place seems to work ok
  if newton_data.use_inexact_nk
    norm_i = newton_data.res_norm_i
    norm_i_1 = newton_data.res_norm_i_1
    gamma = newton_data.krylov_gamma
    #reltol = newton_data.ls.reltol*(norm_i/norm_i_1)^gamma
    reltol = newton_data.ls.reltol
    setTolerances(newton_data.ls, reltol, -1, -1, -1)
    #println(BSTDOUT, "updating krylov reltol to ", reltol)
  end

  return nothing
end


#------------------------------------------------------------------------------
# Psuedo-Transient Continuation (aka. Implicit Euler)
#------------------------------------------------------------------------------

"""
  Reinitailize the underlying data.  This function should be called every
  time newtonInner() is entered
"""
function clearEulerConstants(ls::LinearSolver)

  lo = getInnerLO(ls.lo, NewtonLinearObject)
  clearEulerConstants(lo)

  if !(typeof(ls.pc) <: PCNone)
    pc = getInnerPC(ls.pc, NewtonLinearObject)
    clearEulerConstants(pc)
  end

  return nothing
end

function clearEulerConstants(lo::NewtonLinearObject)

  lo.idata.res_norm_i = 0
  lo.idata.res_norm_i_1 = 0

  return nothing
end

"""
  Record the most recent residual norm

  **Inputs**

   * ls: a [`LinearSolver`](@ref) object
   * res_norm_i: must be a real number
"""
function recordEulerResidual(ls::LinearSolver, res_norm_i)

  lo = getInnerLO(ls.lo, NewtonLinearObject)
  lo.idata.res_norm_i = res_norm_i

  if !(typeof(ls.pc) <: PCNone)
    pc = getInnerPC(ls.pc, NewtonLinearObject)
    pc.idata.res_norm_i = res_norm_i
  end
end


"""
  Get the most recent residual norm and the residual norm the last time
  the PC was recalculated.  This function should be used to inspect the
  values only.  [`useEulerConstants`](@ref) should be used when updating the
  diagonal term in the Jacobian.

  **Inputs**

   * pc or lo: a PC or LO object

  **Outputs**

   * res_norm_i_1: second most recent residual norm
   * res_norm_i: most recent residual norm
"""
function getEulerConstants(pc::AbstractPC)

  pc2 = getInnerPC(pc, NewtonLinearObject)
  return pc2.idata.res_norm_i_1, pc2.idata.res_norm_i
end


function getEulerConstants(lo::AbstractLO)

  lo = getInnerLO(lo, NewtonLinearObject)
  return lo.idata.res_norm_i_1, lo.idata.res_norm_i
end


"""
  Similar to [`getEulerConstants`](@ref), this function both returns the
  Euler constants and shifts them so that the next relative residual will
  be correct.
"""
function useEulerConstants(pc::AbstractPC)
  pc2 = getInnerPC(pc, NewtonLinearObject)
  t1 = pc2.idata.res_norm_i_1
  t2 = pc2.idata.res_norm_i
  pc2.idata.res_norm_i_1 = t2
  return t1, t2
end

function useEulerConstants(lo::AbstractLO)
  lo2 = getInnerLO(lo, NewtonLinearObject)
  t1 = lo2.idata.res_norm_i_1
  t2 = lo2.idata.res_norm_i
  lo2.idata.res_norm_i_1 = t2
  return t1, t2
end

"""
  Check if Euler globalization is initialized
"""
function isEulerInitialized(lo::AbstractLO)
  lo2 = getInnerLO(lo, NewtonLinearObject)
  return !(lo2.idata.res_norm_i == 0 && lo2.idata.res_norm_i_1 == 0)
end

function isEulerInitialized(pc::AbstractPC)
  pc2 = getInnerPC(pc, NewtonLinearObject)
  return !(pc2.idata.res_norm_i == 0 && pc2.idata.res_norm_i_1 == 0)
end

function isEulerInitialized(pc::PCNone)
  return false  # can't globalize what doesn't exist
end




@doc """
### NonlinearSolver.calcTauVec

  This function calculates the spatially varying vector for the timestep.

  The timestep varies according to tau/(1 + |J|^(1/d) ), where `|J|` is the
  determinant of the mapping jacobian and `d` is the dimensionality (2D ord
  3D).

  **Inputs**

  *  mesh
  *  opts
  *  tau:  timestep factor

  **Inputs/Outputs**

   * tau_vec: vector (of length numDof) populated with tau*spatially varying
             factor

  Aliasing restrictions:  none
"""->
function calcTauVec(mesh, opts, tau, tau_vec)
# calculate the spatially varying pseudo-timestep
  #TODO: make tau_vec = 1/tau_vec, so we don't have to do fp division when
  #      applying it
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        dof = mesh.dofs[k, j, i]
        tau_vec[dof] = tau/(1 + sqrt(real(mesh.jac[j, i])))  #TODO: should be |J|^1/dim
#        tau_vec[dof] = tau
      end
    end
  end

  if mesh.myrank == 0
    println(BSTDOUT, "average tau value = ", mean(tau_vec))
  end

  return nothing

end


"""
  Reinitilize an `ImplicitEulerData` object (for example, when calling
  `newtonInner` a second time).

  Uses opts["euler_tau"] to reset the pseudo-time step

  **Inputs**

   * mesh
   * opts
   * idata: `ImplicitEulerData`
"""
function reinitImplicitEuler(mesh, opts, idata::ImplicitEulerData)

  idata.res_norm_i = 0
  idata.res_norm_i_1 = 0
  idata.tau_l = opts["euler_tau"]
  calcTauVec(mesh, opts, idata.tau_l, idata.tau_vec)

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

  tau_l_old = lo.idata.tau_l
  res_norm_i_1, res_norm_i = useEulerConstants(lo)

  # on the first Jacobian calculate, don't update
  if res_norm_i_1 == 0
    return nothing
  end
  
  tau_update = res_norm_i_1/res_norm_i

  # update tau
  lo.idata.tau_l = lo.idata.tau_l * tau_update
  
  if lo.myrank == 0
    println(BSTDOUT, "tau_update factor = ", tau_update)
  end

  scale!(lo.idata.tau_vec, tau_update)

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

  if !isEulerInitialized(lo)
    return nothing
  end

  if mesh.myrank == 0
    println(BSTDOUT, "applying Implicit Euler globalization")
    println(BSTDOUT, "average tau value = ", mean(lo.idata.tau_vec))
  end

  
  lo2 = getBaseObject(lo)
#  println("euler globalization tau = ", lo.tau_l)
  # create the indices

  val = [1/lo.idata.tau_l]
  idx = PetscInt[0]
  idy = PetscInt[0]
  #TODO: replace this with MatDiagonalSet
  for i=1:mesh.numDof
    idx[1] = i + mesh.dof_offset
    idy[1] = i + mesh.dof_offset
    val[1] = -eqn.M[i]/lo.idata.tau_vec[i]
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
                    lo::NewtonLinearObject, b::AbstractArray)
# apply the diagonal update term to the jacobian vector product

  if !isEulerInitialized(lo)
    return nothing
  end


  for i=1:mesh.numDof
    b[i] -= eqn.M[i]*(1/lo.idata.tau_vec[i])*vec[i]
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
