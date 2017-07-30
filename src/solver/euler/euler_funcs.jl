# euler_funcs.jl
# this file contains all the functions that calculate values at a node (except
# for boundary related functions) as well as some additional helper functions.

#------------------------------------------------------------------------------
# function to calculate the Euler flux over the entire mesh
#------------------------------------------------------------------------------
@doc """
### EulerEquationMod.getEulerFlux

  This function calculates the Euler flux across the entire mesh by passing
  pieces of the eqn.q, eqn.aux_vars, eqn.f_xi and eqn.params to a low level
  function.  The flux is calculated in the xi and eta directions,
  scaled (mulitiplied) by the mapping jacobian (so that when performing the
  integral we don't have to explictly divide by the jacobian, it just cancels
  out with the jacobian factor introduced here.

  Calls writeFlux to do any requested output.

  This is a mid level function
"""->
# mid level function
function getEulerFlux{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                        sbp::AbstractSBP,
                                        eqn::EulerData{Tsol, Tres, Tdim}, opts)
# calculate Euler flux in parametric coordinate directions, stores it in eqn.flux_parametric

  nrm = zeros(Tmsh, Tdim)
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement  # loop over nodes on current element
      q_vals = ro_sview(eqn.q, :, j, i)
      aux_vars = ro_sview(eqn.aux_vars, :, j, i)
      # put this loop here (rather than outside) to we don't have to fetch
      # q_vals twice, even though writing to the flux vector is slower
      # it might be worth copying the normal vector rather than
      # doing an view

      for k=1:Tdim  # loop over dimensions
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        # don't do an array view because strided views are type-unstable
#        nrm[1] = mesh.dxidx[k, 1, j, i]
#        nrm[2] = mesh.dxidx[k, 2, j, i]
        flux = sview(eqn.flux_parametric, :, j, i, k)

      	# this will dispatch to the proper calcEulerFlux
        calcEulerFlux(eqn.params, q_vals, aux_vars, nrm, flux)
      end

    end
  end


  writeFlux(mesh, sbp, eqn, opts)

  return nothing
end


@doc """
### EulerEquationMod.writeFlux

  This function writes the real part of Euler flux to a file named Fxi.dat,
  space delimited, controlled by the input options 'writeflux', of type Bool.

  This is a high level function.
"""->
function writeFlux(mesh, sbp, eqn, opts)

   if !eqn.params.writeflux
     return nothing
   end

   fname = "Fxi.dat"
   rmfile(fname)
   writedlm(fname, real(eqn.flux_parametric))

   return nothing
end


@doc """
### EulerEquationMod.getEulerFlux2

  This function calcules the euler flux over the entire mesh directly (ie.
  does not call a low level function.  This function is deprecated, although
  useful for benchmarking purposes.  2D only.

  This is a mid level function
"""->
# this function is deprecated in favor of getEulerFlux()
# useful for benchmarking purposes
function getEulerFlux2{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                                   eqn::EulerData{Tsol}, opts)
# calculates the Euler flux for every node in the xi and eta directions
# eqn is the equation type
# q is the 3D array (4 by nnodes per element by nel), of the conservative variables
# dxidx is the 4D array (2 by 2 x nnodes per element by nel) that specifies the direction of xi and eta in each element (output from mappingjacobian!)
# flux_parametric is populated with the flux in xi direction (same shape as q)
# F_eta is populated with flux in eta direction

# once the Julia developers fix slice notation and speed up subarrays, we won't have to
# vectorize like this (can calculate flux one node at a time inside a dedicated function

q = eqn.q
dxidx = mesh.dxidx
flux_parametric = sview(eqn.flux_parametric, :, :, :, 1)
F_eta = sview(eqn.flux_parametric, :, :, :, 2)

(ncomp, nnodes, nel) = size(q)  # get sizes of things

  for i=1:nel  # loop over elements
    for j=1:nnodes  # loop over nodes within element
      # get direction vector components (xi direction)
      nx = dxidx[1, 1, j, i]
      ny = dxidx[1, 2, j, i]
      # calculate pressure
      press = (eqn.params.gamma-1)*(q[4, j, i] - 0.5*(q[2, j, i]^2 + q[3, j, i]^2)/q[1, j, i])

      # calculate flux in xi direction
      # hopefully elements of q get stored in a register for reuse in eta direction
      U = (q[2, j, i]*nx + q[3, j, i]*ny)/q[1, j, i]
      flux_parametric[1, j, i] = q[1, j, i]*U
      flux_parametric[2, j, i] = q[2, j, i]*U + nx*press
      flux_parametric[3, j, i] = q[3, j, i]*U + ny*press
      flux_parametric[4, j, i] = (q[4, j, i] + press)*U

      # get direction vector components (eta direction)
      nx = dxidx[2, 1, j, i]
      ny = dxidx[2, 2, j, i]

      # calculate xi flux
      U = (q[2, j, i]*nx + q[3, j, i]*ny)/q[1, j, i]
      F_eta[1, j, i] = q[1, j, i]*U
      F_eta[2, j, i] = q[2, j, i]*U + nx*press
      F_eta[3, j, i] = q[3, j, i]*U + ny*press
      F_eta[4, j, i] = (q[4, j, i] + press)*U
    end
  end



  return nothing

end

"""
  Calculates the volume integrals for the weak form, computing the Euler flux
  as needed, rather than using eqn.flux_parametric

  Inputs:
    mesh
    sbp
    eqn: eqn.res is updated with the result
    opts
"""
function calcVolumeIntegrals_nopre{Tmsh, Tsol, Tres, Tdim}(
                                   mesh::AbstractMesh{Tmsh},
                                   sbp::AbstractSBP,
                                   eqn::EulerData{Tsol, Tres, Tdim},
                                   opts)


  # flux in the parametric directions for a given element
  flux_el = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, Tdim)

  nrm = eqn.params.nrm  # vector in parametric direction

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      for k=1:Tdim
        flux_k = sview(flux_el, :, j, k)

        # get the direction vector
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        # consider calculating all directions at once
        # not sure if that will help because the data dependencies are 
        # really simple
        calcEulerFlux(eqn.params, q_j, aux_vars_j, nrm, flux_k)
      end  # end loop k
    end  # end loop j

    res_i = sview(eqn.res, :, :, i)
    for k=1:Tdim
      weakDifferentiateElement!(sbp, k, sview(flux_el, :, :, k), res_i, SummationByParts.Add(), true)
    end

  end  # end loop i

  return nothing
end  # end function





"""
  Calculate (S .*F)1 where S is the skew symmetric part of sbp.Q and F
  is a symmetric numerical flux function.  eqn.res is updated with the result.
  Methods are available for curvilinear and non-curvilinear meshes

  Inputs:
    mesh
    sbp
    eqn
    opts
    functor: the numerical flux function F, of type FluxType
"""
function calcVolumeIntegralsSplitForm{Tmsh, Tsol, Tres, Tdim}(
                                        mesh::AbstractMesh{Tmsh}, 
                                        sbp::AbstractSBP,  
                                        eqn::EulerData{Tsol, Tres, Tdim}, opts,
                                        functor::FluxType)

  if opts["use_staggered_grid"]
    # not planning on implementing the non-curvilinear version of this
    calcVolumeIntegralsSplitFormStaggered(mesh, mesh,mesh2, sbp, mesh.sbp2,
                                          eqn, opts, functor)
  else
    if mesh.coord_order == 1
      calcVolumeIntegralsSplitFormLinear(mesh, sbp, eqn, opts, functor)
    else
      calcVolumeIntegralsSplitFormCurvilinear(mesh, sbp, eqn, opts, functor)
    end
  end

  return nothing
end

@doc """
  Calculate (S .* F)1, where S is the skew-symmetric part of sbp.Q 
  and F is a symmetric numerical flux function.  eqn.res is updated 
  with the result.  Linear (non-curvilinear) meshes only
"""
function calcVolumeIntegralsSplitFormLinear{Tmsh, Tsol, Tres, Tdim}(
                                        mesh::AbstractMesh{Tmsh}, 
                                        sbp::AbstractSBP,  
                                        eqn::EulerData{Tsol, Tres, Tdim}, opts,
                                        functor::FluxType)

#  println("----- entered calcVolumeIntegralsSplitForm -----")
  dxidx = mesh.dxidx
  res = eqn.res
  q = eqn.q
  nrm = eqn.params.nrmD
  aux_vars = eqn.aux_vars
  F_d = eqn.params.flux_valsD
  S = eqn.params.S
  params = eqn.params
  
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q, :, j, i)
      aux_vars_j = ro_sview(aux_vars, :, j, i)
      for k=1:(j-1)  # loop over lower triangle of S
        q_k = ro_sview(q, :, k, i)
        # calcaulate the normal vector in each parametric directions
        for d=1:Tdim
          # get the normal vector
          for p=1:Tdim
            nrm[p, d] = dxidx[d, p, j, i] 
          end
        end

        # calculate the numerical flux functions in all Tdim
        # directions at once
        functor(params, q_j, q_k, aux_vars_j, nrm, F_d)

        @simd for d=1:Tdim
          # update residual
          @simd for p=1:(Tdim+2)
            res[p, j, i] -= 2*S[j, k, d]*F_d[p, d]
            res[p, k, i] += 2*S[j, k, d]*F_d[p, d]
          end

        end  # end d loop
      end  # end k loop


    end  # end j loop
  end  # end i loop

  return nothing
end

@doc """
  Calculate (S .* F)1, where S is the skew-symmetric part of sbp.Q 
  and F is a symmetric numerical flux function.  eqn.res is updated 
  with the result.  This function is used for curvilinear meshes.
"""
function calcVolumeIntegralsSplitFormCurvilinear{Tmsh, Tsol, Tres, Tdim}(
                                        mesh::AbstractMesh{Tmsh}, 
                                        sbp::AbstractSBP,  
                                        eqn::EulerData{Tsol, Tres, Tdim}, opts,
                                        functor::FluxType)

#  println("----- entered calcVolumeIntegralsSplitForm -----")
  dxidx = mesh.dxidx
  res = eqn.res
  q = eqn.q
  nrm = eqn.params.nrmD
  aux_vars = eqn.aux_vars
  F_d = eqn.params.flux_valsD
#  S = eqn.params.S
  S = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)
  params = eqn.params

  # S is calculated in x-y-z, so the normal vectors should be the unit normals
  fill!(nrm, 0.0)
  for d=1:Tdim
    nrm[d, d] = 1
  end
  
  for i=1:mesh.numEl
    # get S for this element
    dxidx_i = ro_sview(dxidx, :, :, :, i)
    calcSCurvilinear(sbp, dxidx_i, S)

    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q, :, j, i)
      aux_vars_j = ro_sview(aux_vars, :, j, i)
      for k=1:(j-1)  # loop over lower triangle of S
        q_k = ro_sview(q, :, k, i)

        # calculate the numerical flux functions in all Tdim
        # directions at once
        functor(params, q_j, q_k, aux_vars_j, nrm, F_d)

        @simd for d=1:Tdim
          # update residual
          @simd for p=1:(Tdim+2)
            res[p, j, i] -= 2*S[j, k, d]*F_d[p, d]
            res[p, k, i] += 2*S[j, k, d]*F_d[p, d]
          end

        end  # end d loop
      end  # end k loop


    end  # end j loop
  end  # end i loop

  return nothing
end


"""
  This function calculates I_S2F.'*(S .* F( uk(I_S2F wk), uk(I_S2F wk)))1.
  This is similar to
  [`calcVolumeIntegralsSplitFormCurvilinear](@ref), but for the staggered grid
  algorithm.

  Inputs:
    mesh_s: the mesh object for the solution grid
    mesh_f: the mesh object for the flux grid
    sbp_s: SBP operator for solution grid
    sbp_f: SBP operator for flux grid
    functor: the numerical flux functor ([`FluxType`](@ref)) used to compute F

  Inputs/Outputs:
    eqn: the equation object (implicitly on the solution grid).  eqn.res is
         updated with the result

  Aliasing Restrictions: none (the meshes and the sbp operators could alias,
                         in which case this algorithm reduces to the 
                         non-staggered version
"""
function calcVolumeIntegralsSplitFormCurvilinear{Tmsh, Tsol, Tres, Tdim}(
                                        mesh_s::AbstractMesh{Tmsh},
                                        mesh_f::AbstractMesh{Tmsh},
                                        sbp_s::AbstractSBP,
                                        sbp_f::AbstractSBP
                                        eqn::EulerData{Tsol, Tres, Tdim}, opts,
                                        functor::FluxType)


  error("calcVolumeIntegralsSplitForm not implemented for staggered grid")

  return nothing
end


# calculating the Euler flux at a node
#------------------------------------------------------------------------------
@doc """
### EulerEquationMod.calcEulerFlux

   This function calculates the Euler flux from the conservative variables at
   a single node in a particular direction.  2D only.

   Inputs:
   params  : ParamaterType{2, :conservative}
   q  : vector of conservative variables
   aux_vars : vector of auxiliary variables
   dir :  vector in direction to calculate the flux

   Inputs/Outputs:
   F  : vector to populate with the flux

   The Tdim paramater of params determine whether this method or the 3D
   version is called.

   This is a low level function
"""->
# low level function
function calcEulerFlux{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# aux_vars is the vector of auxiliary variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux (is a vector of length 4)
# 2D  only


  press = calcPressure(params, q)
#  press = getPressure(aux_vars)
#  press = @getPressure(aux_vars)
  U = (q[2]*dir[1] + q[3]*dir[2])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = (q[4] + press)*U

  return nothing

end

@doc """
###EulerEquationMod.calcEulerFlux_revm

Compute the derivative of the euler flux in reverse mode w.r.t to unit vector
flux direction.

"""->
function calcEulerFlux_revm{Tmsh, Tsol}(params::ParamType{2, :conservative},
                            q::AbstractArray{Tsol,1}, aux_vars,
                            dir::AbstractArray{Tmsh,1}, F_bar, dir_bar)

  # Compute the reverse mode
  # Differentiate euler flux product with F_bar in reverse mode w.r.t dir to get
  # q_bar

  press = calcPressure(params, q)
  U = (q[2]*dir[1] + q[3]*dir[2])/q[1]

  # intermediate function that is only used in computing F so has to be reverse
  # diffed only in F_bar

  # dir_bar has dependence on both F and U
  # Reverse mode using F_bar
  dir_bar[1] += F_bar[2]*press
  dir_bar[2] += F_bar[3]*press
  U_bar = 0.0
  U_bar += F_bar[1]*q[1] + F_bar[2]*q[2] + F_bar[3]*q[3] + F_bar[4]*(q[4] + press)

  # Reverse mode using U_bar
  dir_bar[1] += U_bar*q[2]/q[1]
  dir_bar[2] += U_bar*q[3]/q[1]

  return nothing
end

@doc """
###EulerEquationMod.calcEulerFlux_revq

Compute the derivative of the Euler flux in reverse mode w.r.t q

"""->

function calcEulerFlux_revq{Tmsh, Tsol}(params::ParamType{2, :conservative},
                            q::AbstractArray{Tsol,1}, aux_vars,
                            dir::AbstractArray{Tmsh,1}, F_bar, q_bar)

  press = calcPressure(params, q)
  U = (q[2]*dir[1] + q[3]*dir[2])/q[1]

  U_bar = zero(Tsol)     # Initialize
  press_bar = zero(Tsol) #
  # Reverse diff F[4] = (q[4] + press)*U
  q_bar[4] += F_bar[4]*U
  U_bar += F_bar[4]*(q[4] + press)
  press_bar += F_bar[4]*U

  # Reverse diff F[3] = q[3]*U + dir[2]*press
  q_bar[3] += F_bar[3]*U
  U_bar += F_bar[3]*q[3]
  press_bar += F_bar[3]*dir[2]

  # Reverse diff F[2] =  q[2]*U + dir[1]*press
  q_bar[2] += F_bar[2]*U
  U_bar += F_bar[2]*q[2]
  press_bar += F_bar[2]*dir[1]

  # Reverse diff  F[1] = q[1]*U
  q_bar[1] += F_bar[1]*U
  U_bar += F_bar[1]*q[1]

  # Reverse diff U = (q[2]*dir[1] + q[3]*dir[2])/q[1]
  q_bar[2] += U_bar*dir[1]/q[1]
  q_bar[3] += U_bar*dir[2]/q[1]
  q_bar[1] -= U_bar*(q[2]*dir[1] + q[3]*dir[2])/(q[1]*q[1])

  # Reverse diff press = calcPressure(params, q)
  calcPressure_revq(params, q, press_bar, q_bar)

  return nothing
end


@doc """
# low level function
    Calculates the Euler flux from entropy variables

    Inputs:
    params : ParameterType{Tdim, :entropy}
    q : vector of entropy variables
    aux_vars : vector of auxiliary variables
    dir : vector specifying the direction to caculate the flux

    Inputs/Outputs:
    F  : vector to populate with the flux

    This is a low level function.  The static parameters of
    the ParameterType are used to dispatch to the right method for any
    combination of variable type or equation dimension.
"""->
function calcEulerFlux{Tmsh, Tsol, Tres}(params::ParamType{2, :entropy},
                       q::AbstractArray{Tsol,1},
                       aux_vars::AbstractArray{Tres, 1},
                       dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  # calculate some intermediate quantities
  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  s = gamma - q[1] + k1    # entropy
    # internal energy (rho*i in Hughes) - not specific internal energy e
  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)
  U = q[2]*dir[1] + q[3]*dir[2]
  fac = rho_int/q[4]

  # now we can actually calculate the flux
  F[1] = q[4]*U*fac
  F[2] = (dir[1]*gamma_1*q[4] - q[2]*U)*fac
  F[3] = (dir[2]*gamma_1*q[4] - q[3]*U)*fac
  F[4] = U*(k1 - gamma)*fac

  return nothing
end

@doc """
### EulerEquationMod.calcEulerFlux
  This is the 3D method.  All arguments are same as the 2D version.
"""->
# low level function
function calcEulerFlux{Tmsh, Tsol, Tres}(params::ParamType{3, :conservative},
                       q::AbstractArray{Tsol,1},
                       aux_vars::AbstractArray{Tres, 1},
                       dir::AbstractArray{Tmsh},
                       F::AbstractArray{Tsol,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 5), of the conservative variables at the point
# dir is a vector of length 3 that specifies the direction
# F is populated with the flux (is a vector of length 5)
# 3D  only

# once the Julia developers fix slice notation and speed up subarrays, we can make a faster
# vectorized version of this

  press = calcPressure(params, q)
  U = (q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = q[4]*U + dir[3]*press
  F[5] = (q[5] + press)*U

  return nothing

end

function calcEulerFlux_revm{Tmsh,Tsol,Tres}(q::AbstractArray{Tsol,1},
                            dir::AbstractArray{Tmsh,1},
                            dir_bar::AbstractArray{Tmsh,1},
  F_bar::AbstractArray{Tres,1})
  press = gami*(q[4] - 0.5*(q[2]^2 + q[3]^2 + q[4]^2)/q[1])
  U_bar = zero(Tres)
  # F[5] = (q[5] + press)*U
  U_bar += (q[5] + press)*F_bar[5]
  #F[4] = q[4]*U + dir[3]*press
  U_bar += q[4]*F_bar[4]
  dir_bar[3] += press*F_bar[4]
  # F[3] = q[3]*U + dir[2]*press
  U_bar += q[3]*F_bar[3]
  dir_bar[2] += press*F_bar[3]
  # F[2] = q[2]*U + dir[1]*press
  U_bar += q[2]*F_bar[2]
  dir_bar[1] += press*F_bar[2]
  # F[1] = q[1]*U
  U_bar += q[1]*F_bar[1]
  # U = (q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])/q[1]
  dir_bar[1] += q[2]*U_bar/q[1]
  dir_bar[2] += q[3]*U_bar/q[1]
  dir_bar[3] += q[4]*U_bar/q[1]
end

function calcEulerFlux_revm{Tmsh, Tsol}(params::ParamType{3, :conservative},
                            q::AbstractArray{Tsol,1}, aux_vars,
                            dir::AbstractArray{Tmsh,1}, F_bar, dir_bar)

  # Forward sweep
  press = calcPressure(params, q)
  U = (q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = q[4]*U + dir[3]*press
  F[5] = (q[5] + press)*U

  # Reverse sweep


  return nothing
end


function calcEulerFlux{Tmsh, Tsol, Tres}(params::ParamType{3, :entropy},
                       q::AbstractArray{Tsol,1},
                       aux_vars::AbstractArray{Tres, 1},
                       dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  # calculate some intermediate quantities
  k1 = 0.5*(q[2]^2 + q[3]^2 + q[4]^2 )/q[5]  # a constant from Hughes' paper
  s = gamma - q[1] + k1    # entropy
    # internal energy (rho*i in Hughes) - not specific internal energy e
  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[5])^gamma))^(1/gamma_1)
  U = q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3]
  fac = rho_int/q[5]

  # now we can actually calculate the flux
  F[1] = q[5]*U*fac
  F[2] = (dir[1]*gamma_1*q[5] - q[2]*U)*fac
  F[3] = (dir[2]*gamma_1*q[5] - q[3]*U)*fac
  F[4] = (dir[3]*gamma_1*q[5] - q[4]*U)*fac
  F[5] = U*(k1 - gamma)*fac

  return nothing
end



@doc """
### EulerEquationMod.getAuxVars

  This function calculates any extra variables that are stored across the mesh
  using the conservative variables eqn.q.  Currently only calculates pressure.

  Thi is a mid level function
"""->

#------------------------------------------------------------------------------
# functions for calculating additional quantities: pressure, entropy etc.
#------------------------------------------------------------------------------

# mid level function
function getAuxVars{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                      eqn::EulerData{Tsol, Tres, Tdim})
# calculate all auxiliary variables

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_vals = ro_sview(eqn.q, :, j, i)

      # calculate pressure
      press = calcPressure(eqn.params, q_vals)
      @setPressure(eqn.aux_vars, j, i, press)
    end
  end

  return nothing
end


@doc """
### EulerEquationMod.calcPressure

  This function calculates the pressure from the conservative variables at a
  node in 2D.  It returns a single value.

  Inputs:
    params : ParamType{Tdim, var_type }
    q  : vector of conservative variables


  The parameter of params determines whether the 2D or 3D, conservative
  or entropy method is dispatched.

  This is a low level function.

  Aliasing restrictions: none
"""->
# low level function
function calcPressure{Tsol}(params::ParamType{2, :conservative},
                            q::AbstractArray{Tsol,1} )
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  return  (params.gamma_1)*(q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])

end


@doc """
### EulerEquationMod.calcPressure
  This function calculates pressure using the entropy variables.

  Inputs:
    params : ParamType{2, :entropy}
    q  : vector of entropy varaibles

  returns pressure

  Aliasing restrictions: none
"""->
function calcPressure{Tsol}(params::ParamType{2, :entropy},
                            q::AbstractArray{Tsol,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]*q[2] + q[3]*q[3])/q[4]  # a constant from Hughes' paper
  s = gamma - q[1] + k1    # entropy
  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)
  return gamma_1*rho_int
end

@doc """
###EulerEquationMod.calcPressure_revq

Compute the gradient of pressure w.r.t q in the reverse mode

**Arguments**

* `params` : Parameter object
* `q` : Forward sweep solution variable
* `press_bar` : Reverse pressure gradient
* `q_bar` : Reverse mode solution gradient
`
"""->

function calcPressure_revq{Tsol}(params::ParamType{2, :conservative},
                           q::AbstractArray{Tsol,1}, press_bar, 
                           q_bar)

  gamma_1 = params.gamma_1
  q1_inv = 1.0/q[1]
  q_bar[4] += press_bar*gamma_1
  q_bar[3] -= gamma_1*press_bar*q[3]*q1_inv
  q_bar[2] -= gamma_1*press_bar*q[2]*q1_inv
  q_bar[1] += 0.5*gamma_1*press_bar*(q[2]*q[2] + q[3]*q[3])*q1_inv*q1_inv

  return nothing
end


@doc """
### EulerEquationMod.calcPressure

  3D method.  See 2D method documentation
"""->
# low level function
function calcPressure{Tsol}(params::ParamType{3, :conservative},
                            q::AbstractArray{Tsol,1} )
  # calculate pressure for a node
  # q is a vector of length 5 of the conservative variables
  return  (params.gamma_1)*(q[5] - 0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/q[1])

end

function calcPressure{Tsol}(params::ParamType{3, :entropy},
                            q::AbstractArray{Tsol,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/q[5]  # a constant from Hughes' paper
  s = gamma - q[1] + k1    # entropy
  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[5])^gamma))^(1/gamma_1)
  return gamma_1*rho_int
end



@doc """
### EulerEquationMod.calcSpeedofSound

  This function calculates the speed of sound at a node and returns it.
  Methods are available for both conservative and entropy variables.

  Inputs:
    params:  ParamType{Tdim, var_type}
    q  vector of solution variables at a node

  Returns: speed of sound

  This is a low level function

  Aliasing restrictions: none
"""->

function calcSpeedofSound{Tdim, Tsol}(params::ParamType{Tdim, :conservative},
                                q::AbstractArray{Tsol, 1})
# calculates teh speed of sond at a node
  pressure = calcPressure(params, q)
  return sqrt((params.gamma*pressure)/q[1])

end



function calcSpeedofSound{Tsol}(params::ParamType{2, :entropy},
                                q::AbstractArray{Tsol, 1})
# calculate speed of sound using the same formula as conservative variables,
# just rewriting all variables in entropy variables

#  printbacktrace()
#  println("q = ", q)

  gamma = params.gamma
  gamma_1 = params.gamma_1
  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  pressure = calcPressure(params, q)
  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)
  rho = -q[4]*rho_int

#  println("rho = ", rho)
  return sqrt((params.gamma*pressure)/rho)
end

function calcSpeedofSound{Tsol}(params::ParamType{3, :entropy},
                                q::AbstractArray{Tsol, 1})
# calculate speed of sound using the same formula as conservative variables,
# just rewriting all variables in entropy variables

#  printbacktrace()
#  println("q = ", q)

  gamma = params.gamma
  gamma_1 = params.gamma_1
  k1 = 0.5*(q[2]^2 + q[3]^2 + q[4]^2)/q[5]  # a constant from Hughes' paper
  pressure = calcPressure(params, q)
  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[5])^gamma))^(1/gamma_1)
  rho = -q[5]*rho_int

#  println("rho = ", rho)
  return sqrt((params.gamma*pressure)/rho)
end



@doc """
### EulerEquationMod.calcEntropy

  This function calculates the entropy at a node and returns it.  Method are
  available for conservative and entropy variables, 2D or 3D

  Inputs:
    params: ParamType{Tdim, var_type}, used to dispatch to the right method.
    q: vector of solution variables at a node.

  Returns: entropy

  This is a low level function

  Aliasing Restrictions: none

"""->
function calcEntropy{Tsol}(params::ParamType{2, :conservative},
                           q::AbstractArray{Tsol,1} )

  gamma = params.gamma
  gamma_1 = params.gamma_1

  rho_int = q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1]
  return log(gamma_1*rho_int/(q[1]^gamma))
end

function calcEntropy{Tsol}(params::ParamType{2, :entropy},
                           q::AbstractArray{Tsol,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  return gamma - q[1] + 0.5*(q[2]*q[2] + q[3]*q[3])/q[4]
end

function calcEntropy{Tsol}(params::ParamType{3, :conservative},
                           q::AbstractArray{Tsol,1} )

  gamma = params.gamma
  gamma_1 = params.gamma_1

  rho_int = q[5] - 0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/q[1]
  return log(gamma_1*rho_int/(q[1]^gamma))
end

function calcEntropy{Tsol}(params::ParamType{3, :entropy},
                           q::AbstractArray{Tsol,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  return gamma - q[1] + 0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/q[5]
end

"""
  This function calculates the entropy function U used to define the IR
  variablesat a node and returns it.   It does not return the physical entropy
  s.  This function is agnostic to the dimension of the equation.

  Inputs:
    params: a ParamType
    q: a vector of conservative variables at a node

  Outputs:
    U: the value of the entropy function

  Aliasing restrictions: none.
"""
function calcEntropyIR{Tdim, Tsol}(params::ParamType{Tdim, :conservative},
                           q::AbstractArray{Tsol,1} )
# this calculate the entropy functiion U associated with the IR variables,
# not the physical entropy s

  gamma = params.gamma
  gamma_1 = params.gamma_1

  p = calcPressure(params, q)
  rho = q[1]
  U = -rho*(log(p) - gamma*log(rho))/gamma_1
  return U
end

"""
  This function calculates the vorticity at every node of an element
  (because vorticity requires derivatives of the velocity, it would be
  awkward to compute it at a node level).

  3D, conservative variables only

  Inputs:
    params: a ParamType
    q: numDofPerNode x numNodesPerElement array of conservative variables at 
       each node
    dxidx: the 3 x 3 x numNodesPerElement  scaled  mapping jacobian at the node
    jac: numNodesPerElement vector of the mapping jacobian determinant at
         each node

  Input/Outputs:
    vorticity: a 3 x numNodesPerElement array  containing the vorticity in the 
               x, y, and z directions at each node, overwritten

  Aliasing restrictions: from params: dxidx_element, velocities, velocity_deriv,
                                      velocity_deriv_xy
"""
function calcVorticity{Tsol, Tmsh, Tres}(params::ParamType{3, :conservative}, sbp,
                       q::AbstractMatrix{Tsol}, dxidx::Abstract3DArray{Tmsh}, 
                       jac::AbstractVector{Tmsh}, 
                       vorticity::AbstractMatrix{Tres})

  Tdim = 3
  numDofPerNode = size(q, 1)
  numNodesPerElement = size(q, 2)
  dxidx_unscaled = params.dxidx_element
  velocities = params.velocities
  # velocity derivatives in parametric space
  # first dimension: velocity component, 3rd dimension parametric direction
  velocity_deriv = params.velocity_deriv
  fill!(velocity_deriv, 0.0)
  # velocity derivatives in xy space
  # first dimension: velocity, second dimension derivative direction
  velocity_deriv_xy = params.velocity_deriv_xy
  fill!(velocity_deriv_xy, 0.0)
  fill!(vorticity, 0.0)

  # unscale the mapping jacobian
  for i=1:numNodesPerElement
    jac_i = jac[i]
    for j=1:Tdim
      for k=1:Tdim
        dxidx_unscaled[k, j, i] = dxidx[k, j, i]*jac_i
      end
    end
  end

  # compute velocities
  for i=1:numNodesPerElement
    rho_i = q[1, i]
    for j=1:Tdim
      velocities[j, i] = q[j+1, i]/rho_i
    end
  end

  # differentiate velocities
  for d=1:Tdim
    differentiateElement!(sbp, d, velocities, sview(velocity_deriv, :, :, d))
  end

  for i=1:numNodesPerElement
    for v=1:3 # velocity component
      for cart_dim=1:3  # cartesian direction
        for para_dim=1:3  # parametric direction, summed
          velocity_deriv_xy[v, cart_dim, i] += 
              velocity_deriv[v, i, para_dim]*dxidx_unscaled[para_dim, cart_dim]
        end
      end
    end
  end

  # finally, compute the vorticity
  for i=1:numNodesPerElement
    vorticity[1, i] = velocity_deriv_xy[3, 2, i] - velocity_deriv_xy[2, 3, i]
    vorticity[2, i] = -velocity_deriv_xy[3, 1, i] + velocity_deriv_xy[1, 3, i]
    vorticity[3, i] = velocity_deriv_xy[2, 1, i] - velocity_deriv_xy[1, 2, i]
  end

  return nothing
end

#------------------------------------------------------------------------------
# function to calculate various coefficient matrices
#------------------------------------------------------------------------------

@doc """
### EulerEquationMod.calcA0

  This function calculates the A0 (ie. dq/dv, where q are the conservative
  and v are the entropy variables) for a node, and stores it A0

  The formation of A0 is given in Hughes

  Inputs:
    params : ParamType{Tdim, :entropy}
    q  : vector of entropy variables, length Tdim + 2


  Inputs/Outputs:
  A0 : (Tdim+2)x(Tdim+2) matrix populated with A0.  Overwritten

  Aliasing restrictions: none
"""->
function calcA0{Tsol}(params::ParamType{2, :entropy}, q::AbstractArray{Tsol,1},
                      A0::AbstractArray{Tsol, 2})


  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  k2 = k1 - gamma
  k3 = k1*k1 - 2*gamma*k1 + gamma
#    k4 = k2 - gamma_1
  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)

  fac = rho_int/(gamma_1*q[4])

  # calculate the variables used in Hughes A.1
  c1 = gamma_1*q[4] - q[2]*q[2]
  c2 = gamma_1*q[4] - q[3]*q[3]

  d1 = -q[2]*q[3]

  e1 = q[2]*q[4]
  e2 = q[3]*q[4]

  # populate the matrix
  # the matrix is symmetric, but we don't use it because I think populating
  # the matrix will be faster if the matrix is write-only
  A0[1,1] = -q[4]*q[4]*fac
  A0[2,1] = e1*fac
  A0[3,1] = e2*fac
  A0[4,1] = q[4]*(1-k1)*fac
  A0[1,2] = e1*fac  # symmetric
  A0[2,2] = c1*fac
  A0[3,2] = d1*fac
  A0[4,2] = q[2]*k2*fac
  A0[1,3] = e2*fac  # symmetric
  A0[2,3] = d1*fac  # symmetric
  A0[3,3] = c2*fac
  A0[4,3] = q[3]*k2*fac
  A0[1,4] = q[4]*(1-k1)*fac  # symmetric
  A0[2,4] = q[2]*k2*fac   # symmetric
  A0[3,4] = q[3]*k2*fac  # symmetric
  A0[4,4] = -k3*fac

    return nothing
end

function calcA0{Tsol}(params::ParamType{3, :entropy}, q::AbstractArray{Tsol,1},
                      A0::AbstractArray{Tsol, 2})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2 + q[4]^2)/q[5]  # a constant from Hughes' paper
  k2 = k1 - gamma
  k3 = k1*k1 - 2*gamma*k1 + gamma
  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[5])^gamma))^(1/gamma_1)
  fac = rho_int/(gamma_1*q[5])

  # calculate the variables used in Hughes A.1
  c1 = gamma_1*q[5] - q[2]*q[2]
  c2 = gamma_1*q[5] - q[3]*q[3]
  c3 = gamma_1*q[5] - q[4]*q[4]
  d1 = -q[2]*q[3]
  d2 = -q[2]*q[4]
  d3 = -q[3]*q[4]

  e1 = q[2]*q[5]
  e2 = q[3]*q[5]
  e3 = q[4]*q[5]

  # populate the matrix
  # the matrix is symmetric, but we don't use it because I think populating
  # the matrix will be faster if the matrix is write-only
  A0[1,1] = -q[5]*q[5]*fac
  A0[2,1] = e1*fac
  A0[3,1] = e2*fac
  A0[4,1] = e3*fac
  A0[5,1] = q[5]*(1-k1)*fac
  A0[1,2] = e1*fac  # symmetric
  A0[2,2] = c1*fac
  A0[3,2] = d1*fac
  A0[4,2] = d2*fac
  A0[5,2] = q[2]*k2*fac
  A0[1,3] = e2*fac  # symmetric
  A0[2,3] = d1*fac  # symmetric
  A0[3,3] = c2*fac
  A0[4,3] = d3*fac
  A0[5,3] = q[3]*k2*fac
  A0[1,4] = e3*fac  # symmetric
  A0[2,4] = d2*fac  # symmetric
  A0[3,4] = d3*fac  # symmetric
  A0[4,4] = c3*fac
  A0[5,4] = q[4]*k2*fac
  A0[1,5] = q[5]*(1-k1)*fac  # symmetric
  A0[2,5] = q[2]*k2*fac  # symmetric
  A0[3,5] = q[3]*k2*fac  # symmetric
  A0[4,5] = q[4]*k2*fac  # symmetric
  A0[5,5] = -k3*fac

  return nothing
end


@doc """
# EulerEquationMod.calcA0Inv

  Calculates inv(A0), where A0 = dq/dv, where q are the conservative variables
  at a node and v are the entropy varaibles at a node, using the entropy
  variables.

  Inputs:
    params : ParamType{Tdim, :entropy}
    q  : vector (length 4 or 5) of entropy variables at a node

  Inputs/Outputs:
    A0inv : matrix to be populated with inv(A0).  Overwritten.

  Aliasing restrictions: none
"""->
function calcA0Inv{Tsol}(params::ParamType{2, :entropy},
                   q::AbstractArray{Tsol,1},
                   A0inv::AbstractArray{Tsol, 2})
  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
   s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)

  fac = -1/(rho_int*q[4])

  d1 = -q[2]*q[3]
  e1 = q[2]*q[4]
  e2 = q[3]*q[4]


  A0inv[1,1] = (k1*k1 + gamma)*fac
  A0inv[2,1] = k1*q[2]*fac
  A0inv[3,1] = k1*q[3]*fac
  A0inv[4,1] = (k1 + 1)*q[4]*fac
  A0inv[1,2] = k1*q[2]*fac  # symmetry
  A0inv[2,2] = (q[2]*q[2] - q[4])*fac
  A0inv[3,2] = -d1*fac
  A0inv[4,2] = e1*fac
  A0inv[1,3] = k1*q[3]*fac  # symmetry
  A0inv[2,3] = -d1*fac  # symmetry
  A0inv[3,3] = (q[3]*q[3] - q[4])*fac
  A0inv[4,3] = e2*fac
  A0inv[1,4] = (k1 + 1)*q[4]*fac  # symmetry
  A0inv[2,4] = e1*fac  # symmetry
  A0inv[3,4] = e2*fac  # symmetry
  A0inv[4,4] = q[4]*q[4]*fac

    return nothing
end

function calcA0Inv{Tsol}(params::ParamType{3, :entropy},
                   q::AbstractArray{Tsol,1},
                   A0inv::AbstractArray{Tsol, 2})
  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2 + q[4]^2)/q[5]  # a constant from Hughes' paper
   s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[5])^gamma))^(1/gamma_1)

  fac = -1/(rho_int*q[5])

  d1 = -q[2]*q[3]
  d2 = -q[2]*q[4]
  d3 = -q[3]*q[4]

  e1 = q[2]*q[5]
  e2 = q[3]*q[5]
  e3 = q[4]*q[5]


  A0inv[1,1] = (k1*k1 + gamma)*fac
  A0inv[2,1] = k1*q[2]*fac
  A0inv[3,1] = k1*q[3]*fac
  A0inv[4,1] = k1*q[4]*fac
  A0inv[5,1] = (k1 + 1)*q[5]*fac
  A0inv[1,2] = k1*q[2]*fac  # symmetry
  A0inv[2,2] = (q[2]*q[2] - q[5])*fac
  A0inv[3,2] = -d1*fac
  A0inv[4,2] = -d2*fac
  A0inv[5,2] = e1*fac
  A0inv[1,3] = k1*q[3]*fac  # symmetry
  A0inv[2,3] = -d1*fac  # symmetry
  A0inv[3,3] = (q[3]*q[3] - q[5])*fac
  A0inv[4,3] = -d3*fac
  A0inv[5,3] = e2*fac
  A0inv[1,4] = k1*q[4]*fac  # symmetry
  A0inv[2,4] = -d2*fac  # symmetry
  A0inv[3,4] = -d3*fac  # symmetry
  A0inv[4,4] = (q[4]*q[4] - q[5])*fac
  A0inv[5,4] = e3*fac
  A0inv[1,5] = (k1 + 1)*q[5]*fac  # symmetry
  A0inv[2,5] = e1*fac  # symmetry
  A0inv[3,5] = e2*fac  # symmetry
  A0inv[4,5] = e3*fac  # symmetry
  A0inv[5,5] = q[5]*q[5]*fac

    return nothing
end



@doc """
### EulerEquationMod.calcA0

  This function calculates the A0 (ie. dq/dq, where q are the conservative
  variables at a node), and stores it A0.  This function is provided for
  compatability purposes


  Inputs:
    params : ParamType{2, :entropy}
    q  : vector of conservative variables, length 4


  Inputs/Outputs:
  A0 : 4x4 (or 5x5 in 3D)  matrix populated with A0.  Overwritten

  Aliasing restrictions: none
"""->
function calcA0{Tdim, Tsol}(params::ParamType{Tdim, :conservative}, q::AbstractArray{Tsol,1},
                      A0::AbstractArray{Tsol, 2})


  for i=1:length(q)
    A0[i,i] = one(Tsol)
  end

  return nothing
end



@doc """
# EulerEquationMod.calcA0Inv

  Calculates inv(A0), where A0 = dq/dq, where q are the conservative variables
  at a node.  This function is provided for compatability purposes

  Inputs:
    params : ParamType{2, :entropy}
    q  : vector (length Tdim + 2) of conservative variables at a node

  Inputs/Outputs:
    A0inv : matrix to be populated with inv(A0).  Overwritten.

  Aliasing restrictions: none
"""->
function calcA0Inv{Tdim, Tsol}(params::ParamType{Tdim, :conservative},
                   q::AbstractArray{Tsol,1},
                   A0inv::AbstractArray{Tsol, 2})

  for i=1:length(q)
    A0inv[i,i] = one(Tsol)
  end

  return nothing
end


@doc """
### EulerEquationMod.matVecA0inv

  This function takes a 3D array and multiplies it in place by the inv(A0)
  matrix (calculated at each node), inplace, (where A0 =dq/dv, where q are the
  conservative variables and v are some other variables), inplace.
  Methods are available for conservative and entropy variables.

  For conservative variables, A0 is the identity matrix and this is a no-op

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs:
    res_arr: the array to multiply

  Aliasing restrictions: none
"""->
function matVecA0inv{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP,
                     eqn::EulerData{Tsol, Tres, Tdim, :entropy}, opts,
                     res_arr::AbstractArray{Tsol, 3})
# multiply a 3D array by inv(A0) in-place, useful for explicit time stepping
# res_arr *can* alias eqn.q safely
  A0inv = Array(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  res_vals = Array(Tsol, mesh.numDofPerNode)
  q_vals = Array(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      # copy values into workvec
      for k=1:mesh.numDofPerNode
	q_vals[k] = eqn.q[k, j, i]
	res_vals[k] = res_arr[k, j, i]
      end

      res_view = sview(res_arr, :, j, i)
      # get A0Inv for this node
      calcA0Inv(eqn.params, q_vals, A0inv)

      smallmatvec!(A0inv, res_vals, res_view)
    end
  end

  return nothing
end

# no-op, because for conservative variables this is A0inv is the identity matrix
function matVecA0inv{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP,
                     eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                     opts, res_arr::AbstractArray{Tsol, 3})

  return nothing
end

@doc """
### EulerEquationMod.matVecA0

  This function is the same as matVecA0inv, except it multiplies by A0 not
  A0inv.  See its documention.

"""->
function matVecA0{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                  sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim, :entropy},
                  opts, res_arr::AbstractArray{Tsol, 3})
# multiply a 3D array by inv(A0) in-place, useful for explicit time stepping
# res_arr *can* alias eqn.q safely
# a non-alias tolerant implimention wold avoid copying q_vals
  A0 = Array(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  res_vals = Array(Tsol, mesh.numDofPerNode)
  q_vals = Array(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      # copy values into workvec
      for k=1:mesh.numDofPerNode
	q_vals[k] = eqn.q[k, j, i]
	res_vals[k] = res_arr[k, j, i]
      end

      res_view = sview(res_arr, :, j, i)
      # get A0Inv for this node
      calcA0(eqn.params, q_vals, A0)

      smallmatvec!(A0, res_vals, res_view)
    end
  end

  return nothing
end

#no-op
function matVecA0{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim, :conservative}, opts, res_arr::AbstractArray{Tsol, 3})

  return nothing
end






@doc """
### EulerEquationMod.calcA1

  This function calculates the A1 (ie. dF1/dq, where F1 is the first column of
  the Euler flux) for a node, aka the flux
  Jacobian of the Euler flux in the x direction.  Methods are available for
  both conservative and entropy variables.

  The formation of A1 is given in Hughes paper

  Inputs:
    params : ParamType{2, :entropy}
    q  : vector of entropy variables, length 4

  Inputs/Outputs:
  A1 : 4x4 matrix to be populated.  Overwritten


"""->

function calcA1{Tsol}(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1}, A1::AbstractArray{Tsol, 2})
  gamma_1 = params.gamma_1
  gamma = params.gamma
  u = q[2]/q[1] # Get velocity in the x-direction
  v = q[3]/q[1] # Get velocity in the x-direction

  intvar = gamma_1*(u*u + v*v)/2
  a1 = intvar*2 - gamma*q[4]/q[1]
  A1[1,1] = 0
  A1[2,1] = intvar - u*u
  A1[3,1] = -u*v
  A1[4,1] = a1*u

  A1[1,2] = 1
  A1[2,2] = (3 - gamma)*u
  A1[3,2] = v
  A1[4,2] = gamma*q[4]/q[1] - intvar - gamma_1*u*u

  A1[1,3] = 0
  A1[2,3] = -gamma_1*v
  A1[3,3] = u
  A1[4,3] = -gamma_1*u*v

  A1[1,4] = 0
  A1[2,4] = gamma_1
  A1[3,4] = 0
  A1[4,4] = gamma*u

  return nothing
end


function calcA1{Tsol}(params::ParamType{2, :entropy}, q::AbstractArray{Tsol,1},
                      A1::AbstractArray{Tsol, 2})


  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  k2 = k1 - gamma
  k3 = k1*k1 - 2*gamma*k1 + gamma
  k4 = k2 - gamma_1
  k5 = k2*k2 - gamma_1*(k1 + k2)

  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)

  fac = rho_int/(gamma_1*q[4]*q[4])

  # calculate the variables used in Hughes A.1
  c1 = gamma_1*q[4] - q[2]*q[2]
  c2 = gamma_1*q[4] - q[3]*q[3]

  d1 = -q[2]*q[3]

  e1 = q[2]*q[4]
#  e2 = q[3]*q[4]

  # populate the matrix
  # the matrix is symmetric, but we don't use it because I think populating
  # the matrix will be faster if the matrix is write-only
  A1[1,1] = e1*q[4]*fac
  A1[2,1] = c1*q[4]*fac
  A1[3,1] = d1*q[4]*fac
  A1[4,1] = k2*e1*fac
  A1[1,2] = c1*q[4]*fac  # symmetric
  A1[2,2] = -(c1 + 2*gamma_1*q[4])*q[2]*fac
  A1[3,2] = -c1*q[3]*fac
  A1[4,2] = (c1*k2 + gamma_1*q[2]*q[2])*fac
  A1[1,3] = d1*q[4]*fac  # symmetric
  A1[2,3] = -c1*q[3]*fac  # symmetric
  A1[3,3] = -c2*q[2]*fac
  A1[4,3] = k4*d1*fac
  A1[1,4] = k2*e1*fac  # symmetric
  A1[2,4] = (c1*k2 + gamma_1*q[2]*q[2])*fac   # symmetric
  A1[3,4] = k4*d1*fac  # symmetric
  A1[4,4] = k5*q[2]*fac

    return nothing
end


@doc """
### EulerEquationMod.calcA2

  This function calculates A2 (ie. dF2/dq, where F2 is the second column of the
  Euler flux, aka the flux jacobian in the y direction.
  Methods are available for both conservative and entropy variables.

  The formation of A2 is given in Hughes

  Inputs:
    params : ParamType{2, :entropy},
    q  : vector of entropy variables, length 4
  Inputs/Outputs:
  A2 : 4x4 matrix to be populated.  Overwritten


"""->
function calcA2{Tsol}(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1}, A2::AbstractArray{Tsol, 2})
  gamma_1 = params.gamma_1
  gamma = params.gamma
  u = q[2]/q[1] # Get velocity in the x-direction
  v = q[3]/q[1] # Get velocity in the x-direction

  intvar = gamma_1*(u*u + v*v)/2
  a1 = intvar*2 - gamma*q[4]/q[1]

  A2[1,1] = 0
  A2[2,1] = -u*v
  A2[3,1] = intvar - v*v
  A2[4,1] = a1*v
  A2[1,2] = 0
  A2[2,2] = v
  A2[3,2] = -gamma_1*u
  A2[4,2] = -gamma_1*u*v
  A2[1,3] = 1
  A2[2,3] = u
  A2[3,3] = (3 - gamma)*v
  A2[4,3] = gamma*q[4]/q[1] - intvar - gamma_1*v*v
  A2[1,4] = 0
  A2[2,4] = 0
  A2[3,4] = gamma_1
  A2[4,4] = gamma*v


  return nothing
end

function calcA2{Tsol}(params::ParamType{2, :entropy}, q::AbstractArray{Tsol,1},
                      A2::AbstractArray{Tsol, 2})


  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  k2 = k1 - gamma
  k3 = k1*k1 - 2*gamma*k1 + gamma
  k4 = k2 - gamma_1
  k5 = k2*k2 - gamma_1*(k1 + k2)

  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)

  fac = rho_int/(gamma_1*q[4]*q[4])

  # calculate the variables used in Hughes A.1
  c1 = gamma_1*q[4] - q[2]*q[2]
  c2 = gamma_1*q[4] - q[3]*q[3]

  d1 = -q[2]*q[3]

#  e1 = q[2]*q[4]
  e2 = q[3]*q[4]

  # populate the matrix
  # the matrix is symmetric, but we don't use it because I think populating
  # the matrix will be faster if the matrix is write-only
  A2[1,1] = e2*q[4]*fac
  A2[2,1] = d1*q[4]*fac
  A2[3,1] = c2*q[4]*fac
  A2[4,1] = k2*e2*fac
  A2[1,2] = d1*q[4]*fac  # symmetric
  A2[2,2] = -c1*q[3]*fac
  A2[3,2] = -c2*q[2]*fac
  A2[4,2] = k4*d1*fac
  A2[1,3] = c2*q[4]*fac  # symmetric
  A2[2,3] = -c2*q[2]*fac  # symmetric
  A2[3,3] = -(c2 + 2*gamma_1*q[4])*q[3]*fac
  A2[4,3] = (c2*k2 + gamma_1*q[3]*q[3])*fac
  A2[1,4] = k2*e2*fac  # symmetric
  A2[2,4] = k4*d1*fac   # symmetric
  A2[3,4] = (c2*k2 + gamma_1*q[3]*q[3])*fac  # symmetric
  A2[4,4] = k5*q[3]*fac

    return nothing
end

@doc """

"""->

function calcSteadyFluxJacobian(params, A, q, nrm)

  u = q[2]/q[1]
  v = q[3]/q[1]

  a1 = params.gamma*q[4]/q[1]
  theta = nrm[1]*u + nrm[2]*v
  phi_2 = 0.5*params.gamma_1*(u*u + v*v)

  A[1,1] = 0.0
  A[2,1] = -u*theta + nrm[1]*phi_2
  A[3,1] = -v*theta + nrm[2]*phi_2
  A[4,1] = theta*(phi_2-a1)

  A[1,2] = nrm[1]
  A[2,2] = theta - (params.gamma - 2)*nrm[1]*u
  A[3,2] = nrm[1]*v - params.gamma_1*nrm[2]*u
  A[4,2] = nrm[1]*a1 - params.gamma_1*u*theta

  A[1,3] = nrm[2]
  A[2,3] = nrm[2]*u - params.gamma_1*nrm[1]*v
  A[3,3] = theta - (params.gamma - 2)*nrm[2]*v
  A[4,3] = nrm[2]*a1 - params.gamma_1*v*theta

  A[1,4] = 0.0
  A[2,4] = params.gamma_1*nrm[1]
  A[3,4] = params.gamma_1*nrm[2]
  A[4,4] = params.gamma*theta

  return nothing
end


@doc """
### EulerEquationMod.calcMaxWaveSpeed

  This function calculates the maximum wave speed (ie. acoustic wave speed)
  present in the domain and returns it.  Methods are available for conservative   and entropy variables.  This function uses eqn.q_vec, not eqn.q to do the
  calculation.

  This is a mid level function
"""->
function calcMaxWaveSpeed{Tsol, Tdim, Tres}(mesh, sbp,
                          eqn::EulerData{Tsol, Tres, Tdim, :conservative}, opts)
# calculate the maximum wave speed (ie. characteristic speed) on the mesh
# uses solution vector q, not array
  q = eqn.q_vec
  max_speed = zero(Float64)
  for i=1:mesh.numDofPerNode:length(q)
    q_i = ro_sview(q, i:(i+mesh.numDofPerNode - 1))
    a = calcSpeedofSound(eqn.params, q_i)
    u_nrm = zero(Tsol)
    for j=1:Tdim
      u_j = q_i[j+1]/q_i[1]
      u_nrm += u_j*u_j
    end
    u_norm = sqrt(real(u_nrm))
    speed = real(a) + u_norm

    if speed > max_speed
      max_speed = speed
    end  # end if statement
  end  # end loop over vector q

  max_speed = MPI.Allreduce(max_speed, MPI.MAX, mesh.comm)
  return max_speed
end  # end function


function calcMaxWaveSpeed{Tsol, Tres}(mesh, sbp,
                          eqn::EulerData{Tsol, Tres, 2, :entropy}, opts)
# calculate the maximum wave speed (ie. characteristic speed) on the mesh
# uses solution vector q, not array
  q = eqn.q_vec
  gamma = eqn.params.gamma
  gamma_1 = eqn.params.gamma_1

  max_speed = zero(eltype(q))
  for i=1:mesh.numDofPerNode:length(q)
    q_i = ro_sview(q, i:(i+mesh.numDofPerNode - 1))
    a = calcSpeedofSound(eqn.params, q_i)

    # this is dimension specific
    k1 = 0.5*(q_i[2]*q_i[2] + q_i[3]*q_i[3])/q_i[4]  # a constant from Hughes' paper
    s = gamma - q_i[1] + k1    # entropy
    rho_int = exp(-s/gamma_1)*(gamma_1/((-q_i[4])^gamma))^(1/gamma_1)
    rho = -rho_int*q_i[4]
    u_norm = zero(Tsol)
    for j=1:2
      u_j = rho_int*q[j+1]/rho
      u_norm += u_j*u_j
    end
    u_norm = sqrt(u_norm)
    speed = a + u_norm

    if speed > max_speed
      max_speed = speed
    end  # end if statement
  end  # end loop over vector q

  return max_speed
end  # end function


#------------------------------------------------------------------------------
# some old experimental functions
#------------------------------------------------------------------------------
#=
# this is a test for an algorithmic differentiation package
function getEulerJac_wrapper{T}(q::AbstractArray{T,1}, F::AbstractArray{T,1})

  dir = [1.0, 0.0]
#  F = zeros(T, 4)
  sbp = TriSBP{Float64}()
  mesh = PumiMesh2{Float64}(".null", "../../mesh_files/quarter_vortex3l.smb", 1, sbp; dofpernode=4)
  eqn = EulerData1{T, T}(mesh, sbp)

  @time getEulerFlux(eqn, q, dir, F)

  return F

end

fluxJac = forwarddiff_jacobian!(getEulerJac_wrapper, Float64, fadtype=:dual; n=4, m=4)
=#

function calcMomentContribution!{Tsbp,Tmsh,Tres
  }(sbpface::AbstractFace{Tsbp}, xsbp::AbstractArray{Tmsh,3},
    dforce::AbstractArray{Tres,3}, xyz_about::AbstractArray{Tmsh,1})
  @assert( sbpface.numnodes == size(xsbp,2) == size(dforce,2) )
  @assert( size(xsbp,3) == size(dforce,3) )
  @assert( size(xsbp,1) == size(dforce,1) == size(xyz_about,1) )
  rvec = zeros(Tmsh,3)
  dM = zeros(Tres,3)
  moment = zeros(Tres, length(xyz_about))
  for f = 1:size(dforce,3)
    for i = 1:sbpface.numnodes
      for di = 1:3
        rvec[di] = xsbp[di,i,f] - xyz_about[di]
      end
      crossProd(rvec, view(dforce,:,i,f), dM)
      for di = 1:3
        moment[di] += dM[di]*sbpface.wface[i]
      end
    end
  end

  return moment
end

function calcMomentContribution!{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                 eqn::AbstractSolutionData{Tsol, Tres},  
                                 bndry_nums::Array{Int, 1}, 
                                 xyz_about::AbstractArray{Tmsh, 1})

  moment = zeros(Tres, mesh.dim)
  for i=1:length(bndry_nums)
    start_idx = mesh.bndry_offsets[ bndry_nums[i] ]
    end_idx = mesh.bndry_offsets[ bndry_nums[i] + 1 ] - 1
    face_range = start_idx:end_idx
    bndry_faces = sview(mesh.bndryfaces, face_range)
    coords = ro_sview(mesh.coords_bndry, :, :, face_range)
    nrm = ro_sview(mesh.nrm_bndry, :, :, face_range)
    
    # compute dforce
    dforce = computeDForce(mesh, eqn, bndry_faces, nrm)

    # compute moment
    moment += calcMomentContribution!(mesh.sbpface, coords, dforce, xyz_about)

  end

  return moment
end

function calcMomentContribution_revm!{Tmsh, Tres}(mesh::AbstractMesh, eqn::AbstractSolutionData, bndry_nums::Array{Int, 1}, xyz_about::AbstractArray{Tmsh, 1}, moment_bar::AbstractArray{Tres, 1})

  nfaces = length(mesh.bndryfaces)
  moment = zeros(Tres, mesh.dim)
  for i=1:length(bndry_nums)
    start_idx = mesh.bndry_offsets[ bndry_nums[i] ]
    end_idx = mesh.bndry_offsets[ bndry_nums[i] + 1 ] - 1
    face_range = start_idx:end_idx
    bndry_faces = ro_sview(mesh.bndryfaces, face_range)
    coords = ro_sview(mesh.coords_bndry, :, :, face_range)
    coords_bar = zeros(coords)

    nrm = sview(mesh.nrm_bndry, :, :, face_range)
    nrm_bar = sview(mesh.nrm_bndry_bar, :, :, face_range)
    
    # compute dforce
    dforce = computeDForce(mesh, eqn, bndry_faces, nrm)
    dforce_bar = zeros(dforce)

    #--------------------------------------------------------------------------
    # start reverse sweep
    calcMomentContribution_rev!(mesh.sbpface, coords, coords_bar, dforce, dforce_bar, xyz_about, moment_bar)

    computeDForce_revm!(mesh, eqn, bndry_faces, nrm_bar, dforce_bar)
  end

  return nothing
end


function computeDForce{Tmsh, Tsol, Tres}(mesh::AbstractMesh, 
                                         eqn::AbstractSolutionData{Tsol, Tres},
                                         bndryfaces::AbstractArray{Boundary, 1},
                                         nrm::Abstract3DArray{Tmsh})

  nfaces = length(bndryfaces)
  dforce = zeros(Tres, mesh.dim, mesh.numNodesPerFace, nfaces)
  for i=1:nfaces
    for j=1:mesh.numNodesPerFace
      q_j = ro_sview(eqn.q_bndry, :, j, i)  # q_bndry must already have been populated
      p = calcPressure(eqn.params, q_j)
      for k=1:mesh.dim
        dforce[k, j, i] = p*nrm[k, j, i]
      end
    end
  end

  return dforce
end


function computeDForce_revm!{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                             eqn::AbstractSolutionData{Tsol, Tres},
                             bndryfaces::AbstractArray{Boundary, 1},
                             nrm_bar::Abstract3DArray,
                             dforce_bar::Abstract3DArray)

   nfaces = length(bndryfaces)
   for i=1:nfaces
    for j=1:mesh.numNodesPerFace
      q_j = ro_sview(eqn.q_bndry, :, j, i)
      p = calcPressure(eqn.params, q_j)
      for k=1:mesh.dim
        nrm_bar[k, j, i] = dforce_bar[k, j, i]*p
      end
    end
  end

  return nothing
end

@doc """
### calcMomentContribution_rev!

This is the reverse differentiated version of calcMomentContribution!.  See docs
of calcMomentContribution! for further details of the primal method.  This
function is differentiated with respect to the primal version's `xsbp` and
`dforce` variables.

**Inputs**

* `sbpface`: an SBP face operator type
* `xsbp`: SBP-face nodes in physical space; [coord, sbp node, face]
* `dforce`: scaled force at the sbpface nodes; [coord, sbp node, face]
* `xyz_about`: point about which the moment is taken
* `moment_bar`: left multiplies d(moment)/d(xsbp) and d(moment)/d(dforce)

**InOuts**

* `xsbp_bar`: result of vector Jacobian product; [coord, sbp node, face]
* `dforce_bar`: result of vector Jacobian product; [coord, sbp node, face]

"""->
function calcMomentContribution_rev!{Tsbp,Tmsh,Tsol,Tres
  }(sbpface::AbstractFace{Tsbp}, xsbp::AbstractArray{Tmsh,3},
    xsbp_bar::AbstractArray{Tmsh,3}, dforce::AbstractArray{Tsol,3},
    dforce_bar::AbstractArray{Tsol,3}, xyz_about::AbstractArray{Tmsh,1},
    moment_bar::AbstractArray{Tres,1})
  @assert( sbpface.numnodes == size(xsbp,2) == size(xsbp_bar,2)
           == size(dforce,2) == size(dforce_bar,2) )
  @assert( size(xsbp,3) == size(xsbp,3) == size(dforce,3) == size(dforce_bar,3) )
  @assert( size(xsbp,1) == size(xsbp_bar,1) == size(dforce,1)
           == size(dforce_bar,1) == size(xyz_about,1) )
  rvec = zeros(Tmsh,3)
  rvec_bar = zeros(Tres,3)
  dM_bar = zeros(Tres,3)
  for f = 1:size(dforce,3)
    for i = 1:sbpface.numnodes
      for di = 1:3
        rvec[di] = xsbp[di,i,f] - xyz_about[di]
        rvec_bar[di] = 0
        dM_bar[di] = 0  # ???
      end
      # start reverse sweep
      for di = 1:3
        # moment[di] += dM[di]*sbpface.wface[i]
        dM_bar[di] += moment_bar[di]*sbpface.wface[i]
      end
      # crossProd(rvec, view(dforce,:,i,f), dM)
      crossProd_rev(rvec, rvec_bar, view(dforce,:,i,f),
                       view(dforce_bar,:,i,f), dM_bar)
      for di = 1:3
        # rvec[di] = xsbp[di,i,f] - xyz_about[di]
        xsbp_bar[di,i,f] += rvec_bar[di]
      end
    end
  end
end

function calcMomentContribution!{Tsbp,Tmsh,Tres}(sbpface::AbstractFace{Tsbp}, 
    xsbp::AbstractArray{Tmsh,3},
    dforce::AbstractArray{Tres,3}, xyz_about::AbstractArray{Tmsh,1})
  @assert( sbpface.numnodes == size(xsbp,2) == size(dforce,2) )
  @assert( size(xsbp,3) == size(dforce,3) )
  @assert( size(xsbp,1) == size(dforce,1) == size(xyz_about,1) )
  rvec = zeros(Tmsh,3)
  dM = zeros(Tres,3)
  moment = zeros(Tres, length(xyz_about))
  for f = 1:size(dforce,3)
    for i = 1:sbpface.numnodes
      for di = 1:3
        rvec[di] = xsbp[di,i,f] - xyz_about[di]
      end
      crossProd(rvec, view(dforce,:,i,f), dM)
      for di = 1:3
        moment[di] += dM[di]*sbpface.wface[i]
      end
    end
  end

  return moment
end


