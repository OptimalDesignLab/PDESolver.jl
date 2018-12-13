abstract type AbstractDirichletBC <: BCType end
abstract type AbstractNeumannBC <: BCType end
export isDirichlet, isNeumann

mutable struct DirichletAllZero <: AbstractDirichletBC
end
function (obj::DirichletAllZero)(
                          xy::AbstractArray{Tmsh},
                          gD::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}  # (numDofPerNode))
  gD[:] = 0.0
  return nothing
end

mutable struct DirichletAllOne <: AbstractDirichletBC
end
function (obj::DirichletAllOne)(
                          xy::AbstractArray{Tmsh},
                          gD::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}  # (numDofPerNode)
  gD[:] = 1.0
  return nothing
end

mutable struct DirichletTrig <: AbstractDirichletBC
end
function (obj::DirichletTrig)(
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}  # (numDofPerNode))
  k = 2.0
  q[:] = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  if abs(q[1]) > 1.0E-10
    error("some problem!")
  end
  return nothing
end

mutable struct DirichletPolynial2nd <: AbstractDirichletBC
end
function (obj::DirichletPolynial2nd)(
                          xy::AbstractArray{Tmsh},
                          gD::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}  # (numDofPerNode)
  a = 1.0
  b = 1.0
  c = 1.0
  gD[:] =  a*xy[1]*xy[1] + b*xy[1]*xy[2] + c*xy[2]*xy[2]
  # gD[:] =  -0.25*xy[1]*xy[1] - 0.25*xy[2]*xy[2]
end


mutable struct NeumannAllZero <: AbstractNeumannBC
end

function (obj::NeumannAllZero)(
                          xy::AbstractArray{Tmsh, 1},
                          gN::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}  # (numDofPerNode)
  gN[:] = 1.0
  return nothing
end

mutable struct NeumannAllOne <: AbstractNeumannBC
end

function (obj::NeumannAllOne)(
                          xy::AbstractArray{Tmsh, 1},
                          gN::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}  # (numDofPerNode)
  gN[:] = 1.0
  return nothing
end

global const BCDict = Dict{String, BCType}(
                                                "DirichletAllZero" => DirichletAllZero(),
                                                "DirichletAllOne" => DirichletAllOne(),
                                                "NeumannAllZero" => NeumannAllZero(),
                                                "NeumannAllOne" => NeumannAllOne(),
                                                "DirichletPolynial2nd" => DirichletPolynial2nd(),
                                                "DirichletTrig" => DirichletTrig()
                                               )

function isDirichlet(dBC::AbstractDirichletBC)
  return true
end
function isDirichlet(dBC::AbstractNeumannBC)
  return false
end

function isNeumann(dBC::AbstractDirichletBC)
  return false
end
function isNeumann(dBC::AbstractNeumannBC)
  return true
end

function getBCFunctors(mesh::AbstractMesh, sbp::AbstractOperator, eqn::EllipticData, opts)
  for i = 1:mesh.numBC
    key = string("BC", i, "_name")
    val = opts[key]
    mesh.bndry_funcs[i] = BCDict[val]
  end
end
function interpolateBoundary(mesh::AbstractDGMesh,
                             sbp::AbstractOperator,
                             eqn::AbstractEllipticData{Tsol, Tres},
                             opts,
                             q::AbstractArray{Tsol,3},
                             q_bndry::AbstractArray{Tsol,3}) where {Tsol, Tres}
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, q, q_bndry)
end

function interpolateBoundary(mesh::AbstractDGMesh,
                             sbp::AbstractOperator,
                             eqn::AbstractEllipticData{Tsol, Tres},
                             opts,
                             grad::AbstractArray{Tsol,4},
                             grad_bndry::AbstractArray{Tsol,4}) where {Tsol, Tres}
  @assert(size(grad, 4) == size(grad_bndry, 4))   # Tdim
  @assert(size(grad, 1) == size(grad_bndry, 1))   # numDofPerNode

  dim = size(grad, 4)
  for d=1:dim
    g = sview(grad, :, :, :, d)
    g_bndry = sview(grad_bndry, : , : , :, d)
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, g, g_bndry)
  end
end

function getBCFluxes(mesh::AbstractMesh{Tmsh},
                            sbp::AbstractOperator,
                            eqn::EllipticData{Tsol, Tres, Tdim},
                            opts) where {Tmsh, Tsol, Tres, Tdim}
  #
  # The boundary integrals are categorized into 3 classes.
  # The first is numerical and applies to all boundaries (both Dirichlet
  # and Neumann). The second class is `physical` Dirichlet and the third
  # one is `physical` Neumann. These 3 classes are dealed with separatedly
  # in 3 loops.
  #

  #
  # Take into account Direclet and Neumann boundary conditions
  #
  p = opts["order"]
  Cip = opts["Cip"]
  penalty_method = opts["Flux_name"]
  sbpface = mesh.sbpface
  dq = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)  
  penalty_factor_shahbazi = 0.5*Float64(p + 1.0)*Float64(p + Tdim)/Float64(Tdim)
  penalty = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  penalty_shahbazi = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  sbpface = mesh.sbpface
  numFacesPerElem = 3
  nrm = Array{Tmsh}(Tdim, mesh.numNodesPerFace)
  nrm1 = Array{Tmsh}(Tdim, mesh.numNodesPerFace)
  area = Array{Tmsh}(mesh.numNodesPerFace)
  numFacesPerElem = 3

  eigMax = Array{Tmsh}(mesh.numDofPerNode)
  lambda_dqdx = Array{Tsol}(Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)

  Sat = Array{Tmsh}(mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace) 
  RLR = Array{Tmsh}(mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace)
  R = sview(sbpface.interp, :,:)
  RR = R.'*R
  area_sum = sview(eqn.area_sum, :)

  perm = zeros(Tmsh, sbp.numnodes, sbpface.stencilsize)
  Hinv = zeros(Tmsh, sbp.numnodes, sbp.numnodes)
  Bsqrt = zeros(Tmsh, sbpface.numnodes, sbpface.numnodes)
  for s = 1:sbpface.stencilsize
    perm[sbpface.perm[s, 1], s] = 1.0
  end
  for i = 1:sbp.numnodes
    Hinv[i,i] = 1.0/sbp.w[i]
  end
  for i = 1:sbpface.numnodes
    Bsqrt[i,i] = sqrt(sbpface.wface[i])
  end
  BsqrtRHinvRtBsqrt = Bsqrt*R.'*perm.'*Hinv*perm*R*Bsqrt
  sigma = eigmax(BsqrtRHinvRtBsqrt)

  relax_coef = 1.0
  if haskey(opts, "unstable_coef")
    relax_coef = opts["unstable_coef"]
  end

  for iBC = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[iBC]
    indx1 = mesh.bndry_offsets[iBC+1] - 1
    gD = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
    gg = Array{Tsol}(mesh.numDofPerNode)
    gN = Array{Tsol}(mesh.numDofPerNode)

    for f = indx0:indx1
      xflux = sview(eqn.xflux_bndry, :,:,f)
      yflux = sview(eqn.yflux_bndry, :,:,f)
      flux  = sview(eqn.flux_bndry, :,:,f)
      bndry = mesh.bndryfaces[f]
      elem = bndry.element

      #
      # Compute geometric info on face
      #
      nrm = ro_sview(mesh.nrm_bndry, :, :, f)
      for n=1:mesh.numNodesPerFace
        area[n] = norm(ro_sview(nrm, :, n)) 
        nrm1[:,n] = nrm[:,n] / area[n]
      end

      dqdx = sview(eqn.q_grad_bndry, :, :, f, 1)
      dqdy = sview(eqn.q_grad_bndry, :, :, f, 2)

      lambda = sview(eqn.lambda_bndry, :,:,:,:,f)

      for n = 1:mesh.numNodesPerFace
        for var = 1:mesh.numDofPerNode
          lambda_dqdx[1, var, n] = lambda[1, 1, var, n]*dqdx[var, n] + lambda[1, 2, var, n]*dqdy[var, n] 
          lambda_dqdx[2, var, n] = lambda[2, 1, var, n]*dqdx[var, n] + lambda[2, 2, var, n]*dqdy[var, n] 
        end
      end

      for n = 1:mesh.numNodesPerFace
        xy = sview(mesh.coords_bndry, :, n, f)
        q = sview(eqn.q_bndry, :, n, f)
        #
        # term-2
        #
        xflux[:, n] = -q[:]*nrm[1, n]
        yflux[:, n] = -q[:]*nrm[2, n]
        #
        # term-3
        #
        flux[:, n] = -lambda_dqdx[1, :, n]*nrm[1, n] - lambda_dqdx[2, :, n]*nrm[2, n]
      end
    end

    #
    # physical Dirichlet
    #
    bc_func = mesh.bndry_funcs[iBC]
    if isDirichlet(bc_func) # Dirichlet BC, term-5 and term-7
      for f = indx0:indx1
        xflux = sview(eqn.xflux_bndry, :,:,f)
        yflux = sview(eqn.yflux_bndry, :,:,f)
        flux  = sview(eqn.flux_bndry, :,:,f)
        bndry = mesh.bndryfaces[f]
        elem = bndry.element
        face = bndry.face
        perm = sview(sbpface.perm, :, face)
        stencilSize = sbpface.stencilsize
        #
        # Compute geometric info on face
        #
        nrm = ro_sview(mesh.nrm_bndry, :, :, f)
        for n = 1 : mesh.numNodesPerFace
          # area[n] = norm(view(nrm, :, n))
          area[n] = norm(ro_sview(nrm, :, n))
          nrm1[1, n] = nrm[1, n]/area[n]
          nrm1[2, n] = nrm[2, n]/area[n]
        end
        #
        # compute element size `he`
        #
        elem_vol = 0.0
        for n = 1:mesh.numNodesPerElement
          elem_vol += sbp.w[n]/mesh.jac[n, elem]
        end

        face_area = 0.0
        for n = 1:mesh.numNodesPerFace
          face_area += sbpface.wface[n]*area[n]
        end
        # 
        # On boundary, the area is weighted twice
        #
        area_weight = 0.5*area_sum[elem]/face_area

        if penalty_method == "Shahbazi" || penalty_method == "SAT0"
          he = elem_vol/eqn.area_sum[elem]
        end

        lambda = sview(eqn.lambda, :,:,:,:,elem)

        if penalty_method == "SAT0"
          eigMax[:] = -1.0

          for dof = 1 : mesh.numDofPerNode
            # left element
            for n = 1 : mesh.numNodesPerElement 
              b = real(lambda[1,1,dof,n]) + real(lambda[2,2,dof,n])
              ac = real(lambda[1,1,dof,n]) * real(lambda[2,2,dof,n]) - real(lambda[1,2,dof,n]) * real(lambda[2,1,dof,n])
              root = 0.5*(b + sqrt(b*b - 4*ac))*mesh.jac[n, elem]
              eigMax[dof] = max(eigMax[dof], root)
            end

            for n = 1:mesh.numNodesPerFace
              # 0.25 = 0.5*0.5. 
              # One 1/2 comes from double weight for boundary face;
              # the other comes from 1/sum(sbpface.wface[i])
              penalty[dof, n] = eigMax[dof]*sigma*area_weight*area[n]*area[n]
            end
          end
        elseif penalty_method == "Shahbazi" 
          lambda_face = sview(eqn.lambda_bndry, :,:, :, :, f)

          for n = 1:mesh.numNodesPerFace
            for dof = 1:mesh.numDofPerNode
              penalty[dof, n] =  nrm1[1, n]*nrm1[1, n]*lambda_face[1, 1, dof, n]
              penalty[dof, n] += nrm1[1, n]*nrm1[2, n]*lambda_face[1, 2, dof, n]
              penalty[dof, n] += nrm1[2, n]*nrm1[1, n]*lambda_face[2, 1, dof, n]
              penalty[dof, n] += nrm1[2, n]*nrm1[2, n]*lambda_face[2, 2, dof, n]
              penalty[dof, n] *= penalty_factor_shahbazi*area[n]/he
            end
          end 
        elseif penalty_method == "SAT"

          lambda = sview(eqn.lambda, :,:,:,:,elem)
          Sat[:,:,:] = 0.0  

          # RLR_L and RLR_R are reused for each (d1, d2)
          for d1 = 1:Tdim
            for d2 = 1:Tdim
              # Compute R Λ{d1,d2}/H R^T
              RLR[:,:,:] = 0.0
              for n1 = 1:mesh.numNodesPerFace
                for n2 = 1:mesh.numNodesPerFace
                  for k = 1:stencilSize
                    RLR[:, n2, n1] += lambda[d1, d2, :, perm[k]] / eqn.w[perm[k], elem] * R[k, n2] * R[k, n1] 
                  end
                end
              end
              #
              # N_{d1}i} RLR N_{d2}
              #
              for n1 = 1:mesh.numNodesPerFace
                for n2 = 1:mesh.numNodesPerFace
                  Sat[:, n2, n1] += nrm1[d1, n2]*RLR[:, n2, n1]*nrm1[d2, n12]
                end
              end
            end
          end

          #
          # Left multiply by numFacesPerElem*B
          #
          for row = 1:mesh.numNodesPerFace
            for col = 1:mesh.numNodesPerFace
              Sat[:, row, col] *= area_weight*sbpface.wface[col]*area[col]
            end
          end

        else
          for n = 1:mesh.numNodesPerFace
            for dof = 1:mesh.numDofPerNode
              penalty[dof, n] = Cip*p*p/he
            end
          end
        end

        q = sview(eqn.q_bndry, :, :, f)
        for n = 1:mesh.numNodesPerFace
          xy = sview(mesh.coords_bndry, :, n, f)
          # bc_func(xy, sview(gD, :, n))
          bc_func(xy, gg)
          gD[:, n] = gg[:]
          for dof = 1:mesh.numDofPerNode
            dq[dof, n] = q[dof, n] - gD[dof, n]
          end
        end

        if penalty_method == "Shahbazi" || penalty_method == "SAT0" 
          for n = 1:mesh.numNodesPerFace
            for dof = 1 : mesh.numDofPerNode
              flux[dof, n] += penalty[dof, n]*dq[dof, n] * relax_coef
            end
          end
        elseif penalty_method == "SAT" 
          for n = 1:mesh.numNodesPerFace
            for n1 = 1:mesh.numNodesPerFace
              for dof = 1:mesh.numDofPerNode
                flux[dof, n] += area[n]*Sat[dof, n, n1]*dq[dof, n1] * relax_coef
              end
            end
          end
        else
          error("We should never get here") 
        end

        #
        # term-7
        #
        for n = 1:mesh.numNodesPerFace
          xflux[:, n] += gD[:,n]*nrm[1, n]
          yflux[:, n] += gD[:,n]*nrm[2, n]
        end
        #
      end
      # physical Neumann
      #
    elseif isNeumann(bc_func)   # Neumann BC, term-8
      for f = indx0:indx1
        flux  = sview(eqn.flux_bndry, :,:,f)
        bndry = mesh.bndryfaces[f]
        #
        # Compute geometric info on face
        #
        nrm = ro_sview(mesh.nrm_bndry, :, :, f)
        for n = 1 : mesh.numNodesPerFace
          # area[n] = norm(view(nrm, :, n))
          area[n] = norm(ro_sview(nrm, :, n))
          nrm1[1, n] = nrm[1, n]/area[n]
          nrm1[2, n] = nrm[2, n]/area[n]
        end

        for j = 1:mesh.numNodesPerFace
          xy = sview(mesh.coords_bndry, :, j, f)
          bc_func(xy, gN)
          # term-8
          flux[:, j] -= gN*area[n]
        end
      end
    else
      error("We should never get here!")
    end
  end

  return nothing 
end
