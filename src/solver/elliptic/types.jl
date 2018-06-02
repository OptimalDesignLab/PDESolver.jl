
type ParamType{Tdim, Tsol, Tres, Tmsh} <: AbstractParamType
  f::IOStream
  t::Float64  # current time value
  order::Int  # accuracy of elements (p=1,2,3...)

  q_vals::Array{Tsol, 1}  # resuable temporary storage for q variables at a node
  qg::Array{Tsol, 1}  # reusable temporary storage for boundary condition
  v_vals::Array{Tsol, 1}  # reusable storage for convert back to entropy vars.

  res_vals1::Array{Tres, 1}  # reusable residual type storage
  res_vals2::Array{Tres, 1}  # reusable residual type storage

  flux_vals1::Array{Tres, 1}  # reusable storage for flux values
  flux_vals2::Array{Tres, 1}  # reusable storage for flux values

  sat_vals::Array{Tres, 1}  # reusable storage for SAT term

  A0::Array{Tsol, 2}  # reusable storage for the A0 matrix
  A0inv::Array{Tsol, 2}  # reusable storage for inv(A0)
  A1::Array{Tsol, 2}  # reusable storage for a flux jacobian
  A2::Array{Tsol, 2}  # reusable storage for a flux jacobian

  A_mats::Array{Tsol, 3}  # reusable storage for flux jacobians

  Rmat1::Array{Tres, 2}  # reusable storage for a matrix of type Tres
  Rmat2::Array{Tres, 2}

  nrm::Array{Tmsh, 1}  # a normal vectora

  edgestab_gamma::Float64  # edge stabilization parameter
  # debugging options
  writeflux::Bool  # write Euler flux
  writeboundary::Bool  # write boundary data
  writeq::Bool # write solution variables
  use_edgestab::Bool  # use edge stabilization
  use_filter::Bool  # use filtering
  use_res_filter::Bool # use residual filtering

  filter_mat::Array{Float64, 2}  # matrix that performs filtering operation
  # includes transformations to/from modal representation

  use_dissipation::Bool  # use artificial dissipation
  dissipation_const::Float64  # constant used for dissipation filter matrix

  tau_type::Int  # type of tau to use for GLS stabilization

  krylov_itr::Int  # Krylov iteration number for iterative solve
  krylov_type::Int # 1 = explicit jacobian, 2 = jac-vec prod
  time::Timings

  function ParamType(mesh, sbp, opts, order::Integer)
    # create values, apply defaults

    t = 0.0
    myrank = mesh.myrank
    f = open("log_$myrank.dat", "w")
    q_vals = Array(Tsol, 4)
    qg = Array(Tsol, 4)
    v_vals = Array(Tsol, 4)

    res_vals1 = Array(Tres, 4)
    res_vals2 = Array(Tres, 4)

    flux_vals1 = Array(Tres, 4)
    flux_vals2 = Array(Tres, 4)

    sat_vals = Array(Tres, 4)

    A0 = zeros(Tsol, 4, 4)
    A0inv = zeros(Tsol, 4, 4)
    A1 = zeros(Tsol, 4, 4)
    A2 = zeros(Tsol, 4, 4)
    A_mats = zeros(Tsol, 4, 4, Tdim)

    Rmat1 = zeros(Tres, 4, 4)
    Rmat2 = zeros(Tres, 4, 4)

    nrm = zeros(Tmsh, 2)

    edgestab_gamma = opts["edgestab_gamma"]

    # debugging options
    writeflux = opts[ "writeflux"]
    writeboundary = opts[ "writeboundary"]
    writeq = opts["writeq"]
    use_edgestab = opts["use_edgestab"]
    if use_edgestab println("edge stabilization enabled") end

    use_filter = opts["use_filter"]
    if use_filter println("solution variables filter enabled") end


    use_res_filter = opts["use_res_filter"]
    if use_res_filter println("residual filter enabled") end

    if use_filter || use_res_filter || opts["use_filter_prec"]
      filter_fname = opts["filter_name"]
      filter_mat = calcFilter(sbp, filter_fname, opts)
    else
      filter_mat = Array(Float64, 0,0)
    end

    use_dissipation = opts["use_dissipation"]
    if use_dissipation println("artificial dissipation enabled") end

    dissipation_const = opts["dissipation_const"]

    tau_type = opts["tau_type"]


    krylov_itr = 0
    krylov_type = 1 # 1 = explicit jacobian, 2 = jac-vec prod

    time = Timings()
    return new(f, t, order, q_vals, qg, v_vals, res_vals1, res_vals2, sat_vals, flux_vals1,
               flux_vals2, A0, A0inv, A1, A2, A_mats, Rmat1, Rmat2, nrm,
               edgestab_gamma, writeflux, writeboundary,
               writeq, use_edgestab, use_filter, use_res_filter, filter_mat,
               use_dissipation, dissipation_const, tau_type,
               krylov_itr, krylov_type,
               time)

  end   # end of ParamType function

end  # end type declaration

type EllipticData_{Tsol, Tres, Tdim, Tmsh} <: EllipticData{Tsol, Tres, Tdim}
  params::ParamType{Tdim, Tsol, Tres, Tmsh}
  comm::MPI.Comm
  commsize::Int
  myrank::Int
  #
  # Debug Use Only
  #
  xy_face::Array{Tsol, 4} # (x/y, L/R, numnodes, iface)
  #
  # End Debug
  #
  area_sum::Array{Tmsh, 1}
  w::Array{Tmsh, 2} # H/|J|, (numNodesPerElement, numEl)
  q::Array{Tsol, 3}
  q_irk::Array{Tsol, 4}
  q1::Array{Float64, 3}
  q2::Array{Float64, 3}
  q_face::Array{Tsol, 4} # (numDofs, left/right, numNodesPerFace, numFaces)
  # flux_parametric::Array{Tsol, 2}

  xflux_face::Array{Tres, 3}    # stores (u+ - u-)nx*, (numDofs, numNodes, numFaces)
  yflux_face::Array{Tres, 3}    # stores (u+ - u-)ny*
  flux_face ::Array{Tres, 3}    # stores σ⋅norm

  xflux_bndry::Array{Tres, 3}    # stores (u+ - u-)nx*, (numDofs, numNodes, numFaces)
  yflux_bndry::Array{Tres, 3}    # stores (u+ - u-)ny*
  flux_bndry ::Array{Tres, 3}    # stores σ⋅norm

  res::Array{Tres, 3}            # result of computation (numDofs, numNodesPerElem, numElems)
  res_irk::Array{Tres, 4}            # result of computation (numDofs, numNodesPerElem, numElems)
  res1::Array{Float64, 3}            # result of computation (numDofs, numNodesPerElem, numElems)
  res_vec::Array{Tres, 1}         # result of computation in vector form
  q_vec::Array{Tres,1}            # initial condition in vector form
  q_bndry::Array{Tsol, 3}         # store solution variables interpolated to

  shared_data::Array{SharedFaceData{Tsol}, 1}  # MPI send and receive buffers
  #q_face_send::Array{Array{Tsol, 3}, 1}    # send buffers for sending q values
  # to other processes
  #q_face_recv::Array{Array{Tsol, 3}, 1}    # recieve buffers for q values
  #flux_sharedface::Array{Array{Tres, 3}, 1}    # hold shared face flux
  # derivative
  src_func::SRCType    # functor for source term
  src::Array{Tsol, 3} # (numDofPerNode, numNodesPerElement, numEl)
  #
  # Declaration of variables exclusively used in elliptic problems
  #
  # flux_func::FluxType    # functor for the face flux
  diffusion_func::AbstractDiffn
  functional::AbstractFunctional
  calc_energy::AbstractFunctional
  energy::Array{Tsol, 1}
  # diffusion coefficients
  lambda::Array{Tsol, 5} # (d1, d2, ndof, numnodes, n_elemes)
  lambda_face::Array{Tsol, 6} # (d2, d1, ndof, L/R, numNodesPerFace, n_faces)
  lambda_bndry::Array{Tsol, 5} # (d1, d2, ndof, numNodesPerface, n_bndfaces)
  q_grad::Array{Tsol, 4}    # gradient of q, (dim, numDofPerNode, numNodesPerElem, numElems)
  #
  # TBD: Not sure if we store gradient of q at face cubature points
  #
  q_grad_face::Array{Tsol, 5} #(dim, numDofPerNode, left/right, numNodesPerFace, numFaces)
  q_grad_bndry::Array{Tsol, 4} #(dim, numDofPerNode, numNodesPerFace, numBndries)
  #
  # end declaration of variables exclusively used in elliptic problems
  #
  Minv::Array{Float64, 1}
  M::Array{Float64, 1}
  res_edge::Array{Tres, 4}       # edge based residual used for stabilization
  # (numdof, numNodesPerElement, numEl, numEdgesPerEl)
  flux_func::FluxType
  majorIterationCallback::Function
  # assembleSolution::Function
  disassembleSolution::Function
  multiplyA0inv::Function

  q_face_send::Array{Array{Tsol, 3}, 1}  # send buffers for sending q values to other processes
  q_face_recv::Array{Array{Tsol, 3}, 1}  # recieve buffers for q values

  nstages::UInt8
  istage::UInt8

  function EllipticData_(mesh::AbstractMesh, sbp::AbstractSBP, opts)
    println("\nConstructing EllipticData object")
    println("  Tsol = ", Tsol)
    println("  Tres = ", Tres)
    println("  Tdim = ", Tdim)
    println("  Tmsh = ", Tmsh)

    eqn = new()

    eqn.comm = mesh.comm
    eqn.commsize = mesh.commsize
    eqn.myrank = mesh.myrank

    eqn.params = ParamType{Tdim, Tsol, Tres, Tmsh}(mesh, sbp, opts, mesh.order)

    numfacenodes = mesh.numNodesPerFace
    numfaces = mesh.numInterfaces
    numBndFaces = mesh.numBoundaryFaces
    numvars  = mesh.numDofPerNode
    #
    # volume variables
    #
    # TODO: switch the last two dimensions in all gradient variable 
    # in order to keep consistency with lambda
    #
    eqn.w      = zeros(Tmsh, mesh.numNodesPerElement, mesh.numEl)
    eqn.q      = zeros(Tsol, numvars, sbp.numnodes, mesh.numEl)
    if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "BDF2"
      eqn.q1      = zeros(Float64, numvars, sbp.numnodes, mesh.numEl)
      eqn.q2      = zeros(Float64, numvars, sbp.numnodes, mesh.numEl)
    end
    if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "CN"
      eqn.q1      = zeros(Float64, numvars, sbp.numnodes, mesh.numEl)
    end
    eqn.res    = zeros(Tsol, numvars, sbp.numnodes, mesh.numEl)
    if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "SDIRK4"
      # eqn.res_irk = zeros(Tsol, numvars, sbp.numnodes, mesh.numEl, 6)
      eqn.q_irk   = zeros(Tsol, numvars, sbp.numnodes, mesh.numEl, 6)
    end
    eqn.res1    = zeros(Float64, numvars, sbp.numnodes, mesh.numEl)
    eqn.src    = zeros(Tsol, numvars, sbp.numnodes, mesh.numEl)
    eqn.q_grad = zeros(Tsol, numvars, sbp.numnodes, mesh.numEl, Tdim)
    eqn.lambda = zeros(Tsol, Tdim, Tdim, numvars, sbp.numnodes, mesh.numEl)
    #
    # interface variables
    #
    eqn.q_face    = zeros(Tsol, numvars, 2, numfacenodes, numfaces)
    eqn.q_grad_face = zeros(Tsol, numvars, 2, numfacenodes, numfaces, Tdim)
    eqn.flux_face = zeros(Tsol, numvars, numfacenodes, numfaces)
    eqn.xflux_face = zeros(Tsol, numvars, numfacenodes, numfaces)
    eqn.yflux_face = zeros(Tsol, numvars, numfacenodes, numfaces)
    eqn.lambda_face = zeros(Tsol, Tdim, Tdim, numvars, 2, numfacenodes, numfaces)
    #
    # boundary variables
    #
    eqn.q_bndry    = zeros(Tsol, numvars, numfacenodes, numBndFaces)
    eqn.q_grad_bndry = zeros(Tsol, numvars, numfacenodes, numBndFaces, Tdim)
    eqn.flux_bndry = zeros(Tsol, numvars, numfacenodes, numBndFaces)
    eqn.xflux_bndry = zeros(Tsol, numvars, numfacenodes, numBndFaces)
    eqn.yflux_bndry = zeros(Tsol, numvars, numfacenodes, numBndFaces)
    eqn.lambda_bndry = zeros(Tsol, Tdim, Tdim, numvars, numfacenodes, numBndFaces)
    if mesh.isDG
      for i=1:mesh.npeers
        if opts["precompute_face_flux"]
          eqn.flux_sharedface[i] = zeros(Tres, mesh.numDofPerNode, numfacenodes,
                                         mesh.peer_face_counts[i])
        end
        # eqn.aux_vars_sharedface[i] = zeros(Tres, mesh.numDofPerNode,
                                        # numfacenodes, mesh.peer_face_counts[i])
      end
      eqn.shared_data = getSharedFaceData(Tsol, mesh, sbp, opts)
    else
      eqn.shared_data = Array(SharedFaceData, 0)
    end
#
    # reshape of q and res
    #
    eqn.q_vec = reshape(eqn.q, mesh.numDof)
    eqn.res_vec = reshape(eqn.res, mesh.numDof)


    # flux_func = FluxDict["SIPG"]
    eqn.majorIterationCallback = majorIterationCallback
    # eqn.assembleSolution = assembleSolution
    # eqn.disassembleSolution = disassembleSolution
    eqn.multiplyA0inv = multiplyA0inv

    eqn.Minv = calcMassMatrixInverse(mesh, sbp, eqn)
    eqn.M = calcMassMatrix(mesh, sbp, eqn)

    if opts["use_edge_res"]
      eqn.res_edge = zeros(Tres, mesh.numDofPerNode, sbp.numnodes, mesh.numEl,
                           mesh.numTypePerElement[2])
    else
      eqn.res_edge = zeros(Tres, 0, 0, 0, 0)
    end

    eqn.xy_face = zeros(Tsol, 2, 2, numfacenodes, numfaces)

    eqn.area_sum = zeros(Tmsh, mesh.numEl)
    calcWetArea(mesh, sbp, eqn)

    for el = 1 : mesh.numEl
      for n = 1 : mesh.numNodesPerElement
        eqn.w[n, el] = sbp.w[n]/mesh.jac[n, el]
      end
    end
    eqn.energy = zeros(Tsol, mesh.numDofPerNode)
    eqn.calc_energy = FunctionalDict["energy"]

    eqn.istage = 1
    eqn.nstages = 1
    if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "SDIRK4"
      eqn.nstages = 5
    end

    # TODO: parallel related variables
    #
    eqn.q_face_send = Array(Array{Tsol, 3}, 1)
    eqn.q_face_recv = Array(Array{Tsol, 3}, 1)
    return eqn
  end # end function
end # end type

function calcWetArea{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                             sbp::AbstractSBP,
                                             eqn::EllipticData{Tsol, Tres, Tdim})
  nfaces = length(mesh.interfaces)
  nrm = zeros(Tmsh, mesh.numNodesPerFace, Tdim)
  area = zeros(Tmsh, mesh.numNodesPerFace)
  face_area = zero(Tmsh)
  sbpface = mesh.sbpface

  #
  # Compute the wet area of each element
  # 
  for f = 1:nfaces
    face = mesh.interfaces[f]
    eL = face.elementL
    eR = face.elementR
    fL = face.faceL
    fR = face.faceR
    #
    # Compute the size of face
    face_area = 0.0
    for n=1:mesh.numNodesPerFace

      dxidx = sview(mesh.dxidx_face, :, :, n, f)
      # norm vector in reference element
      # nrm_xi = sview(sbp.facenormal, :, fL)
      # nrm_xi = sview(mesh.sbpface.normal, :, fL)
      # nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
      # nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
      # area[n] = sqrt(nrm_xy[n,1]*nrm_xy[n,1] + nrm_xy[n,2]*nrm_xy[n,2])

      nrm_xy = ro_sview(mesh.nrm_face, :, n, f)
      area[n] = norm(nrm_xy)
      face_area += sbpface.wface[n]*area[n]
    end
    eqn.area_sum[eL] += face_area
    eqn.area_sum[eR] += face_area
  end	

  for bc = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[bc]
    indx1 = mesh.bndry_offsets[bc+1] - 1

    for f = indx0:indx1
      face = mesh.bndryfaces[f].face
      elem = mesh.bndryfaces[f].element
      #
      # Compute the size of face
      face_area = 0.0
      for n=1:mesh.numNodesPerFace
        # dxidx = sview(mesh.dxidx_bndry, :, :, n, f)
        # # norm vector in reference element
        # nrm_xi = sview(mesh.sbpface.normal, :, face)
        # nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
        # nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
        # area[n] = sqrt(nrm[n,1]*nrm[n,1] + nrm[n,2]*nrm[n,2])
        nrm_xy = ro_sview(mesh.nrm_bndry, :, :, f)
        area[n] = norm(nrm_xy)
        face_area += sbpface.wface[n]*area[n]
      end
      eqn.area_sum[elem] += 2.0*face_area
    end
  end
  return nothing
end



function majorIterationCallback(itr::Integer,
                                mesh::AbstractMesh,
                                sbp::AbstractSBP,
                                eqn::AbstractEllipticData,
                                opts,
                                f::IO)
  if opts["write_vis"] && ((itr % opts["output_freq"]) == 0 || itr == 1)
    vals = real(eqn.q_vec)
    saveSolutionToMesh(mesh, vals)
    fname = string("solution_", itr)
    writeVisFiles(mesh, fname)
  end

  return nothing
end

function multiplyA0inv{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                                               sbp::AbstractSBP,
                                               eqn::EllipticData{Tsol, Tres, Tdim},
                                               opts,
                                               res::AbstractArray{Tsol, 3})
  return nothing
end

