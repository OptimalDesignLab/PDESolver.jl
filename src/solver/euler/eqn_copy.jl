# ??: EulerData_ or EulerData?
function copy(eqn::EulerData_, mesh::AbstractMesh, sbp::AbstractSBP, opts::Dict)

  # 1: call constructor on eqn_copy

  eqn_copy = EulerData_{Tsol, Tres, Tdim, Tmsh, var_type}(mesh, sbp, opts)

  # should be zero'ed out

  # 2: copy over fields
  eqn_copy.f = copy(eqn.f)                            # f::IOStream
  eqn_copy.t = copy(eqn.t)                            # t::Float64
  eqn_copy.order = copy(eqn.order)                    # order::Int



end

"""
  TODO Document this
"""
function getStaticParams(opts; return_type="dict")

  if return_type == "dict"
    staticParams = Dict()
    staticParams = opts["staticParams"]

    return staticParams
  elseif return_type == "tuple"

    Tsol = opts["staticParams"]["Tsol"]
    Tres = opts["staticParams"]["Tres"]
    Tmsh = opts["staticParams"]["Tmsh"]
    Tdim = opts["staticParams"]["Tdim"]
    var_type = opts["staticParams"]["var_type"]

    staticParamTuple = (Tsol, Tres, Tmsh, Tdim, var_type)

    return staticParamTuple
  end

end


#=
# list of fields
  f::IOStream
  t::Float64  # current time value
  order::Int  # accuracy of elements (p=1,2,3...)

  #TODO: consider making these vectors views of a matrix, to guarantee
  #      spatial locality
  q_vals::Array{Tsol, 1}  # resuable temporary storage for q variables at a node
  q_vals2::Array{Tsol, 1}
  q_vals3::Array{Tsol, 1}
  qg::Array{Tsol, 1}  # reusable temporary storage for boundary condition
  v_vals::Array{Tsol, 1}  # reusable storage for convert back to entropy vars.
  v_vals2::Array{Tsol, 1}
  Lambda::Array{Tsol, 1}  # diagonal matrix of eigenvalues

  # numDofPerNode x stencilsize arrays for entropy variables
  w_vals_stencil::Array{Tsol, 2}
  w_vals2_stencil::Array{Tsol, 2}

  res_vals1::Array{Tres, 1}  # reusable residual type storage
  res_vals2::Array{Tres, 1}  # reusable residual type storage
  res_vals3::Array{Tres, 1}  

  flux_vals1::Array{Tres, 1}  # reusable storage for flux values
  flux_vals2::Array{Tres, 1}  # reusable storage for flux values

  sat_vals::Array{Tres, 1}  # reusable storage for SAT term

  A0::Array{Tsol, 2}  # reusable storage for the A0 matrix
  A0inv::Array{Tsol, 2}  # reusable storage for inv(A0)
  A1::Array{Tsol, 2}  # reusable storage for a flux jacobian
  A2::Array{Tsol, 2}  # reusable storage for a flux jacobian
  S2::Array{Tsol, 1}  # diagonal matrix of eigenvector scaling

  A_mats::Array{Tsol, 3}  # reusable storage for flux jacobians

  Rmat1::Array{Tres, 2}  # reusable storage for a matrix of type Tres
  Rmat2::Array{Tres, 2}

  P::Array{Tmsh, 2}  # projection matrix

  nrm::Array{Tmsh, 1}  # a normal vector
  nrm2::Array{Tmsh, 1}
  nrm3::Array{Tmsh, 1}

  h::Float64 # temporary: mesh size metric

  cv::Float64  # specific heat constant
  R::Float64  # specific gas constant used in ideal gas law (J/(Kg * K))
  gamma::Float64 # ratio of specific heats
  gamma_1::Float64 # = gamma - 1

  Ma::Float64  # free stream Mach number
  Re::Float64  # free stream Reynolds number
  aoa::Float64  # angle of attack
  rho_free::Float64  # free stream density
  E_free::Float64 # free stream energy (4th conservative variable)

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

  vortex_x0::Float64  # vortex center x coordinate at t=0
  vortex_strength::Float64  # strength of the vortex

  krylov_itr::Int  # Krylov iteration number for iterative solve
  krylov_type::Int # 1 = explicit jacobian, 2 = jac-vec prod

  Rprime::Array{Float64, 2}  # numfaceNodes x numNodesPerElement interpolation matrix
                             # this should live in sbpface instead
  # temporary storage for calcECFaceIntegrals
  A::Array{Tres, 2}
  B::Array{Tres, 3}
  iperm::Array{Int, 1}

  time::Timings

=#
