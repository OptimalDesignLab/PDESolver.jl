using ODLCommonTools

export ArtificialViscosityType

type ArtificialViscosityType{Tsol, Tres, Tdim, Tmsh} <: AbstractSolutionData{Tsol}

  res_type::DataType 
  q::Array{Tsol,3}
  aux_vars::Array{Tres, 3}  # storage for auxiliary variables 
  flux_parametric::Array{Tsol,4}  # flux in xi direction
  res::Array{Tres, 3}  # result of computation
  SL::Array{Tres, 1}  # result of computation in vector form
  SL0::Array{Tres,1}  # initial condition in vector form
  edgestab_alpha::Array{Tmsh, 4} # alpha needed by edgestabilization
  bndryflux::Array{Tsol, 3}  # boundary flux
  stabscale::Array{Tsol, 2}  # stabilization scale factor
  Minv::Array{Float64, 1}  # invese mass matrix

  function ArtificialViscosityType(mesh::PumiMesh2, sbp::SBPOperator, opts)
    eqn = new()  # incomplete initilization
    eqn.params = ParamType{Tdim}(opts)

    
  end

end # End Type
