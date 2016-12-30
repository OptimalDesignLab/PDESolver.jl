@doc """
### ArtificialViscosityType

Composite type with an inner constructor for Artificial Viscosity

**Members**

*  `res_type` :
*  `q`        : Solution Array (eqn.q) where Artificial viscosity is applied
*  `aux_vars` : Auxiliary variables
*  `flux_parametric` : flux in xi direction
*  `res`      : Result of computation
*  `res_vec`  : result of computation in vector form
*  `q_vec`    : initial condition in vector form
*  `edgestab_alpha` : alpha needed by edgestabilization
*  `bndryflux` : boundary flux
*  `stabscale` : stabilization scale factor
*  `Minv`      : inverse mass matrix



"""->

type ArtificialViscosityType{Tsol, Tres, Tdim, Tmsh} <: AbstractSolutionData{Tsol}

  res_type::DataType 
  q::Array{Tsol,3}
  aux_vars::Array{Tres, 3}          # storage for auxiliary variables 
  flux_parametric::Array{Tsol,4}    # flux in xi direction
  res::Array{Tres, 3}               # result of computation
  res_vec::Array{Tres, 1}           # result of computation in vector form
  q_vec::Array{Tres,1}              # initial condition in vector form
  edgestab_alpha::Array{Tmsh, 4}    # alpha needed by edgestabilization
  bndryflux::Array{Tsol, 3}         # boundary flux
  stabscale::Array{Tsol, 2}         # stabilization scale factor
  Minv::Array{Float64, 1}           # inverse mass matrix

  function ArtificialViscosityType(mesh::AbstractMesh, sbp::AbstractSBP, opts)
    eqn = new()  # incomplete initilization
    eqn.params = ParamType{Tdim}(opts)
  end

end # End Type
