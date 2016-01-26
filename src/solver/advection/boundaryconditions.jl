# boundaryconditions.jl

@doc """
### AdvectionEquationMod.calcBoundaryFlux

  This function calculates the boundary flux for the portion of the boundary
  with a particular boundary condition.  The eqn.q are converted to 
  conservative variables if needed

  Inputs:
  mesh : AbstractMesh
  sbp : SBPOperator
  eqn : AdvectionEquation
  functor : a callable object that calculates the boundary flux at a node
  bndry_facenums:  An array with elements of type Boundary that tell which
                   element faces have the boundary condition
  Outputs:
  bndryflux : the array to store the boundary flux, corresponds to 
              bndry_facenums

  The functor must have the signature
  functor( q, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
  where q are the *conservative* variables.
  where all arguments (except params and nrm) are vectors of values at a node.

  params is the ParamType associated with the the EulerEquation object
  nrm = sbp.facenormal[:, current_node]

  This is a mid level function.
"""->
# mid level function
function calcBoundaryFlux{Tmsh,  Tsol, Tres}( mesh::AbstractMesh{Tmsh}, 
                          sbp::SBPOperator, eqn::AdvectionData{Tsol}, 
                          functor::BCType, 
                          bndry_facenums::AbstractArray{Boundary,1}, 
                          bndryflux::AbstractArray{Tres, 3})
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("enterted calcBoundaryFlux")
  

  nfaces = length(bndry_facenums)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
#    println("element = ", bndry_i.element, ", face = ", bndry_i.face)
#    println("interface ", i)
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]

      # get components
      q = view(eqn.q, :, k, bndry_i.element)
      alpha_x = view(eqn.alpha_x, 1, k, bndry_i.element)
      alpha_y = view(eqn.alpha_y, 1, k, bndry_i.element)
      # flux_parametric = view(eqn.flux_parametric, :, k, bndry_i.element, :)
      # aux_vars = view(eqn.aux_vars, :, k, bndry_i.element)
      coords = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("eqn.bndryflux = ", eqn.bndryflux)
      bndryflux_i = view(bndryflux, :, j, i)

      # functor(u, flux_parametric, aux_vars, coords, dxidx, nrm, bndryflux_i, eqn.params)
      functor(q, alpha_x, alpha_y, coords, dxidx, nrm, bndryflux_i)
    end
  end

  return nothing
end

type x5plusy5BC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::x5plusy5BC, u::AbstractArray{Tsol,1}, 
              alpha_x::Tsol, alpha_y::Tsol, coords::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1})

  u_bc = 0.0
  calc_x5plusy5(coords, u_bc) # Calculate the actual analytic value of u at the bondary
  RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx, bndryflux)

  return nothing
end

type exp_xplusyBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::exp_xplusyBC, u::AbstractArray{Tsol,1}, 
              alpha_x::Tsol, alpha_y::Tsol, coords::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1})

  u_bc = 0.0
  calc_exp_xplusy(coords, u_bc)
  RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx, bndryflux)

  return nothing
end


@doc """
### AdvectionEquationMod.BCDict

It stores all the possible boundary condition dictionary options. Whenever a 
new boundary condition is created, it should get added to BCDict.

"""->
global const BCDict = Dict{ASCIIString, BCType} (
"x5plusy5BC" => x5plusy5BC(),
"exp_xplusyBC" => exp_xplusyBC()
)


@doc """
### AdvectionEquationMod.getBCFunctors

This function uses the opts dictionary to populate mesh.bndry_funcs with
the the functors

This is a high level function.

**Inputs**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Input dictionary options

**Outputs**

*  None

"""->
# use this function to populate access the needed values in BCDict
function getBCFunctors{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
                       eqn::AdvectionData{Tsol, Tdim}, opts)
# populate the array mesh.bndry_funcs with the functors for the boundary condition types

#  println("Entered getBCFunctors")

  for i=1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
    println("BCDict[val] = ", BCDict[val])
    mesh.bndry_funcs[i] = BCDict[val]
  end

  return nothing
end