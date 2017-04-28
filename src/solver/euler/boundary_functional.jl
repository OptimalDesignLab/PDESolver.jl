# Calculate the analytical force on the inner boundary of the isentropic vortex
#=
function calc_analytical_forces{Tmsh}(mesh::AbstractMesh{Tmsh}, params::ParamType{2},
	                                    coords::AbstractArray{Tmsh})

  q = zeros(mesh.numDofPerNode)
  calcIsentropicVortex(coords, params, q)  # Get analytical q ath the coordinates
  p = calcPressure(params, q) # Get the analytical pressure
  r = sqrt(coords[1]*coords[1] + coords[2]*coords[2]) # get the curcular radius
  force = 0.5*pi*r*p

  return force
end
=#
@doc """
### EulerEquationMod.calcBndryFunctional

This function calculates a functional on a geometric boundary of a the 
computational space. There is no need to call this function withing the 
nonlinear solve while computing eqn.q

**Inputs**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `g_edge_number` : Geometric edge number

**Outputs**

*  `functional_val` : computed numerical force at the boundary.

"""->

function calcBndryFunctional{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},
                                         sbp::AbstractSBP,
                                         eqn::EulerData{Tsol}, 
                                         opts, 
                                         functor, 
                                         functional_edges)

  functional_val = zero(Tsol)

  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)
    q2 = zeros(Tsol, mesh.numDofPerNode)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = sview(eqn.q_bndry, :, j, global_facenum)
        convertToConservative(eqn.params, q, q2)
        aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = sview(mesh.coords_bndry, :, j, global_facenum)
        dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

        boundary_integrand[1,j,i] = functor(eqn.params, q, aux_vars, [nx, ny])
      
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces

  val_per_geom_edge = zeros(Tsol, 1)

  integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], 
                         boundary_integrand, val_per_geom_edge)

  functional_val += val_per_geom_edge[1]
  end  # End for itr = 1:length(functional_edges)

  return functional_val
end

@doc """
### EulerEquationMod.drag

Computes the force in the X-direction.

**Inputs**

*  `params` : Parameter type
*  `q`      : Solution at a node
*  `aux_vars` : Vector of auxiliary variables
*  `nrm`    : Normal vector in the physical space

**Outputs**

*  `val`    : Momentum derivative in the X-direction

"""->

type drag <: FunctionalType
end

function call{Tsol, Tres, Tmsh}(obj::drag, params, q::AbstractArray{Tsol,1}, 
              aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh})

  euler_flux = zeros(Tsol, length(q))
  calcEulerFlux(params, q, aux_vars, nrm, euler_flux)
  val = euler_flux[2]

  return val
end

@doc """
### EulerEquationMod.lift

Computes the force in the Y-direction.

**Inputs**

*  `params` : Parameter type
*  `q`      : Solution at a node
*  `aux_vars` : Vector of auxiliary variables
*  `nrm`    : Normal vector in the physical space

**Outputs**

*  `val`    : Momentum derivative in the Y-direction

"""->

type lift <: FunctionalType
end

function call{Tsol, Tres, Tmsh}(obj::lift, params, q::AbstractArray{Tsol,1}, 
              aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh})

  euler_flux = zeros(Tsol, length(q))
  calcEulerFlux(params, q, aux_vars, nrm, euler_flux)
  val = euler_flux[3]

  return val
end


@doc """
### EulerEquationMod.FunctionalDict

It stores the names of all possible functional options that can be computed. 
Whenever a new functional is created, it should be added to FunctionalDict.

"""->
global const FunctionalDict = Dict{ASCIIString, FunctionalType}(
"drag" => drag(),
"lift" => lift(),
)


@doc """
### EulerEquationMod.getFunctionalName

Gets the name of the functional that needs to be computed at a particular point

**Inputs**

*  `opts`     : Input dictionary
*  `f_number` : Number of the functional in the input dictionary

**Outputs**

*  `functional` : Returns the functional name in the dictionary

"""->
function getFunctionalName(opts, f_number)

  key = string("functional_name", f_number)
  val = opts[key]

  return functional = FunctionalDict[val]
end

