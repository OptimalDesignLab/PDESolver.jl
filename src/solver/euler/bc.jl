# functions for calculating boundary integrals

include("bc_solvers.jl")  # Roe solvers and related things

"""
  Type that identified a particular node on a particular face of
  a particular boundary.  It also contains the index of the face within
  the array of Boundaries with the same boundary condition

  **Fields**

   * element: the element the face is part of
   * face: the local face number
   * faceidx: the index within the array of Boundaries with this BC
   * node: the node on the face

"""
struct BoundaryNode
  element::UInt32
  face::UInt8
  faceidx::Int  # index within array of faces that have this boundary condition
  node::Int
end

function BoundaryNode(bndry::Boundary, faceidx, node)
  return BoundaryNode(bndry.element, bndry.face, faceidx, node)
end

"""
  Null boundary node (all fields zero).  This is useful as a default value
  for a functiona argument (if the argument is unused).
"""
global const NullBoundaryNode = BoundaryNode(Boundary(0, 0), 0, 0)

@doc """
### EulerEquationMod.getBCFluxes

  This function calls other functions to calculate the boundary fluxes, passing
  them pieces of the array needed.  This populates eqn.bndryflux.  It also
  calls writeBoundary() to do any requested output.  If the options dictionary
  specifies not to precompute the boundary flux, this function will do the
  integration as well and update `eqn.res`.

  This is a mid level function.

  
"""
# this is a mid level function
function getBCFluxes(mesh::AbstractMesh, sbp::AbstractOperator, eqn::EulerData, opts)
  #get all the fluxes for all the boundary conditions and save them in eqn.bndryflux

  #println("mesh.bndry_funcs = ", mesh.bndry_funcs)
  for i=1:mesh.numBC
  #  println("computing flux for boundary condition ", i)
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range = start_index:(end_index - 1)
    bndry_facenums_i = sview(mesh.bndryfaces, idx_range)

    if opts["use_staggered_grid"]

      calcBoundaryFlux_nopre(mesh, mesh.mesh2, sbp, mesh.sbp2, eqn, functor_i,
                             idx_range, bndry_facenums_i)
      
    elseif opts["precompute_boundary_flux"]
      bndryflux_i = sview(eqn.bndryflux, :, :, start_index:(end_index - 1))

      # call the function that calculates the flux for this boundary condition
      # passing the functor into another function avoid type instability
      calcBoundaryFlux(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i, bndryflux_i)
    else
      calcBoundaryFlux_nopre(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i)
    end
  end

  writeBoundary(mesh, sbp, eqn, opts)

  return nothing
end

@doc """
### EulerEquationMod.writeBoundary

  This function writes information about the boundary faces and fluxes to files.
  It is controlled by the input argument writeboundary, of type Bool.

  It generates the files:
    * boundaryfaces.dat : writes mesh.bndryfaces, an array with eltype Boundary
                          to a file, one element per line
    * boundaryflux.dat  : writes the element, local node number and boundary
                          flux to a line in a human readable format
    * boundaryflux2.dat : writes the real part ofmesh.bndryflux to space
                          delimited file

   This is a high level function.
"""
function writeBoundary(mesh, sbp, eqn, opts)

  if !eqn.params.writeboundary
    return nothing
  end

    face_name = "boundaryfaces.dat"
    flux_name = "boundaryflux.dat"
    flux_dlm = "boundaryflux2.dat"

    rmfile(face_name)
    rmfile(flux_name)
    rmfile(flux_dlm)


  # write boundaryfaces.dat
  f = open(face_name, "a+")
  for i=1:length(mesh.bndryfaces)
    println(f, mesh.bndryfaces[i])
  end
  close(f)

  # write boundaryflux.dat
  f = open(flux_name, "a+")
  for i=1:mesh.numBoundaryFaces
    el = mesh.bndryfaces[i].element
    face = mesh.bndryfaces[i].face
    for j=1:mesh.numNodesPerFace  # TODO: should be mesh.numNodesPerFace?
      jb = mesh.facenodes[j, face]
      println(f, "el ", el, ", node_index ", jb, ", flux = ",
               real(eqn.bndryflux[:, j, i]))
    end
  end
  close(f)

  # write boundaryflux2.dat
  writedlm(flux_dlm, real(eqn.bndryflux))

  return nothing
end

@doc """
### EulerEquationMod.interpolateBoundary

  Interpolates the solution variables to the exterior boundary of the mesh
  and calculates any additional quantities at the boundary of the mesh.
  DG only

  Inputs:
    mesh: an AbstractDGMesh
    sbp
    eqn
    opts
    q : the 3D array of solution variables for all elements, numdofpernode x
        numNodesPerElement x numEl

  Inputs/Outputs:
    q_bndry: the array to be populated with the solution interpolated to
             the boundary, numdofpernode x numNodesPerFace x num boundary faces

    eqn.aux_vars_bndry is also populated

    Aliasing restrictions: none
"""
function interpolateBoundary(mesh::AbstractDGMesh, sbp, eqn, opts, q::Abstract3DArray, q_bndry::Abstract3DArray, aux_vars_bndry::Abstract3DArray)

  # interpolate solutions
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, q, q_bndry)

  # calculate aux_vars
  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      q_vals = ro_sview(eqn.q_bndry, :, j, i)
      aux_vars_bndry[1, j, i] = calcPressure(eqn.params, q_vals)
    end
  end

end


@doc """
### EulerEquationMod.calcBoundaryFlux

  This function calculates the boundary flux for the portion of the boundary
  with a particular boundary condition.  The eqn.q are converted to
  conservative variables if needed

  Inputs:
  mesh : AbstractMesh
  sbp : AbstractOperator
  eqn : EulerEquation
  functor : a callable object that calculates the boundary flux at a node
  idx_range: the Range describing which Boundaries have the current BC
  bndry_facenums:  An array with elements of type Boundary that tell which
                   element faces have the boundary condition
  Outputs:
  bndryflux : the array to store the boundary flux, corresponds to
              bndry_facenums

  The functor must have the signature
  functor( q, aux_vars, x, nrm_xy, bndryflux_i, eqn.params)
  where q are the *conservative* variables.
  where all arguments (except params) are vectors of values at a node.

  params is the ParamType associated with the the EulerEquation object
  nrm = mesh.sbpface.normal[:, current_node]

  This is a mid level function.
"""
# mid level function
function calcBoundaryFlux( mesh::AbstractCGMesh{Tmsh},
       sbp::AbstractOperator, eqn::EulerData{Tsol},
       functor::BCType, idx_range::UnitRange,
       bndry_facenums::AbstractArray{Boundary,1},
       bndryflux::AbstractArray{Tres, 3}) where {Tmsh,  Tsol, Tres}
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("enterted calcBoundaryFlux")


  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  nrm_xy = zeros(Tmsh, mesh.dim)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
#    println("boundary ", i, "element = ", bndry_i.element, ", face = ", bndry_i.face)
#    println("interface ", i)
    for j = 1:mesh.numNodesPerFace
      k = mesh.facenodes[j, bndry_i.face]

      # get components
      q = ro_sview(eqn.q, :, k, bndry_i.element)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = ro_sview(eqn.aux_vars, :, k, bndry_i.element)
      x = ro_sview(mesh.coords, :, k, bndry_i.element)
      dxidx = ro_sview(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = ro_sview(mesh.sbpface.normal, :, bndry_i.face)
      calcBCNormal(eqn.params, dxidx, nrm, nrm_xy)
      bndryflux_i = sview(bndryflux, :, j, i)

      functor(eqn.params, q2, aux_vars, x, nrm_xy, bndryflux_i)

    end

  end


  return nothing
end


# DG version
function calcBoundaryFlux( mesh::AbstractDGMesh{Tmsh},
       sbp::AbstractOperator, eqn::EulerData{Tsol},
       functor::BCType, idx_range::UnitRange,
       bndry_facenums::AbstractArray{Boundary,1},
       bndryflux::AbstractArray{Tres, 3}) where {Tmsh,  Tsol, Tres}
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  params = eqn.params
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]

    for j = 1:mesh.numNodesPerFace
      # get components
      q = ro_sview(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
      bndryflux_i = sview(bndryflux, :, j, i)

      bndry_node = BoundaryNode(bndry_i, i, j)

      functor(params, q2, aux_vars, coords, nrm_xy, bndryflux_i, bndry_node)

    end
  end

  return nothing
end

"""
  Like calcBoundaryFlux, but performs the integration and updates res rather
  than storing the flux.
"""
function calcBoundaryFlux_nopre( mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractOperator, eqn::EulerData{Tsol, Tres},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1}) where {Tmsh,  Tsol, Tres}
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  params = eqn.params
  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace

      # get components
      #TODO: this doesn't work if precompute_q_bndr == false ?
      q = ro_sview(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
      bndryflux_i = sview(flux_face, :, j)

     
      bndry_node = BoundaryNode(bndry_i, i, j)

      functor(params, q2, aux_vars, coords, nrm_xy, bndryflux_i, bndry_node)
    end

    res_i = sview(eqn.res, :, :, bndry_i.element)
    boundaryFaceIntegrate!(mesh.sbpface, bndry_i.face, flux_face, res_i,
                           SummationByParts.Subtract())
  end

  return nothing
end


"""
  Staggered grid version
"""
function calcBoundaryFlux_nopre(mesh_s::AbstractDGMesh{Tmsh},
                          mesh_f::AbstractDGMesh{Tmsh},
                          sbp_s::AbstractOperator, sbp_f::AbstractOperator,
                          eqn::EulerData{Tsol, Tres},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1}) where {Tmsh,  Tsol, Tres}
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  nfaces = length(bndry_facenums)
  q_face = zeros(Tsol, mesh_f.numDofPerNode, mesh_f.numNodesPerFace)
  params = eqn.params
  flux_face = zeros(Tres, mesh_f.numDofPerNode, mesh_f.numNodesPerFace)
  res_f = zeros(Tres, mesh_f.numDofPerNode, mesh_f.numNodesPerElement)
  res_s = zeros(Tres, mesh_s.numDofPerNode, mesh_s.numNodesPerElement)
  aux_vars = zeros(Tsol, 1)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    q_vol = ro_sview(eqn.q_flux, :, :, bndry_i.element)
    boundaryFaceInterpolate!(mesh_f.sbpface, bndry_i.face, q_vol, q_face)

    # interpolate to face
    for j = 1:mesh_f.numNodesPerFace

      # get components
      q_j = sview(q_face, :, j)

      # convert to conservative variables if needed
      aux_vars[1] = calcPressure(eqn.params, q_j)
      coords = ro_sview(mesh_f.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh_f.nrm_bndry, :, j, global_facenum)
      bndryflux_i = sview(flux_face, :, j)

      functor(params, q_j, aux_vars, coords, nrm_xy, bndryflux_i)
    end  # end loop j

    # integrate
    fill!(res_f, 0.0)
    fill!(res_s, 0.0)
    boundaryFaceIntegrate!(mesh_f.sbpface, bndry_i.face, flux_face, res_f,
                           SummationByParts.Subtract())

    # interpolate back
    smallmatmat!(res_f, mesh_s.I_S2F, res_s)


    # accumulate into res
    @simd for j=1:mesh_s.numNodesPerElement
      @simd for k=1:mesh_s.numDofPerNode
        eqn.res[k, j, bndry_i.element] += res_s[k, j]
      end
    end

  end  # end loop i

  return nothing
end

#-----------------------------------------------------------------------------
# Boundary condition definitions
# TODO: move these to separate file

"""
  This macro facilitates defining a new boundary condition struct and default
  constructor.  This works for regular boundary condition functors as well
  as reverse mode ones.  The default functor for a boundary condition is:

  ```
    mutable struct Foo{Tsol} <: BCType
      qg::Vector{Tsol}  # for boundary state
      v_vals::Vector{Tsol}  # for entropy variables (if needed)
    end
  ```

  **Inputs**

   * fname: name of new boundary condition struct (`Foo` in the above example)
            If the functor does the reverse mode with respect to the metrics,
            the fname *must* end with `_revm`.  Similarly, reverse mode with
            respect to the solution must end with `_revq`.
   * docstring: (optional) docstring for the newly created functor

   **Example usage**

   ```
     @makeBC Foo \"\"\"
Docstring for Foo
\"\"\"
  ```
"""
macro makeBC(fname::Symbol, docstring="")

  # introducing an empty docstring is not the same as not supplying a docstring
  # figure out which to do
  if docstring != ""
    docex = quote
              @doc $docstring $fname
            end
  else
    docex = :()
  end

  # figure out name of supertype based on name of functor
  fname_str = string(fname)
  if endswith(fname_str, "_revq")
    stype = :BCType_revq
  elseif endswith(fname_str, "_revm")
    stype = :BCType_revm
  else
    stype = :BCType
  end


  return esc(quote
               struct $fname{Tsol} <: $stype
                 #=
                 """
                   Storage for boundary state
                 """
                 qg::Vector{Tsol}
                 """
                   Storage for entropy variables (may be unneeded)
                 """
                 v_vals::Vector{Tsol}
                 =#
               end


               function $fname(mesh::AbstractMesh, eqn::EulerData{Tsol, Tres}) where {Tsol, Tres}
                 #=
                 qg = zeros(Tsol, mesh.numDofPerNode)
                 v_vals = zeros(Tsol, mesh.numDofPerNode)
                 =#

                 return $fname{Tsol}()
               end

               $docex
             end  # end quote
            )
end



@makeBC errorBC """
  This BC generates an error if it is called
"""

function (obj::errorBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  error("errorBC has been called")
end


"""
  Returns the boundary state for a given boundary condition.  Non dirchlet-type
  boundary conditions should throw an exception
"""
function getDirichletState(obj::errorBC, params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  error("getDirichletState for errorBC has been called")
end


@makeBC errorBC_revm """
  This BC generates an error if it is called.  Used as a placeholder if no
  BC is specified
"""
function (obj::errorBC_revm)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm::AbstractArray{Tmsh,1},
              nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  error("errorBC_revm has been called")
end

@makeBC errorBC_revq """
  This BC generates an error if it is called.  Used as a placeholder if no
  BC is specified
"""
function (obj::errorBC_revq)(params::ParamType2,
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tsol, 1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  error("errorBC_revm has been called")
end


@makeBC isentropicVortexBC """
### EulerEquationMod.isentropicVortexBC <: BCTypes

  This type and the associated call method define a functor to calculate
  the flux using the Roe Solver using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.
"""

# low level function
function (obj::isentropicVortexBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.bcdata.qg
  getDirichletState(obj, params, q, aux_vars, coords, nrm, qg, bndry)
  RoeSolver(params, q, qg, aux_vars, nrm, bndryflux)

  return nothing

end # ends the function isentropicVortexBC


@makeBC isentropicVortexBC_revm """
###EulerEquationMod.isentropicVortexBC_revm

Reverse mode for isentropicVortexBC.

"""


function (obj::isentropicVortexBC_revm)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm::AbstractArray{Tmsh,1},
              nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


   # Forward Sweep
  @unpack params.bcdata qg q_bar qg_bar
  fill!(q_bar, 0); fill!(qg_bar, 0)

  getDirichletState(obj, params, q, aux_vars, coords, nrm, qg, bndry)

  # Reverse Sweep

  # RoeSolver(params, q, qg, aux_vars, nrm, bndryflux)
  RoeSolver_revm(params, q, qg, aux_vars, nrm, nrm_bar, bndryflux_bar)
  RoeSolver_revq(params, q, q_bar, qg, qg_bar, aux_vars, nrm, bndryflux_bar)

  getDirichletState_revm(obj, params, q, aux_vars, coords, coords_bar, nrm,
                         nrm_bar, qg_bar, bndry)

  return nothing
end


@makeBC isentropicVortexBC_revq """
###EulerEquationMod.isentropicVortexBC_revq

Reverse mode for isentropicVortexBC.

"""


function (obj::isentropicVortexBC_revq)(params::ParamType2,
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tsol, 1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  # Forward Sweep
  @unpack params.bcdata qg qg_bar
  fill!(qg_bar, 0)

  getDirichletState(obj, params, q, aux_vars, coords, nrm, qg)
  calcIsentropicVortex(params, coords, qg)

  # Reverse Sweep

  # RoeSolver(params, q, qg, aux_vars, nrm, bndryflux)
  RoeSolver_revq(params, q, q_bar, qg, qg_bar, aux_vars, nrm, bndryflux_bar)

  return nothing
end

const isentropicVortexBCs = Union{isentropicVortexBC, isentropicVortexBC_revq,
                                  isentropicVortexBC_revm}
function getDirichletState(obj::isentropicVortexBCs, params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  calcIsentropicVortex(params, coords, qg)
end


function getDirichletState_revq(obj::isentropicVortexBC_revq, params::ParamType2,
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tsol, 1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              qg_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # nothing to do here
  return nothing
end


function getDirichletState_revm(obj::isentropicVortexBC_revm, params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm::AbstractArray{Tmsh,1},
              nrm_bar::AbstractVector{Tmsh},
              qg_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  calcIsentropicVortex_rev(params, coords, coords_bar, qg_bar)
end


@makeBC isentropicVortexBC_physical """
### EulerEquationMod.isentropicVortexBC_physical <: BCTypes

  This type and the associated call method define a functor to calculate
  the actual Euler flux  using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.
"""


function (obj::isentropicVortexBC_physical)(
              params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  calcEulerFlux(params, q, aux_vars, nrm_xy, bndryflux)

  return nothing

end # end function isentropicVortexBC_physical


@makeBC noPenetrationBC """
### EulerEquationMod.noPenetrationBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid velocity is projected into the wall.

  Works in 2D, untested in 3D.

  This is a low level functor
"""

# low level function
function (obj::noPenetrationBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.bcdata.qg
  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)
  v_vals = params.bcdata.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)
  # this is a problem: q is in conservative variables even if
  # params says we are using entropy variables


#  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)
#  calcLFFlux(params, q, v_vals, aux_vars, nrm_xy, bndryflux)
  calcEulerFlux(params, v_vals, aux_vars, nrm_xy, bndryflux)

  return nothing
end


function (obj::noPenetrationBC)(params::ParamType3,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector


  qg = params.bcdata.qg

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)

  # call Roe solver
  #RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  v_vals = params.bcdata.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)
  # this is a problem: q is in conservative variables even if
  # params says we are using entropy variables
  calcEulerFlux(params, v_vals, aux_vars, nrm_xy, bndryflux)
#  calcLFFlux(params, q, v_vals, aux_vars, nrm_xy, bndryflux)

  return nothing
end

@makeBC noPenetrationESBC """
  Entropy stable no penetration boundary condition from Chen and Shu
  "Entropy Stbale High Order Discontinuous Galerkin Methods with Suitable
  quadrature rules for Hyperbolic Conservation Laws.

  Entropy stable for diagonal E SBP operators only
"""

# low level function
function (obj::noPenetrationESBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector


  qg = params.bcdata.qg
  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)
  calcLFFlux(params, q, qg, aux_vars,nrm_xy, bndryflux)
  #calcHLLFlux(params, q, qg, aux_vars,nrm_xy, bndryflux)

  return nothing
end


# revm
@makeBC noPenetrationESBC_revm """
Reverse mode for noPenetrationESBC.

"""
function (obj::noPenetrationESBC_revm)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  @unpack params.bcdata qg q_bar qg_bar
  fill!(q_bar, 0); fill!(qg_bar, 0)

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)

  # Reverse sweep
  calcLFFlux_revm(params, q, qg, aux_vars, nrm_xy, nrm_bar, bndryflux_bar)
  calcLFFlux_revq(params, q, q_bar, qg, qg_bar, aux_vars, nrm_xy, 
                  bndryflux_bar)

  getDirichletState_revm(obj, params, q, aux_vars, coords, coords_bar, nrm_xy,
                         nrm_bar, qg_bar, bndry)

  return nothing
end



# revq

@makeBC noPenetrationESBC_revq """
  Reverse mode wrt q of `noPenetrationESBC`
"""

function (obj::noPenetrationESBC_revq)(params::ParamType{2, :conservative},
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  @unpack params.bcdata qg qg_bar
  fill!(qg_bar, 0)

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)

  # Reverse sweep
  calcLFFlux_revq(params, q, q_bar, qg, qg_bar, aux_vars, nrm_xy, 
                  bndryflux_bar)

  getDirichletState_revq(obj, params, q, q_bar, aux_vars, coords, nrm_xy,
                         qg_bar, bndry)

  return nothing
end


# getDirichletState

const noPenetrationESBCs = Union{noPenetrationESBC, noPenetrationESBC_revm,
                                 noPenetrationESBC_revq}

function getDirichletState(obj::noPenetrationESBCs, params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  # Get the normal momentum
  Unrm = nx*q[2] + ny*q[3]

  # this is equivalent to:
  #   1. computing the normal and tangential components
  #   2. negating the normal component
  #   3. combining the negative normal and non-negated tangent component
  qg[1] = q[1]
  qg[2] = -2*Unrm*nx + q[2]
  qg[3] = -2*Unrm*ny + q[3]
  qg[4] = q[4]

  return nothing
end


function getDirichletState(obj::noPenetrationESBCs, params::ParamType3,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  nz = nrm_xy[3]
  fac = 1.0/(sqrt(nx*nx + ny*ny + nz*nz))
  # normalize normal vector
  nx = nx*fac
  ny = ny*fac
  nz = nz*fac

  # this is momentum, not velocity?
  Unrm = nx*q[2] + ny*q[3] + nz*q[4]

  qg[1] = q[1]
  qg[2] = -2*Unrm*nx + q[2]
  qg[3] = -2*Unrm*ny + q[3]
  qg[4] = -2*Unrm*nz + q[4]
  qg[5] = q[5]

  return nothing
end


function getDirichletState_revm(obj::noPenetrationESBC_revm, params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              qg_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  # Get the normal momentum
  Unrm = nx*q[2] + ny*q[3]

  # this is equivalent to:
  #   1. computing the normal and tangential components
  #   2. negating the normal component
  #   3. combining the negative normal and non-negated tangent component
  #qg[1] = q[1]
  #qg[2] = -2*Unrm*nx + q[2]
  #qg[3] = -2*Unrm*ny + q[3]
  #qg[4] = q[4]

  #-----------------
  # reverse sweep
  nx_bar = -2*Unrm*qg_bar[2]
  ny_bar = -2*Unrm*qg_bar[3]
  Unrm_bar = -2*nx*qg_bar[2] + -2*ny*qg_bar[3]

  nx_bar += q[2]*Unrm_bar
  ny_bar += q[3]*Unrm_bar

  # restore primal state
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  fac_bar = nx*nx_bar + ny*ny_bar
  nx_bar *= fac
  ny_bar *= fac

  den = (nx*nx + ny*ny)^1.5
  nx_bar -= fac_bar*nx/den
  ny_bar -= fac_bar*ny/den

  nrm_bar[1] += nx_bar
  nrm_bar[2] += ny_bar

  return nothing
end

function getDirichletState_revq(obj::noPenetrationESBC_revq, 
              params::ParamType{2, :conservative},
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  # Get the normal momentum
  Unrm = nx*q[2] + ny*q[3]

  # this is equivalent to:
  #   1. computing the normal and tangential components
  #   2. negating the normal component
  #   3. combining the negative normal and non-negated tangent component
  #qg[1] = q[1]
  #qg[2] = -2*Unrm*nx + q[2]
  #qg[3] = -2*Unrm*ny + q[3]
  #qg[4] = q[4]

  #-----------------
  # reverse sweep
  q_bar[1] += qg_bar[1]
  q_bar[2] += qg_bar[2]
  q_bar[3] += qg_bar[3]
  q_bar[4] += qg_bar[4]

  Unrm_bar = -2*nx*qg_bar[2] + -2*ny*qg_bar[3]
  q_bar[2] += nx*Unrm_bar
  q_bar[3] += ny*Unrm_bar

  return nothing
end



@makeBC noPenetrationBC_revm """
###EulerEquationMod.noPenetrationBC_revm

Reverse mode for noPenetrationBC.

"""
function (obj::noPenetrationBC_revm)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.bcdata.qg
  getDirichletState(obj, params, q, aux_vars, coords, nrm, qg, bndry)
  v_vals = params.bcdata.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)

  # Reverse sweep
#  nrm2_bar = zeros(Tmsh, 2)
  qg_bar = zeros(Tsol, 4)  #TODO: stop allocating new arrays
  calcEulerFlux_revm(params, v_vals, aux_vars, nrm, nrm_bar, bndryflux_bar)
  calcEulerFlux_revq(params, v_vals, qg_bar, aux_vars, nrm, bndryflux_bar)

  # TODO: reverse mode convertFromNaturalToWorkingVars(params, qg, v_vals)
  getDirichletState_revm(obj, params, q, aux_vars, coords, coords_bar, nrm,
                         nrm_bar, qg_bar, bndry)

  return nothing
end





@makeBC noPenetrationBC_revq """
  Reverse mode wrt q of `noPenetrationBC`
"""

function (obj::noPenetrationBC_revq)(params::ParamType{2, :conservative},
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector

  data = params.bcdata
  @unpack data qg qg_bar
  fill!(qg_bar, 0)

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)
  calcEulerFlux_revq(params, qg, qg_bar, aux_vars, nrm_xy, bndryflux_bar)
  getDirichletState_revq(obj, params, q, q_bar, aux_vars, coords, nrm_xy, qg_bar, bndry)

  return nothing
end

# boundary conditions really need to have a single type for all 3 versions
# so this kind of thing is not needed
const noPenetrationBCs = Union{noPenetrationBC, noPenetrationBC_revq, noPenetrationBC_revm}
function getDirichletState(obj::noPenetrationBCs, params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector

  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  # Get the normal momentum
  Unrm = nx*q[2] + ny*q[3]

  for i=1:length(q)
    qg[i] = q[i]
  end

  # Subtract the normal component of the momentum from \xi & \eta components
  # of the momentum
  qg[2] -= nx*Unrm
  qg[3] -= ny*Unrm

  return nothing
end


function getDirichletState(obj::noPenetrationBCs, params::ParamType3,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector

  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  nz = nrm_xy[3]
  fac = 1.0/(sqrt(nx*nx + ny*ny + nz*nz))
  # normalize normal vector
  nx = nx*fac
  ny = ny*fac
  nz = nz*fac

  # this is momentum, not velocity?
  Unrm = nx*q[2] + ny*q[3] + nz*q[4]

  for i=1:length(q)
    qg[i] = q[i]
  end

  # calculate normal velocity
  qg[2] -= nx*Unrm
  qg[3] -= ny*Unrm
  qg[4] -= nz*Unrm

  return nothing
end




function getDirichletState_revm(obj::noPenetrationBC_revm, params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              qg_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
  # Forward sweep
  n1 = nrm[1]
  n2 = nrm[2]
#  n1 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
#  n2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  fac = 1.0/(sqrt(n1*n1 + n2*n2))
  nx = n1*fac
  ny = n2*fac
  Unrm = nx*q[2] + ny*q[3]

  # Reverse sweep
  n1_bar = nrm_bar[1]
  n2_bar = nrm_bar[2]


  # qg[2] = qg[2] - nx*Unrm
  # qg[3] = qg[3] - ny*Unrm
  ny_bar = -qg_bar[3]*Unrm
  nx_bar = -qg_bar[2]*Unrm
  Unrm_bar = -qg_bar[3]*ny -qg_bar[2]*nx

  # Unrm = nx*q[2] + ny*q[3]
  nx_bar += Unrm_bar*q[2]
  ny_bar += Unrm_bar*q[3]

  # nx = n1*fac
  # ny = n2*fac
  n2_bar += ny_bar*fac
  n1_bar += nx_bar*fac
  fac_bar = ny_bar*n2 + nx_bar*n1

  # fac = 1.0/(sqrt(n1*n1 + n2*n2))
  n1_bar -= fac_bar*n1*((n1*n1 + n2*n2)^(-1.5))
  n2_bar -= fac_bar*n2*((n1*n1 + n2*n2)^(-1.5))

  # because n1_bar and n2_bar came out of nrm_bar earlier, this is an assigment
  nrm_bar[1] = n1_bar
  nrm_bar[2] = n2_bar

  return nothing
end


function getDirichletState_revq(obj::noPenetrationBC_revq, 
              params::ParamType{2, :conservative},
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  # Get the normal momentum
  Unrm = nx*q[2] + ny*q[3]

  # reverse sweep
  Unrm_bar = -nx*qg_bar[2] - ny*qg_bar[3]
  for i=1:length(q)
    q_bar[i] += qg_bar[i]
  end

  q_bar[2] += nx*Unrm_bar
  q_bar[3] += ny*Unrm_bar

  return nothing
end


#=
function (obj::noPenetrationBC_revm)(params::ParamType3,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, dxidx_bar::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward sweep
  nx = zero(Tmsh)
  ny = zero(Tmsh)
  nz = zero(Tmsh)
  nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
  ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
  nz2 = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]
  fac = 1.0/(sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2))
  # normalize normal vector
  nx = nx2*fac
  ny = ny2*fac
  nz = nz2*fac

  # this is momentum, not velocity?
  Unrm = nx*q[2] + ny*q[3] + nz*q[4]

  qg = params.bcdata.qg
  for i=1:length(q)
    qg[i] = q[i]
  end

  # calculate normal velocity
  qg[2] -= nx*Unrm
  qg[3] -= ny*Unrm
  qg[4] -= nz*Unrm

  v_vals = params.bcdata.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)

  # Reverse Sweep
  nrm2_bar = zeros(Tmsh, 3)
  v_vals_bar = zeros(Tsol, 5)
  calcEulerFlux_revm(params, v_vals, aux_vars, [nx2, ny2, nz2], nrm2_bar, bndryflux_bar)
  calcEulerFlux_revq(params, v_vals, v_vals_bar, aux_vars, [nx2, ny2, nz2], bndryflux_bar)

  nz2_bar = nrm2_bar[3]
  ny2_bar = nrm2_bar[2]
  nx2_bar = nrm2_bar[1]

  # TODO: reverse mode convertFromNaturalToWorkingVars(params, qg, v_vals)
  qg_bar = v_vals_bar

  # qg[4] = qg[4] - nz*Unrm
  nz_bar = -qg_bar[4]*Unrm
  Unrm_bar = -qg_bar[4]*nz
  # qg[3] -= ny*Unrm
  ny_bar = -qg_bar[3]*Unrm
  Unrm_bar -= qg_bar[3]*ny
  # qg[2] -= nx*Unrm
  nx_bar = -qg_bar[2]*Unrm
  Unrm_bar -= qg_bar[2]*nx

  # Unrm = nx*q[2] + ny*q[3] + nz*q[4]
  nx_bar += Unrm_bar*q[2]
  ny_bar += Unrm_bar*q[3]
  nz_bar += Unrm_bar*q[4]

  # nz = nz2*fac
  nz2_bar += nz_bar*fac
  fac_bar = nz_bar*nz2
  # ny = ny2*fac
  ny2_bar += ny_bar*fac
  fac_bar += ny_bar*ny2
  # nx = nx2*fac
  nx2_bar += nx_bar*fac
  fac_bar += nx_bar*nx2

  # fac = 1.0/(sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2))
  nx2_bar -= fac_bar*nx2*((nx2*nx2 + ny2*ny2 + nz2*nz2)^(-1.5))
  ny2_bar -= fac_bar*ny2*((nx2*nx2 + ny2*ny2 + nz2*nz2)^(-1.5))
  nz2_bar -= fac_bar*nz2*((nx2*nx2 + ny2*ny2 + nz2*nz2)^(-1.5))

  # nz2 = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]
  dxidx_bar[1,3] += nz2_bar*nrm[1]
  dxidx_bar[2,3] += nz2_bar*nrm[2]
  dxidx_bar[3,3] += nz2_bar*nrm[3]
  # ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
  dxidx_bar[1,2] += ny2_bar*nrm[1]
  dxidx_bar[2,2] += ny2_bar*nrm[2]
  dxidx_bar[3,2] += ny2_bar*nrm[3]
  # nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
  dxidx_bar[1,1] += nx2_bar*nrm[1]
  dxidx_bar[2,1] += nx2_bar*nrm[2]
  dxidx_bar[3,1] += nx2_bar*nrm[3]

  return nothing
end # End noPenetrationBC_revm 3D
=#
@makeBC unsteadyVortexBC """
### EulerEquationMod.unsteadyVortexBC <: BCTypes

  This type and the associated call method define a functor to calculate
  the flux using the Roe Solver using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.
"""
# low level function
function (obj::unsteadyVortexBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)
  # getting qg
  qg = params.bcdata.qg
  calcUnsteadyVortex(params, coords, qg)

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end # ends the function unsteadyVortex BC

@makeBC unsteadyVortex2BC

# low level function
function (obj::unsteadyVortex2BC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)
  # getting qg
  qg = params.bcdata.qg
  calcUnsteadyVortex2(params, coords, qg)

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end # ends the function unsteadyVortex BC


@makeBC Rho1E2U1VW0BC """
### EulerEquationMod.Rho1E2U1VW0BC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid density is 1, energy = 2, and 
    u = 1.0, v = 0.0, and w = 0.0 (if 3D).

  It should work for 2D and 3D meshes.

  This is a low level functor.

"""
function (obj::Rho1E2U1VW0BC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  #println("in Rho1E2BCU1V0W0")
  qg = params.bcdata.qg

  calcRho1Energy2U1VW0(params, coords, qg)

  #println("qg = ", qg)
  # call Roe solver
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end

@makeBC Rho1E2BC """
### EulerEquationMod.Rho1E2BC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid density is 1, energy = 2, and u = v = 0.0

  This is a low level functor

"""
function (obj::Rho1E2BC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  #println("in Rho1E2BC")
  qg = params.bcdata.qg

  calcRho1Energy2(params, coords, qg)

  #println("qg = ", qg)
  # call Roe solver
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end

@makeBC Rho1E2U3BC """
### EulerEquationMod.Rho1E2U3BC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid density is 1, energy = 2, and x velocity = 0.5

  This is a low level functor

"""
function (obj::Rho1E2U3BC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  #println("in Rho1E2U3Bc")
  qg = params.bcdata.qg

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)

  #println("qg = ", qg)
  # call Roe solver
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end

@makeBC Rho1E2U3BC_revm """
  Reverse mode wrt metrics of of`Rho1E2U3BC
"""
function (obj::Rho1E2U3BC_revm)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward sweep
  qg = params.bcdata.qg
  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)

  # Reverse sweep
  RoeSolver_revm(params, q, qg, aux_vars, nrm_xy, nrm_bar, bndryflux_bar)
end


@makeBC Rho1E2U3BC_revq """
  Reverse mode wrt solution of Rho1E2U3BC
"""
function (obj::Rho1E2U3BC_revq)(params::ParamType,
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  data = params.bcdata
  @unpack data qg qg_bar
  fill!(qg_bar, 0)

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)
  RoeSolver_revq(params, q, q_bar,  qg, qg_bar, aux_vars, nrm_xy, bndryflux_bar)

  return nothing
end


const Rho1E2U3BCs = Union{Rho1E2U3BC, Rho1E2U3BC_revm, Rho1E2U3BC_revq}
function getDirichletState(obj::Rho1E2U3BCs, params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  calcRho1Energy2U3(params, coords, qg)
end


function getDirichletState_revq(obj::Rho1E2U3BC_revq, params::ParamType,
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # nothing to do here
end


function getDirichletState_revm(obj::Rho1E2U3BC_revm, params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              qg_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # nothing to do here
end


@makeBC FreeStreamBC """
### EulerEquationMod.FreeStreamBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state corresponding to the free stream velocity, using rho_free, Ma, aoa, and E_free

  Works in 2D and 3D

  This is a low level functor

"""
function (obj::FreeStreamBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.bcdata.qg

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end


@makeBC FreeStreamBC_revm """
###EulerEquationMod.FreeStreamBC_revm

Reverse mode for FreeStreamBC.

"""
function (obj::FreeStreamBC_revm)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward sweep
  qg = params.bcdata.qg

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)
  # Reverse sweep
  RoeSolver_revm(params, q, qg, aux_vars, nrm_xy, nrm_bar, bndryflux_bar)

  # calcFreeStream does not depend on coords, so no need to reverse mode it
  return nothing
end

@makeBC FreeStreamBC_revq """
  Reverse mode wrt q of `FreeStreamBC`
"""

function (obj::FreeStreamBC_revq)(params::ParamType,
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  data = params.bcdata
  @unpack data qg qg_bar
  fill!(qg_bar, 0)

  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg, bndry)
  RoeSolver_revq(params, q, q_bar,  qg, qg_bar, aux_vars, nrm_xy, bndryflux_bar)

  return nothing
end

const FreeStreamBCs = Union{FreeStreamBC, FreeStreamBC_revq, FreeStreamBC_revm}
function getDirichletState(obj::FreeStreamBCs, params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  calcFreeStream(params, coords, qg)

  return nothing
end


function getDirichletState_revm(obj::FreeStreamBC_revm, params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  # calcFreeStream does not depend on coords, so no need to reverse mode it

  return nothing
end


function getDirichletState_revq(obj::FreeStreamBC_revq, params::ParamType,
              q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # no dependence on q
  return nothing
end




@makeBC FreeStreamBC_dAlpha """
### EulerEquationMod.FreeStreamBC_dAlpha <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state corresponding to the free stream velocity, using rho_free, Ma, aoa, and E_free

  This is a low level functor
"""
function (obj::FreeStreamBC_dAlpha)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.bcdata.qg

  calcFreeStream_dAlpha(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm, bndryflux)

  return nothing
end


@makeBC allOnesBC """
### EulerEquationMod.allOnesBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where all the conservative variables have a value 1.0

  This is a low level functor
"""
function (obj::allOnesBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = zeros(Tsol, 4)
  calcOnes(params, coords, qg)

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function call

@makeBC allZerosBC """
### EulerEquationMod.allZerosBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where all the conservative variables have a value 0.0

  This is a low level functor
"""

function (obj::allZerosBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = zeros(Tsol, 4)
  calcZeros(params, coords, qg)

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function


@makeBC ExpBC """
  Dirichlet BC for [`calcExp`](@ref)
"""

function (obj::ExpBC)(params::ParamType, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.bcdata.qg
  calcExp(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function

@makeBC ExpBC_revm """
Reverse mode of ExpBC
"""

function (obj::ExpBC_revm)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1}, coords_bar::AbstractArray{Tmsh, 1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward Sweep
  @unpack params.bcdata qg q_bar qg_bar
  fill!(q_bar, 0); fill!(qg_bar, 0)

  calcExp(params, coords, qg)

  # Reverse Sweep

  # RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)
  RoeSolver_revm(params, q, qg, aux_vars, nrm_xy, nrm_bar, bndryflux_bar)
  RoeSolver_revq(params, q, q_bar, qg, qg_bar, aux_vars, nrm_xy, bndryflux_bar)

  calcExp_rev(params, coords, coords_bar, qg_bar)

  return nothing
end


@makeBC ExpBC_revq """
Reverse mode of ExpBC
"""


function (obj::ExpBC_revq)(params::ParamType,
              q::AbstractArray{Tsol,1},
              q_bar::AbstractArray{Tres, 1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward Sweep
  @unpack params.bcdata qg qg_bar
  fill!(qg_bar, 0)

  calcExp(params, coords, qg)

  # Reverse Sweep

  # RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)
  RoeSolver_revq(params, q, q_bar, qg, qg_bar, aux_vars, nrm_xy, bndryflux_bar)

  return nothing
end



@makeBC PeriodicMMSBC """
  Dirichlet boundary condition for [`calcPeriodicMMS`](@ref).  This BC
  should not normally be needed, the correct thing to do is make a periodic
  mesh
"""

function (obj::PeriodicMMSBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# use the exact solution as the boundary condition for the PeriodicMMS
# solutions

  qg = params.bcdata.qg
  calcPeriodicMMS(params, coords, qg)
  use_efix = 0
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux, use_efix)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function

@makeBC ChannelMMSBC """
  DirichletBC for [`calcChannelMMS`](@ref)
"""

function (obj::ChannelMMSBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# use the exact solution as the boundary condition for the ChannelMMS
# solutions

  qg = params.bcdata.qg
  calcChannelMMS(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end # end function


@makeBC defaultBC """
  Default boundary condition for Euler.  Computes the Euler flux.  This 
  should only be used for supersonic outflow.
"""

function (obj::defaultBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  calcEulerFlux(params, q, aux_vars, nrm_xy, bndryflux)

  return nothing
end

@makeBC SubsonicInflowBC """
  See NASA/TM-2011-217181: Inflow/Outflow Boundary Conditions with Application
                           to FUN3D by Carlson 
  The derivation has some algebraic mistakes, but the approach is correct

  This BC does not work yet.
"""

function (obj::SubsonicInflowBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  pt = 102010.0/params.p_free  # boundary stagnation pressure
  Tt = 288.6/params.T_free  # boundary stagnation temperature
  # i = interior quantity
  # b = boundary state

  # need normalized outward normal vector
  nrm_fac = 1/sqrt(nrm_xy[1]*nrm_xy[1] + nrm_xy[2]*nrm_xy[2])

  gamma = params.gamma
  gamma_1 = params.gamma_1

  pressi = calcPressure(params, q)
  # magnitude of velocity (negative sign because the normal is outward but
  # the velocity should be inward
  Ui = -nrm_fac*(q[2]*nrm_xy[1] + q[3]*nrm_xy[2])/q[1]
#  vi = q[3]/q[1]
  ai2 = gamma*pressi/q[1]  # speed of sound squared

  # stagnation enthalpy (specific)
  hti = ai2/gamma_1 + Ui*Ui

  # Riemann invarient for the characteristic exiting the domain
  Ri = Ui - 2*sqrt(ai2)/gamma_1

  # this step uses the adiabatic assumption + the Riemann invarient Rb to
  # form a quandratic equation for ab
  # pick the larger of the two roots

  a = 0.5 + 2/gamma_1
  b = 2*Ri
  c = 0.5*gamma_1*(Ri*Ri - hti)

  tmp1 = -0.5*b/a
  tmp2 = 0.5*sqrt(b*b - 4*a*c)/a
  ab1 = tmp1 + tmp2
  ab2 = tmp1 - tmp2

  ab = max(ab1, ab2)  # maximum root is the physically correct one

  # use Riemann invarient to find velocity magnitude on the boundary side
  Ub = Ri + 2*ab/gamma_1
  Mb = Ub/ab

  @assert Mb < 1.0

  operand = 1/(1 + 0.5*gamma_1*Mb*Mb)
  pb = pt*(operand)^(gamma/gamma_1)
  Tb = Tt*operand


  # convert back to conservative variables
  qg = params.bcdata.qg
  rho1 = gamma*pb/(ab*ab)  # this is numerically equivalent to rho2 below,
                           # which is weuird because in this case Tb is never
                           # used
  rho2 = pb/(params.R_ND*Tb)
  qg[1] = rho2  # R is not nondimenstionalized
  qg[2] = -Ub*nrm_xy[1]*nrm_fac*qg[1]  # negative sign because the normal is neg
  qg[3] = -Ub*nrm_xy[2]*nrm_fac*qg[1]
  qg[4] = pb/gamma_1 + 0.5*qg[1]*Ub*Ub

#  R_computed = pb/(qg[1]*Tb)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)


  return nothing
end

@makeBC SubsonicOutflowBC

function(obj::SubsonicOutflowBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  pb = 101300.0/params.p_free  # nondimensionalized pressure
  gamma = params.gamma
  gamma_1 = params.gamma_1

  pressi = calcPressure(params, q)
  # verify Mach number < 1
  ai2 = gamma*pressi/q[1]  # speed of sound squared
  # need normalized outward normal vector
  nrm_fac = 1/sqrt(nrm_xy[1]*nrm_xy[1] + nrm_xy[2]*nrm_xy[2])
  Un = (q[2]*nrm_xy[1] + q[3]*nrm_xy[2])*nrm_fac/q[1]

  @assert Un >= 0  # this should be outflow, not inflow
  @assert (Un*Un)/ai2 < 1

  qg = params.bcdata.qg
  qg[1] = q[1]
  qg[2] = q[2]
  qg[3] = q[3]
  # compute energy from the specified pressure
  qg[4] = pb/gamma_1 + 0.5*(q[2]*q[2] + q[3]*q[3])/q[1]
#  qg[4] = pb/gamma_1 + 0.5*q[1]*(q[2]*q[2]

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end

@makeBC inviscidChannelFreeStreamBC

# low level function
function(obj::inviscidChannelFreeStreamBC)(             
              params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)
  # getting qg
  qg = params.bcdata.qg
  calcInvChannelFreeStream(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end # ends the function unsteadyVortex BC

"""
   This is a special boundary condition used for imposing a numerical
   solution as a boundary condition.

   **Fields**

    * bc_vals: an array of size numDofPerNode x numNodesPerFace x numFaces with
               this boundary condition on it.  This array can be accessed
               using the `faceidx` field of [`BoundaryNode`](@ref).


  This BC is special because it has to store the boundary values it is imposing.
  Whenever this boundary condition is needed, the user should construct a
  new object, with the data inside it, and then register it just before
  constructing the `EulerData` object.
"""
mutable struct reanalysisBC{Tsol} <: BCType
  bc_vals::Array{Tsol, 3}
end


function reanalysisBC(mesh::AbstractMesh, eqn::AbstractSolutionData{Tsol, Tres}) where {Tsol, Tres}
  bc_vals = Array{Tsol}(0, 0, 0)
  return reanalysisBC{Tsol}(bc_vals)
end

function (obj::reanalysisBC)(params::AbstractParamType{Tdim},
        q::AbstractArray{Tsol, 1},
        aux_vars::AbstractArray{Tres, 1},
        coords::AbstractVector{Tmsh},
        nrm_xy::AbstractVector{Tmsh},
        bndryflux::AbstractArray{Tres},
        bndry::EulerEquationMod.BoundaryNode) where {Tmsh, Tsol, Tres, Tdim}


  # get numerical data out of the array
  qg = sview(obj.bc_vals, :, bndry.node, bndry.faceidx)
  EulerEquationMod.RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)
  
  return nothing
end


@makeBC ZeroBC """
  Boundary condition that sets q = 0.  This is useful for testing
"""

function getDirichletState(obj::ZeroBC, params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  fill!(qg, 0)

  return nothing
end


function (obj::ZeroBC)(params::ParamType, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.bcdata.qg
  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end # end function


@makeBC LaplaceBC """
  Boundary condition for testing shock capturing diffusion terms
"""

function getDirichletState(obj::LaplaceBC, params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              qg::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  calcLaplaceSolution(params, coords, qg)

  return nothing
end


function (obj::LaplaceBC)(params::ParamType, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.bcdata.qg
  getDirichletState(obj, params, q, aux_vars, coords, nrm_xy, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end # end function




# every time a new boundary condition is created,
# add it to the dictionary

"""
  Maps boundary conditions names to the functor objects.
  Each functor should be callable with the signature
"""
global const BCDict = Dict{String, Type{T} where T <: BCType}(  # BCType
"errorBC" => errorBC,
"isentropicVortexBC" => isentropicVortexBC,
"noPenetrationBC" => noPenetrationBC,
"noPenetrationESBC" => noPenetrationESBC,
"Rho1E2BC" => Rho1E2BC,
"Rho1E2U1VW0BC" => Rho1E2U1VW0BC,
"Rho1E2U3BC" => Rho1E2U3BC,
"isentropicVortexBC_physical" => isentropicVortexBC_physical,
"FreeStreamBC" => FreeStreamBC,
"allOnesBC" => allOnesBC,
"unsteadyVortexBC" => unsteadyVortexBC,
"unsteadyVortex2BC" => unsteadyVortex2BC,
"ExpBC" => ExpBC,
"PeriodicMMSBC" => PeriodicMMSBC,
"ChannelMMSBC" => ChannelMMSBC,
"subsonicInflowBC" => SubsonicInflowBC,
"subsonicOutflowBC" => SubsonicOutflowBC,
"inviscidChannelFreeStreamBC" => inviscidChannelFreeStreamBC,
"reanalysisBC" => reanalysisBC,
"zeroBC" => ZeroBC,
"LaplaceBC" => LaplaceBC,
"defaultBC" => defaultBC,
)

@doc """
### EulerEquationMod.getBCFunctors

  This function uses the opts dictionary to populate mesh.bndry_funcs with
  the functors.

  The functors must be a subtype of [`BCType`](@ref).

  The function must be callable with the signature:

  ```
    func(params::ParamType,
         q::AbstractArray{Tsol,1},
         aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
         nrm_xy::AbstractArray{Tmsh,1},
         bndryflux::AbstractArray{Tres, 1},
         bndry::BoundaryNode=NullBoundaryNode)
  ```

  with inputs

   * q: solution vector at a node
   * aux_vars: auxiliary variables for the node
   * coords: x-y coordinates of the node
   * nrm_xy: normal vector of the node

  and output

   * bndryflux: the flux vector, length `numDofPerNode`


  The functor must also have have an outer constructor

  ```
    Foo(mesh::AbstractMesh, eqn::EulerData)
  ```

  See the [`@makeBC`](@ref) macro creating a new BC type for with the default
  fields.  If a particular boundary condition needs special fields, then it
  must be constructed by hand (ie. not using the macro).


  This is a high level function.
"""

# use this function to populate access the needed values in BCDict
function getBCFunctors(mesh::AbstractMesh, sbp::AbstractOperator, eqn::EulerData, opts)
# populate the array mesh.bndry_funcs with the functors for the boundary condition types

#  println("Entered getBCFunctors")

  for i=1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
    ctor = BCDict[val]
    mesh.bndry_funcs[i] = ctor(mesh, eqn)
  end

  return nothing
end # End function getBCFunctors

global const BCDict_revm = Dict{String, Type{T} where T <: BCType_revm}(
"errorBC" => errorBC_revm,
"noPenetrationESBC" => noPenetrationESBC_revm,
"noPenetrationBC" => noPenetrationBC_revm,
"Rho1E2U3BC" => Rho1E2U3BC_revm,
"FreeStreamBC" => FreeStreamBC_revm,
"ExpBC" => ExpBC_revm,
"isentropicVortexBC" => isentropicVortexBC_revm,
)

"""
  This function uses the options dictionary to populate mesh.bndry_funcs_revm.

  The functors must be of time [`BCType_revm`](@ref)

  If opts["need_adjoint"] is false, the error BC is supplied instead.

  The function must be callable with the signature

  ```
  func(params::ParamType,
        q::AbstractArray{Tsol,1},
        aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
        nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
        bndryflux_bar::AbstractArray{Tres, 1},
        bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
  ```

  with inputs

   * q: solution vector at a node
   * aux_vars: auxiliary variables for the node
   * coords: x-y coordinates of the node
   * nrm_xy: normal vector of the node
   * bndryflux_bar: seed vector for the flux vector, length `numDofPerNode`

  and outputs

   * nrm_bar: output vector, same length of `nrm_xy`


  The functor must also have an outer constructor

  ```
    Foo(mesh::AbstractMesh, eqn::EulerData)
  ```
"""
function getBCFunctors_revm(mesh::AbstractMesh, sbp::AbstractOperator, eqn::EulerData, opts)

  for i = 1:mesh.numBC
    key_i = string("BC", i, "_name")
    if opts["need_adjoint"]
      val = opts[key_i]
    else
      val = "errorBC"  # placeholder, because the BC in use might not have
                       # a reverse-mode version
    end
    ctor = BCDict_revm[val]
#    println("BCDict_revm[$val] = ", BCDict_revm[val])
    mesh.bndry_funcs_revm[i] = ctor(mesh, eqn)
  end # End for i = 1:mesh.numBC

  return nothing
end # End function getBCFunctors_revm


global const BCDict_revq = Dict{String, Type{T} where T <: BCType_revq}(
"errorBC" => errorBC_revq,
"noPenetrationESBC" => noPenetrationESBC_revq,
"noPenetrationBC" => noPenetrationBC_revq,
"Rho1E2U3BC" => Rho1E2U3BC_revq,
"FreeStreamBC" => FreeStreamBC_revq,
"ExpBC" => ExpBC_revq,
"isentropicVortexBC" => isentropicVortexBC_revq,
)


"""
  This function uses the options dictionary to populate mesh.bndry_funcs_revq.

  The functors must be of time [`BCType_revq`](@ref)

  If opts["need_adjoint"] is false, the error BC is supplied instead.

  The function must be callable with the signature

  ```
  func(params::ParamType,
        q::AbstractArray{Tsol,1}, q_bar::AbstractArray{Tres, 1}
        aux_vars::AbstractArray{Tres, 1},
        nrm_xy::AbstractArray{Tmsh,1},
        bndryflux_bar::AbstractArray{Tres, 1},
        bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
  ```

  with inputs

   * q: solution vector at a node
   * aux_vars: auxiliary variables for the node
   * coords: x-y coordinates of the node
   * nrm_xy: normal vector of the node
   * bndryflux_bar: seed vector for the flux vector, length `numDofPerNode`

  and outputs

   * q_bar: vector to be updated with back-propigation of `bndryflux_bar`,
            same shape as `q`


  The functor must also have an outer constructor

  ```
    Foo(mesh::AbstractMesh, eqn::EulerData)
  ```
"""
function getBCFunctors_revq(mesh::AbstractMesh, sbp::AbstractOperator, eqn::EulerData, opts)

  for i = 1:mesh.numBC
    key_i = string("BC", i, "_name")
    if opts["need_adjoint"]
      val = opts[key_i]
    else
      val = "errorBC"  # placeholder, because the BC in use might not have
                       # a reverse-mode version
    end
    ctor = BCDict_revq[val]
    mesh.bndry_funcs_revq[i] = ctor(mesh, eqn)
  end # End for i = 1:mesh.numBC

  return nothing
end # End function getBCFunctors_revq


