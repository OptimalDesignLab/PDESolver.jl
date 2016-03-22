export getBCFunctors

include("bc_solvers.jl")


@doc """
### EulerEquationMod.getBCFluxes

  This function calls other functions to calculate the boundary fluxes, passing
  them pieces of the array needed.  This populates eqn.bndryflux.  It also
  calls writeBoundary() to do any requested output.

  This is a mid level function
"""->
# this is a mid level function
function getBCFluxes(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)
  #get all the fluxes for all the boundary conditions and save them in eqn.bndryflux

  #println("mesh.bndry_funcs = ", mesh.bndry_funcs)

  for i=1:mesh.numBC
  #  println("computing flux for boundary condition ", i)
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range = start_index:end_index
    bndry_facenums_i = view(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = view(eqn.bndryflux, :, :, start_index:(end_index - 1))
 
    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
    calcBoundaryFlux(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i, bndryflux_i)
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
"""->
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
  for i=1:mesh.numBoundaryEdges
    el = mesh.bndryfaces[i].element
    face = mesh.bndryfaces[i].face
    for j=1:sbp.numfacenodes
      jb = sbp.facenodes[j, face]
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
"""->
function interpolateBoundary(mesh::AbstractDGMesh, sbp, eqn, opts, q::Abstract3DArray, q_bndry::Abstract3DArray)

  # interpolate solutions
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)

  # calculate aux_vars
  for i=1:mesh.numBoundaryEdges
    for j=1:mesh.sbpface.numnodes
      q_vals = view(eqn.q_bndry, :, j, i)
      eqn.aux_vars_bndry[1, j, i] = calcPressure(eqn.params, q_vals)
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
  sbp : AbstractSBP
  eqn : EulerEquation
  functor : a callable object that calculates the boundary flux at a node
  idx_range: the Range describing which Boundaries have the current BC
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
function calcBoundaryFlux{Tmsh,  Tsol, Tres}( mesh::AbstractCGMesh{Tmsh}, 
                          sbp::AbstractSBP, eqn::EulerData{Tsol}, 
                          functor::BCType, idx_range::UnitRange, 
                          bndry_facenums::AbstractArray{Boundary,1}, 
                          bndryflux::AbstractArray{Tres, 3})
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("enterted calcBoundaryFlux")
  

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
#    println("element = ", bndry_i.element, ", face = ", bndry_i.face)
#    println("interface ", i)
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]

      # get components
      q = view(eqn.q, :, k, bndry_i.element)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      flux_parametric = view(eqn.flux_parametric, :, k, bndry_i.element, :)
      aux_vars = view(eqn.aux_vars, :, k, bndry_i.element)
      x = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("eqn.bndryflux = ", eqn.bndryflux)
      bndryflux_i = view(bndryflux, :, j, i)

      functor(q2, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)

    end

  end


  return nothing
end

# DG version
function calcBoundaryFlux{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh}, 
                          sbp::AbstractSBP, eqn::EulerData{Tsol}, 
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1}, 
                          bndryflux::AbstractArray{Tres, 3})
# calculate the boundary flux for the boundary condition evaluated by the functor

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
#    println("element = ", bndry_i.element, ", face = ", bndry_i.face)
#    println("interface ", i)
    for j = 1:mesh.sbpface.numnodes

      # get components
      q = view(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = view(eqn.aux_vars_bndry, :, j, global_facenum)
      x = view(mesh.coords_bndry, :, j, global_facenum)
      dxidx = view(mesh.dxidx_bndry, :, :, j, global_facenum)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("eqn.bndryflux = ", eqn.bndryflux)
      bndryflux_i = view(bndryflux, :, j, i)

      functor(q2, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
    end
  end

  return nothing
end


@doc """
### EulerEquationMod.isentropicVortexBC <: BCTypes

  This type and the associated call method define a functor to calculate
  the flux using the Roe Solver using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.
"""->
type isentropicVortexBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC, 
              q::AbstractArray{Tsol,1}, 
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})


#  println("entered isentropicOvrtexBC (low level)")

  # getting qg
  qg = params.qg
  calcIsentropicVortex(x, params, qg)
  RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)

  return nothing

end # ends the function isentropicVortexBC

@doc """
### EulerEquationMod.isentropicVortexBC_physical <: BCTypes

  This type and the associated call method define a functor to calculate
  the actual Euler flux  using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.
"""->
type isentropicVortexBC_physical <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC_physical, 
              q::AbstractArray{Tsol,1}, 
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  calcEulerFlux(params, q, aux_vars, [nx, ny], bndryflux)

  return nothing

end # end function isentropicVortexBC_physical


@doc """
### EulerEquationMod.noPenetrationBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid velocity is projected into the wall.
  

  This is a low level functor
"""
type noPenetrationBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::noPenetrationBC, q::AbstractArray{Tsol,1},  aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector


  # calculate normal vector in xy space
  nx = zero(Tmsh)
  ny = zero(Tmsh)
  tngt = Array(Tmsh, 2)  # tangent vector
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac  
  ny *= fac

  Unrm = nx*q[2] + ny*q[3]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end

  #qg = copy(q)

  # calculate normal velocity
  qg[2] -= nx*Unrm
  qg[3] -= ny*Unrm

  # call Roe solver
  #RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)
  nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  v_vals = params.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)
  # this is a problem: q is in conservative variables even if
  # params says we are using entropy variables
  calcEulerFlux(params, v_vals, aux_vars, [nx2, ny2], bndryflux)

  #TODO: make this a unary minus, not a fp multiplication
  for i=1:4
    bndryflux[i] = -bndryflux[i]
  end

return nothing


end


@doc """
### EulerEquationMod.unsteadyVortexBC <: BCTypes

  This type and the associated call method define a functor to calculate
  the flux using the Roe Solver using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.
"""->
type unsteadyVortexBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::unsteadyVortexBC, q::AbstractArray{Tsol,1},                aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)

  # getting qg
  qg = params.qg
  calcUnsteadyVortex(x, params, qg)

  RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)

  return nothing

end # ends the function unsteadyVortex BC





@doc """
### EulerEquationMod.Rho1E2U3BC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid density is 1, energy = 2, and x velocity = 0.5

  This is a low level functor
"""
type Rho1E2U3BC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::Rho1E2U3BC, q::AbstractArray{Tsol,1},  
              aux_vars::AbstractArray{Tres, 1},  
              x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, bndryflux::AbstractArray{Tres, 1}, 
              params::ParamType{2})



  #println("in Rho1E2U3Bc")
  qg = params.qg

  calcRho1Energy2U3(x, params, qg)

  #println("qg = ", qg)
  # call Roe solver
  RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)

return nothing


end

@doc """
### EulerEquationMod.FreeStreamBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state corresponding to the free stream velocity, using rho_free, Ma, aoa, and E_free

  This is a low level functor
"""
type FreeStreamBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::FreeStreamBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  qg = params.qg

  calcFreeStream(x, params, qg)
  RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)
  
  return nothing
end


@doc """
### EulerEquationMod.allOnesBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where all the conservative variables have a value 1.0

  This is a low level functor
"""

type allOnesBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::allOnesBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  qg = zeros(Tsol, 4)
  calcOnes(x, params, qg)

  RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function call


# every time a new boundary condition is created,
# add it to the dictionary
#const isentropicVortexBC_ = isentropicVortexBC()
#const noPenetrationBC_ = noPenetrationBC()
global const BCDict = Dict{ASCIIString, BCType}(
"isentropicVortexBC" => isentropicVortexBC(),
"noPenetrationBC" => noPenetrationBC(),
"Rho1E2U3BC" => Rho1E2U3BC(),
"isentropicVortexBC_physical" => isentropicVortexBC_physical(),
"FreeStreamBC" => FreeStreamBC(),
"allOnesBC" => allOnesBC(),
"unsteadyVortexBC" => unsteadyVortexBC()
)

@doc """
### EulerEquationMod.getBCFunctors

  This function uses the opts dictionary to populate mesh.bndry_funcs with
  the the functors

  This is a high level function.
"""->
# use this function to populate access the needed values in BCDict
function getBCFunctors(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)
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

