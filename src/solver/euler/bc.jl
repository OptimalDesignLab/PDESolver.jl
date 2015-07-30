export getBCFunctors

include("bc_solvers.jl")

@doc """
### EulerEquationMod.calcBoundaryFlux

  This function calculates the boundary flux for the portion of the boundary
  with a particular boundary condition.

  Inputs:
  mesh : AbstractMesh
  sbp : SBPOperator
  eqn : EulerEquation
  functor : a callable object that calculates the boundary flux at a node
  bndry_facenums:  An array with elements of type Boundary that tell which
                   element faces have the boundary condition
  Outputs:
  bndryflux : the array to store the boundary flux, corresponds to 
              bndry_facenums

  The functor must have the signature
  functor( q, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)

  where all arguments (except params) are vectors of values at a node.

  This is a mid level function.
"""->
# mid level function
function calcBoundaryFlux{Tmsh, Tsbp, Tsol, Tres}( mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, functor::BCType, bndry_facenums::AbstractArray{Boundary,1}, bndryflux::AbstractArray{Tres, 3})
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("enterted calcBoundaryFlux")
  

  nfaces = length(bndry_facenums)
#=  
  for i=1:length(bndry_facenums)
#    println("i = ", i)
    println("element ", bndry_facenums[i].element, " edge ", bndry_facenums[i].face)
  end
=#

  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
#    println("element = ", bndry_i.element, ", face = ", bndry_i.face)
#    println("interface ", i)
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]

      # get components
      q = view(eqn.q, :, k, bndry_i.element)
      aux_vars = view(eqn.aux_vars, :, k, bndry_i.element)
      x = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("eqn.bndryflux = ", eqn.bndryflux)
      bndryflux_i = view(bndryflux, :, j, i)


      functor(q, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
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
function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC, q::AbstractArray{Tsol,1}, aux_vars::AbstractArray{Tsol, 1},  x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)

  # getting qg
  qg = zeros(Tsol, 4)
  calcIsentropicVortex(x, params, qg)
#  calcVortex(x, eqn, qg)

  RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)

  return nothing

end # ends the function isentropicVortex BC


@doc """
### EulerEquationMod.noPenetrationBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid velocity is projected into the wall.

  This is a low level functor
"""
type noPenetrationBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::noPenetrationBC, q::AbstractArray{Tsol,1}, aux_vars::AbstractArray{Tsol, 1},  x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector

# calculate normal vector in xy space
nx = zero(Tmsh)
ny = zero(Tmsh)
#xy_nrm = Array(Tmsh, 2)  # normal vector
tngt = Array(Tmsh, 2)  # tangent vector
nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
fac = 1.0/(sqrt(nx*nx + ny*ny))
# normalize normal vector
nx *= fac  
ny *= fac

Unrm = nx*q[2] + ny*q[3]


#tngt[1] =  1.0 - xy_nrm[1]
#tngt[2] =  1.0 - xy_nrm[2]

qg = copy(q)

# calculate normal velocity
qg[2] -= nx*Unrm
qg[3] -= ny*Unrm

#println("q = ", q)
#println("qg = ", qg)
#diff = q -qg
#println("diff = ", diff)
# call Roe solver
#RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)
nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

#println("q = ", q)
#println("qg = ", qg)
#println("nscl = ", [nx, ny])
#println("normal vector = ", [nx2, ny2])
calcEulerFlux(params, qg, aux_vars, [nx2, ny2], bndryflux)

for i=1:4
  bndryflux[i] *= -1
end

#bndryflux *= -1  # make it negative because all fluxes are negative
#println("bndryflux = ", bndryflux)
#println("bndryflux = ", bndryflux)
#print("\n")

return nothing


end

@doc """
### EulerEquationMod.Rho1E2U3BC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid density is 1, energy = 2, and x velocity = 0.5

  This is a low level functor
"""
type Rho1E2U3BC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::Rho1E2U3BC, q::AbstractArray{Tsol,1}, aux_vars::AbstractArray{Tsol, 1},  x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})



#println("in Rho1E2U3Bc")
qg = zeros(Tsol, 4)

calcRho1Energy2U3(x, params, qg)
#println("qg = ", qg)
# call Roe solver
RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, params)

return nothing


end






# every time a new boundary condition is created,
# add it to the dictionary
#const isentropicVortexBC_ = isentropicVortexBC()
#const noPenetrationBC_ = noPenetrationBC()
global const BCDict = Dict{ASCIIString, BCType} (
"isentropicVortexBC" => isentropicVortexBC(),
"noPenetrationBC" => noPenetrationBC(),
"Rho1E2U3BC" => Rho1E2U3BC()
)

@doc """
### EulerEquationMod.getBCFunctors

  This function uses the opts dictionary to populatemesh.bndry_funcs with
  the the functors

  This is a high level function.
"""->
# use this function to populate access the needed values in BCDict
function getBCFunctors(mesh::PumiMesh, sbp::SBPOperator, eqn::EulerData, opts)
# populate the array mesh.bndry_funcs with the functors for the boundary condition types

println("Entered getBCFunctors")
println("BCDict = ", BCDict)

for i=1:mesh.numBC
  key_i = string("BC", i, "_name")
  val = opts[key_i]
  println("BCDict[val] = ", BCDict[val])
  mesh.bndry_funcs[i] = BCDict[val]
end

return nothing

end

