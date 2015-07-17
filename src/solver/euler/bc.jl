export getBCFunctors

include("bc_solvers.jl")

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


      functor(q, aux_vars, x, dxidx, nrm, bndryflux_i, eqn)
    end
  end


  return nothing
end



type isentropicVortexBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC, q::AbstractArray{Tsol,1}, aux_vars::AbstractArray{Tsol, 1},  x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, bndryflux::AbstractArray{Tres, 1}, eqn::EulerData{Tsol, 2})


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)

  # getting qg
  qg = zeros(Tsol, 4)
  calcIsentropicVortex(x, eqn, qg)
#  calcVortex(x, eqn, qg)

  RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, eqn)

  return nothing

end # ends the function isentropicVortex BC



type noPenetrationBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::noPenetrationBC, q::AbstractArray{Tsol,1}, aux_vars::AbstractArray{Tsol, 1},  x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, bndryflux::AbstractArray{Tres, 1}, eqn::EulerData{Tsol, 2})
# a clever optimizing compiler will clean this up


# calculate normal vector in xy space
xy_nrm = Array(Tmsh, 2)  # normal vector
tngt = Array(Tmsh, 2)  # tangent vector
xy_nrm[1] = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
xy_nrm[2] = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

tngt[1] =  1.0 - xy_nrm[1]
tngt[2] =  1.0 - xy_nrm[2]

qg = copy(q)

# calculate normal velocity
qg[2] *= tngt[1]
qg[3] *= tngt[2]


# call Roe solver
RoeSolver(q, qg, aux_vars, dxidx, nrm, bndryflux, eqn)

return nothing


end






# every time a new boundary condition is created,
# add it to the dictionary
const isentropicVortexBC_ = isentropicVortexBC()
const noPenetrationBC_ = noPenetrationBC()
global const BCDict = Dict{ASCIIString, BCType} (
"isentropicVortexBC" => isentropicVortexBC(),
"noPenetrationBC" => noPenetrationBC()
)

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

