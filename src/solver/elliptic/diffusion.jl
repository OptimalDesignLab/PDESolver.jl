abstract AbstractDiffn

type DiffnPoly2nd<: AbstractDiffn
end
function call{Tmsh, Tsol}(obj::DiffnPoly2nd,
						  xy::AbstractArray{Tmsh},
						  lambda::AbstractArray{Tsol, 3})
	# @assert(size(lambda, 1) == Tdim)
	# @assert(size(lambda, 2) == Tdim)
	# the 3rd dimension should be dof per node	
	for dof = 1:size(lambda, 3)
		lambda[1, 1, dof] = xy[1]*xy[1] + 1.0
		lambda[1, 2, dof] = xy[1]*xy[2]
		lambda[2, 2, dof] = xy[2]*xy[2] + 1.0
		lambda[2, 1, dof] = lambda[1, 2, dof]
	end
	return nothing
end

type DiffnPoly6th<: AbstractDiffn
end
function call{Tmsh, Tsol}(obj::DiffnPoly6th,
						  xy::AbstractArray{Tmsh},
						  lambda::AbstractArray{Tsol, 3})
	# @assert(size(lambda, 1) == Tdim)
	# @assert(size(lambda, 2) == Tdim)
	# the 3rd dimension should be dof per node	
	for dof = 1:size(lambda, 3)
		lambda[1, 1, dof] = xy[1]^6 + 1.0
		lambda[1, 2, dof] = xy[1]^3*xy[2]^3
		lambda[2, 2, dof] = xy[2]^6 + 1.0
		lambda[2, 1, dof] = lambda[1, 2, dof]
	end
	return nothing
end

type DiffnPoly0th<: AbstractDiffn
end

function call{Tmsh, Tsol}(obj::DiffnPoly0th,
						  xy::AbstractArray{Tmsh},
						  lambda::AbstractArray{Tsol, 3})
	# @assert(size(lambda, 1) == Tdim)
	# @assert(size(lambda, 2) == Tdim)
	# the 3rd dimension should be dof per node	
	for dof = 1:size(lambda, 3)
		lambda[1, 1, dof] = 1.0e1
		lambda[2, 2, dof] = 1.0e1 
		lambda[1, 2, dof] = 0.0 
		lambda[2, 1, dof] = 0.0
	end
	return nothing
end

global const DiffnDict = Dict{ASCIIString, AbstractDiffn}(
	"poly0th" => DiffnPoly0th(),
	"poly2nd" => DiffnPoly2nd(),
	"poly6th" => DiffnPoly6th()
)

function getDiffnFunc(mesh::AbstractMesh,
					  sbp::AbstractSBP,
					  eqn::EllipticData,
					  opts)
	eqn.diffusion_func = DiffnDict[opts["Diffusion"]]
end

function calcDiffn{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
										   sbo::AbstractSBP,
										   eqn::EllipticData{Tsol, Tres, Tdim},
										   diffusion_func::AbstractDiffn,
										   lambda::AbstractArray{Tres, 5} )

	numElems = mesh.numEl
	numNodesPerElement = mesh.numNodesPerElement

	for elem = 1:numElems
		# @bp
		for n = 1:numNodesPerElement
			xy = sview(mesh.coords, :, n, elem)

			diffusion_func(xy, sview(lambda, :, :, :, n, elem))
		end
	end
end

