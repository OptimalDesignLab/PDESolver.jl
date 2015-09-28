# Checking for derivative residual

push!(LOAD_PATH, "../../SummationByParts/src/")
push!(LOAD_PATH,"../src/simple_mesh/")

using SummationByParts
# using using PDESolver
using SimpleMesh

lengthx = 1
lengthy = 1
nedx = 1
nedy = 1
numDofPerNode = 1
d = 1

m = simpleMesh{Float64}(lengthx,lengthy,nedx,nedy,3,numDofPerNode)
sbp = TriSBP{Float64}(degree=d)
node_coord = zeros(Float64, (2,sbp.numnodes,m.numEl))
jacobian = zeros(Float64, (2,2,sbp.numnodes,m.numEl))
jacmod = zeros(Float64, (sbp.numnodes,m.numEl))
func = zeros(Float64, (1,sbp.numnodes,m.numEl))
actderiv = zeros(func)
flux = zeros(func)
flux2 = zeros(func)
phi = zeros(func)
dxidx = zeros(func)
detadx = zeros(func)

#=
for nnpe = 2:5

d = nnpe - 1 # Degree of the shape function
p = 1
for p = 1:d
	q = d - p

	m = simpleMesh{Float64}(lengthx,lengthy,nedx,nedy,nnpe)
	# println(m.ien)

	sbp = TriSBP{Float64}(degree=d) 
	node_coord = zeros(Float64, (2,sbp.numnodes,m.nel))
	jacobian = zeros(Float64, (2,2,sbp.numnodes,m.nel))
	jacmod = zeros(Float64, (sbp.numnodes,m.nel))
	func = zeros(Float64, (1,sbp.numnodes,m.nel))
	actderiv = zeros(func)
	flux = zeros(func)
	flux2 = zeros(func)
	phi = zeros(func)
	dxidx = zeros(func)
	detadx = zeros(func)

	# Get the X-Y coordinates of all the nodes
	for i = 1:m.nel
	  vtx_coord = m.vtx_loc[:,m.ien[1:3,i]]
	  vtx_coord = vtx_coord'
	  node_coord[:,:,i] = calcnodes(sbp,vtx_coord)
	end

	mappingjacobian!(sbp, node_coord, jacobian, jacmod) # Get the Jocabian for transformation between actual and iso-parametric space

	# println(round(jacobian,1));
	func[1,:,:] =  (node_coord[1,:,:].^p).*(node_coord[2,:,:].^q)
	actderiv[1,:,:] = p*(node_coord[1,:,:].^(p-1)).*(node_coord[2,:,:].^q); # Actual derivative of the function
	
	for i = 1:m.nel
	  for j = 1:sbp.numnodes
	    actderiv[1,j,i] = actderiv[1,j,i]/jacmod[j,i]
	  end
	end 


	di = 1
	# phi = dxidx.*func
	for i = 1:m.nel
	  for j = 1:sbp.numnodes
	    phi[1,j,i] = jacobian[1,1,j,i]*func[1,j,i]
	  end
	end
	differentiate!(sbp, di, phi, flux)

	di = 2
	# phi = detadx.*func
	for i = 1:m.nel
	  for j = 1:sbp.numnodes
	    phi[1,j,i] = jacobian[1,2,j,i]*func[1,j,i]
	  end
	end
	differentiate!(sbp, di, phi, flux2)

	flux += flux2

	residual = flux - actderiv

	# println("Residual = ", '\n', round(residual,1))
	flag = 0
	error = 0
	for i = 1:m.nel
	  for j = 1:sbp.numnodes
	    error = residual[1,j,i]
	    if(abs(error) > 1e-10)
	      println("Test Failed, aborting.")
	      println("Error = ", error)
	      flag = 1
	      break
	    end
	  end
	end
	if flag == 0
	  println("Test Passed for ", d," order polynomial",'\n')
	end
end

end # end nnpe loop =#