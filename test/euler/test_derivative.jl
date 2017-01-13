# Checking for derivative residual

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))


using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using Utils
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

include(joinpath(Pkg.dir("PDESolver"),"src/solver/euler/output.jl"))  # printing results to files
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

#function runtest(flag::Int)
println("ARGS = ", ARGS)
println("size(ARGS) = ", size(ARGS))
opts = read_input(ARGS[1])
#opts = read_input("input_vals_channel2.jl")

# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
flag = opts["run_type"]

# timestepping parameters
delta_t = opts["delta_t"]   # delta_t: timestep for RK
t_max = opts["t_max"]       # t_max: maximum time for RK
order = opts["order"]       # order of accuracy

# types of the mesh, SBP, Equation objects
if flag == 1 || flag == 8  # normal run
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Float64
  Tres = Float64
elseif flag == 2  # calculate dR/du
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Float64
  Tres = Float64
elseif flag == 3  # calcualte dR/dx using forward mode
  Tmsh = Dual{Float64}
  Tsbp = Float64
  Tsol = Dual{Float64}
  Tres = Dual{Float64}
elseif flag == 4  # use Newton method using finite difference
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Float64
  Tres = Float64
elseif flag == 5  # use complex step dR/du
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Complex128
  Tres = Complex128
elseif flag == 6 || flag == 7  # evaluate residual error and print to paraview
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Complex128
  Tres = Complex128
end

# create SBP object
println("\nConstructing SBP Operator")
sbp = TriSBP{Tsbp}(degree=order)  # create linear sbp operator

dmg_name = opts["dmg_name"]
smb_name = opts["smb_name"]

# create linear mesh with 4 dof per node
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, arg_dict; dofpernode=1)

# println("mesh successfully created")

# function to be tested is x2y2 and is being differentiated using differntiate.sbp
println("mesh coordinates = \n", mesh.coords)

fval = zeros(sbp.numnodes, mesh.numEl)
for i = 1:mesh.numEl
  for j = 1:sbp.numnodes
  	x = mesh.coords[1,j,i]
  	y = mesh.coords[2,j,i]
    fval[j,i] = x + y
  end
end
println("fval = \n", fval)

# Look at Jacobian
# println("Jacobian = \n", mesh.dxidx)

actderiv = zeros(fval) # Actual derivative w.r.t x
for i = 1:mesh.numEl
  for j = 1:sbp.numnodes
  	x = mesh.coords[1,j,i]
  	y = mesh.coords[2,j,i]
    actderiv[j,i] = 1 
  end
end
println("Actual derivative = \n", actderiv)

sbp_deriv = zeros(actderiv)
sbp_deriv2 = zeros(sbp_deriv)

# println("mesh.jac = ", mesh.jac)

#=fval1 = zeros(fval)
fval2 = zeros(fval)

for i = 1:3
	fval1[i] = fval[i]*mesh.dxidx[1,1,i,1]*mesh.jac[i,1]
end
differentiate!(sbp, 1, fval1, sbp_deriv)
for i = 1:3
	fval2[i] = fval[i]*mesh.dxidx[2,1,i,1]*mesh.jac[i,1]
end
differentiate!(sbp, 1, fval1, sbp_deriv)
=#

differentiate!(sbp, 1, fval, sbp_deriv)

for i = 1:sbp.numnodes
	sbp_deriv[i,1] = sbp_deriv[i]*mesh.dxidx[1,1,i,1]*mesh.jac[i,1]
end 
# println(sbp_deriv)

differentiate!(sbp, 2, fval, sbp_deriv2)

for i = 1:sbp.numnodes
	sbp_deriv2[i,1] = sbp_deriv2[i]*mesh.dxidx[2,1,i,1]*mesh.jac[i,1]
end  
sbp_deriv += sbp_deriv2
println(sbp_deriv) 

error = actderiv - sbp_deriv
println("error = \n", error)

#=
push!(LOAD_PATH, "../../SummationByParts/src/")
# push!(LOAD_PATH,"../src/simple_mesh/")

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
=# 
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
