# startup script for solving an equation

push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
include("../equation/Equation.jl")  # equation types
include("../rk4/rk4.jl")  # timestepping
include("./euler/euler.jl")  # solver functions
include("./euler/ic.jl")  # initial conditions functions
include("./euler/adaptfuncs.jl")  # mesh adapt function
include("./euler/output.jl")  # print results to file
#=
# timestepping parameters
delta_t = 0.5
t_max = 1.00
=#

# create operator
sbp = TriSBP{Float64}()  # create linear sbp operator


# create mesh
dmg_name = ".null"
#smb_name = "tri2l.smb"
#smb_name = "tri18l.smb"
smb_name = "adapt_big.smb"
mesh = PumiMesh2(dmg_name, smb_name, 1; dofpernode=4)  #create linear mesh with 1 dof per node

# create euler equation
eqn = EulerEquation(sbp)

# create vectors to hold solution at current, previous timestep
u0 = zeros(mesh.numDof)  # solution at previous timestep
u = zeros(mesh.numDof) # solution at current timestep



# populate u0 with initial condition
ICsmoothHeaviside(mesh, sbp, eqn, u0)
saveSolutionToMesh(mesh, u0)

# create C-callible function
cfunc3 = cfunction(shockRefine2, Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}, Ptr{Void}, Ptr{TriSBP{Float64}}))

# pass the function to C
createAnisoFunc(mesh.m_ptr, cfunc3, mesh.f_ptr, sbp)

writeVtkFiles("output_pre", mesh.m_ptr)

# run mesh adaptation
runAnisoAdapt(mesh.m_ptr)

writeVtkFiles("output_post", mesh.m_ptr)

# reinitialize mesh 
reinitPumiMesh2(mesh)
u =zeros(mesh.numDof)
retrieveSolutionFromMesh(mesh, u);

printSolution(mesh, u)
printCoordinates(mesh)

#println("u = \n", u)

