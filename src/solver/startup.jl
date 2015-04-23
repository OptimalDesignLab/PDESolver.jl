# startup script for solving an equation

push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
include("../equation/Equation.jl")  # equation types
include("../rk4/rk4.jl")  # timestepping
include("./euler/euler.jl")  # solver functions
include("./euler/ic.jl")  # initial conditions functions


# timestepping parameters
delta_t = 0.5
t_max = 1.00

# create operator
sbp = TriSBP{Float64}()  # create linear sbp operator

# create mesh
dmg_name = ".null"
smb_name = "tri8l.smb"
mesh = PumiMesh2(dmg_name, smb_name, 1; dofpernode=4)  #create linear mesh with 1 dof per node

# create euler equation
eqn = EulerEquation(sbp)
println("eqn.bigQT_xi = \n", eqn.bigQT_xi)
println("eqn.bigQT_eta = \n", eqn.bigQT_eta)
println("sbp.QT_xi' = \n", sbp.Q[:,:,1].')
println("sbp.QT_eta' = \n", sbp.Q[:,:,2].')


# create vectors to hold solution at current, previous timestep
u0 = zeros(mesh.numDof)  # solution at previous timestep
u = zeros(mesh.numDof) # solution at current timestep


# populate u0 with initial condition
ICZero(mesh, sbp, eqn, u0)



function evalEuler(t, x)
# this function is called by time stepping algorithm
# t is the current time
# x is the solution value at the previous timestep
# u = output, the function value at the current timestep
# u is declared outside this function to avoid reallocating memory

u[:] = 0.0  # zero out u before starting
evalVolumeIntegrals(mesh, sbp, eqn, u, x)
evalBoundaryIntegrals(mesh, sbp, eqn, u, x)
applyMassMatrixInverse(mesh, sbp, eqn, u, x)

return u

end  # end evalEuler


# call timestepper
rk4(evalEuler, delta_t, u0, t_max)

