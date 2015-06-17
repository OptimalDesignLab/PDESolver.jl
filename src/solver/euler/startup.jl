# startup script for solving an equation

 push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
 push!(LOAD_PATH, "../../equation")
 push!(LOAD_PATH, "../../common")
 push!(LOAD_PATH, "/users/creanj/.julia/v0.4/PDESolver/src/solver/euler")
#push!(LOAD_PATH, "../../../../PUMI")
using PDESolverCommon
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
include("../../rk4/rk4.jl")  # timestepping
include("./output.jl")  # printing results to files


function getResType(Tmsh::DataType, Tsbp::DataType, Tsol::DataType )
# figure out what type eqn.res needs to be, taking into account
# algorithmic differentiation of the different arguments
# to support reverse mode, will need to consider the type of Tres as an input
  if Tsol <: DualNumbers.Dual # differentiating wrt eqn.q
    Tres = Tsol
  elseif  Tmsh <: DualNumbers.Dual  # differentiating wrt mesh.coords
    Tres = Tmsh
    Tsol = Tmsh  # derivative of coordinates will end up in eqn.q
  else  # no algorithmic differntiation
    Tres = Tsol
  end

  return Tres

end




function runtest(flag::Int)
# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
# timestepping parameters
delta_t = 0.005
#t_max = 0.025
t_max = 5.00
#t_max = 1.0
order = 2  # order of accuracy




# types of the mesh, SBP, Equation objects
if flag == 1  # normal run
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Float64
  Tres = Float64
elseif flag == 2  # calculate dR/du
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Dual{Float64}
  Tres = Dual{Float64}
elseif flag == 3  # calcualte dR/dx using forward mode
  Tmsh = Dual{Float64}
  Tsbp = Float64
  Tsol = Dual{Float64}
  Tres = Dual{Float64}
end


# create operator
sbp = TriSBP{Tsbp}(degree=order)  # create linear sbp operator

# create mesh
dmg_name = ".null"
smb_name = "../../mesh_files/quarter_vortex3q.smb"
#smb_name = "../../mesh_files/quarter_vortex8l.smb"
#smb_name = "../../mesh_files/tri30l.smb"
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp; dofpernode=4)  #create linear mesh with 4 dof per node





# create euler equation
eqn = EulerEquation1{Tsol, Tres}(mesh, sbp)
#eqn = EulerEquation{Tsol}(mesh, sbp, Float64)


#SL0 = zeros(Tsol, mesh.numDof)  # solution at previous timestep
#SL = zeros(Tsol, mesh.numDof) # solution at current timestep
SL = eqn.SL
SL0 = eqn.SL0

# populate u0 with initial condition
# ICZero(mesh, sbp, eqn, SL0)
# ICLinear(mesh, sbp, eqn, SL0)
 ICIsentropicVortex(mesh, sbp, eqn, SL0)
#ICRho1E2(mesh, sbp, eqn, SL0)
#ICRho1E2U3(mesh, sbp, eqn, SL0)

#ICVortex(mesh, sbp, eqn, SL0)
#ICIsentropicVortex(mesh, sbp, eqn, SL0)

SL_exact = deepcopy(SL0)


saveSolutionToMesh(mesh, SL0)

writeVtkFiles("solution_ic",mesh.m_ptr)



# call timestepper
if flag == 1 # normal run
 rk4(evalEuler, delta_t, t_max, mesh, sbp, eqn)
# println("rk4 @time printed above")
elseif flag == 2 # forward diff dR/du

  # define nested function
  function dRdu_rk4_wrapper(u_vals::AbstractVector, SL::AbstractVector)
    eqn.SL0 = u_vals
    eqn.SL0 = SL
    rk4(evalEuler, delta_t, t_max, mesh, sbp, eqn)
    return nothing
  end

  # use ForwardDiff package to generate function that calculate jacobian
  calcdRdu! = forwarddiff_jacobian!(drDu_rk4_wrapper, Float64, fadtype=:dual, n = mesh.numDof, m = mesh.numDof)

  jac = zeros(Float64, mesh.numDof, mesh.numDof)  # array to be populated
  calcdRdu!(eqn.SL0, jac)

elseif flag == 3 # calculate dRdx

  # dRdx here

end







if flag == 1


    SL_diff = SL - SL_exact
    SL_norm = norm(SL_diff)/mesh.numDof
    #SL_side_by_side = [SL_exact  SL]

    #=
    println("\n\n\n")
    println("SL_diff: \n")
    for i=1:size(SL_diff)[1]
      println(SL_diff[i,:])
    end
    println("SL_side_by_side: \n")
    for i=1:size(SL_side_by_side)[1]
      println(i, " ", SL_side_by_side[i,:])
    end
    =#
    println("SL_norm: \n",SL_norm,"\n")


    saveSolutionToMesh(mesh, SL0)
    printSolution(mesh, SL0)
    printCoordinates(mesh)
    writeVtkFiles("solution_done",mesh.m_ptr)
end
end


runtest(1)
