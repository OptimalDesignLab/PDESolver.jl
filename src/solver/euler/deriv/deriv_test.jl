# Test for computing the derivative using the design variables. This includes
# using FFD and mesh movement. The disgn variables are angle of attack and the
# coordinates of the control points


resize!(ARGS, 1)
ARGS[1] = "../input_vals_vortex.jl"
#ARGS[1] = "../input_vals_airfoil.jl"

push!(LOAD_PATH, joinpath(Pkg.dir("FFD.jl"), "src"))

using FreeFormDeformation
using MeshMovement
using ODLCommonTools
import ODLCommonTools.sview


include("../startup.jl")
# include("differentiateByMetrics.jl")
objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)

dFluxdM = EulerEquationMod.getdFdm(mesh, sbp, eqn, opts)
#=
# Check against complex step
pert = complex(0, 1e-20)
epsilon = imag(pert)
complex_dFluxdM = zeros(dFluxdM)
Tdim = 2
nrm = zeros(Complex{Float64}, Tdim)
# ctr = 0
for i = 1:mesh.numEl
  for j = 1:mesh.numNodesPerElement
    q_vals = sview(eqn.q, :,j,i)
    aux_vars = sview(eqn.aux_vars, :, j, i)
    for k=1:Tdim  # loop over dimensions
      for p=1:Tdim
        nrm[p] = mesh.dxidx[k, p, j, i]
      end
      for p = 1:Tdim
        flux = eqn.params.flux_vals1
        fill!(flux, 0.0)
        nrm[p] += pert
        EulerEquationMod.calcEulerFlux(eqn.params, q_vals, aux_vars, nrm, flux)
        complex_dFluxdM[:,p,k,j,i] = imag(flux[:])/epsilon
        # flux[:] = imag(flux[:])/pert
        # println("complex_dFluxdM = $(complex_dFluxdM[:,p,k,j,i]), dFluxdM = $(real(dFluxdM[:,k,p,j,i]))")
        nrm[p] -= pert
      end # End for p = 1:Tdim
    end # End for k = 1:Tdim
  end   # End for j = 1:mesh.numNodesPerElement
end     # End for i = 1:mesh.numEl


ctr = 0
for i = 1:length(dFluxdM)
  err = norm(dFluxdM[i] - complex_dFluxdM[i], 2)
  if err > 1e-12
    ctr += 1
  end
end

println("error counter = $ctr")
=#

deriv_bndry_funcs = EulerEquationMod.getBCDerivFunctors(mesh, sbp, eqn, opts)
# println(deriv_bndry_funcs)
dBndryFluxdm = EulerEquationMod.getdBndryFluxdm(mesh, sbp, eqn, opts, deriv_bndry_funcs)
# println("dBndryFluxdm = \n", dBndryFluxdm)

# Check agains complex step
pert = complex(0, 1e-20)
ctr = 0
for i=1:mesh.numBC
  functor_i = mesh.bndry_funcs[i]
  if functor_i == EulerEquationMod.noPenetrationBC()
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range = start_index:end_index  # TODO: should this be start_index:(end_index - 1) ?
    bndry_facenums_i = sview(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = sview(eqn.bndryflux, :, :, start_index:(end_index - 1))

    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
    EulerEquationMod.complex_calcBoundaryFluxdm(mesh, sbp, eqn, functor_i,
                                        idx_range, bndry_facenums_i, bndryflux_i)

    # Check against analytical value
    dbndryflux_i = sview(dBndryFluxdm, :, 1, :, start_index:(end_index -  1))
    println(" size of dbndryflux_i = ", size(dbndryflux_i), " length = ", length(dbndryflux_i))
    println("size of bndryflux_i = ", size(bndryflux_i), " length = ", length(dbndryflux_i))

    for j = 1:length(bndryflux_i)
      error = norm(dbndryflux_i[j] - bndryflux_i[j])
      if error > 1e-12
        println("dbndryflux_i[$j] = $(dbndryflux_i[j]) , bndryflux_i[$j] = $(bndryflux_i[j])")
        ctr += 1
      end
    end
    

  end # End if noPenetrationBC
end
println("ctr = $ctr")

#=
# Get the adjoint vector
include("../startup.jl")
objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
println("objective.lift_val = $(objective.lift_val)")
orig_val = copy(real(objective.lift_val))
adjoint_vec = zeros(Complex128, mesh.numDof)
EulerEquationMod.calcAdjoint(mesh, sbp, eqn, opts, objective, adjoint_vec)

# Write VTK files
PdePumiInterface.saveSolutionToMesh(mesh, real(adjoint_vec))
PdePumiInterface.writeVisFiles(mesh, "adjoint_field")

# Initialize FFD and MeshWarping
# geom_faces = opts["BC2"]
# include("initMeshMotion.jl")

# Get the partial derivative of the functional w.r.t aoa
dJdAlpha = objective.dLiftdAlpha
println("dJdAlpha = $(real(dJdAlpha))")

# Check dJdALpha against the complex step method
eqn.params.aoa += opts["epsilon"]*im
println("aoa = $(eqn.params.aoa)")
EulerEquationMod.calcBndryFunctional(mesh, sbp, eqn, opts, objective)
dJdAlpha_comp = imag(objective.lift_val)/opts["epsilon"]
println("dJdAlpha = $dJdAlpha, dJdAlpha_comp = $dJdAlpha_comp")
println("dJdAlpha error = ", norm(dJdAlpha - dJdAlpha_comp,2))


# Get the partial derivative of the residual vector w.r.t aoa
eqn.params.aoa = opts["aoa"]
eqn.params.aoa += opts["epsilon"]*im # Imaginary perturbation
fill!(eqn.res_vec, 0.0)
fill!(eqn.res, 0.0)
res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalResidual)
dRdAlpha = imag(eqn.res_vec)/opts["epsilon"]
#=
f = open("dRdAlpha.dat", "w")
for i = 1:length(eqn.res_vec)
  println(f, dRdAlpha[i])
end
close(f)
=#
#=
# Verify against finite difference
eqn.params.aoa = opts["aoa"]
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)
res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
orig_res_vec = copy(eqn.res_vec)
eqn.params.aoa += 1e-6
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)
res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
dRdAlpha_FD = (eqn.res_vec - orig_res_vec)/1e-6

f = open("dRdAlpha_error.dat", "w")
for i = 1:length(eqn.res_vec)
  println(f, dRdAlpha_FD[i] - dRdAlpha[i])
end
close(f)
println("dRdAlpha error norm = ", norm(dRdAlpha_FD - dRdAlpha, 2))
=#

dLdx_adjoint = dJdAlpha + dot(adjoint_vec, dRdAlpha)

#=
#----- Finite Differencing -----#
# Get the design variable array
x_dv = reshape(ffd_map.cp_xyz, length(ffd_map.cp_xyz))
x_dv = append!([eqn.params.aoa], x_dv)
=#
# Finite difference derivative of Lagrangian wrt aoa, which is x_dv[1]

pert = 1e-6 # FD perturbation
eqn.params.aoa = opts["aoa"] + pert
fill!(eqn.q, 0.0)
fill!(eqn.q_vec, 0.0)
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)

# Rerun with the perturbed value
ICfunc_name = opts["IC_name"]
ICfunc = EulerEquationMod.ICDict[ICfunc_name]
ICfunc(mesh, sbp, eqn, opts, eqn.q_vec)
pmesh = mesh
EulerEquationMod.init(mesh, sbp, eqn, opts, pmesh)
call_nlsolver(mesh, sbp, eqn, opts, pmesh)

EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
dLdx = (real(objective.lift_val) - orig_val)/pert
println("orig_val = $orig_val, new_val = $(real(objective.lift_val)), dLdx = $dLdx")
println("dLdx_adjoint = $dLdx_adjoint")

errfd_norm = norm(dLdx - dLdx_adjoint,2)
println("err dLdx = $errfd_norm")
println("dJdAlpha = $dJdAlpha, dJdAlpha_comp = $dJdAlpha_comp")
println("dJdAlpha error = ", norm(dJdAlpha - dJdAlpha_comp,2))
=#
