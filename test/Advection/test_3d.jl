push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

include(joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

function clean_dict(collection)
  for i in keys(collection)
    delete!(collection, i)
  end
end

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")
# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_3d.jl"
include(STARTUP_PATH)

facts("----- Testing 3D -----") do

  facts("----- Testing SummationByParts -----") do
    q = ones(1, mesh.numNodesPerElement, 2)
    res = zeros(q)
    weakdifferentiate!(sbp, 1, q, res, trans=false)

    for i=1:2
      for j=1:mesh.numNodesPerElement
        @fact res[1, j, i] --> roughly(0.0, atol=1e-13)
      end
    end

  end
  # test constant IC -> zero residual
  opts["use_src_term"] = false
  opts["BC1_name"] = "constantBC"
  fill!(eqn.res, 0.0)
  fill!(eqn.q_vec, 2.0)
  fill!(eqn.q, 2.0)
  AdvectionEquationMod.ICConstant(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  AdvectionEquationMod.getBCFunctors(mesh, sbp, eqn, opts)
  evalAdvection(mesh, sbp, eqn, opts)

  for i=1:length(eqn.res)
    @fact eqn.res[i] --> roughly(0.0, atol=1e-13)
  end
  println("size(eqn.res) = ", size(eqn.res))
  writedlm("resout.dat", eqn.res)


function test_bc(flux_exp)
# test summing boundary condition
 fill!(eqn.bndryflux, 0.0)
 AdvectionEquationMod.evalBndry(mesh, sbp, eqn)

 flux_in = zero(eltype(eqn.bndryflux))
 flux_out = zero(flux_in)
 for i=1:mesh.numBoundaryFaces
   net_flux_face = zero(flux_in)
   for j=1:mesh.numNodesPerFace
     net_flux_face += mesh.sbpface.wface[j]*eqn.bndryflux[1, j, i]
   end
   if net_flux_face > 0
     flux_in += net_flux_face
   else
     flux_out += net_flux_face
   end
 end

 # flux = area of face * solution * alpha
 @fact abs(flux_in) --> roughly(flux_exp, atol=1e-13)
 @fact abs(flux_out) --> roughly(abs(flux_in), atol=1e-13)
end



  facts("----- Testing BCSolver -----") do
    # check that the solver produces the regular flux when qL = qR
    q2 = rand(1, mesh.numNodesPerElement, mesh.numEl)
    q1 = eqn.q
    eqn.q = q2
    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalSCResidual(mesh, sbp, eqn)
    nrm = [1., 1, 1]  # arbirary normal vector
    alphas_xy = [eqn.params.alpha_x, eqn.params.alpha_y, eqn.params.alpha_z]
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        u = eqn.q[1, j, i]
        dxidx = sview(mesh.dxidx, :, :, j, i)

        # calculate the flux directly
        alphas_parametric = dxidx*alphas_xy
        flux_parametric = alphas_parametric*u
        net_flux = sum(flux_parametric.*nrm)

        # calculate boundary flux
        bndryflux_calc = AdvectionEquationMod.RoeSolver(u, u, eqn.params, nrm, dxidx)
        # calculate flux from evalSCResidual
        bndryflux_weak = zero(eltype(eqn.flux_parametric))
        for d=1:3
           bndryflux_weak += nrm[d]*eqn.flux_parametric[1, j, i, d]
        end

        for d=1:3
          @fact eqn.flux_parametric[1, j, i, d] --> roughly(flux_parametric[d], atol=1e-13)
        end

         @fact bndryflux_calc --> roughly(net_flux, atol=1e-13)
       end
     end

   end

   facts("----- Testing boundary flux calculation -----") do
     # set alpha_x = 1. all others zero, q = constant, check flux
     fill!(eqn.q, 2.0)
     fill!(eqn.res, 0.0)
     eqn.params.alpha_x = 1.0
     eqn.params.alpha_y = 0.0
     eqn.params.alpha_z = 0.0

     opts["BC1_name"] = "constantBC"
     opts["numBC"] = 1
     AdvectionEquationMod.getBCFunctors(mesh, sbp, eqn, opts)

     test_bc(8.0)
     eqn.params.alpha_x = 0.0
     eqn.params.alpha_y = 1.0
     eqn.params.alpha_z = 0.0

     test_bc(8.0)
     eqn.params.alpha_x = 0.0
     eqn.params.alpha_y = 0.0
     eqn.params.alpha_z = 1.0

     test_bc(8.0)
     eqn.params.alpha_x = 1.0
     eqn.params.alpha_y = 1.0
     eqn.params.alpha_z = 1.0
     test_bc(3*8.0)


   end

end
