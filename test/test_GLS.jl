# Test GLS
using PDESolver
using FactCheck
using ODLCommonTools
using ArrayViews

include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_vortex_GLS.jl"
opts = read_input(ARGS[1])
facts("Check input arguments applied correctly") do
  @fact opts["use_GLS"] --> true
  @fact opts["Q_transpose"] --> true
end


include("../src/solver/euler/startup.jl")
include("../src/solver/euler/SUPG.jl")

fill!(eqn.res, 0.0)
fill!(eqn.res_vec,0.0)
res_0 = zeros(eqn.res_vec)
res_0_norm = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_0)

facts("\nCheck flux jacobian for conservative variables") do
  @fact eqn.Axi[:,:,1,1] --> roughly([0.0 + 0.0im 0.4999999999999998 + 0.0im 2.498001805406602e-16 + 0.0im 0.0 + 0.0im
                 0.04512499999999997 + 0.0im -1.6780363152182847e-16 + 0.0im 0.13435028842544397 + 0.0im 0.19999999999999987 + 0.0im
                 -9.017786517517833e-17 + 0.0im -0.3358757210636099 + 0.0im -2.684858104349255e-16 + 0.0im 9.992007221626407e-17 + 0.0im
                 2.3247095601955315e-16 + 0.0im 0.7378124999999998 + 0.0im 3.235224588227226e-16 + 0.0im -2.3492508413055986e-16 + 0.0im])
  @fact eqn.Axi[:,:,2,1] --> roughly([0.0 + 0.0im 0.5 + 0.0im -8.326672684688674e-17 + 0.0im 0.0 + 0.0im
                 0.011281249999999996 + 0.0im 2.796727192030474e-17 + 0.0im 0.06717514421272198 - 0.0im 0.19999999999999996 + 0.0im
                 7.514822097931527e-18 + 0.0im -0.167937860531805 + 0.0im 4.4747635072487585e-17 + 0.0im -3.330669073875469e-17 + 0.0im
                 -4.0638194054697806e-17 + 0.0im 0.7378125000000003 + 0.0im -1.1911305275447156e-16 + 0.0im 3.9154180688426634e-17 + 0.0im])
  @fact eqn.Axi[:,:,3,1] --> roughly([0.0 + 0.0im 0.5 + 0.0im -1.6653345369377348e-16 + 0.0im 0.0 + 0.0im
                 -0.033843749999999916 + 0.0im 0.2687005768508877 + 0.0im 0.06717514421272185 - 0.0im 0.19999999999999996 + 0.0im
                 0.05640624999999984 + 0.0im -0.16793786053180476 + 0.0im 0.16793786053180484 + 0.0im -6.661338147750938e-17 + 0.0im
                 -0.24023510949074672 + 0.0im 0.7152499999999997 + 0.0im 0.022562499999999687 + 0.0im 0.23511300474452673 + 0.0im])
  @fact eqn.Aeta[:,:,1,1] --> roughly([0.0 + 0.0im 8.326672684688674e-17 + 0.0im 0.5 + 0.0im 0.0 + 0.0im
                 7.514822097931527e-18 + 0.0im -0.3358757210636101 + 0.0im 2.2373817536243792e-17 + 0.0im 3.330669073875469e-17 + 0.0im
                 -0.1805 + 0.0im -5.593454384060949e-17 + 0.0im -0.5374011537017761 + 0.0im 0.19999999999999996 + 0.0im
                 0.4653138270684989 + 0.0im 1.2287046380343727e-16 + 0.0im 0.6475625 + 0.0im -0.47022600948905413 + 0.0im])
  @fact eqn.Aeta[:,:,2,1] --> roughly([0.0 + 0.0im 4.163336342344337e-17 + 0.0im 0.4999999999999999 + 0.0im 0.0 + 0.0im
                 9.393527622414408e-19 + 0.0im -0.16793786053180498 + 0.0im 5.5934543840609466e-18 + 0.0im 1.6653345369377344e-17 + 0.0im
                 -0.04512499999999998 + 0.0im -1.398363596015237e-17 + 0.0im -0.268700576850888 + 0.0im 0.1999999999999999 + 0.0im
                 0.2440242074689959 + 0.0im 6.143523190171865e-17 + 0.0im 0.7152500000000003 + 0.0im -0.23511300474452695 + 0.0im])
  @fact eqn.Aeta[:,:,3,1] --> roughly([0.0 + 0.0im 0.0 + 0.0im 0.5000000000000002 + 0.0im 0.0 + 0.0im
                 0.05640624999999985 + 0.0im -0.16793786053180487 + 0.0im 0.16793786053180484 + 0.0im 0.0 + 0.0im
                 -0.03384374999999992 + 0.0im -0.06717514421272192 - 0.0im -0.2687005768508878 + 0.0im 0.20000000000000004 + 0.0im
                 0.24023510949074675 + 0.0im 0.022562499999999937 + 0.0im 0.71525 + 0.0im -0.2351130047445268 + 0.0im])
  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
      Flux_xi = eqn.Axi[:,:,j,i]*eqn.q[:,j,i]
      Flux_eta = eqn.Aeta[:,:,j,i]*eqn.q[:,j,i]
      @fact Flux_xi --> roughly(eqn.flux_parametric[:,j,i,1]) "Xi mismatch at Node $j element $i"
      @fact Flux_eta --> roughly(eqn.flux_parametric[:,j,i,2]) "Eta mismatch at Node $j element $i"
  	end
  end
end

facts("Check residual vector without any stabilization") do
  @fact eqn.res_vec --> roughly([0.026726907259624955,0.48767456657840463,
  	                             0.5783855773525426,-0.23352036974413096,
  	                             -0.04667697365332135,-0.5286862525347082,
  	                             0.056759263290692086,-0.16958045972628238,
  	                             -0.23493789531609155,-0.018798326722505865,
  	                             -0.45002903493239477,0.026981780972407405])
end

fill!(eqn.res, 0.0)
fill!(eqn.res_vec,0.0)

FluxJacobian(mesh, sbp, eqn)
#println("Node1 = \n",eqn.Aeta[:,:,1,1])
#println("Node2 = \n",eqn.Aeta[:,:,2,1])
#println("Node3 = \n",eqn.Aeta[:,:,3,1])


# println(eqn.Aeta)



GLS(mesh, sbp, eqn)
facts("Check if GLS residual is getting added to the residual") do
  @fact eqn.res[:,1,1] --> roughly([0.111689561162,0.507611237471,0.516680308782,-0.1081473535781])
  @fact eqn.res[:,2,1] --> roughly([0.038285680249,-0.50874958164,-0.00494600527,-0.0442074435602])
  @fact eqn.res[:,3,1] --> roughly([-0.149975241412,0.00113834417,-0.51173430350,0.15235479713840])
end

