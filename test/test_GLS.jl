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
include("../src/solver/euler/GLS.jl")

fill!(eqn.res, 0.0)
fill!(eqn.res_vec,0.0)
# res_0 = zeros(eqn.res_vec)
res_0_norm = calcResidual(mesh, sbp, eqn, opts, evalEuler)

# println("eqn.q = \n", eqn.q)
# println("dxidx = \n", mesh.dxidx)
# println("jacobian determinant = \n", mesh.jac)
# println("coordinates = \n", mesh.coords)
# println("eqn.Axi[:,:,1,1] = \n", eqn.Axi[:,:,1,1])

facts("\nCheck flux jacobian for conservative variables") do
  
  @fact eqn.Axi[:,:,1,1] --> roughly([0.0 + 0.0im 0.5 + 0.0im 0.0 + 0.0im 0.0 + 0.0im
                 0.04512499999999999 + 0.0im 0.0 + 0.0im 0.13435028842544403 - 0.0im 0.19999999999999996 + 0.0im
                 0.0 + 0.0im -0.3358757210636101 + 0.0im 0.0 + 0.0im 0.0 + 0.0im
                 0.0 + 0.0im 0.7378125000000001 + 0.0im 0.0 + 0.0im 0.0 + 0.0im])
  @fact eqn.Axi[:,:,2,1] --> roughly([0.0 + 0.0im 0.5 + 0.0im -8.326672684688674e-17 + 0.0im 0.0 + 0.0im
                 0.011281249999999996 + 0.0im 2.796727192030474e-17 + 0.0im 0.06717514421272198 - 0.0im 0.19999999999999996 + 0.0im
                 7.514822097931527e-18 + 0.0im -0.167937860531805 + 0.0im 4.4747635072487585e-17 + 0.0im -3.330669073875469e-17 + 0.0im
                 -4.0638194054697806e-17 + 0.0im 0.7378125000000003 + 0.0im -1.1911305275447156e-16 + 0.0im 3.9154180688426634e-17 + 0.0im])
  @fact eqn.Axi[:,:,3,1] --> roughly([0.0 + 0.0im 0.5 + 0.0im 0.0 + 0.0im 0.0 + 0.0im
                 -0.033843749999999895 + 0.0im 0.26870057685088766 + 0.0im 0.0671751442127219 - 0.0im 0.19999999999999996 + 0.0im
                 0.056406249999999825 + 0.0im -0.1679378605318048 + 0.0im 0.16793786053180476 + 0.0im 0.0 + 0.0im
                 -0.24023510949074664 + 0.0im 0.7152499999999997 + 0.0im 0.022562499999999926 + 0.0im 0.23511300474452665 + 0.0im])
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


fill!(eqn.res, 0.0)
fill!(eqn.res_vec,0.0)

parametricFluxJacobian(mesh, sbp, eqn)

GLS(mesh, sbp, eqn)
facts("Check if GLS residual is getting added to the residual") do
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      for k = 1:mesh.numDofPerNode
        @fact abs(eqn.res[k,j,i]) --> greater_than(0.0) "GLS values are not getting added to eqn.res at dof $k, node $j, element $i."
      end
    end
  end
end


# println("Ax = \n", Ax)
# println("\nAy = \n", Ay)

facts("Check if the Flux Jacobian in the X & Y direction are being computed accrately") do
  q = ones(Tsol, 4)
  Ax = zeros(Tsol, 4, 4)
  Ay = zeros(Tsol, 4, 4)
  calcFluxJacobian(eqn, q, Ax, Ay)

  @fact Ax --> roughly([0.0 1.0 0.0  0.0
                       -0.6 1.6 -0.4 0.4
                       -1.0 1.0 1.0 0.0
                       -0.6 0.6 -0.4 1.4])
  @fact Ay --> roughly([0.0 0.0 1.0 0.0
                       -1.0 1.0 1.0 0.0
                       -0.6 -0.4 1.6 0.4
                       -0.6 -0.4 0.6 1.4])
end

facts("Check if eigen value factorization is being computed correctly") do
  q = Complex{Float64}[2.0,0.0,-1.3435028842544403,2.236964285714286]
  T = zeros(Tsol,4,4)
  Tinv = zeros(T)
  Lambda = eye(T)
  
  # check Axi factorization is being computed correctly
  dxidx = Complex{Float64}[0.5,0.0]
  calcEigenFactorization(eqn, q, dxidx, T, Tinv, Lambda)
  @fact T --> roughly([1.0 0.0 1.0 1.0
                       0.0 0.0 0.7071067811865476 -0.7071067811865476
                      -0.67175144212722 -1.0 -0.6717514421272202 -0.6717514421
                       0.22562500000000002 0.6717514 1.4756250 1.4756250])
  @fact Tinv --> roughly([0.8195 0.0 -0.5374011537017759 -0.8
                         -0.6717514421272202 0.0 -1.0 0.0
                          0.09025 0.7071067811865475 0.26870057685088794 0.4
                          0.09025 -0.7071067811865475 0.26870057685088794 0.4])
  @fact Lambda --> roughly([0.0 0.0 0.0 0.0
                            0.0 0.0 0.0 0.0
                            0.0 0.0 0.3535533905932738 0.0
                            0.0 0.0 0.0 -0.3535533905932738])

  # Check Aeta factorization in being computed correctly
  dxidx = Complex{Float64}[0.0,0.5]
  calcEigenFactorization(eqn, q, dxidx, T, Tinv, Lambda)
  @fact T --> roughly([1.0 0.0 1.0 1.0
                       0.0 1.0 0.0 0.0
                      -0.6717514421272202 -0.0 0.035355339 -1.3788582233
                       0.225625 0.0 1.0006250000000003 1.950625])
  @fact Tinv --> roughly([0.8195 0.0 -0.5374011537017759 -0.8
                         -0.0 1.0 -0.0 0.0
                          0.56525 0.0 0.9758073580374353 0.4
                         -0.38475 -0.0 -0.43840620433565947 0.4])
  @fact Lambda --> roughly([-0.3358757210636101 0.0 0.0 0.0
                             0.0 -0.3358757210636101 0.0 0.0
                             0.0 0.0 0.01767766952966371 0.0
                             0.0 0.0 0.0 -0.6894291116568838])
end