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
res_0_norm = physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (evalResidual,))

# println("eqn.q = \n", eqn.q)
# println("dxidx = \n", mesh.dxidx)
# println("jacobian determinant = \n", mesh.jac)
# println("coordinates = \n", mesh.coords)
# println("eqn.Axi[:,:,1,1] = \n", eqn.Axi[:,:,1,1])

facts("\nCheck flux jacobian for conservative variables") do
  
  @fact eqn.Axi[:,:,1,1] --> roughly([0.0 0.5 0.0 0.0
                                      0.045125 0.0 0.13435028842 0.2
                                      0.0 -0.335875721 0.0 0.0
                                      0.0 0.7378125 0.0 0.0])
  @fact eqn.Axi[:,:,2,1] --> roughly([0.0 0.5 0.0 0.0
                                      0.01128125 0.0 0.06717514421272198 0.2
                                      0.0 -0.167937860531805 0.0 0.0
                                      0.0 0.7378125 0.0 0.0])
  @fact eqn.Axi[:,:,3,1] --> roughly([0.0 0.5 0.0 0.0
                                     -0.03384375 0.26870057685088766 0.0671751442 0.2
                                      0.05640625 -0.16793786053 0.16793786053 0.0
                                     -0.240235109491 0.71525 0.0225625 0.2351130])
  @fact eqn.Aeta[:,:,1,1] --> roughly([0.0 0.0 0.5 0.0
                                       0.0 -0.3358757210636101 0.0 0.0
                                      -0.1805 0.0 -0.5374011537017761 0.2
                                       0.4653138270684989 0.0 0.6475625 -0.47022601])
  @fact eqn.Aeta[:,:,2,1] --> roughly([0.0 0.0 0.5 0.0
                                       0.0 -0.16793786053180498 0.0 0.0
                                      -0.045125 0.0 -0.268700576850888 0.2
                                       0.2440242074689959 0.0 0.71525 -0.23511300474452695])
  @fact eqn.Aeta[:,:,3,1] --> roughly([0.0 0.0 0.5 0.0
                                       0.05640625 -0.16793786053180487 0.16793786053180484 0.0
                                      -0.03384375 -0.06717514421272192 -0.2687005768508878 0.2
                                       0.24023510949074675 0.0225625 0.71525 -0.2351130047445268])
  
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
  dxidx = Tsol[0.5,0.0]
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
  @fact T*Lambda*Tinv --> roughly(eqn.Axi[:,:,1,1])

  # Check Aeta factorization in being computed correctly
  dxidx = Tsol[0.0,0.5]
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
   @fact T*Lambda*Tinv --> roughly(eqn.Aeta[:,:,1,1])
end

facts("Check if element level product of shape function with nodal flux 
jacobian is being computed correctly") do

  ndof = mesh.numDofPerNode
  nnpe = mesh.numNodesPerElement
  endof = nnpe*ndof # degrees of freedom in an element
  Axidxi = zeros(Tsol, endof, endof)
  shapefuncderiv = ones(Tsol, sbp.numnodes, sbp.numnodes, Tdim)
  Axi = ones(Tsol, ndof, ndof, nnpe)
  Aeta = ones(Axi)
  calcAxidxi(Axidxi, shapefuncderiv, Axi, Aeta, ndof, nnpe)

  @fact Axidxi --> roughly(2*ones(Tsol, endof, endof)) 

end
