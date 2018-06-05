# test Utils module
mutable struct TestParams{Tdim} <: AbstractParamType{Tdim}
  time::Timings
end

mutable struct TestData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
  params::TestParams{2}
  M::Array{Float64, 1}
  Minv::Array{Float64, 1}
  comm::MPI.Comm
end

mutable struct TestMesh{Tmsh} <: AbstractMesh{Tmsh}
  jac::AbstractArray{Float64, 2}
  dofs::AbstractArray{Int, 3}
  comm::MPI.Comm
  commsize::Int
  numNodes::Int
  numNodesPerElement::Int
  numEl::Int
  numDofPerNode::Int
  min_node_dist::Float64
  dim::Int
end

mutable struct FakeSBP
end

"""
  Test BufferedIO type
"""
function test_utils_io()
  @testset "---- Testing IO -----" begin
    fname = "iotest.dat"
    f = open(fname, "w")
    f2 = BufferedIO(f)

    println(f2, 1)
    println(f2, 2)

    fstat = stat(fname)
    @test ( fstat.size )== 0
    flush(f2)
  #  close(f)

    data = readdlm(fname)
    @test ( data[1] )== 1
    @test ( data[2] )== 2

    close(f)
  end

  return nothing
end

#test_utils_io()
add_func1!(EulerTests, test_utils_io, [TAG_UTILS, TAG_SHORTTEST])

"""
  Test Timings type
"""
function test_utils_timing()
  @testset "----- Testing Timings -----" begin
    t = Timings()
  #=
    t.t_volume = 1.0
    t.t_face = 2.0
    t.t_source = 3.0
    t.t_sharedface = 4.0
    t.t_bndry = 5.0
    t.t_send = 6.0
    t.t_wait = 7.0
    t.t_allreduce = 8.0
    t.t_jacobian = 9.0
    t.t_solve = 10.0
    t.t_barrier = 11.0
    t.t_barrier2 = 12.0
    t.t_barrier3 = 13.0
    t.t_barriers[1:7] = 14.0:20.0

    write_timings(t, "timing_breakdown")
    vals = readdlm("timing_breakdown.dat")
    for i=1:20
      @test isapprox( vals[i], Float64(i)) 
    end
  =#

    function mysleep(t)
      sleep(t)
      return 42
    end

    @time_all mysleep(0.0001)
    gc()
    val, telapsed, gc_time, alloced = @time_all mysleep(0.1)
    @test ( val )== 42
    @test isapprox( telapsed, 0.1) atol=0.05
    @test isapprox( gc_time, 0.0) 

  end

  return nothing
end

#test_utils_timing()
add_func1!(EulerTests, test_utils_timing, [TAG_UTILS, TAG_SHORTTEST])

"""
  Test various functions: calcNorm, calcMeshH etc.
"""
function test_utils_misc()
  @testset "----- Testing Utility Functions -----" begin
    M = rand(10)
    Minv = 1./M
    params = TestParams{2}(Timings())
    obj = TestData{Float64, Float64}(params, M, Minv, MPI.COMM_WORLD)
    data = rand(10)

    norm1 = calcNorm(obj, data)
    norm2 = sqrt(sum(data.*M.*data))

    @test isapprox( norm1, norm2) 

    norm1 = calcNorm(obj, data, strongres=true)
    norm2 = sqrt(sum(data.*Minv.*data))


    data = complex.(rand(10), rand(10))
    norm1 = calcNorm(obj, data)
    norm2 = sqrt(sum(real(data).*M.*real(data)))

    @test isapprox( norm1, norm2) 

    norm1 = calcNorm(obj, data, strongres=true)
    norm2 = sqrt(sum(real(data).*Minv.*real(data)))

    @test isapprox( norm1, norm2) 

    numEl = 5
    numNodesPerElement = 3
    numNodes = numEl*numNodesPerElement
    numDofPerNode = 1
    jac = zeros(Float64, numNodesPerElement, numEl)
    fill!(jac, 2)
    dofs = zeros(Int, numDofPerNode, numNodesPerElement, numEl)
    for i=1:length(dofs)
      dofs[i] = i
    end
    mesh = TestMesh{Float64}(jac, dofs, MPI.COMM_WORLD, MPI.Comm_size(MPI.COMM_WORLD), numNodes, numNodesPerElement, numEl, numDofPerNode, 0.25, 2)
    opts = Dict{Any, Any}()
    @test isapprox( calcMeshH(mesh, FakeSBP(), obj, opts), mesh.min_node_dist*1./sqrt(2)) 


  end

  return nothing
end

#test_utils_misc()
add_func1!(EulerTests, test_utils_misc, [TAG_UTILS, TAG_SHORTTEST])

"""
  Test permutation matrices, vectors, functions
"""
function test_utils_perm()
  @testset "----- Testing Permutation Functions -----" begin

    # test permuting rows
    A = [1 2 3 10 ; 4 5 6 11; 7 8 9 12]
    B = zeros(A)
    permvec = [2, 3, 1]
    P = permMatrix(permvec)

    A_orig = copy(A)
    applyPermRowInplace(permvec, A, B)
    A2 = P*A_orig
   
    for i=1:3
      @test ( A2[i, :] )== A[i, :]
    end


    A = [1 2 3; 4 5 6; 7 8 9]
    A_orig = copy(A)
    permvec = [3, 1, 2]
    P = permMatrix(permvec)

    A2 = P*A
    applyPermRowInplace(permvec, A, B)
   
    for i=1:3
      @test ( A2[i, :] )== A[i, :]
    end

    #  test inverse permutation
    invperm = zeros(permvec)
    inversePerm(permvec, invperm)
    applyPermRowInplace(invperm, A, B)

    for i=1:length(A)
      @test ( A[i] )== A_orig[i]
    end

    P2 = permMatrix(invperm)

    @test ( P2.' )== P

    # test permuting columns
    A = [1 2 3; 4 5 6; 7 8 9; 10 11 12]
    B = zeros(A)
    A_orig = copy(A)
    P = permMatrix(permvec)

    applyPermColumnInplace(permvec, A, B)
    A2 = A_orig*P

    for i=1:3
      @test ( (A_orig*P)[:, i] )== A[ :, i ]
    end
  end  # end facts block

  return nothing
end

#test_utils_perm()
add_func1!(EulerTests, test_utils_perm, [TAG_UTILS, TAG_SHORTTEST])

function test_utils_projection()
  @testset "----- testing Projection -----" begin
    t = Timings()
    params = TestParams{2}(t)
    params3 = TestParams{3}(t)

    nrm = zeros(2)
    P = zeros(4,4)
    nrm3 = zeros(3)
    P3 = zeros(5,5)

    q = [1.0, 2.0, 3.0, 7.0]
    qprime = zeros(q)
    qprime2 = zeros(q)
    q2 = zeros(q)
    q2a = zeros(q)
    q3 = [1.0, 2.0, 3.0, 4.0, 15.0]
    q3prime = zeros(q3)
    q3prime2 = zeros(q3)
    q32 = zeros(q3)
    q32a = zeros(q3)

    for theta=0.0:0.1:2*pi
      nrm[1] = cos(theta)
      nrm[2] = sin(theta)

      t1, t2 = Utils.getOrthogonalVector(params, nrm)
      dot_val = abs.(t1*nrm[1] + t2*nrm[2])
      # check normal vector calculation
      @test isapprox( sqrt(t1*t1 + t2*t2), 1.0) atol=1e-12
      @test  dot_val  < 1e-13

      # check projection
      getProjectionMatrix(params, nrm, P)

      # check P is unitary
      @test isapprox( norm(inv(P) - P.'), 0.0) atol=1e-12

      # check q1, q4 are same, magnitude of momentum is same
      A_mul_B!(qprime, P, q)
      @test isapprox( q[1], qprime[1]) atol=1e-13
      @test isapprox( q[4], qprime[4]) atol=1e-13
      @test isapprox( q[2]*q[2] + q[3]*q[3], qprime[2]*qprime[2] + qprime[3]*qprime[3]) atol=1e-12

      # check multiplication functions against A_mul_b
      projectToNT(params, P, q, qprime2)
      @test isapprox( qprime, qprime2) atol=1e-12

      # project back
      At_mul_B!(q2, P, qprime)
      @test isapprox( q2, q) atol=1e-12
      projectToXY(params, P, qprime, q2a)
      @test isapprox( q2a, q) atol=1e-12


    end

    # 3D test
    for theta=0:0.1:pi
      for phi=0:0.1:2*pi
        nrm3[1] = sin(theta)*cos(phi)
        nrm3[2] = sin(theta)*sin(phi)
        nrm3[3] = cos(theta)

        t1, t2, t3 = Utils.getOrthogonalVector(params3, nrm3)
        dot_val = abs.(t1*nrm3[1] + t2*nrm3[2] + t3*nrm3[3])
        @test isapprox( t1*t1 + t2*t2 + t3*t3, 1.0) atol=1e-12
        @test isapprox( dot_val, 0.0) atol=1e-13

        getProjectionMatrix(params3, nrm3, P3)

        @test isapprox( norm(inv(P3) - P3.'), 0.0) atol=1e-12

        # check q1, q4 are same, magnitude of momentum is same
        A_mul_B!(q3prime, P3, q3)
        @test isapprox( q3[1], q3prime[1]) atol=1e-13
        @test isapprox( q3[5], q3prime[5]) atol=1e-13
        @test isapprox( q3[2]*q[2] + q3[3]*q3[3] + q3[4]*q3[4], q3prime[2]*q3prime[2] + q3prime[3]*q3prime[3] + q3prime[4]*q3prime[4]) atol=1e-12

        # check multiplication functions against A_mul_b
        projectToNT(params3, P3, q3, q3prime2)
        @test isapprox( q3prime, q3prime2) atol=1e-12

        # project back
        At_mul_B!(q32, P3, q3prime)
        @test isapprox( q32, q3) atol=1e-12
        projectToXY(params3, P3, q3prime, q32a)
        @test isapprox( q32a, q3) atol=1e-12

      end
    end

  end  # end facts block

  return nothing
end

add_func1!(EulerTests, test_utils_projection, [TAG_UTILS, TAG_SHORTTEST])
