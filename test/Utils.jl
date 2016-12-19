#=
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))


using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
include( joinpath(Pkg.dir("PDESolver"), "src/solver/euler/complexify.jl"))
include( joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))


push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using Utils
using FactCheck
using ODLCommonTools
using MPI

if !MPI.Initialized()
  MPI.Init()
end
=#

type TestParams
  time::Timings
end

# test the tools
type TestData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
  params::TestParams
  M::Array{Float64, 1}
  Minv::Array{Float64, 1}
  comm::MPI.Comm
end

type TestMesh{Tmsh} <: AbstractMesh{Tmsh}
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

type FakeSBP
end


function test_utils_io()
  facts("---- Testing IO -----") do
    fname = "iotest.dat"
    f = open(fname, "w")
    f2 = BufferedIO(f)

    println(f2, 1)
    println(f2, 2)

    fstat = stat(fname)
    @fact fstat.size --> 0
    flush(f2)
  #  close(f)

    data = readdlm(fname)
    @fact data[1] --> 1
    @fact data[2] --> 2

    close(f)
  end

  return nothing
end

#test_utils_io()
add_func1!(EulerTests, test_utils_io)

function test_utils_timing()
  facts("----- Testing Timings -----") do
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
      @fact vals[i] --> roughly(Float64(i))
    end
  =#

    function mysleep(t)
      sleep(t)
      return 42
    end

    @time_all mysleep(0.0001)
    gc()
    val, telapsed, gc_time, alloced = @time_all mysleep(0.1)
    @fact val --> 42
    @fact telapsed --> roughly(0.1, atol=0.05)
    @fact gc_time --> roughly( 0.0)

  end

  return nothing
end

#test_utils_timing()
add_func1!(EulerTests, test_utils_timing)

function test_utils_misc()
  facts("----- Testing Utility Functions -----") do
    M = rand(10)
    Minv = 1./M
    params = TestParams(Timings())
    obj = TestData{Float64, Float64}(params, M, Minv, MPI.COMM_WORLD)
    data = rand(10)

    norm1 = calcNorm(obj, data)
    norm2 = sqrt(sum(data.*M.*data))

    @fact norm1 --> roughly(norm2)

    norm1 = calcNorm(obj, data, strongres=true)
    norm2 = sqrt(sum(data.*Minv.*data))


    data = complex(rand(10), rand(10))
    norm1 = calcNorm(obj, data)
    norm2 = sqrt(sum(real(data).*M.*real(data)))

    @fact norm1 --> roughly(norm2)

    norm1 = calcNorm(obj, data, strongres=true)
    norm2 = sqrt(sum(real(data).*Minv.*real(data)))

    @fact norm1 --> roughly(norm2)

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
    @fact calcMeshH(mesh, FakeSBP(), obj, opts) --> roughly(mesh.min_node_dist*1./sqrt(2))


  end

  return nothing
end

#test_utils_misc()
add_func1!(EulerTests, test_utils_misc)

function test_utils_perm()
  facts("----- Testing Permutation Functions -----") do

    # test permuting rows
    A = [1 2 3 10 ; 4 5 6 11; 7 8 9 12]
    B = zeros(A)
    permvec = [2, 3, 1]
    P = permMatrix(permvec)

    A_orig = copy(A)
    applyPermRowInplace(permvec, A, B)
    A2 = P*A_orig
   
    for i=1:3
      @fact A2[i, :] --> A[i, :]
    end


    A = [1 2 3; 4 5 6; 7 8 9]
    A_orig = copy(A)
    permvec = [3, 1, 2]
    P = permMatrix(permvec)

    A2 = P*A
    applyPermRowInplace(permvec, A, B)
   
    for i=1:3
      @fact A2[i, :] --> A[i, :]
    end




    #  test inverse permutation
    invperm = zeros(permvec)
    inversePerm(permvec, invperm)
    applyPermRowInplace(invperm, A, B)

    for i=1:length(A)
      @fact A[i] --> A_orig[i]
    end

    P2 = permMatrix(invperm)

    @fact P2.' --> P

    # test permuting columns
    A = [1 2 3; 4 5 6; 7 8 9; 10 11 12]
    B = zeros(A)
    A_orig = copy(A)
    P = permMatrix(permvec)

    applyPermColumnInplace(permvec, A, B)
    A2 = A_orig*P

    for i=1:3
      @fact (A_orig*P)[:, i] --> A[ :, i ]
    end
  end  # end facts block

  return nothing
end

#test_utils_perm()
add_func1!(EulerTests, test_utils_perm)


