const test_checkpoint_inputfile = "input_vals_channel_dg.jl"
#const test_checkpoint_moddict = Dict{String, Any}("Flux_name" => "RoeFlux", "use_DG" => true, "new_fname" => "input_vals_channel_dg")

mutable struct TestCheckpointData <: AbstractCheckpointData
  a::Int
  b::Array{Float64, 1}
end

function test_checkpoint(mesh, sbp, eqn, opts)
  @testset "----- Testing Checkpointer -----" begin
    myrank = mesh.myrank
    rand!(eqn.q_vec)

    chkpointer = Checkpointer(myrank, 2, "abc")

    chkpoint_data = TestCheckpointData(42, rand(10))

    # check initialization
    @test ( chkpointer.ncheckpoints )== 2
    @test ( chkpointer.paths[2] )== joinpath(pwd(), "abc_checkpoint2")
    @test ( isdir(joinpath(pwd(), chkpointer.paths[2])) )== true
    for i=1:chkpointer.ncheckpoints
      @test ( chkpointer.status[i] )== Utils.CheckpointFree
      @test ( chkpointer.history[i] )== -1
    end

    @test ( countFreeCheckpoints(chkpointer) )== 2
    @test ( getNextFreeCheckpoint(chkpointer) )== 1

    # save a checkpoint
    @test ( saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data) )== 1

    @test ( countFreeCheckpoints(chkpointer) )== 1
    @test ( getNextFreeCheckpoint(chkpointer) )== 2
    @test ( getOldestCheckpoint(chkpointer) )== 1
    @test ( getLastCheckpoint(chkpointer) )== 1
    @test ( chkpointer.status[1] )== Utils.CheckpointUsed
    @test ( chkpointer.history[1] )== 1

    q1 = copy(eqn.q_vec)
    scale!(eqn.q_vec, 2)
    q2 = copy(eqn.q_vec)
    @test ( saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data) )== 2
    @test ( countFreeCheckpoints(chkpointer) )== 0
    @test ( chkpointer.status[2] )== Utils.CheckpointUsed

    @test ( getOldestCheckpoint(chkpointer) )== 1
    @test_throws Exception  saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)

    # try loading a checkpoint
    scale!(eqn.q_vec, 2)
    loadLastCheckpoint(chkpointer, mesh, sbp, eqn, opts)
    @test isapprox( norm(eqn.q_vec - q2), 0.0) atol=1e-14

    chkpointer2 = Checkpointer(opts, mesh.myrank)

    @test ( chkpointer.ncheckpoints )== chkpointer2.ncheckpoints
    @test isapprox( norm(chkpointer.history - chkpointer2.history), 0.0) atol=1e-14
    @test isapprox( norm(chkpointer.status - chkpointer2.status), 0.0) atol=1e-14
    for i=1:chkpointer.ncheckpoints
      @test ( chkpointer.paths[i] )== chkpointer2.paths[i]
    end

    chkpoint_data2 = readCheckpointData(chkpointer, getLastCheckpoint(chkpointer), mesh.myrank)

    @test ( chkpoint_data.a )== chkpoint_data2.a
    @test isapprox( norm(chkpoint_data.b - chkpoint_data2.b ), 0.0) atol=1e-14

    # corrupt the flagfile to force the loading of the earlier checkpoint
    Utils.deleteFlagFile(chkpointer, getLastCheckpoint(chkpointer))
    opts["most_recent_checkpoint"] = 1
    opts["most_recent_checkpoint_path"] = chkpointer.paths[1]

    chkpointer2 = Checkpointer(opts, mesh.myrank)
    @test ( countFreeCheckpoints(chkpointer2) )== 1
    @test ( getLastCheckpoint(chkpointer2) )== 1


    # test with more checkpoints
    chkpointer = Checkpointer(myrank, 5)
    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)
    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)
    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)
    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)

    freeCheckpoint(chkpointer, 2)
    @test ( countFreeCheckpoints(chkpointer) )== 2
    @test ( getNextFreeCheckpoint(chkpointer) )== 2
    @test ( getLastCheckpoint(chkpointer) )== 4
    @test ( getOldestCheckpoint(chkpointer) )== 1

    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)
    @test ( countFreeCheckpoints(chkpointer) )== 1
    @test ( getNextFreeCheckpoint(chkpointer) )== 5
    @test ( getLastCheckpoint(chkpointer) )== 2
    @test ( getOldestCheckpoint(chkpointer) )== 1

  end

  return nothing
end

add_func2!(EulerTests, test_checkpoint, test_dg_inputfile, [TAG_SHORTTEST, TAG_CHECKPOINT])
