const test_checkpoint_inputfile = "input_vals_channel.jl"
const test_checkpoint_moddict = Dict{ASCIIString, Any}("Flux_name" => "RoeFlux", "use_DG" => true, "new_fname" => "input_vals_channel_dg")

type TestCheckpointData <: AbstractCheckpointData
  a::Int
  b::Array{Float64, 1}
end

function test_checkpoint(mesh, sbp, eqn, opts)
  facts("----- Testing Checkpointer -----") do
    rand!(eqn.q_vec)

    chkpointer = Checkpointer(2, "abc")

    chkpoint_data = TestCheckpointData(42, rand(10))

    # check initialization
    @fact chkpointer.ncheckpoints --> 2
    @fact chkpointer.paths[2] --> joinpath(pwd(), "abc_checkpoint2")
    @fact isdir(joinpath(pwd(), chkpointer.paths[2])) --> true
    for i=1:chkpointer.ncheckpoints
      @fact chkpointer.status[i] --> Utils.CheckpointFree
      @fact chkpointer.history[i] --> -1
    end

    @fact countFreeCheckpoints(chkpointer) --> 2
    @fact getNextFreeCheckpoint(chkpointer) --> 1

    # save a checkpoint
    @fact saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data) --> 1

    @fact countFreeCheckpoints(chkpointer) --> 1
    @fact getNextFreeCheckpoint(chkpointer) --> 2
    @fact getOldestCheckpoint(chkpointer) --> 1
    @fact getLastCheckpoint(chkpointer) --> 1
    @fact chkpointer.status[1] --> Utils.CheckpointUsed
    @fact chkpointer.history[1] --> 1

    q1 = copy(eqn.q_vec)
    scale!(eqn.q_vec, 2)
    q2 = copy(eqn.q_vec)
    @fact saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data) --> 2
    @fact countFreeCheckpoints(chkpointer) --> 0
    @fact chkpointer.status[2] --> Utils.CheckpointUsed

    @fact getOldestCheckpoint(chkpointer) --> 1
    @fact_throws saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)

    # try loading a checkpoint
    scale!(eqn.q_vec, 2)
    loadLastCheckpoint(chkpointer, mesh, sbp, eqn, opts)
    @fact norm(eqn.q_vec - q2) --> roughly(0.0, atol=1e-14)

    chkpointer2 = Checkpointer(opts)

    @fact chkpointer.ncheckpoints --> chkpointer2.ncheckpoints
    @fact norm(chkpointer.history - chkpointer2.history) --> roughly(0.0, atol=1e-14)
    @fact norm(chkpointer.status - chkpointer2.status) --> roughly(0.0, atol=1e-14)
    for i=1:chkpointer.ncheckpoints
      @fact chkpointer.paths[i] --> chkpointer2.paths[i]
    end

    chkpoint_data2 = readCheckpointData(chkpointer, getLastCheckpoint(chkpointer), mesh.myrank)

    @fact chkpoint_data.a --> chkpoint_data2.a
    @fact norm(chkpoint_data.b - chkpoint_data2.b ) --> roughly(0.0, atol=1e-14)

    # corrupt the flagfile to force the loading of the earlier checkpoint
    Utils.deleteFlagFile(chkpointer, getLastCheckpoint(chkpointer))
    opts["most_recent_checkpoint"] = 1
    opts["most_recent_checkpoint_path"] = chkpointer.paths[1]

    chkpointer2 = Checkpointer(opts)
    @fact countFreeCheckpoints(chkpointer2) --> 1
    @fact getLastCheckpoint(chkpointer2) --> 1


    # test with more checkpoints
    chkpointer = Checkpointer(5)
    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)
    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)
    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)
    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)

    freeCheckpoint(chkpointer, 2)
    @fact countFreeCheckpoints(chkpointer) --> 2
    @fact getNextFreeCheckpoint(chkpointer) --> 2
    @fact getLastCheckpoint(chkpointer) --> 4
    @fact getOldestCheckpoint(chkpointer) --> 1

    saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpoint_data)
    @fact countFreeCheckpoints(chkpointer) --> 1
    @fact getNextFreeCheckpoint(chkpointer) --> 5
    @fact getLastCheckpoint(chkpointer) --> 2
    @fact getOldestCheckpoint(chkpointer) --> 1

  end

  return nothing
end

add_func3!(EulerTests, test_checkpoint, test_dg_inputfile, test_dg_moddict, [TAG_SHORTTEST, TAG_CHECKPOINT])
