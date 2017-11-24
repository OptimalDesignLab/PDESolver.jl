# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson, cnResidual

#TODO: stop doing this
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

@doc """
crank_nicolson

  This function performs Crank-Nicolson implicit time solution, using a function 
  of the form du/dt = f(u, t)
  
  Arguments:
    * f  : function evaluation, must have signature (ctx..., opts, t)
    * h  : time step size
    * t_max: time value to stop time stepping (time starts at 0)
    * q_vec: vector of the u values
    * res_vec: vector of du/dt values (the output of the function f)
    * pre_func: function to to be called after the new u values are put into
                q_vec but before the function f is evaluated.  Must have
                signature: pre_func(ctx..., opts)
    * post_func: function called immediately after f is called.  The function
                 must have the signature res_norm = post_func(ctx..., opts, 
                 calc_norm=true),
                 where res_norm is a norm of res_vec, and calc_norm determines
                 whether or not to calculate the norm.
    * ctx: a tuple (or any iterable container) of the objects needed by
           f, pre_func, and post func.  The tuple is splatted before being
           passed to the functions.
    * opts : options dictionary

    Keyword Arguments:
    * majorIterationCallback: a callback function called after the first
                              stage, useful to do output or logging
    * res_tol : keyword arg, residual topping tolerance
    * real_time : do actual time marching, not pseudo-time marching

   The eqn.q_vec should hold the whichever variables (conservative or
   entropy) that the simulation should use.

   The idea behind pre_func, post func, and ctx is that they abstract what kind
   of system rk4 is timestepping.  rk4 only needs to know about q_vec and
   res_vec.

   For physics modules, ctx should be (mesh, sbp, eqn) and q_vec and res_vec 
   should be eqn.q_vec and eqn.res_vec.
"""->
function crank_nicolson(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                        opts, res_tol=-1.0, real_time=true)

  myrank = eqn.myrank
  if myrank == 0
    println(BSTDOUT, "\nEntered Crank-Nicolson")
    println(BSTDOUT, "res_tol = ", res_tol)
  end
  flush(BSTDOUT)

  output_freq = opts["output_freq"]::Int
  write_vis = opts["write_vis"]::Bool
  use_itermax = opts["use_itermax"]::Bool
  jac_type = opts["jac_type"]
  if use_itermax
    itermax = opts["itermax"]
  end

  use_checkpointing = opts["use_checkpointing"]::Bool
  chkpoint_freq = opts["checkpoint_freq"]::Int
  ncheckpoints = opts["ncheckpoints"]::Int


#  if jac_type == 4
#    throw(ErrorException("CN not implemented for matrix-free ops. (jac_type cannot be 4)"))
#  end

  @mpi_master if myrank == 0
    _f1 = open("convergence.dat", "a")
    f1 = BufferedIO(_f1)
  end

  t = 0.0
  t_steps = round(Int, t_max/h)

  # eqn_nextstep = deepcopy(eqn)
  eqn_nextstep = eqn_deepcopy(mesh, sbp, eqn, opts)

  # TODO: copyForMultistage! does not give correct values.
  #     deepcopy works for now, but uses more memory than copyForMultistage!, if it worked
  # eqn_nextstep = copyForMultistage!(eqn)
  eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

  @debug1 println("============ In CN ============")

   # setup all the checkpointing related data
  chkpointer, chkpointdata, skip_checkpoint = explicit_checkpoint_setup(opts, myrank)
  istart = chkpointdata.i

  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop

  # NOTE: eqn_nextstep changed to eqn 20161013
  pc, lo = getCNPCandLO(mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm)
  newton_data, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, ls)

  # this loop is 2:(t_steps+1) when not restarting
  for i = istart:(t_steps + 1)

    t = (i-2)*h
    @mpi_master println(BSTDOUT, "\ni = ", i, ", t = ", t)
    @debug1 println(eqn.params.f, "====== CN: at the top of time-stepping loop, t = $t, i = $i")
    @debug1 flush(eqn.params.f)

    if use_checkpointing && i % chkpoint_freq == 0 && !skip_checkpoint
      @mpi_master println(BSTDOUT, "Saving checkpoint at timestep ", i)
      skip_checkpoint = false
      # save all needed variables to the chkpointdata
      chkpointdata.i = i

      if countFreeCheckpoints(chkpointer) == 0
        freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint
      end
      
      # save the checkpoint
      saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpointdata)
    end

#=
    #----------------------------
    # zero out Jac
    #   this works for both PETSc and Julia matrices.
    #   when jac is a Julia matrix, this is effectively wrapping: fill!(jac, 0.0)
    if jac_type != 4
      PetscMatZeroEntries(jac)
    end
=#
    # TODO: Allow for some kind of stage loop: ES-Dirk

    # TODO: output freq

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

    t_nextstep = t + h

    # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
    if opts["cleansheet_CN_newton"]
      cnNewton(mesh, sbp, opts, h, f, eqn, eqn_nextstep, t)
    else

      ctx_residual = (f, eqn, h, newton_data)
      newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, ls, 
                  rhs_vec, ctx_residual, t, itermax=30, 
                  step_tol=opts["step_tol"], 
                  res_abstol=opts["res_abstol"], res_reltol=opts["res_reltol"],                   res_reltol0=opts["res_reltol0"])
    end

    # do the callback using the current eqn object at time t
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts, STDOUT)

    # need to assemble solution into res_vec?
    res_norm = calcNorm(eqn, eqn.res_vec)
    # logging
    @mpi_master if i % 1 == 0
      println(f1, i, " ", res_norm)
    end
    
    @mpi_master if i % output_freq == 0
      println(BSTDOUT, "flushing convergence.dat to disk")
      flush(f1)
    end

    # check stopping conditions
    if (res_norm < res_tol)  # ???
      if myrank == 0
        println(BSTDOUT, "breaking due to res_tol, res norm = $res_norm")
        close(f1)
        flush(BSTDOUT)
      end
      break
    end


    if use_itermax && i > itermax
      if myrank == 0
        println(BSTDOUT, "breaking due to itermax")
        close(f1)
        flush(BSTDOUT)
      end
      break
    end


    # This allows the solution to be updated from _nextstep without a deepcopy.
    # There are two memory locations used by eqn & eqn_nextstep, 
    #   and this flips the location of eqn & eqn_nextstep every time step
    eqn_temp = eqn
    eqn = eqn_nextstep
    eqn_nextstep = eqn_temp

    # Note: we now need to copy the updated q over for the initial newton guess
    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] = eqn.q_vec[i]
    end
    disassembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q, eqn_nextstep.q_vec)

#    t = t_nextstep
    flush(BSTDOUT)

  end   # end of t step loop

  # final time update
  t += h

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copyForMultistage!(dest, src)
  copyForMultistage!(eqn, eqn_nextstep)


  #TODO: return the NewtonData?
  free(newton_data)

  @debug1 println("============= end of CN: t = $t ===============")
  return t

  flush(BSTDOUT)

end   # end of crank_nicolson function

"""
  Construct CN preconditioner and linear operator based on options dictionary

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function getCNPCandLO(mesh, sbp, eqn, opts)

  jac_type = opts["jac_type"]

  if jac_type <= 2
    pc = PCNone(mesh, sbp, eqn, opts)
  elseif opts["use_volume_preconditioner"]
    pc = CNVolumePC(mesh, sbp, eqn, opts)
  else
    pc = CNMatPC(mesh, sbp, eqn, opts)
  end

  if jac_type == 1
    lo = CNDenseLO(pc, mesh, sbp, eqn, opts)
  elseif jac_type == 2
    lo = CNSparseDirectLO(pc, mesh, sbp, eqn, opts)
  elseif jac_type == 3
    lo = CNPetscMatLO(pc, mesh, sbp, eqn, opts)
  elseif jac_type == 4
    lo = CNPetscMatFreeLO(pc, mesh, sbp, eqn, opts)
  end

  return pc, lo
end
#------------------------------------------------------------------------------
# Regular (matrix-based) PC

type CNMatPC <: AbstractPetscMatPC
  pc_inner::NewtonMatPC
end

function CNMatPC(mesh::AbstractMesh, sbp::AbstractSBP,
                    eqn::AbstractSolutionData, opts::Dict)

  pc_inner = NewtonMatPC(mesh, sbp, eqn, opts)

  return CNMatPC(pc_inner)
end

function calcPC(pc::CNMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  calcPC(pc, mesh, sbp, eqn, opts, ctx_residual, t)
  modifyJacCN(pc, mesh, sbp, eqn, opts, ctx_residual, t)

  return nothing
end


#------------------------------------------------------------------------------
# Volume preconditioner

"""
  Volume integral jacobian preconditioner for CN.

"""
type CNVolumePC <: AbstractPetscMatFreePC
  pc_inner::PetscMatFreePC
  other_pc::NewtonVolumePC
end  # end type definition

function CNVolumePC(mesh::AbstractMesh, sbp,
                              eqn::AbstractSolutionData, opts)

  pc_inner = PetscMatFreePC(mesh, sbp, eqn, opts)
  other_pc = NewtonVolumePC(mesh, sbp, eqn, opts)

  return CNVolumePC(pc_inner, other_pc)
end

function calcPC(pc::CNVolumePC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

# eqn should be eqn_nextstep from the CN function

  cnCalcVolumePreconditioner(pc, mesh, sbp, eqn, opts, ctx_residual, t)
  factorVolumePreconditioner(pc, mesh, sbp, eqn, opts)

  # this needs to come last (because calcVolumePreconditioner calls it too)
  setPCCtx(pc, mesh, sbp, eqn, opts, ctx_residual, t)

  return nothing
end

function applyPC(pc::CNVolumePC, mesh::AbstractMesh, sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)

  applyVolumePreconditioner(pc.other_pc, mesh, sbp, eqn, opts, b, x)

  return nothing
end

#------------------------------------------------------------------------------
# Linear operators

#TODO: make a macro to generate these definitions

type CNDenseLO <: AbstractDenseLO
  lo_inner::NewtonDenseLO
end

function CNDenseLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonDenseLO(pc, mesh, sbp, eqn, opts)

  return CNDenseLO(lo_inner)
end

type CNSparseDirectLO <: AbstractSparseDirectLO
  lo_inner::NewtonSparseDirectLO
end

function CNSparseDirectLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonSparseDirectLO(pc, mesh, sbp, eqn, opts)

  return CNSparseDirectLO(lo_inner)
end


type CNPetscMatLO <: AbstractPetscMatLO
  lo_inner::NewtonPetscMatLO
end

function CNPetscMatLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonPetscMatLO(pc, mesh, sbp, eqn, opts)

  return CNPetscMatLO(lo_inner)
end

"""
  All CN LOs that have matrices
"""
typealias CNMatLO Union{CNDenseLO, CNSparseDirectLO, CNPetscMatLO}


function calcLinearOperator(lo::CNMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  h = ctx_residual[3]
  t_nextstep = t + h
  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t_nextstep)

  modifyJacCN(lo, mesh, sbp, eqn, opts, ctx_residual, t_nextstep)

  return nothing
end


# matrix-free
type CNPetscMatFreeLO <: AbstractPetscMatFreeLO
  lo_inner::NewtonPetscMatFreeLO
end

function CNPetscMatFreeLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                      sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonPetscMatFreeLO(pc, mesh, sbp, eqn, opts)

  return CNPetscMatFreeLO(lo_inner, rhs_func)
end

"""
  Any PC or LO that has a matrix in the field `A`
"""
typealias CNHasMat Union{CNMatPC, CNDenseLO, CNSparseDirectLO, CNPetscMatLO}


function calcLinearOperator(lo::NewtonPetscMatFreeLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)
  h = ctx_residual[3]
  t_nextstep = t + h

  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t_nextstep)

  setLOCtx(lo, mesh, sbp, eqn, opts, ctx_residual, t_nextstep)

  return nothing
end

function applyLinearOperator{Tsol}(lo::CNPetscMatFreeLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  h = ctx[3]
  # the CN residual has the form: I - 0.5*delta_t*jac, where jac is the physics
  # Jacobian.
  # thus a matrix vector product is: v - 0.5*delta_t*jac*v
  applyLinearOperator(lo, mesh, sbp, eqn, opts, ctx_residual, t, x, b)

  scale!(b, -0.5*h)
  @simd for i=1:length(b)
    b[i] += x[i]
  end

  return nothing
end

function applyLinearOperatorTranspose{Tsol}(lo::CNPetscMatFreeLO,
                             mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  # see the non-transpose method for an explanation of the math

  h = ctx[3]
  applyLinearOperatorTranspose(lo, mesh, sbp, eqn, opts, ctx_residual, t, x, b)

  scale!(b, -0.5*h)
  @simd for i=1:length(b)
    b[i] += x[i]
  end

  return nothing
end



"""
  Takes the Jacobian of the physics and modifies it to be the Crank-Nicolson
  Jacobian.

  CN_Jac = I - dt/2 * physicsJac


  **Inputs**

   * lo: a CNHasMat
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: delta_t must be the 3rd element
   * t
"""
function modifyJacCN(lo::CNHasMat, mesh, sbp, eqn, opts, ctx_residual, t)


  lo2 = getBaseLO(lo)
  h = ctx_residual[3]

  PetscMatAssemblyBegin(lo2.A, PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(lo2.A, PETSC_MAT_FINAL_ASSEMBLY)

  # scale jac by -delta_t/2
#  scale_factor = h*-0.5
  petsc_scale_factor = PetscScalar(-h*0.5)
  PetscMatScale(lo2.A, petsc_scale_factor)

  # add the identity
  PetscMatShift(lo2.A, PetscScalar(1))

  return nothing
end

#=
@doc """
###NonlinearSolvers.cnJac

  Jac of the CN calculation.
  Effectively a wrapper for physicsJac, because the CN Jac is:
    CN_Jac = I - dt/2 * physicsJac

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element
    newton_data must be the fourth element
"""->
function cnJac(newton_data, mesh::AbstractMesh, sbp::AbstractSBP,
               eqn_nextstep::AbstractSolutionData, opts, jac, ctx, t)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)

  # should t be t_nextstep?
  physics_func = ctx[1]
  # NOTE: eqn instead of eqn_nextstep, 20161013
  eqn = ctx[2] # unused?
  h = ctx[3]
#  newton_data = ctx[4]  # ??? this is an argumetn

  jac_type = opts["jac_type"]

  t_nextstep = t + h

  # don't compute the jacobian
  #TODO: support using a preconditioning matrix
  if jac_type == 4
    return nothing
  end
  # Forming the CN Jacobian:
  #   call physicsJac with eqn_nextstep & t_nextstep
  #   then form CN_Jac = I + dt/2 * physics_Jac

  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t_nextstep)


  # need to flush assembly cache before performing the scale operation.
  #   These are defined for Julia matrices; they're just noops
  PetscMatAssemblyBegin(jac, PETSc.PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(jac, PETSc.PETSC_MAT_FINAL_ASSEMBLY)

  #--------------------------
  # applying dt/2 to jac
  # Jacobian is always 2D
  scale_factor = h*-0.5
  # make this into a petsc_scale_factor
  petsc_scale_factor = PetscScalar(scale_factor)

  # PetscMatScale is defined for all jac types, PETSc and Julia
  # when jac is julia array, this effectively does: scale!(jac, scale_factor)
  PetscMatScale(jac, petsc_scale_factor)

  # PetscMatAssembly___ not necessary here; PetscMatScale is provably local so it doesn't cache its stuff

  #--------------------------
  # adding identity
  ix_petsc_row = zeros(PetscInt, 1)
  ix_petsc_col = zeros(PetscInt, 1)
  value_to_add = zeros(PetscScalar, 1, 1)
  value_to_add[1,1] = 1.0
  flag = PETSc.PETSC_ADD_VALUES

  for i = 1:mesh.numDof
    ix_petsc_row[1,1] = i + mesh.dof_offset       # jac diag index + local to global offset
    ix_petsc_col[1,1] = i + mesh.dof_offset       # jac diag index + local to global offset

    # PETSc function: set_values1!(Jac, [2], [2], Jac[2,2] + 1)
    # args: array, row (as an array of len 1), col (as an array of len 1), new values (as a 2D array of len 1)
    #   set_values1! has different methods for both PETSc matrices and Julia matrices
    #   for Julia dense & sparse arrays, set_values1! effectively does this: jac[i,i] += 1
    #   in serial, mesh.dof_offset is set to 0 automatically
    set_values1!(jac, ix_petsc_row, ix_petsc_col, value_to_add, flag)
  end

  # set_values1! only caches the results; need to be assembled. This happens in petscSolve in petsc_funcs.jl
  #   (The assemble funcs are defined for Julia matrices; they're just noops)

  # jac is now I + dt/2 * physics_jac

  return nothing

end
=#
#=
"""
  This wrapper around [`cnVolumePCSetUp`](@ref) is passed into Petsc as a
  callback

  **Inputs**

   * pc: a Shell PC object

  **Notes**

  the PC ctx must be: (mesh, sbp, eqn, opts, newton_data, func, ctx_residual,
  t), where ctx_residual is the ctx used by all the CN functions and
  func is the right hand side function passed into [`newtonInner`](@ref).

"""
function cnVolumePCSetUp_wrapper(pc::PC)
  ctx_ptr = PCShellGetContext(pc)
  ctx_petsc_pc = unsafe_pointer_to_objref(ctx_ptr)

  # unpack ctx
  # the ctx_petsc_pc is the same as for a shell matrix for convenience
  mesh = ctx_petsc_pc[1]
  sbp = ctx_petsc_pc[2]
  eqn = ctx_petsc_pc[3]
  opts = ctx_petsc_pc[4]
  newton_data = ctx_petsc_pc[5]
  func = ctx_petsc_pc[6]  # rhs_func from newtonInner
  ctx_residual = ctx_petsc_pc[7]
  t = ctx_petsc_pc[8]

  cnVolumePCSetUp(newton_data, mesh, sbp, eqn, opts, ctx_residual, t)

  return PetscErrorCode(0)
end
=#
#=
"""
  This function calculates and factors the volume preconditioner.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: the Crank-Nicolson ctx used by all CN functions
   * t: current time

  **Inputs/Outputs**

   * newton_data: vol_prec is updated with the new preconditioner
"""
function cnVolumePCSetUp(newton_data, mesh::AbstractDGMesh,
                         sbp::AbstractSBP, eqn_nextstep::AbstractSolutionData,
                         opts, ctx_residual, t)

  cnCalcVolumePreconditioner(newton_data, mesh, sbp, eqn_nextstep, opts,
                             ctx_residual, t)
  factorVolumePreconditioner(newton_data, mesh, sbp, eqn_nextstep, opts)

  return nothing
end
=#

"""
  This function uses [`calcVolumePreconditioner`](@ref) to compute a
  a preconditioner for the matrix used by Crank-Nicolson.  This function
  has the same signature as [`cnJac`](@ref) to facilitate matrix-free
  operations.
"""
function cnCalcVolumePreconditioner(pc::CNVolumePC, mesh::AbstractDGMesh,
               sbp::AbstractSBP, eqn_nextstep::AbstractSolutionData, opts,
               ctx_residual, t)
  #TODO: is this ctx_residual?
  physics_func = ctx[1]
#  eqn = ctx[2]
  h = ctx[3]
  t_nextstep = t + h

  epsilon = opts["epsilon"]
  pert = Complex128(0, epsilon)

  other_pc = pc.other_pc

#  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t_nextstep)
  calcVolumePreconditioner(other_pc, mesh, sbp, eqn_nextstep, opts, pert, 
                           physics_func, t_nextstep)

  volume_prec = other_pc.vol_prec
  scale!(other_pc.volume_jac, -0.5*h)

  # add the identity matrix
  jac_size = mesh.numDofPerNode*mesh.numNodesPerElement

  for i=1:mesh.numEl
    for j=1:jac_size
      other_pc.volume_jac[j, j, i] += 1
    end
  end

  return nothing
end




@doc """
###NonlinearSolvers.cnRhs

  RHS of the CN calculation

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element

"""->
function cnRhs(mesh::AbstractMesh, sbp::AbstractSBP, eqn_nextstep::AbstractSolutionData, opts, rhs_vec, ctx, t)

  # eqn comes in through ctx_residual, which is set up in CN before the newtonInner call

  physics_func = ctx[1]
  eqn = ctx[2]
  h = ctx[3]

  t_nextstep = t + h

  # evalute residual at t_nextstep
  # q_vec -> q
  disassembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q_vec)

  # start parallel communication if needed
  time = eqn.params.time
  time.t_send += @elapsed if opts["parallel_type"] == 2 && mesh.npeers > 0
    startSolutionExchange(mesh, sbp, eqn_nextstep, opts)
  end


  physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
  assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, 
                   eqn_nextstep.res_vec)

  # evalute residual at t
  if opts["parallel_type"] == 2 && mesh.npeers > 0
    startSolutionExchange(mesh, sbp, eqn, opts)
  end
  physics_func(mesh, sbp, eqn, opts, t)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  # compute rhs

  #   what this is doing:
  #   u_(n+1) - 0.5*dt* (del dot G_(n+1)) - u_n - 0.5*dt* (del dot G_n)
  for i = 1:mesh.numDof

    temp1 = eqn_nextstep.q_vec[i] - 0.5*h*eqn_nextstep.res_vec[i]
    temp2 = eqn.q_vec[i] + 0.5*h*eqn.res_vec[i]

    rhs_vec[i] = temp1 - temp2 

    # NOTE: question: is there a sign problem here? should rhs_vec = -rhs_vec ?
    #     NO. this negative gets applied in newton.jl, where res_0[i] = -res_0[i]

  end

  # calculate the norm
  rhs_0_norm = calcNorm(eqn, rhs_vec, strongres=true)


  return rhs_0_norm

end     # end of function cnRhs

# the goal is to replace newton.jl.
# this will go into CN in the time-stepping loop
function cnNewton(mesh, sbp, opts, h, physics_func, eqn, eqn_nextstep, t)
  println("++++++++++++++++ clean sheet Newton being run ++++++++++")

  println("---- physics_func: ",physics_func)

  # Jac on eqn or eqn_nextstep?

  epsilon = 1e-8
  t_nextstep = t + h

  jac = zeros(mesh.numDof, mesh.numDof)

  # emulates physicsJac
  # so we need to step through the jacobian column wise.
  #   d res[1]/d q[1]   d res[1]/d q[2]   d res[1]/d q[3] ...
  #   d res[2]/d q[1]   d res[2]/d q[2]   d res[2]/d q[3] ...
  #   d res[3]/d q[1]   d res[3]/d q[2]   d res[3]/d q[3] ...
  #   ...               ...               ...

  newton_itermax = 2
  delta_q_vec = zeros(eqn_nextstep.q_vec)

  # newton_loop starting here?
  for newton_i = 1:newton_itermax

    #--------------------------
    # emulates physicsJac
    unperturbed_q_vec = copy(eqn_nextstep.q_vec)

    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    # needed b/c physics_func only updates eqn.res
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
    # Comment here about mass matrix inv multiplication TODO
    applyMassMatrixInv(mesh, eqn_nextstep, eqn_nextstep.res_vec)
    unperturbed_res_vec = copy(eqn_nextstep.res_vec)

    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] = eqn_nextstep.q_vec[i] + epsilon

      physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
      assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
      applyMassMatrixInv(mesh, eqn, eqn_nextstep.res_vec)

      jac[:,i] = (eqn_nextstep.res_vec - unperturbed_res_vec)/epsilon

      eqn_nextstep.q_vec[i] = unperturbed_q_vec[i]

    end

    #--------------------------
    # emulates cnJac
    scale!(jac, -0.5*h)
    for i = 1:mesh.numDof
      jac[i,i] += 1
    end

    #--------------------------
    # emulates cnRhs
    #   what this is doing:
    #   u_(n+1) - 0.5*dt* (del dot G_(n+1)) - u_n - 0.5*dt* (del dot G_n)
    #=
    physics_func(mesh, sbp, eqn, opts, t)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)

    rhs_vec = zeros(eqn.q_vec)

    for i = 1:mesh.numDof
      rhs_vec[i] = eqn_nextstep.q_vec[i] - h*0.5*eqn_nextstep.res_vec[i] - eqn.q_vec[i] - h*0.5*eqn.res_vec[i]
    end
    =#
    physics_func(mesh, sbp, eqn, opts, t)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    applyMassMatrixInv(mesh, eqn, eqn.res_vec)
    current_t_step_contribution = zeros(eqn.q_vec)
    for i = 1:mesh.numDof
      current_t_step_contribution[i] = - eqn.q_vec[i] - h*0.5*eqn.res_vec[i]
    end

    # Test for 3D Minv results
    # this works!
#     res_vec_control = deepcopy(eqn.res_vec)
#     res_vec_test = deepcopy(eqn.res_vec)
#     res_control = deepcopy(eqn.res)
#     res_test = deepcopy(eqn.res)
#     applyMassMatrixInv3D(mesh, sbp, eqn, res_test)
#     assembleSolution(mesh, sbp, eqn, opts, res_test, res_vec_test)
#     println("=+=+=+ norm of diff btwn res_vec_test & res_vec_control: ", norm(res_vec_test - res_vec_control))

    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
    applyMassMatrixInv(mesh, eqn_nextstep, eqn_nextstep.res_vec)
    next_t_step_contribution = zeros(eqn.q_vec)
    for i = 1:mesh.numDof
      next_t_step_contribution[i] = eqn_nextstep.q_vec[i] - h*0.5*eqn_nextstep.res_vec[i] 
    end

    rhs_vec = zeros(eqn.q_vec)

    for i = 1:mesh.numDof
      rhs_vec[i] = current_t_step_contribution[i] + next_t_step_contribution[i]
    end

    # TODO: check these args
    rhs_norm = calcNorm(eqn, rhs_vec, strongres=true)

    #--------------------------
    # start of actual Newton
    neg_rhs = scale(rhs_vec, -1.0)

    fill!(delta_q_vec, 0.0)
    delta_q_vec = jac\neg_rhs
    fill!(jac, 0.0)

    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] += delta_q_vec[i]
    end

    rhs_norm_tol = 1e-6
    if rhs_norm < rhs_norm_tol
      println("=== cnNewton converged with rhs_norm under $rhs_norm_tol -- newton iters: $newton_i ===")
      return nothing
    end

  end   # end of newton iterations

  println("=== cnNewton did not converge ===")
  return nothing


end


#TODO: why not use the versions in Utils?
# TODO: comment here
function applyMassMatrixInv(mesh, eqn, vec)

  for k = 1:mesh.numDof
    vec[k] = eqn.Minv[k] * vec[k]
  end

  return vec
end

# TODO: comment here
function applyMassMatrixInv3D(mesh, sbp, eqn, arr)

  for i = 1:mesh.numEl
    for j = 1:sbp.numnodes
      for k = 1:mesh.numDofPerNode
        arr[k, j, i] = eqn.Minv3D[k, j, i] * arr[k, j, i]
      end
    end
  end

  return arr
end


