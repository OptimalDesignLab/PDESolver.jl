# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson_ds, cnResidual

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

   This function supported jacobian/preconditioner freezing with the prefix
   "CN", with the default setting to never recalculate either.  newtonInner
   will use its recalculation policy to recalculate the PC and jacobian.
"""->
function crank_nicolson_ds(f::Function, delta_t::AbstractFloat, t_max::AbstractFloat,
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
  t_steps = round(Int, t_max/delta_t)

  # eqn_nextstep = deepcopy(eqn)
  eqn_nextstep = eqn_deepcopy(mesh, sbp, eqn, opts)

  # TODO: copyForMultistage! does not give correct values.
  #     deepcopy works for now, but uses more memory than copyForMultistage!, if it worked
  # eqn_nextstep = copyForMultistage!(eqn)
  eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode,
                           mesh.numNodesPerElement, mesh.numEl)
  eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode,
                             mesh.numNodesPerElement, mesh.numEl)

  @debug1 println("============ In CN ============")

   # setup all the checkpointing related data
  chkpointer, chkpointdata, skip_checkpoint = explicit_checkpoint_setup(opts, myrank)
  istart = chkpointdata.i

  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop

  pc, lo = getCNPCandLO(mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm, opts)
  newton_data, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, ls)
  newton_data.itermax = 30
  recalc_policy = getRecalculationPolicy(opts, "CN")

  #-------------------------------------------------------------------------------
  # DS of Cd wrt M : setup
  if opts["perturb_Ma"]
    term23 = zero(eqn.params.Ma)
    Ma_pert_mag = opts["perturb_Ma_magnitude"]
    Ma_pert = complex(0, Ma_pert_mag)

    # dJdu was formerly named 'term2'
    dJdu = zeros(eqn.q)
    dJdu_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
    @mpi_master println("Direct sensitivity setup done.")
  end   # end if opts["perturb_Ma"]

  #-------------------------------------------------------------------------------
  # capture DS at the IC
  #   v is the direct sensitivity, du/dM
  #   Ma has been perturbed during setup, in types.jl when eqn.params is initialized
  if opts["write_drag"]
    objective = createFunctional(mesh, sbp, eqn, opts, 1)   # 1 is the functional num
    drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))

    # note about drag writing: file_dict populated and file opened in src/solver/euler/types.jl
    @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
    @mpi_master println(f_drag, 1, " ", drag)     # 1 here is the time step index
    @mpi_master flush(f_drag)
  end
  if opts["write_L2vnorm"]
    @mpi_master f_L2vnorm = eqn.file_dict[opts["write_L2vnorm_fname"]]
  end

  if opts["perturb_Ma"]

    # this is the IC, so it gets the first time step's quad_weight
    # i = 1     # note that timestep loop below starts at i=2
    i = istart - 1    # this should work for restarting
    finaliter = calcFinalIter(t_steps, itermax)     # requires that itermax is set by opts.
                                                    # TODO: make it work for terminating on t_steps
    quad_weight = calcQuadWeight(i, delta_t, finaliter)

    #-----------------------------------------------------------------------------
    # allocation of objects for stabilization routine
    if opts["stabilize_v"]
      # initialize DiagJac. block size = num Dof per El. num blocks = num El.
      stab_A = DiagJac(Complex128, mesh.numDofPerNode*mesh.numNodesPerElement, mesh.numEl)    
      stab_assembler = AssembleDiagJacData(mesh, sbp, eqn, opts, stab_A)

      Bv = zeros(Complex128, length(q_vec))     # vector for holding product of B * v = B * du/dM

      tmp_imag = zeros(Float64, length(eqn.q_vec))      # TODO: comment this
      idqimag_vec = zeros(Bv)                           # TODO: comment this
      @mpi_master f_stabilize_v = open("stabilize_v_updates.dat", "w")      # TODO: buffered IO
    end

    # Note: no stabilization of q_vec at the IC
    v_vec = zeros(q_vec)    # direct sensitivity vector
    for v_ix = 1:length(v_vec)
      v_vec[v_ix] = imag(q_vec[v_ix])/Ma_pert_mag
    end

    # dJdu already allocated above
    evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dJdu)    # dJdu is func_deriv_arr

    # do the dot product of the two terms, and save
    #   this dot product is: dJdu*dudM

    # dJdu_vec already allocated above
    array3DTo1D(mesh, sbp, eqn, opts, dJdu, dJdu_vec)

    for v_ix = 1:length(v_vec)
      # this accumulation occurs across all dofs and all time steps
      term23 += quad_weight * term2_vec[v_ix] * v_vec[v_ix]
    end

    if opts["write_L2vnorm"]
      L2_v_norm = calcNorm(eqn, v_vec)
      @mpi_master println(f_L2vnorm, i, "  ", L2_v_norm)
    end

  end   # end if opts["perturb_Ma"]


  #-------------------------------------------------------------------------------
  # Main timestepping loop
  #   this loop is 2:(t_steps+1) when not restarting
  @mpi_master println(BSTDOUT, "---- Ma @ CN start: ", eqn.params.Ma, " ----")
  for i = istart:(t_steps + 1)

    if opts["perturb_Ma"]
      quad_weight = calcQuadWeight(i, delta_t, finaliter)
    end

    t = (i-2)*delta_t
    @mpi_master println(BSTDOUT, "\ni = ", i, ", t = ", t)
    @debug1 println(eqn.params.f, "====== CN: at the top of time-stepping loop, t = $t, i = $i")
    @debug1 flush(eqn.params.f)

    if use_checkpointing && i % chkpoint_freq == 0
      if skip_checkpoint
        skip_checkpoint = false
      else
        @mpi_master println(BSTDOUT, "Saving checkpoint at timestep ", i)
        skip_checkpoint = false
        # save all needed variables to the chkpointdata
        chkpointdata.i = i

        if countFreeCheckpoints(chkpointer) == 0
          freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint
        end

        # save the checkpoint
        saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpointdata)
      end   # end of if skip_checkpoint check
    end   # end of if use_checkpointing check

#=
    #----------------------------
    # zero out Jac
    #   this works for both PETSc and Julia matrices.
    #   when jac is a Julia matrix, this is effectively wrapping: fill!(jac, 0.0)
    if jac_type != 4
      MatZeroEntries(jac)
    end
=#

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

    t_nextstep = t + delta_t

    #------------------------------------------------------------------------------
    # stabilize q_vec: needs to be before q_vec update
    if opts["stabilize_v"]
      # q_vec is obtained from eqn.q_vec
      # TODO: check treal
      calcStabilizedQUpdate!(mesh, sbp, eqn, opts, stab_A, stab_assembler, treal, Bv, tmp_imag)
    end

    # stabilize q_vec (this only affects the imaginary part of q_vec)
    # TODO TODO - figure out why I did a separate time step for the imaginary part





    # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
    if opts["cleansheet_CN_newton"]
      cnNewton(mesh, sbp, opts, delta_t, f, eqn, eqn_nextstep, t_nextstep)
    else

      ctx_residual = (f, eqn, delta_t, newton_data)

      # recalculate PC and LO if needed
      doRecalculation(recalc_policy, i,
                    ls, mesh, sbp, eqn_nextstep, opts, ctx_residual, t_nextstep)

      newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, ls, 
                  rhs_vec, ctx_residual, t_nextstep)
    end

    # do the callback using the current eqn object at time t
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts, BSTDOUT)

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

    #------------------------------------------------------------------------------
    # write drag at end of time step
    if opts["write_drag"]
      drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))
      @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
      @mpi_master println(f_drag, i, "  ", drag)
      @mpi_master if (i % opts["output_freq"]) == 0
        flush(f_drag)
      end
    end

    #------------------------------------------------------------------------------
    # direct sensitivity of Cd wrt M : calculation each time step
    if opts["perturb_Ma"]

      # v is the direct sensitivity, du/dM
      # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
      # TODO: check q_vec or eqn.q_vec
      for v_ix = 1:length(v_vec)
        v_vec[v_ix] = imag(q_vec[v_ix])/Ma_pert_mag
      end

      # dJdu: partial deriv of functional J wrt state (dCd/du)
      fill!(dJdu, 0.0)      # initialized before timestepping loop
      evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dJdu)    # dJdu is func_deriv_arr

      # do the dot product of the two terms, and save
      fill!(dJdu_vec, 0.0)
      array3DTo1D(mesh, sbp, eqn, opts, dJdu, dJdu_vec)   # dJdu -> dJdu_vec

      for v_ix = 1:length(v_vec)
        # this accumulation occurs across all dofs and all time steps.
        term23 += quad_weight * dJdu_vec[v_ix] * v_vec[v_ix]
      end

      # calculate 'energy' of v_vec and write to file
      if opts["write_L2vnorm"]
        L2_v_norm = calcNorm(eqn, v_vec)
        @mpi_master println(f_L2vnorm, i, "  ", L2_v_norm)
        if (i % opts["output_freq"]) == 0
          @mpi_master flush(f_L2vnorm)
        end
      end
    end       # end if opts["perturb_Ma"]

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
    array1DTo3D(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q_vec, eqn_nextstep.q)

#    t = t_nextstep
    flush(BSTDOUT)

  end   # end of t step loop

  # final time update
  t += delta_t

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copyForMultistage!(dest, src)
  copyForMultistage!(eqn, eqn_nextstep)


  #TODO: return the NewtonData?
  free(newton_data)

  @mpi_master begin
    println("------------------------------------------")
    println("   CN: final time step reached. t = $t")
    println("------------------------------------------")
  end

  if opts["perturb_Ma"]
    @mpi_master if opts["stabilize_v"]
      close(f_stabilize_v)
    end

    @mpi_master close(f_drag)

    @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)
    @mpi_master println(" Ma_pert: ", Ma_pert)
    eqn.params.Ma -= Ma_pert      # removing complex perturbation from Ma now
    @mpi_master println(" pert removed from Ma")
    @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)

    finaliter = calcFinalIter(t_steps, itermax)
    Cd, dCddM = calcDragTimeAverage(mesh, sbp, eqn, opts, delta_t, finaliter)     # will use eqn.params.Ma
    term23 = term23 * 1.0/t     # final step of time average: divide by total time
    global_term23 = MPI.Allreduce(term23, MPI.SUM, mesh.comm)
    total_dCddM = dCddM + global_term23

    # Cd calculations
    @mpi_master begin
      f_total_dCddM = open("total_dCddM.dat", "w")
      println(f_total_dCddM, " dCd/dM: ", dCddM)
      println(f_total_dCddM, " global_term23: ", global_term23)
      println(f_total_dCddM, " total dCd/dM: ", total_dCddM)
      flush(f_total_dCddM)
      close(f_total_dCddM)
      println(" dCd/dM: ", dCddM)
      println(" global_term23: ", global_term23)
      println(" total dCd/dM: ", total_dCddM)
    end

  end   # end if opts["perturb_Ma"]

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

mutable struct CNMatPC <: AbstractPetscMatPC
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
mutable struct CNVolumePC <: AbstractPetscMatFreePC
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

mutable struct CNDenseLO <: AbstractDenseLO
  lo_inner::NewtonDenseLO
end

function CNDenseLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonDenseLO(pc, mesh, sbp, eqn, opts)

  return CNDenseLO(lo_inner)
end

mutable struct CNSparseDirectLO <: AbstractSparseDirectLO
  lo_inner::NewtonSparseDirectLO
end

function CNSparseDirectLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonSparseDirectLO(pc, mesh, sbp, eqn, opts)

  return CNSparseDirectLO(lo_inner)
end


mutable struct CNPetscMatLO <: AbstractPetscMatLO
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
const CNMatLO = Union{CNDenseLO, CNSparseDirectLO, CNPetscMatLO}


function calcLinearOperator(lo::CNMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)

  modifyJacCN(lo, mesh, sbp, eqn, opts, ctx_residual, t)

  return nothing
end


# matrix-free
mutable struct CNPetscMatFreeLO <: AbstractPetscMatFreeLO
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
const CNHasMat = Union{CNMatPC, CNDenseLO, CNSparseDirectLO, CNPetscMatLO}


function calcLinearOperator(lo::CNPetscMatFreeLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)
  
  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)

  setLOCtx(lo, mesh, sbp, eqn, opts, ctx_residual, t)

  return nothing
end

function applyLinearOperator(lo::CNPetscMatFreeLO, mesh::AbstractMesh,
                       sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                       opts::Dict, ctx_residual, t, x::AbstractVector, 
                       b::AbstractVector) where Tsol

  delta_t = ctx[3]
  # the CN residual has the form: I - 0.5*delta_t*jac, where jac is the physics
  # Jacobian.
  # thus a matrix vector product is: v - 0.5*delta_t*jac*v
  applyLinearOperator(lo, mesh, sbp, eqn, opts, ctx_residual, t, x, b)

  scale!(b, -0.5*delta_t)
  @simd for i=1:length(b)
    b[i] += x[i]
  end

  return nothing
end

function applyLinearOperatorTranspose(lo::CNPetscMatFreeLO,
                             mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector) where Tsol

  # see the non-transpose method for an explanation of the math

  delta_t = ctx[3]
  applyLinearOperatorTranspose(lo, mesh, sbp, eqn, opts, ctx_residual, t, x, b)

  scale!(b, -0.5*delta_t)
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
  delta_t = ctx_residual[3]

  assembly_begin(lo2.A, MAT_FINAL_ASSEMBLY)
  assembly_end(lo2.A, MAT_FINAL_ASSEMBLY)

  # scale jac by -delta_t/2
#  scale_factor = h*-0.5
  petsc_scale_factor = PetscScalar(-delta_t*0.5)
  scale!(lo2.A, petsc_scale_factor)

  # add the identity
  diagonal_shift!(lo2.A, 1)

  return nothing
end


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
  delta_t = ctx[3]

  epsilon = opts["epsilon"]
  pert = Complex128(0, epsilon)

  other_pc = pc.other_pc

#  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t_nextstep)
  calcVolumePreconditioner(other_pc, mesh, sbp, eqn_nextstep, opts, pert, 
                           physics_func, t)

  volume_prec = other_pc.vol_prec
  scale!(other_pc.volume_jac, -0.5*delta_t)

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
  delta_t = ctx[3]

  # evalute residual at t_nextstep
  # q_vec -> q
  array1DTo3D(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q_vec)

  # start parallel communication if needed
  time = eqn.params.time
  time.t_send += @elapsed if opts["parallel_type"] == 2 && mesh.npeers > 0
    startSolutionExchange(mesh, sbp, eqn_nextstep, opts)
  end


  physics_func(mesh, sbp, eqn_nextstep, opts, t)
  array3DTo1D(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, 
                   eqn_nextstep.res_vec)

  # evalute residual at t - h
  if opts["parallel_type"] == 2 && mesh.npeers > 0
    startSolutionExchange(mesh, sbp, eqn, opts)
  end
  physics_func(mesh, sbp, eqn, opts, t - delta_t)
  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  # compute rhs

  #   what this is doing:
  #   u_(n+1) - 0.5*dt* (del dot G_(n+1)) - u_n - 0.5*dt* (del dot G_n)
  for i = 1:mesh.numDof

    temp1 = eqn_nextstep.q_vec[i] - 0.5*delta_t*eqn_nextstep.res_vec[i]
    temp2 = eqn.q_vec[i] + 0.5*delta_t*eqn.res_vec[i]

    rhs_vec[i] = temp1 - temp2 

    # NOTE: question: is there a sign problem here? should rhs_vec = -rhs_vec ?
    #     NO. this negative gets applied in newton.jl, where res_0[i] = -res_0[i]

  end

  # calculate the norm
  rhs_0_norm = calcNorm(eqn, rhs_vec, strongres=true)


  return rhs_0_norm

end     # end of function cnRhs

# DS: removed cnNewton

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


