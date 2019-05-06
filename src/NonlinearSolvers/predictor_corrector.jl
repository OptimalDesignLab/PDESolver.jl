# Dissipation based homotopy predictor-corrector globalization for solving
# steady problems using Newton's method
# based on Brown and Zingg, "A monolithic Homotopy Continuation Algorithm 
# with applications to Computational Fluid Dynamics"
# Journal of Computational Physics 321, (2016), 55-75
# specifically, Algorithm 2

mutable struct HomotopyData{Tsol, Tjac}

  physics_func::Function  # the function that will be solved at lambda = 0
  g_func::Function        # the function to be blended with physics_func

  time::Timings
  lambda::Float64 # homotopy parameter
  myrank::Int

  # parameters
  lambda_min::Float64
  lambda_cutoff::Float64
  itermax::Int   # maximum number of homotopy iterations
  res_reltol::Float64  # residual tolerance for the final solve
  res_abstol::Float64  # residual tolerance for the final solve
  krylov_reltol0::Float64  # Krylov tolerance that should be set before
                           # entering Newton's method
  orig_newton_globalize_euler::Bool  # value to reset opts afterwards
  use_pc::Bool  # true if PC is not PCNone
  homotopy_globalize_euler::Bool
  tighten_early::Bool  # tighten nonlinear tolerance when lambda reaches
                       # lambda_cutoff

  # working variables
  iter::Int  # current major iteration
  homotopy_tol::Float64  # the residual tolerance for homotopy iterations
                         # when lambda > 0
  delta_max::Float64  # step size limit
  psi_max::Float64  # max angle between tangent vectors (radians)
  h_max::Float64  # maximum h value
  psi::Float64  # current angle between tangent vectors (radians)
  tan_norm_1::Float64  # previous tangent vector norm
  res_norm::Float64   # current norm of physics (not homotopy) residual
  h::Float64  # step size

  # arrays
  q_vec0::Array{Tsol, 1}
  delta_q::Array{Tsol, 1}
  tan_vec::Array{Tjac, 1}
  tan_vec_1::Array{Tjac, 1}
  dHdLambda_real::Array{Tjac, 1}
  rhs_vec::Array{Tsol, 1}

  # composite objects
  recalc_policy::RecalculationPolicy
  ls::LinearSolver
  newton_data::NewtonData
  fconv::IO

  # homotopy shock capturing
  sensor_easy::AbstractShockSensor  # a shock sensor that is easy to solve
  opts_orig::Dict{String, Bool}     # options dictionary used for resetting
                                    # the real options dictionary to is
                                    # original state

  function HomotopyData{Tsol, Tjac}(mesh, sbp, eqn, opts,
                           physics_func::Function, g_func::Function
                          ) where {Tsol, Tjac}

    println("\nConstructing HomotopyData")
    time = eqn.params.time
    lambda = 1.0  # homotopy parameter
    myrank = mesh.myrank

    # some parameters
    lambda_min = 0.0
    lambda_cutoff = 0.000  # was 0.005
    itermax = opts["itermax"]::Int
    res_reltol=opts["res_reltol"]::Float64
    res_abstol=opts["res_abstol"]::Float64
    krylov_reltol0 = opts["krylov_reltol"]::Float64
    orig_newton_globalize_euler = opts["newton_globalize_euler"]  # reset value
    homotopy_globalize_euler = opts["homotopy_globalize_euler"]
    tighten_early = opts["homotopy_tighten_early"]::Bool

    # counters/loop variables
    iter = 1
    homotopy_tol = 1e-2
    delta_max = 1.0  # step size limit, set to 1 initially,
    #psi_max = 10*pi/180  # angle between tangent limiter
    psi_max = 5*pi/180  # angle between tangent limiter
    h_max = 0.1  # maximum h
    psi = 0.0  # angle between tangent vectors
    tan_norm_1 = 0.0  # previous tangent vector norm
    res_norm = 0.0  # norm of residual (not homotopy)
    h = 0.05  # step size
    lambda -= h  # the jacobian is ill-conditioned at lambda=1, so skip it
    #h = lambda/20

    recalc_policy = getRecalculationPolicy(opts, "homotopy")
    # log file
    if myrank == 0
      fconv = BufferedIO("convergence.dat", "a+")
    else
      fconv = DevNull
    end

    # needed arrays
    q_vec0 = zeros(eqn.q_vec)
    delta_q = zeros(eqn.q_vec)
    tan_vec = zeros(Tjac, length(eqn.q_vec))  # tangent vector
    tan_vec_1 = zeros(tan_vec)  # previous tangent vector
    dHdLambda_real = zeros(Tjac, length(eqn.q_vec)) 

    # stuff for newtonInner
    # because the homotopy function is a pseudo-physics, we can reuse the
    # Newton PC and LO stuff, supplying homotopyPhysics in ctx_residual
    pc, lo = getHomotopyPCandLO(mesh, sbp, eqn, opts)
    use_pc = !(typeof(pc) <: PCNone)

    ls = StandardLinearSolver(pc, lo, eqn.comm, opts)


    # configure NewtonData
    newton_data, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, ls)
    newton_data.itermax = 30

    # homotopy shock capturing
    sensor_easy = createShockSensor(mesh, sbp, eqn, opts,
                                    opts["homotopy_shock_sensor"])
    opts_orig = Dict{String, Bool}()


    obj = new(physics_func, g_func, time, lambda, myrank,
              # parameters
              lambda_min, lambda_cutoff, itermax, res_reltol, res_abstol,
              krylov_reltol0, orig_newton_globalize_euler, use_pc,
              homotopy_globalize_euler, tighten_early,
              # working variables
              iter, homotopy_tol, delta_max, psi_max, h_max, psi,
              tan_norm_1, res_norm, h,
              # arrays
              q_vec0, delta_q, tan_vec, tan_vec_1, dHdLambda_real, rhs_vec,
              # composed objects
              recalc_policy, ls, newton_data, fconv,
              # homotopy shock capturing
              sensor_easy, opts_orig)

    lo.hdata = obj
    setLambda(obj, lambda)
   
    return obj
  end
end

function setLambda(hdata::HomotopyData, lambda::Number)

  hdata.lambda = lambda

  return nothing
end


"""
  This function solves steady problems using a dissipation-based
  predictor-correcor homotopy globalization for Newtons method.

  Inputs:
    physics_func: the function to solve, ie. func(q) = 0  mathematically.
                  func must have the signature func(mesh, sbp, eqn, opts)

    g_func: the function that evalutes G(q), the dissipation.
            If explicit jacobian calculation is used, then
            [`evalHomotopyJacobian`](@ref) must evaluate the jacobian of this 
            function wrt. q
    sbp: an SBP operator
    eqn: a AbstractSolutionData.  On entry, eqn.q_vec must be the
         initial condition.  On exit, eqn.q_vec will be the solution to
         func(q) = 0
    opts: options dictionary


  This function uses Newtons method internally, and supports all the
  different jacobian types and jacobian calculation methods that 
  Newton does.

  On entry, eqn.q_vec should contain the initial guess for q.  On exit
  it will contain the solution for func(q) = 0.

  This function is reentrant.

  **Options Keys**

   * calc_jac_explicit
   * itermax
   * res_reltol
   * res_abstol
   * krylov_reltol
   * homotopy_globalize_euler: activate implicit euler globalization when
                               lambda = 0
   * homotopy_tighten_early: tighten the newton solve tolerance and active
                             implicit Euler globalization when 
                             `lambda < lambda_cutoff`, typically 0.005.
                             This may slow down the Newton solve significantly
   * homotopy_shock_capturing: see below
   * homotopy_shock_sensor: shock sensor for `R_easy` with
                            `homotopy_shock_capturing`

  This function supports jacobian/preconditioner freezing using the
  prefix "newton".  Note that this recalcuation policy only affects
  this function, and not newtonInner, and defaults to never recalculating
  (and letting newtonInner update the jacobian/preconditioner according to
  its recalculationPolicy).

  **Homotopy Shock Capturing**

  Does homotopy from one shock capturing scheme to another.
  The homotopy function is:

  R_euler (1-lambda)*R_hard + lambda*R_easy

  where R_euler is the physics residual, excluding shock capturing, R_hard
  is the shock sensor that is difficult to converge, and R_easy is the shock
  sensor that is easy to converge.  It is recommended that the initial
  condition satisfy R_euler + R_easy.

  In this mode, `g_func` is ignored.

"""
function predictorCorrectorHomotopy(physics_func::Function,
                  g_func::Function,
                  mesh::AbstractMesh{Tmsh}, 
                  sbp::AbstractOperator, 
                  eqn::AbstractSolutionData{Tsol, Tres}, 
                  opts) where {Tsol, Tres, Tmsh}

  #----------------------------------------------------------------------------
  # setup

  Tjac = real(Tres)
  hdata = HomotopyData{Tsol, Tjac}(mesh, sbp, eqn, opts, physics_func, g_func)
  myrank = hdata.myrank

  # because the homotopy function is a pseudo-physics, we can reuse the
  # Newton PC and LO stuff, supplying homotopyPhysics in ctx_residual
  homotopyPhysics = (mesh, sbp, eqn, opts, t) -> _homotopyPhysics(hdata, mesh, sbp, eqn, opts, t)
  rhs_func = physicsRhs
  ctx_residual = (homotopyPhysics,)


  # calculate physics residual
  res_norm = real(physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (physics_func,)))
  res_norm_0 = res_norm

  # print to log file
  @mpi_master println(hdata.fconv, 0, " ", res_norm, " ", 0.0)
  @mpi_master flush(hdata.fconv)

  eqn.majorIterationCallback(0, mesh, sbp, eqn, opts, BSTDOUT)

  #----------------------------------------------------------------------------
  # main loop
  while res_norm > res_norm_0*hdata.res_reltol && res_norm > hdata.res_abstol &&
        hdata.iter < hdata.itermax  # predictor loop

    @unpack hdata h lambda iter
    @mpi_master begin
      println(BSTDOUT, "\npredictor iteration ", iter, ", lambda = ", lambda)
      println(BSTDOUT, "res_norm = ", res_norm)
      println(BSTDOUT, "res_norm/res_norm_0 = ", res_norm/res_norm_0)
    end

    # calculate homotopy residual
    #homotopy_norm = physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (homotopyPhysics,))

    # update tolerances (linear and nonlinar)
    updateTolerances(hdata, opts)

   # calculate the PC and LO if needed
    doRecalculation(hdata.recalc_policy, hdata.iter, hdata.newton_data.ls,
                    mesh, sbp, eqn, opts, ctx_residual, 0.0)

    # do corrector steps
    copy!(hdata.q_vec0, eqn.q_vec)  # save initial q to calculate delta_q later
    newtonInner(hdata.newton_data, mesh, sbp, eqn, opts, rhs_func, hdata.ls, 
                hdata.rhs_vec, ctx_residual)

    # predictor step calculation
    if abs(lambda - hdata.lambda_min) > eps()

      # compute tangent vector
      calcTangentVector(hdata, mesh, sbp, eqn, opts)

      takePredictorStep(hdata, mesh, sbp, eqn, opts)

    end  # end if lambda too large

    # calculate physics residual at new state q
    res_norm = real(physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (physics_func,),))

    # print to log file
    @mpi_master println(hdata.fconv, iter, " ", res_norm, " ", h )
    @mpi_master flush(hdata.fconv)

    eqn.majorIterationCallback(iter, mesh, sbp, eqn, opts, BSTDOUT)

    hdata.iter += 1
  end  # end while loop

  print(BSTDOUT, "\n")

  # inform user of final status
  @mpi_master if hdata.iter >= hdata.itermax
    println(BSTDERR, "Warning: predictor-corrector did not converge in $(hdata.iter) iterations")
  
  elseif res_norm <= hdata.res_abstol
    println(BSTDOUT, "predictor-corrector converged with absolute residual norm $res_norm")
  elseif res_norm/res_norm_0 <= hdata.res_reltol
    tmp = res_norm/res_norm_0
    println(BSTDOUT, "predictor-corrector converged with relative residual norm $tmp")
  end

  # reset options dictionary
  opts["newton_globalize_euler"] = hdata.orig_newton_globalize_euler

  free(hdata.newton_data)

  flush(BSTDOUT)
  flush(BSTDERR)

  return nothing
end


"""
  This function can be wrapped inside a lambda function to make it look like
  evalResidual.
"""
function _homotopyPhysics(hdata::HomotopyData, mesh, sbp, eqn, opts, t)

  # this function is only for use with Newton's method, where parallel
  # communication is started outside the physics
  # q_vec -> q

  if opts["homotopy_shock_capturing"]
    evalHomotopyFunction_sc(hdata, mesh, sbp, eqn, opts, t)

  else
    fill!(eqn.res, 0.0)
    res_homotopy = zeros(eqn.res)


    # calculate physics residual
    # call this function before g_func, to receive parallel communication
    hdata.physics_func(mesh, sbp, eqn, opts, t)

    # calculate homotopy function
    hdata.g_func(mesh, sbp, eqn, opts, res_homotopy)

    # combine (use lambda from outer function)
    lambda_c = 1 - hdata.lambda # complement of lambda
    for i=1:length(eqn.res)
      eqn.res[i] =  lambda_c*eqn.res[i] + hdata.lambda*res_homotopy[i]
    end

    #println("homotopy physics exiting with residual norm ", norm(vec(eqn.res)))
  end

  return nothing
end


"""
  Updates the linear and nonlinear tolerances for the Newton solve.
  The tolerances usually only change then lambda is close to or at 0.
  This also has the effect of resetting the tolerances if Newton changed
  them (for example, doing inexact Newton-Krylov)

  **Inputs**

   * hdata: `HomotopyData`
   * opts: options dictonary

"""
function updateTolerances(hdata::HomotopyData, opts)

  @unpack hdata res_reltol res_abstol homotopy_globalize_euler tighten_early
  @unpack hdata lambda_cutoff lambda_min homotopy_tol krylov_reltol0 lambda
  @unpack hdata newton_data myrank

  # if we have finished traversing the homotopy path, solve the 
  # homotopy problem (= the physics problem because lambda is zero)
  # tightly
  changed_tols = false
  if abs(lambda - lambda_min) <= eps()
    @mpi_master println(BSTDOUT, "tightening homotopy tolerance at lambda = 0")
    changed_tols = true
    homotopy_tol = res_reltol

    reltol = res_reltol*1e-3  # smaller than newton tolerance
    abstol = res_abstol*1e-3  # smaller than newton tolerance

    # enable globalization if required
    opts["newton_globalize_euler"] = opts["homotopy_globalize_euler"]
    newton_data.itermax = opts["itermax"]
  elseif lambda < lambda_cutoff && tighten_early

    # turn on implicit euler and tighten tolerances
    # This helps for shock problems when the algorithm takes very small steps
    # in lambda near the end, because lambda may become too small to prevent
    # Newton from diverging
    @mpi_master println(BSTDOUT, "tightening homotopy tolerance at lambda cutoff")
    changed_tols = true
    homotopy_tol = 1e-4

    reltol = res_reltol*1e-3  # smaller than newton tolerance
    abstol = res_abstol*1e-3  # smaller than newton tolerance

    # enable globalization if required
    opts["newton_globalize_euler"] = opts["homotopy_globalize_euler"]
    newton_data.itermax = opts["itermax"]
  else  # normal iteration
    reltol = krylov_reltol0  # reset in case Newton did inexact Newton-Krylov
    abstol = -1.0
     
  end

  if changed_tols
    @mpi_master begin
      println(BSTDOUT, "setting homotopy tolerance to ", homotopy_tol)
      println(BSTDOUT, "ksp reltol = ", reltol)
      println(BSTDOUT, "ksp abstol = ", abstol)
    end
  end

  # record updated tolerances
  newton_data.res_reltol = homotopy_tol
  hdata.homotopy_tol = homotopy_tol

  # reset tolerances in case newton is doing inexact-NK and thus has changed
  # the tolerances
  setTolerances(newton_data.ls, reltol, abstol, -1, -1)

  return nothing
end
 

"""
  Solve for the tangent vector, overwriting `hdata.tan_vec`

  **Inputs**

   * hdata: `HomotopyData`
   * mesh
   * sbp
   * eqn
   * opts
"""
function calcTangentVector(hdata::HomotopyData, mesh, sbp, eqn, opts)

  @unpack hdata rhs_vec dHdLambda_real tan_vec myrank ls

  # calculate dHdLambda at new q value
  calcdHdLambda(hdata, mesh, sbp, eqn, opts, rhs_vec)
  for i=1:length(rhs_vec)
    dHdLambda_real[i] = real(rhs_vec[i])
  end

  # calculate tangent vector dH/dq * t = dH/dLambda
  @mpi_master println(BSTDOUT, "solving for tangent vector")
  # the accuracy of the tangent vector is not critical, so re-use the
  # Jacobian from the last Newton iteration
  #calcPCandLO(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0)

  tsolve = @elapsed linearSolve(ls, dHdLambda_real, tan_vec)
  eqn.params.time.t_solve += tsolve

  return nothing
end



"""
  This function calculates dH/dLambda, where H is the homotopy function
  calculated by homotopyPhysics.  The differentiation is done analytically.

  **Inputs**

   * hdata: `HomotopyData`
   * mesh
   * sbp
   * eqn: eqn.res and eqn.res_vec are overwritten
   * opts

  Inputs/Outputs
    res_vec: vector to store dH/dLambda in

  Aliasing restrictions: res_vec and eqn.res_vec may not alias
"""
function calcdHdLambda(hdata::HomotopyData, mesh, sbp, eqn, opts, res_vec)

  @unpack hdata physics_func g_func

  # it appears this only gets called after parallel communication is done
  # so no need to start communication here

  if opts["homotopy_shock_capturing"]
    calcdHdLambda_sc(hdata, mesh, sbp, eqn, opts, 0.0, res_vec)
  else
    res_homotopy = zeros(eqn.res)

    # calculate physics residual
    physics_func(mesh, sbp, eqn, opts)
    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)


    # calculate homotopy function
    g_func(mesh, sbp, eqn, opts, res_homotopy)
    array3DTo1D(mesh, sbp, eqn, opts, res_homotopy, res_vec)

    # combine them
    for i=1:length(res_vec)
      res_vec[i] -= eqn.res_vec[i]
    end
  end

  return nothing
end


"""
  This function takes the predictor step.  This requires a few actions:

  * normalize tangent vector (overwrites `hdata.tan_vec`)
  * compute `delta_q` (overwrites `hdata.delta_q` and uses `hdata.q_vec0`
  * applies the various heuristics to update the step size (step length, angle,
    etc.)

"""
function takePredictorStep(hdata::HomotopyData, mesh, sbp, eqn, opts)

  @unpack hdata tan_vec tan_vec_1 tan_norm_1 delta_q q_vec0 h lambda
  @unpack hdata psi_max delta_max h_max lambda_min

  # normalize tangent vector
  # t = [z, -1], where z is a delta_q sized vector, so the norm is
  # sqrt(calcNorm(z)^2 + (1)^2).  In the code tan_vec = z, so what is
  # really happenening is the z component of t is being normalize.  The
  # -1 component will be handled below
  tan_norm = calcNorm(eqn, tan_vec)
  tan_norm = sqrt(tan_norm*tan_norm + 1)
  scale!(tan_vec, 1/tan_norm)

  # compute delta_q
  for i=1:length(eqn.q)
    delta_q[i] = eqn.q_vec[i] - q_vec0[i]
  end
  delta = calcNorm(eqn, delta_q)
  println("L2 norm delta = ", delta)

  psi = psi_max
  if hdata.iter > 1
    # compute phi = acos(tangent_i_1 dot tangent_i)
    # however, use the L2 norm for the z part of the tangent vector
    tan_term = calcL2InnerProduct(eqn, tan_vec_1, tan_vec)
    # now add the normalized -1 components of t_i and t_i_1 to complete
    # the inner product between t_i and t_i_1 = inner_product(z_i, z_i_1) 
    # + -1*-1/(||z_i||*||z_i_1||)
    tan_norm_term = (1/tan_norm)*(1/tan_norm_1)
    arg = tan_term + tan_norm_term
    arg = clamp(arg, -1.0, 1.0)
    psi = acos( arg )
  end

  # calculate step size
  fac = max(real(psi/psi_max), sqrt(delta/delta_max))
  println("psi/psi_max = ", psi/psi_max)
  println("delta/delta_max = ", delta/delta_max)
  println("fac = ", fac)
  h /= fac
  h = min(h, h_max)
  lambda = max(lambda_min, lambda - h)
  if lambda < eps()
    lambda = 0.0
  end
  println("new h = ", h)

  # take predictor step
  for i=1:length(eqn.q_vec)
    eqn.q_vec[i] += h*tan_vec[i]
  end


  # save variables
  setLambda(hdata, lambda)
  copy!(tan_vec_1, tan_vec)
  hdata.tan_norm_1 = tan_norm
  hdata.h = h

  return nothing
end


function getHomotopyPCandLO(mesh, sbp, eqn, opts)

  # get PC
  if opts["jac_type"] <= 2
    pc = PCNone(mesh, sbp, eqn, opts)
  else
    pc = HomotopyMatPC(mesh, sbp, eqn, opts)
  end 

  jactype = opts["jac_type"]
  if jactype == 1
    lo = HomotopyDenseLO(pc, mesh, sbp, eqn, opts)
  elseif jactype == 2
    lo = HomotopySparseDirectLO(pc, mesh, sbp, eqn, opts)
  elseif jactype == 3
    lo = HomotopyPetscMatLO(pc, mesh, sbp, eqn, opts)
  elseif jactype == 4
    lo = HomotopyPetscMatFreeLO(pc, mesh, sbp, eqn, opts)
  end

  return pc, lo
end




# Define PC and LO objects
#------------------------------------------------------------------------------
# PC:

mutable struct HomotopyMatPC <: AbstractPetscMatPC
  pc_inner::NewtonMatPC
  lambda::Float64  # homotopy parameter
end

function HomotopyMatPC(mesh::AbstractMesh, sbp::AbstractOperator,
                    eqn::AbstractSolutionData, opts::Dict)


  pc_inner = NewtonMatPC(mesh, sbp, eqn, opts)
  lambda = 1.0

  return HomotopyMatPC(pc_inner, lambda)
end

function calcPC(pc::HomotopyMatPC, mesh::AbstractMesh, sbp::AbstractOperator,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  # compute the Jacobian of the Newton PC
  calcPC(pc.pc_inner, mesh, sbp, eqn, opts, ctx_residual, t)

  if opts["calc_jac_explicit"]
    A = getBasePC(pc).A
    assembly_begin(A, MAT_FINAL_ASSEMBLY)
    assembly_end(A, MAT_FINAL_ASSEMBLY)

    lambda_c = 1 - lambda # complement of lambda
    scale!(A, lambda_c)

    # compute the homotopy contribution to the Jacobian
    assembler = _AssembleElementData(A, mesh, sbp, eqn, opts)
    evalHomotopyJacobian(mesh, sbp, eqn, opts, assembler, lambda)
  end

  return nothing
end

#------------------------------------------------------------------------------
# LO:

mutable struct HomotopyDenseLO <: AbstractDenseLO
  lo_inner::NewtonDenseLO
  hdata::HomotopyData

  function HomotopyDenseLO(lo_inner::NewtonDenseLO)
    return new(lo_inner)
  end
end

function HomotopyDenseLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractOperator, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonDenseLO(pc, mesh, sbp, eqn, opts)
  return HomotopyDenseLO(lo_inner)
end

mutable struct HomotopySparseDirectLO <: AbstractSparseDirectLO
  lo_inner::NewtonSparseDirectLO
  hdata::HomotopyData

  function HomotopySparseDirectLO(lo_inner::NewtonSparseDirectLO)
    return new(lo_inner)
  end

end

function HomotopySparseDirectLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractOperator, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonSparseDirectLO(pc, mesh, sbp, eqn, opts)

  return HomotopySparseDirectLO(lo_inner)
end

mutable struct HomotopyPetscMatLO <: AbstractPetscMatLO
  lo_inner::NewtonPetscMatLO
  hdata::HomotopyData

  function HomotopyPetscMatLO(lo_inner::NewtonPetscMatLO)
    return new(lo_inner)
  end

end


function HomotopyPetscMatLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    sbp::AbstractOperator, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonPetscMatLO(pc, mesh, sbp, eqn, opts)

  return HomotopyPetscMatLO(lo_inner)
end


mutable struct HomotopyPetscMatFreeLO <: AbstractPetscMatFreeLO
  lo_inner::NewtonPetscMatFreeLO
  hdata::HomotopyData

  function HomotopyPetscMatFreeLO(lo_inner::NewtonPetscMatFreeLO)
    return new(lo_inner)
  end

end

"""
  Homotopy mat-free linear operator constructor

  **Inputs**

   * pc
   * mesh
   * sbp
   * eqn
   * opts
   * rhs_func: rhs_func from [`newtonInner`](@ref)
"""
function HomotopyPetscMatFreeLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    sbp::AbstractOperator, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonPetscMatFreeLO(pc, mesh, sbp, eqn, opts)

  return HomotopyPetscMatFreeLO(lo_inner)
end

#------------------------------------------------------------------------------
# Required functions

"""
  Homotopy matrix-explicit linear operators
"""
const HomotopyMatLO = Union{HomotopyDenseLO, HomotopySparseDirectLO, HomotopyPetscMatLO}


function calcLinearOperator(lo::HomotopyMatLO, mesh::AbstractMesh,
                            sbp::AbstractOperator, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  # ctx_residual[1] = is the *homotopy* function H, so for the non-explict
  # case, calling calcLinearOperator(lo.lo_inner) is fine because it will
  # do coloring
  
  if opts["homotopy_shock_capturing"] && opts["calc_jac_explicit"]
    evalJacobian_homotopy_sc(lo, mesh, sbp, eqn, opts, ctx_residual, t)
  else

    calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)

    if opts["calc_jac_explicit"]
      A = getBaseLO(lo).A
      assembly_begin(A, MAT_FINAL_ASSEMBLY)
      assembly_end(A, MAT_FINAL_ASSEMBLY)

      lambda_c = 1 - lo.hdata.lambda # complement of lambda
      scale!(A, lambda_c)

      # compute the homotopy contribution to the Jacobian
      assembler = _AssembleElementData(A, mesh, sbp, eqn, opts)
      evalHomotopyJacobian(mesh, sbp, eqn, opts, assembler, lo.hdata.lambda)
    end
  end

  return nothing
end


function calcLinearOperator(lo::HomotopyPetscMatFreeLO, mesh::AbstractMesh,
                            sbp::AbstractOperator, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  if opts["homotopy_shock_capturing"]
    error("homotopy shock capturing not supported for matrix-free")
  else
    calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)

    setLOCtx(lo, mesh, sbp, eqn, opts, ctx_residual, 0.0)

    # nothing to do here
  end

  return nothing
end


function applyLinearOperator(lo::HomotopyPetscMatFreeLO, mesh::AbstractMesh,
                       sbp::AbstractOperator, eqn::AbstractSolutionData{Tsol},
                       opts::Dict, ctx_residual, t, x::AbstractVector, 
                       b::AbstractVector) where Tsol

  @assert !(Tsol <: AbstractFloat)  # complex step only!

  if opts["homotopy_shock_capturing"]
    error("homotopy shock capturing not supported for matrix-free")
  else
    # ctx_residual[1] = homotopyPhysics, so this computes both the physics
    # and homotopy contribution
    applyLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t, x, b)
  end

  return nothing
end

function applyLinearOperatorTranspose(lo::HomotopyPetscMatFreeLO, 
                             mesh::AbstractMesh,
                             sbp::AbstractOperator, eqn::AbstractSolutionData{Tsol},
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector) where Tsol

  error("applyLinearOperatorTranspose() not supported by HomotopyPetscMatFreeLO")

end

include("shock_homotopy.jl")
