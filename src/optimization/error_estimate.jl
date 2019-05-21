# compute an adjoint-based error estimate

"""
  Options structure for mesh adapation

  **Mesh Adaptation options**

   * strategy: strategy used by [`getTargetSizes`](@ref), default 1
   * free_mesh: if true, the old mesh object will be freed after adaptation
                is complete, Only set this to false if you know what you
                are doing.  Default true.

   * fixed_fraction: a number between 0 and 1 that determines what fraction
                     of the elements are refined for fixed-fraction based
                     strategy.  Default 0.1
  **SolveAdaptive Options**

   * solve_itermax: maximum number of times h-adaptation/solve PDE cycles
                    to do.  Default no limit
   * element_limit: maximum number of elements to generate.  The adaptive
               workflow will abort the first time the mesh has more
               elements than this, so the returned mesh will have
               (possibly many) more elements than this.  Default no
               limit.

  **Strategy**

  The values for the strategy are:

   * 1: attempt to equidistribute the error by refining or coarsening every
        element.  This strategy attempts to reach the error tolerance in a
        single iteration, which can lead to over-refinment
   * 2: refine a fixed fraction of the elements each iteration, reducing
        their size by half.  The elements with the greatest error are selected.
        This strategy does not consider the error tolerance.  Uses the
        `fixed_fraction` field.
"""
mutable struct AdaptOpts
  strategy::Int
  #TODO: wrap more MeshAdapt options
  free_mesh::Bool

  # options used by different strategies
  fixed_fraction::Float64

  # solveAdaptive options
  solve_itermax::Int
  el_limit::Int

  # private options
  free_ls::Bool
end

"""
  Returns object with default options
"""
function AdaptOpts()
  # mesh adaptation options
  strategy = 1
  free_mesh = true

  fixed_fraction = 0.1
  #TODO: need to figure out about load balancing post adaptation

  # solveAdapative options
  solve_itermax = typemax(Int)
  element_limit = typemax(Int)

  # private options
  free_ls = false

  return AdaptOpts(strategy, free_mesh, fixed_fraction, solve_itermax,
                   element_limit, free_ls)
end


"""
  Writes a visualization file with the error estimate.  The output file
  is called "error_estimate".

  **Inputs**

   * mesh
   * el_err: elementwise error estimate
"""
function writeErrorField(mesh, el_err::AbstractVector)

  fshape_ptr = apf.getConstantShapePtr(mesh.dim)
  f = apf.createPackedField(mesh, "error estimate", 1, fshape_ptr)

  val = zeros(Float64, 1)
  for i=1:mesh.numEl
    val[1] = el_err[i]
    apf.setComponents(f, mesh.elements[i], 0, val)
  end

  writeVisFiles(mesh, "error_estimate")

  apf.destroyField(mesh, f)

  return nothing
end


import NonlinearSolvers: NewtonBJacobiPC

"""
  Computes an adjoint-based error estimate and an estimate of each elements's
  contribution to the total error.

  **Inputs**
   * adapt_opts: an [`AdaptOpts`](@ref) object
   * mesh
   * sbp
   * eqn
   * opts
   * func: the functional
   * ls: a LinearSolver used to solve for the adjoint

  **Outputs**

   * err: the estimated error in the functional
   * el_error: the estimated contribution of each element to `err`

  **AdaptOpts Usage**

   * free_ls: if true, free the linear solver after solving the coarse
              space adjoint.  This should free the memory before any
              fine-space data structures are allocated

"""
function calcErrorEstimate(adapt_opts::AdaptOpts, mesh::AbstractMesh,
                           sbp::AbstractOperator{Tsbp},
                           eqn::AbstractSolutionData{Tsol, Tres},
                           opts,
                           func::AbstractFunctional,
                           ls::LinearSolver) where {Tsol, Tres, Tsbp}

  # compute adjoint on coarse space
  psi_H = zeros(Tres, mesh.numDof)
  calcAdjoint(mesh, sbp, eqn, opts, ls, func, psi_H; recalc_jac=true,
              recalc_pc=true, start_comm=true)
  if adapt_opts.free_ls
    free(ls)
    #TODO: gc here to make sure the memory is really freed?
  end

  # construct objects on fine space
  mesh_h, sbp_h, eqn_h, opts_h = createEnrichedObjects(mesh, sbp, eqn, opts)

  
  # interpolate adjoint, solution to fine space
  psi_h = zeros(Tres, mesh_h.numDof)
  interpField(mesh, sbp, psi_H, mesh_h, sbp_h, psi_h)

  # start parallel communication
  setParallelData(eqn_h.shared_data, opts["parallel_data"])
  startSolutionExchange(mesh_h, sbp_h, eqn_h, opts_h)

  # smooth adjoint on fine space
  #TODO: verbose=false, res_tol < 0
  pc = NewtonBJacobiPC(mesh_h, sbp_h, eqn_h, opts_h, itermax=5, verbose=true, res_tol=1e-8) 
  calcPC(pc, mesh_h, sbp_h, eqn_h, opts_h, (evalResidual,), 0.0)

  dJdu_h_arr = zeros(Tsol, mesh_h.numDofPerNode, mesh_h.numNodesPerElement,
                           mesh_h.numEl)
  dJdu_h = zeros(Tsol, mesh_h.numDof)
  evalFunctionalDeriv_q(mesh_h, sbp_h, eqn_h, opts_h, func, dJdu_h_arr,
                        start_comm=false)
  array3DTo1D(mesh_h, sbp_h, eqn_h, opts_h, dJdu_h_arr, dJdu_h)
  scale!(dJdu_h, -1)  # the adjoint problem is dR/dq^T psi = -dJ/du
  applyPCTranspose(pc, mesh_h, sbp_h, eqn_h, opts_h, 0.0, dJdu_h, psi_h)
  free(pc)

  # compute residual on fine space
  # parallel communication was alredy done by calcPC()
  evalResidual(mesh_h, sbp_h, eqn_h, opts_h)
  array3DTo1D(mesh_h, sbp_h ,eqn_h, opts_h, eqn_h.res, eqn_h.res_vec)

  err, el_error = localizeErrorEstimate(mesh_h, sbp_h, eqn_h, opts_h,
                                        eqn_h.res_vec, psi_h)

  finalize(mesh_h)

  if opts["write_error_estimate"]
    writeErrorField(mesh, el_error)
  end

  return err, el_error

end


"""
  Computes the contribution of each element to the total error

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * R_h: fine grid residual evaluated using the interpolated solution
   * psi_h: fine grid adjoint

  **Outputs**

   * err: the estimate of the total error
   * el_error: the estimate of each element's contribution to the total error
"""
function localizeErrorEstimate(mesh::AbstractMesh,
                              sbp::AbstractOperator,
                              eqn::AbstractSolutionData{Tsol, Tres},
                              opts,
                              R_h::AbstractVector,
                              psi_h::AbstractVector) where {Tsol, Tres}
# the arguments are the fine-space objects
  # compute \sum | psi_j R_j |
  err = zero(Float64)
  @simd for i=1:length(R_h)
    err += real(R_h[i]*psi_h[i])
  end
  err = MPI.Allreduce(abs(err), MPI.SUM, eqn.comm)

  el_error = zeros(Float64, mesh.numEl)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        dof = mesh.dofs[k, j, i]
        el_error[i] += abs(R_h[dof]*psi_h[dof])
      end
    end
  end

  return err, el_error
end


"""
  Given the elemnt-local error, this function computes a new size field.

  **Inputs**

   * adapt_opts: an [`AdaptOpts`](@ref) object
   * mesh
   * sbp
   * eqn
   * opts
   * el_error: vector containing the element contribution to the total
               error, length `mesh.numEl`
   * err_target: the target total error (not per-element)

  **Outputs**

   * el_sizes: vector containing the target size of every element.

  **AdaptOpts Usage**

   * strategy: integer indicating what strategy to use for computing
               the sizes, currently:
     * 1: refine everything using theoretical convergence rate of p + 1
"""
function getTargetSizes(adapt_opts::AdaptOpts, mesh::AbstractMesh,
                         sbp::AbstractOperator,
                         eqn::AbstractSolutionData, opts,
                         el_error::AbstractVector,
                         err_target::Number)

  el_sizes = getElementSizes(mesh)

  # currently we refine everything to equi-distribute the error.  This
  # likely leads to over-refinement

  err_target_el = err_target/mesh.numEl
  rate::Int = convert(Int, opts["order"] + 1)  # theoretical convergence rate

  if adapt_opts.strategy == 1
    getTargetSizes_1(adapt_opts, mesh, sbp, eqn, opts, el_error, err_target,
                   el_sizes)
  elseif adapt_opts.strategy == 2
    getTargetSizes_2(adapt_opts, mesh, sbp, eqn, opts, el_error, err_target,
                     el_sizes)
  else
    error("unrecognized sizing strategy: $strategy")
  end

  return el_sizes
end


import PdePumiInterface.adaptMesh


"""
  Wrapper for taking a new size field and adapting the mesh to it.
  Constructs new `AbstractSolutionData` object

  **Inputs**
  
   * adapt_opts: mesh adaptation options
   * mesh
   * sbp
   * eqn
   * opts
   * el_sizes: vector of length `mesh.numEl` containing the new size for
               each element.

  **Outputs**

   * newmesh
   * newsbp
   * neweqn
   * newopts

"""
function adaptMesh(adapt_opts::AdaptOpts, mesh::AbstractMesh,
                   sbp::AbstractOperator,
                   eqn::AbstractSolutionData,
                   opts,
                   el_sizes::AbstractVector)

  #TODO: need options for mesh adaptation
  #      The main thing we want is load balancing

  # perform mesh adaptation (call PumiInterface function)
  newmesh, q_new = adaptMesh(mesh, sbp, opts, el_sizes, eqn.q_vec,
                             free_mesh=adapt_opts.free_mesh)

  # construct new objects
  newopts = deepcopy(opts)
  newmesh, newsbp, neweqn, newopts, pmesh = createObjects(newmesh, sbp,
                                                          newopts)
  copy!(neweqn.q_vec, q_new)

  # return new objects
  return newmesh, newsbp, neweqn, newopts

end


#------------------------------------------------------------------------------
# Interface functions

"""
  Computes an error estimate for the given functional.

  The PDE should have been solved before calling this function.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * func: the `AbstractFunctional`

  **Outputs**

   * err: the estimated functional error
   * el_error: the estimated contribution of each element to `err`
"""
function calcErrorEstimate(mesh::AbstractMesh,
                       sbp::AbstractOperator,
                       eqn::AbstractSolutionData,
                       opts,
                       func::AbstractFunctional)

  # make Newton linear solver
  adapt_opts = AdaptOpts()

  pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm, opts)
  adapt_opts.free_ls = true

  # call other method
  return calcErrorEstimate(adapt_opts, mesh, sbp, eqn, opts, func, ls)
end



"""
  Does one round of h-adaptation, using an adjoint-based error estimate.
  If the returned objects are used to
  solve the PDE, the functional error should be smaller, however it
  may not meet the supplied tolerance.  To continue doing h-adaptation until
  the tolerance is met, see TODO: other function

  On entry, the solution in `eqn.q_vec` should be the solution to the PDE.

  **Inputs**

   * adapt_opts: an [`AdaptOpts`](@ref) object
   * mesh
   * sbp
   * eqn
   * opts
   * func: functional to use for computing adjoint-based error estimate
   * ls: a LinearSolver object.  The linear operator and preconditioner
         will be re-calculated before the linear solve is done
   * err_target: target functional error

  **Outputs**

   * newmesh: new mesh object
   * newsbp: new sbp object
   * neweqn: new eqn object
   * newopts: new options dictonary
   * err: the estimated error at the *prior* to adaptation.  If this
          value is less than the `err_target`, then no adaptation is done
          and the `err` is the estimated error on the current mesh.


  Note that if the original mesh satisfies the error tolerance, the
  input objects are returned and no mesh adaptation is done.


  **Keyword Arguments**

   * free_ls: if true, the linear solver will be freed after use.  This
              can decrease the maximum memory usage of the algorithm.
"""
function doHAdaptation(adapt_opts::AdaptOpts, mesh::AbstractMesh,
                       sbp::AbstractOperator,
                       eqn::AbstractSolutionData,
                       opts,
                       func::AbstractFunctional,
                       ls::LinearSolver,
                       err_target::Number;
                       free_ls::Bool=false)
  # do all the steps
  adapt_opts.free_ls = adapt_opts.free_ls || free_ls
  
  err, el_err = calcErrorEstimate(adapt_opts, mesh, sbp, eqn, opts, func, ls)

  if err > err_target  #TODO: some tolerance here?
    el_sizes = getTargetSizes(adapt_opts, mesh, sbp, eqn, opts, el_err,
                              err_target)
    newmesh, newsbp, neweqn, newopts = adaptMesh(adapt_opts, mesh, sbp, eqn,
                                                 opts, el_sizes)
  else
    newmesh = mesh
    newsbp = sbp
    neweqn = eqn
    newopts = opts
  end


  return newmesh, newsbp, neweqn, newopts, err
end


"""
  This method constructs its own linear solver, rather than using a provided
  one.  All other arguments are the same, as well as the outputs.

  **Inputs**

   * adapt_opts
   * mesh
   * sbp
   * eqn
   * opts
   * func
   * err_target
"""
function doHAdaptation(adapt_opts::AdaptOpts, mesh::AbstractMesh,
                       sbp::AbstractOperator,
                       eqn::AbstractSolutionData,
                       opts,
                       func::AbstractFunctional,
                       err_target::Number)
  # make Newton linear solver
  pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm, opts)

  # call other method
  return doHAdaptation(adapt_opts, mesh, sbp, eqn, opts, func, ls,
                       err_target, free_ls=true)

end


"""
  This function performs repeated h-adaptation until the error estimate
  meets the tolerance (or some other criteria causes adaptation to stop
  early.

  **Inputs**

   * adapt_opts: an [`AdaptOpts`](@ref).  In particular, see the
                 `SolveAdapative` section for options that affects this
                 function
   * mesh
   * sbp
   * eqn
   * opts
   * func: the [`AbstractFunctional`](@ref) to estimate the error of
   * err_target: the target error

  **Outputs**

   * mesh: the adapted mesh
   * sbp
   * eqn: the equation object for the adapted mesh
   * opts; the options dictionary.

"""
function solveAdaptive(adapt_opts::AdaptOpts, mesh::AbstractMesh,
                       sbp::AbstractOperator, eqn::AbstractSolutionData,
                       opts,
                       func::AbstractFunctional,
                       err_target::Number)
  # we assume the PDE was *not* previously solved, unlike doHAdaptation


  # keep solving until tolerance met or some other criteria
  ic_orig = opts["IC_name"]
  exit_code = 1  # 1 = itermax, 0 = anything else
  for i=1:adapt_opts.solve_itermax
    
    solvePDE(mesh, sbp, eqn, opts)
    # for all future solves, use the previous solution
    opts["IC_name"] = "ICPassThrough"

    old_numEl = mesh.numEl
    mesh, sbp, eqn, opts, err = doHAdaptation(adapt_opts, mesh, sbp, eqn,
                                              opts, func, err_target)

    println("\n\nOn iteration ", i, ", error estimate = ", err, ", numEl = ", old_numEl, ", avg mesh size = ", calcMeshH(mesh, sbp, eqn, opts))

    if opts["write_adapt_vis"]
      writeVisFiles(mesh, "adapt_$i")
    end

    #TESTING
    J = evalFunctional(mesh, sbp, eqn, opts, func)
    J_exact = -1/1.4
    println(BSTDOUT, "J = ", J, ", J_exact = ", J_exact, ", err = ", J - J_exact)

    numel_global = MPI.Allreduce(mesh.numEl, MPI.SUM, eqn.comm)
    if err < err_target
      println(BSTDOUT, "solveAdaptive met error tolerance with final error estimate of ", err, " < ", err_target, " with ", numel_global, " elements")
      exit_code = 0
      break
    end

    if numel_global > adapt_opts.el_limit
      println(BSTDOUT, "solveAdaptive exiting due to element limit.  Number of elements = ", numel_global, " > ", adapt_opts.element_limit)
      exit_code = 0
      break
    end

  end  # end loop i

  flush(BSTDOUT)

  opts["IC_name"] = ic_orig  # restore the initial condition
  return mesh, sbp, eqn, opts
end





