# API functions for physics modules that are supplied automatically (and
# therefore do not need to be extended by each physics)

include("interface2_helper.jl")

"""
  This function returns the fully initialized objects needed to solve an
  equations.

  **Inputs**

   * fname: input file name

  **Outputs**

   * mesh: an AbstractMesh of some kind (depending on input file)
   * sbp: an SBP operator of some kind (depending on input file)
   * eqn: an AbstractSolutionData of some kind (depending on input file)
   * opts: options dictionary
   * pmesh: mesh for computing preconditioning, may be the same object as `mesh`
"""
function createObjects(fname::AbstractString)

  # turn input file into dictionary
  #TODO: finish this
  opts = read_input(fname)
  return createObjects(opts)
end


"""
  Method for creating objects from an options dictionary

  **Inputs**

   * opts: the options dictionary

  **Outputs**

    see other method
"""
function createObjects(opts::Dict)

  # retrieve module specific createObjects
  mod, _createObjects, _checkOptions = retrieve_physics(opts["physics"])

  return _createObjects(opts)
end


"""
  This method constructs a new equation object for the given mesh, sbp
  and opts.

  This is useful for solving on a submesh.

  **Inputs**

   * mesh: an `AbstractMesh`
   * sbp: an `AbstractOperator`
   * opts: options dictionary

  **Outputs**

   same as other method
"""
function createObjects(mesh::AbstractMesh, sbp::AbstractOperator, opts::Dict)

  # make sure this options dictionary has all the default values
  read_input(opts)
  mod, _createObjects, _checkOptions = retrieve_physics(opts["physics"])

  return _createObjects(mesh, sbp, opts)
end


"""
  Additional methods for `solvePDE`, allows solving the PDE starting
  with either a file name or options dictionary.

  See the other method for details

  **Inputs**

   * opts: either a file name to load an options dictionary from or the
           options dictionary itself

  **Outputs**

   same as other method
"""
function solvePDE(opts::Union{Dict, AbstractString})


  read_input(opts)
  mesh, sbp, eqn, opts, pmesh = createObjects(opts)
  ls = solvePDE(mesh, sbp, eqn, opts, pmesh)
#  @assert ls != nothing  # catch case where physics modules forgets to return
                         # the LinearSolver
  return ls
end


"""
  Creates a functional object, using the data described in the options
  dictionary

  **Arguments**

  * `mesh` : Abstract PUMI mesh
  * `sbp`  : Summation-by-parts operator
  * `eqn`  : AbstractSolutionData
  * `opts` : Options dictionary
  * `functional_number` : which functional (of the functionals described
                          in the options dictionary) to create


  **Outputs**

   same as other method

"""
function createFunctional(mesh::AbstractMesh, sbp::AbstractOperator,
                                    eqn::AbstractSolutionData, opts,
                                    functional_number::Int=1)

  dict_val = string("functional_name", functional_number)
  key = string("functional_bcs", functional_number)
  functional_name = opts[dict_val]
  functional_bcs = opts[key]

  return createFunctional(mesh, sbp, eqn, opts, functional_name, functional_bcs)
end


"""
  The functional is evaluated at the state in eqn.q_vec.  Parallel
  communication will be started if required by the functional.

  **Inputs**

   * `mesh` :  Abstract mesh object
   * `sbp`  : Summation-By-Parts operator
   * `eqn`  : AbstractSolutionData object
   * `opts` : Options dictionary
   * `func` : `AbstractFunctional` to be evaluated

  **Outputs**

   * val: the functional value

  **Keyword Arguments**

   * start_comm: if true, parallel communication will be started if required
                 by the functional.  If false, it will not be started, even
                 if required by the functional.  This argument gives the caller
                 a way to improve performance if parallel communication has
                 already been done.  As a side effect, `eqn.q_vec` will not
                 be unpacked into `eqn.q` in this case (because parallel
                 communication is based on `eqn.q`, so they must already be
                 consistent).  Default true

"""
function evalFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::AbstractSolutionData{Tsol}, opts,
            func::AbstractFunctional; start_comm=true) where {Tmsh, Tsol}

  start_comm_q = startCommunicationFunctional(mesh, sbp, eqn, opts,
                                              func, start_comm)

  # evaluate functional
  setupFunctional(mesh, sbp, eqn, opts, func)
  val = _evalFunctional(mesh, sbp, eqn, opts, func)

  # verify implementation finished communication
  if start_comm_q
    assertReceivesWaited(eqn.shared_data)
  end

  return val
end


"""
  Performs reverse-mode differentiation of a functional with respect to the
  metrics.  The functional is evaluated in the state at `eqn.q_vec`.
  The `_bar` fields of the mesh are updated with the result
  It is the callers responsiblity to zero out these fields beforehand, if
  required.  Mesh implementation should provide a function `zeroBarArrays`
  to do this.  The fields of the mesh that are updated are:

   * `dxidx_bar`
   * `jac_bar`
   * `nrm_bndry_bar`
   * `nrm_face_bar`
   * `coords_bndry_bar`
   * `coords_bar`
   * `nrm_sharedface_bar`

  **Inputs**

   * mesh: `_bar` fields are updated
   * sbp
   * eqn
   * opts
   * func: [`AbstractFunctional`](@ref) to compute the derivative of

  **Keyword Arguments**

   * start_comm: see [`evalFunctional`](@ref)

"""
function evalFunctionalDeriv_m(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           func::AbstractFunctional, val_bar::Number=1;
                           start_comm=true
                           ) where {Tmsh, Tsol}

  start_comm_q = startCommunicationFunctional(mesh, sbp, eqn, opts, func,
                                              start_comm)

  # evaluate the functional derivative
  setupFunctional(mesh, sbp, eqn, opts, func)
  _evalFunctionalDeriv_m(mesh, sbp, eqn, opts, func, val_bar)

  # verify implementation finished communication
  if start_comm_q
    assertReceivesWaited(eqn.shared_data)
  end

  return nothing
end


"""
  Computes the derivative of the functional with respect to the solution.
  The derivative will be evaluated at the state in `eqn.q_vec`.

  **Inputs**
   
   * mesh
   * sbp
   * eqn
   * opts
   * func: the [`AbstractFunctional`](@ref) to compute the derivative
                     of

  **Inputs/Outputs**

   * func_deriv_arr: array, same shape as `eqn.q` to overwrite with the
                     derivative of the functional wrt the 3D array form
                     of the solution

  **Keyword Arguments**

   * start_comm: see [`evalFunctional`](@ref)
"""
function evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           func::AbstractFunctional,
                           func_deriv_arr::Abstract3DArray;
                           start_comm=true) where {Tmsh, Tsol}
 
  @assert size(func_deriv_arr, 1) == mesh.numDofPerNode
  @assert size(func_deriv_arr, 2) == mesh.numNodesPerElement
  @assert size(func_deriv_arr, 3) == mesh.numEl

  start_comm_q = startCommunicationFunctional(mesh, sbp, eqn, opts, func,
                                              start_comm)

  # evaluate functional derivative
  setupFunctional(mesh, sbp, eqn, opts, func)
  _evalFunctionalDeriv_q(mesh, sbp, eqn, opts, func, func_deriv_arr)

  # verify the implementation finished communication
  if start_comm_q
    assertReceivesWaited(eqn.shared_data)
  end

  return nothing
end


"""
  This method allows the user to supply a 1D vector as the seed values for
  reverse mode with respect to the metrics.  The residual will be evaluated
  at the state in `eqn.q_vec`, and the `_bar` fields of the mesh will be
  updated with the result. The fields of the mesh that will be
  updated are:

   * `dxidx_bar`
   * `jac_bar`
   * `nrm_bndry_bar`
   * `nrm_face_bar`
   * `coords_bndry_bar`
   * `coords_bar`
   * `nrm_sharedface_bar`
   * `remote_metrics_bar`

  It is the caller's responsibility to zero out these arrays before calling
  this function, if required.  Mesh implementations should provide a
  `zeroBarArrays` function to do this.

  Note: the mesh and equation objects must have the `need_adjoint` option set
  to true in the options dictionary used to create them.


  **Inputs**

   * mesh: the `_bar` fields of the mesh are updated
   * sbp
   * eqn: `eqn.res_bar` is overwritten by unpacking `input_array`
   * opts
   * input_array: 1D array that is the seed vector for reverse mode
   * t: optional argument for time value 

  **Keyword Arguments**

   * start_comm: if true, this function will start parallel communication for
                 the solution,
                 otherwise parallel communication will not be started.
                 As a side-effect, if this argument is false then the
                 derivative will be evaluated at the state in `eqn.q` not
                 `eqn.q_vec` (because the parallel communication was done using
                 `eqn.q`.
"""
function evalResidual_revm(mesh::AbstractMesh, sbp::AbstractOperator,
                     eqn::AbstractSolutionData, opts::Dict,
                     input_array::AbstractArray{T, 1}, t::Number=0.0;
                     start_comm=true) where {T}

  array1DTo3D(mesh, sbp, eqn, opts, input_array, eqn.res_bar)
  # do parallel communication
  start_comm_q = start_comm || getSharedData(eqn.shared_data) != PARALLEL_DATA_ELEMENT
  if start_comm
    array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
    eqn.params.time.t_send += @elapsed if eqn.commsize > 1 && start_comm_q
      setParallelData(eqn.shared_data, PARALLEL_DATA_ELEMENT)
      startSolutionExchange(mesh, sbp, eqn, opts)
    end
  end

  evalResidual_revm(mesh, sbp, eqn, opts, t)

  # verify evalResidual_revm finished communication
  if eqn.commsize > 1
    if start_comm_q
      assertReceivesWaited(eqn.shared_data)
    end
  end

  return nothing
end


"""
  This function allows the user to supply 1D input and output vectors for
  doing reverse mode calculations with respect to the solution.  This function
  computes the matrix-free product psi^T dR/dq, where `psi` is the input
  vector.  The product is computed at the solution in `eqn.q_vec`.

  Note: the mesh and equation objects must have the `need_adjoint` option set
  to true in the options dictionary used to create them.


  **Inputs**

   * mesh
   * sbp
   * eqn: `eqn.q_bar` and `eqn.res_bar` are overwritten
   * opts
   * input_array: vector, same shape as `eqn.res_vec` containing the seed
                  values for the reverse mode calculation
  * t: optional argument for time value

  **Inputs/Outputs**

   * output_array: array, same shape as `eqn.q_vec` to be updated with the
                   result

  **Keyword Arguments**

   * zero_output: is true, overwrite the output array, if false, sum into it,
                  default true
   * start_comm: see [`evalResidual_revm`](@ref).  Parallel communication
                 is always started for `eqn.res_bar` (because the `input_array`
                  needs to be scattered).
"""
function evalResidual_revq(mesh::AbstractMesh, sbp::AbstractOperator,
                     eqn::AbstractSolutionData,
                     opts::Dict, input_array::AbstractVector,
                     output_array::AbstractVector, t::Number=0.0;
                     zero_output=true, start_comm=true)


  # start parallel communication
  start_comm_q = start_comm || getParallelData(eqn.shared_data) != PARALLEL_DATA_ELEMENT
  if start_comm
    array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
    if start_comm_q
      setParallelData(eqn.shared_data, PARALLEL_DATA_ELEMENT)
    end
  end

  fill!(eqn.q_bar, 0)
  array1DTo3D(mesh, sbp, eqn, opts, input_array, eqn.res_bar)
  setParallelData(eqn.shared_data_bar, PARALLEL_DATA_ELEMENT)
  startSolutionExchange_rev2(mesh, sbp, eqn, opts, send_q=start_comm_q)

  # do reverse mode calculation
  evalResidual_revq(mesh, sbp, eqn, opts, t)

  # accumulate into output array
  array3DTo1D(mesh, sbp, eqn, opts, eqn.q_bar, output_array, zero_resvec=zero_output)

  # verify evalResidual_revq finished communication
  if eqn.commsize > 1
    if start_comm_q
      assertReceivesWaited(eqn.shared_data)
    end
    assertReceivesWaited(eqn.shared_data_bar)
  end


  return nothing
end

"""
  This function evaluates a matrix-free Jacobian vector product for any
  physics.  The solution variables must be complex for this to work.

  The derivative is evaluated at the state in `eqn.q_vec`.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * input_array: vector to multiply against
   * t: current time value

  **Inputs/Outputs**

   * output_array: vector to store the result in

  **Keyword Arguments**

   * zero_output: if true, `output_array` is overwritten, otherwise it is added
                  to, default true.
"""
function evaldRdqProduct(mesh::AbstractMesh, sbp::AbstractOperator,
                        eqn::AbstractSolutionData,
                        opts::Dict, input_array::AbstractVector, 
                        output_array::AbstractVector,
                        t::Number=0.0; zero_output=true)

  @assert length(input_array) == length(output_array)
  @assert length(input_array) == mesh.numDof

  h = 1e-20
  pert = Complex128(0, h)
  for i=1:mesh.numDof
    eqn.q_vec[i] += pert*input_array[i]
  end

  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # start parallel communication if needed
  time = eqn.params.time
  time.t_send += @elapsed if opts["parallel_type"] == 2
    startSolutionExchange(mesh, sbp, eqn, opts)
  end

  evalResidual(mesh, sbp, eqn, opts, t)

  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  if zero_output
    @simd for i=1:mesh.numDof
      output_array[i] = imag(eqn.res_vec[i])/h
    end
  else
    @simd for i=1:mesh.numDof
      output_array[i] += imag(eqn.res_vec[i])/h
    end
  end

  removeComplex(eqn)

  if eqn.commsize > 1
    assertReceivesWaited(eqn.shared_data)
  end

  return nothing
end
 
"""
  Returns an [`StandardLinarSolver`](@ref) that can be used to solve the system
  `dR/dq * x = b`, where `R` is the spatial residual and `q` is the solution.
  Users should call `free` on this object when they are finished with it.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * jac_type: an integer specifying what kind of Jacobian to create.  See
               the documentation for the input option `jac_type`.  Default
               value `opts["jac_type"]`

  **Outputs**

   * ls: a `StandardLinaerSolver`
"""
function createLinearSolver
  # this method is implemented with the Newton linear operators
end


"""
  User-face function to get the sparsity pattern of the Jacobian.

  **Inputs**

   * mesh: `AbstractMesh` object
   * sbp:  `AbstractOperator` object
   * eqn: `AbstractSolutionData` object, full initialized.  Each physics
          module should specialize this argument
   * opts: options dictonary.

  **Outputs**

   * disc_type: one of the enums

  **Options Keys**

   * if `preallocate_jacobian_coloring` is true, this function returns
     `COLORING`.  Otherwise is returns whatever the physics module specified.
"""
function getSparsityPattern(mesh::AbstractMesh, sbp::AbstractOperator,
                            eqn::AbstractSolutionData, opts)

  
  if opts["preallocate_jacobian_coloring"]
    disctype = COLORING
  else
    disctype = _getSparsityPattern(mesh, sbp, eqn, opts)
  end

  return disctype
end
# functions declared here but defined elsewhere
function createSBPOperator
end


"""
  Construct objects of degree p + 1

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts

  **Outputs**

   * newmesh: new mesh, using the same underlying apf::Mesh as `mesh`
   * newsbp: degree p + 1 SBP operator
   * neweqn: new `eqn` object with solution interpolated from `eqn.q_vec`
   * newopts: new options dictionary (copy of `opts`)
"""
function createEnrichedObjects(mesh::T, sbp, eqn, opts) where {T <: AbstractMesh}

  opts_h = deepcopy(opts)
  opts_h["order"] += 1

  if opts_h["use_staggered_grid"]
    error("staggered grids not supported for error estimation")
  end

  sbp_h, sbpface_h, shape_type, topo = createSBPOperator(opts_h, Float64, mesh.comm)

  # call constructor for mesh
  mesh_h = copy_mesh(mesh, sbp_h, opts_h, sbpface_h)

  # construct eqn object
  mesh_h, sbp_h, eqn_h, opts_h, pmesh_h = createObjects(mesh_h, sbp_h, opts_h)

  # interpolate field
  interpField(mesh, sbp, eqn.q_vec, mesh_h, sbp_h, eqn_h.q_vec)
  array1DTo3D(mesh_h, sbp_h, eqn_h, opts_h, eqn_h.q_vec, eqn_h.q)

  copyParameters(eqn, eqn_h)

  return mesh_h, sbp_h, eqn_h, opts_h
end


