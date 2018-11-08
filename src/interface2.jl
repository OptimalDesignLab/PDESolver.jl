# API functions for physics modules that are supplied automatically (and
# therefore do not need to be extended by each physics)


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
   * sbp: an `AbstractSBP`
   * opts: options dictionary

  **Outputs**

   same as other method
"""
function createObjects(mesh::AbstractMesh, sbp::AbstractSBP, opts::Dict)

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
  return solvePDE(mesh, sbp, eqn, opts, pmesh)
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
function createFunctional(mesh::AbstractMesh, sbp::AbstractSBP,
                                    eqn::AbstractSolutionData, opts,
                                    functional_number::Int=1)

  dict_val = string("functional_name", functional_number)
  key = string("functional_bcs", functional_number)
  functional_name = opts[dict_val]
  functional_bcs = opts[key]

  return createFunctional(mesh, sbp, eqn, opts, functional_name, functional_bcs)
end

"""
  This method allows the user to supply a 1D vector as the seed values for
  reverse mode with respect to the metrics.  See the other method for details
  of the reverse mode calculation.

  **Inputs**

   * mesh: the `_bar` fields of the mesh are updated
   * sbp
   * eqn: `eqn.res_bar` is overwritten by unpacking `input_array`
   * opts
   * input_array: 1D array that is the seed vector for reverse mode
   * t: optional argument for time value 
"""
function evalResidual_revm(mesh::AbstractMesh, sbp::AbstractSBP,
                     eqn::AbstractSolutionData, opts::Dict,
                     input_array::AbstractArray{T, 1}, t::Number=0.0
                    ) where {T}

  @assert eqn.commsize == 1
  #TODO: do parallel communication on input_array
  array1DTo3D(mesh, sbp, eqn, opts, input_array, eqn.res_bar)

  evalResidual_revm(mesh, sbp, eqn, opts, t)

  return nothing
end

"""
  This function allows the user to supply 1D input and output vectors for
  doing reverse mode calculations with respect to the solution.  See the other
  method for the details of the reverse mode calculation

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
"""
function evalResidual_revq(mesh::AbstractMesh, sbp::AbstractSBP,
                     eqn::AbstractSolutionData,
                     opts::Dict, input_array::AbstractVector,
                     output_array::AbstractVector, t::Number=0.0; zero_output=true)


  @assert mesh.commsize == 1

  #TODO: do parallel communication on input_array
  array1DTo3D(mesh, sbp, eqn, opts, input_array, eqn.res_bar)
  fill!(eqn.q_bar, 0)

  evalResidual_revq(mesh, sbp, eqn, opts, t)

  # accumulate into output array
  array3DTo1D(mesh, sbp, eqn, opts, eqn.q_bar, output_array, zero_resvec=zero_output)

  return nothing
end

