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
  mod, _createOjbects, _checkOptions = retrieve_physics(opts["physics"])

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

  mesh, sbp, eqn, opts, pmesh = createObjects(opts)
  println("typeof(eqn) = ", typeof(eqn))
  return solvePDE(mesh, sbp, eqn, opts, pmesh)
end
