# Reduced mesh composed only of the elements that have shocks in them and their
# neigbors

#TODO: move this to Utils?

#------------------------------------------------------------------------------
# for DG-type shock capturing, need to know interfaces

"""
  Interface that also stores the index that interface corresponds to in the
  parent mesh.
"""
struct RelativeInterface
  iface::Interface  # needed for SBP
  idx_orig::Int
  #peeridx::Int32  # index of the peer on the *original* mesh
end

"""
  Boundary time that also contains the index of the boundary in the parent mesh.
"""
struct RelativeBoundary
  bndry::Boundary
  idx_orig::Int
  fac::Int8  # factor to multiply normal vector by, 1 or -1
end

function RelativeBoundary(bndry::Boundary, idx_orig::Integer)
  return RelativeBoundary(bndry, idx_orig, 1)
end


"""
  Stores a limited about of information about the mesh defined by the elements
  with shocks in them and their neighbors.

  Definitions:

    shocked
    neighbor

  The way this data structure is meant to be used is:

   1) at the beginning of a simulation, construct this object
   2) every time the residual/Jacobian is evaluated, call `push!` repeatedly
      to add the elements that have shocks in them
   3) call [`completeShockMesh`](@ref) to finish constructing the data
      structure
   4) evaluate the shock capturing scheme


  Note that steps 2-4 will happen repeatedly.  This data structure is optimized
  to be re-used several times, under the assumption that the shock will not
  change position too much.  The function [`reset`](@ref) should be called
  before step 2.

  Note that the sizes of the arrays are *minimum* sizes.  In practice they
  may be larger.  Do not call `size()` on these arrays, use the appropriate
  integer field.

  **Fields**

   * elnums_all: the element numbers on the original mesh of the elements on
                 the reduced mesh, length `shockmesh.numEl`  This maps element
                 numbers *from* the reduced mesh *to* the full mesh
   * elnum_mesh: the element numbers on the reduced mesh of the elements
                 on the reduced mesh, length `mesh.numEl`.  This maps element
                 numbers *from* the original mesh *to* the reduced mesh.
   * ifaces: array of [`RelativeInterface`](@ref) objects for the faces
            between elements in the reduced mesh (both elements with shocks in
            them and neighbors. length `shockmesh.numInterfaces`
   * bndryfaces: array of [`RelativeBoundary`](@ref) objects.  Contains only
                 boundary faces that are on the boundary of the original mesh
                 and the element has a shock in it (neighbor elements are
                 not included).
   * numShock: number of elements with shocks in them
   * numNeighbor: number of elements that do not have shocks in them, but
                  are adjacent to elements that do
   * numInterfaces: number of interfaces
   * numBoundaryFaces: number of boundary faces
   * numEl: total number of elements
"""
mutable struct ShockedElements{Tres}
  elnums_all::Vector{Int}  # union of elnums_shock and elnums_neighbor,
                           # and shared elements
  elnums_mesh::Vector{Int}  # temporary array, length mesh.numEl
                            # contains the indices of the selected elements
                            # in elnums_all
                            # contains the indices of elements in elnums_all
                            # TODO: can this be a BitArray?
  ee::Vector{Tres}  # the numerical viscoscity of each element in
                    # neighbor_elnums
  ifaces::Vector{RelativeInterface}
  bndryfaces::Vector{RelativeBoundary}

  shared_interfaces::Vector{Vector{RelativeInterface}}

  numShock::Int  # number of elements with shocks in them.  Because we try
                 # to pre-allocate elnums_shock, its length may be greater than
                 # the number of elements with shocks
  numNeighbor::Int  # number of elements that neighbor an element with a shock
                    # in it, but don't have a shock in it
  numShared::Vector{Int}  # number of elements that are owned by other processes
                          # but share a face with a shocked element
  numInterfaces::Int  # number of interfaces in ifaces
  numSharedInterfaces::Vector{Int}  # number of interfaces shared with each
                                    # peer process
  numBoundaryFaces::Int

  npeers::Int  # number of peers
  peer_indices::Vector{Int}  # array containing the *index* of the peer
                             # process in the original mesh list of peer
                             # processes

  numEl::Int  # total number of elements
  isNeumann::Bool  # if true, boundary is Neumann, else Dirichlet

  # useful ranges for iterating
  local_els::UnitRange{Int}
  neighbor_els::UnitRange{Int}
  shared_els::Vector{UnitRange{Int}}
  bndry_offsets::Vector{Int}


  # internal state used for push operations
  idx_all::Int  # current index in elnums_all
  sz_all::Int  # current size of elnums_all

  idx_if::Int  # current index in ifaces
  sz_if::Int  # current size of ifaces

  idx_sf::Vector{Int}  # current index in shared_interfaces
  sz_sf::Vector{Int}   # current size of shared_interfaces

  idx_b::Int  # current index in bndryfaces
  sz_b::Int  # current size of bndryfaces

  function ShockedElements{Tres}(mesh::AbstractMesh) where {Tres}

    # try to guess initial size
    size_guess = max(div(mesh.numEl, 10), 1)
    elnums_all = Array{Int}(size_guess)  # do this later
    elnums_mesh = zeros(mesh.numGlobalEl)
    ee = Array{Tres}(size_guess)
    ifaces = Array{RelativeInterface}(0)
    bndryfaces = Array{RelativeBoundary}(0)
    shared_interfaces = Vector{Vector{RelativeInterface}}(0)

    numShock = 0
    numNeighbor = 0
    numShared = Vector{Int}(0)
    numInterfaces = 0
    numSharedInterfaces = Vector{Int}(0)
    numBoundaryFaces = 0
    npeers = 0
    peer_indices = Vector{Int}(0)
    numEl = 0
    isNeumann = true
    local_els = 0:0
    neighbor_els = 0:0
    shared_els = Vector{UnitRange{Int}}(0)
    bndry_offsets = Vector{Int}(length(mesh.bndry_offsets))

    idx_all = 1
    sz_all = size_guess
    idx_if = 1
    sz_if = 0
    idx_sf = Vector{Int}(0)
    sz_sf = Vector{Int}(0)
    idx_b = 1
    sz_b = 0

    return new(elnums_all, elnums_mesh, ee,
               ifaces, bndryfaces, shared_interfaces, numShock,
               numNeighbor, numShared,
               numInterfaces, numSharedInterfaces, numBoundaryFaces, npeers,
               peer_indices, numEl, isNeumann,
               local_els, neighbor_els, shared_els, bndry_offsets,
               idx_all, sz_all, idx_if, sz_if, idx_sf, sz_sf, idx_b, sz_b)
  end
end

