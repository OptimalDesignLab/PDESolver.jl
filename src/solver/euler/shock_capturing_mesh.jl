# functions that extract the mininum data necessary from the mesh to do a
# DG-type shock capturing scheme.

#TODO: is this unused?
"""
  This function gets the list of elements (and other required data) to do
  shock capturing with a DG-type dissipation.

  **Inputs**
  
   * mesh
   * sbp
   * eqn
   * opts
   * sensor

  **Inputs/Outputs**

   * shockmesh: a [`ShockedElement`](@ref) object.  Will be reset before
                writing a new set of shocked elements
"""
function getShockCapturingMesh(mesh::AbstractMesh, sbp::AbstractOperator,
                               eqn::EulerData, opts, sensor::AbstractShockSensor,
                               shockmesh::ShockedElements)

  reset(shockmesh)
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)

    add_el = isShockElement(eqn.params, sbp, sensor, q_i, jac_i)

    if add_el
      push!(shockmesh, i)
    end
  end

  completeShockElements(shockmesh)

  return nothing
end



"""
  This function resets the `ShockElements` object to its (more-or-less) initial
  state while keeping the sizes of the existing arrays.  After calling reset
  a new set of shocked elements can be recorded, hopfully with fewer array
  resizes that the first time.

  **Inputs**

   * data: `ShockedElements`
"""
function reset(data::ShockedElements)

  fill!(data.elnums_mesh, 0)
  fill!(data.elnums_all, 0)
  data.numShock = 0
  data.numNeighbor = 0
  data.numInterfaces = 0
  data.numEl = 0

  data.idx_all = 1
  data.idx_if = 1
  data.idx_b = 1
  fill!(data.idx_sf, 1)

  return nothing
end

import Base.push!


"""
  Adds a new element to the list of elements that have a shock in them

  **Inputs**

   * data: `ShockedElements`
   * elnum: the element number
"""
@inline function push!(data::ShockedElements, elnum::Int)

  idx = data.idx_all
  sz = data.sz_all
  if idx > sz
    sz = max(sz*2, 8)
    resize!(data.elnums_all, sz)
    data.sz_all = sz
  end

  data.elnums_all[idx] = elnum
  data.elnums_mesh[elnum] = idx

  data.idx_all += 1
end


"""
  Adds a new element to the list of neighboring elements.  Internal function

  **Inputs**

   * data: `ShockedElements`
   * elnum: the new element number
   * idx: the index to insert the new element number
"""
@inline function push_neighbor(data::ShockedElements, elnum::Integer)

  idx = data.idx_all
  sz = data.sz_all

  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.elnums_all, sz)
    data.sz_all = sz
  end

  data.elnums_all[idx] = elnum
  data.elnums_mesh[elnum] = idx
  data.idx_all += 1

  return nothing
end

"""
  Internal function for adding element to list of shared elements (ie.
  elements that live on other processes that share a face with an element
  on this process that has a shock in it.

  **Inputs**

   * data: `ShockedElements`
   * elnum: new element number
"""
@inline function push_shared(data::ShockedElements, elnum::Integer)

  idx = data.idx_all
  sz = data.sz_all

  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.elnums_all, sz)
    data.sz_all = sz
  end

  data.elnums_all[idx] = elnum
  data.elnums_mesh[elnum] = idx
  data.idx_all += 1
  data.numShared[end] += 1

  return nothing
end




"""
  Internal function for pushing a new interface

  **Inputs**

   * data: `ShockElements`
   * iface: a `RelativeInterface`
   * idx: the current index in data.ifaces
   * sz: the current size of data.ifaces
"""
@inline function push_iface(data::ShockedElements, iface::RelativeInterface)

  idx = data.idx_if
  sz = data.sz_if

  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.ifaces, sz)
    data.sz_if = sz
  end

  data.ifaces[idx] = iface
  data.idx_if += 1

  return nothing
end


@inline function push_sharediface(data::ShockedElements, iface::RelativeInterface)

  idx = data.idx_sf[end]
  sz = data.sz_sf[end]
  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.shared_interfaces[end], sz)
    data.sz_sf[end] = sz
  end

  data.shared_interfaces[end][idx] = iface
  data.numSharedInterfaces[end] += 1
  data.idx_sf[end] += 1

  return nothing
end



"""
  Internal function for pushing new boundary

  **Inputs**

   * data: `ShockElements`
   * bndry: a `RelativeBoundary`

"""
@inline function push_bndry(data::ShockedElements, bndry::RelativeBoundary)

  idx = data.idx_b
  sz = data.sz_b

  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.bndryfaces, sz)
    data.sz_b = sz
  end

  data.bndryfaces[idx] = bndry
  data.idx_b += 1

  return nothing
end



"""
  Replaces the element numbers for an interface

  **Inputs**

   * iface: an `Interface`
   * elementL: the new left element number
   * elementR: the new right elenent number

  **Outputs**

   * iface: a new `Interface` object
"""
function replace_interface(iface::Interface, elementL::Integer, elementR::Integer)
  return Interface(elementL, elementR, iface.faceL, iface.faceR, iface.orient)
end

"""
  Similar to [`replace_interface`](@ref), but for `Boundary` objects.
"""
function replace_boundary(bndry::Boundary, elnum::Integer)
  return Boundary(elnum, bndry.face)
end

#TODO: make idx and sz fields of the type
#      make the push_* functions update the counts
#      get rid of shocked and neighbors arrays, use only elnums_all

"""
  After all the elements that have shocks in them have been added, finished
  constructing the mesh data structure.

  **Inputs**

   * mesh
   * shockmesh
"""
function completeShockElements(mesh::AbstractMesh, data::ShockedElements)

  # Use mesh.interfaces to get both data.interfaces and the neighbor elements

  # The purpose of this is to define a reduced element and interface numbering
  # scheme that only numbers the entities needed to compute the shock capturing
  # terms.
  # The element numbers in data.elnums_all are assigned 1:n for the
  # elements with positive shock indicator, and (n+1):m for the neighboring
  # elements

  data.numShock = data.idx_all - 1
  # neighbor stuff
  # iface stuff
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    elementL = data.elnums_mesh[iface_i.elementL]
    elementR = data.elnums_mesh[iface_i.elementR]
#=
    # if an element has not been seen before and its neighbor has a shock in
    # it, add to the list of neighbor elements
    if (elementL > 0 && elementL <= data.numShock) || (elementR > 0 && elementR <= data.numShock)

      if elementL == 0
        elnum_full = iface_i.elementL  # element numbering in the full numbering
        elementL = data.idx_all  # element number in the reduced numbering

        push_neighbor(data, elnum_full)
        
      elseif elementR == 0  # it shouldn't be possible to see the same pair of
                          # elements more than once (even on periodic meshes)
        elnum_full = iface_i.elementR
        elementR = data.idx_all

        push_neighbor(data, elnum_full)
      end
=#

    if (elementL > 0 && elementL <= data.numShock) &&
       (elementR > 0 && elementR <= data.numShock)
      # record the new interface
      iface_new = RelativeInterface(replace_interface(iface_i, elementL, elementR), i)
      push_iface(data, iface_new)
    end  # end if
  end  # end for

  data.numNeighbor = data.idx_all - 1 - data.numShock
  data.numInterfaces = data.idx_if - 1

  # get shared interfaces
  setupShockmeshParallel(mesh, data)

#=
 # get the list of boundary faces
  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    elnum = data.elnums_mesh[bndry_i.element]
    if elnum > 0 && elnum <= data.numShock
      bndry_new = RelativeBoundary(replace_boundary(bndry_i, elnum), i)
      push_bndry(data, bndry_new)
    end
  end
  data.numBoundaryFaces = data.idx_b - 1
=#
  # get boundary faces and Dirichlet info
  fill!(data.bndry_offsets, 0); data.bndry_offsets[1] = 1
  for i=1:mesh.numBC
    bc_range = mesh.bndry_offsets[i]:(mesh.bndry_offsets[i+1]-1)
    nfaces = 0
    for j in bc_range
      bndry_j = mesh.bndryfaces[j]
      elnum = data.elnums_mesh[bndry_j.element]
      # check if the element is in the data
      if elnum > 0 && elnum <= data.numShock
        bndry_new = RelativeBoundary(replace_boundary(bndry_j, elnum), j)
        push_bndry(data, bndry_new)
        nfaces += 1
      end
    end  # end j
    data.bndry_offsets[i+1] = data.bndry_offsets[i] + nfaces
  end  # end i

  # add all boundaries that are not on the boundary of the original domain.
  # Don't include them in data.bndry_offsets
  @assert mesh.coord_order == 1
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    elementL = data.elnums_mesh[iface_i.elementL]
    elementR = data.elnums_mesh[iface_i.elementR]

    # if an element has not been seen before and its neighbor has a shock in
    # it, add to the list of neighbor elements
    elL_nonzero = elementL > 0 && elementL <= data.numShock
    elR_nonzero = elementR > 0 && elementR <= data.numShock

    # only one element non zero
    if elL_nonzero || elR_nonzero && !(elL_nonzero && elR_nonzero)

      if elementL == 0
        # add face of elementR as boundary
        bndry_new = RelativeBoundary(Boundary(elementR, iface_i.faceR), i, -1)
        push_bndry(data, bndry_new)
        
      elseif elementR == 0
        # add face of elementL as boundary
        #TODO: record orientation info needed for normal vector
        bndry_new = RelativeBoundary(Boundary(elementL, iface_i.faceL), i)
        push_bndry(data, bndry_new)
      end
    end
  end

  data.numBoundaryFaces = data.idx_b - 1


  # compute final counts that will be useful for the user
  numEl = data.numShock + data.numNeighbor
  for i=1:data.npeers
    numEl += data.numShared[i]
  end
  data.numEl = numEl

  # setup ranges
  data.local_els = 1:data.numShock
  data.neighbor_els = (data.numShock+1):(data.numShock + data.numNeighbor)
  startidx = data.numShock + data.numNeighbor + 1
  data.shared_els = Vector{UnitRange{Int}}(data.npeers)
  for i=1:data.npeers
    data.shared_els[i] = startidx:(startidx + data.numShared[i] - 1)
    startidx = data.shared_els[i][end] + 1
  end


  return nothing
end


function setupShockmeshParallel(mesh::AbstractMesh, data::ShockedElements)

  # Its possible to re-use the arrays from last iteration, but it would be
  # tricky when the set of peer processes changes from one iteration to the next
  data.peer_indices = Array{Int}(0)
  data.shared_interfaces = Vector{Vector{RelativeInterface}}(0)
  data.numSharedInterfaces = Array{Int}(0)
  data.numShared = Array{Int}(0)
  data.idx_sf = Array{Int}(0)
  data.sz_sf = Array{Int}(0)
  const INITIAL_SIZE = 8
  found_peer = false
  for peer=1:mesh.npeers
    for i=1:length(mesh.shared_interfaces[peer])
      iface_i = mesh.shared_interfaces[peer][i]

      # elementL is always the local one, so only need to check that one
      elementL = data.elnums_mesh[iface_i.elementL]
      elementR = data.elnums_mesh[iface_i.elementR]

      if (elementL > 0) && (elementL <= data.numShock)

        # add the peer
        # need the found_peer flag because peer_indices[0] is an error, so
        # don't evaluate unless length(peer_indices) > 0
        if !found_peer || data.peer_indices[end] != peer
          found_peer = true
          data.npeers += 1
          push!(data.peer_indices, peer)
          # make the new vector non-zero size to save re-allocating it on the
          # first push
          push!(data.shared_interfaces, Vector{RelativeInterface}(INITIAL_SIZE))
          push!(data.numSharedInterfaces, 0)
          push!(data.numShared, 0)
          push!(data.idx_sf, 1)
          push!(data.sz_sf, INITIAL_SIZE)
        end

        # add elementR to the list of elements
        if elementR == 0
          elnum_full = iface_i.elementR
          elementR = data.idx_all

          push_shared(data, elnum_full)
        end

        # add the interface
        iface_new = RelativeInterface(replace_interface(iface_i, elementL, elementR), i)
        push_sharediface(data, iface_new)

      end  # end if elementL
    end  # end i
  end  # end peer

  return nothing
end
#TODO: have push_* functions return the new element number

#TODO: make one of these for regular meshes
"""
  Helper function to get the index of a shared element in the `RemoteMetrics`
  arrays.

  **Inputs**

   * data: a `ShockedElements`
   * mesh: the original mesh
   * peeridx: the peer index on the `ShockedElements`
   * elnum: the element number on the `ShockedElements`

  **Outputs**

   * idx: the index of the shared element in, for example, `RemoteMetrics`
          arrays
"""
function getSharedElementIndex(data::ShockedElements, mesh::AbstractMesh,
                               peeridx::Integer, elnum::Integer)

  peeridx2 = data.peer_indices[peeridx]
  firstnum = mesh.shared_element_offsets[peeridx2]
  elnum2 = data.elnums_all[elnum]
  return elnum2 - firstnum + 1
end
