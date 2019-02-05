# functions that extract the mininum data necessary from the mesh to do a
# DG-type shock capturing scheme.

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

    Se, ee = getShockSensor(eqn.params, sbp, sensor, q_i, jac_i)

    if ee > 0
      push!(shockmesh, i, ee)
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
  data.idx_shock = 1
  data.numShock = 0
  data.numNeighbor = 0
  data.numInterfaces = 0
  data.numEl = 0

  return nothing
end

import Base.push!

"""
  Adds a new element to the list of elements that have a shock in them

  **Inputs**

   * data: `ShockedElements`
   * elnum: the element number
   * ee: the viscoscity value
"""
function push!(data::ShockedElements, elnum::Int, ee::Number)

  idx = data.idx_shock
  
  if idx > length(data.elnums_shock)
    currsize = length(data.elnums_shock)
    newsize = max(currsize*2, 8)
    resize!(data.elnums_shock, newsize)
    resize!(data.ee, newsize)
  end

  data.elnums_shock[idx] = elnum
  data.elnums_mesh[elnum] = idx
  data.ee[idx] = ee

  data.idx_shock += 1
end


"""
  Sets the viscoscity for an element of the shockmesh.  This function
  should only be called *after* [`completeShockElements`](@ref) has been
  called.  Unlike `push!`, this function does not add a new element to the
  shockmesh, it only sets the viscoscity for an existing element.

  **Inputs**

   * data: `ShockedElements`
   * elnum: element number
   * ee: viscoscity value
"""
function setViscoscity(data::ShockedElements, elnum::Integer, ee::Number)

  data.ee[elnum] = ee
end

"""
  Gets the viscsocity for a given element.

  **Inputs**

   * data::`ShockedElements`
   * elnum: the element number

  **Outputs**

   * val: the viscoscity value
"""
function getViscoscity(data::ShockedElements, elnum::Integer)

  return data.ee[elnum]
end


"""
  Adds a new element to the list of neighboring elements.  Internal function

  **Inputs**

   * data: `ShockedElements`
   * elnum: the new element number
   * idx: the index to insert the new element number
   * sz: the current size of data.elnums_neighbor

  **Outputs**

   * newidx: the next index to insert to
   * sz: the size of `data.elnums_neighbor`
"""
@inline function push_neighbor(data::ShockedElements, elnum::Integer, idx::Integer,
                       sz::Integer)

  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.elnums_neighbor, sz)
  end

  data.elnums_neighbor[idx] = elnum
  data.elnums_mesh[elnum] = idx + data.numShock  #TODO: get rid of elnumsNeighbor


  return idx+1, sz
end



"""
  Internal function for pushing a new interface

  **Inputs**

   * data: `ShockElements`
   * iface: a `RelativeInterface`
   * idx: the current index in data.ifaces
   * sz: the current size of data.ifaces
"""
@inline function push_iface(data::ShockedElements, iface::RelativeInterface,
                    idx::Integer, sz::Integer)

  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.ifaces, sz)
  end

  data.ifaces[idx] = iface

  return idx+1, sz
end


@inline function push_sharediface(data::ShockedElements, iface::RelativeInterface,
                    idx::Integer, sz::Integer)

  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.shared_interfaces, sz)
  end

  data.shared_interfaces[end][idx] = iface
  data.numSharedInterfaces[end] += 1

  return idx+1, sz
end



"""
  Internal function for pushing new boundary

  **Inputs**

   * data: `ShockElements`
   * iface: a `RelativeBoundary`
   * idx: the current index in data.bndryfaces
   * sz: the current size of data.bndryfaces

"""
@inline function push_bndry(data::ShockedElements, bndry::RelativeBoundary,
                            idx::Integer, sz::Integer)

  if idx > sz
    sz = max(2*sz, 8)
    resize!(data.bndryfaces, sz)
  end

  data.bndryfaces[idx] = bndry

  return idx+1, sz
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

  data.numShock = data.idx_shock - 1
  # neighbor stuff
  idx_nb = 1
  sz_nb = length(data.elnums_neighbor)
  # iface stuff
  idx_if = 1  # iface stuff
  sz_if = length(data.ifaces)
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    elementL = data.elnums_mesh[iface_i.elementL]
    elementR = data.elnums_mesh[iface_i.elementR]

    # if an element has not been seen before and its neighbor has a shock in
    # it, add to the list of neighbor elements
    if (elementL > 0 && elementL <= data.numShock) || (elementR > 0 && elementR <= data.numShock)
      new_elnum = data.idx_shock

      if elementL == 0
        elnum_full = iface_i.elementL  # element numbering in the full numbering
        elementL = data.idx_shock  # element number in the reduced numbering
        data.idx_shock += 1

        idx_nb, sz_nb = push_neighbor(data, elnum_full, idx_nb, sz_nb)
        
      elseif elementR == 0  # it shouldn't be possible to see the same pair of
                          # elements more than once (even on periodic meshes)
        elnum_full = iface_i.elementR
        elementR = data.idx_shock
        data.idx_shock += 1

        idx_nb, sz_nb = push_neighbor(data, elnum_full, idx_nb, sz_nb)
      end

      # record the new interface
      iface_new = RelativeInterface(replace_interface(iface_i, elementL, elementR), i)
      idx_if, sz_if = push_iface(data, iface_new, idx_if, sz_if)
    end  # end if
  end  # end for

  data.numNeighbor = data.idx_shock - 1 - data.numShock
  data.numInterfaces = idx_if - 1

  # get shared interfaces
  setupShockmeshParallel(mesh, data, idx_nb, sz_nb)


 # get the list of boundary faces
  idx_b = 1
  sz_b = length(data.bndryfaces)
  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    elnum = data.elnums_mesh[bndry_i.element]
    if elnum > 0 && elnum <= data.numShock
      bndry_new = RelativeBoundary(replace_boundary(bndry_i, elnum), i)
      idx_b, sz_b = push_bndry(data, bndry_new, idx_b, sz_b)
    end
  end
  data.numBoundaryFaces = idx_b - 1


  # setup elnums_all
  # TODO: I think the previous step could use the same array for elnums_shock
  #       and elnums_neighbor, making this step unnecessary

  numEl = data.numShock + data.numNeighbor
  for i=1:data.npeers
    numEl += data.numShared[i]
  end

  data.numEl = numEl
  if length(data.elnums_all) < numEl
    resize!(data.elnums_all, numEl)
    resize!(data.ee, numEl)
  end
  idx = 1
  @simd for i=1:data.numShock
    data.elnums_all[idx] = data.elnums_shock[i]
    idx += 1
  end
  @simd for i=1:(data.numNeighbor + sum(data.numShared))
    data.elnums_all[idx] = data.elnums_neighbor[i]
    data.ee[idx] = 0  # for shared faces, this will be set to the correct value
                      # later
    idx += 1
  end

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


function setupShockmeshParallel(mesh::AbstractMesh, data::ShockedElements,
                               idx_nb::Integer, sz_nb::Integer)

  # Its possible to re-use the arrays from last iteration, but it would be
  # tricky when the set of peer processes changes from one iteration to the next
  data.peer_indices = Array{Int}(0)
  data.shared_interfaces = Vector{Vector{RelativeInterface}}(0)
  data.numSharedInterfaces = Array{Int}(0)
  data.numShared = Array{Int}(0)
  const INITIAL_SIZE = 8
  for peer=1:mesh.npeers
    idx_if = 1
    sz_if = INITIAL_SIZE
    for i=1:length(mesh.shared_interfaces[peer])
      iface_i = mesh.shared_interfaces[peer][i]

      # elementL is always the local one, so only need to check that one
      elementL = data.elnums_mesh[iface_i.elementL]
      elementR = data.elnums_mesh[iface_i.elementR]

      if (elementL > 0) && (elementL <= data.numShock)

        # add the peer
        if peer_indices[end] != peer
          data.npeers += 1
          push!(data.peer_indices, peer)
          # make the new vector non-zero size to save re-allocating it on the
          # first push
          push!(data.shared_interfaces, Vector{RelativeInterface}(INITIAL_SIZE))
          push!(data.numSharedInterfaces, 0)
          push!(data.numShared, 0)
        end


        # add elementR to the list of elements
        if elementR == 0
          elnum_full = iface_i.elementR
          elementR = data.idx_shock
          data.idx_shock += 1

          idx_nb, sz_nb = push_neighbor(data, elnum_full, idx_nb, sz_nb)
          data.numShared[end] += 1
        end
    

        # add the interface
        iface_new = RelativeInterface(replace_interface(iface_i, elementL, elementR), i)

        idx_if, sz_if = push_sharediface(data, iface_new, idx_if, sz_if)

      end  # end if elementL

    end  # end i
  end  # end peer

  return nothing
end


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
