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

  fill!(elnums_mesh, 0)
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
   * idx: the current index in data.interfaces
   * sz: the current size of data.interfaces
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

  #TODO: handle parallel interfaces
  @assert mesh.commsize == 1

  data.numNeighbor = data.idx_shock - 1 - data.numShock
  data.numInterfaces = idx_if - 1

  # setup elnums_all
  # TODO: I think the previous step could use the same array for elnums_shock
  #       and elnums_neighbor, making this step unnecessary

  numEl = data.numShock + data.numNeighbor
  data.numEl = numEl
  if length(data.elnums_all) < numEl
    resize!(data.elnums_all, numEl)
  end
  idx = 1
  @simd for i=1:data.numShock
    data.elnums_all[idx] = data.elnums_shock[i]
    idx += 1
  end
  @simd for i=1:data.numNeighbor
    data.elnums_all[idx] = data.elnums_neighbor[i]
    idx += 1
  end

  return nothing
end
