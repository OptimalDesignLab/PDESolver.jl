# Description:  
#   A better interface for preallocation and populating sparse matrices

# is SparseMatrixCSC exported from base?
# Ti is type of indices
# Tv is typeof of values
import Base.SparseMatrixCSC


"""
### ODLCommonTools.SparseMatrixCSC

  Construct a SparseMatrixCSC for a DG mesh, using the pertNeighborEls
  to determine connectivity and dofs to determine the non-zero indices.
  This sparsity structure is exact for SBPOmega but an over-estimate for
  other operators.

  **Inputs**

   * mesh: a DG mesh
   * Tv: element type of the matrix (the indices are Ints because UMFPack is
         picky)

  **Outputs**

   * mat: the matrix
"""
function SparseMatrixCSC(mesh::AbstractDGMesh, ::Type{Tv}) where Tv
# construct a SparseMatrixCSC that preallocates the space needed for values
# tightly
# this should be exact for DG and a slight overestimate for CG

#  println("----- entered SparseMatrixCSC constructor -----")

  # calculate lengths of vectors
  # nel_per_el: number of elements (including self) each element is connected to
  if mesh.coloringDistance == 2
    println("distance-2 coloring")
    nel_per_el = size(mesh.neighbor_nums, 1)
  else
    throw(ErrorException("Unsupported coloring distance"))
  end

  # calculate some useful values
  ndof = mesh.numDof
  nDofPerElement = mesh.numNodesPerElement*mesh.numDofPerNode
  nvals_per_column = mesh.numNodesPerElement*mesh.numDofPerNode*nel_per_el
  nvals = nvals_per_column*mesh.numDof

  # figure out the offset of first dof on the element
  min_dof_per_element = zeros(eltype(mesh.dofs), mesh.numEl)
  for i=1:mesh.numEl
    dofs_i = view(mesh.dofs, :, :, i)
    min_dof, idx = findmin(dofs_i)
    min_dof_per_element[i] = min_dof
  end

  # the permutation vector is the order in which to visit the elements
  perm = sortperm(min_dof_per_element)

  # visit the elements in the perm order to calculate colptrs
  starting_offset = zeros(eltype(mesh.dofs), mesh.numEl)
  starting_offset[1] = 0
  colptr = Array{Int64}(ndof+1)
  rowvals = zeros{Int64}( nvals)
  colptr[1] = 1
  colptr_pos = 2
  nnz_curr = 0  # current element number of neighbors
  # now handle the rest of the mesh
  for i=1:mesh.numEl
    el_i = perm[i]
    # count number of elements
    nnz_curr = 0
    for j=1:size(mesh.pertNeighborEls, 2)
      if mesh.pertNeighborEls[el_i, j] > 0
        nnz_curr += 1
      end
    end
    # set the colptr values for all dofs on the current element (because they
    # are all the same)
    for j=1:nDofPerElement
      colptr[colptr_pos] = colptr[colptr_pos-1] + nnz_curr*nDofPerElement
      colptr_pos += 1
    end

  end

  # set up rowvals
  dofs_i = zeros(eltype(mesh.dofs), nvals_per_column)
  elnums_i = zeros(eltype(mesh.pertNeighborEls), nel_per_el)
  # loop over elements because all nodes on an element have same sparsity
  for i=1:mesh.numEl

    # get the element numbers
    pos = 1
    for j=1:size(mesh.pertNeighborEls, 2)
      val = mesh.pertNeighborEls[i, j]
      if val > 0
        elnums_i[pos] = val
        pos += 1
      end
    end

    # copy dof numbers into array
    for j = 1:(pos-1)
      el_j = elnums_i[j]
      dofs_j = view(mesh.dofs, :, :, el_j)
      src_j = view(dofs_i, ((j-1)*nDofPerElement+1):(j*nDofPerElement))
      copyDofs(dofs_j, src_j)
    end

    dofs_used = view(dofs_i, 1:(pos-1)*nDofPerElement)
    # sort them
    sort!(dofs_used)
    @assert dofs_used[1] != 0

    ndof_used = length(dofs_used)
    min_dof, idx = findmin(view(mesh.dofs, :, :, i))

    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        # compute starting location in rowvals
        mydof = mesh.dofs[k, j, i]
        start_idx = colptr[mydof]
        # set them in rowvals
        for p=1:ndof_used
          idx = start_idx + p - 1
          rowvals[idx] = dofs_used[p]
        end
      end
    end

  end  # end loop over elements

  nzvals = zeros(Tv, nvals)
  return SparseMatrixCSC(ndof, ndof, colptr, rowvals, nzvals)

end

function copyDofs(src::AbstractArray{T, 2}, dest::AbstractArray{T, 1}) where T
  pos = 1
  for i=1:size(src, 2)
    for j=1:size(src, 1)
      dest[pos] = src[j, i]
      pos += 1
    end
  end
end


function SparseMatrixCSC(sparse_bnds::AbstractArray{Ti, 2}, Tv::DataType) where Ti
# TODO: @doc this
# preallocate matrix based on maximum, minimum non zero
# rows in each column
# the type of sparse_bnds is used for the indicies
# the type of val is used for the values
# the value of val itself is never used

  println("creating SparseMatrixCSC")

  (tmp, n) = size(sparse_bnds)
  num_nz = 0  # accumulate number of non zero entries

  m = maximum(sparse_bnds)  # get number of rows
  colptr = Array{Int64}(n+1)  # should be Ti

  if sparse_bnds[1,1] != 0
    colptr[1] = 1
  else
    colptr[1] = 0
  end

  # count number of non zero entries, assign column pointers
  for i=2:(n+1)
    min_row = sparse_bnds[1, i-1]
    max_row = sparse_bnds[2, i-1]

    num_nz += max_row - min_row + 1
    colptr[i] = num_nz + 1
  end

  rowval = zeros(Int64, num_nz)  # should be Ti
  nzval = zeros(Tv, num_nz)

  # populate rowvals
  pos = 1
  for i=1:n
    num_vals_i = colptr[i+1] - colptr[i]
    min_row = sparse_bnds[1, i]

    # write row values to row values
    for j=1:num_vals_i
      rowval[pos] = min_row + j - 1
      pos += 1
    end
  end

  @assert pos == num_nz + 1  # check for sanity

  println("average bandwidth = ", pos/m)
  return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end


"""
  disc_type constant for SparseMatrixCSC constructor, indicating the face nodes
  of an element are connected to the face nodes of the element that shares
  the face.
"""
global const INVISCID=1

"""
  disc_type constant for SparseMatrixCSC constructor, indicating all nodes
  of an element are connected to all face nodes of its neighbours
"""
global const VISCOUS=2

"""
  disc_type constant for SparseMatrixCSC constructor, indicating all nodes
  of an element are connected to all nodes of its neighbors
"""
global const COLORING=3

"""
  disc_type constant for sparsity pattern of viscous terms when the T4 = 0.
  In this case, all nodes in the stencil of the interpolation operator are
  connected to all nodes of the adjacent element.
"""
global const VISCOUSTIGHT=4

"""
  Construct a SparseMatrixCSC using the exact sparsity pattern, taking into
  account the type of SBP face operator.

  **Inputs**

   * mesh: a DG mesh
   * Tv: element type of the matrix (the indices are Ints because UMFPack is
         picky)
   * disc_type: discretization type, 1 = Inviscid, 2 = viscous
   * face_type: This determined what kind of sbpface is used to compute the
                sparsity pattern. 1 = sbpface is a DenseFace, 2 = sbpface is
                a SparseFace

  **Outputs**

   * mat: the matrix

"""
function SparseMatrixCSC(mesh::AbstractDGMesh, ::Type{Tv}, disc_type::Integer, face_type::Integer) where Tv

  sbpface = mesh.sbpface
  dnnz, onnz = getBlockSparsityCounts(mesh, mesh.sbpface, disc_type, face_type)
  bs = mesh.numDofPerNode

  nzval = zeros(Tv, bs*bs*sum(dnnz))
  rowval = zeros(Int, bs*bs*sum(dnnz))
  colptr = zeros(Int, mesh.numDof + 1)

  # this is really a block matrix so each entry in dnnz describes bs rows

  # set up colptr
  # the matrix is structurally symmetric, so we can use dnnz (which are the
  # per-row values) to set up colptr

  colptr[1] = 1
  for i=1:(mesh.numDof)
    block_idx = div(i - 1, bs) + 1
    nnz_i = dnnz[block_idx]*bs
    colptr[i+1]  = colptr[i] + nnz_i
  end

  @assert colptr[end] - 1 == length(rowval)

  # set up rowvals
  # put all the rowvals into the matrix first, sort them later

  # volume terms
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode

        dof = mesh.dofs[k, j, i]
        idx = colptr[dof]  # index in rowval
        # connect this dof to all other dofs on this element
        for p=1:mesh.numNodesPerElement
          for q=1:mesh.numDofPerNode
            rowval[idx] = mesh.dofs[q, p, i]
            idx += 1
          end
        end

      end  # end loop k
    end  # end loop j
  end  # end loop i


  # interface terms
  nonstencil_table = getNonStencilNodes(mesh.numNodesPerElement, mesh.sbpface)
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]

    permL = sview(sbpface.perm, :, iface_i.faceL)
    permR = sview(sbpface.perm, :, iface_i.faceR)
    nonstencilL = sview(nonstencil_table, :, iface_i.faceL)
    nonstencilR = sview(nonstencil_table, :, iface_i.faceR)

    if disc_type == INVISCID
      for j=1:length(permL)
        for k=1:mesh.numDofPerNode
        
          if face_type == 1  # all permL to all permR
            dofL = mesh.dofs[k, permL[j], iface_i.elementL]
            dofR = mesh.dofs[k, permR[j], iface_i.elementR]

            # get the piece of rowvals for the left and right dofs
            idx_rangeL = colptr[dofL]:(colptr[dofL+1]-1)
            idx_rangeR = colptr[dofR]:(colptr[dofR+1]-1)
            rowval_rangeL = sview(rowval, idx_rangeL)
            rowval_rangeR = sview(rowval, idx_rangeR)


            idxL = getStartIdx(rowval_rangeL)
            idxR = getStartIdx(rowval_rangeR)

            # connect to all dofs in perm on the other element
            for p=1:length(permL)
              for q=1:mesh.numDofPerNode
                dofkL = mesh.dofs[q, permL[p], iface_i.elementL]
                dofkR = mesh.dofs[q, permR[p], iface_i.elementR]

                rowval_rangeL[idxL] = dofkR
                idxL += 1
                rowval_rangeR[idxR] = dofkL
                idxR += 1
              end
            end

          else  # face_type == 2, diagonalE
            nodeL = permL[j]
            # get corresponding node on other element
            pnbr = sbpface.nbrperm[j, iface_i.orient]
            nodeR = permR[pnbr]

            dofL = mesh.dofs[k, nodeL, iface_i.elementL]
            dofR = mesh.dofs[k, nodeR, iface_i.elementR]

            # get piece of rowvals for left and right dofs
            idx_rangeL = colptr[dofL]:(colptr[dofL+1]-1)
            idx_rangeR = colptr[dofR]:(colptr[dofR+1]-1)
            rowval_rangeL = sview(rowval, idx_rangeL)
            rowval_rangeR = sview(rowval, idx_rangeR)

            idxL = getStartIdx(rowval_rangeL)
            idxR = getStartIdx(rowval_rangeR)

            # connect to all dofs of corresponding node on other element
            for q=1:mesh.numDofPerNode
              rowval_rangeL[idxL] = mesh.dofs[q, nodeR, iface_i.elementR]
              idxL += 1
              rowval_rangeR[idxR] = mesh.dofs[q, nodeL, iface_i.elementL]
              idxR += 1
            end
          end  # end if face_type
        end  # end loop k
      end  # end loop j

    elseif disc_type == VISCOUS
      # connect all volume nodes of elementL to nodes in perm of elementR
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          dofL = mesh.dofs[k, j, iface_i.elementL]
          dofR = mesh.dofs[k, j, iface_i.elementR]
          idx_rangeL = colptr[dofL]:(colptr[dofL+1]-1)
          idx_rangeR = colptr[dofR]:(colptr[dofR+1]-1)
          rowval_rangeL = sview(rowval, idx_rangeL)
          rowval_rangeR = sview(rowval, idx_rangeR)

          idxL = getStartIdx(rowval_rangeL)
          idxR = getStartIdx(rowval_rangeR)
          for p=1:length(permR)
            for q=1:mesh.numDofPerNode
              rowval_rangeL[idxL] = mesh.dofs[q, permR[p], iface_i.elementR]
              idxL += 1
              rowval_rangeR[idxR] = mesh.dofs[q, permL[p], iface_i.elementL]
              idxR += 1
            end
          end

        end  # end loop k
      end  # end loop j

    elseif disc_type == VISCOUSTIGHT
      # connect nodes on perm of elementL to all volume nodes of elementR
      
      # do permL to all of R
      for j=1:length(permL)
        dofL = mesh.dofs[k, permL[j], iface_i.elementL]
        dofR = mesh.dofs[k, permR[j], iface_i.elementR]
        idx_rangeL = colptr[dofL]:(colptr[dofL+1]-1)
        idx_rangeR = colptr[dofR]:(colptr[dofR+1]-1)
        rowval_rangeL = sview(rowval, idx_rangeL)
        rowval_rangeR = sview(rowval, idx_rangeR)

        idxL = getStartIdx(rowval_rangeL)
        idxR = getStartIdx(rowval_rangeR)
        for p=1:mesh.numNodesPerElement
          for k=1:mesh.numDofPerNode
            rowval_rangeL[idxL] = mesh.dofs[k, p, iface_i.elementR]
            idxL += 1
            rowval_rangeR[idxR] = mesh.dofs[k, p, iface_i.elementL]
            idxR += 1
          end
        end
      end

      # do nonstencilL to permR
      for j=1:length(nonstencilL)
        dofL = mesh.dofs[k, nonstencilL[j], iface_i.elementL]
        dofR = mesh.dofs[k, nonstencilR[j], iface_i.elementR]
        idx_rangeL = colptr[dofL]:(colptr[dofL+1]-1)
        idx_rangeR = colptr[dofR]:(colptr[dofR+1]-1)
        rowval_rangeL = sview(rowval, idx_rangeL)
        rowval_rangeR = sview(rowval, idx_rangeR)

        idxL = getStartIdx(rowval_rangeL)
        idxR = getStartIdx(rowval_rangeR)

        for p=1:length(permR)
          for k=1:mesh.numDofPerNode
            rowval_rangeL[idxL] = mesh.dofs[k, permR[p], iface_i.elementR]
            idxL += 1
            rowval_rangeR[idxR] = mesh.dofs[k, permL[p], iface_i.elementL]
            idxR += 1
          end
        end
      end


    else  # disc_type = COLORING
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          dofL = mesh.dofs[k, j, iface_i.elementL]
          dofR = mesh.dofs[k, j, iface_i.elementR]
          idx_rangeL = colptr[dofL]:(colptr[dofL+1]-1)
          idx_rangeR = colptr[dofR]:(colptr[dofR+1]-1)
          rowval_rangeL = sview(rowval, idx_rangeL)
          rowval_rangeR = sview(rowval, idx_rangeR)

          idxL = getStartIdx(rowval_rangeL)
          idxR = getStartIdx(rowval_rangeR)
          for p=1:mesh.numNodesPerElement
            for q=1:mesh.numDofPerNode
              rowval_rangeL[idxL] = mesh.dofs[q, p, iface_i.elementR]
              idxL += 1
              rowval_rangeR[idxR] = mesh.dofs[q, p, iface_i.elementL]
              idxR += 1
            end
          end

        end  # end loop k
      end  # end loop j

    end  # end if disc_type
  end  # end loop i


  # check the structure is correct and sort the rowvals
  for i=1:mesh.numDof

    rng = colptr[i]:(colptr[i+1]-1)
    # check that no zeros remain
    for j in rng
      @assert rowval[j] != 0
    end

    rowval_i = sview(rowval, rng)
    sort!(rowval_i)
  end

  @assert colptr[mesh.numDof + 1] == length(rowval) + 1

  

  return SparseMatrixCSC(mesh.numDof, mesh.numDof, colptr, rowval, nzval)
end

"""
  Helper function to get the index of the first zero entry of the array.
  Returns zero otherwise

  **Inputs**

   * arr: the array
"""
function getStartIdx(arr::AbstractVector)

  for i=1:length(arr)
    if arr[i] == 0
      return i
    end
  end

  return 0
end

"""
  Compute the number of block (mesh.numDofPerNode x mesh.numDofPerNode) that
  each blocks is connected to in the jacobian matrix

  **Inputs**

   * mesh: a DG mesh
   * sbpface: an AbstractFace
   * disc_type: discretization type (INVISCID, or VISCOUS)
   * face_type: what type of sbpface to use, 1 = DenseFace, 2 = SparseFAce

  **Outputs**

   * dnnz: vector of length div(mesh.numDof, mesh.numDofPerNode) containing
           the number of blocks each block is connected to in the *local*
           part of the mesh.  Element type UInt16
   * onnz: vector, same length as `dnnz`, containing the number of blocks
           each block is connected to in the *non-local* part of the mesh
"""
function getBlockSparsityCounts(mesh::AbstractDGMesh, sbpface,
                                disc_type::Integer, face_type::Integer)
# disc_type: 1 = inviscid: stencil defined by f(uL, uR)a
#           entropy stable code is covered by disc_type == 1 (at least for
#           meshes with only 1 type of element
# disc_type: 2 = viscous: stencil defined by f(D*uL, uR)
  if face_type != 1 && face_type != 2
    error("unsupported AbstractFace type: $(typeof(sbpface))")
  end

  @assert disc_type == INVISCID || disc_type == VISCOUS || disc_type == COLORING

  bs = mesh.numDofPerNode
  @assert mesh.numDof % bs == 0
  nblocks = div(mesh.numDof, bs)

  dnnz = zeros(UInt16, nblocks)
  onnz = zeros(UInt16, nblocks)

  # volume terms
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      blocknum = div(mesh.dofs[1, j, i] - 1, bs) + 1
      dnnz[blocknum] += mesh.numNodesPerElement
    end
  end

  # face terms
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]

    permL = sview(sbpface.perm, :, iface_i.faceL)
    permR = sview(sbpface.perm, :, iface_i.faceR)

    if disc_type == INVISCID
      for j=1:length(permL)
        dofL = mesh.dofs[1, permL[j], iface_i.elementL]
        dofR = mesh.dofs[1, permR[j], iface_i.elementR]
        block_nodeL = div(dofL - 1, bs) + 1
        block_nodeR = div(dofR - 1, bs) + 1

        if face_type == 1  # all permL to all permR
          dnnz[block_nodeL] += length(permR)
          dnnz[block_nodeR] += length(permL)
        else  # face_type == 2, diagonalE
          dnnz[block_nodeL] += 1
          dnnz[block_nodeR] += 1
        end

      end  # end loop j

    elseif disc_type == VISCOUS
      # term if type f(R*D*qL, qR), and the transpose
      # face_type doesn't matter here because the D*qL term connects all
      # volume nodes together
      for j=1:mesh.numNodesPerElement
        dofL = mesh.dofs[1, j, iface_i.elementL]
        dofR = mesh.dofs[1, j, iface_i.elementR]
        block_nodeL = div(dofL - 1, bs) + 1
        block_nodeR = div(dofR - 1, bs) + 1

        dnnz[block_nodeL] += length(permR)
        dnnz[block_nodeR] += length(permL)
      end  # end loop j

    elseif disc_type == VISCOUSTIGHT
      # when T4=0, the important terms are R^T * T3 * D*u and D^T * T2*R*u
      
      for j=1:length(permL)
        dofL = mesh.dofs[1, permL[j], iface_i.elementL]
        dofR = mesh.dofs[1, permR[j], iface_i.elementR]
        block_nodeL = div(dofL - 1, bs) + 1
        block_nodeR = div(dofR - 1, bs) + 1

        dnnz[block_nodeL] += mesh.numNodesPerElement
        dnnz[block_nodeR] += mesh.numNodesPerElement
      end



    else  # disc_type == COLORING
      # all volume nodes to all volume nodes
      for j=1:mesh.numNodesPerElement
        dofL = mesh.dofs[1, j, iface_i.elementL]
        dofR = mesh.dofs[1, j, iface_i.elementR]
        block_nodeL = div(dofL - 1, bs) + 1
        block_nodeR = div(dofR - 1, bs) + 1

        dnnz[block_nodeL] += mesh.numNodesPerElement
        dnnz[block_nodeR] += mesh.numNodesPerElement
      end  # end loop j
    end  # end if disc_type
  end  # end loop i


  # now do shared faces
  # only figure out the sparsity for the local element (elementL)
  for peer=1:mesh.npeers
    interfaces = mesh.shared_interfaces[peer]

    for i=1:length(interfaces)
      iface_i = interfaces[i]

      permL = sview(sbpface.perm, :, iface_i.faceL)
      permR = sview(sbpface.perm, :, iface_i.faceR)

      if disc_type == INVISCID
        for j=1:length(permL)
          dofL = mesh.dofs[1, permL[j], iface_i.elementL]
          block_nodeL = div(dofL - 1, bs) + 1

          if face_type == 1
            onnz[block_nodeL] += length(permR)
          else
            onnz[block_nodeL] += 1
          end

        end  # end loop j

      elseif disc_type == VISCOUS
        for j=1:mesh.numNodesPerElement
          dofL = mesh.dofs[1, j, iface_i.elementL]
          block_nodeL = div(dofL - 1, bs) + 1

          onnz[block_nodeL] += length(permR)
        end

      elseif disc_type == VISCOUSTIGHT
        for j=1:length(permL)
          dofL = mesh.dofs[1, permL[j], iface_i.elementL]
          block_nodeL = div(dofL - 1, bs) + 1
          onnz[block_nodeL] += mesh.numNodesPerElement
        end

      else
        for j=1:mesh.numNodesPerElement
          dofL = mesh.dofs[1, j, iface_i.elementL]
          block_nodeL = div(dofL - 1, bs) + 1
          onnz[block_nodeL] += mesh.numNodesPerElement
        end

      end  # end if disc_type

    end  # end loop i
  end  # end loop peer

  # no need to do boundaries, their sparsity is covered by the volume terms

  return dnnz, onnz
end







 
#------------------------------------------------------------------------------
# Access methods
import Base.getindex
import Base.setindex!
import Base.fill!

const band_dense = false

if band_dense
  # setindex for dense within the band matrix
  function setindex!(A::SparseMatrixCSC{T, Ti}, v, i::Integer, j::Integer) where {T, Ti}
  # TODO: @doc this
  # get a nonzero value from A
  # for speed, no bounds checking

  #  println("using custom setindex")

    row_start = A.colptr[j]
    row_end = A.colptr[j+1] - 1
    row_min = A.rowval[row_start]
    row_max = A.rowval[row_end]

    if i < row_min || i > row_max
      println(STDERR, "Warning: Cannot change sparsity pattern of this matrix")
      println(STDERR, "    i = ", i, ", j = ", j, " value = ", v)
      return A
    end

    offset = i - row_min  # offset due to row
    valindex = row_start + offset
    A.nzval[valindex] = v

    return A

  end

  function getindex(A::SparseMatrixCSC{T}, i::Integer, j::Integer) where T
  # TODO: @doc this
  # get a nonzero value from A
  # for speed, no bounds checking

  #  println("using custom getindex")

    row_start = A.colptr[j]
    row_end = A.colptr[j+1] - 1
    row_min = A.rowval[row_start]
    row_max = A.rowval[row_end]

    if i < row_min || i > row_max
      return zero(eltype(A.nzval))
    end

    offset = i - row_min  # offset due to row
    valindex = row_start + offset

    return A.nzval[valindex]

  end


else
  function setindex!(A::SparseMatrixCSC{T, Ti}, v, i::Integer, j::Integer) where {T, Ti}
    row_start = A.colptr[j]
    row_end = A.colptr[j+1] - 1
    rowvals_extract = sview(A.rowval, row_start:row_end)
    val_idx = fastfind(rowvals_extract, i)

    #TODO: comment this out after testing
    if val_idx == 0
      throw(ErrorException("entry $i, $j not found"))
    end

    idx = row_start + val_idx - 1
    A.nzval[idx] = v

    return A


  end
#=
  function getindex{T}(A::SparseMatrixCSC{T}, i::Integer, j::Integer)
    row_start = A.colptr[j]
    row_end = A.colptr[j+1] - 1
    rowvals_extract = sview(A.rowval, row_start:row_end)
    val_idx = fastfind(rowvals_extract, i)
    idx = row_start + val_idx -1
    return A.nzval[idx]
   
  end
=#
end  # end if band_dense
function fill!(A::SparseMatrixCSC, val)
  fill!(A.nzval, val)
  return nothing
end

@doc """
### ODLCommonTools.findfast

  This function searches a sorted array for a given value, returning 0
  if it is not found.  

  The algorithm is nearly branchless and performs well compared to
  standard implementations.  

  The search takes a maximum of log2(n) + 2 iterations when the requested
  value is present and n iteration if it is not found.

  Inputs:
    arr: array of integers
    val: value to find

  Outputs:
    idx: the index of the array containing the value, 0 if not found
"""->
function fastfind(a::AbstractArray{T}, val) where T <: Integer

  foundflag = false
  lbound = 1
  ubound = length(a)
  idx = lbound + div(ubound - lbound, 2)
#  itermax = floor(log2(length(a))) + 2
  itermax = length(a)
  itr = 0


#  println("lbound = ", lbound)
#  println("ubound = ", ubound)

  while ( a[idx] != val && itr <= itermax)
#    println("\ntop of loop, idx = ", idx)
    if a[idx] > val  # the value lies in the left half 
      ubound = idx
      idx = lbound + fld(ubound - lbound, 2)
#      println("updating ubound = ", ubound)
    else  # a[idx] < val  # value lies in the right half
      lbound = idx
      idx = lbound + cld(ubound - lbound, 2)
#      println("updating lbound = ", lbound)
    end
    
    itr += 1
  end

    successflag = (itr <= itermax)
  return idx*successflag
#  return idx
end


