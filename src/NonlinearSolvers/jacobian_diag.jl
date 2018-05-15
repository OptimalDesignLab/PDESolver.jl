# get the diagonal block of the jacobian



#------------------------------------------------------------------------------
# Block diagonal matrix type

"""
  A special array type for holding a block-diagonal matrix
"""
immutable DiagJac{T} <: AbstractArray{T, 2}
  A::Array{T, 3}
end

"""
  Outer constructor for [`DiagJac`](@ref).

  **Inputs**

   * T: element type of the matrix
   * blocksize: the blocksize of the diagonal blocks
   * nblock: number of blocks.

   The size of the matrix is `blocksize`*`nblock`.
   Only square block are supported at this time

"""
function DiagJac{T}(::Type{T}, blocksize::Integer, nblock::Integer)

  A = zeros(T, blocksize, blocksize, nblock)

  return DiagJac{T}(A)
end

import Base: length, size, getindex, setindex!, show

"""
  Returns the length of the array *as if it were a dense matrix*
"""
function length(A::DiagJac)

  blocksize, blocksize, nblock = size(A.A)
  n = blocksize*nblock
  return n*n
end

function size(A::DiagJac)

  blocksize, blocksize, nblock = size(A.A)

  n = blocksize*nblock
  return (n, n)
end

"""
  Internal function for computing the indices in the 3D array

  **Inputs**

   * A : DaigJac
   * i: first index
   * j: second index

  **Outputs**

   * i2: first index in 3D array
   * j2: second index in 3D array
   * i2block: 3rd index in 3D array
"""
function computeIndex(A::DiagJac, i::Integer, j::Integer)
# internal function for computing indices

  blocksize, blocksize, nblock = size(A.A)

  i2 = ((i-1) % blocksize ) + 1
  i2block = div(i-1, blocksize) + 1

  j2 = ((j-1) % blocksize) + 1
  j2block = div(j-1, blocksize) + 1

  @assert i2block <= nblock
  @assert i2block > 0
  @assert j2block <= nblock
  @assert j2block > 0

  @assert i2 > 0
  @assert i2 <= blocksize
  @assert j2 > 0
  @assert j2 <= blocksize

  @assert i2block == j2block

  return i2, j2, i2block
end


"""
  Indexing using non-block indices
"""
function getindex(A::DiagJac, i::Integer, j::Integer)

  i2, j2, i2block = computeIndex(A, i, j)
  return A.A[i2, j2, i2block]
end

"""
  Indexing using non-block indices
"""
function setindex!(A::DiagJac, v, i::Integer, j::Integer)

  i2, j2, i2block = computeIndex(A, i, j)
  return A.A[i2, j2, i2block] = v
end

function show(io::IO, A::DiagJac)

  m, n = size(A)
  println("DiagJac: $m x $n block diagonal matrix with block size: $(size( A.A, 1)), number of blocks: $(size(A.A, 3))")
  for i=1:size(A.A, 3)
    println("block ", i, " = ")
    println(A.A[:, :, i])
  end

  return nothing
end

"""
  Matrix-vector multiplication function for [`DiagJac`](@ref)

  **Inputs**

   * A: a `DiagJac`
   * x: vector to multiply against

  **Inputs/Outputs**

   * b: vector to be overwritten with the result
"""
function diagMatVec(A::DiagJac, x::AbstractVector, b::AbstractVector)


  blocksize, blocksize, nblock = size(A.A)

  for i=1:nblock
    idx = ((i-1)*blocksize + 1):(i*blocksize)

    A_ii = sview(A.A, :, :, i)
    x_ii = sview(x, idx)
    b_ii = sview(b, idx)

    smallmatvec!(A_ii, x_ii, b_ii)
  end

  return nothing
end
#------------------------------------------------------------------------------
# New AssembleElementData type for getting the diagonal only

type AssembleDiagJacData{T} <: AssembleElementData
  A::DiagJac{T}
end

function AssembleDiagJacData{T}(mesh, sbp, eqn, opts, jac::DiagJac{T})

  return AssembleDiagJacData{T}(jac)
end

"""
  Assembles the block diagonal part of the volume integral contribution only.
  The off diagonal part of `jac` is never touched.
"""
function assembleElement{T}(helper::AssembleDiagJacData, mesh::AbstractMesh,
                            elnum::Integer, jac::AbstractArray{T, 4})

  for p=1:mesh.numNodesPerElement
    q = p  # diagonal contribution only

    dof1 = mesh.dofs[1, p, elnum] 
    blocknum = div(dof1 - 1, mesh.numDofPerNode) + 1


    # get values
    @simd for j=1:mesh.numDofPerNode
      @simd for i=1:mesh.numDofPerNode
        helper.A.A[i, j, blocknum] += jac[i, j, p, q]
      end
    end
  end  # end loop p

  return nothing
end

function assembleInterface{T}(helper::AssembleDiagJacData,
                              sbpface::DenseFace,
                              mesh::AbstractMesh, iface::Interface,
                              jacLL::AbstractArray{T, 4},
                              jacLR::AbstractArray{T, 4},
                              jacRL::AbstractArray{T, 4},
                              jacRR::AbstractArray{T, 4})

  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)

  permL = sview(sbpface.perm, :, iface.faceL)
  permR = sview(sbpface.perm, :, iface.faceR)

    for p=1:sbpface.stencilsize
      q = p
      pL = sbpface.perm[p, iface.faceL]
      pR = sbpface.perm[p, iface.faceR]
      qL = pL
      qR = pR

      dof1 = mesh.dofs[1, pL, iface.elementL]
      blockL = div(dof1 - 1, numDofPerNode) + 1
      dof1 = mesh.dofs[1, pR, iface.elementR]
      blockR = div(dof1 - 1, numDofPerNode) + 1

      # put values into 2 x 2 block matrix
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          helper.A.A[i, j, blockL] +=  jacLL[i, j, pL, qL]
          helper.A.A[i, j, blockR] +=  jacRR[i, j, pR, qR]
        end
      end
    end

  return nothing
end

function assembleSharedFace{T}(helper::AssembleDiagJacData, sbpface::DenseFace,
                               mesh::AbstractMesh,
                               iface::Interface,
                               jacLL::AbstractArray{T, 4},
                               jacLR::AbstractArray{T, 4})

#    for p in permR # =1:numNodesPerElement
    for p=1:sbpface.stencilsize
      q = p
      pL = sbpface.perm[p, iface.faceL]
      qL = pL

      dof1 = mesh.dofs[1, pL, iface.elementL]
      blockidx = div(dof1 - 1, mesh.numDofPerNode) + 1

      # put values into 2 x 2 block matrix
      for j=1:mesh.numDofPerNode
        for i=1:mesh.numDofPerNode
          helper.A.A[i, j, blockidx] += jacLL[i, j, pL, qL]
        end
      end

    end  # end loop p


  return nothing
end

function assembleBoundary{T}(helper::AssembleDiagJacData, sbpface::DenseFace,
                               mesh::AbstractMesh,
                               bndry::Boundary,
                               jac::AbstractArray{T, 4})

  for p=1:sbpface.stencilsize
    q = p
    pL = sbpface.perm[p, bndry.face]
    qL = pL

    # get dofs for node p
    dof1 = mesh.dofs[1, pL, bndry.element]
    blockidx = div(dof1 - 1, mesh.numDofPerNode) + 1

    # get values
    for j=1:mesh.numDofPerNode
      for i=1:mesh.numDofPerNode
        helper.A.A[i, j, blockidx] += jac[i, j, pL, qL]
      end
    end

  end  # end loop p

  return nothing
end


