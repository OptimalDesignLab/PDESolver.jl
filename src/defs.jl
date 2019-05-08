# some definitions used throughout TIcon

"""
  Abstract type for _AssembleElementData.  See that type for details.
"""
abstract type AssembleElementData end


#------------------------------------------------------------------------------
# sparsity pattern enums

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


