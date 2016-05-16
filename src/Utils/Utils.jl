@doc """
Module Utils:
  This module holds miscellaneous functions used throughout the code
"""->
module Utils

using ODLCommonTools
using ArrayViews
using MPI
using SummationByParts
include("parallel.jl")
include("io.jl")
export disassembleSolution, writeQ, assembleSolution, assembleArray, sview
export initMPIStructures, exchangeFaceData, verifyCommunication, getSendData
export exchangeElementData

@doc """
### Utils.disassembleSolution

  This takes eqn.q_vec (the initial state), and disassembles it into eqn.q, the
  3 dimensional array.  This function uses mesh.dofs
  to speed the process.

  This function also calls writeQ to do any requested output.

  Inputs:
    mesh
    sbp
    eqn
    opts

  This is a mid level function, and does the right thing regardless of equation
  dimension.

  Aliasing restrictions: none
"""->
# mid level function (although it doesn't need Tdim)
function disassembleSolution{T}(mesh::AbstractCGMesh, sbp,
                             eqn::AbstractSolutionData, opts, 
                             q_arr::AbstractArray{T, 3}, 
                             q_vec::AbstractArray{T, 1})
  # disassemble q_vec into eqn.
  for i=1:mesh.numEl  # loop over elements
    for j = 1:mesh.numNodesPerElement
      for k=1:size(q_arr, 1)
	      dofnum_k = mesh.dofs[k, j, i]
	      q_arr[k, j, i] = q_vec[dofnum_k]
      end
    end
  end

  writeQ(mesh, sbp, eqn, opts)

  return nothing
end


function disassembleSolution{T}(mesh::AbstractDGMesh, sbp,
                             eqn::AbstractSolutionData, opts, 
                             q_arr::AbstractArray{T, 3}, 
                             q_vec::AbstractArray{T, 1})
                             
  # no need to do any disassembly for DG
  writeQ(mesh, sbp, eqn ,opts)

end

@doc """
### Utils.writeQ

  This function writes the real part of the solution variables eqn.q to a space 
  delimited file called q.dat, controlled by the input options 'writeq', of type bool

  This is a high level function.
"""->
function writeQ(mesh, sbp, eqn, opts)

  if !opts["writeq"]
    return nothing
  end

  fname = "q.dat"
  rmfile(fname)
  writedlm(fname, eqn.q)

  return nothing
end



@doc """
### Utils.assembleSolution

  This function takes the 3D array of variables in arr and 
  reassmbles is into the vector res_vec.  Note that
  This is a reduction operation and zeros res_vec before performing the 
  operation, unless zero_res is set to false

  This is a mid level function, and does the right thing regardless of
  equation dimension
"""->
# mid level function (although it doesn't need Tdim)
function assembleSolution{Tmsh, Tsol, Tres}(mesh::AbstractCGMesh{Tmsh}, 
                         sbp, eqn::AbstractSolutionData{Tsol}, opts, 
                         arr::Abstract3DArray, res_vec::AbstractArray{Tres,1}, 
                         zero_resvec=true)
# arr is the array to be assembled into res_vec

#  println("in assembleSolution")
  if mesh.isDG
    return nothing
  end

  if zero_resvec
    fill!(res_vec, 0.0)
  end


  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      for k=1:size(arr, 1)  # loop over dofs on the node
        dofnum_k = mesh.dofs[k, j, i]
        res_vec[dofnum_k] += arr[k,j,i]
      end
    end
  end
  
  return nothing
end

function assembleSolution{Tmsh, Tsol, Tres}(mesh::AbstractDGMesh{Tmsh}, 
                         sbp, eqn::AbstractSolutionData{Tsol}, opts, 
                         arr::Abstract3DArray, res_vec::AbstractArray{Tres,1}, 
                         zero_resvec=true)

  # no need to do anything for DG meshes
end




# mid level function (although it doesn't need Tdim)
function assembleArray{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                         sbp, eqn::AbstractSolutionData{Tsol}, opts, 
                         arr::Abstract3DArray, res_vec::AbstractArray{Tres,1}, 
                         zero_resvec=true)
# arr is the array to be assembled into res_vec, using an assignment reduction

#  println("in assembleSolution")

  if zero_resvec
    fill!(res_vec, 0.0)
  end

  if mesh.numDofPerNode == 1
    offset = 0
  else
    offset = 1
  end

  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      for k=1:size(arr, 1)  # loop over dofs on the node

        dofnum_k = mesh.dofs[k, j, i]
        dofnum_k1 = div(dofnum_k, mesh.numDofPerNode) + offset # get node number

        res_vec[dofnum_k1] = arr[k,j,i]
      end
    end
  end
  
  return nothing
end


# it would be better if this used @boundscheck
@doc """
### Utils.safe_views

  This bool value controls whether the function named sview refers to 
  view or unsafe_view from the ArrayViews package
"""->
global const safe_views = false
if safe_views
  global const sview = ArrayViews.view
else
  global const sview = ArrayViews.unsafe_view
end

import Base.flush
function flush(f::IOBuffer)

end

end  # end module

  
