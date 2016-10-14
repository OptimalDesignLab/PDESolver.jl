@doc """
Module Utils:
  This module holds miscellaneous functions used throughout the code
"""->
module Utils

push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using ODLCommonTools
using ArrayViews
using MPI
using SummationByParts
using PdePumiInterface     # common mesh interface - pumi

include("parallel.jl")
include("io.jl")
include("logging.jl")
include("initialization.jl")
export disassembleSolution, writeQ, assembleSolution, assembleArray, sview
export calcNorm, calcMeshH
export initMPIStructures, exchangeFaceData, verifyCommunication, getSendData
export exchangeElementData
export @mpi_master, @time_all, print_time_all
export Timings, write_timings
export sharedFaceLogging
export createMeshAndOperator
export calcBCNormal
export applyPermRow, applyPermRowInplace, applyPermColumn
export applyPermColumnInplace, inversePerm, permMatrix, permMatrix!

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


@doc """
### Utils.calcNorm

  This function calculates the norm of a vector (of length numDof) using the
    SBP norm.

    Inputs:
      eqn:  an AbstractSolutionData
      res_vec:  vector to calculate the norm of

    Keyword arguments:
      strongres: if res_vec is the residual of the weak form, then
                 strongres=true computes (efficiently) the norm of the strong
                 form residual.  Default false
      globalnrm: compute the norm over all processes or not. Default true

    Returns:
      val:  norm of solution using SBP norm (Float64)

    There are no restrctions on the datatype of res_vec (ie. it can be complex)

    Aliasing restrictions: none

"""->
function calcNorm{T}(eqn::AbstractSolutionData, res_vec::AbstractArray{T}; strongres=false, globalnrm=true)
# calculates the norm of a vector using the mass matrix

  val = zero(real(res_vec[1]))

  if !strongres
    for i=1:length(res_vec)
      val += real(res_vec[i])*eqn.M[i]*real(res_vec[i])   # res^T M res
    end
  else  # strongres
    for i=1:length(res_vec)
      val += real(res_vec[i])*eqn.Minv[i]*real(res_vec[i])   # res^T M res
    end
  end

  eqn.params.time.t_allreduce = @elapsed if globalnrm
    val = MPI.Allreduce(val, MPI.SUM, eqn.comm)
  end

  val = sqrt(val)
  return val
end     # end of calcNorm function


@doc """
### Utils.calcMeshH

  This function calculates the average distance between nodes over the entire
  mesh.  This function allocates a bunch of temporary memory, so don't call
  it too often.  This is, strictly speaking, not quite accurate in parallel
  because the divison by length happens before the allreduce.

  Inputs:
    mesh
    eqn
    opts
"""->
function calcMeshH{Tmsh}(mesh::AbstractMesh{Tmsh}, sbp,  eqn, opts)
  jac_3d = reshape(mesh.jac, 1, mesh.numNodesPerElement, mesh.numEl)
  jac_vec = zeros(Tmsh, mesh.numNodes)
  assembleArray(mesh, sbp, eqn, opts, jac_3d, jac_vec)

  dim = mesh.dim
  # scale by the minimum distance between nodes on a reference element
  # this is a bit of an assumption, because for distorted elements this
  # might not be entirely accurate
  h_avg = sum(1./(jac_vec.^(1/dim)))/length(jac_vec)
  h_avg = MPI.Allreduce(h_avg, MPI.SUM, mesh.comm)/mesh.commsize
  h_avg *= mesh.min_node_dist
  return h_avg
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

#=
import Base.flush
function flush(f::IOBuffer)

end
=#

@doc """
### Utils.Timings

  This type accumulates the time spent in each part of the code.
"""->
type Timings
  # timings
  t_volume::Float64  # time for volume integrals
  t_face::Float64 # time for surface integrals (interior)
  t_source::Float64  # time spent doing source term
  t_sharedface::Float64  # time for shared face integrals
  t_bndry::Float64  # time spent doing boundary integrals
  t_dataprep::Float64  # time spent preparing data
  t_stab::Float64  # time spent adding stabilization
  t_send::Float64  # time spent sending data
  t_wait::Float64  # time spent in MPI_Wait
  t_allreduce::Float64 # time spent in allreduce
  t_pert::Float64  # time spent applying peraturbation
  t_alloc::Float64  # time spent allocating the Jacobian
  t_insert::Float64  # time spent inserting values into matrix
  t_func::Float64  # time spent evaluating the residual
  t_color::Float64  # time spent evaluating the colors
  t_jacobian::Float64  # time spent computing jacobian
  t_solve::Float64  # linear solve time
  t_newton::Float64  # time spent in newton loop
  t_barrier::Float64  # time spent in MPI_Barrier
  t_barrier2::Float64
  t_barrier3::Float64
  t_barriers::Array{Float64, 1}

  function Timings()
    nbarriers = 7
    barriers = zeros(Float64, nbarriers)
    return new(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, barriers)
  end
end

@doc """
### Utils.write_timings

  Write the values in a Timings object to a file.  Also writes the names of
  fields to a separate file.

  Inputs:
    t: a  Timings object
    fname: the file name, without extension

  The values are written to the file fname.dat, and the names are written to
  fname_names.dat
"""->
function write_timings(t::Timings, fname::AbstractString)
  timing_names = fieldnames(t)
  nbarriers = length(t.t_barriers)
  nvals = length(timing_names) + nbarriers - 1
  vals = Array(Float64, nvals)
  val_names = Array(ASCIIString, nvals)

  # put all values except those from t_barriers into array
  pos = 1
  for i=1:length(timing_names)
    tname_i = timing_names[i]
    tname_i_str = string(tname_i)
    if tname_i_str != "t_barriers"
      vals[pos] = getfield(t, tname_i)
      val_names[pos] = tname_i_str
      pos += 1
    end
  end

  for i=1:length(t.t_barriers)
    vals[pos] = t.t_barriers[i]
    val_names[pos] = string("t_barriers_", i)
    pos += 1
  end

  fname_ext = string(fname, ".dat")
  writedlm(fname_ext, vals)

  fname2_ext = string(fname, "_names.dat")
  writedlm(fname2_ext, val_names)
end


"""
  Calculate the scaled normal vector in parametric coordinates from the
  face normal and scaled mapping jacobian.  `nrm2` is overwritten with
  the result.
"""
function calcBCNormal(params::AbstractParamType{2}, dxidx::AbstractMatrix, 
                    nrm::AbstractVector, nrm2::AbstractVector)

  nrm2[1] = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  nrm2[2] = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  return nothing
end

function calcBCNormal(params::AbstractParamType{3}, dxidx::AbstractMatrix, 
                    nrm::AbstractVector, nrm2::AbstractVector)

  n1 = nrm[1]; n2 = nrm[2]; n3 = nrm[3]
  nrm2[1] = dxidx[1,1]*n1 + dxidx[2,1]*n2 + dxidx[3,1]*n3
  nrm2[2] = dxidx[1,2]*n1 + dxidx[2,2]*n2 + dxidx[3,2]*n3
  nrm2[3] = dxidx[1,3]*n1 + dxidx[2,3]*n2 + dxidx[3,3]*n3

  return nothing
end


#------------------------------------------------------------------------------
# permutation functions
#------------------------------------------------------------------------------
"""
  Permute the rows of A according to the permvec, storing the result in B
  The permvec contains the source indices for each entry in B, ie.
  B[permvec[i]] comes from A[i].  This is consistent with the mathematical
  definition of a permutation that pre-multiplication by a permutation 
  matrix (obtained from permMatrix) is a row permutation, ie.
  B = P*A

  Aliasing: no aliasing allowed
"""
function applyPermRow{T}(permvec::AbstractVector, A::AbstractMatrix{T},
                            B::AbstractMatrix{T})

  m, n = size(A)
  for i=1:m
    idx = permvec[i]
    for j=1:n
      B[ i, j] = A[idx, j]
    end
  end

  return nothing
end

"""
  Like applyPermRow, but the result is returned in A.  Both A and B
  are overwritten.

  Aliasing: no aliasing allowed
"""
function applyPermRowInplace{T}(permvec::AbstractVector, A::AbstractMatrix{T},
                            B::AbstractMatrix{T})

  applyPermRow(permvec, A, B)

  # copy back
  for i=1:length(A)
    A[i] = B[i]
  end

  return nothing
end


"""
  Permute the columns of A according to permvec.  See applyPermRow for the
  definition of a permvec.  Note that for column permutation the 
  interpretation is that destination column permvec[i] comes from 
  source column i.  This is consistent with the notion that post-multiplying
  by a permutation matrix (obtained from permMatrix) is a column 
  permutation, i.e B = A*P

  Aliasing: no aliasing allowed
"""
function applyPermColumn{T}(permvec::AbstractVector, A::AbstractMatrix{T},
                            B::AbstractMatrix{T})

  m, n = size(A)
  for i=1:m
    for j=1:n
      B[i, permvec[j]] = A[i, j]
    end
  end

  return nothing
end

"""
  Like applyPermColumn, but the result is returned in A.  Both A and B
  are overwritten

  Aliasing: no aliasing allowed
"""
function applyPermColumnInplace{T}(permvec::AbstractVector, 
                            A::AbstractMatrix{T}, B::AbstractMatrix{T})

  applyPermColumn(permvec, A, B)

  for i=1:length(A)
    A[i] = B[i]
  end

  return nothing
end

"""
  Create a permutation matrix from a permutation vector.  Only select
  entries of A are overwritten.

"""
function permMatrix!(permvec::AbstractVector, A::AbstractMatrix)

  m, n = size(A)

  for i=1:m
    A[i, permvec[i]] = 1
  end

  return nothing
end

"""
  Create a permutation matrix from a permutation vector.  The element type
  of the returned matrix is Int.
"""
function permMatrix(permvec::AbstractVector)
  m = length(permvec)
  A = zeros(Int, m, m)
  permMatrix!(permvec, A)

  return A
end


"""
  Compute the permutation vector that corresponds to the inverse permutation

  Inputs:
    permvec: the original permutation vector

  Inputs/Outputs:
    invperm: the inverse permutation vector

  Aliasing: no aliasing allowed
"""
function inversePerm(permvec::AbstractVector, invperm::AbstractVector)

  n = length(permvec)
  for i=1:n
    idx = permvec[i]
    invperm[idx] = i
  end

  return nothing
end

# TODO: write functions to apply inverse permutation from permvec, without
#       needing to explicetly compute the inverse permutation vector

end  # end module

  
