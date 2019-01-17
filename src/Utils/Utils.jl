@doc """
Module Utils:
  This module holds miscellaneous functions used throughout the code
"""->
module Utils

push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using ODLCommonTools
import ODLCommonTools.sview
using ArrayViews
import ArrayViews.view
using MPI
using SummationByParts
using PdePumiInterface     # common mesh interface - pumi
using Input

include("output.jl")
include("parallel.jl")
include("io.jl")
include("logging.jl")
include("projections.jl")
include("complexify.jl")
include("mass_matrix.jl")
include("curvilinear.jl")
include("area.jl")
include("checkpoint.jl")
include("interpolation.jl")

export free
export array1DTo3D, writeQ, array3DTo1D, removeComplex, assembleArray
export calcNorm, calcL2InnerProduct, calcMeshH, calcEuclidianNorm
#export initMPIStructures, exchangeFaceData, verifyCommunication, getSendData
#export startDataExchange
#export exchangeElementData
export  Timings, write_timings
export sharedFaceLogging
export calcBCNormal, calcBCNormal_revm, max_deriv_rev
export applyPermRow, applyPermRowInplace, applyPermColumn
export applyPermColumnInplace, inversePerm, permMatrix, permMatrix!
export arrToVecAssign
export fastzero!, fastscale!, removeComplex
export @verbose1, @verbose2, @verbose3, @verbose4, @verbose5, @unpack, @printit,
       assertArraysUnique, assertFieldsConcrete
# projections.jl functions
export getProjectionMatrix, projectToXY, projectToNT, calcLength

# complexify.jl functions
export absvalue, absvalue_deriv, absvalue_rev, atan2_rev

# output.jl
export printSolution, printCoordinates, printMatrix
export print_qvec_coords

# mass_matrix.jl
export calcMassMatrixInverse, calcMassMatrix, calcMassMatrixInverse3D,
       applyMassMatrixInverse

# curvilinear.jl
export calcSCurvilinear, calcSCurvilinear_rev, calcECurvilinear, calcDCurvilinear

# area.jl
export calcVolumeContribution!, calcVolumeContribution_rev!, calcProjectedAreaContribution!, calcProjectedAreaContribution_rev!, crossProd, crossProd_rev

# io.jl
export BufferedIO, BSTDOUT, BSTDERR

# parallel_types.jl
export SharedFaceData, getSharedFaceData, setParallelData, getParallelData,
       setNewTag

# parallel.jl
export startSolutionExchange, startSolutionExchange_rev,
       startSolutionExchange_rev2, exchangeData,
       exchangeData_rev,
       finishExchangeData, finishExchangeData_rev, finishExchangeData_rev2,
       finishSolutionBarExchange,
       @mpi_master, @time_all, print_time_all, verifyReceiveCommunication,
       MPITagManager, getNextTag, markTagUsed, freeTag, assertReceivesWaited,
       assertSendsWaited

# checkpoint.jl
export Checkpointer, AbstractCheckpointData, readCheckpointData,
       readLastCheckpointData, saveNextFreeCheckpoint, loadLastCheckpoint,
       countFreeCheckpoints,
       getLastCheckpoint, getOldestCheckpoint, freeOldestCheckpoint,
       freeCheckpoint, getNextFreeCheckpoint

# interpolation.jl
export interpField

"""
  Generic function to free any memory belonging to other libraries
"""
function free(x)

  error("generic free() reached")

  return nothing
end

@doc """
### Utils.array1DTo3D

  This takes a solution vector, and disassembles it into the 3D array form.
  This function uses mesh.dofs to speed the process.

  This function also calls writeQ to do any requested output.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * q_vec: vector, same shape as `eqn.q_vec

  **Inputs/Outputs**

   * q_arr: array, same shape as `eqn.q`, to be overwritten

  This is a mid level function, and does the right thing regardless of equation
  dimension.

  This function is optimized for DG meshes where `q_vec` and `q` might
  refer to the same array.

  Aliasing restrictions: none
"""->
# mid level function (although it doesn't need Tdim)
function array1DTo3D(mesh::AbstractCGMesh, sbp,
                          eqn::AbstractSolutionData, opts,
                          q_vec::AbstractVector,
                          q_arr::AbstractArray{T, 3}) where T
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


function array1DTo3D(mesh::AbstractDGMesh, sbp,
                          eqn::AbstractSolutionData, opts,
                          q_vec::AbstractVector,
                          q_arr::AbstractArray{T, 3}) where T

  # we assume the memory layouts of q_arr and q_vec are the same
  if pointer(q_arr) != pointer(q_vec)
    for i = 1:length(q_vec)
      q_arr[i] = q_vec[i]
    end
  end

  # no need to do any disassembly for DG
  writeQ(mesh, sbp, eqn, opts)

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
### Utils.array3DTo1D

  This function takes the 3D array of variables and
  reassembles it into the vector form by summing into the output array.
  Note that
  this is a reduction operation and zeros res_vec before performing the
  operation, unless zero_resvec is set to false

  This is a mid level function, and does the right thing regardless of
  equation dimension

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * res_array: 3D array, same shape as `eqn.res`

  **Inputs/Outputs**

   * res_vec: vector, same shape as `eqn.res_vec`, to sum into

  **Keyword Arguments**

   * zero_resvec: if true (default), zeros `res_vec` before summing into it


"""->
# mid level function (although it doesn't need Tdim)
function array3DTo1D(mesh::AbstractCGMesh{Tmsh},
       sbp, eqn::AbstractSolutionData{Tsol}, opts,
       res_arr::Abstract3DArray, res_vec::AbstractVector;
       zero_resvec=true) where {Tmsh, Tsol}
# arr is the array to be assembled into res_vec

#  println("in array3DTo1D")

  if zero_resvec
    @assert pointer(res_vec) != pointer(res_arr)
    fill!(res_vec, 0)
  end

  # sum into output
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      @simd for k=1:size(res_arr, 1)  # loop over dofs on the node
        dofnum_k = mesh.dofs[k, j, i]
        res_vec[dofnum_k] += res_arr[k,j,i]
      end
    end
  end

  return nothing
end

function array3DTo1D(mesh::AbstractDGMesh{Tmsh},
       sbp, eqn::AbstractSolutionData{Tsol}, opts,
       res_arr::Abstract3DArray, res_vec::AbstractVector;
       zero_resvec=true) where {Tmsh, Tsol}

  # we assume the memory layouts of q_arr and q_vec are the same
  if pointer(res_arr) != pointer(res_vec)
    if zero_resvec
      @simd for i = 1:length(res_vec)
        res_vec[i] = res_arr[i]
      end
    else
      @simd for i = 1:length(res_vec)
        res_vec[i] += res_arr[i]
      end
    end

  # else: nothing to do
  end


  return nothing

end


function arrToVecAssign(mesh::AbstractMesh{Tmsh},
       sbp, eqn::AbstractSolutionData{Tsol}, opts,
       arr::Abstract3DArray, dest_vec::AbstractArray{Tres,1},
       zero_resvec=true) where {Tmsh, Tsol, Tres}

  # This was created so a q -> q_vec operation could be performed, but it is sufficiently
  #   generic to operate on other things.
  # It is the inverse function of array1DTo3D

  # arr is the array to be assembled into dest_vec, using an assignment reduction

  # removed zeroing out of dest vec, as all of it will be overwritten below

  if pointer(arr) == pointer(dest_vec)
    return nothing
  end

  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode  # loop over dofs on the node

        # mesh.dofs indexing:
        #   [dof ix on the node, node ix on the el, el ix]
        dofnum_k = mesh.dofs[k, j, i]

        dest_vec[dofnum_k] = arr[k,j,i]

      end
    end
  end

  return nothing

end


"""
  This function removes the complex part of the following fields of `eqn`:
  `q`, `q_vec`, `res`, `res_vec`

  **Inputs**

   * eqn: an `AbstractSolutionData`
"""
function removeComplex(eqn::AbstractSolutionData)

  @simd for i=1:length(eqn.q_vec)
    eqn.q_vec[i] = real(eqn.q_vec[i])
    eqn.res_vec[i] = real(eqn.res_vec[i])
  end

  if pointer(eqn.q_vec) != pointer(eqn.q)
    @simd for i=1:length(eqn.q)
      eqn.q[i] = real(eqn.q[i])
    end
  end

  if pointer(eqn.res_vec) != pointer(eqn.res)
    @simd for i=1:length(eqn.res)
      eqn.res[i] = real(eqn.res[i])
    end
  end

  # send and receive buffers
  for i=1:length(eqn.shared_data)
    arr = eqn.shared_data[i].q_send
    @simd for j=1:length(arr)
      arr[j] = real(arr[j])
    end

    arr = eqn.shared_data[i].q_recv
    @simd for j=1:length(arr)
      arr[j] = real(arr[j])
    end
  end

  return nothing
end




# mid level function (although it doesn't need Tdim)
function assembleArray(mesh::AbstractMesh{Tmsh},
       sbp, eqn::AbstractSolutionData{Tsol}, opts,
       arr::Abstract3DArray, res_vec::AbstractArray{Tres,1},
       zero_resvec=true) where {Tmsh, Tsol, Tres}
# arr is the array to be assembled into res_vec, using an assignment reduction
# the length of res_vec is mesh.numDof/mesh.numDofPerNode, only the last
# dof on the node is placed into res_vec

#  println("in array3DTo1D")

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
function calcNorm(eqn::AbstractSolutionData, res_vec::AbstractArray{T}; strongres=false, globalnrm=true) where T
# calculates the norm of a vector using the mass matrix

  val = zero(real(res_vec[1]))

  #TODO: make this better for complex number: conj(z)Hz
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

"""
  This function calculate the SBP approximation to the L2 inner product of
  two vectors of length mesh.numDof.

  L2_product = u^T H conj(v)

  Note that this function does not take a square root like [`calcNorm`](@ref)
  does.

  **Inputs**

   * eqn: AbstractSolutionData
   * u: first vector, of length mesh.numDof
   * v: second vector, of length mesh.numDof.  If complex, this vector gets
        conjugated

  **Outputs**

   * val: the value of the inner product

  **Keyword Arguments**

   * globalnrm: if true, computes the norm of the entire (parallel) vector,
                if false, computes the norm of only the local part
"""
function calcL2InnerProduct(eqn::AbstractSolutionData, u::AbstractArray{T}, v::AbstractArray{T2}; globalnrm=true) where {T, T2}

  # TODO: generalize this to 1, mesh.numDofPerNode vectors
  @assert length(u) == length(eqn.M)
  @assert length(v) == length(eqn.M)

  val = zero(real(promote_type(T, T2)))

  for i=1:length(u)
    val += real( u[i] * eqn.M[i] * conj(v[i]) )
  end

  eqn.params.time.t_allreduce = @elapsed if globalnrm
    val = MPI.Allreduce(val, MPI.SUM, eqn.comm)
  end

  return val
end


"""
  This function computes the Euclidian norm of a vector where each MPI
  process owns part of the vector

  Inputs:
    comm: an MPI communicator
    vec: the local part of the vector

  Outputs:
    val: the Euclidian norm of the entire vector

  Note that unlike calcNorm, the time spent in the Allreduce is not logged
  for this function.

"""
function calcEuclidianNorm(comm::MPI.Comm, vec::AbstractVector{T}) where T

  Tnorm = real(T)
  val = zero(Tnorm)
  for i=1:length(vec)
    val += real(conj(vec[i])*vec[i])
  end

  val = MPI.Allreduce(val, MPI.SUM, comm)

  return sqrt(val)
end



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
function calcMeshH(mesh::AbstractMesh{Tmsh}, sbp,  eqn, opts) where Tmsh
  jac_3d = reshape(mesh.jac, 1, mesh.numNodesPerElement, mesh.numEl)
  jac_vec = zeros(Tmsh, mesh.numNodes)
  assembleArray(mesh, sbp, eqn, opts, jac_3d, jac_vec)
  dim = mesh.dim
  # scale by the minimum distance between nodes on a reference element
  # this is a bit of an assumption, because for distorted elements this
  # might not be entirely accurate
  h_avg = sum(1./(jac_vec.^(1/dim)))
  h_avg = MPI.Allreduce(h_avg, MPI.SUM, mesh.comm)
  # caution: overflow
  numnodes = MPI.Allreduce(length(jac_vec), MPI.SUM, mesh.comm)
  h_avg *= mesh.min_node_dist/numnodes
  return h_avg
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
mutable struct Timings
  # timings
  t_volume::Float64  # time for volume integrals
  t_face::Float64 # time for surface integrals (interior)
  t_source::Float64  # time spent doing source term
  t_shock::Float64   # time spent doing shock capturing
  t_sharedface::Float64  # time for shared face integrals
  t_bndry::Float64  # time spent doing boundary integrals
  t_dataprep::Float64  # time spent preparing data
  t_stab::Float64  # time spent adding stabilization
  t_send::Float64  # time spent sending data

  t_volume_diff::Float64  # time for volume integrals
  t_face_diff::Float64 # time for surface integrals (interior)
  t_shock::Float64  # time for shock capturing
  t_source_diff::Float64  # time spent doing source term
  t_sharedface_diff::Float64  # time for shared face integrals
  t_bndry_diff::Float64  # time spent doing boundary integrals
  t_dataprep_diff::Float64  # time spent preparing data
  t_stab_diff::Float64  # time spent adding stabilization


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
  t_timemarch::Float64 # time spent in time marching loop
  t_callback::Float64  # time spent performing callbacks
  t_nlsolve::Float64  # time spent in the nonlinear solver
  t_meshinit::Float64  # time spent initializing the mesh object
  t_barrier::Float64  # time spent in MPI_Barrier
  t_barrier2::Float64
  t_barrier3::Float64
  t_barriers::Array{Float64, 1}

  function Timings()
    nbarriers = 7
    barriers = zeros(Float64, nbarriers)
    return new(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, barriers)
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
  vals = Array{Float64}(nvals)
  val_names = Array{String}(nvals)

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

function calcBCNormal_revm(params::AbstractParamType{2}, dxidx::AbstractMatrix,
                           nrm::AbstractVector, nrm2_bar::AbstractVector,
                           dxidx_bar::AbstractMatrix)

  dxidx_bar[1,1] += nrm2_bar[1]*nrm[1]
  dxidx_bar[2,1] += nrm2_bar[1]*nrm[2]
  dxidx_bar[1,2] += nrm2_bar[2]*nrm[1]
  dxidx_bar[2,2] += nrm2_bar[2]*nrm[2]

  return nothing
end

function calcBCNormal_revm(params::AbstractParamType{3}, dxidx::AbstractMatrix,
                           nrm::AbstractVector, nrm2_bar::AbstractVector, 
                           dxidx_bar::AbstractMatrix)

  dxidx_bar[1,1] += nrm2_bar[1]*nrm[1]
  dxidx_bar[2,1] += nrm2_bar[1]*nrm[2]
  dxidx_bar[3,1] += nrm2_bar[1]*nrm[3]
  dxidx_bar[1,2] += nrm2_bar[2]*nrm[1]
  dxidx_bar[2,2] += nrm2_bar[2]*nrm[2]
  dxidx_bar[3,2] += nrm2_bar[2]*nrm[3]
  dxidx_bar[1,3] += nrm2_bar[3]*nrm[1]
  dxidx_bar[2,3] += nrm2_bar[3]*nrm[2]
  dxidx_bar[3,3] += nrm2_bar[3]*nrm[3]

  return nothing
end

#------------------------------------------------------------------------------
# permutation functions
#------------------------------------------------------------------------------
"""
  Permute the rows of A according to the permvec, storing the result in B
  The permvec contains the source indices for each entry in B, ie.
  B[i] comes from A[permvec[i]].  This is consistent with the mathematical
  definition of a permutation that pre-multiplication by a permutation 
  matrix (obtained from permMatrix) is a row permutation, ie.
  B = P*A

  Aliasing: no aliasing allowed
"""
function applyPermRow(permvec::AbstractVector, A::AbstractMatrix,
                            B::AbstractMatrix)

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
function applyPermRowInplace(permvec::AbstractVector, A::AbstractMatrix,
                            B::AbstractMatrix)

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
function applyPermColumn(permvec::AbstractVector, A::AbstractMatrix,
                            B::AbstractMatrix)

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
function applyPermColumnInplace(permvec::AbstractVector,
                            A::AbstractMatrix, B::AbstractMatrix)

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

"""
  This function zeros out an array, and should be faster than fill! (branch
  free)

  Inputs/Outputs:
    x: an array
"""
@inline function fastzero!(x::AbstractArray)

  @inbounds @simd for i=1:length(x)
    x[i] = 0
  end

  return nothing
end

"""
  This function scales an array by a constant, and should be faster than scale!
  because it is branch free
"""
@inline function fastscale!(x::AbstractArray, c::Number)

  @inbounds @simd for i=1:length(x)
    x[i] *= c
  end

  return nothing
end



@doc """
### Utils.verbose1

  This macro introduces an if statement that causes the expression to be 
  executed only if the variable verbose is greater than or equal to 1.  
  verbose must exist in the scope of the caller

"""->
macro verbose1(ex)
  return quote
#    println("myrank = ", esc(myrank))
    if $(esc(:(verbose >= 1)))
      $(esc(ex))
    end
  end
end

@doc """
### Utils.verbose2

  This macro introduces an if statement that causes the expression to be 
  executed only if the variable verbose is greater than or equal to 2.  
  verbose must exist in the scope of the caller

"""->
macro verbose2(ex)
  return quote
#    println("myrank = ", esc(myrank))
    if $(esc(:(verbose >= 2)))
      $(esc(ex))
    end
  end
end

@doc """
### Utils.verbose3

  This macro introduces an if statement that causes the expression to be 
  executed only if the variable verbose is greater than or equal to 3.  
  verbose must exist in the scope of the caller

"""->
macro verbose3(ex)
  return quote
#    println("myrank = ", esc(myrank))
    if $(esc(:(verbose >= 3)))
      $(esc(ex))
    end
  end
end

@doc """
### Utils.verbose4

  This macro introduces an if statement that causes the expression to be 
  executed only if the variable verbose is greater than or equal to 4.  
  verbose must exist in the scope of the caller

"""->
macro verbose4(ex)
  return quote
#    println("myrank = ", esc(myrank))
    if $(esc(:(verbose >= 4)))
      $(esc(ex))
    end
  end
end

@doc """
### Utils.verbose5

  This macro introduces an if statement that causes the expression to be 
  executed only if the variable verbose is greater than or equal to 5.  
  verbose must exist in the scope of the caller

"""->
macro verbose5(ex)
  return quote
#    println("myrank = ", esc(myrank))
    if $(esc(:(verbose >= 5)))
      $(esc(ex))
    end
  end
end

"""
  This macro enables easy unpacking of the members of a type to local variables
  of the same name.  For example:

  ```
    @unpack foo a b c
  ```

  expands into

  ```
    a = foo.a; b = foo.b; c = foo.c
  ```

  **Inputs**

   * obj: the name of a variable (must be a type with fields)
   * varnames...: 1 or more names to unpack.  An error is thrown if a name is
                  duplicated, presuming that is not what the user intended
"""
macro unpack(obj, varnames...)

  @assert length(varnames) > 0

  # check varnames are unique
  for i=1:length(varnames)
    for j=1:(i-1)
      if varnames[i] == varnames[j]
        error("attempting to unpack duplicate variable name: $(varnames[i])")
      end
    end
  end

  # build the expression
  ex = :($(varnames[1]) = $(obj).$(varnames[1]))
  for i=2:length(varnames)
    ex = :($ex; $(varnames[i])= $(obj).$(varnames[i]))
  end

  return esc(ex)
end

"""
  Given a list of variables, prints them out, one per line

  Ex.  The code

  ```
    a = 2
    b = 3
    @printit a b
  ```
    will print
  ```
    a = 2
    b = 3
  ```

  Variables can also be more general expressions:

  ```
    a = [1 2 3]
    @printit a[1:2]
  ```

  will print

  ```
    a[1:2] = [1, 2]
  ```
"""
macro printit(v...)

  ex = :()  # empty expression
  for s in v
    str = string(s)
#    println("s = ", s)
    ex = :($ex, println($str, " = ", $(esc(s))); )
  end

  return ex
end


"""
  This function asserts that every array field of the given object is a
  unique array.  Note that this is less strict than a non-aliasing check,
  becuase views of an array are not the same as the array.
"""
function assertArraysUnique(obj::T) where {T}
# verify a struct does not have any fields that alias each other

  fnames = fieldnames(T)

  for i=1:length(fnames)
    f_i = getfield(obj, fnames[i])
    if typeof(f_i) <: AbstractArray
      for j=(i+1):length(fnames)
        f_j = getfield(obj, fnames[j])
        if typeof(f_j) <: AbstractArray
          @assert !(f_i === f_j)
        end
      end
    end
  end

  return nothing
end

"""
  Asserts that all fields of the given object have declared field types that
  are concrete.  Note that this only considers the types of the field as
  written in the type declaration, not the runtime type of the fields.

  **Inputs**

   * obj: an object
"""
function assertFieldsConcrete(obj::T) where {T}

  for i=1:nfields(T)
    type_i = fieldtype(T, i)
    @assert isleaftype(type_i)
  end

  return nothing
end


end  # end module
