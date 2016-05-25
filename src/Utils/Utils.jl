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
include("logging.jl")
export disassembleSolution, writeQ, assembleSolution, assembleArray, sview
export initMPIStructures, exchangeFaceData, verifyCommunication, getSendData
export exchangeElementData
export @mpi_master, @time_all, print_time_all
export Timings, write_timings
export sharedFaceLogging

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
    return new(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, barriers)
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



end  # end module

  
