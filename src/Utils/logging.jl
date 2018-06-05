# functions to log data

@doc """
### Utils.sharedFaceLogging

  This function writes files for the function that evalutes the shared face
  integrals.  This function should only be called in debug mode.
  
  If opts["writeqface"] is true,  writes the q values at the face to a file 
  q_sharedface_i_myrank.dat, where i is the peer process index (not MPI rank)
  and myrank is the MPI rank of the current process. 

  If opts["write_fluxface"] is true, it writes eqn.flux_sharedface to a file
  flux_sharedface_i_myrank.dat

  Inputs:
    mesh
    sbp
    eqn: and AbstractSolutionData
    opts: options dictonary
    data: SharedFaceData object
    qL_arr: an array holding solution values at face nodes on the 
            left side of the interface,  numDofPerNode x numfacenodes x 
            numsharedfaces on this partition boundary
    qR_arr: solution values at face nodes on the right side of the interface.
            same shape as qL_arr

  also, the eqn.flux_shared is used to write the face flux.  It is the same
  shape as qL_arr

  Aliasing restrictions: qL_arr, qR_arr, and eqn.flux_sharedface must not alias.
"""->
function sharedFaceLogging(mesh, sbp, eqn::AbstractSolutionData{Tsol}, 
                           opts, data::SharedFaceData, qL_arr, qR_arr) where Tsol

  if opts["writeqface"]
    myrank = data.myrank
    i = data.peeridx
    tmp_arr = zeros(Tsol, mesh.numDofPerNode, 2, mesh.numNodesPerFace, mesh.peer_face_counts[i])
    for j = 1:mesh.peer_face_counts[i]
      for k=1:mesh.numNodesPerFace
        tmp_arr[:, 1, k, j] = qL_arr[:, k, j]
        tmp_arr[:, 2, k, j] = qR_arr[:, k, j]
      end
    end
    println(eqn.params.f, "q_sharedface $i = \n", tmp_arr)
    fname = string("qsharedface_", i, "_", myrank, ".dat")
    writedlm(fname, tmp_arr)

  end  # end if

  if opts["write_fluxface"]
    for i=1:mesh.npeers
      fname = string("fluxsharedface_", i, "_", myrank, ".dat")
      writedlm(fname, eqn.flux_sharedface[i])
    end
  end

  flush(eqn.params.f)
  return nothing
end


