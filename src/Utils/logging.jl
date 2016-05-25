@doc """
### Utils.sharedFaceLogging

  This function writes files for the function that evalutes the shared face
  integrals.  This function should only be called in debug mode

  Inputs:
    mesh
    sbp
    eqn: and AbstractSolutionData
    opts: options dictonary
    qL_arr: an array of arrays holding solution values at face nodes on the 
            left side of the interface, ie. there are mesh.npeer arrays, and 
            each array is numDofPerNode x numfacenodes x numsharedfaces on this
            partition boundary
    qR_arr: solution values at face nodes on the right side of the interface.
            same shape as qL_arr

  also, the eqn.flux_shared is used to write the face flux.  It is the same
  shape as qL_arr

  Aliasing restrictions: qL_arr, qR_arr, and eqn.flux_sharedface must not alias.
"""->
function sharedFaceLogging{Tsol}(mesh, sbp, eqn::AbstractSolutionData{Tsol}, 
                                 opts, qL_arr, qR_arr)

  if opts["writeqface"]
    myrank = mesh.myrank
    for i=1:mesh.npeers
      tmp_arr = zeros(Tsol, mesh.numDofPerNode, 2, mesh.numNodesPerFace, mesh.peer_face_counts[i])
      qL_arr_i = qL_arr[i]
      qR_arr_i = qR_arr[i]
      for j = 1:mesh.peer_face_counts[i]
        for k=1:mesh.numNodesPerFace
          tmp_arr[:, 1, k, j] = qL_arr_i[:, k, j]
          tmp_arr[:, 2, k, j] = qR_arr_i[:, k, j]
        end
      end
      println(eqn.params.f, "q_sharedface $i = \n", tmp_arr)
      fname = string("qsharedface_", i, "_", myrank, ".dat")
      writedlm(fname, tmp_arr)
    end  # end loop over peers

  end  # end if

  if opts["write_fluxface"]
    for i=1:mesh.npeers
      fname = string("fluxsharedface_", i, "_", myrank, ".dat")
      writedlm(fname, eqn.flux_sharedface[i])
    end
  end

  flush(params.f)
  return nothing
end


