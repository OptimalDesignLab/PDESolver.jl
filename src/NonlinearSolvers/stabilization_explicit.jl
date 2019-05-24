"""
  calcStabilizedQUpdate!:
    Provides a stabilized version of the imaginary component of q_vec.

  Input:
    mesh, sbp, eqn, opts: standard
    stab_A: DiagJac type
    stab_assembler: AssembleDiagJacData type
    clipJacData: type ClipJacData from clipJac; storage for operations required for clipJac routines
    t: Time. Make sure to use treal here.

  Old input:
    q_vec: the q_vec to be stabilized. This should be the actual q_vec - only the imaginary component
          of q_vec is used here, because we don't want to affect the solution in any way,
          just the direct sensitivity. This will not be modified here

    Now q_vec is obtained straight from eqn.q_vec

  In/Output:
    Bv: just an array that is prealloc'd. size of q_vec. Real.
        In the formulation, this is B*v = B*imag(q_vec)
"""
function calcStabilizedQUpdate!(mesh, sbp, eqn, opts, 
                                stab_A, stab_assembler, clipJacData,
                                t, Bv, tmp_imag)
                         # q_vec::AbstractArray{Tsol, 1},
                         # Bv::AbstractArray{Tsol, 1})
                         # stab_A::DiagJac,
                         # stab_assembler::AssembleDiagJacData,
  MatZeroEntries(stab_A)

  # Taking the imaginary component of q out -> Storing it in tmp_imag
  fill!(tmp_imag, 0.0)
  for j = 1:length(eqn.q_vec)
    tmp_imag[j] = imag(eqn.q_vec[j])

    # commenting out did not fix negative density
    # eqn.q_vec[j] = complex(real(eqn.q_vec[j]), 0.0)   

    # for array assignment, julia will elevate the real to complex, so this is cleaner
    eqn.q_vec[j] = real(eqn.q_vec[j])     
  end

  evalJacobianStrong(mesh, sbp, eqn, opts, stab_assembler, t)     # calcs Qx*Ax + Qy*Ay     # TODO: is stab_assembler.A complex or real

  filterDiagJac(mesh, opts, real(tmp_imag), clipJacData, stab_A, eigs_to_remove="neg")        # stab_A is now B in the derivation

  # Putting the imaginary component back
  for j = 1:length(eqn.q_vec)
    # println(" j: $j, real(eqn.q_vec[j]): ", real(eqn.q_vec[j]), " tmp_imag[j]: ", tmp_imag[j])

    # commenting out did not fix negative density
    eqn.q_vec[j] = complex(real(eqn.q_vec[j]), tmp_imag[j])
  end

  # Bv fill! to 0's is not required here. See diagMatVec code. Bv is assigned straight into, no += or anything

  # TODO: in filterDiagJac, ublock is taken as q_vec[mesh.dofs[i,j,k]]. 
  #       Is this equivalent to the diagMatVec below?

  # does Bv = B*imag(q_vec)
  diagMatVec(stab_A, mesh, imag(eqn.q_vec), Bv)     # Prof H thinks stab_A needs to be real       # TODO TODO imag(q_vec) needs to have cplx perturbation applied to it???
  # println(" vecnorm(Bv): ", vecnorm(Bv))

  # Note: here is where you can scale Bv to increase its effect, say by a factor of 4

  # application of Bv to q_vec happens outside of this function.

  # The mass matrix needs to be applied to Bv, as it is part of the residual.
  #   see pde_post_func in rk4.jl
  # for j = 1:length(Bv)
    # Bv[j] = eqn.Minv[j]*Bv[j]
  # end

  # Testing Bv


  return nothing

end
