# residual_evaluation.jl: some functions used by NonlinearSolvers relating to 
# residual evaluation

export calcResidual
@doc """
### NonlinearSolvers.calcResidual

  This function takes the eqn object with the solution varaibles stored in
    q_vec, scatters them into q, evaluates the residual, and then gathers the 
    residual values into res_vec.

    Effectively, this is a wrapper around the physics module residual evaluation
      function (which performs eqn.q -> eqn.res) that performs eqn.q_vec ->
      eqn.res_vec.

    The norm of the strong residual (using the SBP norm) is calculated and returned
    (even though the weak form residual is stores in eq.res_vec).

    Inputs:
      mesh:  an AbstractMesh object
      sbp:  an SBP operator
      eqn:  an AbstractSolutionData object
      opts:  options dictonary
      func: residual evaluation function

    Outputs:
      res_0_norm:  norm of residual

    Aliasing restrictions: none
"""->
function calcResidual(mesh, sbp, eqn, opts, func)
# calculate the residual and its norm

  disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)

  if opts["parallel_type"] == 2
    exchangeElementData(mesh, opts, eqn.q, eqn.q_face_send, eqn.q_face_recv, eqn.params.f)
  end

  func(mesh, sbp, eqn, opts)

  assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec, assemble_edgeres=false)

  res_0_norm = calcNorm(eqn, eqn.res_vec, strongres=true)
  res_norm_global = MPI.Allreduce(res_0_norm*res_0_norm, MPI.SUM, mesh.comm)
#  println("residual norm = ", res_0_norm)

 return sqrt(res_norm_global)
end


@doc """
### NonlinearSolvers.assembleResidual

  This function takes the residual in eqn.res and performs an additive reduction
  into res_vec.
  This function wraps assembleSolution.

  Inputs:
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    res_vec: residual vector to put the residual into

  Outputs:
    none

  Aliasing restrictions: none

"""->
function assembleResidual{T}(mesh, sbp, eqn, opts, res_vec::AbstractArray{T, 1}; 
                             assemble_edgeres=true, zero_resvec=true)
# assembles all of the residuals into res_vec
# no aliaising concerns

  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, res_vec)

  if assemble_edgeres

    for i=1:size(eqn.res_edge, 4)
      eqn.assembleSolution(mesh, sbp, eqn, opts, sview(eqn.res_edge, :, :, :, i),
                           res_vec, zero_resvec=zero_resvec)
    end
  end

  return nothing
end


@doc """
### NonlinearSolvers.disassembleSolution

  This function performs the scatter q_vec -> eqn.q

  Inputs:
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    q_vec: vector containing solution variables 

  Outputs:
    none

  Aliasing Restrictions: none

"""->
function disassembleSolution{T}(mesh, sbp, eqn, opts, q_vec::AbstractArray{T, 1})
# scatters the q_vec to the 3d array eqn.q
# no aliasing concerns here
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, q_vec)

  return nothing
end
