# helper functions for interface2.jl.  These functions are not interface
# functions themselves, hence this file

"""
  This function starts parallel communication, if necessary, for
  functionals.  This function also unpacks `eqn.q_vec` into `eqn.q`
  if communication is started.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * func: the functional
   * start_comm: bool specifying if communication must be started.
                 Communication may also be started when `start_comm` is false,
                 if the parallel data of `eqn.shared_data` is different.

  **Outputs**

   * start_comm_q: bool describing if communication was started or not
"""
function startCommunicationFunctional(mesh, sbp, eqn, opts,
                                      func::AbstractFunctional, start_comm::Bool)
  # start communication
  pdata = getParallelDataString(getParallelData(func))
  start_comm_q = start_comm || (pdata != getParallelData(eqn.shared_data) &&
                                pdata != "none")

  if start_comm 
    array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
    eqn.params.time.t_send += @elapsed if eqn.commsize > 1 && start_comm_q
      setParallelData(eqn.shared_data, pdata)
      startSolutionExchange(mesh, sbp, eqn, opts)
    end
  end

  return start_comm_q
end


