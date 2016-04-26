# Adjoint for Euler Equations

push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
using NonlinearSolvers

export calcAdjoint

@doc """
### EulerEquationMod.calcAdjoint

Calculated the adjoint vector for a single functional

**Inputs**



**Outputs**

*  None

"""->

function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGmesh{Tmsh},
	                sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts,
	                functional_number, adjoint_vec::Array{Tsol, 1})



  return nothing
end