# Abstract interface for diffusion penalties


"""
  Abstract type for the different diffusion penalties that can be expressed
  in the framework of Yan et.al. "Interior Penalties for Summation-by-Parts
  Discretizations of Linear Second-Order Differential Equations".

  Each type must implement a method of:

  ```
   applyPenalty(penalty::AbstractDiffusionPenalty, sbp, sbpface,
                diffusion::AbstractDiffusion, iface::Interface,
                delta_w::AbstractMatrix, theta::AbstractMatrix,
                wL::AbstractMatrix, wR::AbstractMatrix,
                nrm_face::AbstractMatrix,
                alphas::AbstractVector,
                jacL::AbstractVector, jacR::AbstractVector,
                res1L::AbstractMatrix, res1R::AbstractMatrix,
                res2L::AbstractMatrix, res2R::AbstractMatrix
               ) where {Tsol, Tres}


  ```
  specializing the `penalty` argument.  See the docstrings for these
  functions for details.
  
"""
abstract type AbstractDiffusionPenalty end


"""
  Applies the interface penalty for 2nd derivative terms.
  More specifically, it applies the penalty matrix
  for both sides of the face at the same time:

  [res1L   = [T1 T3  [delta_w   and   [res1R   = [T1 T3  [-delta_w
   res2L]     T2 T4]  theta]           res2R]     T2 T4]  theta]



  and similarly for the right element.  Note that delta_w is negated for
  the right element, and this function performs the required negation.

  **Inputs**

   * penalty: the [`AbstractDiffusionPenalty`](@ref) object
   * sbp
   * sbpface
   * diffusion: the [`AbstractDiffusion`](@ref) object
   * iface: `Interface` object
   * delta_w: difference of entropy variables at the face, as compuated by
              [`getFaceVariables`](@ref)
   * theta: Dgk * wL + Dgn * wR, as compuated by `getFaceVariables`
   * wL: entropy variables for the left element, `numDofPerNode` x
         `numNodesPerElement`
   * wR: entropy variables for the right element, same size as `wR`
   * nrm_face: (scaled) normal vectors at the face, `dim` x `numNodesPerFace`
   * alphas: vector of length 2 containing alpha_gk and alpha_gn
   * jacL: mapping jacobian determinant for the left element,
           `numNodesPerElement`
   * jacR: mapping jacobian determinant for the right element, same size as
           `jacR`

  **Inputs/Outputs**
  
   * res1L: `numDofPerNode` x `numNodesPerFace`
   * res1R: same size as above
   * res2L: same size as above
   * res2R: same size as above

  All output arrays are overwritten
"""
function applyPenalty(penalty::AbstractDiffusionPenalty, sbp, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix, theta::AbstractMatrix,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      alphas::AbstractVector,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L::AbstractMatrix, res1R::AbstractMatrix,
                      res2L::AbstractMatrix, res2R::AbstractMatrix)

  error("generic fallback for applyPenalty() reached: did you forget to extend it with a new method for your AbstractDiffusionPenalty?")
end

"""
  **Inputs**

   * penalty: the [`AbstractDiffusionPenalty`](@ref) object
   * sbp
   * sbpface
   * diffusion: the [`AbstractDiffusion`](@ref) object
   * iface: `Interface` object
   * delta_w: difference of entropy variables at the face, as compuated by
              [`getFaceVariables`](@ref)
   * delta_w_dotL: array containing derivative of `delta_w` wrt `qL`,
                   `numDofPerNode` x `numDofPerNode` x `numNodesPerFace`
                   x `numNodesPerElement`
   * delta_w_dotR: derivative of `delta_w` wrt `qR, same shape as `delta_w_dotR`
   * theta: Dgk * wL + Dgn * wR, as compuated by `getFaceVariables`
   * theta_dotL derivative of `theta` wrt `qL`, same shape as `delta_w_dotL`
   * theta_dotR derivative of `theta` wrt `qR`, same shape as `delta_w_dotL`
   * wL: entropy variables for the left element, `numDofPerNode` x
         `numNodesPerElement`
   * wR: entropy variables for the right element, same size as `wR`
   * nrm_face: (scaled) normal vectors at the face, `dim` x `numNodesPerFace`
   * alphas: vector of length 2 containing alpha_gk and alpha_gn
   * jacL: mapping jacobian determinant for the left element,
           `numNodesPerElement`
   * jacR: mapping jacobian determinant for the right element, same size as
           `jacR`

  **Inputs/Outputs**
  
   * res1L_dotL: derivative of `res1L` wrt `qL`, same size as `delta_w_dotL`
   * res1L_dotR: derivative of `res1L` wrt `qR`, same size as above`
   * res1R_dotL: derivative of `res1R` wrt `qL`, same size as `delta_w_dotL`
   * res1R_dotR: derivative of `res1R` wrt `qR`, same size as above`
   * res2L: see `applyPenalty` 
   * res2R: see `applyPenalty
   * res2L_dotL: derivative of `res2L` wrt `qL`, same size as `delta_w_dotL`
   * res2L_dotR: derivative of `res2L` wrt `qR`, same size as above`
   * res2R_dotL: derivative of `res2R` wrt `qL`, same size as `delta_w_dotL`
   * res2R_dotR: derivative of `res2R` wrt `qR`, same size as above`

   All of these array are overwritten

  **Outputs**

   * T4_nonzero: true if the T4 penalty is zero, false otherwise. This is
                 used to determine the sparsity pattern of the Jacobian.
"""
function applyPenalty_diff(penalty::AbstractDiffusion, sbp, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix,
                      delta_w_dotL::Abstract4DArray, delta_w_dotR::Abstract4DArray,
                      theta::AbstractMatrix,
                      theta_dotL::Abstract4DArray, theta_dotR::Abstract4DArray,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      alphas::AbstractVector,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L_dotL::Abstract4DArray, res1L_dotR::Abstract4DArray,
                      res1R_dotL::Abstract4DArray, res1R_dotR::Abstract4DArray,
                      res2L::AbstractMatrix, res2R::AbstractMatrix,
                      res2L_dotL::Abstract4DArray, res2L_dotR::Abstract4DArray,
                      res2R_dotL::Abstract4DArray, res2R_dotR::Abstract4DArray)


  error("generic fallback for applyPenalty_diff() reached: did you forget to extend it with a new method for your AbstractDiffusionPenalty?")
end


#------------------------------------------------------------------------------
# Diffusion penalties

include("br2_penalty.jl")

global const DiffusionPenaltyDict = Dict{String, Type{T} where T <: AbstractDiffusionPenalty}(
"BR2" => BR2Penalty,
)

"""
  Returns a fully constructed [`AbstractDiffusionPenalty`](@ref).

  **Inputs**
  
   * mesh
   * sbp
   * eqn
   * opts
   * name: name of penalty to get, defaults to opts["DiffusionPenalty"]

  **Outputs**

   * obj: `AbstractDiffusionPenalty` object
"""
function getDiffusionPenalty(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                             name::String=opts["DiffusionPenalty"]
                            ) where {Tsol, Tres}

  obj = DiffusionPenaltyDict[name]{Tsol, Tres}(mesh, sbp, opts)
  assertArraysUnique(obj)

  return obj
end

 
