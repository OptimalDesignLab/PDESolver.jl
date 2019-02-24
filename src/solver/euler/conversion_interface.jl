"""
  Abstract type that defines the conversion between conservative variables
  and some other set of variables, typically the entropy variables.

  Required functions:

  ```
    convertToEntropy
    convertToConservative
    getA0
    getA0inv
  ```

  See those functions docstrings for details
"""
abstract type AbstractVariables end

"""
  Convert to entropy variables from conservative variables

  **Inputs**

   * obj: the [`AbstractVariables`](@ref) object specifying what variables
          to convert to
   * params: AbstractParamType object
   * qc: the conservative variables at the node, length `numDofPerNode`

  **Inputs/Outputs**

   * qe: vector to be overwritten with the entropy variables, same length
         as `qc`

  Aliasing: `qc` and `qe` may alias.  Implementations must ensure correct
            behavior in this case
"""
function convertToEntropy(obj::AbstractVariables, params::AbstractParamType,
                          qc::AbstractVector, qe::AbstractVector)

  error("abstract method for convertToEntropy() reached")
end

"""
  Convert to conservative variables from entropy variables

  **Inputs**

   * obj: the [`AbstractVariables`](@ref) object specifying what variables
          to convert to
   * params: AbstractParamType object
   * qe: the entropy variables at the node, length `numDofPerNode`

  **Inputs/Outputs**

   * qc: vector to be overwritten with the conservative variables, same length
         as `qe`

  Aliasing: `qc` and `qe` may alias.  Implementations must ensure correct
            behavior in this case
"""
function convertToConservative(obj::AbstractVariables, params::AbstractParamType,
                          qe::AbstractVector, qc::AbstractVector)

  error("abstract method for convertToConservative() reached")
end


"""
  Computes the Jacobian dq/dw, where `q` are the conservative variables and
  `w` are the entropy variables`.

  **Inputs**

   * obj: an [`AbstractVariables`](@ref) object
   * params: AbstractParamType object
   * qc: vector of conservative variables specifying the state at which to
         compute `A0`, length `numDofPerNode`

  **Inputs/Outputs**

   * A0: matrix to be overwritten with `dq/dw, `numDofPerNode` x `numDofPerNode`
"""
function getA0(obj::AbstractVariables, params::AbstractParamType,
               qc::AbstractVector, A0::AbstractMatrix)

  error("abstract method for getA0() reached")
end


"""
  Computes the Jacobian dw/dq = inv(dq/dw), where `q` are the conservative
  variables and `w` are the entropy variables`.

  **Inputs**

   * obj: an [`AbstractVariables`](@ref) object
   * params: AbstractParamType object
   * qc: vector of conservative variables specifying the state at which to
         compute `A0`, length `numDofPerNode`

  **Inputs/Outputs**

   * A0inv: matrix to be overwritten with `dw/dq, `numDofPerNode` x
           `numDofPerNode`
"""
function getA0inv(obj::AbstractVariables, params::AbstractParamType,
                  qc::AbstractVector, A0inv::AbstractMatrix)

  error("abstract method for getA0inv() reached")
end

#TODO: add: entropy function, entropy flux


#------------------------------------------------------------------------------
# Variables used by Ismail and Roe

"""
  Variables used by Ismail and Roe.  See "Affordable, entropy-consistent Euler
  Flux Functions II: Entropy Production at Shocks"
"""
struct IRVariables <: AbstractVariables
end

function convertToEntropy(obj::IRVariables, params::AbstractParamType,
                          qc::AbstractVector, qe::AbstractVector)

  convertToIR_(params, qc, qe)
end


function convertToConservative(obj::IRVariables, params::AbstractParamType,
                          qe::AbstractVector, qc::AbstractVector)

  convertToConservativeFromIR(params, qe, qc)
end

function getA0(obj::IRVariables, params::AbstractParamType, qc::AbstractVector,
               A0::AbstractMatrix)

  getIRA0(params, qc, A0)
end

function getA0inv(obj::IRVariables, params::AbstractParamType, qc::AbstractVector,
                  A0inv::AbstractMatrix)

  getIRA0inv(params, qc, A0inv)
end

#------------------------------------------------------------------------------
# Hughes entropy variables variables

"""
  Variables used by Hughes.  See "A New Finite Element Formulation For
  Computational Fluid Dynamics: Part I"
"""
struct HughesVariables <: AbstractVariables
end

function convertToEntropy(obj::HughesVariables, params::AbstractParamType,
                          qc::AbstractVector, qe::AbstractVector)

  convertToEntropy_(params, qc, qe)
end


function convertToConservative(obj::HughesVariables, params::AbstractParamType,
                          qe::AbstractVector, qc::AbstractVector)

  #TODO: cleanup the naming scheme
  convertToConservative(params, qe, qc)
end

function getA0(obj::HughesVariables, params::AbstractParamType, qc::AbstractVector,
               A0::AbstractMatrix)

  #TODO: cleanup the naming scheme
  calcA0(params, qc, A0)
end

function getA0inv(obj::HughesVariables, params::AbstractParamType, qc::AbstractVector,
                  A0inv::AbstractMatrix)

  calcA0inv(params, qc, A0inv)
end


#------------------------------------------------------------------------------
# Conservative variables 
# Useful for testing to make conservative = entropy

"""
  Trick the solver into thinking the conservative and entropy variables are
  the same.
"""
struct ConservativeVariables <: AbstractVariables
end

function convertToEntropy(obj::ConservativeVariables, params::AbstractParamType,
                          qc::AbstractVector, qe::AbstractVector)

  for i=1:length(qc)
    qe[i] = qc[i]
  end
end


function convertToConservative(obj::ConservativeVariables, params::AbstractParamType,
                          qe::AbstractVector, qc::AbstractVector)

  for i=1:length(qc)
    qc[i] = qe[i]
  end
end

function getA0(obj::ConservativeVariables, params::AbstractParamType, qc::AbstractVector,
               A0::AbstractMatrix)

  fill!(A0, 0)
  for i=1:size(A0, 1)
    A0[i, i] = 1
  end
end

function getA0inv(obj::ConservativeVariables, params::AbstractParamType, qc::AbstractVector,
                  A0inv::AbstractMatrix)

  fill!(A0inv, 0)
  for i=1:size(A0inv, 1)
    A0inv[i, i] = 1
  end
end



