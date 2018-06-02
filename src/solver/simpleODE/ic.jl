# ic.jl
# Needed to initialize a problem.

function ICConstant(mesh::AbstractMesh, sbp::AbstractSBP, eqn::SimpleODEData, opts,
                    u0::AbstractArray)

  for i=1:length(u0)
    u0[i] = 2
  end

  return nothing
end

@doc """
### SimpleODEMod.ICx2plust4

Computes the initial conditions for the state variable
    u = x^2 + t^4

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Simple ODE equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->
function ICx2plust4(mesh::AbstractMesh{Tmsh}, 
sbp::AbstractSBP{Tsbp}, eqn::SimpleODEData{Tsol, Tres, 2}, 
opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}

  t = 0.0

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
#       x = mesh.coords[1, j, i]
#       y = mesh.coords[2, j, i]
#       u0[dofnums_j] = x^2 + t^2
      u0[dofnums_j] = calc_x2_4t3(mesh.coords[:, j, i], eqn.params, t)
    end
  end

  return nothing
end # end function ICx2plust4

@doc """
### SimpleODEMod.ICx2plust3

Computes the initial conditions for the state variable
    u = x^2 + t^3

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Simple ODE equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->
function ICx2plust3(mesh::AbstractMesh{Tmsh}, 
sbp::AbstractSBP{Tsbp}, eqn::SimpleODEData{Tsol, Tres, 2}, 
opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}

  t = 0.0

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
#       x = mesh.coords[1, j, i]
#       y = mesh.coords[2, j, i]
#       u0[dofnums_j] = x^2 + t^2
      u0[dofnums_j] = calc_x2_3t2(mesh.coords[:, j, i], eqn.params, t)
    end
  end

  return nothing
end # end function ICx2plust3

@doc """
### SimpleODEMod.ICx2plust2

Computes the initial conditions for the state variable
    u = x^2 + t^2

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Simple ODE equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICx2plust2(mesh::AbstractMesh{Tmsh}, 
sbp::AbstractSBP{Tsbp}, eqn::SimpleODEData{Tsol, Tres, 2}, 
opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}

  t = 0.0

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
#       x = mesh.coords[1, j, i]
#       y = mesh.coords[2, j, i]
#       u0[dofnums_j] = x^2 + t^2
      u0[dofnums_j] = calc_x2_t2(mesh.coords[:, j, i], eqn.params, t)
    end
  end

  return nothing
end # end function ICx2plust2

@doc """
### SimpleODEMod.ICx2plust2

Computes the initial conditions for the state variable
    u = x^2 + t^2

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Simple ODE equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICallzero(mesh::AbstractMesh{Tmsh}, 
sbp::AbstractSBP{Tsbp}, eqn::SimpleODEData{Tsol, Tres, 2}, 
opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}
  t = 0.0

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
#       x = mesh.coords[1, j, i]
#       y = mesh.coords[2, j, i]
#       u0[dofnums_j] = x^2 + t^2
      u0[dofnums_j] = 0.0
    end
  end

  return nothing
end # end function ICallzero

function ICallOnes(mesh::AbstractMesh{Tmsh}, 
sbp::AbstractSBP{Tsbp}, eqn::SimpleODEData{Tsol, Tres, 2}, 
opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}
  t = 0.0

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
#       x = mesh.coords[1, j, i]
#       y = mesh.coords[2, j, i]
#       u0[dofnums_j] = x^2 + t^2
      u0[dofnums_j] = 1.0
    end
  end

  return nothing
end # end function ICallOnes

# TODO: rename ICallzero to ICallZeros

global const ICDict = Dict{Any, Function}(
"ICConstant" => ICConstant,
"ICx2plust4" => ICx2plust4,
"ICx2plust3" => ICx2plust3,
"ICx2plust2" => ICx2plust2,
"ICallzero" => ICallzero,
"ICallOnes" => ICallOnes,
)
