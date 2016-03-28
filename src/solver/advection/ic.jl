# ic.jl
# Needed to initialize a problem.

@doc """
### AdvectionEquationMod.ICx5plusy5

Computes the initial conditions for the state variable
    u = x^5 + y^5

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICx5plusy5{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
  	  u0[dofnums_j] = x^5 + y^5
  	end
  end

  return nothing
end # end function ICx5plusy5

@doc """
### AdvectionEquationMod.ICexp_xplusy

Computes the initial conditions for the state variable
    u = exp(x + y)

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICexp_xplusy{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
  	  u0[dofnums_j] = exp(x + y)
  	end
  end

  return nothing
end # end function exp_xplusy


function ICsinwave{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y
  	  u0[dofnums_j] = calc_sinwave(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy


function ICsinwavey{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_sinwavey(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy


function ICsinwavey_pert{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_sinwavey_pert(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy

function ICmms1{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_mms1(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy

function ICx4{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_x4(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy


function ICp1{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_p1(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy

function ICp2{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_p2(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy

function ICp3{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_p3(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy

function ICp4{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_p4(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy


function ICp5{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          alpha_x = eqn.alpha_x
          alpha_y, = eqn.alpha_y

  	  u0[dofnums_j] = calc_p5(mesh.coords[:, j, i], alpha_x, alpha_y, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy



@doc """
### AdvectionEquationMod.ICFile

This function reads a vector from a file on disk and set the solution to it.
The vector must contain the same number of entries as there are degrees of 
freedom in the mesh. 

This function is useful for things like restarting from a checkpoint.
In this case, the file should be the output of writedlm(eqn.q).  The degree 
of freedom number must be the same for both simulation for this to work (the 
file contains no degree of freedom number information).


**Inputs**

*  `mesh`
*  `sbp`
*  `eqn`
*  `opts`

**Inputs/Outputs**

*  `u0`: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICFile{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                operator::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, opts, 
                u0::AbstractVector{Tsol})
# populate u0 with initial values from a disk file
# the file name comes from opts["ICfname"]

  fname = opts["ICfname"]
  vals = readdlm(fname)

  @assert length(vals) == mesh.numDof

  for i=1:mesh.numDof
    u0[i] = vals[i]
  end

end

global const ICDict = Dict{Any, Function}(
"ICx5plusy5" => ICx5plusy5,
"ICexp_xplusy" => ICexp_xplusy,
"ICsinwave" => ICsinwave,
"ICsinwavey" => ICsinwavey,
"ICsinwavey_pert" => ICsinwavey_pert,
"ICFile" => ICFile,
"ICmms1" => ICmms1,
"ICx4" => ICx4,
"ICp1" => ICp1,
"ICp2" => ICp2,
"ICp3" => ICp3,
"ICp4" => ICp4,
"ICp5" => ICp5,
)
