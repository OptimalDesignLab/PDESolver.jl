# ic.jl
# Needed to initialize a problem.

function ICConstant(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AdvectionData, opts,
                    u0::AbstractArray)

  for i=1:length(u0)
    u0[i] =  2
  end

  return nothing
end

@doc """
### AdvectionEquationMod.ICx5plusy5

Computes the initial conditions for the state variable
    u = x^5 + y^5 + z^5; z is omitted it 2d

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICx5plusy5{Tmsh, Tsbp, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol, Tres, 2}, 
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

function ICx5plusy5{Tmsh, Tsbp, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol, Tres, 3}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
  	  u0[dofnums_j] = x^5 + y^5 + z^5
  	end
  end

  return nothing
end # end function ICx5plusy5


@doc """
### AdvectionEquationMod.ICexp_xplusy

Computes the initial conditions for the state variable
    u = exp(x + y + z), z is omitted in 2d

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICexp_xplusy{Tmsh, Tsbp, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol, Tres, 2}, 
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

function ICexp_xplusy{Tmsh, Tsbp, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol, Tres, 3}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
          z = mesh.coords[3,j,i]
  	  u0[dofnums_j] = exp(x + y + z)
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
  	  u0[dofnums_j] = calc_sinwave(mesh.coords[:, j, i], eqn.params, eqn.t)
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
  	  u0[dofnums_j] = calc_sinwavey(mesh.coords[:, j, i], eqn.params, eqn.t)
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
  	  u0[dofnums_j] = calc_sinwavey_pert(mesh.coords[:, j, i], eqn.params, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy


function ICsinwavexy{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = mesh.dofs[1, j, i]
      coords = sview(mesh.coords, :, j, i)
      u0[dofnums_j] = sin(2*pi*coords[1]) + sin(2*pi*coords[2])
      if mesh.dim == 3
        u0[dofnums_j] += sin(2*pi*coords[3])
      end
    end
  end

  return nothing
end


function ICmms1{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  u0[dofnums_j] = calc_mms1(mesh.coords[:, j, i], eqn.params, eqn.t)
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
  	  u0[dofnums_j] = calc_x4(mesh.coords[:, j, i], eqn.params, eqn.t)
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
          alpha_x = eqn.params.alpha_x
          alpha_y, = eqn.params.alpha_y

  	  u0[dofnums_j] = calc_p1(mesh.coords[:, j, i], eqn.params, eqn.t)
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
  	  u0[dofnums_j] = calc_p2(mesh.coords[:, j, i], eqn.params, eqn.t)
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
  	  u0[dofnums_j] = calc_p3(mesh.coords[:, j, i], eqn.params, eqn.t)
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
  	  u0[dofnums_j] = calc_p4(mesh.coords[:, j, i], eqn.params, eqn.t)
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
  	  u0[dofnums_j] = calc_p5(mesh.coords[:, j, i], eqn.params, eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy

function ICexp5xplus4yplus2{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                            sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol},
                            opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp5xplus4yplus2(mesh.coords[:, j, i], eqn.params, eqn.t)
    end
  end

  return nothing
end

function ICexp5xplusy{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                            sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol},
                            opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp5xplusy(mesh.coords[:, j, i], eqn.params,
                                      eqn.t)
    end
  end

  return nothing
end

function ICexp3xplusy{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                            sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol},
                            opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp3xplusy(mesh.coords[:, j, i], eqn.params,
                                      eqn.t)
    end
  end

  return nothing
end

function ICexp2xplus2y{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                            sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol},
                            opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp2xplus2y(mesh.coords[:, j, i], eqn.params,
                                       eqn.t)
    end
  end

  return nothing
end

function ICexp_xy{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                            sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol},
                            opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp_xy(mesh.coords[:, j, i], eqn.params, eqn.t)
    end
  end

  return nothing
end

function ICxplusy{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                            sbp::AbstractSBP{Tsbp}, eqn::AdvectionData{Tsol},
                            opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_xplusy(mesh.coords[:, j, i], eqn.params, eqn.t)
    end
  end

  return nothing
end

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
"ICConstant" => ICConstant,
"ICx5plusy5" => ICx5plusy5,
"ICexp_xplusy" => ICexp_xplusy,
"ICsinwave" => ICsinwave,
"ICsinwavey" => ICsinwavey,
"ICsinwavey_pert" => ICsinwavey_pert,
"ICsinwavexy" => ICsinwavexy,
"ICFile" => ICFile,
"ICmms1" => ICmms1,
"ICx4" => ICx4,
"ICp1" => ICp1,
"ICp2" => ICp2,
"ICp3" => ICp3,
"ICp4" => ICp4,
"ICp5" => ICp5,
"ICexp5xplus4yplus2" => ICexp5xplus4yplus2,
"ICexp5xplusy" => ICexp5xplusy,
"ICexp3xplusy" => ICexp3xplusy,
"ICexp2xplus2y" => ICexp2xplus2y,
"ICexp_xy" => ICexp_xy,
"ICxplusy" => ICxplusy,
)
