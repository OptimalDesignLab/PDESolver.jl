# ic.jl
# Needed to initialize a problem.
# TODO: doc these
function ICConstant(mesh::AbstractMesh, sbp::AbstractOperator, eqn::AdvectionData, opts,
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

function ICx5plusy5(mesh::AbstractMesh{Tmsh}, 
sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol, Tres, 2}, 
opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}

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

function ICx5plusy5(mesh::AbstractMesh{Tmsh}, 
sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol, Tres, 3}, 
opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}

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

function ICexp_xplusy(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol, Tres, 2}, 
                    opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}

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

function ICexp_xplusy(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol, Tres, 3}, 
                    opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol, Tres}

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



function ICsinwave(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_sinwave(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy


function ICsinwavey(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_sinwavey(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy


function ICsinwavey_pert(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_sinwavey_pert(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy


function ICsinwavexy(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  n = 1  # number of sin waves
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = mesh.dofs[1, j, i]
      coords = sview(mesh.coords, :, j, i)
      u0[dofnums_j] = sin(2*pi*n*coords[1]) + sin(2*pi*n*coords[2])
      if mesh.dim == 3
        u0[dofnums_j] += sin(2*pi*n*coords[3])
      end
    end
  end

  return nothing
end

function ICsinwave_ampl(mesh::AbstractMesh{Tmsh},
                    sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
                    opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_sinwave_ampl(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function ICsinwave_ampl

function ICmms1(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_mms1(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy

function ICx4(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_x4(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy


function ICp0(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
  	  u0[dofnums_j] = calc_p0(eqn.params, mesh.coords[:, j, i], eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy


function ICp1(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_p1(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy

function ICp2(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_p2(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy

function ICp3(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_p3(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy

function ICp4(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_p4(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy


function ICp5(mesh::AbstractMesh{Tmsh}, 
  sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
  opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_p5(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end # end function exp_xplusy

function ICexp5xplus4yplus2(mesh::AbstractMesh{Tmsh}, 
          sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
          opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp5xplus4yplus2(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end

function ICexp5xplusy(mesh::AbstractMesh{Tmsh}, 
          sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
          opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp5xplusy(eqn.params, mesh.coords[:, j, i],
                                      eqn.t)
    end
  end

  return nothing
end

function ICexp3xplusy(mesh::AbstractMesh{Tmsh}, 
          sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
          opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp3xplusy(eqn.params, mesh.coords[:, j, i],
                                      eqn.t)
    end
  end

  return nothing
end

function ICexp2xplus2y(mesh::AbstractMesh{Tmsh}, 
          sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
          opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp2xplus2y(eqn.params, mesh.coords[:, j, i],
                                       eqn.t)
    end
  end

  return nothing
end

function ICexp_xy(mesh::AbstractMesh{Tmsh}, 
          sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
          opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_exp_xy(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end

function ICxplusy(mesh::AbstractMesh{Tmsh}, 
          sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
          opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_xplusy(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end


function ICunsteadymms(mesh::AbstractMesh{Tmsh}, 
          sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
          opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_unsteadymms(eqn.params, mesh.coords[:, j, i], eqn.t)
    end
  end

  return nothing
end

function ICunsteadypoly(mesh::AbstractMesh{Tmsh}, 
          sbp::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol},
          opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      u0[dofnums_j] = calc_unsteadypoly(eqn.params, mesh.coords[:, j, i], eqn.t)
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
function ICFile(mesh::AbstractMesh{Tmsh}, 
operator::AbstractOperator{Tsbp}, eqn::AdvectionData{Tsol}, opts, 
u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}
# populate u0 with initial values from a disk file
# the file name comes from opts["ICfname"]

  fname = get_parallel_fname(opts["ICfname"], mesh.myrank)
  vals = readdlm(fname)

  @assert length(vals) == mesh.numDof

  for i=1:mesh.numDof
    u0[i] = vals[i]
  end

end

"""
  Dictionary that maps IC names to functions.  Every new IC should be added
  to the list
"""
global const ICDict = Dict{Any, Function}(
"ICConstant" => ICConstant,
"ICx5plusy5" => ICx5plusy5,
"ICexp_xplusy" => ICexp_xplusy,
"ICsinwave" => ICsinwave,
"ICsinwavey" => ICsinwavey,
"ICsinwavey_pert" => ICsinwavey_pert,
"ICsinwavexy" => ICsinwavexy,
"ICsinwave_ampl" => ICsinwave_ampl,
"ICmms1" => ICmms1,
"ICx4" => ICx4,
"ICp0" => ICp0,
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
"ICunsteadymms" => ICunsteadymms,
"ICunsteadypoly" => ICunsteadypoly,
"ICFile" => ICFile,
)
