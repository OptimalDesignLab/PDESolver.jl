# functions that populate the initial conditions
# List of functions:

"""
  This macro facilitates make initial condition functions.  Specifically,
  IC functions that use a function in common_funcs.jl to populate the initial
  condition can be created with this macro

  **Inputs**

   * namestub: the name stub.  If the function from common_funcs.jl is named
               `calcFoo`, then `namestub` should be `Foo` and the new initial
               condition function will be `ICFoo`.
   * docstring: docstring for the newly created function `ICFoo`.  This argument
                is optional, but good programming practices dictate that it
                should be provided in all cases.

  **Example Usage**

   ```
     @makeIC Foo \"\"\"
Docstring for ICFoo
\"\"\"
  ```
"""
macro makeIC(namestub, docstring="""""")

  # names of IC function, common_funcs.jl function
  fname = Symbol(string("IC", namestub))
  calcname = Symbol(string("calc", namestub))


  # introducing an empty docstring is not the same as not supplying a docstring
  # figure out which to do
  if docstring != ""
    docex = quote
              @doc $docstring $fname
            end
  else
    docex = :()
  end

  return esc(
    quote


      function $fname(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol}) where {Tmsh, Tsol,}

        q = zeros(Tsol, mesh.numDofPerNode)
        for i=1:mesh.numEl
          for j=1:mesh.numNodesPerElement
            dofs = sview(mesh.dofs, :, j, i)
            coords = sview(mesh.coords, :, j, i)
            $calcname(eqn.params, coords, q)
            for k=1:mesh.numDofPerNode
              u0[dofs[k]] = q[k]
            end
          end
        end

        return nothing
      end

    end  # end quote
  )

end

@doc """
### EulerEquationMod.ICZero

  Sets all components of the solution to zero

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICZero(mesh::AbstractMesh{Tmsh}, 
operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts, 
u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = operator.numnodes
dofpernode = mesh.numDofPerNode
for i=1:numEl

  for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = 0.0
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
  end
end

return nothing

end  # end function


@doc """
### EulerEquationMod.ICOnes

  Sets all components of the solution to 1.0

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->

function ICOnes(mesh::AbstractMesh{Tmsh}, 
operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts,
u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = 1.0
      u0[dofnum_rhou] = 1.0
      u0[dofnum_rhov] = 1.0
      u0[dofnum_e] = 1.0

      # u0 = 2*u0
    end
  end

  return nothing
end # end function ICOnes


@doc """
### EulerEquationMod.ICRho1E2

  Sets all density values to 1.0 and energy values to 2.0

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->


function ICRho1E2(mesh::AbstractMesh{Tmsh}, 
operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts, 
u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = operator.numnodes
dofpernode = mesh.numDofPerNode
for i=1:numEl
  for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = 1.0
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 2.0
  end
end

return nothing

end  # end function

@makeIC Rho1E2U1VW0 """
### EulerEquationMod.ICRho1E2U1VW0

  Sets the density values 1.0, x momentum to 1.0, 
  v & w momenta to 0.0, and energy to 2.0 at a node.

  It should work for 2D and 3D meshes.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""

@makeIC Rho1E2U3 """
### EulerEquationMod.ICRho1E2U3

  Sets all components density values to 1.0, x and y momenta to 0.35355, and
  energy to 2.0

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""

@makeIC FreeStream """
### EulerEquationMod.ICFreeStream

  Sets all components of the solution to the free stream condition according
  to the angle of attack and and Mach number.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""

@makeIC FreeStream0 """
  Like [`calcFreeStream`](@ref), but uses [`calcFreeStream0`](@ref) instead of
  [`calcFreeStream`](@ref).
"""




@makeIC Vortex """
 what is this? how is it different than ICIsentropic Vortex?
"""

@doc """
### EulerEquationMod.ICsmoothHeavisideder

  Sets the density to the derivative of the smooth Heaviside function, all 
  other components to zero.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICsmoothHeavisideder(mesh::AbstractMesh{Tmsh}, 
            operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol},
            opts, u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}
# calculate the value of the smooth heaviside function derivative at a location x
# x0 is specified within this function

# smooth heaviside  parameters
  x0 = 0
  L = 5
  k = 5

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  for i=1:numEl
  #  dofnums_i = sview(mesh, i)  # get dof nums for this element
  #  coords = sview(mesh, [i])

    for j=1:nnodes
        coords_j = sview(mesh.coords, :, j, i)
        dofnums_j = sview(mesh.dofs, :, j, i)
   
        # get dof numbers for each variable
        dofnum_rho = dofnums_j[1]
        dofnum_rhou = dofnums_j[2]
        dofnum_rhov = dofnums_j[3]
        dofnum_e = dofnums_j[4]

        x = coords_j[1]
        y = coords_j[2]

        # apply initial conditions here
        u0[dofnum_rho] = L*(2*k*e^(-2*k*x))/(e^(-2*k*x) +1 )^2
        u0[dofnum_rhou] = 0.0
        u0[dofnum_rhov] = 0.0
        u0[dofnum_e] = 0.0
    end
  end

  return nothing
end


@doc """
### EulerEquationMod.ICsmoothHeaviside

  Sets the density to the smooth Heaviside function, all other components to
  zero.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICsmoothHeaviside(mesh::AbstractMesh{Tmsh}, 
                           operator::AbstractSBP{Tsbp}, 
                           eqn::EulerData{Tsol}, 
                           opts, u0::AbstractArray{Tsol, 1}) where {Tmsh, Tsbp, Tsol}
  # calculate the value of the smooth heaviside function at a location x
  # x0 is specified within this function

  # smooth heaviside  parameters
  x0 = 0
  L = 5
  k = 5



  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  for i=1:numEl
    for j=1:nnodes
      coords = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = L/(1 + e^(-k*(x-x0)))
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
    end
  end

  return nothing
end

@makeIC IsentropicVortex """

### EulerEquationMod.ICIsentropicVortex

Sets the solution to the steady isentropic vortex solution.

Inputs:
mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""

@doc """
### EulerEquationMod.ICIsentropicVortexWithNoise

  Sets the solution to the steady isentropic vortex solution plus 
  a small random noise component.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICIsentropicVortexWithNoise(mesh::AbstractMesh{Tmsh},
                   operator::AbstractSBP{Tsbp}, 
                   eqn::EulerData{Tsol}, 
                   opts, u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}
  # populate u0 with initial values
  # this is a template for all other initial conditions

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, 4)
  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
      calcIsentropicVortex(eqn.params, coords_j, sol)

      # apply initial conditions here
      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k] + 0.1*rand()
      end

    end
  end

  return nothing

end  # end function



@makeIC UnsteadyVortex """
### EulerEquationMod.ICUnsteadyVortex

  Sets the solution to the unsteady vortex problem.  eqn.params.t is used to
  determine what time to use for the solution.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""
@makeIC UnsteadyVortex2 """
  Vortex travelling at an angle
"""

@makeIC InvChannelIC """
  Initial condition for SU2 bump in inviscid channel case.  See also
  The subsonic inflow and subsonic outflow boundary conditions.
"""
function ICInvChannel(mesh::AbstractMesh{Tmsh}, 
        operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, 
        opts, u0::AbstractArray{Tsol}) where {Tmsh, Tsbp, Tsol}

  sol = zeros(Tsol, mesh.numDofPerNode)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofnums_j = sview(mesh.dofs, :, j, i)
      coords_j = sview(mesh.coords, :, j, i)
      calcInvChannelIC(eqn.params, coords_j, sol)

      for k=1:mesh.numDofPerNode
        u0[dofnums_j[k]] = sol[k]
      end
    end
  end

  return nothing
end

@doc """
### EulerEquationMod.ICFile

  This function reads a vector from a file on disk and set the solution to it.
  The vector must contain the same number of entries as there are degrees of 
  freedom in the mesh. 

  This function is useful for things like restarting from a checkpoint.
  In this case, the file should be the output of writedlm(eqn.q).  The degree 
  of freedom number must be the same for both simulation for this to work (the 
  file contains no degree of freedom number information).


  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICFile(mesh::AbstractMesh{Tmsh}, 
operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts, 
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

@makeIC Exp """
  Assigns exp(k*x*y*z) as the initial condition, of each node, where k is 
  the index of the degree of freedom of the node
"""

@makeIC PeriodicMMS """
  Writes calcPeriodicMMS to the initial condition vector u0
"""

"""
  This function applies the initial condition for the Taylor Green vortex,
  using the constants in Gassner, Winters, and Kopriva's Split form Nodal
  DG paper
"""
function ICTaylorGreen(mesh::AbstractMesh{Tmsh}, sbp, 
          eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol}) where {Tmsh, Tsol,}

  # parameters
  M = 1.0  # Mach number
  gamma_1 = eqn.params.gamma_1
  gamma = eqn.params.gamma

  q = zeros(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofs = sview(mesh.dofs, :, j, i)
      coords = sview(mesh.coords, :, j, i)
      x = coords[1]
      y = coords[2]
      z = coords[3]

      p = 100/gamma + (1/16)*( cos(2*x)*cos(2*z) + 2*cos(2*y) + 2*cos(2*x) + cos(2*y)*cos(2*z) )
      q[1] = 1
      q[2] = q[1]*M*sin(x)*cos(y)*cos(z)
      q[3] = -q[1]*M*cos(x)*sin(y)*cos(z)
      q[4] = 0
      q[5] = p/gamma_1 + (q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/(2*q[1])

      for k=1:mesh.numDofPerNode
        u0[dofs[k]] = q[k]
      end
    end
  end

  return nothing
end

@makeIC ChannelMMS """
  Initial condition of channel MMS
"""

@makeIC Square1D """
  Initial for square wave in 1D
"""

@makeIC Square2D """
  Initial for square wave in 2D
"""

@makeIC SedovExplosion """
  Initial for Sedov Explosion
"""
 

# declare a const dictionary here that maps strings to function (used for input arguments)
"""
  Map IC names to functions.  Generally the name is the same as the function
  name

  Every function in this dictionary must have the signature:

  ```
    Foo(mesh::AbstractMesh{Tmsh}, sbp, 
          eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol}) where {Tmsh, Tsol,}
  ```

  where `u0` is the vector to be overwritten with the initial condition, 
  `mesh` is the mesh object`, `sbp` is the SBP operator, `eqn` is the equation
  object, and `opts` is the optios dictionary.

  Many of the initial condition functions call a function in common_funcs.jl
  to compute the state.  See the macro [`makeIC`](@ref) for an easy way
  to make such functions.  If these functions depend on time, they use
  `eqn.params.t` to get the current time value.  This means IC functions can
  be used to get analytical solutions at times other than `t = 0`.

"""
global const ICDict = Dict{Any, Function}(
"ICZero" => ICZero,
"ICOnes" => ICOnes,
"ICRho1E2" => ICRho1E2,
"ICRho1E2U1VW0" => ICRho1E2U1VW0,
"ICRho1E2U3" => ICRho1E2U3,
"ICFreeStream" => ICFreeStream,
"ICFreeStream0" => ICFreeStream0,
"ICVortex" => ICVortex,
#"ICLinear" => ICLinear,
"ICsmoothHeavisideder" => ICsmoothHeavisideder,
"ICsmoothHeaviside" => ICsmoothHeaviside,
"ICIsentropicVortex" => ICIsentropicVortex,
"ICUnsteadyVortex" => ICUnsteadyVortex,
"ICUnsteadyVortex2" => ICUnsteadyVortex2,
"ICIsentropicVortexWithNoise" => ICIsentropicVortexWithNoise,
"ICInvChannel" => ICInvChannelIC,
"ICFile" => ICFile,
"ICExp" => ICExp,
"ICPeriodicMMS" => ICPeriodicMMS,
"ICTaylorGreen" => ICTaylorGreen,
"ICChannelMMS" => ICChannelMMS,
"ICSquare1D" => ICSquare1D,
"ICSquare2D" => ICSquare2D,
"ICSedovExplosion" => ICSedovExplosion,
)


