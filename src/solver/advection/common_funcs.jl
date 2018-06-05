# common_funcs.jl

@doc """
### AdvectionEquationMod.calc_x5plusy5

  Calculates and returns x^5 + y^5
"""->
function calc_x5plusy5(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh

  x = coords[1]
  y = coords[2]
  u = x^5 + y^5
  
  return u
end # end function calc_x5plusy5

@doc """
### AdvectionEquationMod.calc_exp_xplusy

  Calculates and returns e^(x + y)
"""->
function calc_exp_xplusy(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh

  x = coords[1]
  y = coords[2]
  u = exp(x+y)

  return u
end

@doc """
### AdvectionEquationMod.calc_sinwave

  Calculates and returns sin(-x + t)
"""->
function calc_sinwave(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
 
  x = coords[1]
  y = coords[2]

  return sin(-x + t)
#  return sin( (-x + -y)/2 + t)
end

@doc """
### AdvectionEquationMod.calc_sinwavey

  Calculates and returns sin(y)^2 + 5*sin(y) + 3/sin(y)
"""->
function calc_sinwavey(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  y = coords[2]
  return sin(y)^2 + 5*sin(y) + 3/sin(y)
end

@doc """
### AdvectionEquationMod.calc_sinwavey_pert

  Calculates and returns 1000*sin(x)*calc_sinwavey
"""->
function calc_sinwavey_pert(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  return calc_sinwavey(coords, t)*1000*sin(x)
end

@doc """
### AdvectionEquationMod.calc_sinwave_ampl

  Calculates and returns A*sin(-x + omega*t)
"""->
function calc_sinwave_ampl(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh

  x = coords[1]
  y = coords[2]

  # omega = 1.0
  # A = 2.0
  omega = params.omega
  A = params.sin_amplitude

  return A*sin(-x + omega*t)
end

@doc """
### AdvectionEquationMod.calc_mms1

  Calculates and returns the value of the solution for doing Method of 
  Manufactured solutions.  This is for debugging only, and could change
  at any time.
"""->

function calc_mms1(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  px = 1.5
  py = 1.5
  return sin(px*x)*cos(py*y) + sin(px*x)*cos(py*y)
end

@doc """
### AdvectionEquationMod.calc_mms1dx

  Calculates and returns the x derivative of calc_mms1
"""->
function calc_mms1dx(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  px = 1.5
  py = 1.5
  return px*cos(px*x)*cos(py*y) + px*cos(px*x)*cos(py*y)
end

@doc """
### AdvectionEquationMod.calc_x4

  Calculates and returns a 4th order polynomial in x
"""->
function calc_x4(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  return x.^4 + x.^3 + x.^2 + x + 1
end

@doc """
### AdvectionEquationMod.calc_x5plusy5

  Calculates and returns the x derivative of calc_x4
"""->
function calc_x4der(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  return 4*x.^3 + 3*x.^2 + 2*x + 1
end

"""
  calculates and returns a degree 0 polynomial (ie. a constant)
"""
function calc_p0(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh

  return 1
end

@doc """
### AdvectionEquationMod.calc_p1

  Calculates and returns a 1st order polynomial of x and y
"""->
function calc_p1(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return  x + 1 + y
end

function calc_p1(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return x + y + z + 1
end
  
@doc """
### AdvectionEquationMod.calc_p1dx

  Calculates and returns the x derivative of calc_p1
"""->
function calc_p1dx(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  return 1
end

@doc """
### AdvectionEquationMod.calc_p1dy

  Calculates and returns the y derivative of calc_p1
"""->
function calc_p1dy(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  y = coords[2]
  return 1
end

@doc """
### AdvectionEquationMod.calc_p1dz

  Calculates and returns the z derivative of calc_p1
"""->
function calc_p1dz(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  return 1
end

@doc """
### AdvectionEquationMod.calc_p2

  Calculates and returns a 2nd order polynomial in x and y
"""->
function calc_p2(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return x^2 + x + 1 + y^2 + y
end

function calc_p2(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return x^2 + x + 1 + y^2 + y + z^2 + z
end


@doc """
### AdvectionEquationMod.calc_p2dx

  Calculates and returns a the x derivative of calc_p2
"""->
function calc_p2dx(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  return 2*x + 1
end

@doc """
### AdvectionEquationMod.calc_p2dy

  Calculates and returns the y derivative of calc_p2
"""->
function calc_p2dy(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  y = coords[2]
  return 2*y + 1
end

@doc """
### AdvectionEquationMod.calc_p2dz

  Calculates and returns the z derivative of calc_p2
"""->

function calc_p2dz(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  z = coords[3]
  return 2*z + 1
end


@doc """
### AdvectionEquationMod.calc_p3

  Calculates and returns a 3rd order polynomial in x and y (and z in 3d)
"""->
function calc_p3(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return x^3 + x^2 + x + 1 + y^3 + y^2 + y
end

function calc_p3(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return x^3 + x^2 + x + 1 + y^3 + y^2 + y + z^3 + z^2 + z
end


@doc """
### AdvectionEquationMod.calc_p3dx

  Calculates and returns the x derivataive of calc_p3
"""->
function calc_p3dx(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  return 3*x^2 + 2*x + 1
end

@doc """
### AdvectionEquationMod.calc_p3dy

  Calculates and returns the y derivative of calc_p3
"""->
function calc_p3dy(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  y = coords[2]
  return 3*y^2 + 2*y + 1
end

@doc """
### AdvectionEquationMod.calc_p3dz

  Calculates and returns the z derivative of calc_p3
"""->
function calc_p3dz(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  z = coords[3]
  return 3*z^2 + 2*z + 1
end


@doc """
### AdvectionEquationMod.calc_p4

  Calculates and returns a 4th order polynomial in x and y
"""->
function calc_p4(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return (x^4 + x^3 + x^2 + x + 1) + (y^4 + y^3 + y^2 + y)
end

function calc_p4(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return (x^4 + x^3 + x^2 + x + 1) + (y^4 + y^3 + y^2 + y) + (z^4 + z^3 + z^2 + z)
end

@doc """
### AdvectionEquationMod.calc_p4x

  Calculates and returns the x derivative of calc_p4
"""->

function calc_p4dx(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  return 4*x^3 + 3*x^2 + 2*x + 1
end

@doc """
### AdvectionEquationMod.calc_p4dy

  Calculates and returns the y derivative of calc_p4
"""->
function calc_p4dy(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return 4*y^3 + 3*y^2 + 2*y + 1
end

@doc """
### AdvectionEquationMod.calc_p4dz

  Calculates and returns the z derivative of calc_p4
"""->
function calc_p4dz(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  z = coords[3]
  return 4*z^3 + 3*z^2 + 2*z + 1
end


@doc """
### AdvectionEquationMod.calc_p5

  Calculates and returns a 5th order polynomial in x and y (and z in 3d)
"""->
function calc_p5(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return (x.^5 + x^4 + x^3 + x^2 + x + 1) + (y^5 + y^4 + y^3 + y^2 + y)
end

function calc_p5(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return (x.^5 + x^4 + x^3 + x^2 + x + 1) + (y^5 + y^4 + y^3 + y^2 + y) + (z^5 + z^4 + z^3 + z^2 + z)
end


@doc """
### AdvectionEquationMod.calc_p5dx

  Calculates and returns the x derivative of calc_p5
"""->
function calc_p5dx(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  return 5*x^4 + 4*x^3 + 3*x^2 + 2*x + 1
end

@doc """
### AdvectionEquationMod.calc_p5y

  Calculates and returns the y derivative of calc_p5
"""->

function calc_p5dy(params::ParamTypes, coords::AbstractArray{Tmsh}, t) where Tmsh
  y = coords[2]
  return 5*y^4 + 4*y^3 + 3*y^2 + 2*y + 1
end

@doc """
### AdvectionEquationMod.calc_p5z

  Calculates and returns the z derivative of calc_p5
"""->

function calc_p5dz(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  z = coords[3]
  return 5*z^4 + 4*z^3 + 3*z^2 + 2*z + 1
end


@doc """
### AdvectionEquationMod.calc_exp5xplus4yplus2

Calculates and returns the expression u = exp(5*x + 4*y +2)
"""->

function calc_exp5xplus4yplus2(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return exp(5*x + 4*y +2)
end

@doc """
### AdvectionEquationMod.calc_exp5xplusy

Calculates and return the expression u = exp(5*x + y)
"""->

function calc_exp5xplusy(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return exp(5*x + y)
end

function calc_exp5xplusy(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return exp(5*x + y + z)
end


@doc """
### AdvectionEquationMod.calc_exp3xplusy

Calculates and return the expression u = exp(3*x + y)
"""->

function calc_exp3xplusy(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return exp(3*x + y)
end

function calc_exp3xplusy(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return exp(3*x + y + z)
end


@doc """
### AdvectionEquationMod.calc_exp2xplus2y

Calculates and return the expression u = exp(2*x + 2*y)
"""->

function calc_exp2xplus2y(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return exp(2*x + 2*y)
end

function calc_exp2xplus2y(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return exp(2*x + 2*y + z*z)
end


@doc """
### AdvectionEquationMod.calc_exp_xy

Calculates and returns u = exp(x*y)
"""->

function calc_exp_xy(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return exp(x*y)
end

function calc_exp_xy(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return exp(x*y*z)
end


@doc """
### AdvectionEquationMod.calc_x5plusy5

Calculates and returns u = x+y
"""->

function calc_xplusy(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  return x + y
end

function calc_xplusy(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]
  return x + y + z
end


"""
  u = exp(x + y + z + t) in 3d (z = 0 in 2d)
"""
function calc_unsteadymms(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]

  return exp(x + y + t)
end

function calc_unsteadymms(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]
  z = coords[3]

  return exp(x + y + z + t)
end


function calc_unsteadymmsdx(params::ParamType, coords::AbstractArray{Tmsh}, t) where Tmsh

  calc_unsteadymms(params, coords, t)
end

function calc_unsteadymmsdy(params::ParamType, coords::AbstractArray{Tmsh}, t) where Tmsh

  calc_unsteadymms(params, coords, t)
end

function calc_unsteadymmsdz(params::ParamType3, coords::AbstractArray{Tmsh}, t) where Tmsh

  calc_unsteadymms(params, coords, t)
end


function calc_unsteadymmsdt(params::ParamType, coords::AbstractArray{Tmsh}, t) where Tmsh

  calc_unsteadymms(params, coords, t)
end

function calc_unsteadypoly(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]

  return x + 1 + y + t 
end

function calc_unsteadypolydx(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]

  return 1
end

function calc_unsteadypolydy(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]

  return 1
end

function calc_unsteadypolydt(params::ParamType2, coords::AbstractArray{Tmsh}, t) where Tmsh
  x = coords[1]
  y = coords[2]

  return 1
end


