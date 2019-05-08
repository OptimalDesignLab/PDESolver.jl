# Solves subsonic airfoil
#   ARGS[1] = path to input file


# create temporary directory to do the solve in
const dirname = "./work"
if isdir(dirname)
  error("$dirname already exists, please delete it before proceeding")
end

mkdir(dirname)
cd(dirname)

# load PDESolver and solve the PDE
using PDESolver
using EulerEquationMod

solvePDE(joinpath("../", ARGS[1]))
