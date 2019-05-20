# Solves subsonic airfoil
#   ARGS[1] = path to input file

# load PDESolver and solve the PDE
using PDESolver
using EulerEquationMod


const dirname = "./work"
# create temporary directory to do the solve in
if MPI.Comm_rank(MPI.COMM_WORLD) == 1
  if isdir(dirname)
    error("$dirname already exists, please delete it before proceeding")
  end
  mkdir(dirname)
end
MPI.Barrier(MPI.COMM_WORLD)
cd(dirname)

solvePDE(joinpath("../", ARGS[1]))
