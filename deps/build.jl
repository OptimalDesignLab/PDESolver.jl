# download and build all non metadata dependences

sbp_git = "https://github.com/OptimalDesignLab/SummationByParts.jl.git"
pdesolvercommon_git = "https://github.com/OptimalDesignLab/ODLCommonTools.jl.git"
pumiinterface_git = "https://github.com/OptimalDesignLab/PumiInterface.jl.git"
petsc_git = "https://github.com/JaredCrean2/PETSc.jl.git"
pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "ODLCommonTools")
  Pkg.clone(pdesolvercommon_git)
  Pkg.build("ODLCommonTools")
end

if !haskey(pkg_dict, "SummationByParts")
  Pkg.clone(sbp_git)
  Pkg.build("SummationByParts")
end

if !haskey(pkg_dict, "PumiInterface")
  Pkg.clone(pumiinterface_git)
  Pkg.build("PumiInterface")
end

if !haskey(pkg_dict, "PETSc")
  Pkg.clone(petsc_git)
  petsc_dir = Pkg.dir("PETSc")
  start_dir = pwd()
  cd(petsc_dir)
  run(`git checkout 58999c0`)
  cd(start_dir)
end


