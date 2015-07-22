# download and build all non metadata dependences

sbp_git = "https://github.com/OptimalDesignLab/SummationByParts.jl.git"
pdesolvercommon_git = "https://github.com/OptimalDesignLab/PDESolverCommon.jl.git"
pumiinterface_git = "https://github.com/OptimalDesignLab/PumiInterface.jl.git"
pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "PDESolverCommon")
  Pkg.clone(pdesolvercommon_git)
  Pkg.build("PDESolverCommon")
end

if !haskey(pkg_dict, "SummationByParts")
  Pkg.clone(sbp_git)
  Pkg.build("SummationByParts")
end

if !haskey(pkg_dict, "PumiInterface")
  Pkg.clone(pumiinterface_git)
  Pkg.build("PumiInterface")
end


