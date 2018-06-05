# download and build all non metadata dependences

# install PkgFix if not present
if !isdir(joinpath(Pkg.dir(), "PkgFix"))
  Pkg.clone("https://github.com/OptimalDesignLab/PkgFix.jl.git")
end
Pkg.checkout("PkgFix", "upgrade_0.6")

using PkgFix



include("build_funcs.jl")

"""
  This function install all of PDESolvers dependencies (currently PDESolver
  itself does not need to be built).

  Environmental variables:
    PDESOLVER_INSTALL_DEPS_MANUAL: install all dependencies manually (ie. using 
                                   the particular commit stored in this file.
                                   Note that if the package is already
                                   installed, this does not do anything

    PDESOLVER_FORCE_DEP_INSTALL_ALL: download and install all dependencies,
          checking out the version stored in this file, even if the package
          is already installed

    PDESOLVER_FORCE_DEP_INSTALL_pkg_name: where pkg_name is replaced with the
          name of a dependency, forces the installation of only that dependency.

    PDESOLVER_BUNDLE_DEPS: download packages to a specified directory, do not
                           install

    PDESOLVER_PKGDIR: the specified directory for PDESOLVER_BUNDLE_DEPS

    PDESOLVER_UNBUNDLE_DEPS: install PDESolver, using the packages in
                             PDESOLVER_PKGDIR rather than downloading from the
                             internet

    TODO: add a bundling mode that copies the currently installed packages to
          a directory, rather than downloading them from the internet
"""
function installPDESolver()

  checkInput()

  f = open("install.log", "w")
  #------------------------------------------------------------------------------
  # these packages are always installed manually
  #-----------------------------------------------------------------------------
  # [pkg_name, git url, commit identified]
  std_pkgs = [
              "PkgFix" "https://github.com/OptimalDesignLab/PkgFix.jl.git" "upgrade_0.6";
              "ArrayViews"   "https://github.com/JaredCrean2/ArrayViews.jl.git" "work"
              "ODLCommonTools" "https://github.com/OptimalDesignLab/ODLCommonTools.jl.git" "v0.4";
              "SummationByParts" "https://github.com/OptimalDesignLab/SummationByParts.jl.git" "jc_v0.3";
              "PumiInterface" "https://github.com/OptimalDesignLab/PumiInterface.jl.git" "v0.8";
              "PETSc2" "https://github.com/OptimalDesignLab/PETSc2.jl.git" "v0.2"
              ]


  pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

  # force installation to a specific commit hash even if package is already present
  # or would normally be install via the REQUIRE file

  global const FORCE_INSTALL_ALL = haskey(ENV, "PDESOLVER_FORCE_DEP_INSTALL_ALL")

  # figure out the package directory
  if haskey(ENV, "PDESOLVER_BUNDLE_DEPS")
    pkgdir = ENV["PDESOLVER_PKGDIR"]
  else  # unbundling deps
    pkgdir = Pkg.dir()
  end

  println(f, "PDESOLVER_BUNDLE_DEPS = ", haskey(ENV, "PDESOLVER_BUNDLE_DEPS"))
  println(f, "PDESOLVER_PKGDIR = ", haskey(ENV, "PDESOLVER_PKGDIR"))

  println(f, "\n---Installing non METADATA packages---\n")
  for i=1:size(std_pkgs, 1)
    pkg_name = std_pkgs[i, 1]
    git_url = std_pkgs[i, 2]
    git_commit = std_pkgs[i, 3]
    if haskey(ENV, "PDESOLVER_BUNDLE_DEPS")
      bundle_pkg(pkgdir, pkg_name, git_url, git_commit, f)
    elseif haskey(ENV, "PDESOLVER_UNBUNDLE_DEPS")
      unbundle_pkg(pkgdir, pkg_name, f)
    else  # do regular install
      install_pkg(pkgdir, pkg_name, git_url, git_commit, pkg_dict, f)
    end
  end

  println(f, "\n---Finished installing non METADATA packages---\n")
  #------------------------------------------------------------------------------
  # these packages are not always installed manually
  # this section provides a backup mechanism if the regular version resolution
  # mechanism fails
  #------------------------------------------------------------------------------

  # array of [pkg_name, git_url,  commit_hash/tag/branch]
  # these must be in dependency order for unbundling to work
  pkg_list = [
              "Compat" "https://github.com/JuliaLang/Compat.jl.git" "v0.66.0";
              "SHA" "https://github.com/staticfloat/SHA.jl.git" "v0.5.7";
              "URIParser" "https://github.com/JuliaWeb/URIParser.jl.git" "v0.3.1";
              "BinDeps" "https://github.com/JuliaLang/BinDeps.jl.git" "v0.8.8";
              "FactCheck" "https://github.com/JuliaLang/FactCheck.jl.git" "v0.4.3";
              "MPI" "https://github.com/JuliaParallel/MPI.jl.git" "v0.6.0";
              ]
              
  
    println(f, "\n---Considering manual package installations---\n")
    for i=1:size(pkg_list, 1)

      pkg_name = pkg_list[i, 1]
      git_url = pkg_list[i, 2]
      git_commit = pkg_list[i, 3]

      if haskey(ENV, "PDESOLVER_INSTALL_DEPS_MANUAL") || haskey(ENV, "PDESOLVER_FORCE_DEP_INSTALL_$pkg_name") || FORCE_INSTALL_ALL

        install_pkg(pkgdir, pkg_name, git_url, git_commit, pkg_dict, f, force=true)
      end
      if haskey(ENV, "PDESOLVER_BUNDLE_DEPS")
        bundle_pkg(pkgdir, pkg_name, git_url, git_commit, f)
      elseif haskey(ENV, "PDESOLVER_UNBUNDLE_DEPS")
        unbundle_pkg(pkgdir, pkg_name, f)
      end

    end

    println(f, "\n---Finished manual package installations---\n")


  close(f)

end  # end function


# run the installation
installPDESolver()

# generate the known_keys dictonary
if !haskey(ENV, "PDESOLVER_BUNDLE_DEPS")
  start_dir = pwd()
  input_path = joinpath(Pkg.dir("PDESolver"), "src/input")
  cd(input_path)
  include(joinpath(input_path, "extract_keys.jl"))
  cd(start_dir)
end


