# download and build all non metadata dependences
include("build_funcs.jl")

function installPDESolver()

  f = open("install.log", "a+")
  #------------------------------------------------------------------------------
  # these packages are always installed manually
  #-----------------------------------------------------------------------------
  # [pkg_name, git url, commit identified]
  std_pkgs = [ "ArrayViews" "https://github.com/JaredCrean2/ArrayViews.jl.git" "204cc23138cacbd5fc19e861c8038937541952b6";
              "ODLCommonTools" "https://github.com/OptimalDesignLab/ODLCommonTools.jl.git" "HEAD";
              "SummationByParts" "https://github.com/OptimalDesignLab/SummationByParts.jl.git" "bbbfe37e14a89432dd311c803391de0717d312e3";
              "PumiInterface" "https://github.com/OptimalDesignLab/PumiInterface.jl.git" "eec844306c99ab9dbf10bee673020b479e24dc65";
              "MPI" "https://github.com/JaredCrean2/MPI.jl.git" "d68859248ddeed79abe9005a5c75caa35f10885e"
              ]


  # manually install MPI until the package matures
  # handle Petsc specially until we are using the new version
  petsc_git = "https://github.com/JaredCrean2/Petsc.git"

  pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

  # force installation to a specific commit hash even if package is already present
  # or would normally be install via the REQUIRE file

  global const FORCE_INSTALL_ALL = haskey(ENV, "PDESOLVER_FORCE_DEP_INSTALL_ALL")

  # figure out the package directory
  if haskey(ENV, "PDESOLVER_PKGDIR")
    pkgdir = ENV["PDESOLVER_PKGDIR"]
  else
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


  # now install PETSc
  pkg_name = "PETSc"
  git_url = petsc_git
  already_installed = haskey(pkg_dict, pkg_name) 
  force_specific = haskey(ENV, "PDESOLVER_FORCE_DEP_INSTALL_$pkg_name") 

  if !haskey(ENV, "PDESOLVER_BUNDLE_DEPS") && !already_installed  || force_specific  || FORCE_INSTALL_ALL
    println(f, "Installing package $pkg_name")

    try 
      start_dir = pwd()
      cd(pkgdir)
      run(`git clone $petsc_git`)
      run(`mv -v ./Petsc ./PETSc`)
      Pkg.build("PETSc")
      cd(start_dir)
      println(f, "  Installation appears to have completed sucessfully")
    catch x
      println(f, "Error installing package $pkg_name")
      println(f, "Error is $x")
    end

  elseif haskey(ENV, "PDESOLVER_BUNDLE_DEPS")
    println(f, "Bundling package PETSc")
    start_dir = pwd()
    try
      cd(pkgdir )
      println(f, "current directory is now ", pwd() )
      run(`git clone $petsc_git`)
      run(`mv -v ./Petsc ./PETSc`)
      
      cd("./PETSc/deps")
      println(f, "pwd = ", pwd())
      run(`./download.sh`)
      println(f, "Bundling of package PETSc appears to have completed sucessfully")
    catch x
      print(f, "Bunding of backapge PETSc failed")
      println(f, "Error is $x")
    end
    cd(start_dir)
  elseif haskey(ENV, "PDESOLVER_UNBUNDLE_DEPS")
    unbundle_pkg(pkgdir, pkg_name, f)
  else
    println(f, "bundle deps = ", haskey(ENV, "PDESOLVER_BUNDLE_DEPS"))
    println(f, "Skipping installation of package $pkg_name")
    println(f, "  already installed: ", already_installed, ", specifically forced: ", force_specific, ", force general: ", FORCE_INSTALL_ALL)
  end


  println(f, "\n---Finished installing non METADATA packages---\n")
  #------------------------------------------------------------------------------
  # these packages are not always installed manually
  # this section provides a backup mechanism if the regular version resolution
  # mechanism fails
  #------------------------------------------------------------------------------

  # array of [pkg_name, commit_hash]
  pkg_list = ["Compat" "https://github.com/JuliaLang/Compat.jl.git" "aa23cb98dbdba0eab82d465f87b1a936685548c0";
              "URIParser" "https://github.com/JuliaWeb/URIParser.jl.git" "1c4c5f2af17e57617c018ad060f0ec3c9dc5946b";
              "FactCheck" "https://github.com/JuliaLang/FactCheck.jl.git" "e3739d5fdf0e54bc1e74957c060c693cd8ce9cd6";
#              "ArrayViews" "https://github.com/JuliaLang/ArrayViews.jl.git" "4ec55697fc4f9cba522a5137b96d502230269910";
              "SHA" "https://github.com/staticfloat/SHA.jl.git" "90144b2c9e6dd41582901ca0b311215b6bfb3f10";
              "BinDeps" "https://github.com/JuliaLang/BinDeps.jl.git" "ce03a36a969eedc5641aff1c6d7f8f886a17cc98";
              "NaNMath" "https://github.com/mlubin/NaNMath.jl.git" "969151c5ff8022487379279ebac4a239a400dd44";
              "Calculus" "https://github.com/johnmyleswhite/Calculus.jl.git" "cb42f3699177449a42bdc3461c8aea8777aa8c39";
              "DualNumbers" "https://github.com/JuliaDiff/DualNumbers.jl.git" "34ae2b236f4853028fc60a5efed42bd17a33230f";
              "Debug" "https://github.com/toivoh/Debug.jl.git" "0e733093cd71c67bd40ac1295e54153c3db0c751";
#              "MPI" "https://github.com/JuliaParallel/MPI.jl.git" "c546ee896f314340dc61e8bf7ab71f979c57d73c";
              "ForwardDiff" "https://github.com/JuliaDiff/ForwardDiff.jl.git" "d6714170e667027e9e53aa5daf941c3ef5252e7b"]
  
#=
  pkg_list = ["Compat"  "5b03745a6a948781329c444f08ad67cff63f91f7";
  "URIParser" "1c4c5f2af17e57617c018ad060f0ec3c9dc5946b";
  "FactCheck" "e3caf3143b13d4ce9bc576579d5ed9cd4aaefb11";
  "ArrayViews" "4ec55697fc4f9cba522a5137b96d502230269910"
  "SHA" "6f249c3e3f8ec0503f11e95ac5dca8edf210ab2c";
  "BinDeps" "0ff7ad492955bca87155f62511ac37d662616de3";
  "NaNMath" "adb02521f6fd7dfa7c64bbb38611d1f017e026f1";
  "Calculus" "cb42f3699177449a42bdc3461c8aea8777aa8c39";
  "DualNumbers" "7d4eedfa1bd6c0af7d844d085d0c2ce6b2357aa7";
  "Debug" "0e733093cd71c67bd40ac1295e54153c3db0c751";
  "MPI" "c546ee896f314340dc61e8bf7ab71f979c57d73c";
  "ForwardDiff" "d6714170e667027e9e53aa5daf941c3ef5252e7b"]
=#

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

  # generate the known_keys dictonary
  if !haskey(ENV, "PDESOLVER_BUNDLE_DEPS")
    start_dir = pwd()
    input_path = joinpath(Pkg.dir("PDESolver"), "src/input")
    cd(input_path)
    include(joinpath(input_path, "extract_keys.jl"))
    cd(start_dir)
  end

end  # end function


# run the installation
installPDESolver()
