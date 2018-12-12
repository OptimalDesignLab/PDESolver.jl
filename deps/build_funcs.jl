function checkInput()

  if haskey(ENV, "PDESOLVER_BUNDLE_DEPS")

    if !haskey(ENV, "PDESOLVER_PKGDIR")
      error("PDESOLVER_PKGDIR environment variable must be specified when bundling dependencies")
    end

    pkgdir = ENV["PDESOLVER_PKGDIR"]
    if !isdir(pkgdir)
      error("PDESOLVER_PKGDIR $pkgdir must exist")
    end

    if length(readdir(pkgdir)) != 0
      error("PDESOLVER_PKGDIR $pkgdir must be empty directory")
    end

  end

  if haskey(ENV, "PDESOLVER_UNBUNDLE_DEPS")

    if haskey(ENV, "PDESOLVER_PKGDIR")
      error("Environment variable PDESOLVER_PKGDIR should not be set when unbundling dependencies (PDESOLVER_UNBUNDLE_DEPS).  Copy the contents of the directory directory PDESOLVER_PKGDIR to the Julia package directory (eg. ~/.julia/v0.6) on the new machine")
    end
  end

  return nothing
end


"""
  This function determines if a git identifier is a commit or something else
  (ie. a branch or tag)

  **Inputs**

   * pkg_dir: path to directory containing the repo
   * hash: git object identifier

  **Outputs**

   * bool, true if hash represents a commit, false otherwise
"""
function is_commit(pkg_dir::AbstractString, hash::AbstractString)

  start_dir = pwd()
  if pkg_dir != ""
    cd(pkg_dir)
  end
  
  iscommit = false
  try
    vals = readstring(`git show-ref \| grep $hash`)
    # make sure this is an exact match
    name = split(vals, '/')[end]
    if name[1:end-1] == hash
      iscommit = false
    else
      iscommit = true
    end
  catch x
    iscommit = true
  end

  cd(start_dir)
  return iscommit
end

"""
  Checks out a git object of some kind.  If the object is a tag or a branch,
  it is checkout out as normal.  If it is a commit, a new branch is created
  called detatched_from_commit, where commit is the commit hash

  **Inputs**

   * pkg_dir: directory where the git repo is located
   * hash: git object identifier
"""
function set_hash(pkg_dir::AbstractString, hash::AbstractString)
# this function checks out a particular git hash of a particular package

  new_branchname = "detached_from_$hash"

  start_dir = pwd()
  cd(pkg_dir)
#  branch_name = readstring(`git rev-parse --abbrev-ref HEAD`)
#  if branch_name[1:end-1] != new_branchname
#    run(`git checkout -b pdesolver_version $hash`)
#  end
  iscommit = is_commit("", hash)
  run(`git checkout $hash`)

  if iscommit
    run(`git checkout -b $new_branchname HEAD`)
  end
  cd(start_dir)

  println("finished setting hash")
end

"""
  This function installs a specified package.  Unless overridded, if the package is
  already installed, it will not be modified.

  **Inputs**

  * dir: directory in which to install package
  * pkg_name: name of the package (not including .jl)
  * git_url: URL of git repository
  * git_commit: anything that git checkout can use to identify what to check out
  * pkg_dict: dictionary of installed package
  * f: an IO object to write info to
  * force: force package to install even if already present, (keyword arg), default
           false
"""
function install_pkg(dir::AbstractString, pkg_name::AbstractString, git_url::AbstractString, git_commit::AbstractString, pkg_dict, f; force=false)
# pkg_name = name of package
# git_url: value to be passed to Pkg.clone()
# git_commit: commit has to check out
# pkg_dict: output of Pkg.installed, passed as an argument because it takes so long to get
# f: file handle to write progress to
# force: install the package regardless of other flags, default false

  already_installed = haskey(pkg_dict, pkg_name) 
  pkg_path = joinpath(dir, pkg_name)
  force_specific = haskey(ENV, "PDESOLVER_FORCE_DEP_INSTALL_$pkg_name") 
  if !already_installed || force_specific || FORCE_INSTALL_ALL || force
    println(f, "Installing package $pkg_name")
    try 
      if !isdir(pkg_path)  # if package not already present
        start_dir = pwd()
        cd(dir)
        run(`git clone $git_url`)
        name_ext = string(pkg_name, ".jl")
        run(`mv -v ./$name_ext ./$pkg_name`)
        set_hash(pkg_path, git_commit)
        cd(start_dir)
      end
#      else
#        reset_repo(pkg_name)
#      end
      println(f, "Resolving dependencies")
      Pkg.resolve()
      println(f, "Building package")
      println(f, "  Note: if this step fails, there is no way to report it in this log")
      Pkg.build(pkg_name)
      Libc.flush_cstdio()
      println(f, "  Installation appears to have completed sucessfully")
    catch x
      Libc.flush_cstdio()
      println(f, "Error installing package $pkg_name")
      println(f, "Error is $x")
    end

  else
    println(f, "Skipping installation of package $pkg_name")
    println(f, "already installed: ", already_installed, ", specifically forced: ", force_specific,
               ", force general: ", FORCE_INSTALL_ALL, ", force override: ", force)
  end

  flush(f)
end

function reset_repo(pkg_name::AbstractString)
  start_dir = pwd()
  cd(Pkg.dir(pkg_name))
  run(`git reset --hard`)  # remove changes to files under version control
  run(`git clean -x -f -d`) # remove any files not under version control
end

"""
  This function downloads packages to a specified directory, in preparation for
  copying them to a machine without internet access.

  **Inputs**

   * dir: the directory to put the packages into
   * pkg_name: the name of the package
   * git_url: the git URL of the package
   * git_commit: the name of a commit/tag/branch to checkout
   * f: an IO to write logging messages to
"""
function bundle_pkg(dir::AbstractString, pkg_name::AbstractString, git_url::AbstractString, git_commit::AbstractString, f::IO)
# clone all packages into a specified directory, without installing any of them
  println(f, "Bundling package $pkg_name")

  start_dir = pwd()
  cd(dir)

  println(f, "now in directory ", pwd())
  # check if package already exists
  if isdir(pkg_name)
    println(f, "package $pkg_name is already present, skipping...")
    return
  end

  try
    run(`git clone $git_url`)
    println(f, "ls of pwd() = ", readstring(`ls`))
    name_ext = string(pkg_name, ".jl")
    run(`mv -v ./$name_ext ./$pkg_name`)
    set_hash("./$pkg_name", git_commit)
    println(f, "set hash to ", git_commit)

    if isdir("./$pkg_name/deps")
      cd("./$pkg_name/deps")

      # perform any needed downloads
      if isfile("./download.sh")
        run(`./download.sh`)
      end
    end

    println(f, "Bundling of package $pkg_name appears to have completed successfully")
  catch x
    println(f, "Error bundling package $pkg_name")
    println(f, "error is $x")
  end

  flush(f)
  cd(start_dir)
end

"""
  Unbundle packages and install them, without accessing the internet.
  This is the inverse function to `bundle_pkg()`.

  For this to work, the contents of the directory created by 
  PDESOLVER_BUNDLE_DEPS should be copied to the new machine and placed in the
  Julia package directory (ie. ~/.julia/v0.6 or JULIA_PKGDIR).  When done,
  the contents of this directory should look like

  ~/.julia/v0.6:
    - PDESolver
    - ODLCommonTools
    - PumiInterface
    - etc...


  Note that the order in which the dependencies is installed is important.
  They must be installed in an order such that all the dependencies of a
  package are installed before the package itself

  The order in which packages are installed is determined by the order in
  which they appear listed in `build.jl`.

  **Inputs**

   * dir: the julia package directory
   * pkg_name: the name of the package
   * f: an IO object to write logging messages to
"""
function unbundle_pkg(dir::AbstractString, pkg_name::AbstractString, f::IO)
# install a package by directly calling the build.jl script
# this bypasses the REQUIRE file, and therefore any internet access
# the dependencies must be installed in order for this to work

  println(f, "unbundling package $pkg_name")
  start_dir = pwd()
  pkgdir = joinpath(dir, pkg_name)
  cd(pkgdir)
  if isdir("./deps")
    cd("./deps")

    if isfile("./build.jl")
      try
        buildpath = joinpath(dir, pkg_name, "deps", "build.jl")
        include(buildpath)
        println(f, "unbundling package $pkg_name appears to have completed successfully")
      catch x
        println(f, "Error unbundling package $pkg_name")
        println(f, "Error is $x")
      end

    else
      println(f, "No build.jl file for this package, skipping...")

    end   # end if isfile()

  else  
    println(f, "No deps directory for this package, skipping...")
  end  # end isdir()

  flush(f)
end

