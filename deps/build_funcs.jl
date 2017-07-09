
function set_hash(pkg_dir::AbstractString, hash::AbstractString)
# this function checks out a particular git hash of a particular package

  new_branchname = "pdesolver_version"

  start_dir = pwd()
  cd(pkg_dir)
#  branch_name = readall(`git rev-parse --abbrev-ref HEAD`)
#  if branch_name[1:end-1] != new_branchname
#    run(`git checkout -b pdesolver_version $hash`)
#  end
  run(`git checkout $hash`)
  run(`git checkout -b pdesolver_version HEAD`)
  cd(start_dir)

end

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
       Pkg.resolve()
      Pkg.build(pkg_name)
      println(f, "  Installation appears to have completed sucessfully")
    catch x
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


function bundle_pkg(dir::AbstractString, pkg_name::AbstractString, git_url::AbstractString, git_commit::AbstractString, f)
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
    println(f, "ls of pwd() = ", readall(`ls`))
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

function unbundle_pkg(dir::AbstractString, pkg_name::AbstractString, f)
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

