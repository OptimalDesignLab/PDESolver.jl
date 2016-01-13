
function set_hash(pkg_name::AbstractString, hash::AbstractString)
# this function checks out a particular git hash of a particular package

  start_dir = pwd()
  cd(Pkg.dir(pkg_name))
  run(`git checkout -b pdesolver_version $hash`)
  cd(start_dir)

end

function install_pkg(pkg_name::AbstractString, git_url::AbstractString, git_commit::AbstractString, pkg_dict, f; force=false)
# pkg_name = name of package
# git_url: value to be passed to Pkg.clone()
# git_commit: commit has to check out
# pkg_dict: output of Pkg.installed, passed as an argument because it takes so long to get
# f: file handle to write progress to
# force: install the package regardless of other flags, default false

  already_installed = haskey(pkg_dict, pkg_name) 
  force_specific = haskey(ENV, "PDESOLVER_FORCE_DEP_INSTALL_$pkg_name") 
  if !already_installed || force_specific || FORCE_INSTALL_ALL || force
    println(f, "Installing package $pkg_name")
    try 
      Pkg.clone(git_url)
      set_hash(pkg_name, git_commit)
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
