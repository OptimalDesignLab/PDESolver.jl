
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

  if !haskey(pkg_dict, pkg_name)  || haskey(ENV, "PDESOLVER_FORCE_DEP_INSTALL_$pkg_name")  || FORCE_INSTALL_ALL || force
    println(f, "Installing package $pkg_name")
    try 
      Pkg.clone(git_url)
      set_has(git_url)
      Pkg.build(pkg_name)
      println(f, "  Installation appears to have completed sucessfully")
    catch
      println(f, "Error installing package $pkg_name")
    end

  else
    println(f, "Skipping installation of package $pkg_name")
  end

  flush(f)
end
