find . -name '*.dat' -exec rm '{}' \;
#find . -name 'arg_dict_output.jl' -print
#find . -name 'newmesh_linear' -print
#find . -name 'mesh_complete' -print
#find . -name 'output_init'
#find . -name 

dirnames=("arg_dict_output.jl" "newmesh_linear" "mesh_complete" "output_init" "solution_done" "solution_error" "solution_ic" "solution_relfunc")

for dirname in "${dirnames[@]}";
do
  find . -name $dirname -type d -exec rm -r '{}' \;
done

#rm -v ./*vtu
#rm -v ./convergence.dat
#rm -v ./coords_output.dat
#rm -v ./edge_vertnums.dat
#rm -v ./face_vertnums.dat
#rm -v ./IC.dat
#rm -v ./jacobian*.dat
#rm -v ./newton_solution.dat
#rm -v ./solution.dat
#rm -v ./sparsity_bnds.dat
