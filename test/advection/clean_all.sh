rm -rfv ./solution*
rm -rfv ./mesh_complete
rm -rfv ./output_init
rm -rfv ./adjoint_field
rm -v ./convergence.dat
rm -v ./coords_output.dat
rm -v ./edge_vertnums_*.dat
rm -v ./face_vertnums*.dat
rm -v ./error_calc.dat
rm -v ./IC.dat
rm -v ./arg_dict_output.jl
rm -v ./newton_solution.dat
rm -v ./timing_breakdown_*.dat
rm -v ./rk4_solution.dat
rm -v ./load_balance_*.dat
rm -v ./meshlog_*.dat

# Clean .test/Advection/newton
cd ./newton/parallel
rm -rfv ./mesh_complete
rm -rfv ./output_init
rm -rfv ./solution_done
rm -rfv ./solution_ic
rm -v ./*.dat
rm -v ./*.jl
cd ../serial
rm -rfv ./mesh_complete
rm -rfv ./output_init
rm -rfv ./solution_done
rm -rfv ./solution_ic
rm -v ./*.dat
rm -v ./arg_dict_output.jl
cd ../../





