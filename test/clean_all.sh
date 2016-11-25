#!/bin/bash

rm -v ./arg_dict_output.jl
rm -v ./convergence.dat
rm -v ./coords_output_*.dat
rm -v ./error_calc.dat
rm -v ./edge_vertnums.dat
rm -v ./face_vertnums.dat
rm -v ./IC_*.dat
rm -v ./iotest.dat
rm -v ./jacobian*.dat
rm -v ./load_balance_*.dat
rm -v ./log_*.dat
rm -v ./meshlog_*dat
rm -v ./newton_solution.dat
rm -v ./solution.dat
rm -v ./solution_done.dat
rm -v ./solution_ic.dat
rm -v ./sparsity_bnds.dat
rm -v ./timing_breakdown_*.dat

# Clear directories
rm -rfv ./mesh_complete
rm -rfv ./newmesh_linear
rm -rfv ./output_init
rm -rfv ./solution_done
rm -rfv ./solution_error
rm -rfv ./solution_ic

# Clean ./convergence
cd convergence/
for i in p1 #p2 p3 p4
do
	cd $i
	echo $i
	pwd
	cd ./conservative
	rm err_data.dat
	for j in m1 m2
	do
		cd $j
		pwd
		: '
		rm arg_dict_output.jl
		rm *vtu
		rm *.dat
		rm counts.txt
		'
		cd ../
	done
	cd ../

	cd ./conservative_dg
	pwd
	rm err_data.dat
	for j in m1 m2
	do
		cd $j
		pwd
		: '
		rm arg_dict_output.jl
		rm *vtu
		rm *.dat
		rm counts.txt
		'
		cd ../
	done
	cd ../

	cd ./entropy
	pwd
	rm err_data.dat
	for j in m1 m2
	do
		cd $j
		pwd
		: '
		rm arg_dict_output.jl
		rm *vtu
		rm *.dat
		rm counts.txt
		'
		cd ../
	done
	cd ../

	if [ "$i" = "p1" ] || [ "$i" = "p2" ]; then
		cd ./conservative_3dg
		pwd
		rm err_data.dat	
		for j in m1 m2
		do
			cd $j
			pwd
			: '
			rm arg_dict_output.jl
			rm *vtu
			rm *.dat
			rm counts.txt
			'
			cd ../
		done
		cd ../
	fi
	pwd

	if [ "$i" = "p1" ]; then
		cd ./entropy
		pwd
		rm err_data.dat	
		for j in m1 m2
		do
			cd $j
			pwd
			: '
			rm arg_dict_output.jl
			rm *vtu
			rm *.dat
			rm counts.txt
			'
			cd ../
		done
		cd ../
	fi
	cd ../
	pwd
done


