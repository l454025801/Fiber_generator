#!/bin/bash

STEP1=1
STEP2=1
STEP3=1
STEP4=0
STEP5=0

#######################################################################################################################
################ First create a new directory and a new sub-directory called "prep" ###################################
################ copy all the monomer pdbs into the new_dire/prep/.                 ###################################
################ Then copy the whole "auto_prep" directory into new_dire/           ###################################
#######################################################################################################################

# Write inp file, then generate pdb file with packmols
if (( $STEP1 == 1 )); then
	python -c 'import gen_fibers; gen_fibers.write_inp(density=9, vacancy=0.35, height=5, filler="filler_h", protonation=1)' #, ligand="Z33_GG_h", ratio=50, layers=16)'
	mv build_fiber.inp ../prep/.
	cd ../prep/.
	~/Documents/packmol/packmol < build_fiber.inp
	cd ../
fi


if (( $STEP2 == 1 )); then

	#check working directory
	#if basename $(pwd)=="auto_prep"; then
	#	cd ../
	#	echo $(pwd)
	#fi

	cp prep/fiber.pdb .
	cp system_EM.top system.top

	gmx editconf -f fiber.pdb -o box.gro -d 0
	x=$(awk 'END{print $1}' box.gro)
	y=$(awk 'END{print $2}' box.gro)
	z=$(awk 'END{print $3}' box.gro)
	x_new=$(python -c "print($x+6)")
	y_new=$(python -c "print($y+6)")
	# if 16 layers, the height should be 8	
	gmx editconf -f box.gro -o box.gro -box $x_new $y_new 8
	
	#place the tube in the center
	echo -e "r C12\nq\n" | gmx make_ndx -f box.gro -o tmp.ndx 
        echo -e "10\n0\n" | gmx trjconv -f box.gro -o box_1.gro -s box.gro -center -n tmp.ndx 
	rm tmp.ndx

	gmx solvate -cp box.gro -cs spc216 -o sol -p system
	x_center=$(python -c "print(($x+6)*5)")
	y_center=$(python -c "print(($y+6)*5)")
	echo -e "source auto_prep/remove_water.tcl\nremove_water $x_center $y_center\nexit" | vmd -dispdev text sol.gro
	
	#edit system.top to minus the number of removed water molecules
	sol_line=$(wc -l < sol.gro)
	remo_line=$(wc -l < removed_water.pdb)
	removed_water=$(python -c "print(int(($sol_line-$remo_line-1)/3))")
	
	ori_water=$(grep 'SOL' system.top | awk '{print $2}')
	new_water=$(python -c "print(int($ori_water-$removed_water))")
	sed -i "s/SOL $PARTITION_COLUMN.*/SOL                $new_water/g" system.top
fi


if (( $STEP3 == 1 )); then

	#check working directory
	#if basename $(pwd)=="auto_prep"; then
	#	cd ../
	#fi

	gmx grompp -f MDP/em.mdp -c removed_water -r removed_water -p system -o ions -maxwarn 1
	echo -e "SOL" | gmx genion -s ions -neutral -conc 0.1 -o ions -p system

	gmx grompp -f MDP/em.mdp -c ions -r ions -p system -o em -maxwarn 1
	gmx mdrun -deffnm em -v 

	echo -e "!18\nname 19 fibers\nq\n" | gmx make_ndx -f em.gro
	
fi


if (( $STEP4 == 1 )); then

	#check working directory
	#if basename $(pwd)=="auto_prep"; then
	#	cd ../
	#fi

	gmx grompp -f MDP/md1.mdp -c em -r em -p system -o md1 -n -maxwarn 1
	gmx mdrun -deffnm md1 -v
	
	#small steps npt
	gmx grompp -f MDP/md2.mdp -c md1 -r md1 -p system -o md2 -n -maxwarn 2
        gmx mdrun -deffnm md2 -v
	
	# position restraint fc = 1000
	gmx grompp -f MDP/md2_1.mdp -c md2 -r md2 -p system -o md2_1 -n -maxwarn 2
        gmx mdrun -deffnm md2_1 -v
	
	# fc = 600
	gmx grompp -f MDP/md3.mdp -c md2_1 -r md2_1 -p system -o md3 -n -maxwarn 2
        gmx mdrun -deffnm md3 -v
	
	# fc = 300
	gmx grompp -f MDP/md4.mdp -c md3 -r md3 -p system -o md4 -n -maxwarn 2
        gmx mdrun -deffnm md4 -v
fi

if (( $STEP5 == 1 )); then

	#check working directory
	#if basename $(pwd)=="auto_prep"; then
	#	cd ../
	#fi
	
	#small steps MD
	gmx grompp -f MDP/md5.mdp -c md4 -p system -o md5 -n -maxwarn 2
        gmx mdrun -deffnm md5 -v

	gmx grompp -f MDP/T298.mdp -c md5 -p system -o T298 -n -maxwarn 2
	gmx mdrun -deffnm T298 -v
fi

