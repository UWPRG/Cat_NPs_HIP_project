#!/bin/bash

# load gromacs
source /gscratch/pfaendtner/cnyambura/gmx2020.5_gpu.sh

#make directories in polymer/water only

counter=0

for i in 1; do	
	echo $PWD
	if [ ! -d em_ions/ ] ; then 	
		mkdir -p BSA_pH3_7_only
		cd BSA_pH3_7_only
		echo $PWD
		mkdir -p em_ions
		mkdir -p production
		mkdir -p nvt_eq
		mkdir -p npt_eq	
		mkdir -p scaling
		cp ../new.txt em_ions
		cp ../new_2.txt em_ions
		cp ../j_slurm.sh scaling
		cp ../md_scaling.mdp scaling
		cp ../minim.mdp em_ions
		cp ../nvt.mdp nvt_eq
		cp ../npt.mdp npt_eq
		cp ../md_prod.mdp production	
		cp ../bsa_pH3_7.pdb em_ions	
	fi

##generate gro and top files, 
#before running this step, make sure the polymer pdb is in the right folders, ions don't need to be added 
	
	# launch energy minimization	
	if [ ! -f em_ions/conf.gro ] && [ ! -f em_ions/posre.itp ] ; then 
		cd em_ions	
		echo $PWD	
		ln -s /gscratch/pfaendtner/cnyambura/NEE_home_mox/BSA_Nano_Prep/amber99sb-star-ildnp.ff/
		cat new.txt | gmx pdb2gmx -f bsa_pH3_7.pdb -heavyh yes
		gmx editconf -f conf.gro -o box_pH3_7.gro -c -d 1.0 -bt cubic
		gmx solvate -cp box_pH3_7.gro -cs spc216.gro -o boxsolv.gro -p topol.top 
		gmx grompp -f minim.mdp -c boxsolv.gro -p topol.top -maxwarn 1
		cat new_2.txt | gmx genion -s topol.tpr -o bsolv_ions.gro -p topol.top -nname CL -neutral
		gmx grompp -f minim.mdp -c bsolv_ions.gro -p topol.top -maxwarn 1
		gmx mdrun -nt 32 -gpu_id 0,1,2,3 -nb gpu -pme cpu -v -deffnm topol &> log.txt  
                wait ${!}
	fi	
	
	cd ../
	echo $PWD
			
	#launch NVT equilibration 
	if [ ! -f nvt_eq/confout.gro ] ; then 
		cd nvt_eq/
		echo $PWD
		ln -s ../em_ions/topol.top
                ln -s ../em_ions/topol.trr
		ln -s ../em_ions/topol.edr
		ln -s ../em_ions/topol.log
		ln -s ../em_ions/topol.gro
		ln -s ../em_ions/posre.itp
		ln -s /gscratch/pfaendtner/cnyambura/NEE_home_mox/BSA_Nano_Prep/amber99sb-star-ildnp.ff/
		gmx grompp -f nvt.mdp -c topol.gro -r topol.gro -p topol.top -maxwarn 1
		gmx mdrun -nt 32 -gpu_id 0,1,2,3 -nb gpu -pme cpu -cpi restart -cpo restart -cpt 1.0 &> log.txt   
		wait ${!}
	fi
	
	cd ../	
	echo $PWD
	
	#launch NPT equilibration 
	if [ ! -f npt_eq/confout.gro ] ; then 
		cd npt_eq/
		echo $PWD
		ln -s ../nvt_eq/confout.gro
		ln -s ../em_ions/topol.top
		ln -s ../nvt_eq/restart.cpt
		ln -s ../nvt_eq/ener.edr
		ln -s ../nvt_eq/md.log 	
                ln -s ../em_ions/posre.itp 
		ln -s /gscratch/pfaendtner/cnyambura/NEE_home_mox/BSA_Nano_Prep/amber99sb-star-ildnp.ff/	
		gmx grompp -f npt.mdp -c confout.gro -r confout.gro -t restart.cpt -p topol.top -maxwarn 1 
		rm restart.cpt
		gmx mdrun -nt 32 -gpu_id 0,1,2,3 -nb gpu -pme cpu -cpi restart -cpo restart -cpt 1.0 &> log.txt 
		wait ${!}
	fi		
	cd ../
	echo $PWD

	# Perform scaling
	cd scaling/
	echo $PWD
        for j in 1 2 4 6 8 ; do
		n_thr=$((4*${j}))
		mkdir scal_ngpu${j}_nthr${n_thr}	
		cd scal_ngpu${j}_nthr${n_thr}
		echo $PWD
		c_dir=$PWD
		cp ../j_slurm.sh .
		sed -i "s|--ntasks-per-node=XX|--ntasks-per-node=${n_thr}|g" j_slurm.sh
		sed -i "s|name=XXX|name=scal_gpu${j}_nt${n_thr}|g" j_slurm.sh
		sed -i "s|-nt XX|-nt ${n_thr}|g" j_slurm.sh
		sed -i "s|chdir=XXX|chdir=${c_dir}|g" j_slurm.sh
		sed -i "s|time=XX|time=1:00:00|g" j_slurm.sh
		sed -i "s|q6000:XX|q6000:${j}|g" j_slurm.sh
		if [ $n_thr == 4 ] ; then 
		   gp_str="0"
		   sed -i "s|-gpu_id X|-gpu_id $gp_str|g" j_slurm.sh
		elif [ $n_thr == 8 ] ; then
		   gp_str="0,1"
                   sed -i "s|-gpu_id X|-gpu_id $gp_str|g" j_slurm.sh
		elif [ $n_thr == 16 ] ; then
		   gp_str="0,1,2,3"
                   sed -i "s|-gpu_id X|-gpu_id $gp_str|g" j_slurm.sh
		elif [ $n_thr == 24 ] ; then
		   gp_str="0,1,2,3,4,5"
                   sed -i "s|-gpu_id X|-gpu_id $gp_str|g" j_slurm.sh
		elif [ $n_thr == 32 ] ; then
		   gp_str="0,1,2,3,4,5,6,7"
                   sed -i "s|-gpu_id X|-gpu_id $gp_str|g" j_slurm.sh
		fi
		cat j_slurm.sh
		cp ../md_scaling.mdp .	
		ln -s ../../npt_eq/confout.gro 
		ln -s ../../em_ions/topol.top
		ln -s ../../npt_eq/restart.cpt 
		ln -s ../../npt_eq/md.log 
		ln -s ../../npt_eq/ener.edr 	
		ln -s /gscratch/pfaendtner/cnyambura/NEE_home_mox/BSA_Nano_Prep/amber99sb-star-ildnp.ff/	
		gmx grompp -f md_scaling.mdp -t restart.cpt -c confout.gro -r confout.gro -p topol.top -maxwarn 1
		rm restart.cpt
		sbatch -p ckpt -A pfaendtner-ckpt j_slurm.sh
		wait
                cd ../
        done
        cd ../..
	counter=$((counter+1))
	echo $PWD	       
done 


