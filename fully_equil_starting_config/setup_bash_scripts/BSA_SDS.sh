#!/bin/bash

# load gromacs
source /gscratch/pfaendtner/cnyambura/gmx2020.5_gpu.sh

#make directories in polymer/water only

counter=0

# Box size in nm
n_size=15

#Half box length for packmol
nb_hlen=74

declare -a pos_ions
pos_ions=(0 0 0 50)

declare -a neg_ions
neg_ions=(62 46 14 0)

mkdir -p BSA_SDS_ions
cd BSA_SDS_ions

for i in 16 32 64 128 ; do	
	echo $PWD
	mkdir BSA_SDS_MR_${i}
	cd BSA_SDS_MR_${i}
	for j in 1 2 3 ; do 
		if [ ! -d em_ions/ ] ; then 	
			mkdir trial_${j}
			cd trial_${j}
			mkdir -p em_ions
			mkdir -p production
			mkdir -p nvt_eq
			mkdir -p npt_eq	
			mkdir -p scaling
			cp ../../../new.txt em_ions
			cp ../../../j_slurm.sh scaling
			cp ../../../md_scaling.mdp scaling
			cp ../../../minim.mdp em_ions
			cp ../../../sds_posre.itp em_ions
			cp ../../../HMR_sds.itp em_ions
			cp ../../../SDS.pdb em_ions		
			cp ../../../K.pdb em_ions
			cp ../../../CL.pdb em_ions
			cp ../../../nvt.mdp nvt_eq
			cp ../../../npt.mdp npt_eq
			cp ../../../md_prod.mdp production	
			cp ../../../bsa_pH3_7.pdb em_ions	
		fi

	##generate gro and top files, 
	#before running this step, make sure the polymer pdb is in the right folders, ions don't need to be added 
	
		# launch energy minimization	
		if [ ! -f em_ions/conf.gro ] && [ ! -f em_ions/posre.itp ] ; then 
			cd em_ions	
			echo $PWD	
			ln -s /gscratch/pfaendtner/cnyambura/NEE_home_mox/BSA_Nano_Prep/amber99sb-star-ildnp.ff/
			cat new.txt | gmx pdb2gmx -f bsa_pH3_7.pdb -heavyh yes
			rm conf.gro
			if [ ${pos_ions[counter]} -ne 0 ] ; then 
				cat > packmol.inp << EOF1
#
# BSA/SDS/K/CL in water system

seed -1

tolerance 2.0
filetype pdb
output BSA_SDS_MR${i}t${j}.pdb

structure bsa_pH3_7.pdb
  number 1
  center
  fixed 0. 0. 0. 0. 0. 0.
end structure

structure SDS.pdb
  number ${i}
  inside box -${nb_hlen} -${nb_hlen} -${nb_hlen} ${nb_hlen} ${nb_hlen} ${nb_hlen}
  outside sphere 0. 0. 0. 50.
end structure

structure K.pdb
  number ${pos_ions[counter]}
  inside box -${nb_hlen} -${nb_hlen} -${nb_hlen} ${nb_hlen} ${nb_hlen} ${nb_hlen}
  outside sphere 0. 0. 0. 50.
end structure
EOF1

                	        packmol < packmol.inp
                	        gmx editconf -f BSA_SDS_MR${i}t${j}.pdb -c -box ${n_size} ${n_size} ${n_size} -o MR${i}_bsaSDSt${j}.gro
                	        sed -i '24i\; Include SDS parameters' topol.top 
				sed -i '25i\#include "HMR_sds.itp"' topol.top
				echo "SDS            ${i}" >> topol.top
				echo "K            ${pos_ions[counter]}" >> topol.top	
 			       	gmx solvate -cp MR${i}_bsaSDSt${j}.gro -cs spc216.gro -o boxsolvSDS_MR${i}t${j}.gro -p topol.top 
 			       	gmx grompp -f minim.mdp -c boxsolvSDS_MR${i}t${j}.gro -p topol.top -maxwarn 1
				gmx mdrun -nt 40 -gpu_id 0,1,2,3 -nb gpu -pme cpu -v -deffnm topol &> log.txt
				wait ${!}
			elif [ ${pos_ions[counter]} -eq 0 ] ; then
				cat > packmol.inp << EOF1
#
# BSA/SDS/K/CL in water system

seed -1

tolerance 2.0
filetype pdb
output BSA_SDS_MR${i}t${j}.pdb

structure bsa_pH3_7.pdb
  number 1
  center
  fixed 0. 0. 0. 0. 0. 0.
end structure

structure SDS.pdb
  number ${i}
  inside box -${nb_hlen} -${nb_hlen} -${nb_hlen} ${nb_hlen} ${nb_hlen} ${nb_hlen}
  outside sphere 0. 0. 0. 50.
end structure

structure CL.pdb
  number ${neg_ions[counter]}
  inside box -${nb_hlen} -${nb_hlen} -${nb_hlen} ${nb_hlen} ${nb_hlen} ${nb_hlen}
  outside sphere 0. 0. 0. 50.
end structure
EOF1

                                packmol < packmol.inp
                                gmx editconf -f BSA_SDS_MR${i}t${j}.pdb -c -box ${n_size} ${n_size} ${n_size} -o MR${i}_bsaSDSt${j}.gro
                                sed -i '24i\; Include SDS parameters' topol.top
                                sed -i '25i\#include "HMR_sds.itp"' topol.top
                                echo "SDS            ${i}" >> topol.top
                                echo "CL            ${neg_ions[counter]}" >> topol.top
                                gmx solvate -cp MR${i}_bsaSDSt${j}.gro -cs spc216.gro -o boxsolvSDS_MR${i}t${j}.gro -p topol.top
                                gmx grompp -f minim.mdp -c boxsolvSDS_MR${i}t${j}.gro -p topol.top -maxwarn 1
                                gmx mdrun -nt 40 -gpu_id 0,1,2,3 -nb gpu -pme cpu -v -deffnm topol &> log.txt
                                wait ${!}
			fi
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
			ln -s ../em_ions/sds_posre.itp
			ln -s ../em_ions/HMR_sds.itp
			ln -s /gscratch/pfaendtner/cnyambura/NEE_home_mox/BSA_Nano_Prep/amber99sb-star-ildnp.ff/
			gmx grompp -f nvt.mdp -c topol.gro -r topol.gro -p topol.top -maxwarn 1
			gmx mdrun -nt 40 -gpu_id 0,1,2,3 -nb gpu -pme cpu -cpi restart -cpo restart -cpt 1.0 &> log.txt   
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
			ln -s ../em_ions/sds_posre.itp
			ln -s ../em_ions/HMR_sds.itp 
			ln -s /gscratch/pfaendtner/cnyambura/NEE_home_mox/BSA_Nano_Prep/amber99sb-star-ildnp.ff/	
			gmx grompp -f npt.mdp -c confout.gro -r confout.gro -t restart.cpt -p topol.top -maxwarn 1 
			rm restart.cpt
			gmx mdrun -nt 40 -gpu_id 0,1,2,3 -nb gpu -pme cpu -cpi restart -cpo restart -cpt 1.0 &> log.txt 
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
			sed -i "s|time=XX|time=0:15:00|g" j_slurm.sh
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
			ln -s ../../em_ions/HMR_sds.itp	
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
        	cd ../../
	
	echo $PWD	       
	done 
	cd ../
	echo $counter
	counter=$((counter+1))
	echo $PWD	
done

