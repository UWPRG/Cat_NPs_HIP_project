#!/bin/bash 

## Job Name 
#SBATCH --job-name=XXX

## Allocation Definition
#SBATCH --account=pfaendtner-ckpt
#SBATCH --partition=ckpt

## Resources 
## Nodes 
#SBATCH --nodes=1

## GPUs 
#SBATCH --gres=gpu:q6000:XX

## Tasks per node (Slurm assumes you want to run 28 tasks, remove 2x # and adjust parameter if needed)
#SBATCH --ntasks-per-node=XX 

## Walltime (ten minutes) 
#SBATCH --time=XX 

# E-mail Notification, see man sbatch for options
#SBATCH --mail-type=NONE

## Memory per node 
#SBATCH --mem=40G 

## Specify the working directory for this job 
#SBATCH --chdir=XXX

# Load cuda 11.1
module load cuda/11.1.1-1

# Load c++ compiler
module load gcc/10.1.0

# Load cmake 
module load cmake/3.11.2

# Load gromacs 2020.5
source /gscratch/pfaendtner/jpfaendt/codes/gmx2020.5/bin/GMXRC

# Execute mdrun 
gmx mdrun -nt XX -gpu_id X -nb gpu -pme cpu -cpi restart -cpo restart -cpt 1.0 &> log.txt

