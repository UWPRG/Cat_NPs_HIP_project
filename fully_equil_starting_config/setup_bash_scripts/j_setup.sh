#!/bin/bash 

## Job Name 
#SBATCH --job-name=BSA_DS_setup

## Allocation Definition
#SBATCH --account=pfaendtner
#SBATCH --partition=pfaendtner-gpu

## Resources 
## Nodes 
#SBATCH --nodes=1

## GPUs 
#SBATCH --gres=gpu:2080ti:4

## Tasks per node (Slurm assumes you want to run 28 tasks, remove 2x # and adjust parameter if needed)
#SBATCH --ntasks-per-node=40

## Walltime (ten minutes) 
#SBATCH --time=24:00:00 

# E-mail Notification, see man sbatch for options
#SBATCH --mail-type=NONE

## Memory per node 
#SBATCH --mem=100G 

## Specify the working directory for this job 
#SBATCH --chdir=/gscratch/pfaendtner/cnyambura/NEE_home_mox/HIP_complexation/BSA_IP_sims

./BSA_DS.sh &> log_BSA_DS.txt 
