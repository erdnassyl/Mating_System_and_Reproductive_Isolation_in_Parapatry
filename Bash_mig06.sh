#!/bin/bash
#SBATCH --job-name=slim
#SBATCH --account=irbi
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --array=1-400
#SBATCH --output=slurm_%A_%a.out

module load anaconda/2024.10-1
source /shared/apps/anaconda/2024.10-1/anaconda3/etc/profile.d/conda.sh
conda activate slim_env

n_iter=100


sr_values=(0 0.3 0.6 0.9)


nsr=${#sr_values[@]}

index=$((SLURM_ARRAY_TASK_ID - 1))


sr_index=$(( index % nsr ))
sr=${sr_values[$sr_index]}

slim -d seed=$SLURM_ARRAY_TASK_ID \
     -d L=1e5 \
     -d N=1000 \
     -d gamma=0.04 \
     -d mig=0.06 \
     -d sr=$sr \
     -d mu_m1=1e-9 \
     -d mu_m2=1e-9 \
     -d mu_m3=1e-9 \
     -d s_m1=0 \
     -d s_m2=0 \
     -d s_m3=0 \
     -d h_m1=0.5 \
     -d h_m2=0.5 \
     -d h_m3=0.5 \
     -d rec=1e-5 \
     SLiM_BDMI_MainTest06.txt

