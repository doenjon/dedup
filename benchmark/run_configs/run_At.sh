#!/bin/bash
#SBATCH --job-name=nf_At
#SBATCH --time=7-0
#SBATCH -p ellenyeh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

source ~/.bash_profile
conda activate env2/

nextflow main.nf -c run_configs/At.config 