#!/usr/bin/env bash

#SBATCH -D /group/jrigrp9/stetter_projects/BGS_sims/
#SBATCH -o /group/jrigrp9/stetter_projects/BGS_sims/logs/simulation_out-%j.txt
#SBATCH -e /group/jrigrp9/stetter_projects/BGS_sims/logs/simulation_stderr-%j.txt
#SBATCH -t 14-00:00
#SBATCH -J bgsSim
#SBATCH --array=1-99
#SBATCH --nodes=1-1
#SBATCH --ntasks 1
#SBATCH --exclude bigmem1

echo replicate $SLURM_ARRAY_TASK_ID

#mkdir -p results/tennessen/burnins
#srun python sims/simulate_humans.py demographies/tennessen.csv results/tennessen/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/torres/burnins
#srun python sims/simulate_humans.py demographies/torres.csv results/torres/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/maize/burnins
#srun python sims/simulate_maize.py demographies/maize.csv results/maize/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/fixed_n2/
#srun python sims/fixed_n_sim.py results/fixed_n2/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/model13_long/
#srun python sims/simulate_model13_to_equilibrium.py demographies/model13_for_long_sim.csv results/model13_long/ $SLURM_ARRAY_TASK_ID

mkdir -p results/constantS/gamma2/burnins
srun python sims/simulate_generic_fixed_s.py demographies/generic_models.csv results/constantS/gamma2/ $SLURM_ARRAY_TASK_ID 2
