#!/bin/bash -l

#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -J Debuggin
#SBATCH --nodelist=node06,node07

cd /home/reidpc/GRChomboTom/Examples/ScalarFieldCosmo
srun Main_ScalarFieldCosmo3d.Linux.64.g++.gfortran.DEBUG.OPENMPCC.ex params56.txt
