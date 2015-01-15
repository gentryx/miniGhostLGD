#! /bin/bash
#PBS -q debug
#PBS -l mppwidth=96
#PBS -l walltime=00:10:00
#PBS -N miniGhost-small
#PBS -V

cd $PBS_O_WORKDIR

# Small problem, sized for approximately 90 GB total memory usage.
# The script below is based on a 4 node, 96 MPI rank run on NERSC's Hopper
# To change to a different number of MPI ranks, change --npx, --npy and --npz.
# All other paramters would remain the same.
# E.g., an 8 MPI rank run would use: --npx 2 --npy 2 --npz 2. 

export OMP_NUM_THREADS=1

aprun -ss -n 96 ./miniGhost.x --scaling 1 --nx 672 --ny 672 --nz 672 --num_vars 40 --num_spikes 1 --debug_grid 1 --report_diffusion 21 --percent_sum 100 --num_tsteps 20 --stencil 24 --comm_method 10 --report_perf 1 --npx 4 --npy 4 --npz 6 --error_tol 8

mv results.yaml results-small.yaml
mv results.txt results-small.txt

