#! /bin/bash
#PBS -q regular
#PBS -l mppwidth=49152
#PBS -l walltime=00:30:00
#PBS -N miniGhost-large
#PBS -V

# Large problem, sized for approximately 46 TB total memory usage.
# The script below is based on a 2048 node, 49152 MPI rank run on NERSC's Hopper
# To change to a different number of MPI ranks, change --npx, --npy and --npz.
# All other paramters would remain the same.
# E.g., an 4096 MPI rank run would use: --npx 16 --npy 16 --npz 16. 

export OMP_NUM_THREADS=1

aprun -ss -n 49152 ./miniGhost.x --scaling 1 --nx 5376 --ny 5376 --nz 5376 --num_vars 40 --num_spikes 1 --debug_grid 1 --report_diffusion 21 --percent_sum 100 --num_tsteps 20 --stencil 24 --comm_method 10 --report_perf 1 --npx 32 --npy 32 --npz 48 --error_tol 8

mv results.yaml results-large.yaml

