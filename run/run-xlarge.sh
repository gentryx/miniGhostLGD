#! /bin/bash
#
#PBS -l mppwidth=??
#PBS -l walltime=??

cd $PBS_O_WORKDIR

# Extra large problem, sized for approximately 900 TB total memory usage.
# To change to a different number of MPI ranks, change --npx, --npy and --npz.
# All other paramters would remain the same.
# E.g., a 80,000 MPI rank run would use: --npx 40 --npy 40 --npz 50. 
export OMP_NUM_THREADS=1
aprun -ss -n 96 ./miniGhost.x --scaling 1 --nx 14478 --ny 14478 --nz 14478 --num_vars 40 --num_spikes 1 --debug_grid 1 --report_diffusion 21 --percent_sum 100 --num_tsteps 20 --stencil 24 --comm_method 10 --report_perf 1 --npx ? --npy ? --npz ? --error_tol 8
