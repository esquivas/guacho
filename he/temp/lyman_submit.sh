#!/bin/sh 

#SBATCH --time=00:60:00
#SBATCH --nodes=2
#SBATCH -A tcast009c
#SBATCH -p DevQ
#SBATCH --mail-user=villarrc@tcd.ie
#SBATCH --mail-type=END
#SBATCH --output=ly-outfile
#SBATCH --error=ly-errfile 

module load intel/2018u4

workdir="$SLURM_SUBMIT_DIR"
cd $workdir



JOBID=`echo $PBS_JOBID | sed -e 's/\..*$//'`

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
echo "numnodes: $SLURM_NNODES"
echo "nproc: $SLURM_NPROCS"


mpirun -n 80 ./lyman_alpha_tau_angle


