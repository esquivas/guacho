#!/bin/sh 

#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH -A tcast009c
#SBATCH -p ProdQ
#SBATCH --mail-user=villarrc@tcd.ie
#SBATCH --mail-type=END
#SBATCH --output=outfile
#SBATCH --error=errfile 

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


mpirun -n 80 ./GJ436_H3


