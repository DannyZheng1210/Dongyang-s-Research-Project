#!/bin/bash -l
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job.
#SBATCH --export=ALL
# Define job name
#SBATCH -J frame_0142_frames1_frame0089
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o castep.%u.%N.%j.out
# Define a standard error file
#SBATCH -e castep.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes
# Request the number of nodes
#SBATCH -N 1
# Specify the number of tasks per node 
#SBATCH --ntasks-per-node=64
# Request the number of cpu per task
#SBATCH --cpus-per-task=1
# This asks for 1 hrs
#SBATCH -t 50:00:00
# Insert your own username to get e-mail notifications
#SBATCH --mail-user=sgdzheng@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL
#
export OMP_NUM_THREADS=1
module purge
module load castep/24.1-openmpi5.0.5-gcc14.2.0
module list
echo =========================================================   
echo SLURM job: submitted  date = `date`
date_start=`date +%s`

export TMPDIR=$HOME/localscratch

hostname
echo Current directory: `pwd`

echo "Print the following environmetal variables:"
echo "Job name                     : $SLURM_JOB_NAME"
echo "Job ID                       : $SLURM_JOB_ID"
echo "Job user                     : $SLURM_JOB_USER"
echo "Job array index              : $SLURM_ARRAY_TASK_ID"
echo "Submit directory             : $SLURM_SUBMIT_DIR"
echo "Temporary directory          : $TMPDIR"
echo "Submit host                  : $SLURM_SUBMIT_HOST"
echo "Queue/Partition name         : $SLURM_JOB_PARTITION"
echo "Node list                    : $SLURM_JOB_NODELIST"
echo "Hostname of 1st node         : $HOSTNAME"
echo "Number of nodes allocated    : $SLURM_JOB_NUM_NODES or $SLURM_NNODES"
echo "Number of tasks              : $SLURM_NTASKS"
echo "Number of tasks per node     : $SLURM_TASKS_PER_NODE"
echo "Initiated tasks per node     : $SLURM_NTASKS_PER_NODE"
echo "Requested CPUs per task      : $SLURM_CPUS_PER_TASK"
echo "Requested CPUs on the node   : $SLURM_CPUS_ON_NODE"
echo "Scheduling priority          : $SLURM_PRIO_PROCESS"

echo "Running parallel job:"
echo Job output begins                                           
echo -----------------                                           
echo "Running castep:"

## The next few environment variables ensure good OpenMP performance
## in the MPI job

export OMP_PROC_BIND=true
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp

mpirun  -np $SLURM_NTASKS castep.mpi  frame_0142_frames1_frame0089

ret=$?
echo   
echo ---------------                                           
echo Job output ends                                           
date_end=`date +%s`
seconds=$((date_end-date_start))
minutes=$((seconds/60))
seconds=$((seconds-60*minutes))
hours=$((minutes/60))
minutes=$((minutes-60*hours))
echo =========================================================   
echo SLURM job: finished   date = `date`
echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
echo =========================================================   
exit $ret
