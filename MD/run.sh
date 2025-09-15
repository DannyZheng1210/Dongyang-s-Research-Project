#!/bin/bash
#SBATCH --job-name=MD
#SBATCH --nodes=1
#SBATCH --ntasks=1              
#SBATCH --cpus-per-task=4     
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --output=slurm-%j.out

set -euo pipefail

BASEDIR="${SLURM_SUBMIT_DIR:-$(pwd -P)}"

module purge
module load gromacs/2024.3-openmpi5.0.8-gcc14.2.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export GMX=gmx_mpi

echo "BASE: $BASEDIR"
echo "MPI ranks = $SLURM_NTASKS"
echo "OMP threads per rank = $OMP_NUM_THREADS"

# Input files
TOPOL="$BASEDIR/topol.top"
ASD="$BASEDIR/asd.gro"
MINIM="$BASEDIR/minim.mdp"
NVTMDP="$BASEDIR/nvt.mdp"
NPTMDP="$BASEDIR/npt.mdp"
MDMDP="$BASEDIR/md.mdp"

#creat a new folder
mkdir -p "$BASEDIR/em" "$BASEDIR/nvt" "$BASEDIR/npt" "$BASEDIR/md"

# 1) EM
cd "$BASEDIR/em"
$GMX grompp -f "$MINIM" -c "$ASD" -p "$TOPOL" -o em.tpr -maxwarn 1
mpirun -np 1 $GMX mdrun -s em.tpr -deffnm em -ntomp $OMP_NUM_THREADS -pin on

# 2) NVT
cd "$BASEDIR/nvt"
$GMX grompp -f "$NVTMDP" -c "$BASEDIR/em/em.gro" -r "$BASEDIR/em/em.gro" -p "$TOPOL" -o nvt.tpr -maxwarn 1
mpirun -np 1 $GMX mdrun -s nvt.tpr -deffnm nvt -ntomp $OMP_NUM_THREADS -pin on

# 3) NPT
cd "$BASEDIR/npt"
$GMX grompp -f "$NPTMDP" -c "$BASEDIR/nvt/nvt.gro" -t "$BASEDIR/nvt/nvt.cpt" -r "$BASEDIR/nvt/nvt.gro" -p "$TOPOL" -o npt.tpr -maxwarn 1
mpirun -np 1 $GMX mdrun -s npt.tpr -deffnm npt -ntomp $OMP_NUM_THREADS -pin on

# 4) MD
cd "$BASEDIR/md"
$GMX grompp -f "$MDMDP" -c "$BASEDIR/npt/npt.gro" -t "$BASEDIR/npt/npt.cpt" -p "$TOPOL" -o md.tpr -maxwarn 1
mpirun -np 1 $GMX mdrun -s md.tpr -deffnm md -ntomp $OMP_NUM_THREADS -pin on
