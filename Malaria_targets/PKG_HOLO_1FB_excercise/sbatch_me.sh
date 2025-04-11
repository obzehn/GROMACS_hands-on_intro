#!/bin/bash 

#SBATCH --account=gervasio_teach_19h330
#SBATCH --partition=private-gervasio-gpu
#SBATCH --time 144:00:00
#SBATCH --job-name PKG_1FB
#SBATCH --error jobname-error.e%j
#SBATCH --output jobname-out.o%j
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --nodes 1
#SBATCH --gpus-per-node=nvidia_geforce_rtx_3080:1
#SBATCH --hint=nomultithread

module load GCC/11.3.0
module load OpenMPI/4.1.4
module load CUDA/11.7.0
module load GROMACS/2023.1-CUDA-11.7.0
export OMP_NUM_THREADS=16

runcommand="mdrun -ntmpi 1 -ntomp 16 -gpu_id 0 -pinoffset 0 -pin on -deffnm"

# energy minimization
cd ./step1_em/
if [ ! -e ./em.gro ]; then
  srun gmx grompp -f em.mdp -c ../start.gro -p ../topol.top -o em.tpr -po em_mdout.mdp
  srun gmx ${runcommand} em
fi
cd ../

# NVT equilibration
cd ./step2_nvt/
if [ ! -e ./nvt_1.gro ]; then
  srun gmx grompp -f nvt_1.mdp -c ../step1_em/em.gro -r ../step1_em/em.gro -p ../topol.top -o nvt_1.tpr -po nvt_1_mdout.mdp
  srun gmx ${runcommand} nvt_1
fi

for((i=2;i<7;i++))
do
	if [ ! -e ./nvt_${i}.gro ]; then
	  j=$((${i}-1))
	  srun gmx grompp -f nvt_${i}.mdp -c nvt_${j}.gro -r ../step1_em/em.gro -p ../topol.top -o nvt_${i}.tpr -po nvt_${i}_mdout.mdp
	  srun gmx ${runcommand} nvt_${i}
	fi
done
cd ../

# NPT equilibration
cd ./step3_npt/
if [ ! -e ./npt_1.gro ]; then
  srun gmx grompp -f npt_1.mdp -c ../step2_nvt/nvt_6.gro -r ../step1_em/em.gro -t ../step2_nvt/nvt_6.cpt -p ../topol.top -o npt_1.tpr -po npt_1_mdout.mdp
  srun gmx ${runcommand} npt_1
fi

for((i=2;i<7;i++))
do
	if [ ! -e ./npt_${i}.gro ]; then
	  j=$((${i}-1))
    srun gmx grompp -f npt_${i}.mdp -c npt_${j}.gro -r ../step1_em/em.gro -t npt_${j}.cpt -p ../topol.top -o npt_${i}.tpr -po npt_${i}_mdout.mdp
    srun gmx ${runcommand} npt_${i}
  fi
done
cd ../

# Production
cd ./step4_prod/
if [ ! -e ./prod.tpr ]; then
  srun gmx grompp -f prod.mdp -c ../step3_npt/npt_6.gro -t ../step3_npt/npt_6.cpt -p ../topol.top -o prod.tpr -po prod_mdout.mdp
fi
if [ ! -e ./prod.gro ]; then
  srun gmx ${runcommand} prod -cpi prod.cpt
fi
exit;
