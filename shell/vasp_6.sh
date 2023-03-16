#!/bin/sh
# control options #
#PBS -N vasp_6 
#PBS -l nodes=1:ppn=36:t2
########
#PBS -q t2
#PBS -o out.log
#PBS -j oe

# PATH & EXE
EXE='/usr/local/vasp.6.3.2/bin_i18_impi_avx512/vasp_std'
#EXE='/usr/local/vasp.6.3.2/bin_i18_impi_avx512/vasp_ncl'
#EXE='/usr/local/vasp.6.3.2/bin_i18_impi_avx512/vasp_gam'

#
NUMBER=`cat $PBS_NODEFILE | wc -l`
cd $PBS_O_WORKDIR

# run 
echo job started at `date` >> time
mpirun -np $NUMBER -machinefile $PBS_NODEFILE $EXE  > $PBS_JOBNAME.out
echo job ended at `date` >> time
