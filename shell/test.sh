#!/bin/sh
##### control options #####
#PBS -N 8slab
#PBS -l nodes=4:ppn=20:full
#!!PBS -l nodes=dirac003:ppn=8:d001-004

### MAX of nodes = 8; MAX of ppn = 12 #
##########
#PBS -q full
#PBS -o out.log
#PBS -j oe

## PATH & EXE
EXE="/usr/local/vasp.5.4.4/bin_vtst_i18_ompi_noscal/vasp_std"
#EXE="/usr/local/vasp.5.4.4/bin_vtst_i18_ompi_noscal/vasp_gam"
#EXE="/usr/local/vasp.5.4.4/bin_vtst_i18_ompi_noscal/vasp_ncl"

##
NUMBER=`cat $PBS_NODEFILE | wc -l`
cd $PBS_O_WORKDIR

## run scripts
echo job started at `date` >> time
#!!---------------------------------------------------------------------------------
mpirun -np $NUMBER -machinefile $PBS_NODEFILE $EXE  > $PBS_JOBNAME.out

#!!----------------------------------------------------------------------------------
echo job ended at `date` >> time

