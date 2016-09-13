#!/bin/bash

#              A line beginning with #PBS is a PBS directive.
#              PBS directives must come first; any directives
#                 after the first executable statement are ignored.
#PBS -N stud5_sparseMatrixProdParallel

#          Specify the number of nodes requested and the
#          number of processors per node.
#PBS -l nodes=2:ppn=1
#
###PBS -o stdout_file
###PBS -e stderr_file

#          Specify the maximum cpu and wall clock time. The wall
#          clock time should take possible queue waiting time into
#          account.  Format:   hhhh:mm:ss   hours:minutes:seconds
#PBS -l     cput=2:00:00
#PBS -l walltime=4:00:00
#          Specify the maximum amount of physical memory required per process.
#          kb for kilobytes, mb for megabytes, gb for gigabytes.
###PBS -l pmem=8192mb
##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################
NCPU=`wc -l < $PBS_NODEFILE`

cat $PBS_NODEFILE
WDIR=/scratch/$USER.$PBS_JOBID
cd $WDIR
pwd
mpirun_rsh -rsh -np $PBS_NP -hostfile $PBS_NODEFILE /home/stud5/src/sparseMatrixProdParallel
