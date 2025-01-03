#!/bin/bash -x

#PBS -N NaCl3m
#PBS -l nodes=1:ppn=192
#PBS -j oe
#PBS -q amd192q

#
#define variables
#
cd $PBS_O_WORKDIR
rm -rf hostfile*
touch hostfile0 hostfile
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
cat $PBS_NODEFILE | sort | uniq >> hostfile0
sort hostfile0 > hostfile

#set temp dir
export TempDir=/scratch/$PBS_JOBID
mkdir -p $TempDir
if [ $? -ne 0 ];then
    echo "TempDir Create Failed!";exit
fi
cd $TempDir

export NPROCS=`wc -l hostfile |gawk '//{print $1}'`

source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh
OMP_NUM_THREADS=48
export OMP_NUM_THREADS
export g16root=/home/jinfeng
source $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR=$TempDir

cp $PBS_O_WORKDIR/* .
./perform_calc.x
cp -rf $TempDir/out.data $TempDir/force_all.dat $PBS_O_WORKDIR
#rm -rf $TempDir

