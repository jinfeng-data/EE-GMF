#!/bin/bash -x

#PBS -N nvt
#PBS -l nodes=1:ppn=1:gpus=1
#PBS -j oe
#PBS -q gpu

#
#define variables
#

cd $PBS_O_WORKDIR

#cuda-10.1
export CUDA_HOME=/public/software/compiler/cuda-11.3
export PATH=$PATH:$CUDA_HOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib64

export DP_HOME=/public/home/xiaohe/soft/deepmd-kit-2.1.1-gpu
export PATH=$DP_HOME/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DP_HOME/lib

lmp -in in.nvt.lammps > lmp.nvt.log

