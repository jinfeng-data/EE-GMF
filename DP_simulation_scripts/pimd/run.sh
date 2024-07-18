#!/bin/bash -x

#PBS -N 290k-cmd
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

export PATH="/public/home/xiaohe/soft/amber18/miniconda/bin:$PATH"

export DP_HOME=/public/home/xiaohe/soft/deepmd-kit-2.1.1-gpu
export PATH=$DP_HOME/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DP_HOME/lib
export PATH=/public/home/xiaohe/jfmao/soft/miniconda3/bin:$PATH
source /public/home/xiaohe/jinfeng/i-pi-master/env.sh

rm -rf /tmp/ipi_localhost11
i-pi input.xml > ipi.log & sleep 10
lmp -in in.lmp > lmp.log
