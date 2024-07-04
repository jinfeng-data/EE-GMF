export NPROCS=`wc -l hostfile |gawk '//{print $1}'`
mpirun -machinefile hostfile -np $NPROCS ./frag.x
