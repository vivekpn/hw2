#!/bin/bash
#$ -N Mandelbrot
#$ -q class8-intel
#$ -pe mpi 16
#$ -R y

# Grid Engine Notes:
# -----------------
# 1) Use "-R y" to request job reservation otherwise single 1-core jobs
#    may prevent this multicore MPI job from running.   This is called
#    job starvation.

# Module load boost
module load boost/1.57.0

# Module load OpenMPI
module load openmpi-1.8.3/gcc-4.9.2

# Run the program 
# mpirun -np 11 ./mandelbrot_ms 1000 1000


echo "Script began:" `date`
echo "Node:" `hostname`
echo "Current directory: ${PWD}"

echo ""
echo "=== Running 5 trials of Mandelbrot Master and save... ==="
for trial in 1 2 3 4 5 ; do
  echo "*** Trial ${trial} ***"
  mpirun -np 11 ./mandelbrot_ms 1000 1000
done

echo ""
echo "=== Done! ==="

# eof


