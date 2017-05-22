#!/bin/bash
#$ -N Mandelbrot
#$ -q pub8i
#$ -pe mpi 64
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
#for trial in 1 2 3 4 5 ; do
echo "*** np =1 ***"
mpirun -np 1 ./mandelbrot_joe 4096 4096
echo "*** np =2 ***"
mpirun -np 2 ./mandelbrot_joe 4096 4096
echo "*** np =4 ***"
mpirun -np 4 ./mandelbrot_joe 4096 4096
echo "*** np =8 ***"
mpirun -np 8 ./mandelbrot_joe 4096 4096
echo "*** np =16 ***"
mpirun -np 16 ./mandelbrot_joe 4096 4096
echo "*** np =32 ***"
mpirun -np 32 ./mandelbrot_joe 4096 4096
echo "*** np =64 ***"
mpirun -np 64 ./mandelbrot_joe 4096 4096


#done

echo ""
echo "=== Done! ==="

# eof


