#!/bin/bash -f
#SBATCH --time=4:00:00
#SBATCH -A ccsi
#SBATCH -p batch
#SBATCH --mem=32G
#SBATCH --ntasks-per-node 32
#SBATCH --job-name=run.makezones
#SBATCH --mail-user yuanf@ornl.gov
#SBATCH -o ./%j-output.txt
#SBATCH -e ./%j-error.txt
#SBATCH -N 1
#SBATCH -n 7
#SBATCH -c 1

cd ./

if mpirun -np 7 ./makezones
then
  wait
  echo 'successfully DONE!'

else
  exit &?
fi

