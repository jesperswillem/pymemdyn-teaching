#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 0-08:00:00

ml anaconda

python /home/jespers/software/pymemdyn/pymemdyn -p protein.pdb -l E3R --lpg l
