#!/bin/bash
cd eqProd
/home/apps/apps/.spack/sandybridge/gcc-10.2.0/gromacs-2021-lx52hldze4odq56frfkrf4vdhifqb7do/bin/gmx mdrun -s topol.tpr -o traj.trr -e ener.edr -c confout.gro -g md_eqBW.log -x traj.xtc
