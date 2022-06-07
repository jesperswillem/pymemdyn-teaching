#!/bin/bash
cd eq
gmx mdrun -s topol.tpr -o traj.trr -e ener.edr -c confout.gro -g md_eq1000.log -x traj.xtc
