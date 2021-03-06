CONGRATULATIONS!!
=================
You have  performed an  MD equilibration  of your  receptor, including
lipids, water molecules and counterions. For more  details on the methods  
followed please take some time to read
reference [1].


The performed equilibration includes the following stages:
----------------------------------------------------------

|   STAGE    | RESTRAINED ATOMS        | FORCE CONSTANT       | TIME           |
|:----------:|:-----------------------:|:--------------------:|:--------------:|
|  -         |   -                     |kJ/(mol?nm^2)        | ns             |
|Minimization|   -                     | -                    |(Max. 500 steps)|
|Equil. 1    |Protein Heavy Atoms      | 1000                 | 0.5            |
|Equil. 2    |Protein Heavy Atoms      | 800                  | 0.5            |
|Equil. 3    |Protein Heavy Atoms      | 600                  | 0.5            |
|Equil. 4    |Protein Heavy Atoms      | 400                  | 0.5            |
|Equil. 5    |Protein Heavy Atoms      | 200                  | 0.5            |
|Equil. 6    |Venkatakrishnan Pairs /  | 200 /                | 2.5 /          |
|            |C-alpha Atoms            | 200                  | 2.5            |

In this folder you will find several files related to this simulation:


INPUT:
------
    - popc.itp              # Topology of the lipids  
    - ffoplsaa_mod.itp      # Modified OPLSAA-FF, to account for lipid modifications  
    - ffoplsaabon_mod.itp   # Modified OPLSAA-FF(bonded), to account for lipid modifications   
    - ffoplsaanb_mod.itp    # Modified OPLSAA-FF(non-bonded), to account for lipid modifications
    - topol.tpr             # Input for the first equilibration stage
    - topol.top             # Topology of the system
    - protein.itp           # Topology of the protein
    - index.ndx             # Index file with appropriate groups for GROMACS
    - prod.mdp              # Example of a parameter file to configure a production run (see TIPS)


STRUCTURES:
-----------
    - hexagon.pdb           # Initial structure of the system, with the receptor centered in the box 
    - confout.gro           # Final structure of the system (see TIPS)
    - load_gpcr.pml         # Loads the initial structure and the trajectory in pymol


TRAJECTORY FILES:
-----------------
    - traj_pymol.xtc        # Trajectory of the whole system for visualization in pymol. 1 snapshot/100 ps
    - traj_EQ.xtc           # Trajectory of the whole system in .xtc format: 1 snapshot/50 ps 
    - ener_EQ.edr           # Energy file of the trajectory
    - load_gpcr.pml         # Script to load the equilibration trajectory in pymol.


REPORTS:
--------
In the "reports" subfolder, you will find the following files:  
    - tot_ener.xvg, tot_ener.log    # System total energy plot and log
    - temp.xvg, temp.log            # System temperature plot and log
    - pressure.xvg, pressure.log    # System pressure plot and log
    - volume.xvg, volume.log        # System volume plot and log
    - rmsd-all-atom-vs-start	    # All atoms RMSD plot
    - rmsd-backbone-vs-start.xvg    # Backbone RMSD plot
    - rmsd-calpha-vs-start.xvg      # C-Alpha RMSD plot
    - rmsf-per-residue.xvg          # Residue RMSF plot

LOGS:
-----
In the "logs" subfolder, you will find the log files of mdrun:  
    - eq_{force_constant}.log       # log of stages with restrained heavy atoms of the receptor
    - eqCA.log                      # log of the stage with restrained C-alfa atoms of the receptor


**NOTE ON GROMACS METHODS**
To integrate  the equations of  motion we have selected  the leap-frog
integrator with  a 2 femtosecond timestep.   Longe-range electrostatic
interactions  in periodic  boundary  conditions are  treated with  the
particle mesh  Ewald method.  We  use a Nose-Hoover thermostat  with a
tau_t of 0.5 picoseconds and  a Parinello-Rahman barostat with a tau_p
of 2.0.   The pressure  coupling is  semiisotropic, meaning  that it's
isotropic in the x and y  directions but different in the z direction.
Since  we are  using  pressure coupling  we are  working  with an  NPT
ensemble. This  is done both in  the all-atom restrained steps  and in
the alpha-carbon atom restrained part.   All of these details are more
explicitly stated in the Rodriguez et al. [1] publication.


**TIPS**  

NOTE: these tips work for GROMACS version >= 4.5 and < 5.0. For later 
versions, adjustments are required, but the principle remains the same.

- If you want to configure a .tpr input file for a **production** run, you
can use the template 'prod.mdp' file by introducing the number of 
steps (nsteps), and thus the simulation time, you want to run.  

After that, you just have to type:  

    grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr  
    mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c confout.gro -g production.log -x traj_prod.xtc

- If  you  want  to  create  a  PDB file  of  your  system  after  the
equilibration, with the receptor centered in the box, type:  

    echo 1 0 | trjconv -pbc mol -center -ur compact -f confout.gro -o confout.pdb

- If you want to create an xmgrace graph of the root mean square
  deviation for c-alpha atoms in the 5.0 ns of simulation you can use:  

    echo 3 3 | g_rms -f traj_EQ.xtc -s topol.tpr -o rmsd-calpha-vs-start.xvg

- You may want to get a pdb file of your last frame. You can first
check the total time of your trajectory and then use this time to
request the last frame with:

    gmxcheck -f traj_pymol.xtc
    echo 1 | trjconv -b 5000 -e 5000 -f traj_pymol.xtc -o last51.pdb


References
----------

[1] Rodr?guez D., Pi?eiro ?. and Guti?rrez-de-Ter?n H.   
Molecular Dynamics Simulations Reveal Insights into Key Structural Elements of Adenosine Receptors   
Biochemistry (2011), 50, 4194-208.   
