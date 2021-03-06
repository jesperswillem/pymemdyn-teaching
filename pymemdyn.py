#!/usr/bin/env python3.7
import argparse
import os
import shutil
import glob

import complex
import gromacs
import membrane
import protein
import pymemqueue
import settings


class Run(object):
    def __init__(self, pdb, *args, **kwargs):
        """
        A molecular dynamics *Run()*  MUST be given a *pdb* file.

        This class tries to initialize a full complex to send to simulation.
        Given a set of molecules (protein, ligand, other ligand, waters, ...),
        this class would try to build a full embedded-in-membrane complex.

        The complex is stored in self.g (a *Gromacs* object), and thus
        can be **run** through g.recipe and g.run_recipe procedure. See
        gromacs.py for more information.

        The queueing system is also created here to be used in certain steps.
        """

        self.pdb = pdb
        self.own_dir = kwargs.get("own_dir") or ""
        self.repo_dir = kwargs.get("repo_dir") or ""
        self.ligand = kwargs.get("ligand") or ""
        self.alosteric = kwargs.get("alosteric") or ""
        self.waters = kwargs.get("waters") or ""
        self.ions = kwargs.get("ions") or ""
        self.cho = kwargs.get("cho") or ""
        self.ligpargen = kwargs.get("ligpargen") or ""
        self.library = kwargs.get("library") or ""
        self.restraint = kwargs.get("restraint") or ""
        self.queue = kwargs.get("queue") or ""
        self.debug = kwargs.get("debug") or False

        if self.pdb:
            self.pdb = protein.Protein(pdb=self.pdb).check_number_of_chains()

        sugars = {"ligand": "Ligand",
                  "alosteric": "Alosteric",
                  "waters": "CrystalWaters",
                  "ions": "Ions",
                  "cho": "Cholesterol"}

        if self.ligpargen or self.library:
            protein.Sugar_prep.__init__(self)

        for sugar_type, class_name in sugars.items():
            if getattr(self, sugar_type):
                base_name = getattr(self, sugar_type)
                setattr(self,
                        sugar_type,
                        getattr(protein, class_name)(
                            pdb=base_name + ".pdb",
                            itp=base_name + ".itp",
                            ff=base_name + ".ff"))

        self.membr = membrane.Membrane()

        prot_complex = protein.ProteinComplex(
            monomer=self.pdb,
            ligand=self.ligand or None,
            alosteric=self.alosteric or None,
            waters=self.waters or None,
            ions=self.ions or None,
            cho=self.cho or None)

        full_complex = complex.MembraneComplex()

        if self.pdb.__class__.__name__ == "Dimer":
            '''The box for the dimers is slightly bigger'''
            full_complex.box_height = 3.5
            full_complex.box_width = 1.2
        full_complex.complex = prot_complex
        full_complex.membrane = self.membr

        self.g = gromacs.Gromacs(membrane_complex=full_complex)

        # NOTE: If not provided in command line, self.queue is set to
        # NoQueue
        if self.queue:
            if self.queue == "slurm":
                my_queue = queue.Slurm()
            elif self.queue == "pbs":
                my_queue = queue.PBS()
            elif self.queue == "pbs_ib":
                my_queue = queue.PBS_IB()
            elif self.queue == "svgd":
                my_queue = queue.Svgd()
        else:
            my_queue = pymemqueue.NoQueue()

        self.g.queue = my_queue

    def clean(self):
        """
        Removes all previously generated files
        """
        to_unlink = ["#index.ndx.1#", "#index.ndx.2#", "#index.ndx.3#", 
                     "#index.ndx.4#", "#index.ndx.5#", "#output.pdb.1#", 
                     "#proteinopls.pdb.1#", "#proteinopls.pdb.2#", 
                     "#proteinopls.pdb.3#", "#proteinopls.pdb.4#", 
                     "#protpopc.pdb.1#", "#tmp.pdb.1#", "#topol.top.1#", 
                     "#topol.top.2#", "#topol.top.3#", "#topol.tpr.1#", 
                     "#topol.tpr.2#", "#topol.tpr.3#", "#topol.tpr.4#", 
                     "disre.itp", "ener_EQ.edr", "ffoplsaa_mod.itp", 
                     "ffoplsaabon_mod.itp", "ffoplsaanb_mod.itp", 
                     "hexagon.pdb", "index.ndx", "ions.itp", "ligand_ha.ndx", 
                     "MD_output.tgz", "mdout.mdp", "mdrun.sh", "min.pdb", 
                     "output.pdb", "popc.gro", "popc.itp", "popc.pdb", 
                     "posre.itp", "posre_lig.itp", "pressure.log", 
                     "pressure.xvg", "protein.itp", "protein.top", 
                     "protein_ca200.itp", "proteinopls.fasta", 
                     "proteinopls.pdb", "proteinopls_bw.aln", 
                     "proteinopls_CA.pdb", "protpopc.pdb", 
                     "rmsd-all-atom-vs-start.xvg", 
                     "rmsd-backbone-vs-start.xvg", "rmsd-calpha-vs-start.xvg", 
                     "rmsf-per-residue.xvg", "spc.itp", "steep.mdp", 
                     "temp.log", "temp.xvg", "tmp.pdb", "tmp_proteinopls.pdb",
                     "topol.top", "topol.tpr", "tot_ener.log", "tot_ener.xvg",
                     "traj_EQ.xtc", "traj_pymol.xtc", "volume.log", 
                     "volume.xvg", "water.gro", "water.pdb","GROMACS.log"]

        dirs_to_unlink = ["Rmin", "eq", "eqProd"]

        for target in to_unlink:
            if os.path.isfile(target): os.unlink(target)

        for target in dirs_to_unlink:
            if os.path.isdir(target): shutil.rmtree(target)

        for backup in glob.glob('#*#'):
            os.unlink(backup)

        #for backup in glob.glob('*_backup*'):
        #    os.unlink(backup)

        return True

    def moldyn(self):
        """
        Run all steps in a molecular dynamics simulation of a membrane protein
        """
        if self.restraint == "bw":
            steps = ["Init", "Minimization", "Equilibration", "Relax", 
                     "BWRelax", "BWCollectResults"]
        elif self.restraint == "ca":
            steps = ["Init", "Minimization", "Equilibration", "Relax", 
                     "CARelax", "CACollectResults"]

        for step in steps:
            self.g.select_recipe(stage=step, debug=self.debug)
            self.g.run_recipe(debug=self.debug)

    def moldyn_notebook_info(self,stage="Init"):
        """
        Run all steps in a molecular dynamics simulation of a membrane protein
        """
        self.g.run_recipe_info(debug=self.debug)

    def moldyn_notebook_run(self, command_name, stage="Init"):
        self.g.run_recipe_notebook(command_name)#, debug=self.debug)

    def light_moldyn(self):
        """
        This is a function to debug a run in steps
        """
        steps = ["BWRelax", "CollectResults"]
#        steps = ["CollectResults"]

        for step in steps:
            self.g.select_recipe(stage=step, debug=self.debug)
            self.g.run_recipe(debug = self.debug)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='pymemdyn',
        description=' == Setup Molecular Dynamics for Membrane Proteins given a PDB. == ')

    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s 1.5.1')

    parser.add_argument('-b',
                        dest = "own_dir",
                        help = "Working dir if different from actual dir",
                        default = os.getcwd())

    parser.add_argument('-r',
                        dest = "repo_dir",
                        help = "Path to templates of fixed files. If not \
            provided, take the value from settings.TEMPLATES_DIR.",
                        default = settings.TEMPLATES_DIR)

    parser.add_argument('-p',
                        dest = "pdb",
                        required = True,
                        help = "Name of the pdb to insert into membrane for MD (mandatory). \
            Use the pdb extension. (e.g. -p myprot.pdb)")

    parser.add_argument('-l', '--lig',
                        dest = "ligand",
                        help = "Name of the ligand, without extension. Three files must be \
            present along with the molecule pdb: the ligand, its itp and \
            its force field. (See --lpg and --lib for more information)")

    parser.add_argument('-a', "--alo",
                        dest = "alosteric",
                        help = "Name of the alosteric interaction, without extension. Three \
            files must be present along with the molecule pdb: the alosteric, \
            its itp and its force field. (See --lpg and --lib for more information)")

    parser.add_argument('-w','--waters',
                        dest = "waters",
                        help = "Crystalized water molecules. File name without extension.")

    parser.add_argument('-i', '--ions',
                        dest = "ions",
                        help = "Crystalized ions file name without extension.")

    parser.add_argument('-c', '--cho',
                        dest = "cho",
                        help = "Crystalized cholesterol molecules file name\
            without extension.")
            
    parser.add_argument('--lpg',
                        dest = "ligpargen",
                        help = "Indicate which ligand or cofactor pdb and topology files \
            are generated using LigParGen. Options: any combination of l (ligand), a \
            (allosteric), and c (cholesterol)")
    
    parser.add_argument('--lib',
                        dest = "library",
                        help = "Indicate which ligand or cofactor itp and ff you wish to \
            retrieve from the PyMemDyn library. Options: any combination of l (ligand), a \
            (allosteric), w (waters), i (ions), and c (cholesterol)")

    parser.add_argument('--res',
                        dest = "restraint",
                        help = "Position restraints during MD production run. Options: bw \
            (Ballesteros-Weinstein Restrained Relaxation - default), ca (C-Alpha Restrained \
            Relaxation)",
                        default = "bw")

    parser.add_argument('-q', '--queue',
                        dest = "queue",
                        help = "Queueing system to use (slurm, pbs, pbs_ib and svgd supported)",
                        default = "")

    parser.add_argument('-d', '--debug',
                        action="store_true")

    args = parser.parse_args()

    if not (os.path.isdir(args.own_dir)):
        os.makedirs(args.own_dir)
        print ("Created working dir {0}".format(args.own_dir))
    os.chdir(args.own_dir)

    # Removes files previously generated in Run()
    to_unlink = [''.join([args.pdb[:-4], "-his.pdb"]), 
                 ''.join([args.pdb, "~"])]
   
    for target in to_unlink:
        if os.path.isfile(target): os.unlink(target)

    run = Run(own_dir = args.own_dir,
              repo_dir = args.repo_dir,
              pdb = args.pdb,
              ligand = args.ligand,
              alosteric = args.alosteric,
              waters = args.waters,
              ions = args.ions,
              cho = args.cho,
              ligpargen = args.ligpargen,
              library = args.library,
              restraint = args.restraint,
              queue = args.queue,
              debug = args.debug)


    run.clean()
    # Rewrite old GROMACS.log file if on a re-run
    f = open("GROMACS.log", "w")
    f.close()


    run.moldyn()
#    run.light_moldyn()