#!/usr/bin/env python

# Written by: Victoria Lim and Daisy Kyu @ Mobley Lab UCI

# Purpose: Generate Tripos/GAFF/GAFF2 mol2 files for purposes of
#  comparing minimizations with different force fields. This script
#  generates input files for all the molecules listed in the single input file.
# Usage: python file.py -i /path/to/input.sdf -l /the/output/dir

import os
import time
import commands
import multiprocessing
import openmoltools
import openeye.oechem as oechem

def make_path(filename):
    """
    This function makes subdirectories.
    
    Parameters
    ----------
    filename: String name of directory
    """
    path = os.path.split(filename)[0]
    if not os.path.exists(filename):
        os.makedirs(path)

def GenTriposGAFF(mol):
    """
    This function reads in an OEChem molecule, charges it using the AM1-BCC scheme,
    and generates a Tripos mol2 file. The Tripos mol2 file is used to create GAFF 
    and GAFF2 *.mol2 files, topology files (*.prmtop), and coordinate files (*.inpcrd).
    
    Parameters
    ----------
    mol: OEChem molecule
    
    Returns
    ----------
    True if the function was successful. 
    """
    molName = mol.GetTitle()

    # set filenames for tripos files (mol2)
    tripos_filename = os.path.join("tripos_mol2", molName+".mol2" )

    # set filenames for gaff files (mol2, frcmod, prmtop, inpcrd)
    gaff_mol2_filename = os.path.join("./gaff_mol2/", "%s.mol2" % molName)
    frcmod_filename = os.path.join("./gaff_mol2/", "%s.frcmod" % molName)
    prmtop_filename = os.path.join("./gaff_mol2/", "%s.prmtop" % molName)
    inpcrd_filename = os.path.join("./gaff_mol2/", "%s.inpcrd" % molName)
    # if files already generated, then skip
    if not (os.path.exists(inpcrd_filename) and os.path.exists(prmtop_filename)):
        # calculate partial charges
        try:
            chgmol = openmoltools.openeye.get_charges(mol, keep_confs=1)
            #chgmol = openmoltools.openeye.get_charges(mol, keep_confs=1)
        except RuntimeError, e: # catch specified warning or Runtime error
            print("Skipping %s in %s due to RuntimeError." % (mol.GetTitle(), inputfile))
            return False
        except ValueError, e:
            print("Skipping %s in %s due to ValueError." % (mol.GetTitle(), inputfile))
            return False
        # reset the title to original title (get_charges renames mol)
        chgmol.SetTitle(mol.GetTitle())

        # generate tripos mol2 files
        openmoltools.openeye.molecule_to_mol2(chgmol, tripos_filename)

        # generate the gaff mol2 file and the .frcmod file
        _, _ = openmoltools.amber.run_antechamber(molName, tripos_filename,charge_method=None,\
         gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)
        # generate the gaff .inpcrd and .prmtop files
        openmoltools.amber.run_tleap(molName, gaff_mol2_filename, frcmod_filename,\
         prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)
    else:
        print("Already found gaff and tripos files for %s. Skipping." % molName)
    gaff2_inpcrd=('./gaff2_mol2/%s.inpcrd' % molName)
    gaff2_prmtop=('./gaff2_mol2/%s.prmtop' % molName)
    #generate gaff2 mol2 and frcmod files
    if not (os.path.isfile(gaff2_inpcrd) and os.path.isfile(gaff2_prmtop)):
        os.system('antechamber -i %s/%s.mol2 -fi mol2 -o %s/%s.mol2 -fo mol2 -at gaff2'\
         % ('tripos_mol2', molName, 'gaff2_mol2', molName))
        os.system('parmchk2 -i %s/%s.mol2 -f mol2 -s gaff2 -o %s/%s.frcmod'\
         % ('gaff_mol2', molName, 'gaff2_mol2', molName))

        # generates gaff2 inpcrd and prmtop files
        try:
            openmoltools.amber.run_tleap(molName, 'gaff2_mol2/%s.mol2' % molName,\
             'gaff2_mol2/%s.frcmod' % molName, prmtop_filename = 'gaff2_mol2/%s.prmtop'\
             % molName, inpcrd_filename = 'gaff2_mol2/%s.inpcrd' % molName, leaprc = 'leaprc.gaff2')
        except IOError, e:
            print ("Skipping %s in %s due to IOError." % (mol.GetTitle(), inputfile))
            return False
        except RuntimeError, e:
            print ("Skipping %s in %s due to RuntimeError for GAFF2." % (mol.GetTitle(), inputfile))
            return False
    else:
        print("Already found GAFF2 files for %s. Skipping." % molName)
    # generate gromacs .top and .gro files
    #openmoltools.utils.convert_via_acpype( molName, prmtop_filename, inpcrd_filename)
    return True


# ------------------------- Parse Inputs ----------------------- #

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-i','--input',
            help = "Name of the full path file with all molecules to generate Tripos mol2 files, GAFF mol2 files, GAFF prmtop files, and GAFF inpcrd files.",
            type = "string",
            dest = 'inputfile')

    parser.add_option('-l','--directory',
            help = "Name of the full path directory which will contain all generated Tripos/GAFF files.",
            type = "string",
            dest = 'workdir')
    (opt, args) = parser.parse_args()
    inputfile = opt.inputfile

    # change directory to working directory
    os.chdir(opt.workdir)

    # make subdirectories for tripos and gaff files
    make_path('tripos_mol2/')
    make_path('gaff_mol2/')
    make_path('gaff2_mol2/')
    # read input file with all molecules
    ifs = oechem.oemolistream(inputfile)

    ff = open("timer.dat","a")
    delta=[]
    for mol in ifs.GetOEMols():
        start = time.clock()
        try:
            tstart = time.clock()
            p = multiprocessing.Process(target=GenTriposGAFF, args=(mol,))
            p.start()
            p.join(25)
            if p.is_alive():
                print("Skipping %s in %s due to timeout after 25 sec" % (mol.GetTitle(), inputfile))
                p.terminate()
                p.join()
        except Exception, e:
            print e, mol.GetTitle()
            continue
        end = time.clock()
        delta.append( end - start)
    delta2=sum(delta)
    # Set up a timer.dat file that keeps track of how long each molecule takein seconds
    ff.write(str(inputfile)+('\t%.2f'%delta2)+'\n')
    ff.close()
    # delete extra generated files
    commands.getoutput('rm ./gaff_mol2/*.frcmod')
    commands.getoutput('rm ./gaff2_mol2/*.frcmod')
    commands.getoutput('rm ./ANTECHAMBER*')
    commands.getoutput('rm ./ATOMTYPE.INF')
    commands.getoutput('rm ./leap.log')

    # close files
    ifs.close()
