#!/usr/bin/env python

import os
import time
import commands
import multiprocessing
import openmoltools
import openeye.oechem as oechem
import openeye.oeomega as oeomega


def make_path(filename):
    path = os.path.split(filename)[0]
    if not os.path.exists(filename):
        os.makedirs(path)


def GenTriposGAFF(mol):
    molName = mol.GetTitle()

    # set filenames for tripos files (mol2)
    tripos_filename = os.path.join("tripos_mol2", molName+".mol2" )

    # set filenames for gaff files (mol2, frcmod, prmtop, inpcrd)
    gaff_mol2_filename = os.path.join("./gaff_mol2/", "%s.mol2" % molName)
    frcmod_filename = os.path.join("./gaff_mol2/", "%s.frcmod" % molName)
    prmtop_filename = os.path.join("./gaff_mol2/", "%s.prmtop" % molName)
    inpcrd_filename = os.path.join("./gaff_mol2/", "%s.inpcrd" % molName)

    # if files already generated, then skip
    if os.path.exists(inpcrd_filename) and os.path.exists(prmtop_filename):
        print("Already found files for %s. Skipping." % molName)
        return True    

    # calculate partial charges
    try:
        chgmol = openmoltools.openeye.get_charges(mol)
    except RuntimeError, e: # catch specified warning or Runtime error
        print("Skipping %s in %s due to RuntimeError." % (mol.GetTitle(), inputfile))
        #errfile.write("Skipping %s in %s due to RuntimeError." % (mol.GetTitle(), inputfile))
        return False
    except ValueError, e:
        print("Skipping %s in %s due to ValueError." % (mol.GetTitle(), inputfile))
        #errfile.write("Skipping %s in %s due to ValueError." % (mol.GetTitle(), inputfile))
        return False

    # reset the title to original title (get_charges renames mol)
    chgmol.SetTitle(mol.GetTitle())

    # generate tripos mol2 files
    openmoltools.openeye.molecule_to_mol2(chgmol, tripos_filename)

    # generate the gaff mol2 file and the .frcmod file
    _, _ = openmoltools.amber.run_antechamber(molName, tripos_filename,charge_method=None, gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)

    # generate the gaff .inpcrd and .prmtop files
    openmoltools.amber.run_tleap(molName, gaff_mol2_filename, frcmod_filename, prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)

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

    parser.add_option('-d','--directory',
            help = "Name of the full path directory which will contain all generated Tripos/GAFF files.",
            type = "string",
            dest = 'workdir')
    (opt, args) = parser.parse_args()
    inputfile = opt.inputfile

    # change directory to working directory
    os.chdir(opt.workdir)

    # make subdirectories for tripos and gaff files
    make_path( 'tripos_mol2/' )
    make_path( 'gaff_mol2/')

    # read input file with all molecules
    ifs = oechem.oemolistream(inputfile)
    
    # Set up an error file to record skipped troublesome molecules
    #errfile = open('errMols.txt', 'a')

    for mol in ifs.GetOEMols():

        try:
            tstart = time.clock()
            p = multiprocessing.Process(target=GenTriposGAFF, args=(mol,))
            p.start()
            p.join(25)
            if p.is_alive():
                print("Skipping %s in %s due to timeout after 25 sec" % (mol.GetTitle(), inputfile))
                #errfile.write("Skipping %s in %s due to timeout after 25 sec" % (mol.GetTitle(), inputfile))
                p.terminate()
                p.join()
        except Exception, e:
            print e, mol.GetTitle()
            continue

    # delete extra generated files
    commands.getoutput('rm ./gaff_mol2/*.frcmod')
    commands.getoutput('rm ./ANTECHAMBER*')
    commands.getoutput('rm ./ATOMTYPE.INF')
    commands.getoutput('rm ./leap.log')
    
    # close files
    ifs.close()
    #errfile.close()

