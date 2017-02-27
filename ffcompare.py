#!/usr/bin/env python

"""

### Description: This Python script loops over a directory of mol2 files, and
#      minimizes each with the specified force field (see supported FFs below).
#      In the directory where this script is called, a subdirectory is created
#      for the force field (if it doesn't already exist), and output mol2 files
#      are sorted into their respective subdirectories.
#   Minimizations completed using OpenEye (MMFF94*), and OpenMM (GAFF*, SMIRFF).
#   For N input *.mol2 files, there should be N output *.mol2 files, unless a
#      molecule file was corrupted or could not otherwise be minimized.
    For minimizations of different force field types, use the SAME input mol2 files.
       E.g. both GAFF and MMFF94 minimizations using Tripos mol2 input files.
       This will retain the same atom types / file format, with only coordinates
       being updated.

### Example usage:
# - python ffcompare.py --fftype smirff --ffxml smirff99Frosst.ffxml --inmols /path/to/mol2s > output.dat
# - python ffcompare.py -f gaff -i /path/to/mol2s -g /path/to/gaffFiles > output.dat

### Dependencies:
# - For all:    needs fftype and location of mol2 files. (-f, -i)
# - For GAFFx:  needs directory with both *.prmtop and *.inpcrd files (-g)
# - For SMIRFF: needs smirff99.ffxml file. (-x)


### Supported force fields:
# - MMFF94 and MMFF94S (calling 'mmff' will return both)
# - GAFF
# - GAFF2
# - Smirff99Frosst

### To do / Ideas:
# Update documentation!!
# Instead of having a different boolean variable for every FF (since right now
#  the user can only specify one at a time anyway), we could use index args.
#  I.e., 1=MMFF94*, 2=GAFF, 3=GAFF2, 4=SMIRFF. Pro: less variables, Con: script
#  would be a bit less clear (e.g., 'if dogaff' would turn into 'if fftype==1')
#
# Should we separate MMFF94 and MMFF94S?
# Should we unseparate GAFF and GAFF2?
"""

import os, sys, glob
import numpy as np

import openeye.oechem as oechem
import openeye.oeomega as oeomega
import openeye.oeszybki as oeszybki

import parmed
from parmed.openmm import topsystem
from parmed import unit as u

import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import Platform
from simtk.unit import *

from smarty import forcefield
from smarty import forcefield_utils as ff_utils




# -------------------------- Functions ------------------------- #


def writeUpdatedMol(Mol, fname):

    # Open output file to write molecule.
    ofs = oechem.oemolostream()
    if os.path.exists(fname):
        print("Output .mol2 file already exists. Skipping.\n")
        return
    if not ofs.open(fname):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % fname)

    # write out mol2 file and close the output filestream.
    oechem.OEWriteConstMolecule(ofs, Mol)
    ofs.close()
    return True

def optMMFF(Mol, FF, fname):
    """

    Take an OEMol, conduct an energy minimization, and write
       this molecule's output to a .mol2 file.
    Note: the optimization type is BFGS.

    Parameters
    ----------
    Mol: an OEChem molecule
    FF: string for OEForceFieldType to use. Either "MMFF94" or "MMFF94S"
    fname: string name of the output .mol2 file to save molecule.


    Returns
    -------
    boolean True if the function successfully completed, False otherwise

    """

    tmpmol = oechem.OEMol( Mol)  # work on a copy of the molecule

    # Open output file to write molecule.
    ofs = oechem.oemolostream()
    if os.path.exists(fname):
        print("Output .mol2 file already exists. Skipping.\n")
        return
    if not ofs.open(fname):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % fname)

    # set general energy options along with the run type specification
    optSzybki = oeszybki.OESzybkiOptions()
    optSzybki.SetSolventModel(oeszybki.OESolventModel_Sheffield)
    optSzybki.SetOptimizerType(oeszybki.OEOptType_BFGS)

    # set the particular force field
    if FF == "MMFF94":
        optSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94)
    elif FF == "MMFF94S":
        optSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94S)
    else:
        print( 'optMMFF failed for %s' %  tmpmol.GetTitle() )
        return False

    # create additional dependencies, then perform opt (in if statement)
    szOpt = oeszybki.OESzybki( optSzybki)  # generate minimization engine
    szResults = oeszybki.OESzybkiResults() # make object to hold szybki results
    if not szOpt(tmpmol, szResults):
        print( 'optMMFF failed for %s' %  tmpmol.GetTitle() )
        return False

    # write out mol2 file and close the output filestream.
    oechem.OEWriteConstMolecule(ofs, tmpmol)
    ofs.close()

    return True




def prepSMIRFF(Mol, FF_file):
    """
    Creates OpenMM Topology, System, and initial positions of given molecule.

    Parameters
    ----------
    Mol: an OEChem molecule
    FF_file: string name of *.ffxml file with path


    Returns
    -------
    Topology:  OpenMM topology for this Mol
    System:    OpenMM system for this Mol
    Positions: OpenMM positions for this Mol

    """
    # check that the *.ffxml file exists
    if not os.path.exists(FF_file):
        print("Cannot find the .ffxml file for prepSMIRFF!")
        return

    ff = forcefield.ForceField(FF_file)
    Topology, System, Positions = ff_utils.create_system_from_molecule(ff, Mol, verbose = False)
    return Topology, System, Positions


def prepGAFFx(parm):
    """
    Creates topology, system, and coordinates for AMBER
       Prmtop and Inpcrd input files. This function should work
       with either GAFF or GAFF2 files.

    Parameters
    ----------
    parm


    Returns
    -------
    Topology:  OpenMM topology for this mol's Prmtop and Inpcrd files
    System:    OpenMM system for this mol's Prmtop and Inpcrd files
    Positions: OpenMM positions for this mol's Prmtop and Inpcrd files

    """

    Topology = parm.topology
    Positions = parm.positions
    System = parm.createSystem(nonbondedMethod=app.NoCutoff)

    return Topology, System, Positions


def minimizeOpenMM(Topology, System, Positions):
    """
    Minimize molecule with specified topology, system, and positions
       using OpenMM. Return the positions of the optimized moelcule.

    Parameters
    ----------
    Topology
    System
    Positions

    Returns
    -------

    """

    # need to create integrator but don't think it's used
    integrator = mm.LangevinIntegrator(
            300.0 * u.kelvin,
            1.0 / u.picosecond,
            2.0 * u.femtosecond)

    simulation = app.Simulation(Topology, System, integrator)
    simulation.context.setPositions(Positions)
    #simulation.minimizeEnergy()
    simulation.minimizeEnergy(tolerance=5.0E-9, maxIterations=1500)

    Topology.positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)


    return Topology.positions




# --------------------------- Main Function ------------------------- #


def load_and_minimize(infiles, dommff, dosmirff, ffxml, dogaff, dogaff2, gaffdir):


    molfiles = glob.glob(os.path.join(infiles, '*.mol2'))

    ### Loop over mols.
    for f in molfiles:

        ### Load molecule from .mol2 file.
        ifs = oechem.oemolistream()
        if not ifs.open(f):
             oechem.OEThrow.Warning("Unable to open %s for reading" % f)
        try:
            mol = next(ifs.GetOEMols())
        except StopIteration:
            print('No mol loaded for %s (StopIteration exception)' % mol.GetTitle())
        ifs.close()

        if dommff:  # Optimize with MMFF94 and MMFF94S. Produce output mol2 files for each.

            mf1mol = oechem.OEGraphMol( mol )  # make a copy of mol for MMFF94
            mf2mol = oechem.OEGraphMol( mol )  # make a copy of mol for MMFF94S

            # create subdirectory for MMFF94 and MMFF94S jobs
            if not os.path.exists(os.path.join(os.getcwd(), 'MMFF94')):
                os.makedirs(os.path.join(os.getcwd(), 'MMFF94'))
            if not os.path.exists(os.path.join(os.getcwd(), 'MMFF94S')):
                os.makedirs(os.path.join(os.getcwd(), 'MMFF94S'))

            ### PART ONE: MMFF94
            print('Starting on MMFF94 optimization for %s' % mf1mol.GetTitle())
            fname = mf1mol.GetTitle()+'.mol2'
            fulln = os.path.join(os.getcwd()+'/MMFF94', fname)
            if os.path.exists(fulln):
                print('Optimization file %s for MMFF94 already exists' % (mol.GetTitle()))
                continue
            if not optMMFF(mf1mol, 'MMFF94', fulln):
                print('MMFF94 minimization failed for molecule %s:'\
                        %  (mf1mol.GetTitle()) )
                continue

            ### PART TWO: MMFF94S
            print('Starting on MMFF94S optimization for %s' % mf2mol.GetTitle())
            fname = mf2mol.GetTitle()+'.mol2'
            fulln = os.path.join(os.getcwd()+'/MMFF94S', fname)
            if os.path.exists(fulln):
                print('Optimization file %s for MMFF94S already exists' % (mol.GetTitle()))
                continue
            if not optMMFF(mf2mol, 'MMFF94S', fulln):
                print('MMFF94S minimization failed for molecule %s:'\
                        %  (mf2mol.GetTitle()) )
                continue



        if dosmirff: # Optimize with smirff99frosst. Produce output mol2 files.

            ### Set output file name and make a copy of molecule for opt.
            print('Starting on SMIRFF optimization for %s' % mol.GetTitle())
            if not os.path.exists(os.path.join(os.getcwd(), 'SMIRFF')):
                os.makedirs(os.path.join(os.getcwd(), 'SMIRFF'))
            fname = mol.GetTitle()+'.mol2'
            fulln = os.path.join(os.getcwd()+'/SMIRFF', fname)

            ### Check that optimized file doesn't already exist.
            if os.path.exists(fulln):
                print('Optimization file for %s already exists' % mol.GetTitle())
                continue

            sfmol = oechem.OEGraphMol( mol)

            ### Generate topology, system, and position.
            top, syst, pos = prepSMIRFF(sfmol, ffxml)
            parm = topsystem.load_topology(top, system=syst, xyz=pos)

            ### Use parmed to write out the mol2 file from optimized coordinates.
            parm.positions = minimizeOpenMM(top, syst, pos)

            ### Format OpenMM positions to a form readable to oemol
            concat_coords = []
            for atomi in parm.positions:
                concat_coords += [float(i) for i in atomi._value]
            sfmol.SetCoords(oechem.OEFloatArray(concat_coords))
            writeUpdatedMol(sfmol, fulln)

        if dogaff or dogaff2: # Optimize with GAFF or GAFF2. Produce output mol2 files.

            ### Set output file name and make a copy of molecule for opt.
            if not os.path.exists(os.path.join(os.getcwd(), opt.fftype.upper())):
                os.makedirs(os.path.join(os.getcwd(), opt.fftype.upper()))
            fname = "%s.mol2" % (mol.GetTitle())
            fulln = os.path.join(os.getcwd()+'/'+opt.fftype.upper(), fname)
            tmpmol = oechem.OEGraphMol( mol)
            if dogaff: print('Starting on GAFF optimization for %s' % tmpmol.GetTitle())
            if dogaff2: print('Starting on GAFF2 optimization for %s' % tmpmol.GetTitle())

            ### Get relevant files and file names.
            ffield = opt.fftype.upper()
            prmFile = os.path.join(gaffdir,mol.GetTitle()+'.prmtop')
            inpFile = os.path.join(gaffdir,mol.GetTitle()+'.inpcrd')

            ### Check: optimized file not existing, gaff dependencies present, gaff files have content
            if os.path.exists(fulln):
                print('Optimization file %s already exists' % (mol.GetTitle()))
                continue
            if not (os.path.exists(prmFile) and os.path.exists(inpFile)):
                print('%s.inpcrd or %s.prmtop files do not exist' \
 % (mol.GetTitle(), mol.GetTitle()))
                continue
            if not (os.path.getsize(prmFile) and os.path.getsize(inpFile) > 0):
                print('%s.inpcrd or %s.prmtop files are empty' \
 % (mol.GetTitle(), mol.GetTitle()))
                continue

            parm = parmed.load_file(prmFile, inpFile)

            ### Generate topology, system, and position.
            top, syst, pos = prepGAFFx(parm)

            ### Use parmed to write out mol2 file from optimized coordinates.
            parm.positions = minimizeOpenMM(top, syst, pos)

            ### Format OpenMM positions to a form readable to oemol
            concat_coords = []
            for atomi in parm.positions:
                concat_coords += [float(i) for i in atomi._value]
            tmpmol.SetCoords(oechem.OEFloatArray(concat_coords))
            writeUpdatedMol(tmpmol, fulln)


# ------------------------- Parse Inputs ----------------------- #



if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage = " ... update this from description at top of python script ... ")

    parser.add_option('-f','--fftype',
            help = "REQUIRED! Force field type. Supported options:\
 'mmff', 'gaff', 'gaff2', 'smirff'.",
            type = "string",
            dest = 'fftype')

    parser.add_option('-i', '--inmols',
            help = "REQUIRED! Path to directory containing all .mol2 files\
 to be minimized.",
            type = "string",
            dest = 'inmols')

    parser.add_option('-x', '--ffxml',
            default = None,
            help = "Force field file for smirff99frosst.\
 Required when force field type is 'smirff'.",
            type = "string",
            dest = 'ffxml')

    parser.add_option('-g', '--gaffdir',
            help = "Directory containing GAFFx parameter/topology files\
 (*.prmtop) and coordinates files (*.inpcrd). Required if and only if --fftype\
 is 'gaff' or 'gaff2'.",
            type = "string",
            default = None,
            dest = 'gaffdir')

    (opt, args) = parser.parse_args()

    ### Check required fields.
    if opt.fftype == None:
        parser.error("ERROR: No force field was specified.")
    if opt.inmols == None:
        parser.error("ERROR: No directory for input mol2 files provided.")


    ### Check that all dependencies are present.
    # For SMIRFF: needs *.ffxml (currently smirff99Frosst.ffxml)
    if opt.fftype == 'smirff':
        if opt.ffxml == None:
            parser.error("ERROR: No ffxml file provided for smirff99frosst force field.")

    # For GAFFx:  needs *.prmtop and *.inpcrd files
    if opt.fftype == 'gaff' or opt.fftype == 'gaff2':
        if opt.gaffdir == None:
            parser.error("ERROR: No GAFF* files' directory provided for GAFF* force field.")


    ### Process command line inputs.

    # Force field type.
    dommff = opt.fftype == 'mmff'      # is True if opt.type == 'mmff'
    dosmirff = opt.fftype == 'smirff'  # is True if opt.type == 'smirff'
    dogaff = opt.fftype == 'gaff'      # is True if opt.type == 'gaff'
    dogaff2 = opt.fftype == 'gaff2'    # is True if opt.type == 'gaff2'

    # File names.
    infiles = opt.inmols
    load_and_minimize(infiles, dommff, dosmirff, opt.ffxml, dogaff, dogaff2, opt.gaffdir)

