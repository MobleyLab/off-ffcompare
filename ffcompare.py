#!/usr/bin/env python


### Description: This Python script takes one SDF file (containing
#      a bunch of molecules) with specified input parameters and 
#      runs energy minimzations on each of those molecules. 
#      A mol2 with optimized coordinates is output for each
#      molecule in the SDF file. (N molecules --> N *.mol2 files)
#      Relevant FF subdirectories are created in the directory
#      that houses this python script.


### Use this script with command line inputs in slurm script (bash).
#     Ex. python ffcompare.py --fftype smirff --ffxml smirff99Frosst.ffxml --inmols /path/to/mol2s



### Dependencies
# For all:    needs fftype and location of mol2 files.
# For GAFFx:  needs *.prmtop and *.inpcrd files
# For SMIRFF: needs smirff99.ffxml file.


### Supported force fields
# MMFF94 and MMFF94S (calling 'mmff' will return both)
# GAFF
# GAFF2
# Smirff99Frosst




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
    szOpt = oeszybki.OESzybki( optSzybki)   # generate minimization engine
    szResults = oeszybki.OESzybkiResults()  # make object to hold szybki results
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
    boolean True if the function successfully completed, False otherwise
    
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
    
    Topology.positions = simulation.context.getState(getPositions=True).getPositions()
    
    return Topology.positions




# --------------------------- Main Function ------------------------- #


def load_and_minimize(infiles, dommff, dosmirff, ffxml, dogaff, dogaff2, prmFile, inpFile):


    molfiles = glob.glob(os.path.join(infiles, '*.mol2'))
    
    ### Loop over mols.
    for f in molfiles:

        ### Load molecule from .mol2 file.
        ifs = oechem.oemolistream()
        if not ifs.open(f):
             oechem.OEThrow.Warning("Unable to open %s for reading" % f)
        mol = next(ifs.GetOEMols())
        ifs.close()
    
        if dommff:  # Optimize with MMFF94 and MMFF94S. Produce output mol2 files for each.
    
            mf1mol = oechem.OEGraphMol( mol )    # make a copy of mol for MMFF94
            mf2mol = oechem.OEGraphMol( mol )    # make a copy of mol for MMFF94S

            # create subdirectory for MMFF94 and MMFF94S jobs
            if not os.path.exists(os.path.join(os.getcwd(), 'MMFF94')):    
                os.makedirs(os.path.join(os.getcwd(), 'MMFF94'))
            if not os.path.exists(os.path.join(os.getcwd(), 'MMFF94S')):
                os.makedirs(os.path.join(os.getcwd(), 'MMFF94S'))
    
            ### PART ONE: MMFF94
            print('Starting on MMFF94 optimization for %s' % mf1mol.GetTitle())
            fname = mf1mol.GetTitle()+'.mol2'
            fulln = os.path.join(os.getcwd()+'/MMFF94', fname)
            if not optMMFF(mf1mol, 'MMFF94', fulln):
                print('MMFF94 minimization failed for molecule %s:'\
                        %  (mfmol.GetTitle()) )
                continue
    
            ### PART TWO: MMFF94S
            print('Starting on MMFF94S optimization for %s' % mf2mol.GetTitle())
            fname = mf2mol.GetTitle()+'.mol2'
            fulln = os.path.join(os.getcwd()+'/MMFF94S', fname)
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
            sfmol = oechem.OEGraphMol( mol)
    
            ### Generate topology, system, and position.
            top, sys, pos = prepSMIRFF(sfmol, ffxml)
            parm = topsystem.load_topology(top, system=sys, xyz=pos)
    
            ### Use parmed to write out the mol2 file from optimized coordinates.
            parm.positions = minimizeOpenMM(top, sys, pos)
            mol2f = parmed.formats.Mol2File
            mol2f.write(parm,fulln)
        
        if dogaff: # Optimize with GAFF or GAFF2. Produce output mol2 files.
    
            ### Set output file name and make a copy of molecule for opt.
            if not os.path.exists(os.path.join(os.getcwd(), opt.fftype.upper())):    
                os.makedirs(os.path.join(os.getcwd(), opt.fftype.upper())) 
            fname = "%s.mol2" % (mol.GetTitle()) 
            fulln = os.path.join(os.getcwd()+'/'+opt.fftype.upper(), fname)
            tmpmol = oechem.OEGraphMol( mol)
            print('Starting on GAFFx optimization for %s' % tmpmol.GetTitle())
    
            ### Get relevant files and file names.
            ffield = opt.fftype.upper()
            parm = parmed.load_file(prmFile, inpFile)
    
            ### Generate topology, system, and position.
            top, sys, pos = prepGAFFx(pmfFile, inpFile)
    
            ### Use parmed to write out the mol2 file from optimized coordinates.
            parm.positions = minimizeOpenMM(top, sys, pos)
            mol2f = parmed.formats.Mol2File
            mol2f.write(parm,fulln)
    
    

# ------------------------- Parse Inputs ----------------------- #



if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage = "This script sets up and optimizes molecules in OpenMM applicable to a variety of several force field types.\nReads in a single SDF file containing all molecules,\nand writes out individual mol2 files per molecule.\nExample: python ffcompare.py --fftype mmff --inmols /path/to/mol2s")

    parser.add_option('-f','--fftype',
            help = "REQUIRED! Force field type. Supported options are 'mmff', 'gaff', 'gaff2', 'smirff'.",
            type = "string",
            dest = 'fftype')

    parser.add_option('-i', '--inmols',
            help = "REQUIRED! Path to all directory containing all (no more, no less) .mol2 files.",
            type = "string",
            dest = 'inmols')

    parser.add_option('-x', '--ffxml',
            default = None,
            help = "Force field file for smirff99frosst. Needed when force field type is 'smirff'.",
            type = "string",
            dest = 'ffxml')

    parser.add_option('-p', '--prmtop',
            help = "GAFFx parameter/topology file. Needed when force field type is 'gaff' or 'gaff2'.",
            type = "string",
            default = None,
            dest = 'prmtop')

    parser.add_option('-c','--inpcrd',
            help = "GAFFx coordinate file. Needed when force type is 'gaff' or 'gaff2'.",
            type = "string",
            default = None,
            dest = 'inpcrd')

    (opt, args) = parser.parse_args()

    ### Check required fields.
    if opt.fftype == None:
        parser.error("ERROR: No force field was specified.")
    if opt.inmols == None:
        parser.error("ERROR: No directory location for mol2 files was provided.")


    ### Check that all dependencies are present.
    # For SMIRFF: needs *.ffxml (currently smirff99Frosst.ffxml) file and *.sdf file.
    if opt.fftype == 'smirff':
        if opt.ffxml == None:
            parser.error("ERROR: No ffxml file provided for smirff99frosst force field.")

    # For GAFFx:  needs *.prmtop and *.inpcrd files
    if opt.fftype == 'gaff' or opt.fftype == 'gaff2':
        if opt.prmtop == None:
            parser.error("ERROR: No parameter/topology file provided for GAFF or GAFF2 force field.")
        if opt.inpcrd == None:
            parser.error("ERROR: No coordinate file provided for GAFF or GAFF2 force field.")
    
    
    ### Process command line inputs.
    
    # Force field type.
    dommff = opt.fftype == 'mmff'      # is True if opt.type == 'mmff'
    dosmirff = opt.fftype == 'smirff'  # is True if opt.type == 'smirff'
    dogaff = opt.fftype == 'gaff'      # is True if opt.type == 'gaff'
    dogaff2 = opt.fftype == 'gaff2'    # is True if opt.type == 'gaff2'
    
    # File names.
    infiles = opt.inmols

    load_and_minimize(infiles, dommff, dosmirff, opt.ffxml, dogaff, dogaff2, opt.prmtop, opt.inpcrd)

