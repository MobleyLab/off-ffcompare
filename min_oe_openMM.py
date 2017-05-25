#!/usr/bin/env python

"""

Written by Victoria Lim and Daisy Kyu with contributions from Caitlin Bannan
Mobley Lab UCI

### Description: This Python script loops over a directory of mol2 files, and
#      minimizes each with the specified force field (see supported FFs below).
#      In the directory where this script is called, a subdirectory is created
#      for the force field (if it doesn't already exist), and output mol2 files
#      are sorted into their respective subdirectories.
#   Minimizations completed using OpenEye (MMFF94*), and OpenMM (GAFF*, SMIRNOFF).
#   For N input *.mol2 files, there should be N output *.mol2 files, unless a
#      molecule file was corrupted or could not otherwise be minimized.
    For minimizations of different force field types, use the SAME input mol2 files.
       E.g. both GAFF and MMFF94 minimizations using Tripos mol2 input files.
       This will retain the same atom types / file format, with only coordinates
       being updated.

### Example usage:
# - python min_oe_openMM.py --inmols [directory with mol2] \
        --ffxml [SMIRNOFF ffxml] --gaffdir [directory with GAFF input files] \
        --gaff2dir [directory with GAFF2 input files] \
        --dommff [True or False] --log [log file name] \
        --outdir [directory to store output]

### Dependencies:
# - For all:    location of mol2 files. (-f)
# - For GAFFx:  needs directory with both *.prmtop and *.inpcrd files (-g, -G)
# - For SMIRNOFF: needs smirnoff*.ffxml file. (-x)
# - For MMFF94(S): set dommff=True
# - Optional:   specify output directory and name of a log file

### Supported force fields:
# - MMFF94 and MMFF94S (setting dommff=True will return both)
# - GAFF
# - GAFF2
# - SMIRNOFF formated force fields

### To do / Ideas:
"""

import os, glob
import numpy as np

import openeye.oechem as oechem
import openeye.oeszybki as oeszybki

import parmed
from parmed import unit as u

import simtk.openmm as mm
from simtk.openmm import app

from openforcefield.typing.engines.smirnoff import forcefield
from openforcefield.typing.engines.smirnoff import forcefield_utils as ff_utils


# -------------------------- Functions ------------------------- #

def writeUpdatedMol(Mol, fname, log):
    """

    Parameters
    ----------
    Mol: an OEChem molecule
    fname: str - name of the output mol2 file
    log: open file to write output to

    Returns
    -------
    True if the function completed successfully

    """

    # Open output file to write molecule.
    ofs = oechem.oemolostream()
    if os.path.exists(fname):
        log.write("Output .mol2 file already exists. Skipping.\n")
        return False
    if not ofs.open(fname):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % fname)

    # write out mol2 file and close the output filestream.
    oechem.OEWriteConstMolecule(ofs, Mol)
    ofs.close()
    return True

def optMMFF(Mol, FF, fname, log):
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

    if os.path.exists(fname):
        log.write("Output .mol2 file already exists. Skipping.\n")
        return False

    # set general energy options along with the run type specification
    optSzybki = oeszybki.OESzybkiOptions()
    optSzybki.SetSolventModel(oeszybki.OESolventModel_NoSolv)
    optSzybki.SetOptimizerType(oeszybki.OEOptType_BFGS)

    # set the particular force field
    if FF == "MMFF94":
        optSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94)
    elif FF == "MMFF94S":
        optSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94S)
    else:
        log.write( 'optMMFF failed for %s\n' %  tmpmol.GetTitle() )
        return False

    # create additional dependencies, then perform opt (in if statement)
    szOpt = oeszybki.OESzybki( optSzybki)  # generate minimization engine
    szResults = oeszybki.OESzybkiResults() # make object to hold szybki results
    if not szOpt(tmpmol, szResults):
        log.write( 'optMMFF failed for %s\n' %  tmpmol.GetTitle() )
        return False

    # try writing output file and return result
    return writeUpdatedMol(tmpmol, fname, log)





def optSMIRNOFF(input_Mol, FF_file, output_mol2, log ):
    """

    Creates OpenMM Topology, System, and initial positions of given molecule.

    Parameters
    ----------
    Mol: an OEChem molecule
    FF_file: string name of *.ffxml file with path
    output_mol2: string path/filename.mol2 to save output structure
    log: open file to write output data to

    Returns
    -------
    Boolean: True if output is successfully created

    """
    # make copy of the input_Mol
    Mol = oechem.OEMol(input_Mol)

    # check that the *.ffxml file exists
    if not os.path.exists(FF_file):
        log.write("Cannot find the ffxml file (%s)!\n"\
                % FF_file)
        return

    ff = forcefield.ForceField(FF_file)
    top, syst, pos  = ff_utils.create_system_from_molecule(ff, Mol, verbose = False)

    ### minimize with OpenMM and return the final coordinates in the OEform
    min_pos = minimizeOpenMM(top, syst, pos)
    Mol.SetCoords(oechem.OEFloatArray(min_pos))

    return writeUpdatedMol(Mol, output_mol2, log)


def optGAFFx(mol, gaffdir, output_mol2, log):
    """
    Uses OEMol and GAFF input files to minimize the system and write
    output mol2

    Parameters
    ----------
    mol: OEMol
    gaffdir: directory with GAFF input files
    output_mol2: string, path/to/mol2 where output structure should be saved
    log: open file to write log data

    Returns
    -------
    boolean: True minimization succeeded, False otherwise

    """
    ### Check: optimized file not existing, gaff dependencies present, gaff files have content
    if os.path.exists(output_mol2):
        log.write('Optimization file %s already exists\n' % (output_mol2))
        return False

    ### Locate GAFF input files
    prmFile = os.path.join(gaffdir,mol.GetTitle()+'.prmtop')
    inpFile = os.path.join(gaffdir,mol.GetTitle()+'.inpcrd')

    # Check prmtop file
    if not os.path.exists(prmFile):
        log.write("%s does not exist, skipping minimization\n" % prmFile)
        return False
    elif not os.path.getsize(prmFile) > 0:
        log.write("%s is empty, skipping minimization\n" % prmFile)
        return False

    # Check inpcrd file
    if not os.path.exists(inpFile):
        log.write("%s does not exist, skipping minimization\n" % inpFile)
        return False
    elif not os.path.getsize(inpFile):
        log.write("%s is empty, skipping minimization\n" % inpFile)
        return False

    # make copy of mol2 file
    tmpmol = oechem.OEMol(mol)

    # load input files and create parmed system
    parm = parmed.load_file(prmFile, inpFile)
    Topology = parm.topology
    Positions = parm.positions
    System = parm.createSystem(nonbondedMethod=app.NoCutoff)

    ### Use parmed to write out mol2 file from optimized coordinates.
    minimized_positions = minimizeOpenMM(Topology, System, Positions)
    tmpmol.SetCoords(oechem.OEFloatArray(minimized_positions))

    return writeUpdatedMol(tmpmol, output_mol2, log)



def minimizeOpenMM(Topology, System, Positions):
    """

    Minimize molecule with specified topology, system, and positions
       using OpenMM. Return the positions of the optimized moelcule.

    Parameters
    ----------
    Topology:  OpenMM topology
    System:    OpenMM system
    Positions: OpenMM positions

    Returns
    -------
    concat_coords: list of positions ready to add to an OEMol

    """

    # need to create integrator but don't think it's used
    integrator = mm.LangevinIntegrator(
            300.0 * u.kelvin,
            1.0 / u.picosecond,
            2.0 * u.femtosecond)

    simulation = app.Simulation(Topology, System, integrator)
    simulation.context.setPositions(Positions)
    simulation.minimizeEnergy(tolerance=5.0E-9, maxIterations=1500)

    positions =  simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    positions = positions/u.angstroms
    coordlist = list()
    for atom_coords in positions:
        coordlist += [i for i in atom_coords]
    return coordlist


# --------------------------- Main Function ------------------------- #

def load_and_minimize(infiles, log, output_dir,
        dommff = False, ffxml = None, gaff = None, gaff2 = None):
    """

    This function finds all mol2 file in the input directory,
    loads the molecules and minimizes them
    (if it does not already exist) based on the specified fftype.

    Parameters
    ----------
    infiles: Path to directory containing all mol2 files
    log: an open file to write output data
    output_dir: string, directory for output files to be written
    dommff: Boolean - True to minimize with MMFF94 or MMFF94S
    ffxml: String - path and name of *.ffxml file for SMIRNOFF minimization
    gaff: String - path to directory containing GAFF Inpcrd and Prmtop files
    gaff2: String - path to directory containing GAFF2 Inpcrd and Prmtop files

    """
    molfiles = glob.glob(os.path.join(infiles, '*.mol2'))

    ### Loop over mols.
    for f in molfiles:
        mol_title = f.split('.')[0]
        ### Load molecule from .mol2 file.
        ifs = oechem.oemolistream()
        if not ifs.open(f):
             oechem.OEThrow.Warning("Unable to open %s for reading" % f)
        try:
            mol = next(ifs.GetOEMols())
        except StopIteration: # No molecule in this mol2, move on
            log.write('No mol loaded for %s (StopIteration exception)\n'\
                    % mol_title)
            continue
        ifs.close()

        # Output name is always [molTitle].mol2
        fname = mol.GetTitle()+'.mol2'

        if dommff:  # Optimize with MMFF94 and MMFF94S. Produce output mol2 files for each.

            # create subdirectory for MMFF94 and MMFF94S jobs
            if not os.path.exists(os.path.join(output_dir, 'MMFF94')):
                os.makedirs(os.path.join(output_dir, 'MMFF94'))
            if not os.path.exists(os.path.join(output_dir, 'MMFF94S')):
                os.makedirs(os.path.join(output_dir, 'MMFF94S'))

            ### PART ONE: MMFF94
            print('Starting on MMFF94 optimization for %s' % mol.GetTitle())

            mmff94_file = os.path.join(output_dir+'/MMFF94', fname)
            if not optMMFF(mol, 'MMFF94', mmff94_file, log):
                log.write('MMFF94 minimization failed for molecule %s:\n'\
                        %  (mol.GetTitle()) )

            ### PART TWO: MMFF94S
            print('Starting on MMFF94S optimization for %s'\
                    % mol.GetTitle())
            mmff94s_file = os.path.join(output_dir+'/MMFF94S', fname)
            if not optMMFF(mol, 'MMFF94S', mmff94s_file, log):
                log.write('MMFF94S minimization failed for molecule %s:\n'\
                        %  (mol.GetTitle()) )

        # Optimize with smirnoff ffxml. Produce output mol2 files.
        if ffxml is not None:

            ### Set output file name and make a copy of molecule for opt.
            print('Starting on SMIRNOFF optimization for %s\n' % mol.GetTitle())
            if not os.path.exists(os.path.join(output_dir, 'SMIRNOFF')):
                os.makedirs(os.path.join(output_dir, 'SMIRNOFF'))

            fulln = os.path.join(output_dir+'/SMIRNOFF', fname)

            ### Check that optimized file doesn't already exist.
            if not os.path.exists(fulln):
                optSMIRNOFF(mol, ffxml, fulln, log)
            else:
                log.write('SMIRNOFF minimized molecule file for %s already \
                        exists\n' % mol.GetTitle())


        # Optimize with GAFF. Produce output mol2 files.
        if gaff is not None:

            ### Set output file name and make a copy of molecule for opt.
            if not os.path.exists(os.path.join(output_dir, 'GAFF')):
                os.makedirs(os.path.join(output_dir, 'GAFF'))

            fulln = os.path.join(output_dir+'/GAFF', fname)
            print('Starting on GAFF optimization for %s\n'\
                    % mol.GetTitle())
            optGAFFx(mol, gaff, fulln, log)

        # Optimize with GAFF2. Produce output mol2 files.
        if gaff2 is not None:
            ### Set output file name and make a copy of molecule for opt.
            if not os.path.exists(os.path.join(output_dir, 'GAFF2')):
                os.makedirs(os.path.join(output_dir, 'GAFF2'))

            fulln = os.path.join(output_dir+'/'+'GAFF2', fname)
            print('Starting on GAFF2 optimization for %s\n'\
                    % mol.GetTitle())
            optGAFFx(mol, gaff2, fulln, log)


# ------------------------- Parse Inputs ----------------------- #



if __name__ == '__main__':
    from optparse import OptionParser
    usage_string="""
    Minimize molecules using a variety of forcefields.
    The molecules to be minimized are specified in an input
    "inmols" directory which is expected to contain Tripos mol2 files
    for all molecules you with to minimize.

    usage:
    python min_oe_openMM.py --inmols mol2Directory \
        [--ffxml SMIRNOFF.ffxml] --gaffdir gaffDirectory \
        --gaff2dir gaff2Directory --dommff True --log output.dat \
        --outdir outputDirectory]
    """

    parser = OptionParser(usage = usage_string)

    parser.add_option('-i', '--inmols',
            default = None,
            help = "REQUIRED! Path to directory containing all Tripos .mol2 \
files to be minimized.",
            type = "string",
            dest = 'inmols')

    parser.add_option('-x', '--ffxml',
            default = None,
            help = "OPTIONAL! Force field file for the SMIRNOFF being used.\
SMIRNOFF minimization is only performed when this is specified.",
            type = "string",
            dest = 'ffxml')

    parser.add_option('-g', '--gaffdir',
            help = "OPTIONAL! Directory containing GAFF parameter/topology \
 files (*.prmtop) and coordinates files (*.inpcrd). GAFF minimization only \
 performed if this directory is specified. ",
            type = "string",
            default = None,
            dest = 'gaffdir')

    parser.add_option('-G', '--gaff2dir',
            help = "OPTIONAL! Directory containing GAFF parameter/topology \
 files (*.prmtop) and coordinates files (*.inpcrd). GAFF minimization only \
 performed if this directory is specified. ",
            type = "string",
            default = None,
            dest = 'gaff2dir')

    parser.add_option('-f', '--dommff',
            help = "OPTIONAL! If True, minimizations will be performed \
with the MMFF94 and MMFF94S forcefields. Results will be in the form of \
mol2 files in the respective directory.",
            type = 'choice', choices = ['True', 'False'],
            default = 'False',
            dest = 'dommff')

    parser.add_option('-o', '--outdir',
            help = "OPTIONAL! location for output mol2 files. In this \
directory, sub-directories will be created for each force field type. \
The default is the current working directory.",
            type = 'string',
            default = None,
            dest = 'outdir')

    parser.add_option('-l', '--log',
            help = "OPTIONAL! Name for output file, this will be place \
in the output directory. Default is output.dat.",
            type = "string",
            default = "output.dat",
            dest = 'log')

    (opt, args) = parser.parse_args()
    ### Check required fields.
    if opt.inmols is None:
        parser.error("ERROR: No input directory was specified.")
    ### Check inmols exists!
    if not os.path.isdir(opt.inmols):
        parser.error("ERROR: input directory (%s) was not found" % opt.inmols)

    dommff = opt.dommff == 'True'

    ### Check that at least one force field input was provided
    if (opt.gaffdir is None) and (opt.gaff2dir is None) and (opt.ffxml is None) and (not dommff):
        parser.error("ERROR: You must specify information for at least \
one force field with the option --gaff, --gaff2, --dommff, --ffxml")

    ### Check for output directory
    if opt.outdir is None:
        outdir = os.getcwd()
    else:
        if not os.path.isdir(opt.outdir):
            parser.error("ERROR: output directory %s does not exist" % opt.outdir)
        outdir = opt.outdir

    ### Open log file
    log = open(os.path.join(outdir, opt.log), 'w')

    ### Check if input directories exist
    if (opt.gaffdir is not None) and (not os.path.isdir(opt.gaffdir)):
        parser.error("ERROR: GAFF directory %s does not exist" % opt.gaffdir)
    if (opt.gaff2dir is not None) and (not os.path.isdir(opt.gaff2dir)):
        parser.error("ERROR: GAFF directory %s does not exist" % opt.gaff2dir)
    if (opt.ffxml is not None) and (not os.path.isfile(opt.ffxml)):
        parser.error("ERROR: FFXML file %s does not exist" % opt.ffxml)

    load_and_minimize(opt.inmols, log, outdir, dommff, opt.ffxml, opt.gaffdir, opt.gaff2dir)
    log.close()
