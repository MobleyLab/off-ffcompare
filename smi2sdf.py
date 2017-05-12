#!/usr/bin/env python

## This script takes a list of smiles as "filename.smi" and for each smiles
##   molecule, 1 omega conformer is generated. All mols are written into 1 sdf
##   file. Individual molecule sdf files can also be written.

## Edit wdir and smiles in 'Variables' section, then "python file.py"
## OR, import and call smi2sdf.smi2sdf(wdir, smiles)

import os, sys
import openeye.oechem as oechem
import openeye.oeomega as oeomega
import openeye.oeszybki as oeszybki

def GenerateConfs(Mol):
    """
    Generate conformers of molecule from its SMILES string.

    Parameters
    ----------
    Mol:          OEChem molecule

    Returns
    -------
    molWithConfs: OEChem molecule with omega-generated conformers

    """
    molWithConfs = oechem.OEMol(Mol)
    omega = oeomega.OEOmega()
    maxConfs = 1
    omega.SetMaxConfs(maxConfs)
    omega.SetStrictStereo(False)
    omega.SetEnumNitrogen( oeomega.OENitrogenEnumeration_All)
    if not omega(molWithConfs):
        print( "omega failed on %s" % mol.GetTitle() )
        return None
    else:
        return molWithConfs


def smi2sdf(sdfout, smiles):
    """
    From a file containing smiles strings, generate omega conformers,
       resolve steric clashes, do a quick MM opt, and write SDF output.

    Parameters
    ----------
    sdfout: str - output sdf file. E.g. "name.sdf"
    smiles: str - name of the smiles file. E.g. "name.smi"

    """
    ### Read in smiles file.
    ifs = oechem.oemolistream()
    if not ifs.open(smiles):
        oechem.OEThrow.Warning("Unable to open %s for reading" % smiles)

    ### Open output file to write molecules.
    ofs = oechem.oemolostream()
    if os.path.exists(sdfout):
        #sys.exit("Output .sdf file already exists. Exiting.\n")
        print("Output .sdf file already exists. Exiting.\n")
        return
    if not ofs.open(sdfout):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % sdfout)

    ### For each molecule: label atoms, generate 1 conf, write output.
    for smimol in ifs.GetOEMols():
        oechem.OETriposAtomNames(smimol)
        mol = GenerateConfs(smimol)
        oechem.OEWriteConstMolecule(ofs, mol)

    ### Close files.
    ifs.close()
    ofs.close()


def smi2indivSdf(wdir, smiles):
    """
    From a file containing smiles strings, generate omega conformers,
       resolve steric clashes, do a quick MM opt, and write SDF output
       as individual files.

    Parameters
    ----------
    wdir: str - working directory containing .smi file
    smiles: str - name of the smiles file. E.g. "name.smi"

    """
    os.chdir(wdir)

    ### Read in smiles file.
    ifs = oechem.oemolistream()
    if not ifs.open(smiles):
        oechem.OEThrow.Warning("Unable to open %s for reading" % smiles)


    ### For each molecule: label atoms, generate 1 conf, write output.
    for i, smimol in enumerate(ifs.GetOEMols()):
        oechem.OETriposAtomNames(smimol)
        mol = GenerateConfs(smimol)

        ### Open output file to write molecules.
        sdfout = smiles.split('.')[0] + '-%d.sdf' %(i)
        ofs = oechem.oemolostream()
        if os.path.exists(sdfout):
            #sys.exit("Output .sdf file already exists. Exiting.\n")
            print("Output .sdf file already exists. Exiting.\n")
            return
        if not ofs.open(sdfout):
            oechem.OEThrow.Fatal("Unable to open %s for writing" % sdfout)
        oechem.OEWriteConstMolecule(ofs, mol)
        ofs.close()

    ### Close files.
    ifs.close()

# ----------------- Run from commandline  ----------------------------------

if __name__ == '__main__':
    from optparse import OptionParser

    usage_string = """
    This script is used to convert a list of SMILES string
    to an sdf file with all molecules.
    smi files have the format:
    [SMILES] [Name]

    usage:
    python smi2sdf.py --smiles smilesFile.smi --wdir working/directory/
    """
    parser = OptionParser(usage = usage_string)
    parser.add_option('-s', '--smiles',
            help = "REQUIRED! SMILES file",
            type = "string",
            dest = "smiles")

    parser.add_option('-o', '--sdf',
            help = "REQUIRED! output sdf file",
            type = "string",
            dest = "sdf")

    (opt, args) = parser.parse_args()

    if opt.smiles is None:
        parser.error("Must provide smiles file")
    if opt.sdf is None:
        parser.error("Must provide output SDF filename")

    smi2sdf(opt.sdf, opt.smiles)
