# This is the dataextractor script. 
# It was made for use in the Open ForceField forcefield-compare project.
# It was written primarily by Jordan Ehrman and Caitlin Bannan, with 
# help from Victoria Lim. 
# The two user-inputs this script needs is a list of forcefields, 
# separated by commas, and the path to a directory. 
# This directory must have a specific format: It must contain at least
# two subdirectories with the names of the specified force fields. Each 
# of these subdirectories should then include the .mol2 files that should 
# be compared by TFD and TanimotoCombo. 
# Last updated April 15, 2019. 

# Importing

import os, sys
import pandas as pd
import numpy as np
import itertools
import oeneye.oechem as oechem
from openeye import oeshape
from rdkit import Chem, RDConfig, Geometry
from rdkit.Chem import AllChem, TorsionFingerprints

def tanimotocombo(ref_mol, query_mol):
    """
    This is the TanimotoCombo function. It takes in two OEMols. 
    It does not matter (for our purposes) which mol is the reference mol and 
    which one is the query mol. The TanimotoCombo distance between the two
    mols is the same no matter which is the reference and which is the query.
    The GetTanimotoCombo() function comes from openeye.
     
    Args:
        ref_mol (oemol) An oemol that has already been read in. 
        query_mol (oemol) Another oemol that has already been read in. 
 
    Returns: 
        res.GetTanimotoCombo() (float) The TanimotoCombo value.
    """
    # Prepare reference molecule for calculation
    # With default options this will remove any explicit
    # hydrogens present and add color atoms
    prep = oeshape.OEOverlapPrep()
    prep.Prep(ref_mol)
    # Get appropriate function to calculate both shape and color
    # By default the OEOverlapFunc contains OEGridShapeFunc for shape
    # and OEExactColorFunc for color
    func = oeshape.OEOverlapFunc()
    func.SetupRef(ref_mol)
    res = oeshape.OEOverlapResults()
    prep.Prep(query_mol)
    func.Overlap(query_mol, res)
    return(res.GetTanimotoCombo())


def TFD_for_oemols(ref_mol, query_mol):
    """
    This is the TFD_for_oemols script. 
    It makes use of RDKit's TFD calculation and the function rdmol_from_oemol.
    TFD_for_oemols takes in two OEMOLs. 
    It does not matter which mol is the ref mol and which is the querymol. 
    TFD metric is the same no matter which is the ref and which is the query.
    First, OEmols are made RDKit compatible. Then,
    TFD is computed and returned using RDKit's TorsionFingerprints
    Module. Takes one input reference mol2 and one input query mol2.
    
    Args: 
        ref_mol (oemol) An oemol that has already been read in. 
        query_mol (oemol) An oemol that has already been read in. 
    
    Returns: 
        tfd (float) The torsion fingerprint deviation between ref and query.
    """
    # converts refmol to one readable by RDKit
    rrdmol2 = rdmol_from_oemol(ref_mol)
    # converts querymol to one readable by RDKit
    qrdmol2 = rdmol_from_oemol(query_mol)
    # If there was a mistake in the conversion process, return -1
    if (Chem.MolToSmiles(qrdmol2) != Chem.MolToSmiles(rrdmol2)):
        tfd = -1
    else:
        # calculates the TFD
        tfd = TorsionFingerprints.GetTFDBetweenMolecules(rrdmol2, qrdmol2)
    return tfd



def rdmol_from_oemol(oemol):
    """
    Converts one oemol to a format recognizable by RDKit. 
    Written by Caitlin Bannan.

    Args: 
        oemol (oemol) an oemol that has already been read in. 

    Returns: 
        rdmol.getmol() (rdmol) The same molecule, now made RDKit compatible. 
    """
    print("Starting molecule")

    # start function
    rdmol = Chem.RWMol()

    # RDKit keeps bond order as a type instead using these values, I don't really understand 7,
    # I took them from Shuzhe's example linked above
    _bondtypes = {
        1: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        4: Chem.BondType.QUADRUPLE,
        5: Chem.BondType.QUINTUPLE,
        6: Chem.BondType.HEXTUPLE,
        7: Chem.BondType.ONEANDAHALF, }

    # atom map lets you find atoms again
    map_atoms = dict()  # {oe_idx: rd_idx}
    for oea in oemol.GetAtoms():
        oe_idx = oea.GetIdx()
        rda = Chem.Atom(oea.GetAtomicNum())
        rda.SetFormalCharge(oea.GetFormalCharge())
        rda.SetIsAromatic(oea.IsAromatic())

        # unlike OE, RDK lets you set chirality directly
        cip = oechem.OEPerceiveCIPStereo(oemol, oea)
        if cip == oechem.OECIPAtomStereo_S:
            rda.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
        if cip == oechem.OECIPAtomStereo_R:
            rda.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)

        map_atoms[oe_idx] = rdmol.AddAtom(rda)

    # As discussed above, setting bond stereochemistry requires neighboring bonds
    # so we will store that information by atom index in this list
    stereo_bonds = list()
    # stereo_bonds will have tuples with the form (rda1, rda2, rda3, rda4, is_cis)
    # where rda[n] is an atom index for a double bond of form 1-2=3-4
    # and is_cis is a Boolean is True then onds 1-2 and 3-4 are cis to each other

    aro_bond = 0
    for oeb in oemol.GetBonds():
        # get neighboring rd atoms
        rd_a1 = map_atoms[oeb.GetBgnIdx()]
        rd_a2 = map_atoms[oeb.GetEndIdx()]

        # AddBond returns the total number of bonds, so addbond and then get it 
        rdmol.AddBond(rd_a1, rd_a2)
        rdbond = rdmol.GetBondBetweenAtoms(rd_a1, rd_a2)

        # Assign bond type, which is based on order unless it is aromatic
        order = oeb.GetOrder()
        if oeb.IsAromatic():
            rdbond.SetBondType(_bondtypes[1.5])
            rdbond.SetIsAromatic(True)
        else:
            rdbond.SetBondType(_bondtypes[order])
            rdbond.SetIsAromatic(False)

        # If the bond has specified stereo add the required information to stereo_bonds
        if oeb.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
            # OpenEye determined stereo based on neighboring atoms so get two outside atoms
            n1 = [n for n in oeb.GetBgn().GetAtoms() if n != oeb.GetEnd()][0]
            n2 = [n for n in oeb.GetEnd().GetAtoms() if n != oeb.GetBgn()][0]

            rd_n1 = map_atoms[n1.GetIdx()]
            rd_n2 = map_atoms[n2.GetIdx()]

            stereo = oeb.GetStereo([n1, n2], oechem.OEBondStereo_CisTrans)
            if stereo == oechem.OEBondStereo_Cis:
                print('cis')
                stereo_bonds.append((rd_n1, rd_a1, rd_a2, rd_n2, True))
            elif stereo == oechem.OEBondStereo_Trans:
                print('trans')
                stereo_bonds.append((rd_n1, rd_a1, rd_a2, rd_n2, False))

                # add bond stereochemistry:
    for (rda1, rda2, rda3, rda4, is_cis) in stereo_bonds:
        # get neighbor bonds
        bond1 = rdmol.GetBondBetweenAtoms(rda1, rda2)
        bond2 = rdmol.GetBondBetweenAtoms(rda3, rda4)

        # Since this is relative, the first bond always goes up
        # as explained above these names come from SMILES slashes so UP/UP is Trans and Up/Down is cis
        bond1.SetBondDir(Chem.BondDir.ENDUPRIGHT)
        if is_cis:
            bond2.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
        else:
            bond2.SetBondDir(Chem.BondDir.ENDUPRIGHT)

    # if oemol has coordinates (The dimension is non-zero)
    # add those coordinates to the rdmol
    if oechem.OEGetDimensionFromCoords(oemol) > 0:
        conformer = Chem.Conformer()
        oecoords = oemol.GetCoords()
        for oe_idx, rd_idx in map_atoms.items():
            (x,y,z) = oecoords[oe_idx]
            conformer.SetAtomPosition(rd_idx, Geometry.Point3D(x,y,z))
        rdmol.AddConformer(conformer)

    # Save the molecule title
    rdmol.SetProp("_Name", oemol.GetTitle())

    # Cleanup the rdmol
    # Note I copied UpdatePropertyCache and GetSSSR from Shuzhe's code to convert oemol to rdmol here:
    rdmol.UpdatePropertyCache(strict=False)
    Chem.GetSSSR(rdmol)
    # I added AssignStereochemistry which takes the directions of the bond set
    # and assigns the stereochemistry tags on the double bonds
    Chem.AssignStereochemistry(rdmol, force=False)

    print("Final Molecule: ", Chem.MolToSmiles(Chem.RemoveHs(rdmol), isomericSmiles=True))
    return rdmol.GetMol()


def make_molname_df(directory,ffdirectorylist):
    """
    This is the make_molname_df function. It finds the names of mol2 files
    that are present in all forcefield-minimized directories in the directory.
    Its output is a pandas dataframe. 

    Args: 
        directory (str) path to directory containing FFsubdirectories
        ffdirectorylist (list) list of force fields

    Returns: 
        all_ff_df (dataframe) dataframe containing names of molecules
            present in all FF subdirectories
    """
    dname = {}
    # for every forcefield, create dataframe of molnames
    for ffdirect in ffdirectorylist:
        tempdirectory = '%s/%s/' % (directory, ffdirect)
        list_of_mol2 = sorted(os.listdir(tempdirectory))
        list_of_molnames = list()
        for mol2 in list_of_mol2:
            molName = mol2.split('.')[0]
            list_of_molnames.append(molName)
        tempnamedf = pd.DataFrame(np.array(list_of_molnames))
        dname["%s" % ffdirect] = tempnamedf
    # merges all dataframes generated above into one dataframe
    # only keeps molecules present in all dataframes
    twonamedf = dname[ffdirectorylist[0]].merge(dname[ffdirectorylist[1]])
    for i in range(2, len(ffdirectorylist)):
        all_ff_df = dname[ffdirectorylist[i]].merge(twonamedf)
    all_ff_df = all_ff_df.rename({0: "MolNames"},axis='columns')
    return all_ff_df


def all_info_df(ffdirectorylist, all_ff_df):
    """
    This is the all_info_df function. It takes in the list of forcefields,
    as well as the dataframe of all molecule names, and runs TFD and Tanimoto
    Combo on all molecules. Its output is a dataframe of all this data.

    Args: 
        ffdirectorylist (list) list of ff to compare 
        all_ff_df (dataframe) dataframe created by make_molname_df func above.

    Returns: 
        all_ff_df (dataframe) same dataframe with appended columns. 
    """
    # Creating empty dictionaries that TFD and TANI scores will go in later,
    # As well as a heavyatomlist for putting heavy atoms in
    heavyatomlist = list()
    TFDdict = {}
    TANIdict = {}
    # Creates combinations of forcefields and puts them into dictionaries
    for i,j in list(itertools.combinations(ffdirectorylist,2)):
        TFDdict['%s %s' % (i, j)] = {}
        TANIdict['%s %s' % (i, j)] = {}
    # Generates all the data
    for molname in all_ff_df['MolNames']:
        print(molname)
        mol_file = '%s' % molname + '.mol2'
        refmolin = oechem.oemolistream('%s/%s/%s' % (directory,ffdirectorylist[0],mol_file))
        refmolhev = oechem.OEGraphMol()
        oechem.OEReadMolecule(refmolin,refmolhev)
        heavyvalue = oechem.OECount(refmolhev, oechem.OEIsHeavy())
        heavyatomlist.append(heavyvalue)
        refmolin.close()
        # Gets TanimotoCombo and TFD values
        for i,j in list(itertools.combinations(ffdirectorylist,2)):
            refmolin = oechem.oemolistream('%s/%s/%s' % (directory,i, mol_file))
            refmol = oechem.OEGraphMol()
            oechem.OEReadMolecule(refmolin,refmol)
            qmolin = oechem.oemolistream('%s/%s/%s' % (directory,j, mol_file))
            qmol = oechem.OEGraphMol()
            oechem.OEReadMolecule(qmolin,qmol)
            # Getting TFD
            TFDvalue = TFD_for_oemols(refmol,qmol)
            TFDdict['%s %s' % (i, j)]['%s' % molname]=TFDvalue
            # Getting TanimotoCombo
            TANIvalue = tanimotocombo(refmol,qmol)
            TANIdict['%s %s' % (i, j)][molname]=TANIvalue
            qmolin.close()
            refmolin.close()
    # Loads data into dataframe
    for key in TFDdict:
        tempdf = pd.DataFrame.from_dict(TFDdict['%s' % key],'index')
        tempdf = tempdf.rename({0:'TFD %s' % key}, axis = 'columns')
        tempdf['MolNames'] = tempdf.index
        all_ff_df = all_ff_df.merge(tempdf, on='MolNames')
    for key in TANIdict:
        tempdf = pd.DataFrame.from_dict(TANIdict['%s' % key],'index')
        tempdf = tempdf.rename({0:'TANI %s' % key}, axis = 'columns')
        tempdf['MolNames'] = tempdf.index
        all_ff_df = all_ff_df.merge(tempdf, on='MolNames')
    tempdf = pd.DataFrame(np.array(heavyatomlist))
    tempdf = tempdf.rename({0:'HeavyAtomCount'}, axis = 'columns')
    all_ff_df['HeavyAtomCount'] = tempdf['HeavyAtomCount']
    return all_ff_df

# Performs above functions.
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-f','--ff',
            help = "REQUIRED: Name(s) of the force field directory(s) separated by commas. At least three required.",
            type = "string",
            dest = 'ff')

    parser.add_option('-d','--directory',
            help = "REQUIRED: Name of the full path directory which contains the force field directories listed in --ff. If None, then the current working directory is used",
            type = "string",
            dest = 'directory')
    (opt, args) = parser.parse_args()
    if opt.ff == None:
        parser.error("ERROR: No force fields were specified.")
    if opt.directory == None: 
        parser.error("ERROR: No working directory provided.")
    

    fflist = opt.ff.split(',')
    directory = opt.directory
    molnamedf = make_molname_df(directory,fflist)
    endmoldf = all_info_df(fflist,molnamedf)
    #Exports all data as csv 
    endmoldf.to_csv('%s/alldata.csv' % directory)

