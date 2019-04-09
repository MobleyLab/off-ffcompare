<<<<<<< HEAD
#Importing
=======
#importing
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
from rdkit import Chem, RDConfig, Geometry
from rdkit.Chem import AllChem, TorsionFingerprints
import os, sys
import openeye.oechem as oechem
from openeye import oeshape
import pandas as pd
import numpy as np
import itertools


<<<<<<< HEAD
def tanimotocombo(ref_mol, query_mol):
    """This is the TanimotoCombo function. It takes in two OEMols. 
    It does not matter (for our purposes) which mol is the reference mol and 
    which one is the query mol. The TanimotoCombo distance between the two
    mols is the same no matter which is the reference and which is the query.
    The GetTanimotoCombo() function comes from openeye. 
    """
=======
# In[2]:


# This is the tanimotocombo script. It takes in two OEMOLs.
# It does not matter which mol is the reference mol and which one is the query mol. 
# The TanimotoCombo distance between the two mols is the same no matter the order. 
#tanimotocombo function comes from openeye.
def tanimotocombo(ref_mol, query_mol):
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
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
    return("%s" % res.GetTanimotoCombo())


<<<<<<< HEAD
def TFDr(ref_mol, query_mol):
    """
    This is the TFDr script. 
    It makes use of RDKit's TFD calculation and the function rdmol_from_oemol.
    TFDr takes in two OEMOLs. 
    It does not matter which mol is the ref mol and which is the querymol. 
    TFD metric is the same no matter which is the ref and which is the query.
=======
# In[3]:


#This is the TFDr script. 
#It makes use of RDKit's TFD calculation and the function rdmol_from_oemol, defined below
#TFDr takes in two OEMOLs
#It does not matter which mol is the ref mol and which is the querymol
#TFD distance is the same no matter the order
def TFDr(ref_mol, query_mol):
    """
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
    First, OEmols are made RDKit compatible. Then,
    TFD is computed and returned using RDKit's TorsionFingerprints
    Module. Takes one input reference mol2 and one input query mol2.
    """
<<<<<<< HEAD
    #converts refmol to one readable by RDKit
    rrdmol2 = rdmol_from_oemol(ref_mol)
    #converts querymol to one readable by RDKit
    qrdmol2 = rdmol_from_oemol(query_mol)
    #If there was a mistake in the conversion process, return -1
    if (Chem.MolToSmiles(qrdmol2) != Chem.MolToSmiles(rrdmol2)):
        tfd = -1
    else:
        #calculates the TFD
=======
    #actually convert
    rrdmol2 = rdmol_from_oemol(ref_mol)
    #actually convert
    qrdmol2 = rdmol_from_oemol(query_mol)
    if (Chem.MolToSmiles(qrdmol2) != Chem.MolToSmiles(rrdmol2)):
        tfd = -1
    else:
        #get TFD
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
        tfd = TorsionFingerprints.GetTFDBetweenMolecules(rrdmol2, qrdmol2)
    return tfd


<<<<<<< HEAD

def rdmol_from_oemol(oemol):
    """
    Converts one oemol to a format recognizable by RDKit. 
    Written by Caitlin Bannan.
=======
# In[4]:


#This turns oemol into rdmol-compatible ones, for use with RDKit's TFD script. 
def rdmol_from_oemol(oemol):
    """
    Creates an openeye molecule object that is identical to the input rdkit molecule
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
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


<<<<<<< HEAD
def make_molname_df(directory,ffdirectorylist):
    """This is the make_molname_df function. It finds the names of mol2 files
    that are present in all forcefield-minimized directories in the directory.
    Its output is a pandas dataframe. 
    """
    dname = {}
    # for every forcefield, create dataframe of molnames
    for ffdirect in ffdirectorylist:
=======
# In[5]:


# This finds the molnames that are present in all three forcefield-minimized directories
# In the future, I'll change the first line to make it able to take in different force fields. 
def make_molname_df(directory,ffdirectorylist):
    allnamedf = pd.DataFrame(columns=['allMolName'])
    dname = {}
    for ffdirect in ffdirectorylist :
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
        tempdirectory = '%s/%s/' % (directory, ffdirect)
        list_of_mol2 = sorted(os.listdir(tempdirectory))
        list_of_molnames = list()
        for mol2 in list_of_mol2:
            molName = mol2.split('.')[0]
            list_of_molnames.append(molName)
        tempnamedf = pd.DataFrame(np.array(list_of_molnames))
        dname["%s" % ffdirect] = tempnamedf
<<<<<<< HEAD
    # merges all dataframes generated above into one dataframe
    # only keeps molecules present in all dataframes
=======
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
    twonamedf = dname[ffdirectorylist[0]].merge(dname[ffdirectorylist[1]])
    for i in range(2, len(ffdirectorylist)):
        threenamedf = dname[ffdirectorylist[i]].merge(twonamedf)
    threenamedf = threenamedf.rename({0: "MolNames"},axis='columns')
    return threenamedf


<<<<<<< HEAD
def all_info_df(ffdirectorylist,threenamedf):
    """
    This is the all_info_df function. It takes in the list of forcefields,
    as well as the dataframe of all molecule names, and runs TFD and Tanimoto
    Combo on all molecules. Its output is a dataframe of all this data. 
    """
    #Creating empty dictionaries that TFD and TANI scores will go in later,
    #as well as a heavyatomlist for putting heavy atoms in
    heavyatomlist = list()
    TFDdict = {}
    TANIdict = {}
    #creates combinations of forcefields and puts them into dictionaries
    for i,j in list(itertools.combinations(ffdirectorylist,2)):
        TFDdict['%s %s' % (i, j)] = {}
        TANIdict['%s %s' % (i, j)] = {}
    #generates all the data
=======
# In[6]:


#This makes dictionaries of TFD and TANI values for each combination of forcefields 
def all_info_df(ffdirectorylist,threenamedf):
    heavyatomlist = list()
    TFDdict = {}
    TANIdict = {}
    for i,j in list(itertools.combinations(ffdirectorylist,2)):
        TFDdict['%s %s' % (i, j)] = {}
        TANIdict['%s %s' % (i, j)] = {}

>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
    for molname in threenamedf['MolNames']:
        print(molname)
        mol_file = '%s' % molname + '.mol2'
        refmolin = oechem.oemolistream('%s/%s/%s' % (directory,ffdirectorylist[0],mol_file))
        refmolhev = oechem.OEGraphMol()
        oechem.OEReadMolecule(refmolin,refmolhev)
<<<<<<< HEAD
        #tries to get heavyatom count -- if it can't, returns -1. 
=======
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
        try:
            heavyvalue = oechem.OECount(refmolhev, oechem.OEIsHeavy())
        except (Exception, ValueError):
            heavyvalue = -1
        heavyatomlist.append(heavyvalue)
        refmolin.close()
<<<<<<< HEAD
        #Gets TanimotoCombo and TFD values
=======
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
        for i,j in list(itertools.combinations(ffdirectorylist,2)):
            refmolin = oechem.oemolistream('%s/%s/%s' % (directory,i, mol_file))
            refmol = oechem.OEGraphMol()
            oechem.OEReadMolecule(refmolin,refmol)
            qmolin = oechem.oemolistream('%s/%s/%s' % (directory,j, mol_file))
            qmol = oechem.OEGraphMol()
            oechem.OEReadMolecule(qmolin,qmol)
<<<<<<< HEAD
            #Tries to set TFD value. If it cannot, sets as -1. 
=======
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
            try:
                TFDvalue = TFDr(refmol,qmol)
            except (Exception, ValueError):
                TFDvalue = -1
            TFDdict['%s %s' % (i, j)]['%s' % molname]=TFDvalue
<<<<<<< HEAD
            #Tries to set TanimotoCombo value. If it cannot, sets as -1. 
=======
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
            try:
                TANIvalue = tanimotocombo(refmol,qmol)
            except (Exception, ValueError):
                TANIvalue = -1
            TANIdict['%s %s' % (i, j)]['%s' % molname]=TANIvalue
            qmolin.close()
            refmolin.close()
<<<<<<< HEAD
    #Loads data into dataframe
=======
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
    for key in TFDdict:
        tempdf = pd.DataFrame.from_dict(TFDdict['%s' % key],'index')
        tempdf = tempdf.rename({0:'TFD %s' % key}, axis = 'columns')
        tempdf['MolNames'] = tempdf.index
        threenamedf = threenamedf.merge(tempdf, on='MolNames')
    for key in TANIdict:
        tempdf = pd.DataFrame.from_dict(TANIdict['%s' % key],'index')
        tempdf = tempdf.rename({0:'TANI %s' % key}, axis = 'columns')
        tempdf['MolNames'] = tempdf.index
        threenamedf = threenamedf.merge(tempdf, on='MolNames')
    tempdf = pd.DataFrame(np.array(heavyatomlist))
    tempdf = tempdf.rename({0:'HeavyAtomCount'}, axis = 'columns')
    threenamedf['HeavyAtomCount'] = tempdf['HeavyAtomCount']
    return threenamedf

<<<<<<< HEAD
#Performs above functions.
=======

# In[8]:


>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
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
        print("ERROR: No working directory provided.")
    

    fflist = opt.ff.split(',')
    directory = opt.directory
    molnamedf = make_molname_df(directory,fflist)
    endmoldf = all_info_df(fflist,molnamedf)
<<<<<<< HEAD
    #Exports all data as csv 
=======
>>>>>>> d46fa0678973c22ff2e2d4155e82beff69ea24ef
    endmoldf.to_csv('%s/alldata.csv' % directory)

