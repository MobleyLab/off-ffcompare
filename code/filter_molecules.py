"""
This script was created to filter the DrugBank and eMolecule
databases to gaurentee the following for each molecule in our
final set

* No more than 200 heavy atoms
* No metals
* Gaurentee appropriate valency (no hypervalent carbon for example)

Authors
* Caitlin C. Bannan
    UC Irvine, Mobley Group, bannanc@uci.edu
"""

from openeye import oechem
import string
from openeye import oeomega
import smarty

def check_valence(mol):
    """
    Checks for hypervalency
    Parameter
    ---------
    mol - OEMol()

    Return
    ------
    boolean - True (no inappropriate valency)
              False (an atom with atomic number < 10 has > 4 Valence)
    """
    for atom in mol.GetAtoms():
        atomNum = atom.GetAtomicNum()
        # find number of neighbors to this atom
        valence = atom.GetValence()
        if atomNum <= 10: # first row elements
            if valence > 4:
                print("Found a #%i atom with valence %i in molecule %s" % (atomNum, valence, oechem.OECreateIsoSmiString(mol)))
                return False
    return True

def keep_molecule(mol, remove_smirks = list()):
    """
    Determines if the molecule will be stored.

    Parameters
    ----------
    mol - OEMol
    remove_smirks - list of SMIRKS strings you don't want in your molecules

    Returns
    -------
    boolean - True (molecule meets the requirements below)
            - has no metal atoms
            - no more than 200 heavy atoms
            - has none of the SMIRKS in remove_smirks list
            - molecule has appropriate valency
    """
    # Check number of metal atoms
    if oechem.OECount(mol, oechem.OEIsMetal()) > 0:
        return False
    # Check number of heavy atoms
    if oechem.OECount(mol, oechem.OEIsHeavy()) > 200:
        return False
    # Check for patterns in remove smirks list
    for smirks in remove_smirks:
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smirks):
            continue
        ss = oechem.OESubSearch(qmol)
        matches = [match for match in ss.Match(mol, False)]
        if len(matches) > 0:
            return False
    # check valency
    return check_valence(mol)

def parse_smile(smile_f):
    """
    Parses file with the form
    "molecule_name    SMILES"

    Parameter
    ---------
    smile_f: str, "path/to/smiles_file"

    Return
    ------
    list of smiles in file
    """
    f = open(smile_f,'r')
    lines = f.readlines()
    f.close()
    # get the SMILES string from each line in the file
    smiles = [l.split()[1] for l in lines]
    return smiles

def eMolecules_filtering(input_f, current_smiles = list()):
    """
    This function was used to filter eMolecules database
    and the eMolecules_incremental database.
    It creates all the filtered output files with 1000 molecules
    in each sdf file and 1,000,000 molecule-ID to smiles strings in each
    text file

    Parameter
    ---------
    input_f : string "path/to/inputfile.sdf"
    current_smiles : list of strings; smiles already in your molecule sets
    """
    set_name = input_f.split('.')[0]
    output_f = set_name+"_%i.sdf"
    smiles_base = set_name+"_%i.txt"
    molecule_name = set_name+"_%i_%i"

    # Load and check input file
    ifs = oechem.oemolistream(input_f)
    if not ifs.IsValid():
        raise Exception("Error: input_file (%s) was not valid" % input_f)

    errs = oechem.oeosstream()
    oechem.OEThrow.SetOutputStream(errs)

    molecule = oechem.OECreateOEGraphMol()
    count = 0
    smile_count = 0
    saved = 0
    switch = False

    # first output file
    current_letter = 1000
    ofs_file = output_f%current_letter
    ofs = oechem.oemolostream(ofs_file)
    if not ofs.IsValid():
        raise Exception("output file %s is not valid" % ofs_file)
    add_smiles = open(smiles_base % current_letter, 'a')

    while oechem.OEReadMolecule(ifs, molecule):
        # count input file molecules
        count +=1

        if switch: # If True create new output file
            switch = False
            ofs.close()
            current_letter += 1
            ofs_file = output_f % current_letter
            # Load and check output file
            ofs = oechem.oemolostream(ofs_file)
            if not ofs.IsValid():
                raise Exception("output file %s is not valid" % ofs_file)
            print("Switching to file %s, currently saved %i molecules" % (ofs_file, saved))
            if current_letter%100 == 0:
                add_smiles.close()
                add_smiles = open(smiles_base % current_letter, 'a')

        # IF smiles in current list skip the molecule
        smi = oechem.OECreateIsoSmiString(molecule)
        if smi in current_smiles:
            smile_count += 1
            continue

        # Make copy of molecule before making changes
        mol_copy = oechem.OEMol(molecule)
        oechem.OEAddExplicitHydrogens(mol_copy)
        # if the molecule meets our requirements save to current output
        if keep_molecule(mol_copy):
            mol_title = molecule_name % (current_letter,count)
            mol_copy.SetTitle(mol_title)
            add_smiles.writelines("%s\t\t%s\n" % (mol_title, smi))
            oechem.OEWriteMolecule(ofs, mol_copy)
            saved += 1
            if saved%1000 == 0:
                switch = True

    print("%i molecules in input file" % (count))
    print("%i molecules were had repeated isomeric SMILES" % smile_count)
    print("%i molecules saved to output files" % (saved))

    ifs.close()
    ofs.close()

if __name__=="__main__":
    # =========DrugBank=======================
    # DrugBank was done a bit differently so it could not
    # use the eMolecules_filtering method written above
    input_f = "DrugBank.sdf"
    output_f = "DrugBank_%s.sdf"
    molecule_name = "DrugBank_%s%i"

    # get letters to diferentiate output
    letters = string.ascii_letters
    letters = [l for l in letters]

    # get current smiles
    smiles_f = "smiles_to_ID_off-compare.txt"
    current_smiles = parse_smile(smiles_f)
    add_smiles = open(smiles_f, 'a')

    # Load and check input file
    ifs = oechem.oemolistream(input_f)
    ifs.SetFormat(oechem.OEFormat_SDF)
    if not ifs.IsValid():
        raise Exception("Error: input_file (%s) was not valid" % input_f)

    molecule = oechem.OECreateOEGraphMol()
    count = 0
    smile_count = 0
    saved = 0
    switch = False

    # first output file
    current_letter = letters.pop(0)
    ofs_file = output_f%current_letter
    ofs = oechem.oemolostream(ofs_file)
    ifs.SetFormat(oechem.OEFormat_SDF)
    if not ofs.IsValid():
        raise Exception("output file %s is not valid" % ofs_file)

    while oechem.OEReadMolecule(ifs, molecule):
        # count input file
        count +=1

        if switch: # If True, open new output file
            switch = False
            ofs.close()
            current_letter = letters.pop(0)
            ofs_file = output_f % current_letter
            # Load and check output file
            ofs = oechem.oemolostream(ofs_file)
            if not ofs.IsValid():
                raise Exception("output file %s is not valid" % ofs_file)
            print("Switching to file %s, currently saved %i molecules" % (ofs_file, saved))

        # get isomeric smiles string
        smi = oechem.OECreateIsoSmiString(molecule)
        # if it isn't a new molecule skip it
        if smi in current_smiles:
            smile_count += 1
            continue

        # create and save molecule name in form DrugBank_[letter][number]
        mol_title = molecule_name % (current_letter,count)

        # Make copy before making changes to molecule
        mol_copy = oechem.OEMol(molecule)
        mol_copy.SetTitle(mol_title)
        oechem.OEAddExplicitHydrogens(mol_copy)
        # Determine if molecule meets requirements
        keep = keep_molecule(mol_copy)
        if keep:
            # Add smile to current list and file
            current_smiles.append(smi)
            add_smiles.writelines("%s\t\t%s\n" % (mol_title, smi))
            # write molecule to output file
            oechem.OEWriteMolecule(ofs, mol_copy)
            saved += 1
            # for every 1000 molecules make new output file or 'switch'
            if saved%1000 == 0:
                switch = True

    print("%i molecules in input file" % (count))
    print("%i molecules were had repeated isomeric SMILES" % smile_count)
    print("%i molecules saved to output files" % (saved))

    ofs.close()
    ifs.close()
    add_smiles.close()

    # ========== filter eMolecules and eMolecules_incremental ======
    if os.path.isfile("eMolecules.sdf.gz"):
        eMolecules_filtering("eMolecules.sdf.gz", current_smiles)
    if os.path.isfile("eMolecules_incremental.sdf.gz"):
        eMolecules_filtering("eMolecules_incremental.sdf.gz", current_smiles)


