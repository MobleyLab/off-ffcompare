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
        valence = atom.GetValence()
        if atomNum <= 10:
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
    if oechem.OECount(mol, oechem.OEIsMetal()) > 0:
        return False
    if oechem.OECount(mol, oechem.OEIsHeavy()) > 200:
        return False
    for smirks in remove_smirks:
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smirks):
            continue
        ss = oechem.OESubSearch(qmol)
        matches = [match for match in ss.Match(mol, False)]
        if len(matches) > 0:
            return False
    return check_valence(mol)

def parse_smile(smile_f):
    """
    Parses file with the form "molecule_name    SMILES"
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
    smiles = [l.split()[1] for l in lines]
    return smiles

