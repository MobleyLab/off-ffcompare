from openeye import oechem
import string
from openeye import oeomega
import smarty
from filter_molecule_methods import *

if __name__=="__main__":
    # Check input files
    input_f = "eMolecules.sdf.gz"
    output_f = "eMolecules_%i.sdf"
    smiles_base = "eMolecules_%i.txt"
    molecule_name = "eMolecules_%i_%i"

    # get current smiles
    smiles_f = "smiles_to_ID_off-compare.txt"
    current_smiles = parse_smile(smiles_f)

    # Load and check input file
    ifs = oechem.oemolistream(input_f)
    if not ifs.IsValid():
        raise Exception("Error: input_file (%s) was not valid" % input_f)

    smirks_file = "remove_smirks.smarts"
    smirks = smarty.utils.parse_odds_file(smirks_file, False)
    remove_smirks = smirks[0]

    errs = oechem.oeosstream()
    oechem.OEThrow.SetOutputStream(errs)

    molecule = oechem.OECreateOEGraphMol()
    count = 0
    warnings = 0
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
        # count input file
        count +=1
        if ("warning" in errs.str().lower()):
            warnings += 1
            errs.clear()
            continue

        if switch:
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

        smi = oechem.OECreateIsoSmiString(molecule)
        mol_copy = oechem.OEMol(molecule)
        oechem.OEAddExplicitHydrogens(mol_copy)
        if smi in current_smiles:
            smile_count += 1
            continue

        keep = keep_molecule(mol_copy, remove_smirks)
        if keep:
            mol_title = molecule_name % (current_letter,count)
            mol_copy.SetTitle(mol_title)
            add_smiles.writelines("%s\t\t%s\n" % (mol_title, smi))
            oechem.OEWriteMolecule(ofs, mol_copy)
            saved += 1
            if saved%1000 == 0:
                switch = True
        errs.clear()

    print("%i molecules in input file" % (count))
    print("%i molecules resulted in warnings when parsing" % warnings)
    print("%i molecules were had repeated isomeric SMILES" % smile_count)
    print("%i molecules saved to output files" % (saved))

    ifs.close()

