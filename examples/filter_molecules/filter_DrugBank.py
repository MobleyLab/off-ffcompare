from openeye import oechem
import string
from openeye import oeomega
import smarty
from filter_molecule_methods import *

if __name__=="__main__":
    # Check input files
    input_f = "DrugBank.sdf"
    output_f = "DrugBank_%s.sdf"
    molecule_name = "DrugBank_%s%i"

    # get letters for output file from , 'a', 'aa', 'aaa', 'aaaa' and uppercase
    letters = string.ascii_letters
    letters = [l for l in letters]
    letters += [l+l for l in letters]
    letters += [l+l for l in letters]

    # get current smiles
    smiles_f = "smiles_to_ID_off-compare.txt"
    current_smiles = parse_smile(smiles_f)
    add_smiles = open(smiles_f, 'a')

    # Load and check input file
    ifs = oechem.oemolistream(input_f)
    ifs.SetFormat(oechem.OEFormat_SDF)
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
    current_letter = letters.pop(0)
    ofs_file = output_f%current_letter
    ofs = oechem.oemolostream(ofs_file)
    ifs.SetFormat(oechem.OEFormat_SDF)
    if not ofs.IsValid():
        raise Exception("output file %s is not valid" % ofs_file)

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
            current_letter = letters.pop(0)
            ofs_file = output_f % current_letter
            # Load and check output file
            ofs = oechem.oemolostream(ofs_file)
            if not ofs.IsValid():
                raise Exception("output file %s is not valid" % ofs_file)
            print("Switching to file %s, currently saved %i molecules" % (ofs_file, saved))

        smi = oechem.OECreateIsoSmiString(molecule)
        mol_title = molecule_name % (current_letter,count)
        molecule.SetTitle(mol_title)
        mol_copy = oechem.OEMol(molecule)
        oechem.OEAddExplicitHydrogens(mol_copy)
        if smi in current_smiles:
            smile_count += 1
            continue

        keep = keep_molecule(mol_copy, remove_smirks)
        print("Before Title = %s" % mol_copy.GetTitle())
        print("Before MCMolTitle = %s" % mol_copy.GetMCMolTitle())
        if keep:
            molWithConfs = oechem.OEMol()
            oechem.OEParseSmiles(molWithConfs, smi)
            print("Begin defining omega...")
            omega = oeomega.OEOmega()
            omega.SetMaxConfs( 1)
            omega.SetStrictStereo( False)
            omega.SetSampleHydrogens( True)
            print("run omega on 2d structure...")
            status = omega(molWithConfs)
            print("status is %s" % status)
            if status:
                current_smiles.append(smi)
                print("After Title = %s" % mol_copy.GetTitle())
                print("After MCMolTitle = %s" % mol_copy.GetMCMolTitle())
                add_smiles.writelines("%s\t\t%s\n" % (mol_title, smi))
                oechem.OEWriteMolecule(ofs, mol_copy)
                print("Saved molecule...")
                saved += 1
                if saved%1000 == 0:
                    switch = True
        else:
            print("Molecule not saved %s" % mol_copy.GetTitle())
        errs.clear()

    print("%i molecules in input file" % (count))
    print("%i molecules resulted in warnings when parsing" % warnings)
    print("%i molecules were had repeated isomeric SMILES" % smile_count)
    print("%i molecules saved to output files" % (saved))

    ifs.close()

