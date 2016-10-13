import os
import openmoltools
import commands
import openeye.oechem as oechem
import openeye.oeomega as oeomega

def make_path(filename):
    try:
        path = os.path.split(filename)[0]
        os.makedirs(path)
    except OSError:
        pass


# read in SDF file with all molecules
sdf_File = "./DrugBank_d.sdf"
ifs = oechem.oemolistream(sdf_File)

# Set up an error file to record skipped troublesome molecules
errfile = open('errMols.txt', 'a')

for xmol in ifs.GetOEMols():
    
    # generate a single conformer with any stereochem
    omega = oeomega.OEOmega()
    omega.SetStrictStereo(False)
    omega.SetMaxConfs( 1 )
    if omega(xmol): 
        mol = xmol
    else:
        errfile.write("%s\t%s" % (sdf_File, xmol.GetTitle()))
        continue

    molName = mol.GetTitle()

    # set filenames for tripos files (mol2)
    tripos_filename = os.path.join("tripos_mol2", molName+".mol2" )

    # set filenames for gaff files (mol2, frcmod, prmtop, inpcrd)
    gaff_mol2_filename = os.path.join("./gaff_mol2/", "%s.mol2" % molName)
    frcmod_filename = os.path.join("./gaff_mol2/", "%s.frcmod" % molName)
    prmtop_filename = os.path.join("./gaff_mol2/", "%s.prmtop" % molName)
    inpcrd_filename = os.path.join("./gaff_mol2/", "%s.inpcrd" % molName)

    # make subdirectories for tripos and gaff files
    make_path( 'tripos_mol2/' )
    make_path( 'gaff_mol2/')

    # calculate charge
    chgmol = openmoltools.openeye.get_charges(mol) 

    # reset the title to original title (get_charges renames mol)
    chgmol.SetTitle(mol.GetTitle())

    # generate tripos mol2 files
    openmoltools.openeye.molecule_to_mol2(chgmol, tripos_filename)

    # generate the gaff mol2 file and the .frcmod file
    _, _ = openmoltools.amber.run_antechamber(molName, tripos_filename,charge_method=None, gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)

    # generate the gaff .inpcrd and .prmtop files
    openmoltools.amber.run_tleap(molName, gaff_mol2_filename, frcmod_filename, prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)

    # generate gromacs .top and .gro files
    #openmoltools.utils.convert_via_acpype( molName, prmtop_filename, inpcrd_filename)

    commands.getoutput('rm ./gaff_mol2/*.frcmod')
ifs.close()
errfile.close()
