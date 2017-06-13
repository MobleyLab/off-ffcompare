#!/usr/bin/env python
### Authors:

# Daisy Y. Kyu
# Nam Thi
# Victoria Lim (limvt@uci.edu)
# Caitlin C. Bannan (bannanc@uci.edu)

### Description: This Python script calculates the RMSD values between
#     specified force fields. The values are output into a 'RMSD.txt'
#     or some user-specified output file.

import os
import openeye.oechem as oechem

#-------------------- Functions---------------------#


def RMSD(ref_mol2, query_mol2):
    """
    From one input reference molecule and one input query molecule,
    the RMSD is computed and returned.

    Parameters
    ---------
    ref_mol2: str - mol2 file of the reference force field
    query_mol2: str - mol2 file of the query force field
    
    Returns
    -------
    rms: float - RMSD in Angstroms

    """


    # open reference molecule
    ifsRef = oechem.oemolistream(ref_mol2)
    print ("Opening reference molecule:", ref_mol2)

    # check if the file exist
    if not ifsRef.open(ref_mol2):
        print ("Unable to locate %s. Skipping." % ref_mol2.split('/')[-1])

    # open query molecule
    #queryFile = ("/work/cluster/nthi/ForceField-Comparison/%s/%s/%s" % (homeDir, ffList, fName) )
    print ("Opening query molecule: ", query_mol2  )
    if os.path.exists(query_mol2):
        ifsQuery = oechem.oemolistream(query_mol2)
    # set flavor for the input file
        flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
        ifsRef.SetFlavor(oechem.OEFormat_MOL2, flavor)
        ifsQuery.SetFlavor(oechem.OEFormat_MOL2, flavor)

        # create "blank" object
        rmol = oechem.OEGraphMol()
        qmol = oechem.OEGraphMol()

        # load molecule from files
        oechem.OEReadMolecule(ifsRef,rmol)
        oechem.OEReadMolecule(ifsQuery,qmol)

        # calculate rmsd setting automorph, heavyOnly, and overlay to True
        rms = oechem.OERMSD(rmol,qmol, True, True, True)
    else:
        rms = -2
    return rms

#----------------------------Parse Inputs-----------------------------------#

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-r','--ref',
            help = "REQUIRED: Name(s) of the reference force field(s) separated by commas",
            type = "string",
            dest = 'ref')

    parser.add_option('-c','--compare',
            help = "REQUIRED: Name(s) of the compared force field(s) separated by commas",
            type = "string",
            dest = 'compare')

    parser.add_option('-d','--directory',
            help = "OPTIONAL: Name of the full path directory which contains the force field directories listed in --ref and --compare. If None, then the current working directory is used",
            type = "string",
            dest = 'directory')

    parser.add_option('-o', '--output',
            help = "OPTIONAL: Name of results output file. This will be saved to the directory with force fields. If no output file is specified, RMSD.txt is used.",
            type = 'string',
            dest = 'output',
            default = 'RMSD.txt')

    (opt, args) = parser.parse_args()

    ### Check required fields.
    if opt.ref == None:
        parser.error("ERROR: No force field was specified.")
    if opt.compare == None:
        parser.error("ERROR: No force field was specified.")
    if opt.directory == None:
        print("No working directory provided. Using current directory.")
        directory = os.getcwd()
    else:
        directory = opt.directory


    # set up an error file to record molecule does not exist
    errFile = open('%s/rmsd_errfile.txt' % directory,'a')
    # set up a file to catch negative value
    nValue = open('%s/rmsd_negative_value.txt' % directory,'a')
    # set up a log file for RMSD
    logFile = open('%s/%s' % (directory, opt.output) ,'a')

    # Split up reference and compare force fields
    refFFs = opt.ref.split(',')
    
    for ref in refFFs:
        logFile.write("# Reference Force Field: %s \n" % ref)
        logFile.write("# Molecule Set Directory: %s \n" % directory)

        # get the list of comparison force fields and remove reference
        listFFs = opt.compare.split(',')
        #listFFs.remove(ref)

        refMols = os.listdir(directory + '/' + ref + '/')
        ff_string = "\t".join(['%-9s' % f for f in listFFs])
        logFile.write("%-20s\t%s\n" % ("# MolName", ff_string))

        # loop through each file in the directory and feed them into the function
        for mol2_file in refMols:
            # skip non-mol2 files
            if mol2_file.split('.')[-1] != 'mol2':
                continue
            rms_list = list()
            molName = mol2_file.split('.')[0]
            for queryFF in listFFs:
                ref_file = directory + '/' + ref + '/' + mol2_file
                query_file = directory + '/' + queryFF + '/' + mol2_file
                value =  RMSD(ref_file,query_file)

                # different SMILE strings detected
                if value == -1:
                    oechem.OEThrow.Warning("Negative RMSD value detected for %s: %s" % (molName, value) )
                    nValue.write("%s\t%s\t%s\t%.3e\n" % (molName, ref, query_file, value) )
                    rms_list.append("Neg\t")
                # write mol2 file that does not exist into a file
                elif value == -2:
                    oechem.OEThrow.Warning("Unable to locate %s. Skipping." % query_file)
                    errFile.write("This queryMol does not exist: %s\n" % query_file )
                    rms_list.append("NaN\t")
                else:
                    rms_list.append("%.3e" % value)

            #for each query mol2 file that match reference mol2 file, write out the rms value to the list
            rms_string = "\t".join(rms_list)
            logFile.write("%-20s\t%s\n" % (molName,rms_string))

        logFile.write('#\n')
    errFile.close()
    nValue.close()
    logFile.close()
