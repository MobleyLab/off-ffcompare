#!/usr/bin/env python
### Authors:

# Daisy Y. Kyu
# Nam Thi
# Victoria Lim (limvt@uci.edu)
# Caitlin C. Bannan (bannanc@uci.edu)

### Description: This Python script calculate the RMSD valuesbetween
#     each force field. The rms values is output into a txt files as
#     "name of molecule/referece force field/query force field/rms"
### TODO:
#     1. Add Parse Input to call the script from python i.e python RMSD.py --opls2005...
#     2. Move the text file out of the function?

import os
import openeye.oechem as oechem

#-------------------- Functions---------------------#


def RMSD(ref_mol2, query_mol2):
    """
    From a set of reference molecules and query molecules, the RMSD values are
    computed and wrote into a file.

    Parameters
    ---------
    ref_mol2: str - mol2 file of the reference force field
    query_mol2: str - mol2 file of the query force field


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

        # rmsd calculate
        rms = oechem.OERMSD(rmol,qmol, True, True , True)
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
            help = "OPTIONAL: Name of output file with results, it will be saved to the specified directory. If None provided RMSD.txt is used",
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
        print("no working directory provided using current directory")
        directory = os.getcwd()
    else:
        directory = opt.directory


    # set up an error file to record molecule does not exist
    errFile = open('%s/rmsd_logfile.txt' % directory,'a')
    # set up a file to catch negative value
    nValue = open('%s/negative_value.txt' % directory,'a')
    # set up a log file for RMSD
    logFile = open('%s/%s' % (directory, opt.output) ,'a')

    # Split up reference and compare force fields
    refFFs = opt.ref.split(',')
    listFFs = opt.compare.split(',')

    for ref in refFFs:
        logFile.write("# Reference Force Field: %s \n" % ref)
        logFile.write("# Molecue Set Directory: %s \n" % directory)

        refMols = os.listdir(directory + '/' + ref + '/')
        ff_string = "\t".join(listFFs)
        logFile.write("# MolName \t %30s \n" % ff_string)

        # loop through each file in the directory and feed them into the funciton
        for mol2_file in refMols:
            rms_list = list()
            molName = mol2_file.split('.')[0]
            for queryMol in listFFs:
                ref_file = directory + '/' + ref + '/' + mol2_file
                query_file = directory + '/' + queryMol + '/' + mol2_file
                value =  RMSD(ref_file,query_file)

                # different SMILE strings detected
                if value == -1:
                    oechem.OEThrow.Warning("Negative RMSD's value detected for %s-%s" % (molName, value) )
                    nValue.write("%s\t%s\t%s\t%.3e\n" % (molName, ref, queryMol, value) )
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
            logFile.write("%5s\t%s\n" % (molName,rms_string))

        logFile.write('#\n')
    errFile.close()
    nValue.close()
    logFile.close()
