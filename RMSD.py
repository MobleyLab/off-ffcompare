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


def RMSD(ffRef, ffList, molSet, molName):
    """
    From a set of reference molecules and query molecules, the RMSD values are
     computed and wrote into a text file.  
   
    Parameters
    ---------
    ffRef: str - the reference forcefield
    ffList: str - the query forcefield
    molSet: str - set of molecule
    molName: str - name of the molecules
 
    """
    # Set up an error file to record molecule does not exist
    errFile = open('rmd_logfile.txt','a')

    # Set up a file that catch negative values
    nValue = open('negative_value.txt','a') 
  
    # Set up a log file for RMSD
    logFile = open('RMSD.txt', 'a')
   
 
    # open reference molecule
    refFile = ("/work/cluster/nthi/ForceField-Comparison/%s/%s/%s" % (molSet, ffRef, molName) )
    ifsRef = oechem.oemolistream(refFile)
    print ("Opening reference molecule:", refFile)

    # check if the file exist
    if not ifsRef.open(refFile):
        print ("Unable to locate %s. Skipping." % refFile)

    # open query molecule
    queryFile = ("/work/cluster/nthi/ForceField-Comparison/%s/%s/%s" % (molSet, ffList, molName) )
    print ("Opening query molecule: ", queryFile)   
    if os.path.exists(queryFile):
        ifsQuery = oechem.oemolistream(queryFile)
   
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
        molNames = rmol.GetTitle()
        logFile.write("%s\t%s\t%s\t%.3e\n" % (molNames, ffRef,ffList, rms) )
       
        # different SMILE strings detected
        if rms == -1:
            oechem.OEThrow.Warning("Negative RMSD's value detected for %s-%s" % (molNames, ffList) )  
            nValue.write("%s\t%s\t%s\t%.3e\n" % (molNames, ffRef, ffList, rms) ) 

    # write mol2 file that does not exist into a file
    else: 
        oechem.OEThrow.Warning("Unable to locate %s. Skipping." % queryFile)
        print ("Unable to locate %s. Skipping." % molName)
        errFile.write("This queryMol does not exist: %s\t%s\t%s\tNaN\n" % (ffList,molSet,molName) )
    errFile.close()
    nValue.close()
    logFile.close()

#------------------------------Script----------------------------------------#

refFile = '/work/cluster/nthi/ForceField-Comparison/DrugBank/SMIRFF'
refMols = os.listdir(refFile)

listFFs = ['OPLS2005','OPLS3','GAFF2','GAFF']

# loop through each file in the directory and feed them into the funciton
for fName in refMols:
    for queryFF in listFFs:
        value = RMSD('SMIRFF',queryMol,'DrugBank',fName)
        print (value)
        
#----------------------------Parse Inputs-----------------------------------#

