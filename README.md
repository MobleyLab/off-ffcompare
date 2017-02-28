# Comparing single molecule minimizations with different all atom force fields

**NOTE: only Python 2 is supported for the time being.**

### Project Goals

This project relates to the larger goals of the 
[Open Force Field Group's](https://github.com/open-forcefield-group)
effort to automate force field parameterization. 

Here we have created scripts to minize small molecules using a variety of force fields with two goals in mind. 
1. Find places where the SMIRNOFF99Frosst force field behaves differently from currently accepted force fields.
2. Discover molecules that are minized to different conformations by different force fields. These molecules will likely be used in future SMIRNOFF parameterizations. 

### Force fields 

The follow is a list of the force fields being considered here:

* [SMIRNOFF99Frosst](https://github.com/open-forcefield-group/smirff99Frosst)
* GAFF
* GAFF2
* MMFF94
* MMFF94S
* OPLS3
* OPLS2005

### Contents

* *ffcompare.py*: Script to read in mol2 files and run minimization using OpenMM (or oechem.OESzybki).

* *genTriposGAFFandGAFF2.py*: Script to generate mol2 files with Tripos, GAFF, and GAFF2 atom types.

* *smirff99Frosst.ffxml*: FFXML file for SMIRFF

* *smi2sdf.py*: Script to generate a collection SDF file from a list of SMILES strings for use in genTriposGAFFandGAFF2.py

___________________________________________________________________________

### Instructions

1. Generate initial SDF files.
    a. [Caitlin will include stuff on filtering DrugBank && eMolecules]
       * example command that Caitlin used
    b. If starting from a list of SMILES strings, can generate a single conformer with oechem
       * [VTL will fill this in]
   
2. Generate Tripos, GAFF, and GAFF2 mol2 files.
   * python genTriposGAFFandGAFF2.py -i /path/to/sdf/files -l /path/to/output/files
   
3. Perform minimization.
   a. [Daisy fill this in]
      * python ffcompare.py --fftype smirff --ffxml smirff99Frosst.ffxml --inmols /location/to/mol2Files
   b. [Nam fill this in with description]
      * [Nam fill this in with python command]
      
4. Evaluate RMSD.
   * [Nam fill this in with the python command]
