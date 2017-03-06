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

* *filter_molecule.py* : This script was used to filter the eMolecules and DrugBank databases to meet the requirements for this project.  

___

### Instructions

1. Generate initial SDF files.
    1. filter_molecule.py was used to filter DrugBank and eMolecules databases available online. The method `eMolecule_filtering` in that script could be used to filter any other large database of molecules to meet these requirements for all molecules:
        * < 200 heavy atoms
        * no metals
        * proper valency, that is no first row elements with > 5 bonds
    2. If starting from a list of SMILES strings, can generate a single conformer with oechem
    * [VTL will fill this in]

2. Generate Tripos, GAFF, and GAFF2 mol2 files.
    * `python genTriposGAFFandGAFF2.py -i /path/to/sdf/files -l /path/to/output/files` 
3. Perform minimization.
   1. [Daisy fill this in]
   * python ffcompare.py --fftype smirff --ffxml smirff99Frosst.ffxml --inmols /location/to/mol2Files
   2. [Nam fill this in with description]
   * [Nam fill this in with python command]
      
4. Evaluate RMSD.
   * [Nam fill this in with the python command]
