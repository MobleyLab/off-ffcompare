# Comparing single molecule minimizations with different all atom force fields

**NOTE: only Python 2 is supported for the time being.**

### Project Goals

This project relates to the larger goals of the 
[Open Force Field Group's](https://github.com/open-forcefield-group)
effort to automate force field parameterization. 

Here we have created scripts to minimize small molecules using a variety of force fields with two goals in mind. 
1. Find places where the SMIRNOFF99Frosst force field behaves differently from currently accepted force fields.
2. Discover molecules that are minimized to different conformations by different force fields. These molecules will likely be used in future SMIRNOFF parameterizations. 

### Force fields 

The follow is a list of the force fields being considered here:

* [SMIRNOFF99Frosst](https://github.com/open-forcefield-group/smirff99Frosst)
* [GAFF](http://ambermd.org/antechamber/gaff.html)
* [GAFF2](https://mulan.swmed.edu/group/gaff.php)
* [MMFF94](http://open-babel.readthedocs.io/en/latest/Forcefields/mmff94.html)
* [MMFF94S](http://open-babel.readthedocs.io/en/latest/Forcefields/mmff94.html)
* [OPLS3](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00864)
* [OPLS2005](http://dx.doi.org/10.1002/jcc.20292)

### Contents

* *ffcompare.py*: Script to read in mol2 files and run minimization using OpenMM (or oechem.OESzybki).

* *genTriposGAFFandGAFF2.py*: Script to generate mol2 files with Tripos, GAFF, and GAFF2 atom types.

* *smirff99Frosst.ffxml*: FFXML file for SMIRFF

* *smi2sdf.py*: Script to generate a collection SDF file from a list of SMILES strings for use in genTriposGAFFandGAFF2.py

* *filter_molecules.py* : This script was used to filter the eMolecules and DrugBank databases to meet the requirements for this project.  

___

### Instructions

1. Generate initial SDF files.
    1. filter_molecules.py was used to filter DrugBank and eMolecules databases available online. The method `eMolecule_filtering` in that script could be used to filter any other large database of molecules to meet these requirements for all molecules:
        * < 200 heavy atoms
        * no metals
        * proper valency, that is no first row elements with > 5 bonds
        * `python filter_molecules.py` assuming DrugBank.sdf, eMolecules.sdf.gz and eMolecules_incremental.sdf.gz are in the current directory. 
    2. If starting from a list of SMILES strings, can generate a single conformer with oechem
    * [VTL will fill this in]

2. Generate Tripos, GAFF, and GAFF2 mol2 files.
    * `python genTriposGAFFandGAFF2.py -i /path/to/sdf/files -l /path/to/output/files` 
3. Perform minimization.
   1. ffcompare.py was used to minimize mol2 files with a specified forcefield type (fftype). Supported fftypes include GAFF, GAFF2, MMFF94, MMFF94S and SMIRFF. Inmols should specify path to tripos mol2.
   * SMIRFF: `python ffcompare.py --fftype smirff --ffxml smirff99Frosst.ffxml --inmols /location/to/triposmol2Files`
   * GAFF: `python ffcompare.py --fftyle gaff --inmols /location/to/triposmol2Files -g /location/to/gaff_inpcrd_prmtopFiles`
   * GAFF2: `python ffcompare.py --fftype gaff2 --inmols /location/to/triposmol2Files -g /location/to/gaff2_inpcrd_prmtopFiles`
   * MMFF94 & MMFF94S: `python ffcompare.py --fftype mmff --inmols /location/to/triposmol2Files`
   2. [Nam fill this in with description]
   * [Nam fill this in with python command]
      
4. Evaluate RMSD.
   * [Nam fill this in with the python command]
