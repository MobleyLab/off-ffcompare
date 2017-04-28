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

The following is a list of the force fields being considered here:

* [SMIRNOFF99Frosst](https://github.com/open-forcefield-group/smirnoff99Frosst)
* [GAFF](http://ambermd.org/antechamber/gaff.html)
* [GAFF2](https://mulan.swmed.edu/group/gaff.php)
* [MMFF94](http://open-babel.readthedocs.io/en/latest/Forcefields/mmff94.html)
* [MMFF94S](http://open-babel.readthedocs.io/en/latest/Forcefields/mmff94.html)
* [OPLS3](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00864)
* [OPLS2005](http://dx.doi.org/10.1002/jcc.20292)

### Contents

* `min_oe_openMM.py`: Script to read in mol2 files and run minimization using OpenMM (or oechem.OESzybki).

* *genTriposGAFFandGAFF2.py*: Script to generate mol2 files with Tripos, GAFF, and GAFF2 atom types.

* *smirnoff99Frosst.ffxml*: FFXML file for SMIRNOFF

* *smi2sdf.py*: Script to generate a collection SDF file from a list of SMILES strings for use in genTriposGAFFandGAFF2.py

* *filter_molecules.py* : This script was used to filter the eMolecules and DrugBank databases to meet the requirements for this project.  

___

### Instructions

1. Generate initial SDF files.
    1. Use filter_molecules.py to filter DrugBank and eMolecules databases available online. The method `eMolecule_filtering` in that script could be used to filter any other large database of molecules to meet these requirements for all molecules:
        * < 200 heavy atoms
        * no metals
        * proper valency, that is no first row elements with > 5 bonds
        * `python filter_molecules.py` assuming DrugBank.sdf, eMolecules.sdf.gz and eMolecules_incremental.sdf.gz are in the current directory. 
    2. If starting from a list of SMILES strings, can generate a single conformer with oechem using smi2sdf.py
        * Change variables, then call `python smi2.sdf.py`.

2. Generate Tripos, GAFF, and GAFF2 mol2 files.
    * `python genTriposGAFFandGAFF2.py -i /path/to/sdf/files -l /path/to/output/files` 
3. Perform minimization.
   1. Use oe_min_openMM.py to minimize mol2 files with a specified forcefield type (fftype). Supported fftypes include GAFF, GAFF2, MMFF94, MMFF94S and SMIRNOFF. The inmols flag should specify path to tripos \*.mol2 files (not the mol2 files themselves).

    * This script can be used for one or more of the accepted forcefields
    * All: `python oe_min_openMM.py --inmols /path/triposMol2Files/ --dommff True --ffxml smirnoff99Frosst.ffxml --gaffdir /path/GAFF/Files/ --gaff2dir /path/GAFF2/Files` 
    * GAFF: `python oe_min_openMM.py --inmols /path/triposMol2Files/ --gaffdir /path/GAFF/Files/`  

   2. OPLS.py was used to minimize mol2 files with a specific forcefield type (fftype). Supported fftypes inclue OPLS3, OPLS2005. Input should specifiy path to tripos mol2 files and optimizetype should specify the forcefield type.
   * OPLS3: `python OPLS.py --input /path/to/triposmol2Files --optimizetype "OPLS3" --outdir /path/output_directory/`
   * OPLS2005: `python OPLS.py --input /path/to/triposmol2Files --optimizetype "OPLS2005" --outdir /path/output_directory/`
      
4. Evaluate RMSD.
   * `python --ref [name of reference force-field] --compare [names of compared force fields] --directory /path/containing/directories`
