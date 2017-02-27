
    ==============
    |  CONTENTS  |
    ==============
    
Last updated: Feb 27 2017
NOTE: only Python 2 is supported for the time being.


### Contents

ffcompare.py
    Script to read in mol2 files and run minimization using OpenMM (or oechem.OESzybki).

genTriposGAFFandGAFF2.py
    Script to generate mol2 files with Tripos, GAFF, and GAFF2 atom types.
    TODO: remove anything related to 'commands' module.

smirff99Frosst.ffxml
    FFXML file for SMIRFF
    We talked about everyone pointing to a master copy of the .ffxml file though. (tbd)

smi2sdf.py
    Script to generate a collection SDF file from a list of SMILES strings.
    This can be used to feed into genTriposGAFFandGAFF2.py

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
