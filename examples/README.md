# Examples

Here we show an example of how to use the provided scripts
on a set of 8 Alkanes, Ethers, and Alcohols (AlkEthOH\_chain\_tiny).

Below is a list of command line calls that were used
and a list of output files/directories for each.

* `AlkEthOH_chain_tiny.smi` - smiles strings for this molecules

### SMILES to SDF
```
python smi2sdf.py --smiles examples/AlkEthOH_chain_tiny.smi \
--sdf examples/AlkEthOH_chain_tiny.sdf
```

* `AlkEthOH_chain_tiny.sdf` - SDF file with all molecules in input smiles file

### Generating input files
```
python genMOL2.py -i examples/AlkEthOH_chain_tiny.sdf -l examples/
```

* `tripos\_mol2/` - directory of Tripos mol2 files
* `gaff\_mol2/` - directory of input GAFF files
* `gaff2\_mol2/` - directory of input GAFF2 files
* `timer.dat` - timer data for creating mol2 files

### Minimizing Molecules
```
python min_oe_openMM.py --inmols examples/tripos_mol2/ \
    --ffxml smirnoff99Frosst.ffxml \
    --gaffdir examples/gaff_mol2/  \
    --gaff2dir examples/gaff2_mol2/ \
    --dommff True \
    --outdir examples/
```

* `GAFF(2)/` - minimized mol2 files with GAFF(2)
* `MMFF94(S)/` - minimized mol2 files with MMFF94(S)
* `SMIRNOFF/` - minimized mol2 files with SMIRNOFF
* `output.dat` - output file reporting any errors with min_oe_openMM.py

```
python OPLS.py --idir examples/tripos_mol2/ \
    --dir2005 examples/OPLS2005/ \
    --dir3 examples/OPLS3/
```

* `OPLS3` - minimized mol2 files with OPLS3
* `OPLS2005` - minimized mol2 files with OPLS2005

### RMSD calculations

```
python RMSD.py \
    --ref SMIRNOFF,GAFF,GAFF2,MMFF94,MMFF94S,OPLS2005,OPLS3 \
    --compare SMIRNOFF,GAFF,GAFF2,MMFF94,MMFF94S,OPLS2005,OPLS3 \
    --directory examples/ 
```

* `RMSD.txt` - complete data table for RMSDs comparing each pair of forcefields for each molecule
* `negavite_value.txt` - file that would enumerate molecules that did not match from different force field
* `rmsd_errfile.txt` - file that would enumerate errors with opening minimized mol2 files if they existed 
