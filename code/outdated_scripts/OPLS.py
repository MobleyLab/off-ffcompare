#!/usr/bin/env python

### Author:
# Nam Thi
# Caitlin C. Bannan (bannanc@uci.edu)
# Victoria Lim (limvt@uci.edu)

### Description: This Python script minimizes mol2 files in the
#     given directory with the specified force field. Minimizations
#     are completed with the ffld_server utility from Schrodinger.
#   For N input *.mol2 files, there should be N output
#     *.mol2 files for which ever force field(s) are specified.


### Example usuage:
# - python OPLS.py --idir /input_directory/with/mol2s --dir2005 output/directory/OPLS2005 --dir3 output/directory/OPLS3 > output.dat
# - python OPLS.py -i /input_directory/with/mol2s -d output/directory/OPLS2005 -D output/directory/OPLS3 > output.dat

### Note:
# - Currently compatible with mol2 file types.
#   Call `$SCHRODINGER/utilities/ffld_server -hh` for more options.

### Supported force fields:
# - OPLS3
# - OPLS2005

import os
import glob
import sys

  # ---------------------------- Functions ------------------------- #

def OPLS3_minimize(in_mol2, out_mol2):
    """
    Performs a minimization using OPLS3 on the input mol2 file.
    Assumes the input and output files have already been checked.

    Parameters
    ----------
    in_mol2: string, input mol2 file location
    out_mol2: string, output mol2 file location
    """
    cmd = '$SCHRODINGER/utilities/ffld_server  -version 16 -charges_from_ct \
            -virt -no_cm1a_bcc -opt -BFGS -imol2 %s -omol2 %s'\
            % (in_mol2, out_mol2)
    os.system(cmd)
    return

def OPLS2005_minimize(in_mol2, out_mol2):
    """
    Performs a minimization using OPLS2005 on the input mol2 file.
    Assumes the input and output files have already been checked.

    Parameters
    ----------
    in_mol2: string, input mol2 file location
    out_mol2: string, output mol2 file location
    """
    cmd = '$SCHRODINGER/utilities/ffld_server -version 14 -charges_from_ct \
            -virt -no_cm1a_bcc -opt -BFGS -imol2 %s -omol2 %s'\
            % (in_mol2, out_mol2)

    os.system(cmd)

  # ------------------------- Parse Inputs ----------------------- #


if __name__ == '__main__':
    from optparse import OptionParser

    usage_string="""\
            This script is used to minimize molecules in mol2 files
            with OPLS3 or OPLS2005. You must have an environment variable
            $SCHRODINGER in order for the minimizations to work.

            usage: python OPLS.py --idir [path to mol2 directory]
            --dir2005 [path to OPLS2005 output directory]
            --dir3 [path to OPLS3 output directory]

            if both dir2005 and dir3 are None then no minimizations occur.
            If dir2005 or dir3 are not None, but the directory doesn't exist
            then it is created.
            """

    parser = OptionParser(usage=usage_string)

    parser.add_option('-i','--idir',
            help = "REQUIRED: Path to directory containing all mol2 files to be minimized.",
            type = "string",
            dest = 'idir')

    parser.add_option('-d','--dir2005',
            help = "OPTIONAL: Directory for OPLS2005 minimization results, required for OPLS2005 minimization",
            type = "string",
            dest = 'dir2005')

    parser.add_option('-D', '--dir3',
            help = "OPTIONAL: Directory for OPLS3 minimization results, required for OPLS3 minimization",
            type = "string",
            dest = 'dir3')

    (opt, args) = parser.parse_args()
    # Check input directory
    if opt.idir is None:
        parser.print_help()
        parser.error("ERROR: you must provide an input directory")
    if not os.path.isdir(opt.idir):
        parser.print_help()
        parser.error("ERROR: input directory (%s) does not exist" % opt.idir)

    # Check that at least one output directory exists
    if opt.dir3 is None and opt.dir2005 is None:
        parser.print_help()
        parser.error("ERROR: must provide at least one output directory for minimiztion to occur")

    # check the the not None directory exists
    if opt.dir3 is not None and (not os.path.isdir(opt.dir3)):
        os.mkdir(opt.dir3)
    # if it doesn't exist create it
    if opt.dir2005 is not None and (not os.path.isdir(opt.dir2005)):
        os.mkdir(opt.dir2005)

    # check schrodinger utilities tool directory at least exists
    try:
        os.environ.get('SCHRODINGER')
    except:
        parser.print_help()
        parser.error("ERROR: cannont find environment variable $SCHRODINGER. Please add it before continuing.")

    # Get Mol2 Files!
    mol2_files = glob.glob("%s/*.mol2" % opt.idir)
    # check that you have some mol2 files
    if len(mol2_files) == 0:
        print("No mol2 files found in input directory (%s)" % opt.idir)
        sys.exit()

    # Loop through all input mol2 files
    for mol2 in mol2_files:
        mol2_basename = mol2.split('/')[-1]
        print("Working on %s" % mol2_basename)

        if opt.dir2005 is not None:
            out_mol2 = "%s/%s" % (opt.dir2005, mol2_basename)
            # check if out_mol2 exists:
            if os.path.isfile(out_mol2) and os.path.getsize(out_mol2) > 0:
                print("Skipping OPLS2005 output (%s) already exists and is not empty." % out_mol2)
            else:
                print("Minimizing with OPLS2005...")
                OPLS2005_minimize(mol2, out_mol2)

        if opt.dir3 is not None:
            out_mol2 = "%s/%s" % (opt.dir3, mol2_basename)
            # check if out_mol2 exists:
            if os.path.isfile(out_mol2) and os.path.getsize(out_mol2) > 0:
                print("Skipping OPLS3 output (%s) already exists and is not empty." % out_mol2)
            else:
                print("Minimizing with OPLS3...")
                OPLS3_minimize(mol2, out_mol2)


