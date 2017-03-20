#!/usr/bin/env python
### Author:
# Daisy Y. Kyu
# Nam Thi
# Victoria Lim (limvt@uci.edu)
# Caitlin C. Bannan (bannanc@uci.edu)

### Description: This Python script minimize a mol2 files 
#     with specific forcefield type and location provided.
#     For N input *.mol2 files, there should be N output
#     *.mol2 files, unless the input file is missing


### Example usuage:
# - python OPLS.py --input /input_directory/mol2 --optimizetype fftype > output.dat
# - python OPLS.py -i /input_directory/OPLS2005/mol2 -o OPLS2005  > output.dat

### Dependencies: 
# - For all: needs fftype and location of mol2 files. (-i, -o)
# - cmd command: refer to the help message in utilities on 
#   specifying structures from other file formats.


### Supported force field:
# - OPLS3
# - OPLS2005

import os

  # ---------------------------- Functions ------------------------- #

def OPLSMin(ifile,opttype):
    """
   
    Take one mol2 file and do the minimization, then output into
    a specific location
   
    Parameters:
    ifile: input file directory
    opttpe: force field type
    Note: the output directory is specified inside the command line.
          Refer to the help message in utilities for more information.
    
   
    """

    if opttype == 'OPLS3':  #specify which fftype to be run

        cmd = '$SCHRODINGER/utilities/ffld_server  -version 16 -charges_from_ct -virt \
               -no_cm1a_bcc -opt -BFGS -imol2 /input_location/%s \ 
               -omol2 /output_location/%s'% (ifile,ifile.split('/')[-1])

        os.system(cmd)      #execute the command from the shell

    if opttype == 'OPLS2005':

        cmd = '$SCHRODINGER/utilities/ffld_server -version 14 -charges_from_ct -virt \
               -no_cm1a_bcc -opt -BFGS -imol2 /input_location/%s \ 
               -omol2 /output_location/%s' % (ifile,ifile.split('/')[-1])

        os.system(cmd)

  # ------------------------- Parse Inputs ----------------------- #


if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-i','--input',
            help = "Path to directory containing all mol2 files to be minimized.",
            type = "string",
            dest = 'ifile')

    parser.add_option('-o','--optimizetype',
            help = "Name of the force field type to be used for minimization",
            type = "string",
            dest = 'opttype')
    (opt, args) = parser.parse_args()
    OPLSMin(opt.ifile,opt.opttype)


