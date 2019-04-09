"""This script takes a directory full of mol2 files, and generates a dataframe
of quantities of functional groups present within all mol2 files in 
the directory. It outputs a dataframe enumerating single checkmol descriptors,
called onefuncdict.csv. Then, it outputs a dataframe enumerating combinations 
of two checkmol descriptors, called twofuncdict.csv. 
"""

# Importing 
from openeye import oechem  
import os
import pandas as pd  
import numpy as np
from openmoltools.utils import get_checkmol_descriptors
import itertools
import re

if __name__ == '__main__':
    from optparse import OptionParser
    # user input: directory that includes mol2 files of interest
    parser = OptionParser()
    parser.add_option('-d','--directory',
            help = "REQUIRED: directory",
            type = "string",
            dest = 'directory')
    (opt,args) = parser.parse_args()
    if opt.directory == None: 
        print("ERROR: No working directory provided.")
    directory = opt.directory
    # creating empty dictionary
    desdict = {}
    # loading dictionary with descriptors for every molecule in directory
    for molfile in os.listdir('%s/' % directory):
        tempdes = get_checkmol_descriptors('%s/%s' % (directory,molfile))
        desdict[molfile] = tempdes
    # loading dictionary into dataframe, and formatting dataframe
    descdf = pd.DataFrame.from_dict(desdict, orient = "index")
    descdf[0] = descdf[0].apply(tuple)
    descdf['MolNames'] = descdf.index
    descdf = descdf.reset_index(drop=True)
    descdf['Descriptors'] = descdf[0]
    descdf = descdf.drop(0, axis = 'columns')
    # creating new dataframe to start enunmerating values
    totdf = descdf
    listoflist = totdf["Descriptors"].values
    mergedlist = list(itertools.chain.from_iterable(listoflist))
    noreplist = list(set(mergedlist))
    # formatting new dataframe to make descriptors easily countable
    totdf['liststring'] = totdf['Descriptors'].apply(lambda x:\
     ','.join(map(str, x)))
    totdf['liststring'] = totdf['liststring'].apply(lambda x:''.join(e\
    for e in x if e.isalnum()))
    # counting combinations of descriptors
    twofuncdict = {}
    for func in noreplist:
        print(func)
        for func2 in noreplist:
            func = ''.join(e for e in func if e.isalnum())
            func2 = ''.join(i for i in func2 if i.isalnum())
            twofuncdict[func + '_' + func2] =\
            len(totdf[(totdf['liststring'].str.contains(func)) &\
            (totdf['liststring'].str.contains(func2))])
    # exporting csv
    twofuncdictdf = pd.DataFrame.from_dict(twofuncdict, orient ='index')
    twofuncdictdf.to_csv('twofuncdict.csv')

    # counting single descriptors 
    onefuncdict = {}
    for func in noreplist:
        func = ''.join(e for e in func if e.isalnum())
        onefuncdict[func] = len(totdf[totdf['liststring'].str.contains(func)])
    # exporting csv
    onefuncdictdf = pd.DataFrame.from_dict(onefuncdict, orient ='index')
    onefuncdictdf.to_csv('onefuncdict.csv')

