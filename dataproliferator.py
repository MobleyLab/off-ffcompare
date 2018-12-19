
# coding: utf-8

# In[1]:


# importing
import pandas as pd
import numpy as np
import collections
from itertools import combinations


# In[2]:


# this function reads in the alldata.csv. It expects the 
# format of the alldata.csv produced by the differencemeasure.py
# make heavyatomcount user-input
def read_in(directory,heavyatom):
    alldatadf = pd.read_csv("%s/alldata.csv" % directory, index_col = 0)
    #removing errors
    for column in alldatadf:
        if alldatadf[column].dtype == float:
            alldatadf = alldatadf[alldatadf[column] >= 0]
    # removing molecules above some max HeavyAtomCount
    alldatadf = alldatadf[alldatadf["HeavyAtomCount"] <= heavyatom]
    return alldatadf


# In[1]:


# prints a human-readable statistical summary of each column of data
# make title dependent on heavyatomcount user input
def print_stat(directory,alldataframe,heavyatom):
    statdict = {}
    for column in alldataframe:
        tempstat = alldataframe[column].describe(percentiles = [.25,.5,.75,.95])
        statdict[column] = tempstat
    with open('%s/statistics%d.txt' % (directory, heavyatom),'w') as data:
        data.write(str(statdict))


# In[4]:


# creates all combinations of force fields for use in flagging
def get_ff_combos(alldataframe):
    fflistlist = []
    space = " "
    for columns in alldataframe:
        fflist = columns.split(' ')[-2:]
        if len(fflist) == 2:
            fflist = fflist[0] + space + fflist[1]
        else:
            fflist = fflist[0]
        fflistlist.append(fflist)
    fflist_but_not_ffs = ['MolNames','HeavyAtomCount']
    ffs = set(fflistlist) - set(fflist_but_not_ffs)
    return ffs


# In[2]:


# parses through user inputs
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-i','--tfddc',
            help = "REQUIRED: TFD difference cutoff",
            type = "float",
            dest = 'tfddc')
    
    parser.add_option('-j','--tanisc',
            help = "REQUIRED: TANI same cutoff",
            type = "float",
            dest = 'tanisc')
    
    parser.add_option('-k','--tfdsc',
            help = "REQUIRED: TFD same cutoff",
            type = "float",
            dest = 'tfdsc')

    parser.add_option('-l','--heavy',
            help = "REQUIRED: Heavy Atom Limit",
            type = "float",
            dest = 'heavylim')
    
    parser.add_option('-d','--directory',
            help = "REQUIRED: Name of the full path directory which contains the alldata.csv file",
            type = "string",
            dest = 'directory')
    
    (opt, args) = parser.parse_args()
    if opt.tfddc == None:
        parser.error("ERROR: No TFD difference cutoff was specified.")
    if opt.tanisc == None:
        parser.error("ERROR: No TANI same cutoff was specified.")
    if opt.tfdsc == None:
        parser.error("ERROR: No TFD same cutoff was specified.")
    if opt.directory == None: 
        print("ERROR: No working directory provided.")
    
    directory = opt.directory
    tfddifcutoff = opt.tfddc
    tanisamecutoff = opt.tanisc
    tfdsamecutoff = opt.tfdsc
    heavyatomlimit = opt.heavylim
    
    #runs above functions
    alldatadf = read_in(directory,heavyatomlimit)
    print_stat(directory,alldatadf,heavyatomlimit)
    ffs = get_ff_combos(alldatadf)
    
    #flags molecules as similar or different
    #different molecules are defined different by TFD, same by TaniCo
    #same molecules defined as similar by TFD
    flaggerdicts = {}
    antiflaggerdicts = {}
    for ffcombo in ffs:
        tempflagdatadf = alldatadf[(alldatadf['TFD %s' % ffcombo] > tfddifcutoff) & (alldatadf['TANI %s' % ffcombo] > tanisamecutoff)]
        flaggerdicts[ffcombo] = tempflagdatadf
        tempantiflagdatadf = alldatadf[(alldatadf['TFD %s' % ffcombo] < tfdsamecutoff)]
        antiflaggerdicts[ffcombo] = tempantiflagdatadf
        
    #counts same flags and difference flags for each flagged molecule
    #'extend' appends list-like column to list
    #'Counter' counts frequency of object in list 
    flagcountmolname = []
    for i in flaggerdicts:
        flagcountmolname.extend((flaggerdicts[i]["MolNames"]))
    counter = collections.Counter(flagcountmolname)
    antiflagcountmolname = []
    for i in flaggerdicts:
        antiflagcountmolname.extend((antiflaggerdicts[i]["MolNames"]))
    anticounter = collections.Counter(antiflagcountmolname)
    
    #organizes and exports csvs of data
    antiflagsetdf = pd.DataFrame({'same_flag freq':anticounter})
    flagsetdf = pd.DataFrame({'dif_flag freq':counter})
    antiflagsetdf["MolNames"] = antiflagsetdf.index
    flagsetdf["MolNames"] = flagsetdf.index
    flagsetdf = flagsetdf.merge(alldatadf)
    flagsetdf = flagsetdf.sort_values('dif_flag freq',axis=0,ascending=False)
    flagsetdf = flagsetdf.merge(antiflagsetdf)
    flagsetdf.to_csv('%s/flagset%d.csv' % (directory,heavyatomlimit))
    flagsetlitedf = flagsetdf.select_dtypes(include=['O','int'])
    flagsetlitedf.to_csv('%s/flagsetlite%d.csv' % (directory,heavyatomlimit))

