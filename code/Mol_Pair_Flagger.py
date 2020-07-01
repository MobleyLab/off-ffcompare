# This is the dataproliferator script. It takes the alldata.csv produced by
# the dataextractor script and generates useful data from it. User inputs 
# include:
# --tfddc, the TFD difference cutoff after which the script will consider 
# conformers different
# --tfdsc, the TFD similarity cutoff before which the script will consider
# conformers the same 
# --tanisc, the TanimotoCombo similarity cutoff before which the script 
# will consider molecules similar in atom coordinates (and thus eligible
# for flagging) 
# --heavy, the max number of heavy atoms eligible for analysis 
# --directory, the directory that includes the alldata.csv file and where
# new csv's will be written. 
# Written by Jordan Ehrman, from MobleyLab, UCI. 
# Last updated April 16, 2019. 

# Importing
import pandas as pd
import numpy as np
import collections
from itertools import combinations


def read_in(directory,heavyatom):
    """
    This function reads in the alldata.csv. It expects the format of the 
    alldata.csv produced by the dataextractor.py script. It uses the user-
    specified directory and heavy atom count as inputs. 

    Args: 
        directory (str) path that contains csv file 
        heavyatom (int) max number of heavy atoms to consider for analysis

    Returns: 
        alldatadf (dataframe) dataframe of molecules to consider for analysis
    """
    # Reading in csv
    alldatadf = pd.read_csv("%s/alldata.csv" % directory, index_col = 0)
    # Removing errors, which were coded as -1 by dataextractor.py
    for column in alldatadf:
        if alldatadf[column].dtype == float:
            alldatadf = alldatadf[alldatadf[column] >= 0]
    # Removing molecules above some max HeavyAtomCount
    alldatadf = alldatadf[alldatadf["HeavyAtomCount"] <= heavyatom]
    return alldatadf


def print_stat(directory,alldataframe,heavyatom):
    """
    This prints a human-readable statistical summary of each column
    of data. Title of output is dependent on heavyatomcount. 

    Args: 
        directory (str) path to write output to 
        alldataframe (dataframe) dataframe output from read_in function
        heavyatom (int) max heavy atom count for analysis

    Returns: 
        statistics(heavyatom).txt Human readable text file describing 
        statistics of each column of alldataframe. 
    """
    # Creating empty dictionary to populate later
    statdict = {}
    # Finding the percentiles that correspond to each column of dataframe
    for column in alldataframe:
        tempstat = alldataframe[column].describe(percentiles\
        = [.25,.5,.75,.95])
        statdict[column] = tempstat
    # Exporting statistics(heavyatomcount).txt
    with open('%s/statistics%d.txt' % (directory, heavyatom),'w') as data:
        data.write(str(statdict))


def get_ff_combos(alldataframe):
    """This function finds all combinations of forcefields using only the
    titles of the columns of the imported csv file.

    Args: 
        alldataframe (dataframe) dataframe produced by read_in function. 

    Returns: 
        ffs (set) set of all combinations of force fields 
    """
    fflistlist = []
    space = " "
    # Generating combinations of force fields
    for columns in alldataframe:
        fflist = columns.split(' ')[-2:]
        if len(fflist) == 2:
            fflist = fflist[0] + space + fflist[1]
        else:
            fflist = fflist[0]
        fflistlist.append(fflist)
    # Removing objects that are not forcefields from ffs
    fflist_but_not_ffs = ['MolNames','HeavyAtomCount']
    # And removing repeats
    ffs = set(fflistlist) - set(fflist_but_not_ffs)
    return ffs


# Parses through user inputs
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
    
    # Runs above functions
    alldatadf = read_in(directory,heavyatomlimit)
    print_stat(directory,alldatadf,heavyatomlimit)
    ffs = get_ff_combos(alldatadf)
    
    # Flags molecules as similar or different
    # Different molecules are defined different by TFD, same by TaniCo
    # Same molecules defined as similar by TFD
    flaggerdicts = {}
    antiflaggerdicts = {}
    for ffcombo in ffs:
        tempflagdatadf = alldatadf[(alldatadf['TFD %s' % ffcombo]\
        > tfddifcutoff) & (alldatadf['TANI %s' % ffcombo] > tanisamecutoff)]
        flaggerdicts[ffcombo] = tempflagdatadf
        tempantiflagdatadf = alldatadf[(alldatadf['TFD %s' % ffcombo]\
        < tfdsamecutoff)]
        antiflaggerdicts[ffcombo] = tempantiflagdatadf
        
    # Counts same flags and difference flags for each flagged molecule
    # 'Extend' appends list-like column to list
    # 'Counter' counts frequency of object in list 
    flagcountmolname = []
    for i in flaggerdicts:
        flagcountmolname.extend((flaggerdicts[i]["MolNames"]))
    counter = collections.Counter(flagcountmolname)
    antiflagcountmolname = []
    for i in flaggerdicts:
        antiflagcountmolname.extend((antiflaggerdicts[i]["MolNames"]))
    anticounter = collections.Counter(antiflagcountmolname)
    
    # Organizes and exports csvs of data
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
    
    # Finds conformers that are different and exports dataframes of them 
    ffindivlist = list()
    indivdict = {}
    for ff in ffs:
        templist = ff.split(' ')
        ffindivlist.extend(templist)
    ffindivset = set(ffindivlist)
    for forcef in ffindivset:
        columndict = {}
        for columntitle in flagsetdf:
            if forcef in columntitle:
                if 'TFD' in columntitle:
                    firstdf = flagsetdf[flagsetdf[columntitle] > tfddifcutoff]
                elif 'TANI' in columntitle:
                    firstdf = flagsetdf[flagsetdf[columntitle] > tanisamecutoff]
                else:
                    continue
            else:
                if 'TFD' in columntitle:
                    firstdf = flagsetdf[flagsetdf[columntitle] < tfdsamecutoff]
                else: 
                    continue
            columndict[columntitle] = firstdf
        forcedifdf = columndict[list(columndict.keys())[0]].merge\
        (columndict[list(columndict.keys())[1]])
        for i in range(1,len(list(columndict.keys()))):
            forcedifdf = forcedifdf.merge(columndict[list(columndict.keys())\
            [i]])
        indivdict[forcef] = forcedifdf
    for difforce in indivdict:
        indivdict[difforce].to_csv('%s/indivdif%s%d.csv'\
        % (directory,difforce,heavyatomlimit))

