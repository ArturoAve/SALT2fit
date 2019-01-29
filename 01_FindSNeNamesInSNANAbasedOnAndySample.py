#!/usr/bin/env python
# coding: utf-8

# # SNANA file names comparisons
#
# - Compare 2 columns of names looking for differences.
# - Create a text file with the repeated files.

# # Find and match the names between SNANA database and the SNeIa in my Hubble diagrams
#
#
# I'm given with two table files containing each one a column of SNe names. This script read that column and create a list of those names that are and are not in both lists.
#
# The names in the list may have (or not) prefixes or suffixes; this script is able to read and compare that kind of lists, for instance, in one list the name to compare is "sn2006bh__u_CSP_21_CSP_J" while in the other list is written as "JLA2014_CSP_2006bh.dat": this script is able to output: "Yes, 2006bh is in both lists".
#

# ### USER

import numpy as np
import matplotlib as plt

#--------------------------------------------------
#     Table 1

#      Any YJHK GP Hubble diagram
# DirTables_1 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
# 10Compute/TheTemplates/AllBands/Plots/HubbleDiagram/GaussianProcess/'
# DataFile_1 = 'Table_TotalMu_AllBands_Notes_.txt'

#      Y-band template-method Hubble diagram
DirTables_1 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/Y_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi3_EBVh0.4_Method7_MinData3_vpec150_ok/plots_HD/'
DataFile_1 = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

table_1 = np.genfromtxt(DirTables_1+DataFile_1, usecols=(0), dtype=['S28'])

print "# %s SNe in this list."%len(table_1)

# Ignore the first characters in the name (e.g., the "sn" in "sn2006bo")
AppendNameCharacters_1 = 2   #  2

#--------------------------------------------------
#  Output directory

#      Any YJHK GP Hubble diagram
# DirSaveOutput = '/Users/arturo/Dropbox/Research/SoftwareResearch/\
# SNANA/Odyssey/home_snana/lowz/2_AndySampleOnly/'

#      Y-band template-method Hubble diagram
DirSaveOutput = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/Y_band/snana/'

# Copy the SNANA LC files to the output directory too?
copyLCfiles = True

#--------------------------------------------------
NotebookName = "01_FindSNeNamesInSNANAbasedOnAndySample.ipynb"

# Samples in SNANA database
# /n/panlfs3/aavelino/SNDATA_ROOT/lcmerge/2017_05_12/

#     Table 2

# Location of the ".LIST" files.
DirTables_2 = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/1_lcmerge_LISTS_SNDATA_ROOT_2018_01_02/LIST_files_Originals/'

# Names of the ".LIST" files in SNANA database.
# 1st argument: The name of the ".LIST" file, that it is equal to
# the folder's name in SNANA's "lcmerge" LC folder.
# 2nd argument: The number correponds to the number of characters to ignore before
# the actual SN name in the LC-file name.
# 3rd argument: The subsample associated to that folder, i.e., CfA, CSP, Others.

SamplesSNANA = [
('CFA3_4SHOOTER2', 15, 'CfA'),
('CFA3_KEPLERCAM', 15, 'CfA'),
('CFA4', 5, 'CfA'),
('CSPDR2', 7, 'CSP'),
('ESSENCE_WV07', 13, 'Oth'),
('HST_Riess06_GOLD', 12, 'Oth'),
('HST_Riess06', 12, 'Oth'),
('JLA2014_CfAIII_4SHOOTER2', 25, 'CfA'),
('JLA2014_CfAIII_KEPLERCAM', 25, 'CfA'),
('JLA2014_CSP', 12, 'CSP'),
('JLA2014_HST', 12, 'Oth'),
('JLA2014_LOWZ_LANDOLT_a', 20, 'Oth'),
('JLA2014_LOWZ_LANDOLT_a', 20, 'CSP'),
('JLA2014_LOWZ_LANDOLT_b', 13, 'CfA'),
('JLA2014_LOWZ_LANDOLT_c', 14, 'CfA'),
('JLA2014_LOWZ_LANDOLT_d', 18, 'Oth'),
('JLA2014_SNLS', 13, 'Oth'),
# ('LOWZ_CFA3', 10, 'CfA'), # SNANA files with wrong format.
('LOWZ_JRK07', 11, 'CfA'),
('LOWZ_JRK07', 11, 'CSP'),
('LOWZ_JRK07', 11, 'Oth'),
('PS1s_CFA1_JRK07_RS14', 11, 'CfA'),
('PS1s_CFA2_JRK07_RS14', 11, 'CfA'),
('PS1s_CFA3_4SHOOTER2_RS14', 15, 'CfA'),
('PS1s_CFA3_KEPLERCAM_RS14', 15, 'CfA'),
('PS1s_CFA4_p1_RS14', 8, 'CfA'),
('PS1s_CFA4_p2_RS14', 8, 'CfA'),
('PS1s_CSPDR2_V_RS14', 7, 'CSP'),
('PS1s_PS1_RS14', 0, 'Oth'),
# ('SDSS_allCandidates+BOSS', 0, 'CfA'),
('SDSS_HOLTZ08', 13, 'Oth'),
('SNLS_Ast06', 11, 'Oth'),
('SNLS3year_CFA3_4SHOOTER2', 25, 'CfA'),
('SNLS3year_CFA3_KEPLERCAM', 25, 'CfA'),

# Version SNDATA_ROOT_2018_01_02 don't have any LC file.
# ('SNLS3year_CSP', 7, 'CSP'),

('SNLS3year_JRK07', 16, 'CfA'),
('SNLS3year_JRK07', 16, 'Oth'),
('SNLS3year_MEGACAM', 18, 'Oth'),
('SNLS3year+METADATA', 19, 'Oth')
]

#--------------------------------------------------
# If "copyLCfiles = True" then look for the SNANA LC files in the
# following directory (usually the "lcmerge" dir in SNANA):

if copyLCfiles:
    DirSNANAlcmerge = '/Users/arturo/Documents/Research/SNANA/SNDATA_ROOT_2018_05_04/lcmerge/'

# ## Automatic

# Function to identify if a string is an integer number or a letter.
# This will be used in the dictionary construction to properly read some SN names.

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

# Tests
print is_number('5'), is_number('e')
# True False

# #### Get the name of this ipython notebook
# To print it in the output text files as reference

get_ipython().run_cell_magic('javascript', '', 'var kernel = IPython.notebook.kernel;\nvar thename = window.document.getElementById("notebook_name").innerHTML;\nvar command = "NotebookName = " + "\'"+thename+".ipynb"+"\'";\nkernel.execute(command);')

print '#', (NotebookName)
# Update_zcmb_in_SNANA_datafiles_v1_0.ipynb

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

# #### Main loop
# Find the elements of table 1 that ARE in table 2

# Find the elements of table 1 that ARE in table 2

from shutil import copyfile # To copy/paste data files.

File_snana = open(DirSaveOutput+'SNe_inGPHD_AndySample_snanaNamesOnly.txt','w')
File_MyFormat = open(DirSaveOutput+'SNe_inGPHD_AndySample.txt','w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_line = '#'+'-'*50 + '\n'

File_MyFormat.write('# Data table created by: Arturo Avelino \n')
File_MyFormat.write('# On date: %s \n'%text_timenow)
File_MyFormat.write('# Script used: %s \n'%NotebookName)
File_MyFormat.write(text_line)

for k in range(len(SamplesSNANA)):

    FolderSample = SamplesSNANA[k][0]
    DataFile_2 = FolderSample+'.LIST'
    # print DataFile_2
    table_2 = np.genfromtxt(DirTables_2+DataFile_2, usecols=(0), dtype=['S50'])

    # Ignore the first characters in the name (e.g., the "sn" in "sn2006bo")
    AppendNameCharacters_2 = SamplesSNANA[k][1]
    SuffixNameCharacters_2 = -4

    # Identify the sample of table 2 according to the labels used in table 1.
    # With this I can select those SN names that are in both tables AND that correspond
    # to the same subsample.
    sampleTable_2 = SamplesSNANA[k][2]  # (CfA, CSP, Oth). "Oth" stands for "Others"
    # print sampleTable_2

    #------------------------

    text_01 = '#'+'-'*50
    text_02 = '# The following names are common in both tables.'
    text_03 = '# table 1: %s '%DataFile_1
    text_04 = '# table 2: %s '%DataFile_2
    text_05 = '# Sample type of table 2:'+sampleTable_2

    print text_01
    print text_02
    print text_03
    print text_04
    print text_05

    File_snana.write(text_01+'\n'); File_snana.write(text_05+'\n');
    File_snana.write(text_02+'\n'); File_snana.write(text_03+'\n');
    File_snana.write(text_04+'\n');
    File_MyFormat.write(text_01+'\n'); File_MyFormat.write(text_05+'\n');
    File_MyFormat.write(text_02+'\n'); File_MyFormat.write(text_03+'\n');
    File_MyFormat.write(text_04+'\n');

    #------------------------

    count = 0 # Reset counter

    for i in range(len(table_1)):
        name_1_int = table_1['f0'][i]
        if   name_1_int[7] == '_': name_1 = name_1_int[AppendNameCharacters_1:7] # To read correctly, "sn2011B"
        elif name_1_int[7] != '_':
            if is_number(name_1_int[7]): name_1 = name_1_int[AppendNameCharacters_1:15] # To read correctly, "snf20080514-002"
            else: name_1 = name_1_int[AppendNameCharacters_1:8]  # To read correctly, "sn1998bu"

        sampleFlat_1 = name_1_int[19:22]

        for j in range(len(table_2)):

            name_2_int = table_2['f0'][j]
            name_2 = name_2_int[AppendNameCharacters_2:SuffixNameCharacters_2]   # Remove the appendix and suffix in the name
            # print name_2

            # Write down the SN's names
            if name_1 == name_2 and sampleFlat_1 == sampleTable_2:
                File_MyFormat.write('%s    %s   %s \n'%(name_1_int, name_2_int, FolderSample))
                File_snana.write('%s \n'%name_2_int)
                print name_1_int, '|', name_2_int
                # print name_2_int

                if copyLCfiles:

                    # First give the actual name for'JLA2014_LOWZ_LANDOLT' folder.
                    if (FolderSample=='JLA2014_LOWZ_LANDOLT_a' or
                        FolderSample=='JLA2014_LOWZ_LANDOLT_b' or
                        FolderSample=='JLA2014_LOWZ_LANDOLT_c' or
                        FolderSample=='JLA2014_LOWZ_LANDOLT_d'):
                        FolderSample = 'JLA2014_LOWZ_LANDOLT'

                    #- Force the creation of the directory to save the outputs.
                    #- "If the subdirectory does not exist then create it"
                    import os # To use command line like instructions
                    if not os.path.exists(DirSaveOutput+FolderSample+'/data/'):
                        os.makedirs(DirSaveOutput+FolderSample+'/data/')

                    # Copy file from Dir to DirDestination directory.
                    copyfile(DirSNANAlcmerge+FolderSample+'/'+name_2_int,
                             DirSaveOutput+FolderSample+'/data/'+name_2_int)

                count = count + 1

    text_10 = '# %s names are equal in both tables.'%count
    text_empty = '#'
    print text_10
    print text_empty
    File_snana.write(text_10+'\n'); File_MyFormat.write(text_10+'\n')
    File_snana.write(text_empty+'\n'); File_MyFormat.write(text_empty+'\n')

File_snana.close()
File_MyFormat.close()

# print ''
print '# All done.'

# # Find repeated names in a list
#
# Find the repeated names in a given list

# USER

# DirList = '/Users/arturo/Dropbox/Research/SoftwareResearch/\
# SNANA/Odyssey/panlfs3_aavelino/lcmerge/0_lcmerge_LISTS_2018_05_04/'

DirList = DirSaveOutput

# DirSaveOutput = DirList

FileName = 'SNe_inGPHD_AndySample.txt'

#-------------------
# Preferred sample folder to consider to use that photometric data.
# This is ordered based on the number of SNe I selected manually from each folder:
# the first folder's name corresponds to the one with more SNe, and so on.
PreferredFolders = ['PS1s_CSPDR2_V_RS14', 'LOWZ_JRK07', 'PS1s_CFA3_KEPLERCAM_RS14',
                    'JLA2014_LOWZ_LANDOLT', 'PS1s_CFA4_p1_RS14','PS1s_CFA4_p2_RS14',
                    'CSPDR2', 'JLA2014_CSP', 'CFA3_KEPLERCAM',
                    'SNLS3year_JRK07', 'SNLS3year_CFA3_KEPLERCAM', 'CFA4'
                   ]

# #### Automatic

import numpy as np
import matplotlib as plt
#-------------------
FileData = np.genfromtxt(DirList+FileName, usecols=(0, 1, 2), dtype=['S32', 'S45', 'S45'])
print '# FileData lenght = %s rows'%len(FileData)

# FileData lenght = 99 rows

# ##### Main loop

# Find repeated names in a given list

file_1 = open(DirSaveOutput+FileName[:-4]+'_Repeated_.txt', 'w')
file_1.write('#    Repeated \n')
file_1.write('# SN names that are repeated AT LEAST one time, or are unique. \n')
file_1.write("# NOTE: I've selected automatically the SNe photometric data files that \n\
# I want to finally use, and their corresponding sample folder. \n")

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_line = '#'+'-'*50 + '\n'

file_1.write(text_line)
file_1.write('# Data table created by: Arturo Avelino \n')
file_1.write('# On date: %s \n'%text_timenow)
file_1.write('# Script used: %s \n'%NotebookName)
file_1.write(text_line)

file_1.write('# SN names                      SNANA names                                \
Folder sample               Notes  \n')
# file_1.write('# \n')

ii = 0; count1 = 0; count2 = 0;
listRepeatedNames = []
listRepeatedNames_int = []
listUniqueNames = []
listNamesAlreadyWritten =[]

# print '# Repeated names:'

for i in range(len(FileData)):
    name_a = FileData['f0'][i]
    nameSnana_a  = FileData['f1'][i]
    folder_a = FileData['f2'][i]

    # List to store the different folders and LC files in SNANA
    # about the SNe "name_a".
    listRepeated_OneName = []
    listRepeated_OneName += [(name_a, nameSnana_a, folder_a)]

    # Reset:
    flag_repeated = 0
    ii = i + 1

    #--------------------------------------
    for j in range(len(FileData)-ii):

        name_b     = FileData['f0'][j+ii]
        nameSnana  = FileData['f1'][j+ii]
        folder_b   = FileData['f2'][j+ii]

        # Find repeated SNe. It could detect several times the same SNe.
        if name_b == name_a:
            flag_repeated = 1 # Flag that this name is repeated.
            count1 = count1 + 1
            listRepeatedNames += [(name_b, nameSnana, folder_b)]
            listRepeated_OneName += [(name_b, nameSnana, folder_b)]

            # Write down just one time the name of each SNe:
            if name_b not in listRepeatedNames_int:
                listRepeatedNames_int += [name_b]

        # Put in a list the SNe that were not repeated, i.e., are unique.
        # Check at the end of loop with indez "j".
        elif (j == (len(FileData)-ii-1) and flag_repeated==0 and
              name_a not in listRepeatedNames_int):
            print "%s: no repeated."%name_a
            listUniqueNames += [(name_a, nameSnana_a, folder_a)]

    #--------------------------------------

    # Put in a list the SNe that were not repeated, i.e., are unique.
    # Check at the end of loop with indez "i".
    if (i == (len(FileData)-1) and
          name_a not in listRepeatedNames_int):
        print "%s: no repeated."%name_a
        listUniqueNames += [(name_a, nameSnana_a, folder_a)]

    #--------------------------------------

    # Write in the text file the SNE from my preferred SNANA folders only.

    flag_PreferredFolders = 0 # reset.

    for j2 in range(len(PreferredFolders)):

        for j3 in range(len(listRepeated_OneName)):

            if (listRepeated_OneName[j3][2] == PreferredFolders[j2]
                and flag_PreferredFolders == 0
                and name_a not in listNamesAlreadyWritten):

                file_1.write('%-29s  %-42s  %-35s # \n'%(
                    listRepeated_OneName[j3][0], listRepeated_OneName[j3][1],
                    listRepeated_OneName[j3][2] ))

                print listRepeated_OneName[j3][0]
                # Put the name in the list of SNe already written in the text file,
                # it implies that I've already chosen a preferred SNANA folder to pick
                # the LC file.
                listNamesAlreadyWritten += [listRepeated_OneName[j3][0]]

                # flag that I'm choosing a preferred folder for this SN.
                flag_PreferredFolders = 1

                count2 = count2 + 1

                break

    #--------------------------------------

text_10 = '# Number of names repeated at least one time = %s names. \n'%count2
text_11 = '# Approximated number of extra times some names are repeated =< %s times. \n'%count1

file_1.write(text_10); file_1.write(text_11); file_1.write(text_line);

print text_10
print text_11

#-------------------
# Write the SNe that were not repeated

file_1.write('#   SNe that were not repeated but must be in Hubble diagram \n');
file_1.write('# \n');

if len(listUniqueNames) > 0:
    for i3 in range(len(listUniqueNames)):
        file_1.write('%-29s  %-42s  %-35s # \n'%(listUniqueNames[i3][0],
                        listUniqueNames[i3][1], listUniqueNames[i3][2]  ))

file_1.write(text_line);

#-------------------
# Write a scratch list with all the repeated events.
# Print and write to a text file the names that are repeated at least 1 time.

file_1.write('# \n');
file_1.write('#        Just in case I need this \n');
text_14 = '# %s times I found repeated SNe names. They are: \n'%len(listRepeatedNames)
file_1.write(text_14);
# print text_14

# Sort the list based on the sample-folder name (i.e., the third column in the list).
sortedList = sorted(listRepeatedNames, key=lambda folderName: folderName[2])

for i2 in range(len(sortedList)):
    file_1.write('# %-29s  %-42s  %-35s  # \n'%(sortedList[i2][0],
                                                sortedList[i2][1],
                                                sortedList[i2][2]))
    # print listRepeatedNames[i2]


#-------------------
file_1.close()

# # Find SNe in my list that are NOT in SNANA database

# User

import numpy as np
import matplotlib as plt

#--------------------------------------------------
# DirTables_1 = "/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
# 10Compute/TheTemplates/AllBands/Plots/HubbleDiagram/GaussianProcess/"

# DataFile_1 = 'Table_TotalMu_AllBands_Notes_.txt'
# table_1 = np.genfromtxt(DirTables_1+DataFile_1, usecols=(0), dtype=['S28'])

# Ignore the first characters in the name (e.g., the "sn" in "sn2006bo")
# AppendNameCharacters_1 = 0

#-----------------------------------------------------------------------------80
# DirTables_2 = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/\
# Odyssey/panlfs3_aavelino/lcmerge/0_lcmerge_LISTS_2018_05_04/'

DirTables_2 = DirSaveOutput

#DataFile_2 = 'SNe_inGPHD_AndySample.txt'
DataFile_2 = FileName

table_2 = np.genfromtxt(DirTables_2+DataFile_2, usecols=(0), dtype=['S30'])

# Ignore the first characters in the name (e.g., the "sn" in "sn2006bo")
AppendNameCharacters_2 = 0

#--------------------------------------------------
print '# Size of table_1:', len(table_1), len(table_1.T)
print '# Size of table_2:', len(table_2), len(table_2.T)
# Size of table_1: 59 59
# Size of table_2: 120 120

print ''
print "# Checking that I'm reading table 2 correctly:"
print '#     ', table_2['f0'][0][AppendNameCharacters_2:]
print '#     ', table_2['f0'][1][AppendNameCharacters_2:]
print '#     ', table_2['f0'][2][AppendNameCharacters_2:]

# sn2005cf__U_36_B_1_CfA_J
# sn2005el__U_55_B_3_CfA_J
# sn2005eq__U_22_B_2_CfA_J

# ### Automatic

# Find the elements of table 1 that are NOT in table 2

count = 0

print '#', '#'*50
print '# '
print '#      Missing SNe in SNANA database'
print '# '
print '# The following names in table 1 are NOT in table 2:'
print '# table 1: %s '%DataFile_1
print '# table 2: %s '%DataFile_2
print '# Dir table 1: %s'%DirTables_1
print '# Dir table 2: %s'%DirTables_2

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_line = '#'+'-'*50

print text_line
print '# Data table created by: Arturo Avelino'
print '# On date: %s '%text_timenow
print '# Script used: %s '%NotebookName
print text_line

for i in range(len(table_1)):
    for j in range(len(table_2)):
        if table_1['f0'][i] == table_2['f0'][j]:
            # print '-'
            break
    else:
        if j==(len(table_2)-1):
            print "# %s"%table_1['f0'][i]
            count = count + 1

# print ''
print '#', count, 'names in table 1 are NOT in table 2.'
print '#', '#'*50

