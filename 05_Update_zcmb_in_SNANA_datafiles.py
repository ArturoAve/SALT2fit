#!/usr/bin/env python
# coding: utf-8

# # Code to replace the (zcmb, zhelio) values to the updated values in Andy's compilation in the SNANA photometry data files
#
# - The output files will be saved in a folder called updated_zcmb

# User

# Type of sample. lowz = low-z sample of SNe Ia
# ps1 = panstarrs, des = DES. This definition is simply to include
# or not the prefix 'sn' in the name of the SNe when defining the
# names to be used as entries in the dictionary in the second loop
# of the main loop.
# The lowz sample needs the sufix 'sn' while ps1 and des don't.
sampletype = 'lowz'    # (lowz, ps1, des)

zcmb_update = True  # Update zcmb?
zhelio_update = True  # Update zhelio?

NotebookName = '02_Update_zcmb_in_SNANA_datafiles.ipynb'

#--------------------------------------------------
#     LOW-Z

"""
FoldersList = [
"CFA3_KEPLERCAM",
"CFA4",
"CSPDR2",
"JLA2014_CSP",
"JLA2014_LOWZ_LANDOLT",
# "LOWZ_CFA3", # Bad SNANA file format.
"LOWZ_FromAndyFriedman",
"LOWZ_JRK07",
"PS1s_CFA3_KEPLERCAM_RS14",
"PS1s_CFA4_p1_RS14",
"PS1s_CFA4_p2_RS14",
"PS1s_CSPDR2_V_RS14",
"SNLS3year_CFA3_KEPLERCAM",
"SNLS3year_JRK07" ]

dirSnanaFolders = '/Users/arturo/Dropbox/Research/\
SoftwareResearch/SNANA/Odyssey/home_snana/lowz/2_AndySampleOnly/'
"""

#--------------------------------------------------
#   LOW-Z. INVESTIGATE Y-BAND TREND IN HD

FoldersList = [
"CSPDR2",
"JLA2014_CSP",
"LOWZ_JRK07",
"PS1s_CFA3_KEPLERCAM_RS14",
"PS1s_CFA4_p1_RS14",
"PS1s_CSPDR2_V_RS14" ]

dirSnanaFolders = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/Yband_trend/snana/'

#--------------------------------------------------
#     RAISIN
"""
# FolderSample = 'Data/'
FolderSample = ''

FoldersList = []

# dirSnanaFolders = '/Users/arturo/Dropbox/Research/Articulos/12_RAISINs/\
# Data/RAISIN_2/Data/DES/2017_11_21/Photometry/fit_SomeLinesCommented_ok/'

dirSnanaFolders = '/Users/arturo/Dropbox/Research/Articulos/12_RAISINs/Data/\
RAISIN_2/Data/DES/2017_11_21/Photometry/fit_SomeLinesCommented_ok/\
version_02/Data/1_zcmb_no_updated/'

"""

#--------------------------------------------------
cc = 299792.458  # Speed of light (km/s)

# # Automatic

# #### Function to identify string or number

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

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

# #### Get the name of this ipython notebook
# To print it in the output text files as reference.

get_ipython().run_cell_magic('javascript', '', 'var kernel = IPython.notebook.kernel;\nvar thename = window.document.getElementById("notebook_name").innerHTML;\nvar command = "NotebookName = " + "\'"+thename+".ipynb"+"\'";\nkernel.execute(command);')

print '#', (NotebookName)

# Update_zcmb_in_SNANA_datafiles_v1_0.ipynb

# ### Metadata LOW-Z SAMPLE
#
# I don't need to run these cells if I'm updating the RAISINs SNANA files

import numpy as np

DirMetadata = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/'

# Reading the metadata file
infoSNe_data = np.genfromtxt(DirMetadata+
                             'carrick_Flow_corrections_snnames_v1.txt',
                            dtype=['S17', float,float, 'S40',float,float,
                                   float,float,float,float,'S16', int ])

# Create a dictionary:
# {snname: ra, dec, zhelio, e_zhel, zcmb, e_zcmb, zcmbFlow, e_zcmbFlow, code}

infoSNe_dict = {infoSNe_data['f0'][i]: np.array( [ infoSNe_data['f1'][i],
                infoSNe_data['f2'][i],
                infoSNe_data['f4'][i]/cc, infoSNe_data['f5'][i]/cc,
                infoSNe_data['f6'][i]/cc, infoSNe_data['f7'][i]/cc,
                infoSNe_data['f8'][i]/cc, infoSNe_data['f9'][i]/cc,
                infoSNe_data['f11'][i]] )
                for i in range(len(infoSNe_data)) }

print infoSNe_dict['sn1991T']
# [  1.88542500e+02   2.66556000e+00   5.79067269e-03   3.33564095e-06
#    3.19220839e-03   2.00138457e-05   6.60456908e-03   5.00346143e-04
#    2.00000000e+00]

"""
# (z_helio, error_z_helio)
# From WoodVasey & Andy Friedman metadata file.

DirMetadata = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp/MyNotesAbout/'

InfoSN_zHelio = np.genfromtxt(DirMetadata+'WoodVasey_Andy/WoodVasey_AndyFriedman_zCMB_2017_08_11_Converted_.txt',
                            usecols=(0,1,2), dtype=['S24', float, float])

# Create a dictionary:
InfoSN_zHelio_dict ={}

for i in range(len(InfoSN_zHelio)):
    snname_int1 = InfoSN_zHelio[i][0]
    zHelio_int1 = InfoSN_zHelio[i][1]/cc
    err_zHelio_int1 = InfoSN_zHelio[i][2]/cc

    InfoSN_zHelio_dict[snname_int1] = [ zHelio_int1, err_zHelio_int1  ]

InfoSN_zHelio_dict['sn2006D']
# [0.0085258982732647672, 1.6678204759907602e-05]
"""
0

"""
# (z_helio, error_z_helio, z_CMB, error_z_CMB)
# From WoodVasey & Andy Friedman metadata file for the SPECIAL CASES

InfoSN_zCMB_Special = np.genfromtxt(DirMetadata+'WoodVasey_Andy/WoodVasey_AndyFriedman_zCMB_2017_08_11_SpecialCases_Converted_Ho70.txt',
                            usecols=(0,1,2, 14, 15),
                            dtype=['S10', float, float, float, float])

# Create a dictionary:
InfoSN_zCMB_Special_dict ={}

for i in range(len(InfoSN_zCMB_Special)):
    snname_int2 = InfoSN_zCMB_Special[i][0]
    zHelio_int2 = InfoSN_zCMB_Special[i][1]/cc
    err_zHelio_int2 = InfoSN_zCMB_Special[i][2]/cc
    zCMB_int2 = InfoSN_zCMB_Special[i][3]/cc
    err_zCMB_int2 = InfoSN_zCMB_Special[i][4]/cc

    InfoSN_zCMB_Special_dict[snname_int2] = [ zHelio_int2, err_zHelio_int2,
                                              zCMB_int2, err_zCMB_int2  ]

InfoSN_zCMB_Special_dict['sn1999cl']
# [0.0076085970114698484,
# 1.0006922855944561e-05,
# 0.0031922083910463153,
# 2.0013845711889123e-05]
"""
0

"""
# Using Michael Foley flow-corrected z_CMB + Cepheid distances + special cases
# compiled by Andy Friedman

InfoSN_zCMB_MFoley = np.genfromtxt(DirMetadata+'zCMB_FlowCorrected_MichaelFoley_original.txt',
                            usecols=(0, 1, 2, 5), dtype=['S18', float, float, float])

# Create a final dictionary: (snname: zhelio, err_zhelio, zcmb, err_zcmb, RA, DEC)
InfoSN_dict ={}

for i in range(len(InfoSN_zCMB_MFoley)):
    snname_int3     = InfoSN_zCMB_MFoley[i][0]
    zHelio_int3     = InfoSN_zHelio_dict[snname_int3][0]
    err_zHelio_int3 = InfoSN_zHelio_dict[snname_int3][1]
    RA_int  = InfoSN_zCMB_MFoley[i][1]
    DEC_int = InfoSN_zCMB_MFoley[i][2]


    #---- Define z_CMB and err_z_CMB ------

    if snname_int3 in list(InfoSN_zCMB_Special['f0']): # Special cases
        zCMB_int3 = InfoSN_zCMB_Special_dict[snname_int3][2]
        err_zCMB_int3 = InfoSN_zCMB_Special_dict[snname_int3][3]
    else: # Michael Foley's zcmb values
        zCMB_int3 = InfoSN_zCMB_MFoley[i][3]
        err_zCMB_int3 = 150/cc

    InfoSN_dict[snname_int3] = [ zHelio_int3, err_zHelio_int3,
                                 zCMB_int3, err_zCMB_int3, RA_int, DEC_int ]

InfoSN_dict['sn1998bu']
# [0.0029920699339274241,
# 1.3342563807926082e-05,
# 0.0021141927350000001,
# 0.0005003461427972281,
# 161.69179,
# 11.83531]
"""
0

# ### Metadata RAISINs
#
# Don't need to run these cells if I'm updating the LOW-Z SNANA files.

import numpy as np

# Metadata file

"""
DirMetadata = '/Users/arturo/Dropbox/Research/Articulos/12_RAISINs/\
Metadata/useful_but_tmp/'

InfoSN_z = np.genfromtxt(DirMetadata+'table_zcmb_EBVmw_.txt',
                            dtype=['S15', float, float, float,
                                  float, float, float])
# Create a dictionary:
InfoSN_z_dict ={}

for i in range(len(InfoSN_z)):
    snname_int1 = InfoSN_z[i][0]
    zHelio_int1 = InfoSN_z[i][1]
    e_zhelio_int1 = 0 # It is just a placeholder for now
    zcmb_int1 = InfoSN_z[i][4]
    e_zcmb_int1 = 0 # It is just a placeholder for now
    ebv_int1   = InfoSN_z[i][5]
    e_ebv_int1 = InfoSN_z[i][6]

    InfoSN_z_dict[snname_int1] = [ zHelio_int1, e_zhelio_int1,
                                   zcmb_int1,   e_zcmb_int1,
                                   ebv_int1,    e_ebv_int1]

InfoSN_z_dict['DES15C1nhv']
# [0.42047000000000001,
#  0,
#  0.42011100000000001,
#  0,
#  0.0091999999999999998,
#  0.00020000000000000001]

"""
0

# # Read/replace zcmb, zhelio in SNANA photometry datafile

# ### Main loop

# Read all the photometry file names in a given sample folder:
import os # To use command line like instructions
import glob # To read the files in my directory
# Read the time and date right now
now = datetime.datetime.now()

# Create a log file:
log_file = open(dirSnanaFolders+'Log_Updated_zcmb_files.log','w')

log_file.write("# SNANA file with updated zCMB (written in REDSHIFT_FINAL) based \n\
# on Michael Foley + Cepheid distances + special cases compiled by Andy Friedman. \n")

text_line = '#'+'-'*60+'\n'

log_file.write(text_line)
log_file.write('# Data table created by: Arturo Avelino \n')
text_01 = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
log_file.write('# On date: %s \n'%text_01)
log_file.write('# Script used: %s \n'%NotebookName)
log_file.write(text_line)

#--------------------------------------------------
# Reset
countSN = 0 # Count number of SNe with updated values
countNoAndy = 0 # SNe that are NOT in Andy's sample
FolderSample = ''

# Loop over folder samples
for i in range(len(FoldersList)):

    FolderSample = FoldersList[i]
    print "#-------- %s ---------"%FolderSample

    # Change the working directory where the data files are located
    os.chdir(dirSnanaFolders+FolderSample+'/data/')

    #--- Reading the data files in the folder
    list_SNe = glob.glob('*.DAT')
    if len(list_SNe) == 0: list_SNe = glob.glob('*.dat')

    # Create the folder where the output will be saved
    #- Force the creation of the directory to save the outputs.
    #- "If the subdirectory does not exist then create it"
    DirSaveOutput = dirSnanaFolders+FolderSample+'/updated_zcmb/'
    import os # To use command line like instructions
    if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

    # Create the ".LIST" text file for this folder sample.
    LIST_textfile = open(DirSaveOutput+FolderSample+'.LIST', 'w')
    for file2 in list_SNe:
        # LIST_textfile.write("AndySampleOnly/%s/%s\n"%(FolderSample,file2))
        LIST_textfile.write("Yband_trend/%s/%s\n"%(FolderSample,file2)) # TEMPORAL
    LIST_textfile.close()

    IGNORE_textfile = open(DirSaveOutput+FolderSample+'.IGNORE', 'w')
    IGNORE_textfile.write("\n"); IGNORE_textfile.close()

    README_textfile = open(DirSaveOutput+FolderSample+'.README', 'w')
    README_textfile.write("zCMB, zHELIO updated from Andy Friedman compilation\n");
    README_textfile.close()

    #-----------------------------------

    # Reset
    LIST_list = []

    # Loop 0: over all the SNANA files in a given folder sample.
    for file in list_SNe:

        snanafile = file

        # resetting this variable first.
        snname = ''
        snname_int1 = ''
        snname_int2 = ''

        # Loop over each line of a given SNANA file
        for line_1 in open(dirSnanaFolders+FolderSample+'/data/'+snanafile):

            # Find the line_1 with the SN name and define a temporal variable
            if line_1[:5] == 'SNID:':
                snname_int1 = line_1[6:]

                # Remove the blank spaces before and after the temporal SN name:
                for ss in snname_int1:
                    if ss != ' ': snname_int2 = snname_int2 + ss

                if   sampletype == 'lowz' : snname = 'sn'+snname_int2[:-1]
                elif sampletype == 'ps1' or sampletype == 'des':
                    snname = snname_int2[:-1]

                # print snname

        #-----------------------------------

        if sampletype == 'lowz':

            # Loop 2: # Write down the SNANA file but with the
            # updated (z_CMB, zhelio)

            # Check if this SN is in the metadata file where I have
            # updated (z_CMB, zhelio) values for certain SNe.
            if snname in infoSNe_data['f0']:

                log_file.write('   %s \n'%snanafile)
                print '   %s'%snanafile

                # Open the new text file:
                newSNANAtext = open(DirSaveOutput+snanafile, 'w')

                newSNANAtext.write("# SNANA file with updated zCMB (written in REDSHIFT_FINAL and \
REDSHIFT_CMB) based on \n# Carrick flow model compiled by Andy Friedman. \n")
                # newSNANAtext.write(text_line)
                newSNANAtext.write('# Data table created by: Arturo Avelino \n')
                # text_01 = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
                # newSNANAtext.write('# On date: %s \n'%text_01)
                # newSNANAtext.write('# Script used: %s \n'%NotebookName)
                newSNANAtext.write(text_line)

                # Loop over each line of a given SNANA file
                for line_2 in open(dirSnanaFolders+FolderSample+'/data/'+snanafile):

                    # Find the line_2 with 'REDSHIFT_FINAL' and replace the value of zcmb.
                    if line_2[:14] == 'REDSHIFT_FINAL':

                        if infoSNe_dict[snname][8] == 0:
                            zcmb_new     = infoSNe_dict[snname][6]
                            err_zcmb_new = infoSNe_dict[snname][7]
                        else:
                            zcmb_new     = infoSNe_dict[snname][4]
                            err_zcmb_new = infoSNe_dict[snname][5]

                        # Write the updated zcmb line_2 to the text file
                        text_int1 = 'REDSHIFT_FINAL: %1.5f +- %1.6f (CMB)'%(zcmb_new, err_zcmb_new)
                        newSNANAtext.write(text_int1+'\n')

                        log_file.write('old: '+ line_2)
                        log_file.write('new: '+ text_int1 + '\n')
                        log_file.write(' \n')
                        print 'old:', line_2, 'new:', text_int1

                    elif line_2[:12] == 'REDSHIFT_CMB' and zcmb_update:

                        if infoSNe_dict[snname][8] == 0:
                            zcmb_new     = infoSNe_dict[snname][6]
                            err_zcmb_new = infoSNe_dict[snname][7]
                        else:
                            zcmb_new     = infoSNe_dict[snname][4]
                            err_zcmb_new = infoSNe_dict[snname][5]

                        # Write the updated zcmb line_2 to the text file
                        text_int1 = 'REDSHIFT_CMB: %1.5f +- %1.6f (CMB)'%(zcmb_new, err_zcmb_new)
                        newSNANAtext.write(text_int1+'\n')

                        log_file.write('old: '+ line_2)
                        log_file.write('new: '+ text_int1 + '\n')
                        log_file.write(' \n')
                        print 'old:', line_2, 'new:', text_int1

                    elif line_2[:14] == 'REDSHIFT_HELIO' and zhelio_update:

                        zhel_new     = infoSNe_dict[snname][2]
                        err_zhel_new = infoSNe_dict[snname][3]

                        # Write the updated zhelio line_2 to the text file
                        text_int1 = 'REDSHIFT_HELIO: %1.5f +- %1.6f (HEL)'%(zhel_new, err_zhel_new)
                        newSNANAtext.write(text_int1+'\n')

                        log_file.write('old: '+ line_2)
                        log_file.write('new: '+ text_int1 + '\n')
                        log_file.write(' \n')
                        print 'old:', line_2, 'new:', text_int1

                    # Write all the other "line_2" to the text file
                    else: newSNANAtext.write(line_2)

                newSNANAtext.close()
                countSN = countSN + 1

            else:
                text_log_02 = '   %s: No in Andy s sample. \n'%snanafile
                log_file.write(text_log_02)
                log_file.write(' \n')
                print text_log_02
                countNoAndy = countNoAndy + 1

        #-----------------------------------

        elif sampletype == 'ps1' or sampletype == 'des':

            # Loop 2: # Write down the SNANA file but with the
            # updated (z_CMB, zhelio)

            # Check if this SN is in the metadata file where I have
            # updated (z_CMB, zhelio) values for certain SNe.
            if snname in InfoSN_z['f0']:

                log_file.write('   %s \n'%snanafile)
                print '   %s'%snanafile

                # Open the new text file:
                newSNANAtext = open(DirSaveOutput+snanafile, 'w')

                newSNANAtext.write("# SNANA file with updated zhelio and/or zCMB \
(written in REDSHIFT_FINAL). \n")
                newSNANAtext.write('# Data table created by: Arturo Avelino \n')
                text_01 = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
                newSNANAtext.write('# On date: %s \n'%text_01)
                newSNANAtext.write('# Script used: %s \n'%NotebookName)
                text_line = '#'+'-'*60+'\n'
                newSNANAtext.write(text_line)

                zhel_new     = InfoSN_z_dict[snname][0]
                err_zhel_new = InfoSN_z_dict[snname][1]
                zcmb_new     = InfoSN_z_dict[snname][2]
                err_zcmb_new = InfoSN_z_dict[snname][3]

                # Loop over each line of a given SNANA file
                for line_2 in open(dirSnanaFolders+FolderSample+'/'+snanafile):

                    # Find the line_2 with 'REDSHIFT_xxx' and replace the value of zxxx.
                    # if  zcmb_update == True:
                    if line_2[:15] == 'REDSHIFT_FINAL:':
                        text_int0 = line_2[25:]
                        text_int1 = 'REDSHIFT_FINAL:  %1.5f  %s'%(zcmb_new, text_int0)
                        text_log1 = 'old: '+ line_2
                        text_log2 = 'new: '+ text_int1
                        log_file.write(text_log1); log_file.write(text_log2);
                        print text_log1, text_log2

                    elif line_2[:15] == 'REDSHIFT_HELIO:':
                        text_int0 = line_2[25:]
                        text_int1 = 'REDSHIFT_HELIO:  %1.5f  %s'%(zhel_new, text_int0)
                        text_log1 = 'old: '+ line_2
                        text_log2 = 'new: '+ text_int1
                        log_file.write(text_log1); log_file.write(text_log2);
                        print text_log1, text_log2

                    # Write all the other "line_2's" lines to the text file
                    else: text_int1 = line_2

                    # Write the updated zcmb line_2 to the text file
                    newSNANAtext.write(text_int1)

                newSNANAtext.close()
                countSN = countSN + 1

            else:
                text_log_02 = '   %s: No in Andy s sample. \n'%snanafile
                log_file.write(text_log_02)
                log_file.write(' \n')
                print text_log_02
                countNoAndy = countNoAndy + 1

        #-----------------------------------

    text_10 = '%s SNe were updated their zcmb values peacefully :) \n'%countSN
    text_11 = "%s SNe are not in the metadata file. \n"%countNoAndy

    log_file.write(text_line); log_file.write(text_10); log_file.write(text_11)
    print
    print text_10, text_11

#-------------------
log_file.close()

