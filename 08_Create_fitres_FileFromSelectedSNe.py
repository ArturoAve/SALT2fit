#!/usr/bin/env python
# coding: utf-8

# # Create a SNANA ".fitres" file from selected SNe
#
# The list then will be used as input in "SNANA_SALT2_DistanceMu_v0_1.ipynb" to compute the
# distance moduli, but over all, the RESIDUAL distance moduli from this particular subsample.

import numpy as np

# To read arguments in command line
# Used in the ".py" version of this notebook.
import sys

# ## User

#-----------------------------------------------------------------------------80

## Terminal or notebook version of this script?
ScriptVersion = 'notebook' # ( terminal , notebook )

#-------------------
if ScriptVersion == 'notebook':

    ## Dir where are the SNANA -folders- containing the *.FITRES files
    DirSNANAFittedFolders = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/2_AndySampleOnly/BVr_FUDGEMAGERROR_0_005_NEW_done/FUDGEALL_ITER1_MAXFRAC_0_20/'

    # The suffix I'll suffix to the output file.
    suffixOutfile = ''

    # Write also the CSV version of the file?
    # This very useful to easily open the file in any spreadsheet app.
    out_csv = False

    # SNANA version used to compute?
    # This is useful because the output ".fitres" file is different
    # between different SNANA versions.
    snana_version = 'v10_58g'     # ('v10_58g', 'v10_35g')

#-------------------
elif ScriptVersion == 'terminal':

    ## Dir where are the SNANA -folders- containing the *.FITRES files
    DirSNANAFittedFolders = sys.argv[1]

    # The suffix I'll suffix to the output file.
    suffixOutfile = sys.argv[2]

    # Write also the CSV version of the file?
    # This very useful to easily open the file in any spreadsheet app.
    out_csv = sys.argv[3] == 'True'

    # SNANA version used to compute?
    # This is useful because the output ".fitres" file is different
    # between different SNANA versions.
    snana_version = sys.argv[4]     # ('v10_58g', 'v10_35g')

#-------------------
# Directory where are located the list "SNe_inGPHD_AndySample_Repeated_Notes_.txt"
# with the names of the SNe that I want to create the ".fitres" file from.

DirLists = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/'

ListSNeToBeUsed = 'SNe_inGPHD_AndySample_Repeated_Notes_.txt'
ListSNe = np.genfromtxt(DirLists+ListSNeToBeUsed, dtype=['S30', 'S43', 'S40'])

#######################################################

#    CUTOFFS

# Apply additional cutoffs to the sample?:
ApplyCutoffs = False

if ApplyCutoffs == False:
    c_limits = -0.3, 0.3   # color
    x1_limits = -3, 3   # light-curve parameter
    x1ERR_limits = 1    # error in the light-curve parameter
    PKMJDERR_limits = 2  # error in the determination of MJD at T_Bmax
    FITPROB_limits = 0.001   # probability of the fit

#######################################################

NotebookName= '08_Create_fitres_FileFromSelectedSNe.ipynb'

#--------------------------------------------------
print '%s SNe in the list (regardless the cutoffs).'%len(ListSNe)
# 50 SNe in the list.

# Dir Save Output

DirSaveOutput = DirSNANAFittedFolders

"""
if ApplyCutoffs == True:
    DirSaveOutput = DirLists+'3_GaussianProcessSubsample/cutoffs_yes/'
elif ApplyCutoffs == False:
    DirSaveOutput = DirLists+'3_GaussianProcessSubsample/cutoffs_no/'
"""

#- Force the creation of the directory to save the outputs.
#- "If the subdirectory does not exist then create it"
import os # To use command line like instructions
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

# ## Useful functions

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
if ScriptVersion == 'notebook':
    print is_number('5'), is_number('e')
    # True False

# #### Get the name of this ipython notebook
# To print it in the output text files as reference

# %%javascript
# var kernel = IPython.notebook.kernel;
# var thename = window.document.getElementById("notebook_name").innerHTML;
# var command = "NotebookName = " + "'"+thename+".ipynb"+"'";
# kernel.execute(command);

if ScriptVersion == 'notebook': print '#', (NotebookName)
# Update_zcmb_in_SNANA_datafiles_v1_0.ipynb

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

# ## Automatic

# #### Main loop

#-----------------------------------------------------------------------------80
file_1 = open(DirSaveOutput+'GP_subsample_%s.fitres'%suffixOutfile,'w')
if out_csv:
    file_2 = open(DirSaveOutput+'GP_subsample_%s.csv'%suffixOutfile,'w')

# Write down the SNANA headers:

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_line = '#'+'-'*50 + '\n'

# These two lines MUST be at the first ones on the text file in order to be
# easily readed by "SNANA_SALT2_DistanceMu_v1_1.ipynb"
if snana_version == 'v10_35g':
    file_1.write('NVAR:  25 \n')
    file_1.write('VARNAMES: CID       z        zERR     x0           x0ERR        c      cERR     x1      x1ERR   PKMJD   PKMJDERR mB      mBERR    COVx0x1     COVx0c      COVx1c     CHI2   NDOF  FITPROB      SNRMAX1    SNRMAX2    SNRMAX3 IDSURVEY TYPE    Subsample \n')

elif snana_version == 'v10_58g':
    file_1.write('NVAR:  35 \n')
    file_1.write('VARNAMES: CID IDSURVEY TYPE FIELD CUTFLAG_SNANA   zHEL   zHELERR     zCMB     zCMBERR       zHD     zHDERR   VPEC VPECERR HOST_LOGMASS HOST_LOGMASS_ERR   SNRMAX1     SNRMAX2      SNRMAX3    PKMJD    PKMJDERR     x1            x1ERR              c         cERR         mB      mBERR           x0           x0ERR       COV_x1_c      COV_x1_x0    COV_c_x0     NDOF    FITCHI2      FITPROB   Subsample \n')

    if out_csv:
        file_2.write('VARNAMES:, CID, IDSURVEY, TYPE, FIELD, CUTFLAG_SNANA,   zHEL,   zHELERR,     zCMB,     zCMBERR,       zHD,     zHDERR,   VPEC, VPECERR, HOST_LOGMASS, HOST_LOGMASS_ERR,   SNRMAX1,     SNRMAX2,      SNRMAX3,    PKMJD,    PKMJDERR,     x1,            x1ERR,              c,         cERR,         mB,      mBERR,           x0,           x0ERR,       COV_x1_c,      COV_x1_x0,    COV_c_x0,     NDOF,    FITCHI2,      FITPROB,   Subsample,\n')


file_1.write(text_line)

file_1.write('# SNANA-like "fitres" file created from selected SNe from the list: \n')
file_1.write('# %s \n'%ListSNeToBeUsed)
file_1.write('# %s = SNANA version used to for the computations. \n'%snana_version)

file_1.write('# Data table created by: Arturo Avelino \n')
file_1.write('# On date: %s \n'%text_timenow)
file_1.write('# Script used: %s \n'%NotebookName)
file_1.write(text_line)

if ApplyCutoffs == True:
    file_1.write("# Cutoff applied: %s < c < %s | %s < x1 < %s | x1ERR < %s,  \n"%
                 (c_limits[0], c_limits[1], x1_limits[0], x1_limits[1], x1ERR_limits) )
    file_1.write("# PKMJDERR < %s | FITPROB > %s. \n"%(PKMJDERR_limits, FITPROB_limits))
    file_1.write(text_line)

#-----------------------------------------------------------------------------80

# Reset
countSN = 0; countSNeNoPassCuts=0; count_byhand = 0
flag_ok = 0

for i in range(len(ListSNe)):
    name          = ListSNe['f0'][i]
    namefile        = ListSNe['f1'][i]
    snSnanaFolder = ListSNe['f2'][i]

    # Determine the subsample (CfA, CSP, Others) and create a
    # flag accordingly.
    subsample_append = name[-3:]
    flag_subsample = 0 # Reset
    if   subsample_append == 'CfA': flag_subsample = 1
    elif subsample_append == 'CSP': flag_subsample = 2
    elif subsample_append == 'ers': flag_subsample = 3

    # print "%s | %s "%(name, snSnanaFolder)

    # Read correctly the name of the SNe.
    if   name[7] == '_': snName = name[2:7] # To read correctly, e.g., "sn2011B"
    elif name[7] != '_':
        # To read correctly, e.g., "snf20080514-002"
        if is_number(name[7]): snName = name[2:15]
        else: snName = name[2:8]  # To read correctly, e.g., "sn1998bu"

    if ScriptVersion == 'notebook':
        print "%s | %-10s | %s "%(name, snName, snSnanaFolder)

    Dir_int1 = DirSNANAFittedFolders+snSnanaFolder+'/'

    if snana_version == 'v10_35g':
        fitresFile = np.genfromtxt(Dir_int1+snSnanaFolder+'.fitres',
                    skip_header=2, dtype=['S3', 'S15',
                      float, float, float, float, float, float, float,
                      float, float, float, float, float, float, float, float,
                      float, float, float, float, float, float, float, float ] )
    elif snana_version == 'v10_58g':
        fitresFile = np.genfromtxt(Dir_int1+snSnanaFolder+'.FITRES.TEXT',
                    skip_header=8, dtype=['S3', 'S15', float, float, 'S4',
                                  float, float, float, float, float,
                                  float, float, float, float, float,
                                  float, float, float, float, float,
                                  float, float, float, float, float,
                                  float, float, float, float, float,
                                  float, float, float, float, float])

    # Find the index where the SN is located in the fitres file.
    index_int = np.where(fitresFile['f1'] == snName)[0][0]

    # When there is only one SNe in the fitres files then index_int = 0. However,
    # this produces an error in the following part of the code because for an array
    # of dimension 1 it is not needed to indicate the row in the array! To prevent
    # the error, I do not include the fitres files with only 1 SN, then I have to
    # copy/paste its information by hand.
    try:
        len(fitresFile['f1'])
        # flag "Ok" if I don't get any error when "len(fitresFile['f1'])". This
        # means that the fitres files has at least 2 rows or more.
        flag_ok = 1
    except:
        print "%s | %s. The info in the fitres file for this SN has to be copied/pasted by hand :("%(namefile, snSnanaFolder)
        flag_ok = 0
        count_byhand += 1

    #-----------

    if snana_version == 'v10_35g' and flag_ok:
        c_par  = fitresFile['f6'][index_int]
        x1_par = fitresFile['f8'][index_int]
        x1ERR_par = fitresFile['f9'][index_int]
        PKMJDERR_par = fitresFile['f11'][index_int]
        FITPROB_par  = fitresFile['f19'][index_int]

    elif snana_version == 'v10_58g' and flag_ok:
        c_par  = fitresFile['f23'][index_int]
        x1_par = fitresFile['f21'][index_int]
        x1ERR_par = fitresFile['f22'][index_int]
        PKMJDERR_par = fitresFile['f20'][index_int]
        FITPROB_par  = fitresFile['f34'][index_int]

    if ApplyCutoffs == True  and flag_ok: # SNe that pass the cutoffs:
        if (c_par > c_limits[0] and c_par < c_limits[1] and
            x1_par > x1_limits[0] and x1_par < x1_limits[1] and
            x1ERR_par < x1ERR_limits and PKMJDERR_par < PKMJDERR_limits and
            FITPROB_par > FITPROB_limits):

            snName_print = 'SN: %-13s'%fitresFile['f1'][index_int]
            countSN = countSN + 1

        else:
            snName_print = '## SN: %-13s'%fitresFile['f1'][index_int]
            countSNeNoPassCuts += 1

    elif flag_ok:
        snName_print = 'SN: %-13s'%fitresFile['f1'][index_int]
        countSN = countSN + 1
        if out_csv:
            snName_print2 = 'SN:, %-13s'%fitresFile['f1'][index_int]

    if snana_version == 'v10_35g' and flag_ok:
    # Write down a line in the text list with all the fitres information for a given SN.
        file_1.write('%s  %.5f  %.5f  %.5e  %.3e  %7.4f  %.4f  %7.4f  %.4f  %.3f  %.2f %.4f  %.4f  %10.3e  %10.3e  %10.3e  %5.1f  %4.0f  %.3e  %9.3f  %9.3f  %9.3f  %2.0f     %2.0f %.0f \n'%(
        snName_print, fitresFile['f2'][index_int], fitresFile['f3'][index_int],
        fitresFile['f4'][index_int], fitresFile['f5'][index_int], fitresFile['f6'][index_int],
        fitresFile['f7'][index_int], fitresFile['f8'][index_int], fitresFile['f9'][index_int],
        fitresFile['f10'][index_int], fitresFile['f11'][index_int],
        fitresFile['f12'][index_int], fitresFile['f13'][index_int],
        fitresFile['f14'][index_int], fitresFile['f15'][index_int], fitresFile['f16'][index_int],
        fitresFile['f17'][index_int], fitresFile['f18'][index_int], fitresFile['f19'][index_int],
        fitresFile['f20'][index_int], fitresFile['f21'][index_int], fitresFile['f22'][index_int],
        fitresFile['f23'][index_int], fitresFile['f24'][index_int], flag_subsample  ))

    if snana_version == 'v10_58g' and flag_ok:
    # Write down a line in the text list with all the fitres information for a given SN.
        file_1.write('%s  %3.0f  %3.0f  %s     %3.0f     %9.5f  %.6f  %9.5f  %.6f  %9.5f  %.6f   %3.0f    %3.0f    %10.5f    %10.5f    %10.4f  %10.4f  %10.4f  %10.3f  %6.3f  %12.4e   %12.4e  %12.4e  %12.4e %9.4f  %10.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %4.0f  %10.4f  %12.4e   %.0f \n'%(
        snName_print, fitresFile['f2'][index_int], fitresFile['f3'][index_int],
        fitresFile['f4'][index_int], fitresFile['f5'][index_int], fitresFile['f6'][index_int],
        fitresFile['f7'][index_int], fitresFile['f8'][index_int], fitresFile['f9'][index_int],
        fitresFile['f10'][index_int], fitresFile['f11'][index_int], fitresFile['f12'][index_int],
        fitresFile['f13'][index_int], fitresFile['f14'][index_int], fitresFile['f15'][index_int],
        fitresFile['f16'][index_int], fitresFile['f17'][index_int], fitresFile['f18'][index_int],
        fitresFile['f19'][index_int], fitresFile['f20'][index_int], fitresFile['f21'][index_int],
        fitresFile['f22'][index_int], fitresFile['f23'][index_int], fitresFile['f24'][index_int],
        fitresFile['f25'][index_int], fitresFile['f26'][index_int], fitresFile['f27'][index_int],
        fitresFile['f28'][index_int], fitresFile['f29'][index_int], fitresFile['f30'][index_int],
        fitresFile['f31'][index_int], fitresFile['f32'][index_int], fitresFile['f33'][index_int],
        fitresFile['f34'][index_int], flag_subsample))

        # Copy/paste of above but with separated by commas.
        if out_csv:
            file_2.write('%s  ,%3.0f  ,%3.0f  ,%s     ,%3.0f     ,%9.5f  ,%.6f  ,%9.5f  ,%.6f  ,%9.5f  ,%.6f   ,%3.0f    ,%3.0f    ,%10.5f    ,%10.5f    ,%10.4f  ,%10.4f  ,%10.4f  ,%10.3f  ,%6.3f  ,%12.4e   ,%12.4e  ,%12.4e  ,%12.4e ,%9.4f  ,%10.4e  ,%12.4e  ,%12.4e  ,%12.4e  ,%12.4e  ,%12.4e  ,%4.0f  ,%10.4f  ,%12.4e   ,%.0f, \n'%(
        snName_print2, fitresFile['f2'][index_int], fitresFile['f3'][index_int],
        fitresFile['f4'][index_int], fitresFile['f5'][index_int], fitresFile['f6'][index_int],
        fitresFile['f7'][index_int], fitresFile['f8'][index_int], fitresFile['f9'][index_int],
        fitresFile['f10'][index_int], fitresFile['f11'][index_int], fitresFile['f12'][index_int],
        fitresFile['f13'][index_int], fitresFile['f14'][index_int], fitresFile['f15'][index_int],
        fitresFile['f16'][index_int], fitresFile['f17'][index_int], fitresFile['f18'][index_int],
        fitresFile['f19'][index_int], fitresFile['f20'][index_int], fitresFile['f21'][index_int],
        fitresFile['f22'][index_int], fitresFile['f23'][index_int], fitresFile['f24'][index_int],
        fitresFile['f25'][index_int], fitresFile['f26'][index_int], fitresFile['f27'][index_int],
        fitresFile['f28'][index_int], fitresFile['f29'][index_int], fitresFile['f30'][index_int],
        fitresFile['f31'][index_int], fitresFile['f32'][index_int], fitresFile['f33'][index_int],
        fitresFile['f34'][index_int], flag_subsample))

file_1.write(text_line)
text_10 = '# %s SNe Ia passed the cutoffs (if applied). \n'%countSN
text_11 = "# %s SNe Ia didn't pass the cutoffs (##). \n"%countSNeNoPassCuts
text_12 = "# %s SNe Ia I need to copy/paste their info from the fitres files by hand. \n"%count_byhand
file_1.write(text_10); file_1.write(text_11); file_1.write(text_12)

file_1.close()
if out_csv: file_2.close()

print text_line, text_10, text_11, text_12
#-----------------------------------------------------------------------------80

file_1.close();file_1.close();file_1.close();
file_1.close();file_1.close();file_1.close();
if out_csv:
    file_2.close();file_2.close();file_2.close();

print '# All done smoothly.'

