#!/usr/bin/env python
# coding: utf-8

# # Code to convert from Andy's "wstd" LC files to SNANA format
#
# ## Notes
# - The output files will be saved in a subfolder called 'snana'

import numpy as np
import glob # To read the files in my directory
import os # To use command line like instructions
import datetime # Get the current date and time

#--------------------------------------------------------60
code_created_by = 'Arturo_Avelino'
# On date: 2017.01.10 (yyyy.mm.dd)
code_name = '03_WstdAndy_to_SNANA.ipynb'
code_version = '0.1.7'
code_last_update = '2019.05.10'
code_location = 'https://github.com/ArturoAvelino/SALT2fit/blob/master/03_WstdAndy_to_SNANA.py'

##############################################################################80

# # User

# Directory where Andy's wstd files to be converted are located:

dirwstd = '/Users/arturo/Downloads/tmp/Research/wstd_snana/Wstd2/Others/input_Wstd/'

DirSaveOutput = dirwstd+'snana/'

# String to be printed in the 'SURVEY' line in SNANA files.
# NOTE: For low-z CSP data set "Survey = 'CSP' ". This allows to fit the data in
# SNANA/SALT2 with no issues.
# Options:
#      'LOWZ' for low-z CfA and Others LCs.
#      'CSP' for low-z CSP LCs.
Survey = 'LOWZ'

#--------------------------------------------------------60
# Filter's name matching between Andy and SNANA

# Write -all- the filter's name in the line "FILTERS: " in the
# SNANA-format LC file? If False, then it will be written the
# filters listed in the list 'filterListToConvert' below. If
# true, then it will write down -all- the filters names. For those
# that are not in the list 'filterListToConvert' then it will
# write the filter's name using Andy's names.
write_all_filters = True

# List of filter names that will be converted their names from Andy's
# to SNANA convention's names.
# Any other filter not listed here will not be written in the
# line "FILTERS: " in the SNANA-like output file,
# unless "write_all_filters = True".
filterListToConvert = ['r_prime', 'i_prime',
                       'B_CTIO1p3m','V_CTIO1p3m','R_CTIO1p3m','I_CTIO1p3m',
                       'u_CSP', 'g_CSP', 'r_CSP', 'i_CSP', 'B_CSP', 'V_CSP',
                        'B', 'V', 'R', 'I']

# Create a dictionary with the conversion names from above
# The dictionary structure is: (filter_Andy: filter_SNANA)
FilterNameConversion_dict = {}
FilterNameConversion_dict['r_prime'] = ['r']
FilterNameConversion_dict['i_prime'] = ['i']
FilterNameConversion_dict['B_CTIO1p3m'] = ['B']
FilterNameConversion_dict['V_CTIO1p3m'] = ['V']
FilterNameConversion_dict['R_CTIO1p3m'] = ['R']
FilterNameConversion_dict['I_CTIO1p3m'] = ['I']
FilterNameConversion_dict['u_CSP'] = ['u']
FilterNameConversion_dict['g_CSP'] = ['g']
FilterNameConversion_dict['r_CSP'] = ['r']
FilterNameConversion_dict['i_CSP'] = ['i']
FilterNameConversion_dict['B_CSP'] = ['B']
FilterNameConversion_dict['V_CSP'] = ['o']

# To verify if the SNANA name is the correct one:
FilterNameConversion_dict['B'] = ['B']
FilterNameConversion_dict['V'] = ['V']
FilterNameConversion_dict['R'] = ['R']
FilterNameConversion_dict['I'] = ['I']

#--------------------------------------------------------60
# Rename the output files using the first 'TrimFileName' characters
# of the datafile.

# "-4" = use the full name of the file except the extension characters ".dat"
TrimFileName = -4
# TrimFileName = 18

# Print on the output files the date-time and script used to
# create them?
print_date_scriptName = True

cc = 299792.458  # Speed of light (km/s)

##############################################################################80

# # Automatic

#- Force the creation of the directory to save the outputs.
#- "If the subdirectory does not exist then create it"
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

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
# Read the time and date now
now = datetime.datetime.now()

# ### Metadata

MetadataFile = 'carrick_Flow_corrections_snnames_v1.txt'
DirMetadata = '/Users/arturo/Downloads/tmp/Research/wstd_snana/Wstd2/'

# Reading the metadata file
infoSNe_data = np.genfromtxt(DirMetadata+MetadataFile,
                            dtype=['S17', float,float, 'S40',float,float,
                                   float,float,float,float,'S16', int ])

# Create a dictionary:
# {snname: zhelio, e_zhel, zcmb, e_zcmb, zcmbFlow,
#          e_zcmbFlow, code, ra, dec}

InfoSN_dict = {infoSNe_data['f0'][i]: [
                infoSNe_data['f4'][i]/cc, infoSNe_data['f5'][i]/cc,
                infoSNe_data['f6'][i]/cc, infoSNe_data['f7'][i]/cc,
                infoSNe_data['f8'][i]/cc, infoSNe_data['f9'][i]/cc,
                infoSNe_data['f11'][i],
                infoSNe_data['f1'][i], infoSNe_data['f2'][i] ]
                for i in range(len(infoSNe_data)) }

InfoSN_dict['sn1998bu']
# [0.0029620491653595902,
#  3.3356409519815205e-06,
#  0.0023683050759068795,
#  8.6726664751519533e-05,
#  0.003138838135814611,
#  0.00050034614279722807,
#  1,
#  161.69166999999999,
#  11.835279999999999]

#-----------------------------------------------------------------------------80

# Metadata file with (t_Bmax, EBV_MW)

MetadataFile2 = 'LowzSNe_metadata.txt'
DirMetadata2 = '/Users/arturo/Downloads/tmp/Research/wstd_snana/Wstd2/'

DistMu_np = np.genfromtxt(DirMetadata2+ MetadataFile2,
                             dtype=['S15',float,float,float,float])

#----- Create a dictionary -----
# (snname: TBmax, err_TBmax, EBV_MW, err_EBV_MW)

DistMu_dict ={}
for i in range(len(DistMu_np)):

    # Sn name
    snName_1 = DistMu_np['f0'][i]

    TBmax_int_1 = DistMu_np['f1'][i]

    # add 53000 MJD to tBmax of some CSP SNe.
    if TBmax_int_1 < 4000:
        TBmax_int = TBmax_int_1 + 53000
    else: TBmax_int = TBmax_int_1

    err_TBmax_int = DistMu_np['f2'][i]

    EBV_MW_int     = DistMu_np['f3'][i]
    err_EBV_MW_int = DistMu_np['f4'][i]

    DistMu_dict[snName_1] = [TBmax_int, err_TBmax_int, EBV_MW_int, err_EBV_MW_int]

# print '#', DistMu_dict['sn1998bu']

# #### Convert fluxes to zeropoint = 27.5 (default in SNANA)

def flux_snana(flux_old, zp_Andy):

    zp_snana = 27.5
    flux_new = flux_old * 10**(0.4*(zp_snana - zp_Andy))

    return flux_new

print '# Test:', flux_snana(587.41, 25)
# Test: 5874.1

#-----------------------------------------------------------------------------80

# # Read/convert Andy's wstd file to SNANA-format text file

# Read the names of all the photometry files to be converted

os.chdir(dirwstd)

#- Reading the LC data file names
the_list = glob.glob('*.Wstd.dat')

print '# %s SNe in this list'%len(the_list)

# #### Main loop

zp_Andy = 25  # Andy's zeropoint: 25 mag

count_lcs = 0 # count LCs converted to SNANA format

# Loop over Andy's wstd LC files.
for wstdFile in the_list:

    wstd_np = np.genfromtxt(dirwstd+wstdFile, usecols=[1,2,3,4,5] ,
                        dtype=[float,'S14', float, float, float])

    #====================================================

    # Create the variable "snName" containing the first 8 (or 7)
    # letters of the SNe file name
    try:
        if   wstdFile[7] == '_':
            snName = wstdFile[:7]  # To read correctly, e.g., "sn2011B_"
        elif wstdFile[7] != '_':
            # To read correctly, e.g., "snf20080514-002"
            if is_number(wstdFile[7]): snName = wstdFile[:15]
            else: snName = wstdFile[:8]  # To read correctly, e.g., "sn1998bu"
    except: snName = wstdFile[:6]  # To read correctly, e.g., "sn2011B"

    #-------------------------------
    # Create a list of filters that are in a given photometric file:

    # Reset
    ListFilters = []
    filtername_int_2 = ''; filtersnana_2 = '';

    # Loop over the photometry:
    for j in range(len(wstd_np)):

        # Read the filter's name
        filtername_int_2 = wstd_np['f1'][j]

        # Convert the filter's name to SNANA
        if filtername_int_2 in filterListToConvert:
            filtersnana_2 = FilterNameConversion_dict[filtername_int_2][0]
        else:
            if write_all_filters:
                filtersnana_2 = filtername_int_2

        # Create a list of unique filters
        if filtersnana_2 not in ListFilters:
            ListFilters += [filtersnana_2]

    # print '# Filters in this file: ',ListFilters

    #---- Create a single string with the name of all the filters ----
    # This will be written in the SNANA-format text file
    # in the row "FILTERS: "

    ListFiltersToPrint = ''
    for name in ListFilters:
        ListFiltersToPrint = ListFiltersToPrint+name

    # print '# Text to print in the field FILTERS:', ListFiltersToPrint

    #--------------------------------------
    # Determine the number of observations based on the filters
    # to be considered

    NOBS = 0

    # Loop over the photometry:
    for j2 in range(len(wstd_np)):
        # Read the filter's name
        filtername_int_2 = wstd_np['f1'][j2]

        if write_all_filters:
            NOBS += 1
        elif filtername_int_2 in filterListToConvert:
            NOBS += 1

    #====================================================

    # Read the time and date right now
    now = datetime.datetime.now()

    snana_file = open(DirSaveOutput+wstdFile[:TrimFileName]+'_snana.dat', 'w')

    text_line_1 = '#'+'-'*60+'\n'
    snana_file.write('#    %s \n'%snName)
    snana_file.write("# Andy Friedman's wstd file converted to SNANA-like format. \n")
    snana_file.write('# Source file: %s \n'%wstdFile)
    snana_file.write("# Metadata information from: %s \n"%MetadataFile)
    snana_file.write("# Located at: \n")
    snana_file.write("# %s \n"%DirMetadata)
    snana_file.write("# write_all_filters = %s. \n"%write_all_filters)

    snana_file.write(text_line_1)
    snana_file.write('# Data table created by: %s\n'%code_created_by)
    if print_date_scriptName:
        text_01 = now.strftime("%Y.%m.%d (yyyy.mm.dd); %H:%M hrs (ET).")
        snana_file.write('# On date: %s \n'%text_01)
        snana_file.write('# Script used: %s \n'%code_name)
        snana_file.write('# Script version: %s \n'%code_version)
        snana_file.write('# Script location:\n')
        snana_file.write('# %s\n'%code_location)
    snana_file.write(text_line_1)

    #---------------------

    #   SNANA header

    # Create a final dictionary:
    # (snname: zhelio, err_zhelio, zcmb, err_zcmb)

    zhel = InfoSN_dict[snName][0]
    err_zhel = InfoSN_dict[snName][1]

    # Flag to determine the appropiate z_cmb:
    flag_zcmb  = InfoSN_dict[snName][6]
    if flag_zcmb > 0.1:
        zcmb     = InfoSN_dict[snName][2]
        err_zcmb = InfoSN_dict[snName][3]
    else:
        zcmb     = InfoSN_dict[snName][4]
        err_zcmb = InfoSN_dict[snName][5]

    RA       = InfoSN_dict[snName][7]
    DEC      = InfoSN_dict[snName][8]

    EBV_MW  = DistMu_dict[snName][2]
    PEAKMJD = DistMu_dict[snName][0]

    snana_file.write('SURVEY: %s \n'%Survey)
    snana_file.write('SNID: %s \n'%snName[2:])
    snana_file.write('IAUC: %s \n'%snName[2:])
    snana_file.write('RA: %s deg \n'%RA)
    snana_file.write('DECL: %s deg \n'%DEC)
    snana_file.write('MWEBV: %s MW E(B-V) \n'%EBV_MW)
    snana_file.write('REDSHIFT_HELIO: %1.5f +- %1.6f (HEL)\n'%(zhel, err_zhel))
    snana_file.write('REDSHIFT_CMB: %1.5f +- %1.6f (CMB)\n'%(zcmb, err_zcmb))
    snana_file.write('REDSHIFT_FINAL: %1.5f +- %1.6f (CMB)\n'%(zcmb, err_zcmb))
    snana_file.write('SEARCH_PEAKMJD: %.3f \n'%PEAKMJD)
    snana_file.write('FILTERS: %s \n'%ListFiltersToPrint)

    #---------------------

    snana_file.write(text_line_1)
    snana_file.write('NOBS: %s \n'%NOBS)
    snana_file.write('NVAR: 7 \n')
    snana_file.write('VARLIST:  MJD  FLT           FIELD       FLUXCAL         FLUXCALERR  MAG  MAGERR  \n')

    #---------------------

    # Loop over photometry:
    for i in range(len(wstd_np)):

        filtername_int1 = wstd_np['f1'][i]

        if filtername_int1 in filterListToConvert:
            filtersnana = FilterNameConversion_dict[filtername_int1][0]
        else: filtersnana = filtername_int1

        flux_old         = wstd_np['f2'][i]
        err_flux_old_low = wstd_np['f3'][i]
        err_flux_old_hig = wstd_np['f4'][i]

        flux_new         = flux_snana(flux_old, zp_Andy)
        err_flux_new_low = flux_snana(err_flux_old_low, zp_Andy)
        err_flux_new_hig = flux_snana(err_flux_old_hig, zp_Andy)

        average_errorFlux = (err_flux_new_low + err_flux_new_hig)/2

        # Write the line in the text file
        if write_all_filters:
            snana_file.write('OBS: %.3f  %-12s NULL  %15.4f %15.4f    0     0 \n'%(wstd_np['f0'][i],
                        filtersnana, flux_new, average_errorFlux))

        elif filtername_int1 in filterListToConvert:
            snana_file.write('OBS: %.3f  %-12s NULL  %15.4f %15.4f    0     0 \n'%(wstd_np['f0'][i],
                        filtersnana, flux_new, average_errorFlux))

    snana_file.write('END:')
    snana_file.close();

    count_lcs += 1
    print wstdFile[0:40]

#--------------------------------------------------------60
print '\n# All %s LCs converted smoothly'%count_lcs

snana_file.close();snana_file.close();snana_file.close();

