#!/usr/bin/env python
# coding: utf-8

# # Code to convert from Andy's "wstd" LC files to SNANA format
#
# ## Notes
# - The output files will be saved in a subfolder called 'snana'

# User

# Directory where Andy's wstd file to be converted is located:

dirwstd = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/2_AndySampleOnly/LOWZ_FromAndyFriedman/Andy2SNANA_format/wstd_originals/'

DirSaveOutput = dirwstd+'snana/'

# String to be printed in the 'SURVEY' line in SNANA files.
# NOTE: For low-z CSP data set "Survey = 'CSP' ". This allows to fit the data in
# SNANA/SALT2 with no issues.
# Options:
#      'LOWZ' for low-z CfA LCs.
#      'CSP' for low-z CSP LCs.
Survey = 'LOWZ'

#-------------------
# Filter's name matching between Andy and SNANA

# Write -all- the filter's name in the line "FILTERS: " in the
# SNANA-format LC file? If False, then it will be written the
# filters listed in the list 'filterListToConvert' below. If
# true, then it will write down -all- the filters names, for those
# that are not in the list 'filterListToConvert' then it will
# write the filter's name using Andy's names.
write_all_filters = False

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

#--------------------------------------------------
# Rename the output files using the first 'TrimFileName' characters
# of the datafile.
# "-4" = use the full name of the file except the extension characters ".dat"
TrimFileName = -4
# TrimFileName = 18

# Print on the output files the date-time and script used to
# create them?
print_date_scriptName = False

NotebookName = '03_WstdAndy_to_SNANA.ipynb'

cc = 299792.458  # Speed of light (km/s)

# # Automatic

import numpy as np
#- Force the creation of the directory to save the outputs.
#- "If the subdirectory does not exist then create it"
import os # To use command line like instructions
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

5+6

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

# # WstdAndy_to_SNANA_v1_1.ipynb

# ### Metadata

# (z_helio, error_z_helio)
# From WoodVasey & Andy Friedman metadata file.

"""
import numpy as np

DirMetadata = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/\
AndyLCComp/MyNotesAbout/'

InfoSN_zHelio = np.genfromtxt(DirMetadata+
                'WoodVasey_Andy/WoodVasey_AndyFriedman_zCMB_2017_08_11_Converted_.txt',
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

# (z_helio, error_z_helio, z_CMB, error_z_CMB)
# From WoodVasey & Andy Friedman metadata file for the SPECIAL CASES
"""
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

# Using Michael Foley flow-corrected z_CMB + Cepheid distances + special cases
# compiled by Andy Friedman
"""
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

DirMetadata = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/'

MetadataFile = 'carrick_Flow_corrections_snnames_v1.txt'

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

# #### DistanceMu_All_BeforeCutoffs.txt
# #### (for the J-band Gaussian-process Hubble diagram)
#
# Reading the 'DistanceMu_All_BeforeCutoffs.txt' from GP Hubble diagram for J band: it contains *almost* all the SNe that I need. I use the information in this file to retrieve (t_Bmax, EBV_MW)

DirJband = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/J_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi_1e6_EBVh0.4_Method7_MinData3_vpec150_ok/plots_HD/'

# DirJband = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
# 10Compute/TheTemplates/AllBands/Plots/HubbleDiagram/GaussianProcess/'

DistMu_np = np.genfromtxt(DirJband+
                          'DistanceMu_All_BeforeCutoffs_.txt',
                             dtype=['S30',
       float,float,float,float,float,float,float,float,float,float,
       float,float,float,float,float,float,float,float,float,float,
       float,float,float,float,float,float,float,float,float,float,
       float,float,float])

#----- Create a dictionary -----
# (snname: TBmax, err_TBmax, EBV_MW, err_EBV_MW)

DistMu_dict ={}
for i in range(len(DistMu_np)):

    # Sn name
    # Create the variable "snName" containing the first 8 (or 7)
    # letters of the SNe file name
    snname_int_1 = DistMu_np['f0'][i]
    try:
        if   snname_int_1[7] == '_':
            snName_1 = snname_int_1[:7]  # To read correctly, e.g., "sn2011B_"
        elif snname_int_1[7] != '_':
            # To read correctly, e.g., "snf20080514-002"
            if is_number(snname_int_1[7]): snName_1 = snname_int_1[:15]
            else: snName_1 = snname_int_1[:8]  # To read correctly, e.g., "sn1998bu"
    except: snName_1 = snname_int_1[:6]  # To read correctly, e.g., "sn2011B"

    TBmax_int     = DistMu_np['f14'][i]
    err_TBmax_int = DistMu_np['f15'][i]
    EBV_MW_int     = DistMu_np['f23'][i]
    err_EBV_MW_int = DistMu_np['f24'][i]

    DistMu_dict[snName_1] = [TBmax_int, err_TBmax_int, EBV_MW_int, err_EBV_MW_int]

print '#', DistMu_dict['sn1998bu']

# [50953.113988, 0.081145, 0.0217, 0.0002]

# #### Convert fluxes to zeropoint = 27.5 (default in SNANA)

def flux_snana(flux_old, zp_Andy):

    zp_snana = 27.5
    flux_new = flux_old * 10**(0.4*(zp_snana - zp_Andy))

    return flux_new

print '# Test:', flux_snana(587.41, 25)

# Test: 5874.1

# # Read/convert Andy's wstd file to SNANA-format text file

# Read the names of all the photometry files to be converted

import glob # To read the files in my directory
import os # To use command line like instructions

os.chdir(dirwstd)

#- Reading the LC data file names
the_list = glob.glob('*.Wstd.dat')

print '# %s SNe in this list'%len(the_list)

the_list

# #### Main loop

zp_Andy = 25  # Andy's zeropoint: 25 mag

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
    snana_file.write('# Data table created by: Arturo Avelino \n')
    if print_date_scriptName:
        text_01 = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
        snana_file.write('# On date: %s \n'%text_01)
        snana_file.write('# Script used: %s \n'%NotebookName)
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

    print wstdFile[0:40]

print '\n# All the conversion done smoothly'

