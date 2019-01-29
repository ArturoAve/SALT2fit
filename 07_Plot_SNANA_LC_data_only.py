#!/usr/bin/env python
# coding: utf-8

# # Plot the photometry from the SNANA files
#
# This is useful to have a quick inspection of the data and trim the light-curves to the MJD around is the the B-band maximum.

# # User

# Number of header lines in the SNANA files
# 55 for DES-RAISINs2
SkipHeaderRows = 55

#-------------------
# Dir Snana files
DirSNANA = '/Users/arturo/Dropbox/Research/Articulos/12_RAISINs/Data/RAISIN_2/DES/2017_11_21/Photometry/Around_Bmax/Data/test_2/'

DirSaveOutput = DirSNANA

# # Automatic

import numpy as np
from matplotlib import pyplot as plt
import os # To use command line like instructions
import glob # To read the files in my directory

5+4

import os # To use command line like instructions
import glob # To read the files in my directory

# Change the working directory where the data files are located
os.chdir(DirSNANA)

list_files = glob.glob("*.dat")

print "%s SNANA files found."%len(list_files)

countFiles = 0
for snanafile in list_files:

    # snanafile = 'des_real_01285317.dat'
    print snanafile
    photometry = np.genfromtxt(DirSNANA+snanafile,
                            skip_header=SkipHeaderRows,
                            # usecols = (1,2,4,5,6),
                            # dtype=[float, 'S1', float, float, int]
                            usecols = (1,2,4,5),
                            dtype=[float, 'S1', float, float]
                              )

    #----- Plotting -------
    fig = plt.figure()
    # plt.plot(photometry['f0'], photometry['f2'])
    plt.errorbar(photometry['f0'], photometry['f2'], yerr=photometry['f3'],
                 fmt='.')

    plt.title(snanafile)
    plt.xlabel("MJD")
    plt.ylabel("Flux")
    plt.grid(True)
    plt.savefig(DirSaveOutput+snanafile[:-4]+'.png')
    plt.close()

    countFiles = countFiles + 1

print "%s SNANA ligth curves fitted."%countFiles

