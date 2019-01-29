#!/usr/bin/env python
# coding: utf-8

# # Diverse tools

# # Check if the SALT2 covariance matrices for (mB, x1, c) are positive definite and symmetric
# These are required condition for covariance matrices

import numpy as np

#-----------------------------------------------------------------------------80
# USER

datafile = 'GP_subsample_BVrig_No_FUDGEMAGERROR.fitres'

DirData = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/3_GaussianProcessSubsample/cutoffs_no/'

# Automatic

# Numpy array with the following values:
# (SNname, x1ERR, cERR, mBERR, x0ERR, COV_x1_c, COV_x1_x0, COV_c_x0)

data_np = np.genfromtxt(DirData+datafile, comments='#', skip_header=2,
                usecols=(1, 22, 24, 26, 28, 29, 30, 31),
                dtype=['S20', float, float, float, float, float, float, float])

# Check the first rows of this array
data_np[:5]

# Main loop

for i1 in range(len(data_np)):

    # Name of the SNe.
    snname = data_np[i1][0]

    # Construct the covariance matrix
    covmatrix = np.array([[data_np[i1][1]**2, data_np[i1][5],    data_np[i1][6]],
                          [data_np[i1][5],    data_np[i1][2]**2, data_np[i1][7]],
                          [data_np[i1][6],    data_np[i1][7],    data_np[i1][4]**2] ])

    try:
        np.linalg.cholesky(covmatrix)
        # print "# %s: Matrix -is- positive definite."%snname
    except Exception as e:
        ## If the matrix is NOT positive definite, print it:
        print '# %s:'%snname; print(e)

