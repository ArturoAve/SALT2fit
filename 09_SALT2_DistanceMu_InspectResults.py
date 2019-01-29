#!/usr/bin/env python
# coding: utf-8

# # Inspection: Compare uncertainties and some other values output from SALT2 fit

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad as intquad
# To read arguments in command line
# Used in the ".py" version of this notebook.
import sys

#-----------------------------------------------------------------------------80
# USER

## Terminal or notebook version of this script?
ScriptVersion = 'terminal' # ( terminal , notebook )

#-----------------------------------------------------------------------------80

if ScriptVersion == 'notebook':

    ## Name of the FITRES file
    FitresFile = 'GP_subsample_.fitres'

    ## Define the directory where the input FITRES file is located
    DirSALT2Fit = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/2_AndySampleOnly/BVr_mix_fudge/'

    # Peculiar velocity uncertainty:
    sigma_vPec = 150  # 150 km/s

    ## This corresponds to the SNANA version used
    ## during the computations. This is needed because different versions
    ## of SNANA produce different format of the output fitres files.
    ## Options: ( v10_58g , pantheon )
    snana_version = 'v10_58g'

    # Global parameters in the Tripp formula for: low-z.
    ## low-z values from Table 6 of Pantheon
    alpha = 0.147
    beta =  3.00
    sigmaInt = 0.11

#--------------------------------------------------
elif ScriptVersion == 'terminal':

    ## Name of the FITRES file
    FitresFile = sys.argv[1]

    ## Define the directory where the input FITRES file is located
    DirSALT2Fit = sys.argv[2]

    # Peculiar velocity uncertainty:
    sigma_vPec = int(sys.argv[3]) # km/s

    ## This corresponds to the SNANA version used
    ## during the computations. This is needed because different versions
    ## of SNANA produce different format of the output fitres files.
    ## Options: ( v10_58g , pantheon )
    snana_version = sys.argv[4]

    # Global parameters in the Tripp formula for: low-z.
    alpha = float(sys.argv[5])
    beta =  float(sys.argv[6])
    sigmaInt = float(sys.argv[7])

#-----------------------------------------------------------------------------80

# Define the directory where I'll save the output.
DirSaveOutput = DirSALT2Fit

## The suffix I'll suffix to the output file.
# suffixOutfile = '0_80'
suffixOutfile = ''

NotebookName = '09_SALT2_DistanceMu_InspectResults.ipynb'

#--------------------------------------------------
# The distance modulus is determined from the SNANA-SALT2 fit and the global
# parameters reported in Table 6 of Scolnic et al. 2017
# (https://arxiv.org/abs/1710.00845).
# Global parameters in the Tripp formula for: low-z.

MoFix = 19.36 # Absolute magnitude at T_Bmax.

cc = 299792.458  # Speed of light (km/s)

# ### Automatic

# Get the name of this ipython notebook.
# To print it in the output text files as reference

# %%javascript
# var kernel = IPython.notebook.kernel;
# var thename = window.document.getElementById("notebook_name").innerHTML;
# var command = "NotebookName = " + "'"+thename+".ipynb"+"'";
# kernel.execute(command);

if ScriptVersion == 'notebook':
    print '#', (NotebookName)
    # Update_zcmb_in_SNANA_datafiles_v1_0.ipynb

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

# ### SALT2 to mu function

#-----------------------------------------------------------------------------80
# Define the function to convert from SALT2 parameters to distance mu

def salt2mu(z=None, x0=None, c=None, cerr=None, x1=None, x1err=None,
             mb=None, mberr=None, cov_x1_x0=None, cov_c_x0=None, cov_x1_c=None,
            alpha=None, beta=None, sigint=None, Mo=None, sigma_vPec=None):

    # The distance modulus mean value:
    mu_out = mb + x1*alpha - beta*c + MoFix

    #  Computing the distance-modulus uncertainty (Tripp formula):
    sf = -2.5/(x0*np.log(10.0))
    cov_mb_c = cov_c_x0*sf
    cov_mb_x1 = cov_x1_x0*sf

    # Propagation of uncertainty from Tripp formula and accountig for correlations:
    sigma2_mu_Tripp = (mberr**2. + (alpha**2.)*(x1err**2.) + (beta**2.)*(cerr**2.) +
                         2.0*alpha*cov_mb_x1 - 2.0*beta*cov_mb_c -
                         2.0*alpha*beta*cov_x1_c )

    # Uncertainty in distance modulus due to peculiar velocity uncertainty:
    zerr = (sigma_vPec/cc)*5.0/np.log(10)*(1.0+z)/(z*(1.0+z/2.0))

    # Total distance-modulus uncertainty:
    muerr_out = np.sqrt(sigma2_mu_Tripp + zerr**2. + (0.055**2.)*(z**2.) +
                        sigint**2. )

    return (mu_out, muerr_out, np.sqrt(sigma2_mu_Tripp),
            mberr**2. + (alpha**2.)*(x1err**2.) + (beta**2.)*(cerr**2.),
            2.0*alpha*cov_mb_x1, -2.0*beta*cov_mb_c, -2.0*alpha*beta*cov_x1_c)

# Read the input FITRES file

if snana_version == 'v10_58g':
    snanaSalt2Fit_np = np.genfromtxt(DirSALT2Fit+FitresFile, skip_header=8,
                        dtype=['S3', 'S15', float, float, 'S4',
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float ])

elif snana_version == 'pantheon':
    snanaSalt2Fit_np = np.genfromtxt(DirSALT2Fit+FitresFile, skip_header=0,
                        dtype=['S3', 'S15', float, float, float,
                              'S4',  float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float])

if ScriptVersion == 'notebook':
    print '# %s SNe in this file.'%(len(snanaSalt2Fit_np))
    print
    print '# Checking the first entries of the numpy array:'
    print snanaSalt2Fit_np[0]
    print
    print snanaSalt2Fit_np['f1'][0:8]

# Create some arrays and compute the distance modulus

if snana_version == 'v10_58g':
    z_np =  snanaSalt2Fit_np['f10']
    x0_np = snanaSalt2Fit_np['f27']
    c_np = snanaSalt2Fit_np['f23']
    cerr_np = snanaSalt2Fit_np['f24']
    x1_np = snanaSalt2Fit_np['f21']
    x1err_np = snanaSalt2Fit_np['f22']
    mb_np =  snanaSalt2Fit_np['f25']
    mberr_np =  snanaSalt2Fit_np['f26']
    cov_x1_c_np =  snanaSalt2Fit_np['f29']
    cov_x1_x0_np =  snanaSalt2Fit_np['f30']
    cov_c_x0_np =  snanaSalt2Fit_np['f31']

elif snana_version == 'pantheon':
    z_np =  snanaSalt2Fit_np['f9'] #
    x0_np = snanaSalt2Fit_np['f26'] #
    c_np = snanaSalt2Fit_np['f22'] #
    cerr_np = snanaSalt2Fit_np['f23'] #
    x1_np = snanaSalt2Fit_np['f20'] #
    x1err_np = snanaSalt2Fit_np['f21'] #
    mb_np =  snanaSalt2Fit_np['f24'] #
    mberr_np =  snanaSalt2Fit_np['f25'] #
    cov_x1_c_np =  snanaSalt2Fit_np['f28']
    cov_x1_x0_np =  snanaSalt2Fit_np['f29']
    cov_c_x0_np =  snanaSalt2Fit_np['f30']

#-----------------------------------------------------------------------------80
# Computing the distance modulus

mu_SALT2 = salt2mu(z_np, x0_np, c_np, cerr_np, x1_np, x1err_np,
                   mb_np, mberr_np, cov_x1_x0_np, cov_c_x0_np, cov_x1_c_np,
                  alpha, beta, sigmaInt, MoFix, sigma_vPec)

# ## Compare the uncertainties
# Write a table with the comparisons

txtline = '#%s\n'%('-'*79)

#-------------------
# Define the text for some variables

# Total uncertainty on distance modulus as suggested by David Jones
err_mu_total_txt = 'smu_total'

# Uncertainty on distance modulus due to fitting error only.
# This is what is defined as 'sigma_N'
err_mu_photo_txt  = 'smu_Tripp'
err2_mu_photo_txt = 'smu_Tripp^2'

# Uncertainty on the apparent magnitude at t_Bmax.
err_mB_txt = 'mberr'

#-----------------------------------------------------------------------------80
#           Header text

textfile_1 = open(DirSaveOutput+'Table_inspect_%s.txt'%(suffixOutfile),'w')

textfile_1.write('#        Inspection of some uncertainties\n')
textfile_1.write('# Fitres file: %s\n'%FitresFile)
textfile_1.write('# Located at: %s\n'%DirSALT2Fit)
textfile_1.write('# SNANA version: %s\n'%snana_version)
textfile_1.write(txtline)

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Author = '# Data table created by: Arturo Avelino \n'
text_Date   = '# On date: %s \n'%text_timenow
text_script = '# Script used: %s \n'%NotebookName
textfile_1.write(text_Author); textfile_1.write(text_Date);
textfile_1.write(text_script);
textfile_1.write(txtline)

textfile_1.write('# %s ==  sqrt( mberr^2 + (alpha^2)*(x1err^2) + \
(beta^2)*(cerr^2) + 2.0*alpha*cov_mb_x1 - 2.0*beta*cov_mb_c -\
2.0*alpha*beta*cov_x1_c )\n'%(err_mu_photo_txt))
textfile_1.write('# mberr == uncertainty in the SALT2 apparent magnitude at B max.\n')
textfile_1.write('# diff_1 == %s - %s\n'%(err_mu_photo_txt,err_mB_txt))
textfile_1.write('# cov_1 ==  2.0*alpha*cov_mb_x1\n')
textfile_1.write('# cov_2 == -2.0*beta*cov_mb_c\n')
textfile_1.write('# cov_3 == -2.0*alpha*beta*cov_x1_c\n')
textfile_1.write('# cov123 == cov_1 + cov_2 + cov_3\n')
textfile_1.write('# noCov ==  mberr^2 + (alpha^2)*(x1err^2) + \
(beta^2)*(cerr^2)\n')
textfile_1.write('# sqrt(noCov) ==  squareRoot(noCov)\n')
textfile_1.write('# noCov+cov123 == noCov + cov123\n')
textfile_1.write('#\n# Assumed values of the global parameters in the Tripp formula:\n')
textfile_1.write('# alpha = %s | beta = %s | sigmaInt = %s | \
M (abs mag) = %s\n'%(alpha, beta, sigmaInt, MoFix))
textfile_1.write('# Peculiar velocity = %s km/s\n'%sigma_vPec)
textfile_1.write(txtline)
# textfile_1.write('# \n')

textfile_1.write('# SN name,     NumID,  %s, %s, %s,    %s^2,    \
diff_1,     cov123,      noCov,   sqrt(noCov), noCov+cov123,  \
sqrt(noCov+cov123),   cov_1,      cov_2,      cov_3,  sample,  \n'%(
    err_mu_photo_txt, err2_mu_photo_txt, err_mB_txt, err_mB_txt ))

#-----------------------------------------------------------------------------80
# Write each row of the output table

countSN = 0
for i in range(len(mu_SALT2[0])):

    if snana_version == 'v10_58g':
        snname = snanaSalt2Fit_np['f1'][i]
        col = snanaSalt2Fit_np['f23'][i]
        cerr = snanaSalt2Fit_np['f24'][i]
        x1 = snanaSalt2Fit_np['f21'][i]
        x1err = snanaSalt2Fit_np['f22'][i]
        mb =  snanaSalt2Fit_np['f25'][i]
        mberr =  snanaSalt2Fit_np['f26'][i]
        x0 = snanaSalt2Fit_np['f27'][i]
        sampleflag = snanaSalt2Fit_np['f35'][i]

    elif snana_version == 'pantheon':
        snname = snanaSalt2Fit_np['f1'][i]
        col = snanaSalt2Fit_np['f22'][i]
        cerr = snanaSalt2Fit_np['f23'][i]
        x1 = snanaSalt2Fit_np['f20'][i]
        x1err = snanaSalt2Fit_np['f21'][i]
        mb =  snanaSalt2Fit_np['f24'][i]
        mberr =  snanaSalt2Fit_np['f25'][i]
        x0 = snanaSalt2Fit_np['f26'][i]
        sampleflag = snanaSalt2Fit_np['f4'][i]

    diff_1 = mu_SALT2[2][i]-mberr

    # if diff_1 < 0: # To write down the suspicious SNe only.
    if 1 < 2 :
        countSN += 1
        textfile_1.write('%-14s, %2i,     %.5f,   %.5f,   %.5f,   %.5f,   %8.5f,   %8.5f,   %8.5f,   %8.5f,   %11.5f,    %12.5f,       %8.5f,   %8.5f,   %8.5f, %i, \n'%(
        snname,  countSN,
        mu_SALT2[2][i], mu_SALT2[2][i]**2.0, mberr, mberr**2.,
        diff_1,
        mu_SALT2[4][i] +mu_SALT2[5][i] +mu_SALT2[6][i],
        mu_SALT2[3][i], np.sqrt(mu_SALT2[3][i]),
        mu_SALT2[3][i]+(mu_SALT2[4][i] +mu_SALT2[5][i] +mu_SALT2[6][i]),
        np.sqrt(mu_SALT2[3][i]+(mu_SALT2[4][i] +mu_SALT2[5][i] +mu_SALT2[6][i])),
        mu_SALT2[4][i], mu_SALT2[5][i], mu_SALT2[6][i], sampleflag ) )


textfile_1.write(txtline)
text_050 = '# %s SNe Ia in this list.\n'%countSN
textfile_1.write(text_050)

textfile_1.close()
#-----------------------------------------------------------------------------80

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();
print "# All done smoothly"

