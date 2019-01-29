#!/usr/bin/env python
# coding: utf-8

# # Distance modulus $\mu$ from SNANA-SALT2 fit
#
# The distance modulus is determined from the SNANA-SALT2 fit and the global parameters reported in Table 6 of Scolnic et al. 2017 (https://arxiv.org/abs/1710.00845).
#
# This script outputs:
#     - a text file (with my format) of the distance moduli to create the Hubble diagram with 11.ipynb
#     - a simple Hubble diagram plot.

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad as intquad

# To read arguments in command line
# Used in the ".py" version of this notebook.
import sys

# # USER

#-----------------------------------------------------------------------------80

## Terminal or notebook version of this script?
ScriptVersion = 'notebook' # ( terminal , notebook )

#-----------------------------------------------------------------------------80

if ScriptVersion == 'notebook':
    # FITRES file:
    FitresFile = 'GP_subsample__.fitres'

    # Dir where the FITRES file is located:
    DirSALT2Fit = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/2_AndySampleOnly/BVr_mix_FUDGEMAGERROR/FUDGEALL_ITER1_MAXFRAC_0_20/fudge_0_010_CSP_only_BEST_PrevVersion/'

    # Peculiar velocity uncertainty:
    sigma_vPec = 150  # 150 km/s

    ## This corresponds to the SNANA version used
    ## during the computations. This is needed because different versions
    ## of SNANA produce different format of the output fitres files.
    ## Options: ( v10_58g , v10_35g , v10_45j , pantheon )
    FitresFormat = 'v10_58g'

    # Global parameters in the Tripp formula for: low-z.
    ## low-z values from Table 6 of Pantheon
#     alpha = 0.147
#     beta =  3.00
#     sigmaInt = 0.11

    # high-z values from Table 6 of Pantheon
    alpha = 0.165
    beta =  2.93
    sigmaInt = 0.09

#--------------------------------------------------
elif ScriptVersion == 'terminal':
    # FITRES file:
    FitresFile = sys.argv[1]

    # Dir where the FITRES file is located:
    DirSALT2Fit = sys.argv[2]

    # Peculiar velocity uncertainty:
    sigma_vPec = int(sys.argv[3])  # 150 km/s

    ## This corresponds to the SNANA version used
    ## during the computations. This is needed because different versions
    ## of SNANA produce different format of the output fitres files.
    ## Options: ( v10_58g , v10_35g , v10_45j , pantheon )
    FitresFormat = sys.argv[4]

    # Global parameters in the Tripp formula for: low-z.
    alpha = float(sys.argv[5])
    beta =  float(sys.argv[6])
    sigmaInt = float(sys.argv[7])

#--------------------------------------------------
# Dir to the save the output
DirSaveOutput = DirSALT2Fit

#--------------------------------------------------
# ORIGINAL
HoFix = 73.24  #
OmMFix = 0.28 # Omega_Matter
OmLFix = 0.72 # Omega_Lambda
wFix = -1.0 # Dark energy EoS

# TEMPORAL: SNANA default values
# HoFix = 70.0  # 70 = SNANA default value
# OmMFix = 0.3 # Omega_Matter # 0.3 = SNANA default value
# OmLFix = 0.7 # Omega_Lambda # 0.7 = SNANA default value
# wFix = -1.0 # Dark energy EoS # -1  = SNANA default value

MoFix = 19.36 # Absolute magnitude at T_Bmax. This is just a nuisance parameter
cc = 299792.458  # Speed of light (km/s)

NotebookName = '10_SALT2_DistanceMu.ipynb'

#-----------------------------------------------------------------------------80

# # Automatic
#
# ## Some definitions

# OK
"""
# Distance modulus from modified Tripp formula as in Scolnic+17

def muSALT2_func(mB, M, alpha, x1, beta, c1, gamma, tau, Mass, MassStep, Delta_mB, Delta_c, Delta_x1):

    DeltaM = gamma/(1 + np.exp(-(Mass-MassStep)/tau))
    DeltaB = Delta_mB - beta*Delta_c + alpha*Delta_x1
    distanceMu = mB - M + alpha*x1 - beta*c1 + DeltaM + DeltaB

    return distanceMu, DeltaM, DeltaB

# Test from table 16.
print '#', muSALT2_func(21.86, -18, 0.156, -2.3, 3.02, 0.05, 0.053, 0.001, 10.90, 10.13, 0.04, 0.025, 0.7)
# (39.476899999999993, 0.052999999999999999, 0.07369999999999999)
"""
0

# OK
"""
# Distance modulus from modified Tripp formula as in Rest+14

def muSALT2_Rest14_func(mB, M, alpha, x1, beta, c1):
    distanceMu = mB - M + alpha*x1 - beta*c1
    return distanceMu

#--------------------------------------------------
# Test

Mfix = -19.35
alphaFix = 0.147
betaFix = 3.13

print '# sn1990O (JRK07)', muSALT2_Rest14_func(15.993, Mfix, alphaFix, 0.476, betaFix, -0.054)
# sn1990O (JRK07) 35.581992

print '# sn2005el (CfA3)', muSALT2_Rest14_func(14.671, Mfix, alphaFix, -1.242, betaFix, -0.070)
# sn2005el (CfA3) 34.057526
"""
0

"""
# David Jones: original code

def salt2muDavid(z=None, x0=None, c=None, cerr=None, x1=None, x1err=None,
            mb=None, mberr=None, cov_x1_x0=None, cov_c_x0=None, cov_x1_c=None,
            alpha=None, beta=None, sigint=None, M=None,  peczerr=0.00083):

    sf = -2.5/(x0*np.log(10.0))
    cov_mb_c = cov_c_x0*sf
    cov_mb_x1 = cov_x1_x0*sf
    mu_out = mb + x1*alpha - beta*c + 19.36
    invvars = 1.0 / (mberr**2.+ alpha**2. * x1err**2. + beta**2. * cerr**2. + \
                         2.0 * alpha * (cov_x1_x0*sf) - 2.0 * beta * (cov_c_x0*sf) - \
                         2.0 * alpha*beta * (cov_x1_c) )

    zerr = peczerr*5.0/np.log(10)*(1.0+z)/(z*(1.0+z/2.0))
    muerr_out = np.sqrt(1/invvars + zerr**2. + 0.055**2.*z**2.)
    if sigint: muerr_out = np.sqrt(muerr_out**2. + sigint**2.)
    return(mu_out,muerr_out, np.sqrt(1/invvars))


# 2009ag
salt2muDavid(0.01019, 4.1224e-02, 6.1295e-02, 3.5659e-02, -1.2569e-01, 2.8985e-02,
       14.0971,  1.2470e-01, 4.1521e-05, -1.4239e-04, 4.1521e-05,
       alpha_test, beta_test, sigmaInt_test, MoFix, 0.0005003461427972281)
"""
0

# ### Function to convert from SALT2 parameters to distance modulus

#-----------------------------------------------------------------------------80
# Based on David Jones codes:

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


    # return (mu_out, muerr_out, np.sqrt(sigma2_mu_Tripp)) # original
    # TEMPORAL
    return (mu_out, muerr_out, np.sqrt(sigma2_mu_Tripp),
            mberr**2. + (alpha**2.)*(x1err**2.) + (beta**2.)*(cerr**2.),
            2.0*alpha*cov_mb_x1, -2.0*beta*cov_mb_c, -2.0*alpha*beta*cov_x1_c)

# print '# Test'
alpha_test = 0.147
beta_test =  3.00
sigmaInt_test = 0.11

# print '# 2006le (PS1s_CFA3_KEPLERCAM_RS14):', \
# salt2mu(0.0172, 0.26064E-01, -0.0669,  0.0640,  0.8947,  0.0506,
#         14.5949,  0.2197,   0.895E-04,  -0.349E-03,  -0.104E-02,
#         alpha_test, beta_test, sigmaInt_test, MoFix, sigma_vPec)

# print '# 2006lf (PS1s_CFA3_KEPLERCAM_RS14):', \
# salt2mu(0.0124, 0.78991E-01, -0.2086,  0.1013, -1.3905, 0.0674,
#         13.3911,  0.3646, 0.228E-03,  -0.376E-02,  -0.721E-03,
#         alpha_test, beta_test, sigmaInt_test, MoFix, sigma_vPec)

#-----------------------------------------------------------------------------80
# Inverse of the dimensionless Hubble parameter
def InvEHubblePar(z, OmM, wde):
    "Dimensionless Hubble parameter"
    InvEHubbleParInt = 1.0/(np.sqrt(OmM*(1.0+z)**3.0 + (1.0-OmM)*(1.+z)**(3.*(1.+wde))))
    return InvEHubbleParInt

# ---- The luminosity distance ----
def LumDistanceVec(z, OmM, wde, Ho):
    "Luminosity distance"
    LumDistanceVecInt = 0.
    LumDistanceVecInt = cc*(1.+z)*intquad(InvEHubblePar, 0., z, args=(OmM, wde))[0]/Ho
    return LumDistanceVecInt

# ---- Distance modulus scalar ----
def DistanceMu(z, OmM, wde, Ho):
    "Distance modulus"
    DistanceMuInt = 5.0*np.log10(LumDistanceVec(z, OmM, wde, Ho)) + 25.0
    return DistanceMuInt

# ---- Distance modulus Vector ----
def DistanceMuVector(z, OmM, wde, Ho):
    "Distance modulus"
    DistanceMuInt= []
    for i in range(len(z)):
        DistanceMuInt += [5.0*np.log10(LumDistanceVec(z[i], OmM, wde, Ho)) + 25.0]
    return DistanceMuInt

#--------------------------------------------------
ztest1 = 0.01

print '# Checking that the functions work well:', DistanceMu(ztest1, OmMFix, wFix, HoFix)
# Checking that the functions work well: 33.0773926577

# sigma^2_mu from the peculiar velocity uncertainty
# This function is used to determine in the sections "Intrinsic dispersion"
# and "Optical RMS", to determine the intrinsic dispersion.

def sigma2_pec(zcmb, err_zcmb, vpec):
    sigma2_pecInt = ((5/(zcmb*np.log(10)))*np.sqrt((vpec/cc)**2 + err_zcmb**2))**2
    return sigma2_pecInt

# Test
sigma2_pec(0.0109942726, 0.0010420420, 150)
# 0.052125037354329877

# ### Cepheid distances

# OK but not used
"""
DirSNeWithCepheid = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp/MyNotesAbout/'

ListSNeCepheid = np.genfromtxt(DirSNeWithCepheid+'SNeWithCepheidDistances.txt', dtype=['S10',
                                                float,float,float,float,float,float])
# ListSNeCepheid = np.genfromtxt(DirSNeWithCepheid+'SNeWithCepheidDistances_hack.txt', dtype=['S10',
#                                                 float,float,float,float,float,float])

ListSNeCepheid['f0']
"""
0

# #### Get the name of this ipython notebook
# To print it in the output text files as reference

# %%javascript
# var kernel = IPython.notebook.kernel;
# var thename = window.document.getElementById("notebook_name").innerHTML;
# var command = "NotebookName = " + "'"+thename+".ipynb"+"'";
# kernel.execute(command);

print '#', (NotebookName)
# Update_zcmb_in_SNANA_datafiles_v1_0.ipynb

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

# # Compute $(\mu, \sigma_{\mu})$

#-----------------------------------------------------------------------------80
# Read the FITRES file

#  SNANA_v10_35g fitres file:
if FitresFormat == 'v10_35g':
    print "SNANA file version = v10_35g"
    print "NOTE: suffix a column at the end of each row with the subsample flag"
    try:
        snanaSalt2Fit_np =np.genfromtxt(DirSALT2Fit+FitresFile[:-7]+'Notes_.fitres',
                                        skip_header=2, comments='#',
                                        dtype=['S3', 'S15',
                  float, float, float, float, float, float, float, float,
                  float, float, float, float, float, float, float, float,
                  float, float, float, float, float, float, float, float] )

    except:
        snanaSalt2Fit_np =np.genfromtxt(DirSALT2Fit+FitresFile,
                                        skip_header=2, comments='#',
                                        dtype=['S3', 'S15',
                  float, float, float, float, float, float, float, float,
                  float, float, float, float, float, float, float, float,
                  float, float, float, float, float, float, float, float] )

#  SNANA_v10_58g fitres file:
elif FitresFormat == 'v10_58g':
    print "SNANA file version = v10_58g"
    snanaSalt2Fit_np = np.genfromtxt(DirSALT2Fit+FitresFile, skip_header=8,
                        dtype=['S3', 'S15', float, float, 'S4',
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float ])

#  Pantheon fitres file:
elif FitresFormat == 'pantheon':
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


#  David Jones fitres file:
elif FitresFormat == 'v10_45j':
    snanaSalt2Fit_np =np.genfromtxt(DirSALT2Fit+FitresFile,
                    skip_header=2, comments='#',
                    dtype=['S3', 'S15',float, float, 'S4',
      float, float, float, float, float, float, float, float, float, float,
      float, float, float, float, float, float, float, float, float, float,
      float, float, float, float, float, float, float, float, float, float] )

print '# %s SNe in this file.'%(len(snanaSalt2Fit_np))
# 56 SNe in this file.`

#-----------------------------------------------------------------------------80
# Computing the distance modulus

# Define some variable arrays:

# Reading my own fitres files:
if FitresFormat == 'v10_35g':
    z_np =  snanaSalt2Fit_np['f2'] # z_CMB
    zerr_np = snanaSalt2Fit_np['f3'] # Not needed to compute the distance modulus.
    x0_np = snanaSalt2Fit_np['f4']
    c_np = snanaSalt2Fit_np['f6']
    cerr_np = snanaSalt2Fit_np['f7']
    x1_np = snanaSalt2Fit_np['f8']
    x1err_np = snanaSalt2Fit_np['f9']
    mb_np =  snanaSalt2Fit_np['f12']
    mberr_np =  snanaSalt2Fit_np['f13']
    cov_x1_c_np =  snanaSalt2Fit_np['f16']
    cov_x1_x0_np =  snanaSalt2Fit_np['f14']
    cov_c_x0_np =  snanaSalt2Fit_np['f15']
    chi2_np = snanaSalt2Fit_np['f17'] # Not needed to compute the distance modulus.
    Ndof_np = snanaSalt2Fit_np['f18']  # Not needed to compute the distance modulus.

# Reading my own fitres files:
elif FitresFormat == 'v10_58g':
    z_np =  snanaSalt2Fit_np['f10'] # z_CMB
    zerr_np = snanaSalt2Fit_np['f11'] # Not needed to compute the distance modulus.
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
    Ndof_np = snanaSalt2Fit_np['f32']  # Not needed to compute the distance mu.
    chi2_np = snanaSalt2Fit_np['f33'] # Not needed to compute the distance mu.

elif FitresFormat == 'pantheon':
    z_np =  snanaSalt2Fit_np['f9'] #
    zerr_np = snanaSalt2Fit_np['f10'] # Not needed to compute the distance modulus.
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
    Ndof_np = snanaSalt2Fit_np['f31']  # Not needed to compute the distance mu.
    chi2_np = snanaSalt2Fit_np['f32'] # Not needed to compute the distance mu.

#   David Jones fitres file:
elif FitresFormat == 'v10_45j':
    z_np =  snanaSalt2Fit_np['f6'] # z_CMB
    zerr_np = snanaSalt2Fit_np['f7'] # Not needed to compute the distance modulus.
    x0_np = snanaSalt2Fit_np['f25']
    c_np = snanaSalt2Fit_np['f21']
    cerr_np = snanaSalt2Fit_np['f22']
    x1_np = snanaSalt2Fit_np['f19']
    x1err_np = snanaSalt2Fit_np['f20']
    mb_np =  snanaSalt2Fit_np['f23']
    mberr_np =  snanaSalt2Fit_np['f24']
    cov_x1_c_np =  snanaSalt2Fit_np['f27']
    cov_x1_x0_np =  snanaSalt2Fit_np['f28']
    cov_c_x0_np =  snanaSalt2Fit_np['f29']
    chi2_np = snanaSalt2Fit_np['f31'] # Not needed to compute the distance modulus.
    Ndof_np = snanaSalt2Fit_np['f30']  # Not needed to compute the distance modulus.

#-------------------
# Computing the distance modulus:
mu_SALT2 = salt2mu(z_np, x0_np, c_np, cerr_np, x1_np, x1err_np,
                   mb_np, mberr_np, cov_x1_x0_np, cov_c_x0_np, cov_x1_c_np,
                  alpha, beta, sigmaInt, MoFix, sigma_vPec)

if ScriptVersion == 'notebook':
    # Inspect for 'nan' values in the total distance-modulus uncertanty.
    # It means there is something wrong for that SN. I could 'comment' that SN.

    print "# Total distance-modulus uncertanty:"
    mu_SALT2[1]

if ScriptVersion == 'notebook':
    # Inspect for negative values for the distance-modulus uncertainty that
    # comes from Tripp formula. It means there is something wrong for that SN.
    # I could 'comment' that SN.

    print "# Distance-modulus uncertanty from propagating the Tripp formula:"
    mu_SALT2[2]

# # Determine the mean absolute magnitude

# Define the average in the absolute value of the residual
# These functions are going to be minimized.

# AbsResidual values
res_mu_np_int = mu_SALT2[0] - DistanceMuVector(z_np, OmMFix, wFix, HoFix)

# Inverse-variance weighted average function to be minimized
def WeightedAbsResidual_fun(deltaMo, UpLimit, LowLimit, AbsResidual_Roof):
    residuals_int = res_mu_np_int + deltaMo
    WeightedAbsResidual_int = np.average(residuals_int, weights=mu_SALT2[2] )

    if deltaMo < UpLimit and deltaMo > LowLimit:
        # For some unknown reason, simplex is search the maximum instead of
        # the minimimum, so I had to define the average absolute mag as 1/WeightedAbsResidual
        # WeightedAbsResidual_final = 1/WeightedAbsResidual_int
        WeightedAbsResidual_final = WeightedAbsResidual_int

    else: WeightedAbsResidual_final = AbsResidual_Roof

    return abs(WeightedAbsResidual_final)

#--------------------------------------------------
# Simple average function to be minimized
def AverageAbsResidual_fun(deltaMo, UpLimit, LowLimit, AbsResidual_Roof):
    residuals_int = res_mu_np_int + deltaMo
    # AverageAbsResidual_int = np.mean(residuals_int)
    AverageAbsResidual_int = np.sum(residuals_int)/len(residuals_int)

    if deltaMo < UpLimit and deltaMo > LowLimit:
        # For some unknown reason, simplex is search the maximum instead of
        # the minimimum, so I had to define the average absolute mag as 1/AverageAbsResidual_int
        # AverageAbsResidual_final = 1/AverageAbsResidual_int
        AverageAbsResidual_final = AverageAbsResidual_int

    else: AverageAbsResidual_final = AbsResidual_Roof

    return abs(AverageAbsResidual_int)

if ScriptVersion == 'notebook':
    print '# Tests:'
    print '#', WeightedAbsResidual_fun(-2.3, 3, -3, 1)
    print '#', AverageAbsResidual_fun(-0.0153, 0.2, -0.2, 1)

    # Tests:
    # 2.22003154831089
    # 0.06383338274987903

# Determine the value of the constant deltaMo in order to obtain a
# weighted average value of zero in the Hubble residual

import scipy.optimize

# Search limits:
# UpLimit = 0.2; LowLimit = -0.2
UpLimit = 3; LowLimit = -3

# Assume this value for deltaMo in case of the search is outside the range
# indicated above.
Residual_Roof = 10

WeightedAbsResidual_Out = scipy.optimize.minimize_scalar(WeightedAbsResidual_fun,
                                        args=(UpLimit, LowLimit, Residual_Roof))

AverageAbsResidual_Out = scipy.optimize.minimize_scalar(AverageAbsResidual_fun,
                                        args=(UpLimit, LowLimit, Residual_Roof))

# Redefining the values:
SimplexResult_1 = [WeightedAbsResidual_Out['x'] ]
SimplexResult_2 = [AverageAbsResidual_Out['x'] ]

print WeightedAbsResidual_Out
print '#', SimplexResult_1, ' = value of deltaMo that minimize the Hubble residual.'
print

print AverageAbsResidual_Out
print '#', SimplexResult_2, ' = value of deltaMo that minimize the Hubble residual.'

"""
     fun: 6.2205205415565593e-12
    nfev: 28
     nit: 27
 success: True
       x: 0.0045658128791809327
# [0.0045658128791809327]  = value of deltaMo that minimize the Hubble residual.

     fun: 6.6409884697359133e-12
    nfev: 30
     nit: 29
 success: True
       x: 0.003668901213883702
# [0.003668901213883702]  = value of deltaMo that minimize the Hubble residual.
"""
0

# Verifying the minimim using a different python function to find the minimum.

# Determine the value of the constant deltaMo in order to obtain a
# weighted average value of zero in the Hubble residual

from scipy.optimize import fmin

# Search limits:
# UpLimit_test = 0.1; LowLimit_test = -0.1
UpLimit_test = 3; LowLimit_test = -3

# Assume this value for deltaMo in case of the search is outside the range
# indicated above.
Residual_Roof_test = 10

InitialGuess_test = -0.0423

# Find the value of deltaMo that minimizes the averager functions
SimplexResult_1_test = fmin(WeightedAbsResidual_fun, InitialGuess_test, # xtol=1e-10, ftol=1e-10,
                          # retall=True,
                          args=(UpLimit_test, LowLimit_test, Residual_Roof_test) )
print SimplexResult_1_test, ' = value of deltaMo that minimize the Hubble residual.'
print

SimplexResult_2_test = fmin(AverageAbsResidual_fun, InitialGuess_test,  # xtol=1e-10, ftol=1e-10,
                          # retall=True,
                          args=(UpLimit_test, LowLimit_test, Residual_Roof_test) )
print SimplexResult_2_test, ' = value of deltaMo that minimize the Hubble residual.'

if ScriptVersion == 'notebook':
    # Print results

    print '# deltaMo =', SimplexResult_1, ', for a weighted average abs residual of:',     WeightedAbsResidual_fun(SimplexResult_1[0], UpLimit, LowLimit, Residual_Roof)

    print '# deltaMo =', SimplexResult_2, ', for an average abs residual of:',     AverageAbsResidual_fun(SimplexResult_2[0], UpLimit, LowLimit, Residual_Roof)

# deltaMo min = [-0.03932578] , for a weighted average abs residual of: 3.21942477339e-05
# deltaMo min = [-0.03502969] , for an average abs residual of: 8.94368852025e-06

# ### Compute final residual values and check the Hubble residual

# Inverse-variance weighted average

res_mu_weight_np = (mu_SALT2[0] - DistanceMuVector(z_np, OmMFix, wFix, HoFix)
                    +SimplexResult_1[0])
WeightedResidual = np.average(res_mu_weight_np, weights=mu_SALT2[2] )

print "# NOTE: THE FOLLOWING VALUE SHOULD BE VERY SMALL (< 1e-5), OTHERWISE, IT'S WRONG."
print '#', WeightedResidual
# 3.21942477339e-05

# Simple average

res_mu_mean_np = (mu_SALT2[0] - DistanceMuVector(z_np, OmMFix, wFix, HoFix)
                  +SimplexResult_2[0])

print "# NOTE: THE FOLLOWING VALUE SHOULD BE VERY SMALL (< 1e-5), OTHERWISE, IT'S WRONG."
print '#', np.mean(res_mu_mean_np)
# NOTE: THE FOLLOWING VALUE SHOULD BE VERY SMALL (~ 1e-10), OTHERWISE, IT S WRONG.
# 8.94368852025e-06

# # Write the results to a text files

# ### My format
#
# To be read by '11_DistanceMu_HubbleDiagram.ipynb'

#-----------------------------------------------------------------------------80

# Write down to a text file the results using the same datatable
# format output than my 11_DistanceMu_HubbleDiagram_v2_17.ipynb code.

textfile_1 = open(DirSaveOutput+'DistanceMu_Good_AfterCutoffs_Main_.txt','w')
# To share with collaborators:
textfile_2 = open(DirSaveOutput+'DistanceMu_fromSALT2.txt','w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Author = '# Data table created by: Arturo Avelino \n'
text_Date   = '# On date: %s \n'%text_timenow
text_script = '# Script used: %s \n'%NotebookName
text_line = '#'+'-'*70 + '\n'

text_0010 = '# FITRES files format: %s\n'%FitresFormat

text_01 = '# Distance modulus determined from the combination of SALT2.JLA-B14 \
output fit \n# parameters and the Tripp formula, using the global parameters \
reported by Scolnic etal. 2017. \n'
text_02 = '# Values of the global parameters used: \n'
text_03 = '# alpha = %s | beta = %s |  sigmaInt = %s \n'%(alpha,beta,sigmaInt)
text_04 = '# M = %s. Assumed value for the absolute magnitude in Tripp \
formula. \n'%MoFix

text_05 = "# Delta_M = %s. Additional value to the absolute magnitude in Tripp \n\
# formula to have a zero-value weighted average in the Hubble residual plot. It is, \n\
# the value of 'Delta_M' was -added- to the net distance-modulus value computed directly \n\
# from the SALT2 parameters using David Jones formula. So to recover the actual SALT2 \n\
# distance modulus, I have to subtract Delta_M to 'mu' reported in this table. \
\n"%SimplexResult_1[0]

text_06 = '# Om_matter = %s | w = %s | Ho = %s km/s/Mpc | Flat. \
Cosmological parameters assumed. \n'%(OmMFix, wFix, HoFix)
text_07 = '# %s km/s: Peculiar velocity uncertainty assumed.  \n'%sigma_vPec

text_08_1 = '# SN name           z_CMB   error_zcmb   mu     \
error_mu   mu_residual chi2_dof  \
Sample  appMag_TBmax  error_appMagTBmax  mu_LCDM     sigma_muLCDM_vPec AbsMagTBmax \
error_AbsMagTBmax  TBmax       Error_TBmax    PhaseB        zhelio        err_zhelio   \
dm15         error_dm15      EBVhost      error_EBVhost   EBV_MW        error_EBV_MW     \
Alamb      err_Alamb      R_F           mu_Snoopy     error_muSnoopy   TBmax           \
Error_TBmax  appMag_TBmax  error_appMagTBmax  Notes  \n'

text_08_2 = '# SN name           z_CMB   error_zcmb   mu     error_mu   \
appMag_TBmax  error_appMagTBmax  mu_LCDM     sigma_muLCDM_vPec \n'

textfile_1.write(text_01); textfile_1.write(text_02); textfile_1.write(text_03);
textfile_1.write(text_04); textfile_1.write(text_05); textfile_1.write(text_06);
textfile_1.write(text_07);
textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date);
textfile_1.write(text_script);
textfile_1.write(text_line)
textfile_1.write(text_0010)
textfile_1.write(text_08_1)

textfile_2.write(text_01); textfile_2.write(text_02); textfile_2.write(text_03);
textfile_2.write(text_04); textfile_2.write(text_05); textfile_2.write(text_06);
textfile_2.write(text_07);
textfile_2.write(text_line)
textfile_2.write(text_Author); textfile_2.write(text_Date); textfile_2.write(text_script);
textfile_2.write(text_line)
textfile_2.write(text_0010)
textfile_2.write(text_08_2)

#-----------------------------------------------------------------------------80
# Reset
countSN = 0;
SampleFlag = 0

for i in range(len(snanaSalt2Fit_np)):

    snname = snanaSalt2Fit_np['f1'][i]
    zz_int = z_np[i]
    e_zz_int = zerr_np[i]
    muLCDM_int = DistanceMu(zz_int, OmMFix, wFix, HoFix)
    sigma_muLCDM_vPec_int = np.sqrt(sigma2_pec(zz_int, e_zz_int, sigma_vPec))
    chi2dof_int = chi2_np[i]/Ndof_np[i]
    mu_final = mu_SALT2[0][i] +SimplexResult_1[0]
    muerr_final = mu_SALT2[2][i] # Using the Tripp uncertainty.
    muResidual_int = res_mu_weight_np[i]
    mb_int = mb_np[i]
    mberr_int = mberr_np[i]

    if FitresFormat == 'v10_35g' :
        SampleFlag = snanaSalt2Fit_np['f25'][i]
    elif FitresFormat == 'v10_58g' :
        SampleFlag = snanaSalt2Fit_np['f35'][i]
    elif FitresFormat == 'pantheon' :
        SampleFlag = snanaSalt2Fit_np['f4'][i]
    else: SampleFlag = 1


    textfile_1.write('sn%-16s  %.5f  %.5f  %.5f  %.5f  %10.6f  %10.6f    %.0f     %.4f       %.4f            %.8f  %.8f   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%(
        snname, zz_int, e_zz_int,
        mu_final, muerr_final, muResidual_int,
        chi2dof_int, SampleFlag,
        mb_int, mberr_int,
        muLCDM_int, sigma_muLCDM_vPec_int) )

    #-------

    textfile_2.write('sn%-16s  %.5f  %.5f  %.5f  %.5f     %.4f       %.4f            %.8f  %.8f  \n'%(
        snname, zz_int, e_zz_int,
        mu_final, muerr_final,
        mb_int, mberr_int,
        muLCDM_int, sigma_muLCDM_vPec_int) )

    countSN = countSN + 1

text_12 = '# %s SNeIa in this list.'%countSN
textfile_1.write(text_line);textfile_2.write(text_line)
textfile_1.write(text_12);  textfile_2.write(text_12);

textfile_1.close(); textfile_2.close();
#-----------------------------------------------------------------------------80

textfile_1.close(); textfile_1.close(); textfile_1.close();
textfile_2.close(); textfile_2.close(); textfile_2.close();

# #### SNANA-SALT2mu.exe format
#
# This doesn't work yet for the case of reading David Jones fitres file or any other format apart of "v10_35g"

# Write down to a text file the results using the same datatable
# format output than SALT2mu.exe.

if FitresFormat == 'v10_35g':
    textfile_1 = open(DirSaveOutput+'DistanceMu_SALT2muFormat_.txt','w')

    now = datetime.datetime.now() # Read the time and date right now
    text_timenow = now.strftime("%Y-%m-%d (yyyy/mm/dd); %H:%M hrs.")
    text_line = '#'+'-'*80 + '\n'

    #--- Copy/paste of the cell above -----
    textfile_1.write('# Distance modulus determined from the combination of SALT2.JLA-B14 \
    output fit parameters and the Tripp formula, using the global parameters reported \
    by Scolnic etal. 2017. \n')
    textfile_1.write('# Values of the global parameters used: \n')
    textfile_1.write('# alpha = %s | beta = %s |  sigmaInt = %s \n'%(alpha,beta,sigmaInt))
    textfile_1.write('# M = %s. Assumed value for the absolute magnitude in Tripp \
    formula. \n'%MoFix)

    textfile_1.write(text_05)

    textfile_1.write('# Om_matter = %s | w = %s | Ho = %s km/s/Mpc | Flat. \
    Cosmological parameters assumed. \n'%(OmMFix, wFix, HoFix))
    textfile_1.write('# %s km/s: Peculiar velocity uncertainty assumed.  \n'%sigma_vPec)
    textfile_1.write(text_line)

    textfile_1.write('# Data table created by: Arturo Avelino \n')
    textfile_1.write('# On date: %s \n'%text_timenow)
    textfile_1.write('# Script used: %s \n'%NotebookName)
    textfile_1.write(text_line)

    #------- end copy/paste -------------

    textfile_1.write('# SN name       z        zERR     x0           x0ERR       c       cERR     \
    x1      x1ERR   PKMJD   PKMJDERR mB      mBERR    COVx0x1     COVx0c      COVx1c     CHI2   NDOF  \
    FITPROB      SNRMAX1    SNRMAX2    SNRMAX3 SURVEY TYPE   MU      MUERR    MU_residual   \n')

    for i in range(len(snanaSalt2Fit_np)):

        mu_final = mu_SALT2[0][i] +SimplexResult_1[0]

        textfile_1.write('SN: %-10s  %.5f  %.5f  %.5e  %.3e  %7.4f  %.4f  %7.4f  %.4f  %.3f  %.2f     %.4f  %.4f  %10.3e  %10.3e  %10.3e  %5.1f  %4.0f  %.3e  %9.3f  %9.3f  %9.3f  %2.0f     %2.0f      %.5f  %.5f  %10.6f   \n'%(
            snanaSalt2Fit_np['f1'][i], snanaSalt2Fit_np['f2'][i], snanaSalt2Fit_np['f3'][i],
            snanaSalt2Fit_np['f4'][i], snanaSalt2Fit_np['f5'][i], snanaSalt2Fit_np['f6'][i],
            snanaSalt2Fit_np['f7'][i], snanaSalt2Fit_np['f8'][i], snanaSalt2Fit_np['f9'][i],
            snanaSalt2Fit_np['f10'][i], snanaSalt2Fit_np['f11'][i],
            snanaSalt2Fit_np['f12'][i], snanaSalt2Fit_np['f13'][i],
            snanaSalt2Fit_np['f14'][i], snanaSalt2Fit_np['f15'][i], snanaSalt2Fit_np['f16'][i],
            snanaSalt2Fit_np['f17'][i], snanaSalt2Fit_np['f18'][i], snanaSalt2Fit_np['f19'][i],
            snanaSalt2Fit_np['f20'][i], snanaSalt2Fit_np['f21'][i], snanaSalt2Fit_np['f22'][i],
            snanaSalt2Fit_np['f23'][i], snanaSalt2Fit_np['f24'][i],
            mu_final, mu_SALT2[2][i], res_mu_weight_np[i]  ))

    textfile_1.close()
#-----------------------------------------------------------------------------80

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();

# # Hubble diagram plot

#-----------------------------------------------------------------------------80
# Plot settings:

fontsizePlot = 12

#--------------------------------------------------
#    Theoretical values
# Array plot the -theoretical- (spectroscopic) distance modulus
nbins1= 51
z1 = np.linspace(min(z_np)-0.001, max(z_np)+0.001, nbins1)
mu0 = DistanceMuVector(z1, OmMFix, wFix, HoFix)

#--------------------------------------------------
# PLOTTING

plt.errorbar(z_np, mu_SALT2[0]+SimplexResult_1[0], yerr=mu_SALT2[2], fmt='.')
plt.plot(z1, mu0)

plt.grid(True)
plt.xlabel('Redshift', fontsize=fontsizePlot)
plt.ylabel(r'Distance modulus $\mu$', fontsize=fontsizePlot)
plt.title('SALT2 Hubble diagram')
plt.tight_layout()

plt.savefig(DirSaveOutput+'Plot_Hubble_SALT2_.png')
plt.close()

plt.close();plt.close();plt.close();
plt.close();plt.close();plt.close();

# # Hubble residual plot

#    Theoretical values
# Array plot the -theoretical- (spectroscopic) distance modulus
nbins1= 51
z1 = np.linspace(min(z_np)-0.001, max(z_np)+0.001, nbins1)
mu1 = np.zeros(len(z1))

#--------------------------------------------------
# res_mu_np_int
# res_mu_weight_np

plt.errorbar(z_np, res_mu_weight_np, yerr=mu_SALT2[2], fmt='.')

plt.plot(z1, mu1)
# plt.text(min(z_np)+0.01, min(res_mu_weight_np), r'$\sigma_{\rm pec}$ = %s km/s'%sigma_vPec)
plt.grid(True)
plt.xlabel('Redshift', fontsize=fontsizePlot)
plt.ylabel(r'$\mu - \mu_{\rm \Lambda CDM}$', fontsize=fontsizePlot+3)
plt.title('SALT2 Hubble residual')
plt.tight_layout()

plt.savefig(DirSaveOutput+'Plot_HubbleRes_SALT2_.png')
plt.close()

plt.close();plt.close();plt.close();
plt.close();plt.close();plt.close();

print "# All done smoothly."

