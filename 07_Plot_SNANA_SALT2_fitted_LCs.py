#!/usr/bin/env python
# coding: utf-8

# # Plot the SNANA-SALT2 fits

import numpy as np
from matplotlib import pyplot as plt
import os # To use command line like instructions
import glob # To read the name of the files in a given directory
# To read arguments in command line
# Used in the ".py" version of this notebook.
import sys

#-----------------------------------------------------------------------------80
## USER

## Terminal or notebook version of this script?
ScriptVersion = 'notebook' # ( terminal , notebook )

#-------------------
## Notebook version

if ScriptVersion == 'notebook':

    ## Data file to be plotted:
    filefit = 'PS1s_CSPDR2_V_RS14.LCPLOT.TEXT'

    # Directory where the data file is located:
    DirMain = '/Users/arturo/Downloads/PS1s_CSPDR2_V_RS14_Fudge_0_02/'

#-------------------
# Terminal version

elif ScriptVersion == 'terminal':

    ## Directory where the data file is located:
    DirMain = sys.argv[1]

    ## Data file to be plotted:
    filefit = sys.argv[2]

# # Automatic

#- Force the creation of the directory to save the outputs.
#- "If the subdirectory does not exist then create it"
DirSaveOutput = DirMain+'plots/'
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

# Upload the data

# Reset flag to indicate that the file was not found and skip
# the rest of the code.
flag_no_error = 1

try:
    fit_data = np.genfromtxt(DirMain+filefit, skip_header=2,
                         dtype=['S3', 'S15', float, float,
                                float, float, float, 'S4', float, float])
except:
    print "# SNANA Light curve file not found."
    flag_no_error = 0 # Set to false.

# ## Create a list with the number of rows for each band

#-----------------------------------------------------------------------------80
# Create the list of rows

if flag_no_error:

    # Initialize-reset some values
    count_rows_dataflag = 0;
    count_rows_fit = 0;
    count_rows_band = 0;
    count_rows_SN = 0;
    count_sne = 0;

    snname1 = ''; band1 = '';
    snname2 = ''; band2 = ''; flagdata2 = 1

    rows_list = []; rows_list_int1 = [];

    #-----------------------------------

    for i3 in range(len(fit_data)):

        snname2 = fit_data['f1'][i3]
        band2 = fit_data['f7'][i3]
        flagdata2 = int(fit_data['f6'][i3])

        # Very initial values for these variables.
        if count_rows_SN == 0:
            # Flag to indicate I will write the name of the SN.
            flag_writesnname = 0 # reset
            snname1 = snname2; band1 = band2;
            # print 'Initial values:', snname1, band1

        if (snname2 == snname1 and i3 != (len(fit_data)-1) ):
            count_rows_SN += 1

            if (band2 == band1 and i3 != (len(fit_data)-1) ):
                count_rows_band += 1

                if (flagdata2 == 1 or flagdata2== -1):
                    count_rows_dataflag += 1
                else: count_rows_fit += 1

            elif (band2 != band1 or i3 == (len(fit_data)-1) ):

                # Write the SN name and first filter array:
                if flag_writesnname == 0:

                    # Fix the fact that for the first SNe and filter it is
                    # written -1 instead of 0.
                    if count_sne == 0: row_ini_int = 0
                    else: row_ini_int = i3-count_rows_band-1

                    rows_list_int1 += [snname1, [band1, row_ini_int,
                                       i3-count_rows_fit-1, i3-1] ]
                    flag_writesnname = 1 # update flag

                else:
                    rows_list_int1 += [ [band1, i3-count_rows_band-1,
                                         i3-count_rows_fit-1, i3-1] ]
                # Update
                band1 = band2

                # Reset:
                count_rows_dataflag = 0;
                count_rows_fit = 0;
                count_rows_band = 0;

        elif (snname2 != snname1 or i3 == (len(fit_data)-1) ):

            rows_list_int1 += [ [band1, i3-count_rows_band-1,
                                 i3-count_rows_fit-1, i3-1] ]

            rows_list += [rows_list_int1]
            # print 'rows_list:', rows_list

            snname1 = snname2
            count_sne += 1

            # Reset:
            rows_list_int1 = [];
            count_rows_dataflag = 0;
            count_rows_fit = 0;
            count_rows_band = 0;
            count_rows_SN = 0;

    # print rows_list

# # Plotting

#-----------------------------------------------------------------------------80
if flag_no_error:

    #     Plot settings

    Myfontsize = 12

    # Data dots:
    MarkerSize = 6
    BarWidth = 1
    MyCapeSize = 2
    fmt_int = '.'

    # For the confidence-interval plots (1.=68.3% CL or 1.9600=95% CL )
    CL_size = 1.

    count_SNe_plotted = 0;

    #------------------------

    for i2 in range(len(rows_list)): # loop over SNe.

        # Number of fitted bands per supernova
        nbands = len(rows_list[i2]) - 1

        plt.clf()
        fig, axes = plt.subplots(1, nbands, sharex=False, sharey=True,
                                 figsize=(12, 5))

        # add a big axes, hide frame
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axes
        plt.tick_params(labelcolor='none', top=False, bottom=False,
                        left=False, right=False)

        # plt.tick_params(labelcolor='none', top='off', bottom='off',
        #                 left='off', right='off')

        # Close the space between the subplots
        plt.subplots_adjust(wspace = .2)

        for i1 in range(1,(nbands+1)): # Loop over bands

            snname = rows_list[i2][0]
            band = rows_list[i2][i1][0]
            row_ini = rows_list[i2][i1][1]
            row_end = rows_list[i2][i1][2]
            rowfit_ini = row_end + 1
            rowfit_end = rows_list[i2][i1][3]

            fit_tt = np.atleast_2d(fit_data['f3'][rowfit_ini:rowfit_end]).T

            # The LC data
            axes[i1-1].errorbar(fit_data['f3'][row_ini:row_end],
                                fit_data['f4'][row_ini:row_end],
                         yerr = fit_data['f5'][row_ini:row_end],
                             fmt=fmt_int, color='black',  ms=MarkerSize,
                             elinewidth=BarWidth, capsize=MyCapeSize )

             # Confidence level band of the SALT2 best fit model
            axes[i1-1].fill(np.concatenate([fit_tt, fit_tt[::-1]]),
                            np.concatenate([fit_data['f4'][rowfit_ini:rowfit_end] -
                                  CL_size * fit_data['f5'][rowfit_ini:rowfit_end],
                                           (fit_data['f4'][rowfit_ini:rowfit_end] +
                                  CL_size * fit_data['f5'][rowfit_ini:rowfit_end])[::-1]]),
                                  alpha=0.5, fc='grey')

            # The SALT2 best fit model
            axes[i1-1].plot(fit_data['f3'][rowfit_ini:rowfit_end],
                            fit_data['f4'][rowfit_ini:rowfit_end],
                            color='red')

            axes[i1-1].grid(ls='--', alpha=0.5)
            axes[i1-1].set_xlim(-16,60)
            axes[i1-1].set_xlabel('Phase (days)')
            axes[i1-1].set_ylabel('Flux', fontsize =Myfontsize)
            axes[i1-1].set_title('%s: %s band'%(snname,band))

        #----------------------------------
        # plt.title('Scatter in the Hubble residuals: NIR vs Optical')

        plt.savefig(DirSaveOutput+'%s_SALT2fit_.png'%(snname))
        plt.close()
        count_SNe_plotted += 1
#-----------------------------------------------------------------------------80

if flag_no_error:
    plt.close(); plt.close(); plt.close();plt.close(); plt.close();
    plt.close(); plt.close(); plt.close();plt.close(); plt.close();

if flag_no_error:
    print "# done: %s SNe Ia plotted their SALT2 fit."%count_SNe_plotted
else:
    print 'Please find the location and name of the light curve file to be plotted.'

