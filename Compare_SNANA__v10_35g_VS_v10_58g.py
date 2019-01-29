#!/usr/bin/env python
# coding: utf-8

# # Compare values computed with SNANA versions v10_35g vs v10_58g
#
# This is very useful. I use the version v10_35g to plot the SALT2 fits, because the Odyssey/SNANA_v10_58g was compiled without ROOT.
# The files to compare are the fitres files.

# ## User

import numpy as np

#--------------------------------------------------
#    v10_35g

DirTables_1 = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/tests_ok/LOWZ_JRK07/snana_v10_35g_salt2ext/'

DataFile_1 = 'LOWZ_JRK07.fitres'

fitresFile_1 = np.genfromtxt(DirTables_1+DataFile_1, skip_header=2,
                            dtype=['S3', 'S15',
                                  float, float, float, float, float, float, float,
                                  float, float, float, float, float, float, float, float,
                                  float, float, float, float, float, float, float, float])

print "# %s SNe in this table 1."%len(fitresFile_1)

#-----------------------------------------------------------------------------80
#   v10_58g

DirTables_2 = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/tests_ok/LOWZ_JRK07/snana_v10_58g_salt2ext/'

DataFile_2 = 'snfit_LOWZ_JRK07.FITRES.TEXT'

fitresFile_2 = np.genfromtxt(DirTables_2+DataFile_2, skip_header=7,
                            dtype=['S3', 'S15',
                                  float, float, 'S4',  float, float, float, float, float,
                                  float, float, float, float, float, float, float, float,
                                  float, float, float, float, float, float, float, float,
                                  float, float, float, float, float, float, float, float,
                                  float])

print "# %s SNe in this table 2."%len(fitresFile_2)

#-----------------------------------------------------------------------------80
#    Dir to save the output report

DirSaveOutput = '/Users/arturo/Dropbox/Research/SoftwareResearch/SNANA/Odyssey/home_snana/lowz/tests_ok/LOWZ_JRK07/'

# DirSaveOutput = DirTables_2

NotebookName ='Compare_SNANA__v10_35g_VS_v10_58g__v1_0.ipynb'
#-----------------------------------------------------------------------------80
# 43 SNe in this table 1.
# 43 SNe in this table 2.

# ## Automatic

import numpy as np
from matplotlib import pyplot as plt
5+6

get_ipython().run_cell_magic('javascript', '', 'var kernel = IPython.notebook.kernel;\nvar thename = window.document.getElementById("notebook_name").innerHTML;\nvar command = "NotebookName = " + "\'"+thename+".ipynb"+"\'";\nkernel.execute(command);')

print(NotebookName)
# 11_DistanceMu_HubbleDiagram_v2_17.ipynb

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

# ### Main loop

# Limits in the tolerance to report a difference between both values.

z_tol = 0.001; zERR_tol = 0.001;
x0_tol = 0.001; x0ERR_tol = 0.0005;
c_tol = 0.008; cERR_tol = 0.009;
x1_tol = 0.05; x1ERR_tol = 0.05;
PKMJD_tol = 2; PKMJDERR_tol = 2;
mB_tol = 0.01; mBERR_tol = 0.01;
COV_x0_x1_tol = 1e-6; COV_x0_c_tol = 1e-6; COV_x1_c_tol = 1e-4;
CHI2_tol = 3; NDOF_tol = 1; FITPROB_tol = 0.01;

text_file = open (DirSaveOutput+'Comparison_SALT2params_SNANA_v10_35g_-_v10_58g_.txt','w')

text_01 = '# Comparison between SNANA versions v10_35g vs v10_58g \n'
text_02 = '# Fitres 1 (v10_35g): %s. \n'%DataFile_1
text_03 = '# Fitres 2 (v10_58g): %s. \n'%DataFile_2
text_04 = '# Location table 1: %s. \n'%DirTables_1
text_05 = '# Location table 2: %s. \n'%DirTables_2

text_line = '#'+'-'*50 + '\n'
text_Author = '# Data table created by: Arturo Avelino. \n'
now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Date   = '# On date: %s \n'%text_timenow
text_script = '# Script used: %s \n'%NotebookName

text_file.write(text_01); text_file.write(text_02);
text_file.write(text_03); text_file.write(text_04);
text_file.write(text_05);

text_file.write(text_line);
text_file.write(text_Author);
text_file.write(text_Date);
text_file.write(text_script);
text_file.write(text_line);

text_file.write("# sn           z_1      zHD_2    z1-z2   1   zERR_1  zHDERR_2  \
zE1-ZE2 1   x0_1     x0_2    x01-x02  1   x0ERR_1  x0ERR_2  x0E1-x0E2 1  c_1       \
c_2      c1-c2   1   cERR_1   cERR_2  cE1-cE2  1   x1_1     x1_2    x1_1-x1_2 1   \
x1ERR_1  x1ERR_2 x1E1-x1E2 1  PKMJD_1     PKMJD_2   PKMJD1-PKMJD2 1    PKMJDERR_1 \
PKMJDERR_2  PKE1-PKE2   1  mB_1     mB_2      mB1-mB2 1   mBERR_1  mBERR_2 \
mBE1-mBE_2 1   COV_x0_x1_1     COV_x0_x1_2       C01_1-C01_2     1   COV_x0_c_1      \
COV_x0_c_2        Cx0c1-Cx0c2     1   COV_x1_c_1      COV_x1_c_2        Cx1c1-Cx1c2     \
1    CHI2_1     CHI2_2  CHI21-CHI22 1  NDOF_1 NDOF_2 N1-N2 1  FITPROB_1  FITPROB_2   \
FP1-FP2   #  Notes \n")

print text_line, text_01, text_02, text_03, text_04, text_05,
text_line, text_Author, text_Date, text_script, text_line

count_same = 0
# Loop over SNe in the fitres file
for i1 in range(len(fitresFile_1['f1'])):
    sn_1 = fitresFile_1['f1'][i1];
    z_1 = fitresFile_1['f2'][i1]; zERR_1 = fitresFile_1['f3'][i1]
    x0_1 = fitresFile_1['f4'][i1]; x0ERR_1 = fitresFile_1['f5'][i1];
    c_1 = fitresFile_1['f6'][i1]; cERR_1 = fitresFile_1['f7'][i1];
    x1_1 = fitresFile_1['f8'][i1]; x1ERR_1 = fitresFile_1['f9'][i1];
    PKMJD_1 = fitresFile_1['f10'][i1]; PKMJDERR_1 = fitresFile_1['f11'][i1];
    mB_1 = fitresFile_1['f12'][i1]; mBERR_1 = fitresFile_1['f13'][i1];
    COV_x0_x1_1 = fitresFile_1['f14'][i1]; COV_x0_c_1 = fitresFile_1['f15'][i1];
    COV_x1_c_1 = fitresFile_1['f16'][i1];
    # print COV_x0_x1_1, COV_x0_c_1, COV_x1_c_1
    CHI2_1 = fitresFile_1['f17'][i1]; NDOF_1 = fitresFile_1['f18'][i1];
    FITPROB_1 = fitresFile_1['f19'][i1];
    # print CHI2_1, NDOF_1, FITPROB_1

    #-------------------

    for i2 in range(len(fitresFile_2['f1'])):

        sn_2 = fitresFile_2['f1'][i2];

        if sn_2 == sn_1:

            print "%12s  %s"%(sn_1, sn_2)
            # Reset the parameters
            note_z=''; note_zErr=''; note_x0=''; note_x0ERR='';
            note_c=''; note_cERR=''; note_x1=''; note_x1ERR='';
            note_PKMJD=''; note_PKMJDERR=''; note_mB=''; note_mBERR='';
            note_COV_x0_x1=''; note_COV_x0_c=''; note_COV_x1_c='';
            note_CHI2=''; note_NDOF=''; note_FITPROB='';

            zHD_2 = fitresFile_2['f10'][i2]; zHDERR_2 = fitresFile_2['f11'][i2];
            x0_2 = fitresFile_2['f27'][i2]; x0ERR_2 = fitresFile_2['f28'][i2];
            c_2 = fitresFile_2['f23'][i2]; cERR_2 = fitresFile_2['f24'][i2];
            x1_2 = fitresFile_2['f21'][i2]; x1ERR_2 = fitresFile_2['f22'][i2];
            PKMJD_2 = fitresFile_2['f19'][i2]; PKMJDERR_2 = fitresFile_2['f20'][i2];
            mB_2 = fitresFile_2['f25'][i2]; mBERR_2 = fitresFile_2['f26'][i2];
            COV_x0_x1_2 = fitresFile_2['f30'][i2]; COV_x0_c_2 = fitresFile_2['f31'][i2];
            COV_x1_c_2 = fitresFile_2['f29'][i2];
            # print COV_x0_x1_2, COV_x0_c_2, COV_x1_c_2
            CHI2_2 = fitresFile_2['f33'][i2]; NDOF_2 = fitresFile_2['f32'][i2];
            FITPROB_2 = fitresFile_2['f34'][i2];
            # print CHI2_2, NDOF_2, FITPROB_2

            #---------------------------------------------------

            if abs(z_1-zHD_2) > z_tol: note_z = "z_1-zHD_2=%7.4f | "%(z_1-zHD_2)
            if abs(zERR_1-zHDERR_2) > zERR_tol: note_zErr = "zERR_1-zHDERR_2=%7.4f | "%(zERR_1-zHDERR_2)
            if abs(x0_1-x0_2) > x0_tol: note_x0 = "x0_1-x0_2=%7.4f | "%(x0_1-x0_2)
            if abs(x0ERR_1-x0ERR_2) > x0ERR_tol: note_x0ERR= "x0ERR_1-x0ERR_2=%7.4f | "%(x0ERR_1-x0ERR_2)
            if abs(c_1-c_2) > c_tol: note_c = "c_1-c_2=%7.4f | "%(c_1-c_2)
            if abs(cERR_1-cERR_2) > cERR_tol: note_cERR = "cERR_1-cERR_2=%7.4f | "%(cERR_1-cERR_2)
            if abs(x1_1-x1_2) > x1_tol: note_x1 = "x1_1-x1_2=%7.4f | "%(x1_1-x1_2)
            if abs(x1ERR_1-x1ERR_2) > x1ERR_tol: note_x1ERR = "x1ERR_1-x1ERR_2=%7.4f | "%(x1ERR_1-x1ERR_2)
            if abs(PKMJD_1-PKMJD_2) > PKMJD_tol: note_PKMJD = "PKMJD_1-PKMJD_2=%7.4f | "%(PKMJD_1-PKMJD_2)
            if abs(PKMJDERR_1-PKMJDERR_2) > PKMJDERR_tol: note_PKMJDERR = "PKMJDERR_1-PKMJDERR_2=%7.4f | "%(PKMJDERR_1-PKMJDERR_2)
            if abs(mB_1-mB_2) > mB_tol: note_mB = "mB_1-mB_2=%7.4f | "%(mB_1-mB_2)
            if abs(mBERR_1-mBERR_2) > mBERR_tol: note_mBERR = "mBERR_1-mBERR_2=%7.4f | "%(mBERR_1-mBERR_2)
            if abs(COV_x0_x1_1-COV_x0_x1_2) > COV_x0_x1_tol: note_COV_x0_x1 = "COV_x0_x1_1-COV_x0_x1_2=%12.4e | "%(COV_x0_x1_1-COV_x0_x1_2)
            if abs(COV_x0_c_1-COV_x0_c_2) > COV_x0_c_tol: note_COV_x0_c = "COV_x0_c_1-COV_x0_c_2=%12.4e | "%(COV_x0_c_1-COV_x0_c_2)
            if abs(COV_x1_c_1-COV_x1_c_2) > COV_x1_c_tol: note_COV_x1_c = "COV_x1_c_1-COV_x1_c_2=%12.4e | "%(COV_x1_c_1-COV_x1_c_2)
            if abs(CHI2_1-CHI2_2) > CHI2_tol: note_CHI2 = "CHI2_1-CHI2_2=%7.4f | "%(CHI2_1-CHI2_2)
            if abs(NDOF_1-NDOF_2) > NDOF_tol: note_NDOF = "NDOF_1-NDOF_2=%7.4f | "%(NDOF_1-NDOF_2)
            if abs(FITPROB_1-FITPROB_2) > FITPROB_tol: note_FITPROB = "FITPROB_1-FITPROB_2=%7.4f | "%(FITPROB_1-FITPROB_2)

            text_file.write("%-12s  %7.4f  %7.4f  %7.4f  1  %7.4f  %7.4f  %7.4f  1  %7.4f  %7.4f  %7.4f  1  %7.4f  %7.4f  %7.4f    1  %7.4f  %7.4f  %7.4f  1  %7.4f  %7.4f  %7.4f  1  %7.4f  %7.4f  %7.4f   1  %7.4f  %7.4f  %7.4f   1  %7.4f  %7.4f  %7.4f     1    %7.4f    %7.4f     %7.4f     1  %7.4f  %7.4f  %7.4f  1  %7.4f  %7.4f  %7.4f    1  %12.4e    %12.4e      %12.4e     1  %12.4e    %12.4e      %12.4e     1  %12.4e    %12.4e      %12.4e     1  %9.4f  %9.4f  %9.4f  1  %4.0f  %4.0f  %4.0f    1  %7.4f    %7.4f    %9.6f  # %s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s \n"%(sn_1,
            z_1, zHD_2, z_1-zHD_2, zERR_1, zHDERR_2, zERR_1-zHDERR_2,
            x0_1, x0_2, x0_1-x0_2, x0ERR_1, x0ERR_2, x0ERR_1-x0ERR_2,
            c_1, c_2, c_1-c_2, cERR_1, cERR_2, cERR_1-cERR_2,
            x1_1,x1_2, x1_1-x1_2, x1ERR_1,x1ERR_2,x1ERR_1-x1ERR_2,
            PKMJD_1, PKMJD_2, PKMJD_1-PKMJD_2, PKMJDERR_1,PKMJDERR_2,PKMJDERR_1-PKMJDERR_2,
            mB_1,mB_2, mB_1-mB_2, mBERR_1,mBERR_2, mBERR_1-mBERR_2,
            COV_x0_x1_1, COV_x0_x1_2, COV_x0_x1_1-COV_x0_x1_2,
            COV_x0_c_1,  COV_x0_c_2,  COV_x0_c_1-COV_x0_c_2,
            COV_x1_c_1,  COV_x1_c_2,  COV_x1_c_1-COV_x1_c_2,
            CHI2_1, CHI2_2, CHI2_1-CHI2_2, NDOF_1, NDOF_2, NDOF_1-NDOF_2,
            FITPROB_1, FITPROB_2, FITPROB_1-FITPROB_2,
            note_z, note_zErr, note_x0, note_x0ERR, note_c, note_cERR,
            note_x1, note_x1ERR, note_PKMJD, note_PKMJDERR, note_mB, note_mBERR,
            note_COV_x0_x1, note_COV_x0_c, note_COV_x1_c,
            note_CHI2, note_NDOF, note_FITPROB ))

            count_same += 1


text_10 = '# %s names are in both tables. \n'%count_same
text_file.write(text_line);
text_file.write(text_10)
print text_10

text_file.close();

text_file.close(); text_file.close();

