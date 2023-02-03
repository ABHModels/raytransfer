# Python script to take calculated photons from ray-tracing code, calculate the transfer function, and generate the FITS file
# Modify deformation parameter values, isco.dat file name, and colb column name to use for other deformation parameters

import numpy as np
import math
import itertools as it
import scipy.interpolate
from scipy import interpolate
import subprocess
from astropy.io import fits
from subprocess import Popen, PIPE
import sys
import os
import glob

# Inclination angle values - mu0 is cos(inc)
mu0 = [0.0349447653, 0.09718278, 0.15948, 0.2165542, 0.270481, 0.3221819, 0.3721757, 0.420793, 0.4682622, 0.5147499, 0.5603828, 0.6052601, 0.6494616, 0.6930526, 0.7360878, 0.7786132, 0.8206683, 0.8622873, 0.9035001, 0.9443328, 0.9848086238, 0.9986296296]

# Spin values
a = [-0.998, -0.8, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.5804822, 0.6381477, 0.6800805, 0.7142378, 0.7435757, 0.7695684, 0.7930734, 0.8146393, 0.8346417, 0.8533508, 0.8709681, 0.8876485, 0.9035142, 0.9186634, 0.9331764, 0.9471198, 0.9605495, 0.9735132, 0.9860516, 0.9982]

# Accretion rate - disk thickness
Mdl = 0.0

# Deformation parameter values
epsi3 = 0.0
#a13 = 0.0
a22 = 0.0
a52 = 0.0

# Initialize the FITS file
prihdr = fits.Header()
prihdr['INFO'] = "Transfer functions for a non-Kerr metric"
prihdr['AUTHORS'] = 'A. B. Abdikamalov, D. Ayzenberg, C. Bambi,  S. Nampalliwar, '
prihdr['COMMENTS'] = "Makes a file containing tabulated transfer functions."
prihdu = fits.PrimaryHDU(header=prihdr)
aux=[prihdu]

fits_defpar=np.zeros([30,30])

# Read in isco.dat file containing the spin, defpar grid
spin, defpar_d, isco = np.loadtxt('isco_a13.dat', unpack=True)

# Set 30 by 30 spin, defpar array
k = 0
for i in np.arange(0,30):
    for j in np.arange(0,30):
        fits_defpar[i][j] = defpar_d[k]
        k+=1

fits_defpar.sort()

# Set columns for FITS file
cola = fits.Column(name='a', format='1E',array=a)
colb = fits.Column(name='a13', format='30E',array=fits_defpar)
cols = fits.ColDefs([cola, colb])
tbhdu = fits.BinTableHDU.from_columns(cols)
aux.append(tbhdu)

# Store the inclination values used
colmu0 = fits.Column(name='mu0', format='1E', array=mu0)
cols = fits.ColDefs([colmu0])
tbhdu = fits.BinTableHDU.from_columns(cols)
aux.append(tbhdu)

# Loops over all spins and inclinations

k = 0

for spin in a:
    defpar_d = fits_defpar[k]
    k = k+1
    for defpar in defpar_d:
        for i_obs in mu0:
            spin=np.float(spin)

            ## Read data file and calculate TRF
            filename = 'photons/photons4trf_a%.05Le.i%.02Le.Mdl%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat'%(spin,i_obs,Mdl,epsi3,a13,a22,a52)
            with open(filename, 'rb+') as f:
                f.seek(-1,2)
                lchar = f.read()
                if(lchar!="\n"):
                    print (filename)

            re, g, xobs, yobs, cosem, gminus, xobsminus, yobsminus, gplus, xobsplus, yobsplus = np.loadtxt(filename, unpack=True)

            fits_radii=np.zeros([100,1])
            fits_gmins=np.zeros([100,1])
            fits_gmaxs=np.zeros([100,1])
            fits_trff1=np.zeros([100,40])
            fits_trff2=np.zeros([100,40])
            fits_cosem1=np.zeros([100,40])
            fits_cosem2=np.zeros([100,40])

            #Loop over various radii
            for i in np.arange(1,101):

                gmax = g[(i+1)*84-1];
                gmin = g[i*84];

                trff1=[]
                trff2=[]
                cosem1=[]
                cosem2=[]

                #Branch 1
                for j in range(1, 41):

					# Set positions of current value and values at higher and lower radius
                    poscur = i*84+j;
                    posrplus = poscur-84;
                    posrminus = poscur+84;

					# Calculate transfer function

                    dxdr = (xobs[posrplus]-xobs[posrminus])/(re[posrplus]-re[posrminus]);
                    dydr = (yobs[posrplus]-yobs[posrminus])/(re[posrplus]-re[posrminus]);
                    dxdg = (xobsplus[poscur]-xobsminus[poscur])/(gplus[poscur]-gminus[poscur]);
                    dydg = (yobsplus[poscur]-yobsminus[poscur])/(gplus[poscur]-gminus[poscur]);
                    jac = -(dxdg*dydr-dxdr*dydg);

                    recur = re[poscur];
                    gcur = g[poscur];

					# In case there are redshift values larger than the maximum
                    if(gmax-gcur < 0):
                        gcur = gmax;

                    trff1.append(gcur*math.sqrt((gmax-gcur)*(gcur-gmin))*jac/(math.pi*recur));
                    cosem1.append(cosem[poscur]);

                #Branch 2
                for j in range(43, 83):

					# Set positions of current value and values at higher and lower radius
                    poscur = i*84+j;
                    posrplus = poscur-84;
                    posrminus = poscur+84;

					# Calculate transfer function

                    dxdr = (xobs[posrplus]-xobs[posrminus])/(re[posrplus]-re[posrminus]);
                    dydr = (yobs[posrplus]-yobs[posrminus])/(re[posrplus]-re[posrminus]);
                    dxdg = (xobsplus[poscur]-xobsminus[poscur])/(gplus[poscur]-gminus[poscur]);
                    dydg = (yobsplus[poscur]-yobsminus[poscur])/(gplus[poscur]-gminus[poscur]);
                    jac = (dxdg*dydr-dxdr*dydg);

                    recur = re[poscur];
                    gcur = g[poscur];

                    # In case there are redshift values larger than the maximum
                    if(gmax-gcur < 0):
                        gcur = gmax;

                    trff2.append(gcur*math.sqrt((gmax-gcur)*(gcur-gmin))*jac/(math.pi*recur));
                    cosem2.append(cosem[poscur]);

                #Organization of the FITS file

                fits_radii[i-1]=re[i*84]
                fits_gmins[i-1]=gmin
                fits_gmaxs[i-1]=gmax
                fits_trff1[i-1]=trff1
                fits_trff2[i-1]=trff2
                fits_cosem1[i-1]=cosem1
                fits_cosem2[i-1]=cosem2

            col1 = fits.Column(name='r', format='1E',array=fits_radii)
            col2 = fits.Column(name='gmin', format='1E',array=fits_gmins)
            col3 = fits.Column(name='gmax', format='1E',array=fits_gmaxs)
            col4 = fits.Column(name='trff1', format='40E',array=fits_trff1)
            col5 = fits.Column(name='trff2', format='40E',array=fits_trff2)
            col6 = fits.Column(name='cosne1', format='40E',array=fits_cosem1)
            col7 = fits.Column(name='cosne2', format='40E',array=fits_cosem2)

            cols = fits.ColDefs([col1, col2,col3, col4,col5, col6,col7])
            tbhdu = fits.BinTableHDU.from_columns(cols)

            #Append each case
            aux.append(tbhdu)

# Write the FITS file
thdulist = fits.HDUList(aux)
thdulist.writeto('Trf.fits')
thdulist.close()
thdulist = fits.HDUList(aux)
