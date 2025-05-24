"""
Testing weird nick
"""
import h5py
import pylab as plt
import pandas as pd
import numpy as np


infile = '../../Sting-SED-eagle-rr14_00.hdf5'
file = h5py.File(infile)
mags = file['SED']['ap_dust']['total'][8]
cut = np.where(file['SED']['ap_dust']['total'][11] < 21.2)
infile = "../../mocksky_00.hdf5"
file = h5py.File(infile)
redshift = file['galaxies']['zcos'][:]

plt.scatter(redshift[cut], mags[cut], s=0.01, color='k', alpha=0.2)
plt.ylim(0, 40)
plt.xlim(0, 1)
plt.show()


# Testing the mock catalog generation
infile = '../../waves_wide_gals.parquet'
df = pd.read_parquet(infile)
plt.scatter(df['zcos'], df['total_ap_dust_u_VST'], s=0.01, color='k', alpha=0.2)
plt.ylim(0, 40)
plt.xlim(0, 1)
plt.show()

plt.hist(df['zcos'], bins=500)
plt.yscale('log')
plt.show()
