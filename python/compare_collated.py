import os
import netCDF4
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

if len(sys.argv) < 3:
    print "usage: python diff_reference.py <original_file> <collated> <collated_acspo>"
    print "time as yyyy-mm-dd"
    sys.exit()


original_file_path = sys.argv[1]
collated_file = sys.argv[2]
collated_acspo_file = sys.argv[3]

cdf = netCDF4.Dataset(original_file_path)

l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

sst_original = read_var(cdf,'sea_surface_temperature')

cdf = netCDF4.Dataset(collated_file)
sst_collated = read_var(cdf,'sea_surface_temperature')

cdf = netCDF4.Dataset(collated_acspo_file)
sst_collated_acspo = read_var(cdf,'sea_surface_temperature')

sst_acspo = sst_original.copy()
sst_acspo[cloud_mask] = np.nan

plt.figure()

ax1 = plt.subplot(221)
img1 = ax1.imshow(sst_original, vmin=270, vmax=307)
ax1.imshow(land,interpolation='nearest')
ax1.set_title("Original SST")
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(222, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_collated, vmin=270, vmax=307)
ax2.imshow(land,interpolation='nearest')
ax2.set_title("Collated SST")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(223, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_collated_acspo, vmin=270, vmax=307)
ax2.imshow(land,interpolation='nearest')
ax2.set_title("Collated Acspo SST")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(224, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_acspo, vmin=270, vmax=307)
ax2.imshow(land,interpolation='nearest')
ax2.set_title("Collated - Collated ACSPO")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

plt.show()