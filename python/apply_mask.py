import sys
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import cmocean

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

if len(sys.argv) < 2:
    print "usage: python test_output.py <original_file>"
    sys.exit()

l2p_file_path = sys.argv[1]

l2p_filename = l2p_file_path.split('/')[-1]

clear_folder = "../data/clear/"

l3_clear_path = clear_folder + l2p_filename

cdf = netCDF4.Dataset(l2p_file_path)
sst = read_var(cdf, 'sea_surface_temperature')

l2p_flags = read_var(cdf,'l2p_flags')
cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((land_mask.shape[0],land_mask.shape[1],4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

cdf = netCDF4.Dataset(l3_clear_path)
clear_mask = read_var(cdf, 'brightness_temperature_11um2')
clear_mask = clear_mask != 0

removed = (cloud_mask == 0) & (clear_mask)

sst_masked = sst.copy()
sst_masked[clear_mask] = np.nan
sst[cloud_mask] = np.nan


plt.figure()

ax1 = plt.subplot(131)
img1 = ax1.imshow(sst,vmin=270,vmax=307, cmap='jet')
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(132, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_masked,vmin=270,vmax=307,cmap='jet')
#ax1.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(133, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(removed,cmap='gray')
#ax1.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

plt.show()