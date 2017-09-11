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

original_file_path = sys.argv[1]

cdf = netCDF4.Dataset(original_file_path)

l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
cloud = np.zeros((5500,5500,4))
r = 173/256.0
g = 216/256.0
b = 230/256.0
cloud[cloud_mask] = [r,g,b,1]

sst = read_var(cdf,'sea_surface_temperature')

granule_name = original_file_path.split('/')[-1]
loc_collated = "../data/collated_mat2_ls/"

cdf = netCDF4.Dataset(loc_collated+granule_name)
collated = read_var(cdf,'sea_surface_temperature')

loc_clear = "../data/clear/"
loc_smooth = "../data/scollated_ls/"
loc_reinstated = '../data/reinstated_ls/'

cdf = netCDF4.Dataset(loc_clear+granule_name)
mask = read_var(cdf,'brightness_temperature_11um2')

cdf = netCDF4.Dataset(loc_reinstated+granule_name)
reinstated = read_var(cdf,'sea_surface_temperature')

cdf = netCDF4.Dataset(loc_smooth+granule_name)
smooth = read_var(cdf,'sea_surface_temperature')

#cdf = netCDF4.Dataset(loc_approx+granule_name)
#approx = read_var(cdf,'sea_surface_temperature')

clear_sst = sst.copy()
clear_sst[mask != 0] = np.nan

acspo_sst = sst.copy()
acspo_sst[cloud_mask] =  np.nan

plt.figure()
"""
vmin = 296
vmax = 306
"""

vmin = 273
vmax = 282

inds = np.s_[600:741,2550:2751]

ax1 = plt.subplot(321)
img1 = ax1.imshow(sst[inds],vmin=vmin,vmax=vmax)
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land[inds],interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(322, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(reinstated[inds],vmin=vmin,vmax=vmax)
#ax1.imshow(cloud,interpolation='nearest')
ax2.imshow(land[inds],interpolation='nearest')
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(323, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(acspo_sst[inds],vmin=vmin,vmax=vmax)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land[inds],interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(324, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(smooth[inds],vmin=vmin,vmax=vmax)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land[inds],interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(325, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(clear_sst[inds],vmin=vmin,vmax=vmax)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land[inds],interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(326, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(collated[inds],vmin=vmin,vmax=vmax)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land[inds],interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

plt.show()
