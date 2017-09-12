import os
import sys

import netCDF4
from skimage import measure
from scipy import ndimage, signal
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as cl
from mpl_toolkits.axes_grid1 import make_axes_locatable

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

def gradient(img):
    h = np.array([0.036420, 0.248972, 0.429217, 0.248972, 0.036420])
    hp = np.array([0.108415, 0.280353, 0, -0.280353, -0.108415])

    img = np.pad(img, 2, mode='reflect')
    mode = 'valid'
    dX = signal.convolve(signal.convolve(img, hp[:,np.newaxis], mode=mode),
    h[np.newaxis,:], mode=mode)
    dY = signal.convolve(signal.convolve(img, h[:,np.newaxis], mode=mode),
    hp[np.newaxis,:], mode=mode)
    return dX, dY

if len(sys.argv) < 2:
    print "usage: python compute_gradient.py <granule>"
    sys.exit()

l2p_granule_file = sys.argv[1]

cdf = netCDF4.Dataset(l2p_granule_file)
sst = read_var(cdf,"sea_surface_temperature")
dt_analysis = read_var(cdf,'dt_analysis')
l2p_flags = read_var(cdf,'l2p_flags')

cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((land_mask.shape[0],land_mask.shape[1],4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

acspo_mask = sst.copy()
acspo_mask[cloud_mask] = np.nan

cdf = netCDF4.Dataset('../test_eigen/bt_ratio.nc')
bt3_ratio = read_var(cdf,'data');

cdf = netCDF4.Dataset('../test_eigen/cold_ratio.nc')
bt2_ratio = read_var(cdf,'data')

sst_dX, sst_dY = gradient(sst)
sst_mag = np.sqrt(sst_dX**2 + sst_dY**2)

plt.figure()

ax1 = plt.subplot(221)

img1 = ax1.imshow(bt3_ratio,vmin=-1,vmax=1,cmap='jet')
ax1.set_title("BT3 Ratio")
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(222, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst,vmin=287,vmax=297,cmap='jet')
ax2.set_title("SST")
#ax1.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(223, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(dt_analysis,vmin=-3,vmax=3,cmap='jet')
ax2.set_title("DT Analysis")
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(224, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(acspo_mask,vmin=287,vmax=297,cmap='jet')
ax2.set_title("SST ACSPO Mask")
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

plt.show()