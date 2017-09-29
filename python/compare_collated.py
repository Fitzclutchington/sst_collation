import os
import netCDF4
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import utils


if len(sys.argv) < 2:
    print "usage: python diff_reference.py <original_file>"
    print "time as yyyy-mm-dd"
    sys.exit()


l2p_file_path = sys.argv[1]

l2p_filename = l2p_file_path.split('/')[-1]

clear_folder = "../data/clear/"

l3_clear_path = clear_folder + l2p_filename

collated_folder = '../data/collated_hour/'

l3_collated_path = collated_folder + l2p_filename

reinstated_folder = '../data/reinstated_ls/'

l3_reinstated_path = reinstated_folder + l2p_filename

smooth_folder = '../data/scollated_ls/'

l3_smooth_path = smooth_folder + l2p_filename

l2p_flags = utils.read_var(l2p_file_path,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
land = np.zeros((land_mask.shape[0],land_mask.shape[1],4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]


sst_original = utils.read_var(l2p_file_path,'sea_surface_temperature')
dt_analysis = utils.read_var(l2p_file_path,'dt_analysis')
bt08 = utils.read_var(l2p_file_path,'brightness_temperature_08um6')
bt11 = utils.read_var(l2p_file_path,'brightness_temperature_11um2')
bt12 = utils.read_var(l2p_file_path,'brightness_temperature_12um3')

bt_ratio = bt08 + 0.8*bt11 - 1.8*bt12
sst_acspo = sst_original.copy()
sst_acspo[cloud_mask] = np.nan

pass1_mask = utils.read_var(l3_clear_path,'brightness_temperature_11um2')
sst_pass1 = sst_original.copy()
sst_pass1[pass1_mask!=0] = np.nan

sst_collated = utils.read_var(l3_collated_path,'sea_surface_temperature')
sst_reinstated = utils.read_var(l3_reinstated_path,'sea_surface_temperature')
sst_smooth = utils.read_var(l3_smooth_path,'sea_surface_temperature')

dX,dY = np.gradient(sst_original)
sst_grad = np.sqrt(dX**2 + dY**2)
vmin = 270
vmax = 307
plt.figure()

ax1 = plt.subplot(231)
img1 = ax1.imshow(sst_original, vmin=vmin, vmax=vmax, cmap='jet')
ax1.imshow(land,interpolation='nearest')
ax1.set_title("Original SST")
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(232, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(dt_analysis, vmin=-3, vmax=3, cmap='jet')
ax2.imshow(land,interpolation='nearest')
ax2.set_title("DT Analysis")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(233, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_acspo, vmin=vmin, vmax=vmax, cmap='jet')
ax2.imshow(land,interpolation='nearest')
ax2.set_title("SST ACSPO Mask")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(234, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_reinstated, vmin=vmin, vmax=vmax, cmap='jet')
ax2.imshow(land,interpolation='nearest')
ax2.set_title("SST Reinstated")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(235, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_pass1, vmin=vmin, vmax=vmax, cmap='jet')
ax2.imshow(land,interpolation='nearest')
ax2.set_title("SST Pass1")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(236, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_grad, vmin=0, vmax=0.5, cmap='gray')
ax2.imshow(land,interpolation='nearest')
ax2.set_title("SST Smooth")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

plt.show()