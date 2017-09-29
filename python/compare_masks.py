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


mask_list_path = sys.argv[2]
l2p_file_path = sys.argv[1]

l2p_filename = l2p_file_path.split('/')[-1]

clear_folder = "../data/clear/"

l3_clear_path = clear_folder + l2p_filename

smooth_folder = "../data/scollated_ls/"

l3_smooth_path = smooth_folder + l2p_filename

approx_folder = "../data/approx_ls/"

l3_approx_path = approx_folder + l2p_filename


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


pass1_mask = utils.read_var(l3_clear_path,'brightness_temperature_11um2')
sst_pass1 = sst_original.copy()
sst_pass1[pass1_mask!=0] = np.nan

smooth = utils.read_var(l3_smooth_path,'sea_surface_temperature')

approx = utils.read_var(l3_approx_path,'sea_surface_temperature')

mask_locs = utils.get_mask_list(mask_list_path)
titles = ['Eigen', 'Laplacian', 'Magdiff Bt12', 'Mag Ratio', 'Median diff', 'Standard Deviation Diff']
vmins = [0,   0, -0.1,   0,   0,   0]
vmaxs = [0.5, 0.8,  0.5, 0.5, 0.5, 0.5]
vmin = 270
vmax = 307



plt.figure()

ax1 = plt.subplot(221)
img1 = ax1.imshow(sst_original, vmin=vmin, vmax=vmax, cmap='jet')
ax1.imshow(land,interpolation='nearest')
ax1.set_title("Original SST")
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)



ax2 = plt.subplot(222, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(sst_pass1, vmin=vmin, vmax=vmax, cmap='jet')
ax2.imshow(land,interpolation='nearest')
ax2.set_title("Pass 1 Mask")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(223, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(smooth, vmin=vmin, vmax=vmax, cmap='jet')
ax2.imshow(land,interpolation='nearest')
ax2.set_title("Smooth")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(224, sharex=ax1,sharey=ax1)
img1 = ax2.imshow(approx, vmin=vmin, vmax=vmax, cmap='jet')
ax2.imshow(land,interpolation='nearest')
ax2.set_title("Approx")
div1 = make_axes_locatable(ax2)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

plt.show()