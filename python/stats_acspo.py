#REFERENCE_SST IS DATE OF REFERENCE

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
    print "usage: python stats_clean.py <reference> <time>"
    print "time as yyyy-mm-dd"
    sys.exit()

clear_files = sorted(glob.glob("../data/clear/*.nc"))
approx_files = sorted(glob.glob("../data/approx_ls/*.nc"))
collated2_files = sorted(glob.glob("../data/collated_mat2_ls/*.nc"))
scollate_files = sorted(glob.glob("../data/scollated_ls/*.nc"))
r_collate_files = sorted(glob.glob("../data/above_ls/*.nc"))
reinstated_files = sorted(glob.glob("../data/reinstated_ls/*.nc"))
acspo_files = sorted(glob.glob("../../acspo_test/data/collated_acspo/*.nc"))
original_files = [line.rstrip('\n') for line in open('../ahitest.txt')]
orig_filenames = [x.split('/')[-1] for x in original_files]
approx_filenames = [x.split('/')[-1] for x in approx_files]
clear_filenames = [x.split('/')[-1] for x in clear_files]


collated2_std = []
scollated_std = []
acspo_std = []
pass1_std = []
approx_std = []
reinstated_std = []
col_acspo_std = []

collated2_mean = []
scollated_mean = []
acspo_mean = []
pass1_mean = []
approx_mean = []
reinstated_mean = []
col_acspo_mean = []

acspo_clear = []
approx_clear = []
collated2_clear = []
scollated_clear = []
pass1_clear = []
reinstated_clear = []
col_acspo_clear = []

ref_file = sys.argv[1]
time_file = sys.argv[2]
outfile = 'means_' + time_file
print outfile


cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")
sza = read_var(cdf,"satellite_zenith_angle")
mask = np.abs(sza) < 67

ref_sst[~mask] = np.nan

total_files = len(approx_files)

times =[]
time_stamps = []
hour_ind = []

for i,approx_file in enumerate(approx_filenames):
    time_stamps.append(approx_file[8:12])
    if approx_file[10:12] == '00':
        times.append(approx_file[8:12])
        hour_ind.append(i)

for i in range(total_files):

    base_file = approx_files[i].split('/')[-1]

    cdf = netCDF4.Dataset(approx_files[i])
    vals = read_var(cdf,"sea_surface_temperature")   
    cdf.close()
    diff = vals-ref_sst
    approx_std.append(np.nanstd(diff))
    approx_mean.append(np.nanmean(diff))
    approx_clear.append(np.isfinite(vals).sum())

    cdf = netCDF4.Dataset(reinstated_files[i])
    vals = read_var(cdf,"sea_surface_temperature")   
    cdf.close()
    diff = vals-ref_sst
    reinstated_std.append(np.nanstd(diff))
    reinstated_mean.append(np.nanmean(diff))
    reinstated_clear.append(np.isfinite(vals).sum())

    cdf = netCDF4.Dataset(acspo_files[i])
    vals = read_var(cdf,"sea_surface_temperature")   
    cdf.close()
    diff = vals-ref_sst
    col_acspo_std.append(np.nanstd(diff))
    col_acspo_mean.append(np.nanmean(diff))
    col_acspo_clear.append(np.isfinite(vals).sum())

    sst_filename = original_files[orig_filenames.index(base_file)]
    cdf = netCDF4.Dataset(sst_filename)
    l2p_flags = read_var(cdf,'l2p_flags')
    sst_clear = read_var(cdf, 'sea_surface_temperature')
    
    sst_acspo = sst_clear.copy()
    sst_acspo[~mask] = np.nan
    cdf.close()

    
    cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
    sst_acspo[cloud_mask] = np.nan    
    acspo_std.append(np.nanstd(sst_acspo - ref_sst))
    acspo_mean.append(np.nanmean(sst_acspo-ref_sst))
    acspo_clear.append(np.isfinite(sst_acspo).sum())
    


    clear_filename = clear_files[clear_filenames.index(base_file)]
    cdf = netCDF4.Dataset(clear_filename)
    clear = read_var(cdf,"brightness_temperature_11um2")
    cloud_mask = clear == 0
    
    sst_clear[~mask] = np.nan
    sst_clear[~cloud_mask] = np.nan
    diff = sst_clear-ref_sst
    pass1_std.append(np.nanstd(diff))
    pass1_mean.append(np.nanmean(diff))
    pass1_clear.append(np.isfinite(sst_clear).sum())
   


    cdf = netCDF4.Dataset(collated2_files[i])
    col = read_var(cdf,"sea_surface_temperature")
    col[~mask] = np.nan
    
    diff = col - ref_sst
    collated2_std.append(np.nanstd(diff))
    collated2_mean.append(np.nanmean(diff))
    collated2_clear.append(np.isfinite(col).sum())
    cdf.close()

    cdf = netCDF4.Dataset(scollate_files[i])
    col = read_var(cdf,"sea_surface_temperature")
    col[~mask] = np.nan
    
    diff = col - ref_sst
    scollated_std.append(np.nanstd(diff))
    scollated_mean.append(np.nanmean(diff))
    scollated_clear.append(np.isfinite(col).sum())
    cdf.close() 
    

    print "finished file " + base_file + " " + str(i) + "/" + str(total_files)

"""
outfile = '../data/' +time_file + "/stds"
sio.savemat(outfile, {"sst":sst_std, "pass2":pass2_std, "collated1":collated, "collated2":collated2, "scollated":scollated,"reinstated":reinstated_std, "approx1":approx_std, "approx2":approx2_std, "time_stamps":time_stamps})

outfile = '../data/' +time_file + "/means"
sio.savemat(outfile,{ 'sst':sst_mean,'pass2':pass2_mean,'collated1':collated_mean,'collated2':collated2_mean,'scollated':scollated_mean,'reinstated':reinstated_mean,'approx1':approx_mean,'approx2':approx2_mean,'time_stamps':time_stamps})

outfile = '../data/' +time_file + "/num_obs"
sio.savemat(-outfile,{ "sst":acspo,"pass2":clear_pass2,"collated1":clear_col,"collated2":clear_col2,"scollated":clear_scol,"reinstated":clear_reinstated,"approx1":clear_approx,"approx2":clear_approx2, "time_stamps":time_stamps})
"""

plt.figure()
plt.grid()
plt.plot(collated2_std, c='r', lw=3, label="collated2")
plt.plot(scollated_std, c='m', label="scollated")
plt.plot(acspo_std, c='k', lw=2,label="acspo")
plt.plot(pass1_std, c='c',label="pass1")
plt.plot(approx_std, c="#FFA500", label="approx")
plt.plot(reinstated_std, c="g", label="reinstated")
plt.plot(col_acspo_std, c="g", lw=3, label="collated acspo")

plt.legend()
plt.xticks(hour_ind,times,rotation='vertical')
plt.title("Standard Deviation " + time_file)
plt.show()

#plt.savefig("stds.png")

plt.figure()
plt.grid()
plt.plot(collated2_mean, c='r', lw=3, label="collated2")
plt.plot(scollated_mean, c='m', label="scollated")
plt.plot(acspo_mean, c='k', lw=2,label="acspo")
plt.plot(pass1_mean, c='c',label="pass1")
plt.plot(approx_mean, c="#FFA500", label="approx")
plt.plot(reinstated_mean, c="g", label="reinstated")
plt.plot(col_acspo_mean, c="g", lw=3, label="collated acspo")
plt.legend()
plt.xticks(hour_ind,times,rotation='vertical')
plt.title("Means " + time_file)
plt.show()
#plt.savefig("means.png")


plt.figure()
plt.grid()
plt.plot(collated2_clear, c='r', lw=3, label="collated2")
plt.plot(scollated_clear, c='m', label="scollated")
plt.plot(acspo_clear, c='k', lw=2,label="acspo")
plt.plot(pass1_clear, c='c',label="pass1")
plt.plot(approx_clear, c="#FFA500", label="approx")
plt.plot(reinstated_clear, c="g", label="reinstated")
plt.plot(col_acspo_clear, c="g", lw=3, label="collated acspo")
plt.legend()
plt.title("Number Clear " + time_file)
plt.xticks(hour_ind,times,rotation='vertical')
plt.show()
    
#plt.savefig("num_obs.png")
#outfile = '../data/' +time_file + "/num_obs"
#sio.savemat(outfile,{ "sst":acspo,"pass2":clear_pass2,"collated1":clear_col,"collated2":clear_col2,"scollated":clear_scol,"reinstated":clear_reinstated,"approx1":clear_approx,"approx2":clear_approx2, "time_stamps":time_stamps})


