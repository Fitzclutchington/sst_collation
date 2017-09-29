import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import sys
import utils

if len(sys.argv) < 4:
    print "usage: python time_series.py <point_file> <reference_file> <file_list>"
    sys.exit()

point_file = sys.argv[1]
ref_file = sys.argv[2]
file_list_path = sys.argv[3]

clear_files = sorted(glob.glob("../data/clear/*.nc"))
local_min_files = sorted(glob.glob("../data/local_min/*.nc"))
original_files = [line.rstrip('\n') for line in open(file_list_path)]

orig_filenames = [x.split('/')[-1] for x in original_files]
clear_filenames = [x.split('/')[-1] for x in clear_files]


original_filenames = [x.split('/')[-1] for x in original_files]



def read_points(point_file):
    points = []
    with open(point_file, 'r') as f:
        for line in f:
            point = map(int,line.split(','))
            points.append((point[0],point[1]))
    return points

points = read_points(point_file)
times = []

hour_ind = []
x_lims = []

for i,clear_file in enumerate(clear_filenames):
    if clear_file[10:12] == '00':
        times.append(clear_file[8:12])
        hour_ind.append(i)
    if clear_file[8:12] == '0000':
        x_lims.append(i)

x_lims = tuple(x_lims)

sst = {}
clear = {}
local_min = {}
refs = {}
h = {}


for i in points:    
    sst[i] = []
    clear[i] = []
    local_min[i] = []
    refs[i] = []
    h[i] = []


total_files = len(local_min_files)
for i in range(total_files):

    base_file = local_min_files[i].split('/')[-1]
    
    
    filename = "../data/clear/"+base_file

    clearbinarync=netCDF4.Dataset(filename)
    clear_binary = utils.read_var(filename,'brightness_temperature_11um2')
    
    mask = clear_binary != 0

    filename = original_files[orig_filenames.index(base_file)]
    sst_val = utils.read_var(filename,'sea_surface_temperature')
    
    filename = local_min_files[i]
    local_min_vals = utils.read_var(filename,'sea_surface_temperature')
    local_min_mask = local_min_vals != 0

    print "finished" , filename , str(i+1) , "/", str(total_files)

    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])  


        if(local_min_mask[j[0],j[1]]):
            local_min[j].append(sst[j][i])
        else:
            local_min[j].append(np.nan)

        if(mask[j[0],j[1]]):
            clear[j].append(sst[j][i])
        else:
        	clear[j].append(np.nan)

        
        
       

sst_ref = utils.read_var(ref_file,"sst_reynolds")
for j in points:
    refs[j] = np.full(total_files,sst_ref[j[0],j[1]])


x_inds = range(total_files)
x_inds_approx = range(92)
titles = ["Masked", "Clear"]
for j,i in enumerate(points):
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title("2017 - 04 - 02   Location: " + str(i), fontsize=20)
    plt.grid()

    plt.plot(x_inds,sst[i], 'b', label="Original SST")

    plt.plot(x_inds,local_min[i], 'or', fillstyle='none', markersize=20, label='local min')
    plt.plot(x_inds, clear[i], 'g.', markersize=20, label='clear')


    plt.plot(x_inds,refs[i], '-g', lw=3, label="reference")
    plt.plot(x_inds,sst[i], 'k.')
    plt.legend(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(hour_ind, times, rotation=-45, fontsize=15)
    plt.show()

"""
pixels = [(3791,972),(1604,3290),(1006,1969),(1472,4105),(845,2894)]

for pixel in pixels:
    outfile = "time_series_" + str(pixel[0]) + "_" + str(pixel[1])
    sio.savemat(outfile,{ "sst":sst[pixel],"acspo":acspo[pixel],"eigen":eig[pixel],"nn":NN[pixel],"diag":diag[pixel],"ls":ls[pixel],
        "ls_full":full[pixel],"ls_acspo":ls_acspo[pixel],"ls_acspo2":ls_acspo2[pixel],"times":times })
"""