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
collated_files = sorted(glob.glob("../data/collated_hour/*.nc"))
collated_files = sorted(glob.glob("../data/collated2_hour/*.nc"))
approx_files = sorted(glob.glob("../data/approx_ls/*.nc"))
reinstated_files = sorted(glob.glob("../data/reinstated_ls/*.nc"))
smooth_files = sorted(glob.glob("../data/scollated_ls/*.nc"))
original_files = [line.rstrip('\n') for line in open(file_list_path)]

orig_filenames = [x.split('/')[-1] for x in original_files]
clear_filenames = [x.split('/')[-1] for x in clear_files]
collated_filenames = [x.split('/')[-1] for x in collated_files]

original_filenames = [x.split('/')[-1] for x in original_files]
approx_filenames = [x.split('/')[-1] for x in approx_files]
collated_folder = "../data/collated_hour/"
collated2_folder = "../data/collated2_hour/"

NN_bit = 2;
DIAG_bit = 4;
EIGEN_bit= 8;
LS_bit = 16;
ACSPO_bit = 32;
BT_RATIO_bit = 64

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
smooth = {}
NN = {}
eig = {}
acspo = {}
diag = {}
ls = {}
bt_ratio = {}
polar = {}
collated = {}
collated2 = {}
approx = {}
reinstated = {}
refs = {}
h = {}


for i in points:    
    sst[i] = []
    smooth[i] = []
    NN[i] = []
    eig[i] = []
    acspo[i] = []
    polar[i] = []
    diag[i] = []
    ls[i] = []
    refs[i] = []
    bt_ratio[i] = []
    collated[i] = []
    collated2[i] = []
    reinstated[i] = []
    approx[i] = []
    h[i] = []


total_files = len(clear_files)
for i in range(total_files):

    base_file = clear_files[i].split('/')[-1]
    
    
    filename = clear_files[i]

    clearbinarync=netCDF4.Dataset(filename)
    clear_binary = utils.read_var(filename,'brightness_temperature_11um2')
    


    NN_mask = np.bitwise_and(clear_binary,NN_bit).astype(bool)
    eig_mask = np.bitwise_and(clear_binary,EIGEN_bit).astype(bool)
    diag_mask = np.bitwise_and(clear_binary,DIAG_bit).astype(bool)
    acspo_mask = np.bitwise_and(clear_binary,ACSPO_bit).astype(bool)
    ls_mask = np.bitwise_and(clear_binary,LS_bit).astype(bool)
    bt_ratio_mask = np.bitwise_and(clear_binary,BT_RATIO_bit).astype(bool)
    polar_mask = clear_binary > 32

    filename = original_files[orig_filenames.index(base_file)]
    sst_val = utils.read_var(filename,'sea_surface_temperature')
    
    if base_file in collated_filenames:
        collated_vals = utils.read_var(collated_folder+base_file,'sea_surface_temperature')
        collated2_vals = utils.read_var(collated2_folder+base_file,'sea_surface_temperature')

    if base_file in approx_filenames:
        smooth_vals = utils.read_var(smooth_files[i], 'sea_surface_temperature')
        approx_vals = utils.read_var(approx_files[i], 'sea_surface_temperature')
        reinstated_vals = utils.read_var(reinstated_files[i], 'sea_surface_temperature')

    print "finished" , filename , str(i+1) , "/", str(total_files)

    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])      
        if base_file in approx_filenames:
            approx[j].append(approx_vals[j[0],j[1]])
            reinstated[j].append(reinstated_vals[j[0],j[1]])
            smooth[j].append(smooth_vals[j[0],j[1]])

        if base_file in collated_filenames:
            collated[j].append(collated_vals[j[0],j[1]])
            collated2[j].append(collated2_vals[j[0],j[1]])

        if(NN_mask[j[0],j[1]]):
            NN[j].append(sst[j][i])
        else:
        	NN[j].append(np.nan)

        if(eig_mask[j[0],j[1]]):
            eig[j].append(sst[j][i])
        else:
            eig[j].append(np.nan)

        if(ls_mask[j[0],j[1]]):
            ls[j].append(sst[j][i])
        else:
            ls[j].append(np.nan)

        if(diag_mask[j[0],j[1]]):
            diag[j].append(sst[j][i])
        else:
            diag[j].append(np.nan)

        if(acspo_mask[j[0],j[1]]):
            acspo[j].append(sst[j][i])
        else:
            acspo[j].append(np.nan)

        if(bt_ratio_mask[j[0],j[1]]):
            bt_ratio[j].append(sst[j][i])
        else:
            bt_ratio[j].append(np.nan)
        
        if(polar_mask[j[0],j[1]]):
            polar[j].append(sst[j][i])
        else:
            polar[j].append(np.nan)

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

    plt.plot(x_inds,NN[i], 'r*', markersize=25, label='NN')
    plt.plot(x_inds,eig[i], 'g^', markersize=20, label='eig')
    plt.plot(x_inds,acspo[i], 'bv', markersize=20, label='acspo')

    plt.plot(x_inds,diag[i], 'c.', markersize=20, label='diag')
    plt.plot(x_inds,ls[i], '+', c="g", markersize=20, label='ts2')
    plt.plot(x_inds,bt_ratio[i], 'x', c="r", markersize=20, label='ratio')
    plt.plot(x_inds,polar[i], 'x', c="g", markersize=20, label='polar')

    plt.plot(x_inds,refs[i], '-g', lw=3, label="reference")

    plt.plot(x_inds_approx,approx[i], '-m', lw=3, label='approx')
    plt.plot(x_inds_approx,smooth[i], '-', c='#FFA500', lw=3, label='smooth')
    plt.plot(x_inds_approx,reinstated[i], 'or', markersize=10, fillstyle='none', label='reinstated')
    plt.plot(hour_ind,collated[i], '-k', lw=3, label='collated')
    plt.plot(hour_ind,collated2[i], '-r', lw=3, label='collated above')
    plt.plot(hour_ind,collated2[i], '.k', label='collated above')
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