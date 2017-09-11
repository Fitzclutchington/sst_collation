import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import sys

if len(sys.argv) < 4:
    print "usage: python time_series.py <point_file> <variable_name> <reference_file>"
    sys.exit()

point_file = sys.argv[1]
var = sys.argv[2]
ref_file = sys.argv[3]

original_files = [line.rstrip('\n') for line in open('../ahitest.txt')]
#original_files = sorted(glob.glob("../data/2016-11-15/*.nc"))
clear_files = sorted(glob.glob("../data/clear/*.nc"))
full_files = sorted(glob.glob("../data/ls_full/*.nc"))
acspo_files = sorted(glob.glob("../data/ls_acspo/*.nc"))
acspo2_files = sorted(glob.glob("../data/ls_acspo2/*.nc"))

original_files = [line.rstrip('\n') for line in open('../ahitest.txt')]

orig_filenames = [x.split('/')[-1] for x in original_files]
clear_filenames = [x.split('/')[-1] for x in clear_files]


original_filenames = [x.split('/')[-1] for x in original_files]

filter_window_lag = 1
second_pass_lag = 25/2
smooth_lag = 19/2
collated_inds = []

NN_bit = 2;
DIAG_bit = 4;
EIGEN_bit= 8;
LS_bit = 16;
ACSPO_bit = 32;

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

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
NN = {}
eig = {}
acspo = {}
diag = {}
ls = {}
full = {}
ls_acspo = {}
ls_acspo2 = {}
refs = {}
h = {}


for i in points:    
    sst[i] = []
    NN[i] = []
    eig[i] = []
    acspo[i] = []
    diag[i] = []
    ls[i] = []
    full[i] = []
    ls_acspo[i] = []
    ls_acspo2[i] = []
    refs[i] = []
    h[i] = []


total_files = len(clear_files)
for i in range(total_files):

    base_file = clear_files[i].split('/')[-1]
    
    
    filename = clear_files[i]

    clearbinarync=netCDF4.Dataset(filename)
    clear_binary = read_var(clearbinarync,'brightness_temperature_11um2')
    clearbinarync.close()


    NN_mask = np.bitwise_and(clear_binary,NN_bit).astype(bool)
    eig_mask = np.bitwise_and(clear_binary,EIGEN_bit).astype(bool)
    diag_mask = np.bitwise_and(clear_binary,DIAG_bit).astype(bool)
    acspo_mask = np.bitwise_and(clear_binary,ACSPO_bit).astype(bool)
    ls_mask = np.bitwise_and(clear_binary,LS_bit).astype(bool)
    filename = original_files[orig_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    orignc.close()

    lsnc = netCDF4.Dataset(full_files[i])
    ls_full_val = read_var(lsnc,"least_squares")

    lsnc = netCDF4.Dataset(acspo_files[i])
    acspo_val = read_var(lsnc,"least_squares")

    lsnc = netCDF4.Dataset(acspo2_files[i])
    acspo2_val = read_var(lsnc,"least_squares")

    print "finished " + filename

    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])      
        full[j].append(ls_full_val[j[0],j[1]])
        ls_acspo[j].append(acspo_val[j[0],j[1]])
        ls_acspo2[j].append(acspo2_val[j[0],j[1]])

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

        

cdf = netCDF4.Dataset(ref_file)
sst_ref = read_var(cdf,"sst_reynolds")
for j in points:
    refs[j] = np.full(total_files,sst_ref[j[0],j[1]])


x_inds = range(total_files)

for j,i in enumerate(points):
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title("2017 - 04 - 02   Location: " + str(i) +"\n" + var, fontsize=20)
    plt.grid()

    plt.plot(x_inds,sst[i], 'b', label="Original SST")

    plt.plot(x_inds,NN[i], 'r*', markersize=25, label='NN')
    plt.plot(x_inds,eig[i], 'g^', markersize=20, label='eig')
    plt.plot(x_inds,acspo[i], 'bv', markersize=20, label='acspo')

    plt.plot(x_inds,diag[i], 'c.', markersize=20, label='diag')
    plt.plot(x_inds,ls[i], '+', c="#FF4500", markersize=20, label='ts2')
    
    plt.plot(x_inds,full[i], 'c', label='full')
    plt.plot(x_inds,ls_acspo[i], 'm', label='ls_acspo')
    plt.plot(x_inds,ls_acspo2[i], 'r', label='ls_acspo2')

    plt.plot(x_inds,refs[i], '-g', lw=3, label="reference")

    plt.plot(x_inds,sst[i], 'k.')
    plt.legend(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(hour_ind, times, rotation=-45, fontsize=15)
    plt.show()

pixels = [(3791,972),(1604,3290),(1006,1969),(1472,4105),(845,2894)]

for pixel in pixels:
    outfile = "time_series_" + str(pixel[0]) + "_" + str(pixel[1])
    sio.savemat(outfile,{ "sst":sst[pixel],"acspo":acspo[pixel],"eigen":eig[pixel],"nn":NN[pixel],"diag":diag[pixel],"ls":ls[pixel],
        "ls_full":full[pixel],"ls_acspo":ls_acspo[pixel],"ls_acspo2":ls_acspo2[pixel],"times":times })