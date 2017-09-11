import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import sys

if len(sys.argv) < 2:
    print "usage: python time_series.py <point_file>"
    sys.exit()

point_file = sys.argv[1]
ref_file = sys.argv[3]

original_files = [line.rstrip('\n') for line in open('../ahitest.txt')]
clear_files = sorted(glob.glob("../data/clear/*.nc"))
approx_files = sorted(glob.glob("../data/approx_ls/*.nc"))
collated2_files = sorted(glob.glob("../data/collated_mat2_ls/*.nc"))
scollate_files = sorted(glob.glob("../data/scollated_ls/*.nc"))
above_files = sorted(glob.glob("../data/above_ls/*.nc"))
reinstated_files = sorted(glob.glob("../data/reinstated_ls/*.nc"))
acspo_files = sorted(glob.glob("../../acspo_test/data/collated_acspo/*.nc"))

original_files = [line.rstrip('\n') for line in open('../ahitest.txt')]
orig_filenames = [x.split('/')[-1] for x in original_files]
approx_filenames = [x.split('/')[-1] for x in approx_files]
clear_filenames = [x.split('/')[-1] for x in clear_files]



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

#points = [(1738,2880), (1533,1649), (4255,1488),   (3975,1405), (4443,3479), (4516,3858), (1885,4360), (4763,2138),   (1906,507)]
points = read_points(point_file)
#y_lims = [(299,302),   (290,295),   (292.5,295.5), (296,299),   (294,297),   (290,295),   (298,303),   (284.5,287.5), (298,302) ]
times = []

hour_ind = []
x_lims = []

for i,approx_file in enumerate(approx_filenames):
    if approx_file[10:12] == '00':
        times.append(approx_file[8:12])
        hour_ind.append(i)
    if approx_file[8:12] == '0000':
        x_lims.append(i)

x_lims = tuple(x_lims)

sst = {}
acspo = {}
col_acspo = {}
collated2 = {}
approx = {}
pass1 = {}
scollate = {}
reinstated = {}
above = {}
h = {}

for i in points: 
    reinstated[i] = []   
    sst[i] = []
    acspo[i] = []
    col_acspo[i] = []
    collated2[i] = []
    approx[i] = []
    pass1[i] = []
    scollate[i] = []
    above[i] = []
    h[i] = []

total_files = len(collated2_files)
for i in range(total_files):

    base_file = collated2_files[i].split('/')[-1]
        
    collatednc = netCDF4.Dataset(collated2_files[i])
    collated2_vals = read_var(collatednc,'sea_surface_temperature')
    collatednc.close()

    collatednc = netCDF4.Dataset(acspo_files[i])
    col_acspo_vals = read_var(collatednc,'sea_surface_temperature')
    collatednc.close()

    approxnc = netCDF4.Dataset(approx_files[i])
    approx_vals = read_var(approxnc,'sea_surface_temperature')
    approxnc.close()

    scolnc = netCDF4.Dataset(scollate_files[i])
    scollated_vals = read_var(scolnc,'sea_surface_temperature')
    scolnc.close()

    reinstatednc = netCDF4.Dataset(reinstated_files[i])
    reinstated_vals = read_var(reinstatednc,'sea_surface_temperature')
    reinstatednc.close()

    """
    abovenc = netCDF4.Dataset(above_files[i])
    above_vals = read_var(abovenc,'sea_surface_temperature')
    abovenc.close()
    """

    filename = original_files[orig_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    acspo_val = read_var(orignc, 'l2p_flags')
    orignc.close()

    acspo_val = np.bitwise_and(acspo_val,-16384).astype(bool)

    filename = clear_files[clear_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    clear_val = read_var(orignc,'brightness_temperature_11um2')
    orignc.close()


    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])
        
        collated2[j].append(collated2_vals[j[0],j[1]])

        approx[j].append(approx_vals[j[0],j[1]])

        scollate[j].append(scollated_vals[j[0],j[1]])

        reinstated[j].append(reinstated_vals[j[0],j[1]])
        #above[j].append(above_vals[j[0],j[1]])

        col_acspo[j].append(col_acspo_vals[j[0],j[1]])

        if(not acspo_val[j[0],j[1]]):
            acspo[j].append(sst[j][i])
        else:
        	acspo[j].append(np.nan)

        if clear_val[j[0],j[1]] == 0:
            pass1[j].append(sst[j][i])
        else:
            pass1[j].append(np.nan)
    
    print "finished file " + str(i+1) + "/" + str(total_files)
for j in points:
    for i in hour_ind:
        h[j].append(collated2[j][i])

cdf = netCDF4.Dataset(ref_file)
sst_ref = read_var(cdf,"sst_reynolds")

refs = {}
for j in points:
    refs[j] = np.full(total_files,sst_ref[j[0],j[1]])

pixels = [(1472,4105),(1604, 3290),(4221,1392),(4500,1690)]

for pixel in pixels:
    outfile = "time_series_" + str(pixel[0]) + "_" + str(pixel[1])
    sio.savemat(outfile,{ "sst":sst[pixel],"acspo":acspo[pixel],"collated2":collated2[pixel],"scollated":scollate[pixel],"approx":approx[pixel],"pass1":pass1[pixel],"times":times })

x_inds = range(total_files)
for j,i in enumerate(points):
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title("2017 - 01 -08   Location: " + str(i), fontsize=20)
    plt.grid()
    plt.plot(x_inds,sst[i], 'b', label="Original SST")
    plt.plot(x_inds,sst[i], 'b.')
    plt.plot(x_inds,acspo[i], 'bo', fillstyle='none', markeredgewidth=2, markersize=12, label="ACSPO Mask")
    plt.plot(x_inds,collated2[i], '-k', linewidth=3, markersize=15, label="Collated 2")
    plt.plot(x_inds,approx[i], c='r', label="Approx")
    plt.plot(x_inds,scollate[i], c='g', label="Scollate")
    plt.plot(x_inds,reinstated[i], "kx", markersize=20, label="Reinstated")
    plt.plot(x_inds,pass1[i], 'r.',  markersize=10, label="pass1")
    plt.plot(x_inds,col_acspo[i], 'r', lw=3, markersize=10, label="col acspo")
    plt.plot(hour_ind,h[i], '.', c='c', markersize=20, label="Hour")
    plt.plot(x_inds,refs[i], c="g", lw=2, label="Ref")
    plt.xticks(hour_ind, times, rotation=-45, fontsize=15)
    #plt.xlim(x_lims)
    #plt.ylim(y_lims[j])
    plt.legend(fontsize=15)
    plt.yticks(fontsize=15)
    
    plt.show()
