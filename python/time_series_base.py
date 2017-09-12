"""
This program generates a time series of sst, reference and quality level
flags.  The reference is generated via sst+dt_analysis

Note the ploting functionality rarly changes, look into making a carry over function
"""

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import sys


"""
put these functions into utils file?
"""
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


if len(sys.argv) < 3:
    print "usage: python time_series.py <point_file> <granule_list>"
    sys.exit()

point_file = sys.argv[1]
l2p_granule_list_file = sys.argv[2]

# absolute paths to granules
l2p_granule_files = [line.rstrip('\n') for line in open(l2p_granule_list_file)]

# filenames of granules
l2p_granule_filenames = [x.split('/')[-1] for x in l2p_granule_files]


points = read_points(point_file) #pixels to generate time series for
viirs_points = [(1040,1380),(1025,1415),(1025,1440),(915,1333),(920,1380)]

times = []

# produce x label values for each hour
# create function for these lines
hour_ind = []
hour_labels = []
x_lims = []

for i,filename in enumerate(l2p_granule_filenames):
    times.append(filename[8:12])
    if filename[10:12] == '00':
        hour_labels.append(filename[8:12])
        hour_ind.append(i)
    if filename[8:12] == '0000':
        x_lims.append(i)

x_lims = tuple(x_lims)

#dictionaries where var[key] = time series of value
sst_time_series = {}

quality_flags_time_series = {}
bt3_test_time_series = {}
#bt08_time_series = {}
#bt10_time_series = {}
#bt11_time_series = {}
#bt12_time_series = {}
reference_time_series = {}
hour_time_series = {}
viirs_time_series = {}
viirs_inds = {}

for i in points:    
    sst_time_series[i] = []
    quality_flags_time_series[i] = []
    reference_time_series[i] = []
    hour_time_series[i] = []
    bt3_test_time_series[i] = []
    viirs_time_series[i] = []
    viirs_inds[i] = []
    #bt08_time_series[i] = []
    #bt10_time_series[i] = []
    #bt11_time_series[i] = []
    #bt12_time_series[i] = []

viirs_inds[(1040,1380)] = [36,43]
viirs_inds[(1025,1415)] = [36, 43, 81]
viirs_inds[(1025,1440)] = [36,43,81, 87]
viirs_inds[(915,1333)] = [36,43,81, 87]
viirs_inds[(920,1380)] = [36,43,81, 87]

viirs_time_series[(1040,1380)] = [291.29, 291.43]
viirs_time_series[(1025,1415)] = [290.94, 291.21, 291.16]
viirs_time_series[(1025,1440)] = [291.23, 291.36, 292.13, 293.58]
viirs_time_series[(915,1333)] = [288.33, 288.58, 288.96, 289.38]
viirs_time_series[(920,1380)] = [290.33, 289.30, 291.14, 291.51]

# open each granule and append all given pixel values to respective dictionary list
total_files = len(l2p_granule_files)
for i in range(total_files):
    
    # open granlue and get matrix data
    filename = l2p_granule_files[i] 

    cdf = netCDF4.Dataset(filename)
    sst_mat = read_var(cdf,'sea_surface_temperature')
    bt08_mat = read_var(cdf,'brightness_temperature_08um6')
    #bt10_mat = read_var(cdf,'brightness_temperature_10um4')
    bt11_mat = read_var(cdf,'brightness_temperature_11um2')
    bt12_mat = read_var(cdf,'brightness_temperature_12um3')
    dt_mat = read_var(cdf,'dt_analysis')
    quality_level_mat = read_var(cdf,'quality_level')

    # get reference using sst + dt
    reference_mat = sst_mat - dt_mat
    bt3_mat = bt08_mat + 0.8*bt11_mat - (1+0.8)*bt12_mat

    for point in points:

        row, col = point

        sst_time_series[point].append(sst_mat[row,col])  
     
        #bt08_time_series[point].append(bt08_mat[row,col])
        #bt10_time_series[point].append(bt10_mat[row,col])
        #bt11_time_series[point].append(bt11_mat[row,col])
        #bt12_time_series[point].append(bt12_mat[row,col])

        reference_time_series[point].append(reference_mat[row,col])

        if(quality_level_mat[row,col] == 5):
            quality_flags_time_series[point].append(sst_time_series[point][i])
        else:
        	quality_flags_time_series[point].append(np.NAN)

        if(bt3_mat[row,col] < -1):
            bt3_test_time_series[point].append(sst_time_series[point][i])
        else:
            bt3_test_time_series[point].append(np.NAN)
        

    print "finished " + filename

x_inds = range(total_files)

cdf = netCDF4.Dataset(l2p_granule_files[0])
lons_mat = read_var(cdf,'lon')

#create plot
for i,point in enumerate(points):
    
    row,col = point
    # get difference in utc
    lon = lons_mat[row,col]
    utc = round(lon/15)

    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(point))
    plt.title("2017 - 08 - 31   Location: " + str(point) +"\nUTC " + str(utc) , fontsize=20)
    plt.grid()

    plt.plot(x_inds,sst_time_series[point], 'b', label="Original SST")
    #plt.plot(x_inds,bt08_time_series[point], 'r', label="bt08")
    #plt.plot(x_inds,bt10_time_series[point], 'c', label="bt10")
    #plt.plot(x_inds,bt11_time_series[point], 'k', label="bt11")
    #plt.plot(x_inds,bt12_time_series[point], 'm', label="bt12")
    plt.plot(x_inds,sst_time_series[point], 'b.')
    plt.plot(x_inds,quality_flags_time_series[point], 'ro', markersize=5, label="Clear Sky")
    plt.plot(x_inds,bt3_test_time_series[point], 'go', fillstyle='none', markersize=10, markeredgewidth=2, label="BT3 masked")
    plt.plot(x_inds,reference_time_series[point], '-g', lw=3, label="Reference")

    if point in viirs_points:
        plt.plot(viirs_inds[point],viirs_time_series[point],'ko',markersize=10,label="VIIRS")

    plt.legend(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(hour_ind, hour_labels, rotation=-45, fontsize=15)
    plt.show()

"""
pixels = [(3791,972),(1604,3290),(1006,1969),(1472,4105),(845,2894)]

for pixel in pixels:
    outfile = "time_series_" + str(pixel[0]) + "_" + str(pixel[1])
    sio.savemat(outfile,{ "sst":sst[pixel],"acspo":acspo[pixel],"eigen":eig[pixel],"nn":NN[pixel],"diag":diag[pixel],"ls":ls[pixel],
        "ls_full":full[pixel],"ls_acspo":ls_acspo[pixel],"ls_acspo2":ls_acspo2[pixel],"times":times })
"""