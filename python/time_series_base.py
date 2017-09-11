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
times = []

# produce x label values for each hour
# create function for these lines
hour_ind = []
x_lims = []

for i,filename in enumerate(l2p_granule_filenames):
    if filename[10:12] == '00':
        times.append(filename[8:12])
        hour_ind.append(i)
    if filename[8:12] == '0000':
        x_lims.append(i)

x_lims = tuple(x_lims)

#dictionaries where var[key] = time series of value
sst_time_series = {}
quality_flags_time_series = {}
reference_time_series = {}
hour_time_series = {}

for i in points:    
    sst_time_series[i] = []
    quality_flags_time_series[i] = []
    reference_time_series[i] = []
    hour_time_series[i] = []


# open each granule and append all given pixel values to respective dictionary list
total_files = len(l2p_granule_files)
for i in range(total_files):
    
    # open granlue and get matrix data
    filename = l2p_granule_files[i] 

    cdf = netCDF4.Dataset(filename)
    sst_mat = read_var(cdf,'sea_surface_temperature')
    dt_mat = read_var(cdf,'dt_analysis')
    quality_level_mat = read_var(cdf,'quality_level')

    # get reference using sst + dt
    reference_mat = sst_mat + dt_mat


    for point in points:

        row, col = point

        sst_time_series[point].append(sst_mat[row,col])      
        reference_time_series[point].append(reference_mat[row,col])

        if(quality_level_mat[row,col] == 5):
            quality_flags_time_series[point].append(sst_time_series[point][i])
        else:
        	quality_flags_time_series[point].append(np.NAN)

        

    print "finished " + filename




x_inds = range(total_files)

#create plot
for i,point in enumerate(points):
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title("2017 - 08 - 31   Location: " + str(point) +"\n", fontsize=20)
    plt.grid()

    plt.plot(x_inds,sst_time_series[point], 'b', label="Original SST")
    plt.plot(x_inds,sst_time_series[point], 'b.', label="Original SST")
    plt.plot(x_inds,quality_flags_time_series[point], 'ro', markersize=10, label="Clear Sky")
    plt.plot(x_inds,reference_time_series[point], '-g', lw=3, label="Reference")

    
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