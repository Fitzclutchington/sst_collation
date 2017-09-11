import netCDF4
import numpy as np
import glob
import sys

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

if len(sys.argv) < 3:
	print "usage: python save_crop.py <output_folder> <file_list>"
	print "date as yyyy-mm-dd"
	


min_ys = [  940,  1800, 1300]
max_ys = [  1081, 2001, 1801]
min_xs = [  2000, 4500, 1200]
max_xs = [  2201, 4801, 2201]

min_ys = [  3760  ]
max_ys = [  3860 ]
min_xs = [  900 ]
max_xs = [  1000 ]

output_folder = sys.argv[1].strip('/')

file_list = sys.argv[2]
data_folder = '../data/'
output_folder = '/'.join([data_folder,output_folder])

x_lims = []
original_files = [line.rstrip('\n') for line in open(file_list)]
orig_filenames = [x.split('/')[-1] for x in original_files]

for i,f in enumerate(orig_filenames):
    if f[8:12] == '0000':
        x_lims.append(i)

filenames = orig_filenames[x_lims[0]:x_lims[1]+1]
original_files = original_files[x_lims[0]:x_lims[1]+1]
orig_filenames = orig_filenames[x_lims[0]:x_lims[1]+1]

original_folder = 'original_sst'
acspo_folder = 'acspo'
clear_folder = 'clear'
reinstated_folder = 'reinstated_ls'
collated_folder = 'collated_mat2_ls'
scollated_folder = 'scollated_ls'
approx_folder = 'approx_ls'


for j in range(len(max_xs)):
	max_y = max_ys[j]
	min_y = min_ys[j]
	max_x = max_xs[j]
	min_x = min_xs[j]

	crop_folder = '_'.join(map(str,[min_y,max_y,min_x,max_x]))
	print crop_folder
	w = max_x-min_x
	h = max_y-min_y

	save_folder = '/'.join([output_folder, crop_folder])
	"""
	folders = [clear_folder]
	variable_name = "brightness_temperature_11um2"
	for folder in folders:
		files = ['/'.join([folder, f]) for f in filenames]
		for i,f in enumerate(files):

			f_name = data_folder + f
			cdf = netCDF4.Dataset(f_name)
			data = read_var(cdf,variable_name)[min_y:max_y,min_x:max_x]
			cdf.close()
			
			save_loc = save_folder+ "/" + f
			rootgrp = netCDF4.Dataset(save_loc, "w", format="NETCDF4")
			height = rootgrp.createDimension("height", h)
			width = rootgrp.createDimension("width", w)

			var = rootgrp.createVariable(variable_name,"u1",("height","width"))
			var[:] = data
			rootgrp.close()
			print f

	
	# nothing needs to be done to these data, they are all in a format to be cropped
	folders = [collated_folder, scollated_folder,  approx_folder] 
	variable_name = 'sea_surface_temperature'
	for folder in folders:
		files = ['/'.join([folder, f]) for f in filenames]
		for i,f in enumerate(files):

			f_name = data_folder + f
			cdf = netCDF4.Dataset(f_name)
			data = read_var(cdf,variable_name)[min_y:max_y,min_x:max_x]
			cdf.close()
			
			save_loc = save_folder+ "/" + f
			rootgrp = netCDF4.Dataset(save_loc, "w", format="NETCDF4")
			height = rootgrp.createDimension("height", h)
			width = rootgrp.createDimension("width", w)

			var = rootgrp.createVariable(variable_name,"f4",("height","width"))
			var[:] = data
			rootgrp.close()
			print f

	# open reinstated and convert to boolean flag
	variable_name_output = 'mask'
	files = ['/'.join([reinstated_folder, f]) for f in filenames]
	for i,f in enumerate(files):

		f_name = data_folder + f
		cdf = netCDF4.Dataset(f_name)
		data = read_var(cdf,variable_name)[min_y:max_y,min_x:max_x]
		cdf.close()
		
		mask = np.zeros((h,w)).astype(np.uint8)
		mask[np.isfinite(data)] = 255

		save_loc = save_folder+ "/" + f
		rootgrp = netCDF4.Dataset(save_loc, "w", format="NETCDF4")
		height = rootgrp.createDimension("height", h)
		width = rootgrp.createDimension("width", w)

		var = rootgrp.createVariable(variable_name_output,"u1",("height","width"))
		var[:] = mask
		rootgrp.close()
		print f
	"""
	# save sst and l2p flags
	variable_name = "sea_surface_temperature"
	l2p_name = "l2p_flags"
	for i,f in enumerate(original_files):

		f_name = orig_filenames[i]

		cdf = netCDF4.Dataset(f)
		data = read_var(cdf,variable_name)[min_y:max_y,min_x:max_x]
		l2p_flags = read_var(cdf,l2p_name)[min_y:max_y,min_x:max_x]
		l2p_flags = np.bitwise_and(l2p_flags,-16384).astype(bool)
		cdf.close()


		save_loc = '/'.join([save_folder ,original_folder , f_name])
		print save_loc
		rootgrp = netCDF4.Dataset(save_loc, "w", format="NETCDF4")
		height = rootgrp.createDimension("height", h)
		width = rootgrp.createDimension("width", w)

		var = rootgrp.createVariable(variable_name,"f4",("height","width"))
		var[:] = data
		rootgrp.close()
		

		mask = np.zeros((h,w)).astype(np.uint8)
		mask[~l2p_flags] = 255

		save_loc = '/'.join([save_folder ,acspo_folder, f_name])
		rootgrp = netCDF4.Dataset(save_loc, "w", format="NETCDF4")
		height = rootgrp.createDimension("height", h)
		width = rootgrp.createDimension("width", w)

		var = rootgrp.createVariable("mask","u1",("height","width"))
		var[:] = mask
		rootgrp.close()
		print save_loc
	