import netCDF4
import numpy as np
import glob
import sys
from utils import read_var, get_crop

if len(sys.argv) < 3:
	print "usage: python save_crop.py <output_folder> <file_list> <crops_list>"
	print "date as yyyy-mm-dd"
	



output_folder = sys.argv[1].strip('/')

file_list = sys.argv[2]
crop_path = sys.argv[3]

crops = get_crop(crop_path)

data_folder = '../data'
output_folder = '/'.join([data_folder,output_folder])

x_lims = []
original_files = [line.rstrip('\n') for line in open(file_list)]
orig_filenames = [x.split('/')[-1] for x in original_files]


original_folder = 'original_sst'
acspo_folder = 'acspo'
clear_folder = 'clear'
reinstated_folder = 'reinstated_ls'
collated_folder = 'collated_mat2_ls'
scollated_folder = 'scollated_ls'
approx_folder = 'approx_ls'


for crop in crops:

	ys = crop[0]
	xs = crop[1]
	crop_folder = '_'.join(map(str,[ys.start,ys.stop,xs.start,xs.stop]))
	print crop_folder
	w = xs.stop-xs.start
	h = ys.stop-ys.start

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

		data = read_var(f,variable_name)[crop]
		l2p_flags = read_var(f,l2p_name)[crop]
		l2p_flags = np.bitwise_and(l2p_flags,-16384).astype(bool)


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
	