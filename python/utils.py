import numpy as np
import netCDF4

def read_var(filepath, variable):
    cdf = netCDF4.Dataset(filepath)
    data = np.squeeze(cdf[variable][:])
    cdf.close()
    return data

def get_crop(filepath):
    crops = []
    with open(filepath, 'r') as f:
        for line in f:
            y_crop, x_crop = line.split(',')
            min_y, max_y = map(int, y_crop.split(':'))
            min_x, max_x = map(int, x_crop.split(':'))
            crops.append(np.s_[min_y:max_y+1,min_x:max_x+1])
    return crops

def get_mask_list(filepath):
    mask_paths = []
    with open(filepath, 'r') as f:
        for line in f:
            mask_paths.append(line.strip())
    return mask_paths