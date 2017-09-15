import os
import sys
from utils import get_crop

if len(sys.argv) < 3:
    print "usage: python make_directory.py <output_folder> <brightness_temp> <crop_list>"
    print "brightness temp = 0 run for generated data"
    print "brightness temp = 1 run for brightness temps"
    sys.exit()

save_date = sys.argv[1]
data_folder = "../data"
top_folder = os.path.join(data_folder,save_date)

crop_path = sys.argv[3]
crop_list = get_crop(crop_path)


if sys.argv[2] == '0':
    data = ['original_sst', 'acspo', 'approx_ls', 'clear', 'reinstated_ls','collated_mat2_ls', 'scollated_ls', 'approx_ls']
    data = ['original_sst', 'acspo']
else:
    data = ['brightness_temperature_08um6', 'brightness_temperature_10um4','brightness_temperature_11um2','brightness_temperature_12um3']

data_folder = os.path.join("../data",top_folder)

if not os.path.exists(top_folder):
    os.makedirs(top_folder)


for crop in crop_list:
    xs = crop[1]
    ys = crop[0]

    max_y = ys.stop
    min_y = ys.start
    max_x = xs.stop
    min_x = xs.start

    crop_folder = '_'.join(map(str,[min_y,max_y,min_x,max_x]))
    print crop_folder

    folder2 = '/'.join([top_folder,crop_folder])
    if not os.path.exists(folder2):
        os.makedirs(folder2)

    
    for d in data:
        folder3 = '/'.join([folder2,d])
        if not os.path.exists(folder3):
            os.makedirs(folder3)
    