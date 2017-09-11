import os
import sys

if len(sys.argv) < 3:
    print "usage: python make_directory.py <output_folder> <brightness_temp>"
    print "brightness temp = 0 run for generated data"
    print "brightness temp = 1 run for brightness temps"
    sys.exit()

save_date = sys.argv[1]
data_folder = "../data"
top_folder = os.path.join(data_folder,save_date)

min_ys = [  940,  1800, 1300]
max_ys = [  1081, 2001, 1801]
min_xs = [  2000, 4500, 1200]
max_xs = [  2201, 4801, 2201]

min_ys = [  3760  ]
max_ys = [  3860 ]
min_xs = [  900 ]
max_xs = [  1000 ]

if sys.argv[2] == '0':
    data = ['original_sst', 'acspo', 'approx_ls', 'clear', 'reinstated_ls','collated_mat2_ls', 'scollated_ls', 'approx_ls']
    data = ['original_sst', 'acspo']
else:
    data = ['brightness_temperature_08um6', 'brightness_temperature_10um4','brightness_temperature_11um2','brightness_temperature_12um3']

data_folder = os.path.join("../data",top_folder)

if not os.path.exists(top_folder):
    os.makedirs(top_folder)


for j in range(len(min_xs)):
    max_y = max_ys[j]
    min_y = min_ys[j]
    max_x = max_xs[j]
    min_x = min_xs[j]

    crop_folder = '_'.join(map(str,[min_y,max_y,min_x,max_x]))
    print crop_folder

    folder2 = '/'.join([top_folder,crop_folder])
    if not os.path.exists(folder2):
        os.makedirs(folder2)

    
    for d in data:
        folder3 = '/'.join([folder2,d])
        if not os.path.exists(folder3):
            os.makedirs(folder3)
    