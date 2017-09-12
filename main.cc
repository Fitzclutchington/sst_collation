#include <netcdf.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdarg.h>
#include <cmath>
#include <iostream>
#include <fstream> 
#include <Eigen/Dense>
#include <chrono>

using namespace std;
using namespace cv;
using namespace Eigen;
#include "hermite.h"
#include "parameters.cc"
#include "io.cc"
#include "calc.cc"
#include "mask.cc"
#include "filter.cc"
#include "smooth_no_padding.cc"
#include "interpolate.cc"
#include "generate_histogram.cc"
#include "collated_fit_no_padding.cc"
#include "least_square_mask.cc"
#include "least_square_collate.cc"
#include "lsq_test.cc"

#define NDEBUG


int
main(int argc, char *argv[])
{
    if(argc < 3){
        eprintf("Usage: ./ahil2c <granule_list> <reference_file>");
    }

     // vectors used to pass file names between functions
    vector<string> original_paths;
    vector<string> clear_paths;
    vector<string> hourly_paths;
    
    get_file_names(argv[1],original_paths); // list of paths to original granules
    string reference_file = argv[2];   // path to reference file

    // Initialize arrays to hold various l2p masks
    Mat1b land_mask(HEIGHT,WIDTH);
    Mat1b border_mask(HEIGHT,WIDTH);
    Mat1b l2p_mask(HEIGHT,WIDTH);  

    //function to get land mask and invalid(corners) mask
    get_l2pmask(original_paths[0].c_str(),land_mask,l2p_mask);

    //get_lats(original_paths[0], lats);
    get_landborders(land_mask, border_mask, LAND_KERNEL);
    
    printf("read land mask\n");
    //release memory of unused masks
    land_mask.release();

    
    string filename;
    string clearpath;

    /*
    //for debugging purposes
    int i;
    
    for(i =0;i<203;i++){
        filename = generate_filename(original_paths[i]);
        clearpath = "data/clear" + filename;
        clear_paths.push_back(clearpath.c_str());
    }
    */
    least_square_test(clear_paths, original_paths, l2p_mask);
    
    
    printf("starting cloud filter\n");      
    
    
    //least_square_mask_acspo(clear_paths,original_paths, l2p_mask);
    
    
    filter_clouds(original_paths, l2p_mask, border_mask, clear_paths);
    
    
    int dims[3] = { HEIGHT_HIST, WIDTH_HIST, DEPTH_HIST};
    Mat1f histogram (3,dims);
    Mat1b hist_mask(3,dims);
    
    histogram_3d(clear_paths, original_paths, reference_file,0, histogram);
    
    //open_hist_3d("histogram.nc", histogram);
    //SAVENC(histogram);
    mask_histogram(histogram, hist_mask, false);
    
    
    //Mat1b hist_mask(1,1);
    least_squares_collate(clear_paths, original_paths, l2p_mask, hist_mask, reference_file, hourly_paths);
    histogram_3d(hourly_paths,hourly_paths, reference_file, 2, histogram);
    SAVENC(histogram);
    
}
