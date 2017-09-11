#include "cholesky.cc"

// least_square_mask_acspo
//
// input:
// clear_paths - vector to hold the paths of resulting clear_paths
// original_paths - vector of paths to original granules
// l2p_mask - 2d matrix of l2p flags as masks
//
// This function reads the SST and ACSPO mask from original_paths
// applies the ACSPO mask to the SST
// then computes least square decomposition on the masked SST
// the resulting curve is used as a mask and lsd is run again 
// and the resulting curve is used as a mask 
void
least_square_mask_acspo(vector<string> &clear_paths,const vector<string> &original_paths, const Mat1b &l2p_mask)
{
	int j,y,x,t;
	int time_size = original_paths.size();
	int dims[3] = {HEIGHT, WIDTH, time_size};

	string variable_name = "sea_surface_temperature";
	string filename, clearpath;
	Mat1f sst_samples(3,dims);
	Mat1f ls_approx(3,dims);
	Mat1f sst_a(3,dims); // adjusted sst after first iteration of lsq

	Mat1b clear_mask(3,dims);

	for(j = 0; j < time_size; ++j){

        read_acspo(original_paths[j],clear_mask,j);  // open acspo granule
        readgranule_oneband(original_paths[j], sst_samples, j, variable_name);  // open original granule
        //apply_mask_slice(clear_mask,sst_samples,j,true); // apply mask
        //mask_l2p(clear_mask,l2p_mask,j);

        filename = generate_filename(original_paths[j]);
        clearpath = "data/clear" + filename;
        clear_paths.push_back(clearpath.c_str()); // append save path to clear_paths

        printf("opened file: %s\n",clear_paths[j].c_str());
    }


    // compute nn and diag mask
    // since they use i-1 and i+1 pixels in time we start 
    // the counter from 1 to avoid falling off the matrix
    for(j = DIAG_LAG; j < time_size - DIAG_LAG; ++j ){
        compute_nnmask(sst_samples,clear_mask,j);
        compute_diagmask(sst_samples,clear_mask,j);
        printf("computed NN mask for file %d/%d\n",j,time_size - DIAG_LAG);
    }


    least_square_decomposition(sst_samples, ls_approx, l2p_mask, GAMMA, time_size, true); // compute lsq

    printf("finished first least squares\n");
    
    for(t = 0; t < time_size; ++t){

    	apply_mask_slice(clear_mask,sst_samples,t,true); // apply mask
        mask_l2p(clear_mask,l2p_mask,t);

		for(y = 0; y < HEIGHT; ++y){
			for(x = 0; x < WIDTH; ++ x){

				sst_a(y,x,t) = sst_samples(y,x,t);

				if(l2p_mask(y,x) == 0){
					
					// adjust sst_a 
					if(std::isnan(sst_samples(y,x,t)) || std::isnan(ls_approx(y,x,t)) || sst_samples(y,x,t) < ls_approx(y,x,t) - T_LS){
						sst_a(y,x,t) = NAN; 
					}
					else if(std::isfinite(sst_samples(y,x,t)) && std::isfinite(ls_approx(y,x,t)) && sst_samples(y,x,t) < ls_approx(y,x,t)){
						sst_a(y,x,t) = ls_approx(y,x,t);						
					}
					
					// set LS mask
					if(std::isnan(sst_samples(y,x,t)) || std::isnan(ls_approx(y,x,t)) || sst_samples(y,x,t) < ls_approx(y,x,t) - DELTA_SST){
						clear_mask(y,x,t) = clear_mask(y,x,t) | LS;
					}
					
				}
			}
		}
		//save_test_nc_float(clear,clear_paths[t].c_str());
		printf("masked file: %s\n",clear_paths[t].c_str());
	}
	
	
	
	printf("masked first leastsquares\n");
	least_square_decomposition(sst_a, ls_approx, l2p_mask, GAMMA, time_size, true); //perform second iterations of lsq

	printf("finished second least_square_decomposition\n");
    for(t = 0; t < time_size; ++t){
		for(y = 0; y < HEIGHT; ++y){
			for(x = 0; x < WIDTH; ++ x){
				if(l2p_mask(y,x) == 0){
					// set LS mask 
					if(std::isnan(sst_samples(y,x,t)) || std::isnan(ls_approx(y,x,t)) || sst_samples(y,x,t) < ls_approx(y,x,t) - T_LS){
						//TODO increase histogram for this pixel
						//clear(y,x) = clear(y,x) | LS_APPROX;
						clear_mask(y,x,t) = clear_mask(y,x,t) | LS;
					}
				}

			}
		}
	}
	
	// save masks	
	save_mat(clear_paths, clear_mask, "brightness_temperature_11um2",true,HEIGHT, WIDTH);

}