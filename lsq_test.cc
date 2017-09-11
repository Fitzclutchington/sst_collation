void
least_square_test(vector<string> &clear_paths,const vector<string> &original_paths, const Mat1b &l2p_mask)
{
	int j,y,x,t;
	int time_size = original_paths.size();
	int dims[3] = {HEIGHT, WIDTH, time_size};

	vector<string> ls_full_paths;
	vector<string> ls_acspo_paths;
	vector<string> ls_acspo2_paths;

	string variable_name = "sea_surface_temperature";
	string filename, ls_full_path, ls_acspo_path, ls_acspo2_path;
	Mat1f sst_samples(3,dims);
	Mat1f ls_approx(3,dims);
	Mat1f sst_a(3,dims);

	Mat1b clear_mask(3,dims);

	for(j = 0; j < time_size; ++j){
        // read granule one band
        //read_mask_binary(clear_paths[j],clear_mask,j);S
        //read_mask(clear_paths[j],clear_mask,j);
        read_acspo(original_paths[j],clear_mask,j);
        readgranule_oneband(original_paths[j], sst_samples, j, variable_name);
        //apply_mask_slice(clear_mask,sst_samples,j,true);

        filename = generate_filename(original_paths[j]);
        ls_full_path = "data/ls_full" + filename;
        ls_acspo_path = "data/ls_acspo" + filename;
        ls_acspo2_path = "data/ls_acspo2" + filename;

        ls_full_paths.push_back(ls_full_path);
        ls_acspo_paths.push_back(ls_acspo_path);
        ls_acspo2_paths.push_back(ls_acspo2_path);
        printf("opened file: %s\n",ls_acspo_paths[j].c_str());
    }

    //interpolate_hermite(sst_samples,l2p_mask,time_size,9000,9000);
    least_square_decomposition(sst_samples, ls_approx, l2p_mask, GAMMA, time_size, true);
    save_mat(ls_full_paths, ls_approx, "least_squares",true,HEIGHT, WIDTH);

    for(j = 0; j < time_size; ++j){
    	apply_mask_slice(clear_mask,sst_samples,j,true);
    }

    least_square_decomposition(sst_samples, ls_approx, l2p_mask, GAMMA, time_size, true);
	save_mat(ls_acspo_paths, ls_approx, "least_squares",true,HEIGHT, WIDTH);

    printf("finished first least squares\n");
    
    for(t = 0; t < time_size; ++t){
		for(y = 0; y < HEIGHT; ++y){
			for(x = 0; x < WIDTH; ++ x){
				//clear(y,x) = clear_mask(y,x,t);
				
				//clear(y,x) = clear(y,x) & (~LS_APPROX1);
				sst_a(y,x,t) = sst_samples(y,x,t);

				if(l2p_mask(y,x) == 0){
					
					if(std::isnan(sst_samples(y,x,t)) || std::isnan(ls_approx(y,x,t)) || sst_samples(y,x,t) < ls_approx(y,x,t) - T_LS){
						sst_a(y,x,t) = NAN; 
					}
					else if(std::isfinite(sst_samples(y,x,t)) && std::isfinite(ls_approx(y,x,t)) && sst_samples(y,x,t) < ls_approx(y,x,t)){
						sst_a(y,x,t) = ls_approx(y,x,t);						
					}
					
					if(std::isnan(sst_samples(y,x,t)) || std::isnan(ls_approx(y,x,t)) || sst_samples(y,x,t) < ls_approx(y,x,t) - DELTA_SST){
						clear_mask(y,x,t) = clear_mask(y,x,t) | LS;
					}
					
				}
			}
		}
		//save_test_nc_float(clear,clear_paths[t].c_str());
		printf("masked file: %s\n",ls_full_paths[t].c_str());
	}
	
	
	
	printf("masked first leastsquares\n");
	least_square_decomposition(sst_a, ls_approx, l2p_mask, GAMMA, time_size, true);

	printf("finished second least_square_decomposition\n");
    for(t = 0; t < time_size; ++t){
		for(y = 0; y < HEIGHT; ++y){
			for(x = 0; x < WIDTH; ++ x){
				if(l2p_mask(y,x) == 0){

					if(std::isnan(sst_samples(y,x,t)) || std::isnan(ls_approx(y,x,t)) || sst_samples(y,x,t) < ls_approx(y,x,t) - T_LS){
						//TODO increase histogram for this pixel
						//clear(y,x) = clear(y,x) | LS_APPROX;
						clear_mask(y,x,t) = clear_mask(y,x,t) | LS;
					}
				}
			}
		}
	}
	
	save_mat(ls_acspo2_paths, ls_approx, "least_squares",true,HEIGHT, WIDTH);
}