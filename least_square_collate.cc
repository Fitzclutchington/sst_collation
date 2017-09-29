void
least_squares_collate(const vector<string> &clear_paths, const vector<string> &original_paths, const Mat1b &l2p_mask, 
                      const Mat1b &hist_mask, const Mat1f &reference_sst, const Mat1f &sza, const Mat1b &border_mask, vector<string> &collated_hour_paths)
{

	int time_size = clear_paths.size();
    printf("time_size = %d\n",time_size);
    int i;
    int collated_size, collated_interp_size;
    string filename;

    // All folder_loc and _paths variables are used to save intermediate outputs
    string approx_folder_loc = "data/approx_ls";
    string s_collated_folder_loc = "data/smooth_ls";
    string collated2_folder_loc = "data/collated_mat2_ls";
    string reinstated_folder_loc = "data/reinstated_ls";
    string collated_hour_loc = "data/collated_hour";
    //string collated2_hour_loc = "data/collated2_hour";


    vector<string> approx_paths;
    vector<string> scollated_paths;
    vector<string> collated2_paths;
    //vector<string> collated2_hour_paths;
    vector<string> reinstated_paths;
    

    vector<int> collated_inds;
    vector<int> interp_inds;

    Mat1f clear_samples, approx, collated, collated_interp, sst;
    Mat1w clear_mask;

    // collation is performed on every hour so we find all filenames at time 'xx00'
    // as well as their positions in the file paths arrays
    
    for(i = 0; i < time_size; ++i){
        filename = generate_filename(clear_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            collated_inds.push_back(i);
        }
    }

    collated_size = collated_inds.size(); // how many hours of colaltion we have
    printf("collated size = %d\n",collated_size);
    // filenames and indicies for the rest of the interval
    // we start computation of approximation, smooth, and collation at the first time 'xx00'
    // last granule is ignored to account for the removal of large derivatives
    for(i = collated_inds[0]  ; i <= collated_inds[collated_size - 1]; ++i){
        filename = generate_filename(clear_paths[i]);  
        approx_paths.push_back(approx_folder_loc + filename); 
        collated2_paths.push_back(collated2_folder_loc + filename);
        scollated_paths.push_back(s_collated_folder_loc + filename);
        reinstated_paths.push_back(reinstated_folder_loc + filename);

        if(filename[11] == '0' && filename[12] == '0'){
            collated_hour_paths.push_back(collated_hour_loc + filename);
            //collated2_hour_paths.push_back(collated2_hour_loc + filename);
        }
    }


    // size from first hour to last hour in SAMPLING_RATE steps
    collated_interp_size = collated2_paths.size();
    printf("collated_interp_size = %d\n",collated_interp_size);

    int dims[3] = {HEIGHT, WIDTH, collated_interp_size};
    int dims_collate[3] = { HEIGHT, WIDTH, collated_size};

    // matrix stack to hold sst with masks applied
    clear_samples.create(3,dims);

    // clear mask holds the actual mask data
    clear_mask.create(3,dims);
    
    // collated_interp is used as a temporary variable for data that is not used throughout
    // the entire algorithm, such as the smooth data
    collated_interp.create(3,dims);

    // full sst
    sst.create(3,dims);
    

    // interp_inds hold the indicies from 0 - collated_interp_size
    // of all 'xx00' our granules.
    for(i = 0; i < collated_interp_size; ++i){
        filename = generate_filename(approx_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            interp_inds.push_back(i);
        }
    }


    string variable_name = "sea_surface_temperature";

    //open files and load data
    for(i = 0; i < collated_interp_size; ++i)
    {
        
    	read_mask(clear_paths[i+collated_inds[0]],clear_mask,i);
    	readgranule_oneband(original_paths[i+collated_inds[0]], sst, i, variable_name);
    	readgranule_oneband(original_paths[i+collated_inds[0]], clear_samples, i, variable_name);
        apply_mask_slice(clear_mask, clear_samples, i, true);
        printf("opened file %s\n",clear_paths[i+collated_inds[0]].c_str());
        printf("opened file %s\n",original_paths[i+collated_inds[0]].c_str());
 
    }

    //Remove local min, this may be more appropriate in the masking step
    remove_local_min(clear_samples,l2p_mask, collated_interp_size);
    
    
    ///////////////////////////////
    // calculate smooth     ///////
    ///////////////////////////////

    // 3d moving window on masked sst
    int window[3] = {3, 3, COLLATED_SMOOTH_LAG};
	smooth_samples(clear_samples, collated_interp,l2p_mask, reference_sst, window, collated_interp_size);
	printf("finished smoothing of collated values\n");

    // perform the collated step on the smooth data in order to work out any extra kinks
    collated.create(3,dims_collate);
    collate_samples(collated_interp, collated_interp, collated, interp_inds,l2p_mask,collated_interp_size, reference_sst, border_mask, false, GAMMA2);
    interpolate_collation(collated, collated_interp, l2p_mask, interp_inds, collated_interp_size,T_INTERP,2*INTERP_DIST);
    save_mat(scollated_paths, collated_interp, "sea_surface_temperature",true,HEIGHT, WIDTH);

    /////////////////////////////////////////////////
    // REINSTATE CLEAR WHERE ORIGINAL SST > SMOOTH //
    /////////////////////////////////////////////////
    
    printf("starting reinstated\n");

    Mat1f t_cold(HEIGHT,WIDTH);
    Mat1f t_warm(HEIGHT,WIDTH); 

    // dt thresholds are based on sza
    compute_dt_thresholds(t_cold, t_warm, sza);
    
    // try to reinstate masked out pixels based on how close the pixels are to the smooth curve
    compute_reinstated(clear_samples, sst, collated_interp, clear_mask, t_cold, t_warm, collated_interp_size, reference_sst);
    clear_mask.release();
    distribution_check_3d(hist_mask, clear_samples, l2p_mask, collated_interp_size, reference_sst, sza, reinstated_paths, original_paths[0]);
    remove_local_min(clear_samples,l2p_mask, collated_interp_size);
    save_mat(reinstated_paths, clear_samples, "sea_surface_temperature",true,HEIGHT, WIDTH);
    

    ////////////////////////////////
    // calculate approximation    //
    ////////////////////////////////

    // project smooth time series on clear pixels
    approx.create(3,dims);
	calc_approximate(collated_interp,clear_samples, l2p_mask, approx, reference_sst, collated_interp_size);
	save_mat(approx_paths, approx, "sea_surface_temperature",true, HEIGHT, WIDTH);


    
    ////////////////////////////////
    // Collate samples            //
    ////////////////////////////////
    
    // perform collation algorithm that averages the approximated data with the masked sst at every hour
    // returns a matrix collated of size (HEIGHT, WIDTH, num_hours) -> collated
	printf("Starting Collation on sea_surface_temperature\n");
	collate_samples(clear_samples, approx, collated, interp_inds,l2p_mask,collated_interp_size, reference_sst, border_mask, false, GAMMA2);
	printf("Finished Collation on sea_surface_temperature\n");
	approx.release();
    
    // remove any collated pixels where no clear samples existed, might be able to perform this step in the collate_samples function
    remove_long_interpolation(clear_samples, collated, interp_inds, collated_interp_size);

    // This function was necessary when we interpolated the collated values from hourly into SAMPLING_RATE time
    // Mayno longer be necessary
    remove_last_value(collated,collated_size);
    save_mat(collated_hour_paths, collated, "sea_surface_temperature", true, HEIGHT, WIDTH);

    /*
    Mat1f above(3,dims);
	interpolate_collation(collated, collated_interp, l2p_mask, interp_inds, collated_interp_size,T_INTERP,2*INTERP_DIST);
    
    compute_above_samples(clear_samples, collated_interp, above, l2p_mask, collated_interp_size);

    float sigma = 1/6.0;

    collate_sample_weighted_average(above, above, collated, collated_inds, 0, l2p_mask, collated_interp_size, reference_sst, border_mask,sigma );
	save_mat(collated2_hour_paths, collated, "sea_surface_temperature",true, HEIGHT, WIDTH);
    */
    
    
}