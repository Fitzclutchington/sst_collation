void
least_squares_collate(const vector<string> &clear_paths, const vector<string> &original_paths, const Mat1b &l2p_mask, 
                       Mat1b &hist_mask, const string ref_file, vector<string> &collated_hour_paths)
{

	int time_size = clear_paths.size();
    printf("time_size = %d\n",time_size);
    int i;
    int collated_size, collated_interp_size;
    string filename;

    string approx_folder_loc = "data/approx_ls";
    string s_collated_folder_loc = "data/scollated_ls";
    string collated2_folder_loc = "data/collated_mat2_ls";
    string reinstated_folder_loc = "data/reinstated_ls";
    string collated_hour_loc = "data/collated_hour";

    vector<string> approx_paths;
    vector<string> scollated_paths;
    vector<string> collated2_paths;
    vector<string> reinstated_paths;

    vector<int> collated_inds;
    vector<int> interp_inds;

    Mat1f clear_samples, approx, collated, collated_interp, sst;
    Mat1b clear_mask;

    //generate save paths for collated_values as well as indicies of hour
    for(i = 0; i < time_size; ++i){
        filename = generate_filename(clear_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            collated_inds.push_back(i);
        }
    }

    collated_size = collated_inds.size(); // how many hours of colaltion we have
    printf("collated size = %d\n",collated_size);
    // filenames and indicies for the rest of the interval
    // last granule is ignored to account for the removal of large derivatives
    for(i = collated_inds[0]  ; i <= collated_inds[collated_size - 1]; ++i){
        filename = generate_filename(clear_paths[i]);  
        approx_paths.push_back(approx_folder_loc + filename); 
        collated2_paths.push_back(collated2_folder_loc + filename);
        scollated_paths.push_back(s_collated_folder_loc + filename);
        reinstated_paths.push_back(reinstated_folder_loc + filename);
        if(filename[11] == '0' && filename[12] == '0'){
            collated_hour_paths.push_back(collated_hour_loc + filename);
        }
    }


    collated_interp_size = collated2_paths.size();
    printf("collated_interp_size = %d\n",collated_interp_size);

    int dims[3] = {HEIGHT, WIDTH, collated_interp_size};
    int dims_collate[3] = { HEIGHT, WIDTH, collated_size};

    clear_samples.create(3,dims);
    clear_mask.create(3,dims);
    
    collated_interp.create(3,dims);
    sst.create(3,dims);
    

    //MatrixXf Gamma(collated_size,collated_size);

    for(i = 0; i < collated_interp_size; ++i){
        filename = generate_filename(approx_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            interp_inds.push_back(i);
        }
    }


    string variable_name = "sea_surface_temperature";

    //open files
    for(i = 0; i < collated_interp_size; ++i)
    {
        
    	read_mask(clear_paths[i+collated_inds[0]],clear_mask,i);
    	readgranule_oneband(original_paths[i+collated_inds[0]], sst, i, variable_name);
    	readgranule_oneband(original_paths[i+collated_inds[0]], clear_samples, i, variable_name);
        apply_mask_slice(clear_mask, clear_samples, i, true);
        printf("opened file %s\n",clear_paths[i+collated_inds[0]].c_str());
        printf("opened file %s\n",original_paths[i+collated_inds[0]].c_str());
        /*
        readgranule_oneband(reinstated_paths[i], clear_samples, i, variable_name);
        readgranule_oneband(approx_paths[i], approx, i, variable_name);
        printf("opened file %s\n",reinstated_paths[i].c_str());
        */
    }

    
    ///////////////////////////////
    // calculate smooth     ///////
    ///////////////////////////////
    int window[3] = {3, 3, COLLATED_SMOOTH_LAG};
	smooth_samples_collated(clear_samples, collated_interp,l2p_mask, ref_file, window, collated_interp_size);
	printf("finished smoothing of collated values\n");

	//save_mat(scollated_paths, collated_interp, "sea_surface_temperature",true,HEIGHT, WIDTH);
    
    
    /////////////////////////////////////////////////
    // REINSTATE CLEAR WHERE ORIGINAL SST > SMOOTH //
    /////////////////////////////////////////////////
    
    printf("starting reinstated\n");

    Mat1f t_cold(HEIGHT,WIDTH);
    Mat1f t_warm(HEIGHT,WIDTH);

 

    compute_dt_thresholds(t_cold, t_warm, ref_file);
    //SAVENC(t_cold);SAVENC(t_warm);
    compute_reinstated(clear_samples, sst, collated_interp, clear_mask, t_cold, t_warm, collated_interp_size, ref_file);
    clear_mask.release();
    distribution_check_3d(hist_mask, clear_samples, l2p_mask, collated_interp_size,ref_file, reinstated_paths);

    //save_mat(reinstated_paths, clear_samples, "sea_surface_temperature",true,HEIGHT, WIDTH);
    
    ////////////////////////////////
    // calculate approximation    //
    ////////////////////////////////
    approx.create(3,dims);
	calc_approximate(collated_interp,clear_samples, l2p_mask, approx, ref_file, collated_interp_size);
	//save_mat(approx_paths, approx, "sea_surface_temperature",true, HEIGHT, WIDTH);

	//set_gamma(Gamma, GAMMA_S);
    printf("set Gamma\n");
    
    ////////////////////////////////
    // Collate samples            //
    ////////////////////////////////
    collated.create(3,dims_collate);
	printf("Starting Collation on sea_surface_temperature\n");
	collate_samples(clear_samples, approx, collated, interp_inds,l2p_mask,collated_interp_size, ref_file, false, GAMMA2);
	printf("Finished Collation on sea_surface_temperature\n");
	approx.release();
    remove_long_interpolation(clear_samples, collated, interp_inds, collated_interp_size);
    remove_last_value(collated,collated_size);
    save_mat(collated_hour_paths, collated, "sea_surface_temperature", true, HEIGHT, WIDTH);

    /*
	remove_long_interpolation(clear_samples, collated, interp_inds, collated_interp_size);
	interpolate_collation(collated, collated_interp, l2p_mask, interp_inds, collated_interp_size,T_INTERP,2*INTERP_DIST);
	remove_last_value(collated_interp,collated_interp_size);
    //SAVENC(collated_interp);
	save_mat(collated2_paths, collated_interp, "sea_surface_temperature",true, HEIGHT, WIDTH);
    */

}