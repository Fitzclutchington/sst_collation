void
generate_sst_histogram(const vector<string> &pass2_files, const Mat1f &pass2, 
	                   const Mat1f &reference, const Mat1f &lon, Mat1f &hist)
{
	int i,x,y,m,n;
	int time_size = pass2_files.size();
	printf("time_size = %d\n",time_size);
	int d_range = 10;
	int mod = SAMPLING_RATE * 24;
	float step = 0.2;
	float hour_stamp_n,minute_stamp_n;
	float dt;
	string filename, minute_stamp, hour_stamp;
	printf("intialized+vars\n");
	hist.setTo(0);
	for(i = 0; i < time_size; ++i){
		filename = generate_filename(pass2_files[i]);
		printf("filename = %s\n",filename.c_str());
		minute_stamp = filename.substr(11,2);
		minute_stamp_n = atoi(minute_stamp.c_str()) / 60.0;
		hour_stamp = filename.substr(9,2) ;

		hour_stamp_n = atoi(hour_stamp.c_str())+minute_stamp_n;
		printf("starting loop\n");
		for(y = 0; y < HEIGHT; ++y){
			for(x = 0; x < WIDTH; ++x){
				if(std::isfinite(pass2(y,x,i)) && std::isfinite(reference(y,x))){
					dt = pass2(y,x,i) - reference(y,x);
					m = floor((dt+d_range)/step);
					n = ((((int) floor(SAMPLING_RATE*(hour_stamp_n + lon(y,x)/15)))  % mod) + mod) % mod;

					if(std::isfinite(m) && std::isfinite(n)){
						//printf("val of time = %d, val of dt = %d",n,m);
						hist(n,m)++;
					}
				}
			}
		}
	}
}

void
generate_sst_histogram_3d(const vector<string> &pass2_files, const Mat1f &pass2, 
	                   const Mat1f &reference_sst, const Mat1f &lon, const Mat1f &sza, Mat1f &hist)
{
	int i,x,y,m,n,p;
	int time_size = pass2_files.size();
	int d_range = 10;
	int mod = SAMPLING_RATE * 24;
	float step = 0.2;
	float hour_stamp_n,minute_stamp_n;
	float dt;
	string filename, minute_stamp, hour_stamp;

	hist.setTo(0);
	for(i = 0; i < time_size; ++i){


		
		filename = generate_filename(pass2_files[i]);
	
		minute_stamp = filename.substr(11,2);
	
		minute_stamp_n = atoi(minute_stamp.c_str()) / 60.0;
	
		hour_stamp = filename.substr(9,2) ;
	

		hour_stamp_n = atoi(hour_stamp.c_str())+minute_stamp_n;

		for(y = 0; y < HEIGHT; ++y){
			for(x = 0; x < WIDTH; ++x){
				if(std::isfinite(pass2(y,x,i)) && std::isfinite(reference_sst(y,x))){
					dt = pass2(y,x,i) - reference_sst(y,x);
					m = floor((dt +d_range)/step);
					n = ((((int) floor(SAMPLING_RATE*(hour_stamp_n + lon(y,x)/15)))  % mod) + mod) % mod;
					p = floor(fabs(sza(y,x)));

					if(std::isfinite(m) && std::isfinite(n) && std::isfinite(p)){
						if(m >= 0 && m < WIDTH_HIST && n >= 0 && n < HEIGHT_HIST && p >= 0 && p < DEPTH_HIST){
							//printf("val of time = %d, val of dt = %d",n,m);
							hist(n,m,p)++;
						}

					}
				}
			}
		}
	}

}

void
compute_expected_value(const Mat1f &hist, Mat1f &mean)
{
	int y,x;
	double sum,count;
    float *delta = new float[WIDTH_HIST];
	float step = 0.2;
	int d_range = 10 ;
	for(x = 0; x< WIDTH_HIST; ++x){
		delta[x] = -d_range + step*x;
	}

	for(y = 0; y < HEIGHT_HIST; ++y){
		sum = 0; count = 0;
		for(x = 0; x < WIDTH_HIST; ++x){
			sum += delta[x]*hist(y,x);
			count += hist(y,x);
		}
		if(count != 0){
			mean(y,0) = sum/count;
		}
	}
	delete [] delta;
}

/*
void
compute_enhanced_reference(const Mat1f &pass2, const Mat1f &means, const Mat1f &lons, const Mat1f &ref, Mat1f &rec, int time_size)
{
	int y,x,t;
	int offset = (time_size - TIME_SIZE)/2;
	float positions[TIME_SIZE*3];
	float *p;

	// make indexing for padding easier
	for(t = 0 ; t < TIME_SIZE; ++t){
		positions[t] = means[t];
		positions[t+TIME_SIZE] = means[t];
		positions[t+2*TIME_SIZE] = means[t];
	}
	*p = positions + TIME_SIZE - offset;

	for(y = 0; y < HEIGHT; ++y){
		for(x = 0; x < WIDTH; ++x){

			if(l2p_mask == 0){
				pos = floor((lons(y,x)/15));
				for(t = 0; t < time_size; ++t){
					if(np.isfinite(pass2(y,x,t))){
						m.append(p[pos+t]);
						dt.append(pass2(y,x,t) - ref(y,x));
						ind.append(t);
					}
				}

				len = ind.size();
				MatrixXf Q(len,2);
				MatrixXf
				VectorXf dt_m(len);
				for(t = 0; t < len; ++t){
					Q(t,0) = 1;
					Q(t,1) = m[t];
					dt_m(t) = dt[t];
				}

				z = (Q*Q.transpose()).llt.solve(Q.transpose()*dt_m);
				if(z[1] >= 0){
					for(t = 0; t < len; ++t){
						rec(y,x,ind[t]) = ref(y,x)+z[0]+z[1]*m[t];
					}
					
				}
				else{
					
				}

			}

		}
	}
}
*/



void
distribution_check(const Mat1b &hist_mask,Mat1f &approx,const Mat1b &l2p_mask,
	               int collated_interp_size,string ref_file, const vector<string> &approx1_paths)
{
	int x,y,t,m,n;
	vector<float> hour_stamps;
	string filename, minute_stamp, hour_stamp;
	float minute_stamp_n, hour_stamp_n;
	float dt;
	int d_range = 10;
	float step = 0.2;
	int mod = SAMPLING_RATE * 24;

	Mat1f lon(HEIGHT,WIDTH);
	Mat1f sst_ref(HEIGHT,WIDTH);

	get_var(ref_file,lon,"longitude");
    get_var(ref_file,sst_ref, "sst_reynolds");

	for(t = 0; t < collated_interp_size; ++t){
		filename = generate_filename(approx1_paths[t]);
		printf("filename = %s\n",filename.c_str());
		minute_stamp = filename.substr(11,2);
		printf("minute_stamp = %s\n",minute_stamp.c_str());
		minute_stamp_n = atoi(minute_stamp.c_str()) / 60.0;
		hour_stamp = filename.substr(9,2) ;
		hour_stamp_n = atoi(hour_stamp.c_str())+minute_stamp_n;
		hour_stamps.push_back(hour_stamp_n);
		printf("hour stamp at %d = %f\n",t,hour_stamps[t]);
	}

	//SAVENC(lon);
	for(y = 0; y < HEIGHT; ++y){
		for(x = 0; x < WIDTH; ++x){
			lon(y,x) = lon(y,x)/15.0;
		}
	}
	//save_test_nc_float(lon,"lons2.nc");
	
	for(y = 0; y < HEIGHT; ++y){
		for(x = 0; x < WIDTH; ++x){
			if(l2p_mask(y,x) == 0){
				//printf("pass\n");
				for(t = 0; t < collated_interp_size; ++t){
					if(std::isfinite(approx(y,x,t)) && std::isfinite(sst_ref(y,x)) ){
						//printf("in t loop\n");
						dt = approx(y,x,t) - sst_ref(y,x);
						m = floor((dt+d_range)/step);
						n =  (((((int) floor(SAMPLING_RATE*(hour_stamps[t] + lon(y,x))))  % mod) + mod) % mod);
						//printf("m = %d, n = %d\n",m,n);
						if(std::isnan(m) || std::isnan(n) || hist_mask(n,m) == 0){
							approx(y,x,t) = NAN;
						}
					}
				}
			}
		}
	}
}


void
distribution_check_3d(const Mat1b &hist_mask,Mat1f &approx,const Mat1b &l2p_mask,
	               int collated_interp_size,const Mat1f &reference_sst, const Mat1f &sza, const vector<string> &approx1_paths, const string path)
{
	int x,y,t,m,n,p;
	vector<float> hour_stamps;
	string filename, minute_stamp, hour_stamp;
	float minute_stamp_n, hour_stamp_n;
	float dt;
	int d_range = 10;
	float step = 0.2;
	int mod = SAMPLING_RATE * 24;

	Mat1f lon(HEIGHT,WIDTH);


	get_var(path,lon,"lon");


	for(t = 0; t < collated_interp_size; ++t){
		filename = generate_filename(approx1_paths[t]);
		printf("filename = %s\n",filename.c_str());
		minute_stamp = filename.substr(11,2);
		printf("minute_stamp = %s\n",minute_stamp.c_str());
		minute_stamp_n = atoi(minute_stamp.c_str()) / 60.0;
		hour_stamp = filename.substr(9,2) ;
		hour_stamp_n = atoi(hour_stamp.c_str())+minute_stamp_n;
		hour_stamps.push_back(hour_stamp_n);
		printf("hour stamp at %d = %f\n",t,hour_stamps[t]);
	}

	//SAVENC(lon);
	for(y = 0; y < HEIGHT; ++y){
		for(x = 0; x < WIDTH; ++x){
			lon(y,x) = lon(y,x)/15.0;
		}
	}
	//save_test_nc_float(lon,"lons2.nc");
	
	for(y = 0; y < HEIGHT; ++y){
		for(x = 0; x < WIDTH; ++x){
			if(l2p_mask(y,x) == 0){
				//printf("pass\n");
				for(t = 0; t < collated_interp_size; ++t){
					if(std::isfinite(approx(y,x,t)) && std::isfinite(reference_sst(y,x)) ){
						//printf("in t loop\n");
						dt = approx(y,x,t) - reference_sst(y,x);
						m = floor((dt+d_range)/step);
						n =  (((((int) floor(SAMPLING_RATE*(hour_stamps[t] + lon(y,x))))  % mod) + mod) % mod);
						p = floor(fabs(sza(y,x)));
						//printf("m = %d, n = %d\n",m,n);
						if(std::isnan(m) || std::isnan(n) || std::isnan(p) || hist_mask(n,m,p) == 0){
							approx(y,x,t) = NAN;
						}
					}
				}
			}
		}
	}
}

void
histogram_2d(const vector<string> &mask_files,const vector<string> &sst, string ref_file, int mode, Mat1f &hist)
{
	int i;
	int time_size = mask_files.size();
	int dims[3] = {HEIGHT,WIDTH,time_size};
	//int ORIGINAL_LAG = FILTER_WINDOW_LAG+SECOND_PASS_LAG;

	Mat1b clear_masks(HEIGHT,WIDTH);
	Mat1f clear_samples(3,dims);
	Mat1f lons(HEIGHT,WIDTH);
	Mat1f reference(HEIGHT,WIDTH);
	
	Mat1f sza(HEIGHT,WIDTH);
	//Mat1f means(height,1);
	get_var(ref_file,lons,"longitude");
    get_var(ref_file,reference, "sst_reynolds");
    get_var(ref_file, sza, "satellite_zenith_angle");
    printf("initialization complete\n");

    // pass 2
    if(mode == 0){
    	for(i = 0; i < time_size; ++i){
	        read_mask(mask_files[i],clear_masks,-1);
	        printf("read mask %s\n",sst[i].c_str());
	        readgranule_oneband(sst[i],clear_samples,i,"sea_surface_temperature");
	        printf("read file1 %s\n",sst[i].c_str());
	        apply_mask_slice(clear_masks,clear_samples,i,false);
	        printf("read file1 %s\n",sst[i].c_str());
	        remove_sza(sza,clear_samples,i);
	        printf("read file1 %s\n",sst[i].c_str());
    	}
    }
    //acspo
    else if(mode == 1){
    	for(i = 0; i < time_size; ++i){
	        read_acspo(mask_files[i],clear_masks, -1);
	        readgranule_oneband(mask_files[i],clear_samples,i,"sea_surface_temperature");
	        apply_mask_slice(clear_masks,clear_samples,i,false);
	        printf("read file1 %s\n",mask_files[i].c_str());
	    }
    }
	else{
		for(i = 0; i < time_size; ++i){
	        readgranule_oneband(mask_files[i],clear_samples,i,"sea_surface_temperature");
	        printf("read file1 %s\n",mask_files[i].c_str());
	    }
	}
	//SAVENC(clear_samples);
	printf("starting histogram\n");
    generate_sst_histogram(mask_files, clear_samples, reference,  lons, hist);

    //compute_expected_value(hist, means, height,width);
    //save_and_update("means.nc",means,"means", true, height,1); 
    printf("saving histogram\n");
    //SAVENC(hist);
    if(mode == 0) save_and_update("histogram_pass2.nc", hist,"histogram", true, HEIGHT_HIST, WIDTH_HIST);
    else if(mode == 1) save_and_update("histogram_acspo.nc", hist,"histogram", true, HEIGHT_HIST, WIDTH_HIST);
	else save_and_update("histogram_reinstated.nc", hist,"histogram", true, HEIGHT_HIST, WIDTH_HIST);

}

void
histogram_3d(const vector<string> &mask_files,const vector<string> &sst, Mat1f &reference_sst, Mat1f &sza, int mode, Mat1f &hist)
{
	int i;
	int time_size = mask_files.size();
	int dims[3] = {HEIGHT,WIDTH,time_size};
	//int ORIGINAL_LAG = FILTER_WINDOW_LAG+SECOND_PASS_LAG;

	Mat1b clear_mask(HEIGHT,WIDTH);
	Mat1f clear_samples(3,dims);
	Mat1f lons(HEIGHT,WIDTH);
	
	get_var(sst[0].c_str(),lons,"lon");
    printf("initialization complete\n");

    // pass 2
    if(mode == 0){
    	for(i = 0; i < time_size; ++i){
	        read_mask(mask_files[i],clear_mask,-1);
	        readgranule_oneband(sst[i],clear_samples,i,"sea_surface_temperature");
	        apply_mask_slice(clear_mask,clear_samples,i,false);
	        printf("read file1 %s\n",sst[i].c_str());
    	}
    }
    //acspo
    else if(mode == 1){
    	for(i = 0; i < time_size; ++i){
	        read_acspo(mask_files[i],clear_mask, -1);
	        readgranule_oneband(mask_files[i],clear_samples,i,"sea_surface_temperature");
	        apply_mask_slice(clear_mask,clear_samples,i,false);
	        printf("read file1 %s\n",mask_files[i].c_str());
	    }
    }
	else{
		for(i = 0; i < time_size; ++i){
	        readgranule_oneband(mask_files[i],clear_samples,i,"sea_surface_temperature");
	        printf("read file1 %s\n",mask_files[i].c_str());
	    }
	}


	//SAVENC(clear_samples);
	printf("starting histogram\n");
    generate_sst_histogram_3d(mask_files, clear_samples, reference_sst,  lons, sza, hist);
 
    printf("saving histogram\n");


}

void
distribution_check_io(const vector<string> &pass2_files, const vector<string> &sst, Mat1b &hist_mask, Mat1b &l2p_mask, string ref_file)
{	
	printf("in function\n");
	int i;
	int time_size = pass2_files.size();
	printf("time size = %d\n", time_size);
	int dims[3] = {HEIGHT,WIDTH,time_size};
	Mat1f clear_samples(3,dims);
	Mat1b clear_mask(3,dims);
	Mat1b clear_masks(HEIGHT,WIDTH);
	//int ORIGINAL_LAG = FILTER_WINDOW_LAG+SECOND_PASS_LAG;
	printf("initiated vars\n");

	for(i = 0; i < time_size; ++i){
		//printf("file = %s\n",pass2_files[i].c_str());
        read_mask(pass2_files[i],clear_masks,-1);
        readgranule_oneband(sst[i],clear_samples,i,"sea_surface_temperature");
        apply_mask_slice(clear_masks,clear_samples,i,false);
        //remove_sza(sza,clear_samples,i);
        printf("read file1 %s\n",sst[i].c_str());
    }


	distribution_check(hist_mask, clear_samples, l2p_mask, time_size,ref_file, pass2_files);

	//SAVENC(clear_samples);
	generate_mask_3d(clear_samples, clear_mask, time_size);
	SAVENC(clear_mask);
	save_mat(pass2_files, clear_mask, "brightness_temperature_11um2",true,HEIGHT, WIDTH);
}