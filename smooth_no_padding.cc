void
smooth_samples_collated(const Mat1f &collated_interp, Mat1f &collated_smooth,const Mat1f &l2p_mask,
                        const Mat1f reference_sst,const int* window,const int sample_size)
{
    int i,y,x,t,cur_ind;
    int threshold;
    float val;
    
    int time_lag = window[2];
    int dims[3] = {HEIGHT,WIDTH,sample_size}; // dimensions of array containing masked data


    
    Mat1f smooth_output(HEIGHT, WIDTH); // matrix of smooth vals to save to file

    Mat1f masked_data(3,dims);

    Mat1f time_sum(HEIGHT,WIDTH);
    time_sum.setTo(0);
    
    Mat1f time_count(HEIGHT,WIDTH);
    time_count.setTo(0);


    // always compute on sst
    // open sst and apply mask
    for(i = 0; i < sample_size; ++i){

        //subtract reference from granule
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                masked_data(y,x,i) = collated_interp(y,x,i) - reference_sst(y,x);
            }
        }
    }

    printf("subtracted reference\n");

    //determine first time sum and time count
    for(t = 0; t < time_lag + 1; ++t){
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                val = masked_data(y,x,t);
                if(std::isfinite(val)){
                    time_count(y,x)++;
                    time_sum(y,x) += val;
                }
            }
        }
    }
    printf("generated time sum\n");


    threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * time_lag;

    // compute first smoothb
    windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold); // rewrite this function to only work with time sum
    // add back reference
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            collated_smooth(y,x,0) = smooth_output(y,x) + reference_sst(y,x);
        }
    }

    printf("generated first colalted smooth\n");
    //compute smooth before there arent enough granules to fill entire window
    for(cur_ind = 1; cur_ind < time_lag+1; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "right");
        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (time_lag+cur_ind);
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                collated_smooth(y,x,cur_ind) = smooth_output(y,x) + reference_sst(y,x);
            }
        } 
    }   

    printf("generated first colalted smooth interval\n");

    //compute in interval where entire window will be full
    for(; cur_ind < sample_size - time_lag; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "both");

        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (2*time_lag +1);
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                collated_smooth(y,x,cur_ind) = smooth_output(y,x) + reference_sst(y,x);
            }
        }
    }
    printf("generated middle colalted smooth\n");
    //compute end interval
    for(; cur_ind < sample_size; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "left");

        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (time_lag + sample_size - cur_ind); // i counts how many granules to ignore as window shrinks
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                collated_smooth(y,x,cur_ind) = smooth_output(y,x) + reference_sst(y,x);
            }
        }         
    }
    printf("generated final colalted smooth\n");

    masked_data.release();
}

void
moving_average_time(const Mat1f &collated_interp, Mat1f &collated_smooth,
                    const string ref_file,const int time_lag ,const int sample_size)
{
    int i,y,x,t,cur_ind;
    int threshold;
    float val;
    
    int dims[3] = {HEIGHT,WIDTH,sample_size}; // dimensions of array containing masked data

    
    Mat1f reference(HEIGHT,WIDTH);

    Mat1f masked_data(3,dims);

    Mat1f time_sum(HEIGHT,WIDTH);
    time_sum.setTo(0);
    
    Mat1f time_count(HEIGHT,WIDTH);
    time_count.setTo(0);

    // get the reference data
    get_var(ref_file,reference,"sst_reynolds");

    // always compute on sst
    // open sst and apply mask
    for(i = 0; i < sample_size; ++i){

        //subtract reference from granule
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                masked_data(y,x,i) = collated_interp(y,x,i) - reference(y,x);
            }
        }
    }

    printf("subtracted reference\n");

    //determine first time sum and time count
    for(t = 0; t < time_lag + 1; ++t){
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                val = masked_data(y,x,t);
                if(std::isfinite(val)){
                    time_count(y,x)++;
                    time_sum(y,x) += val;
                }
            }
        }
    }
    printf("generated time sum\n");


    threshold = 0;


    // add back reference
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(time_count(y,x) > threshold){
                collated_smooth(y,x,0) = (time_sum(y,x)/time_count(y,x)) + reference(y,x);
            }
        }
    }

    printf("generated first colalted smooth\n");
    //compute smooth before there arent enough granules to fill entire window
    for(cur_ind = 1; cur_ind < time_lag+1; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "right");
        threshold = 0; 
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                if(time_count(y,x) > threshold){
                    collated_smooth(y,x,cur_ind) = (time_sum(y,x)/time_count(y,x)) + reference(y,x);
                }
            }
        } 
    }   

    printf("generated first colalted smooth interval\n");

    //compute in interval where entire window will be full
    for(; cur_ind < sample_size - time_lag; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "both");

        threshold = 0;
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                if(time_count(y,x) > threshold){
                    collated_smooth(y,x,cur_ind) = (time_sum(y,x)/time_count(y,x)) + reference(y,x);
                }
            }
        }
    }
    printf("generated middle colalted smooth\n");
    //compute end interval
    for(; cur_ind < sample_size; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "left");

        threshold = 0; // i counts how many granules to ignore as window shrinks
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                if(time_count(y,x) > threshold){
                    collated_smooth(y,x,cur_ind) = (time_sum(y,x)/time_count(y,x)) + reference(y,x);
                }
            }
        }         
    }
    printf("generated final colalted smooth\n");

    masked_data.release();
}