void
calc_approximate(const Mat1f &bt11_smooth, const Mat1f &bt11_clear, const Mat1b &l2p_mask, Mat1f &bt11_approx, const Mat1f reference_sst, int time_size)
{
    int x,y,z;
    double d1,d2,coeff;
    //placeholder
    int count=0;
    float val =0;
    vector<float> data_cl2;
    vector<float> smooth;
    vector<int> clear_inds;


    //Mat1f diffs(HEIGHT,WIDTH);
    bt11_approx.setTo(NAN);
    

    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(l2p_mask(y,x) == 0 ){
                count = 0;
                d1=0;d2 =0;
                data_cl2.clear();
                smooth.clear();
                clear_inds.clear();
                for(z = 0; z < time_size; ++z){
                    
                    if(std::isfinite(bt11_smooth(y,x,z)) && std::isfinite(bt11_clear(y,x,z))){

                        d1 += bt11_clear(y,x,z)*bt11_smooth(y,x,z);
                        d2 += bt11_smooth(y,x,z)*bt11_smooth(y,x,z);
                        count++;
                                
                    }                   
                }
            
                //diffs(y,x) = d2-d1;
                if(d2 != 0 && count > MIN_POINTS){
                    coeff = d1/d2;
                    for(z = 0; z < time_size; ++z){
                        val = coeff * bt11_smooth(y,x,z);
                        if((val-reference_sst(y,x)) > T_COLD_DT && (val - reference_sst(y,x)) < T_WARM_DT){
                            bt11_approx(y,x,z) = val;
                        }

                    }
                }
            }
        }
    }
}


void
set_rhs_v_window(const Mat1f &sst, VectorXf &rhs, MatrixXf &v, const vector<int> &collated_inds,
            const int x, const int y, int &first_valid, int &last_valid, const int sample_size, const vector<float> weights)
{
    int i, j, k, t, ind;
    int collated_size = collated_inds.size();
    float sst_sum, sst_count;
    int w = 1;

    first_valid = -1; last_valid = -1;
    //for each hour in collated
    //auto start = std::chrono::system_clock::now();
    for(i = 0; i < collated_size; ++i){
        //for each value in current collated hour
        sst_sum = sst_count = 0;
        ind = collated_inds[i];
        
        // handle case where there are not enough granules in the front of the interval 
        if(ind < COLLATED_LAG){
            for(t = -ind;t < COLLATED_LAG+1; ++t){
                for(k = -w; k < w+1; ++k){
                    for(j = -w; j < w+1; ++j){                                
                        if(std::isfinite(sst(y+k,x+j,ind+t))){
                            sst_sum += weights[t+COLLATED_LAG]*sst(y+k,x+j,ind+t);
                            sst_count += weights[t+COLLATED_LAG];
                        }
                    }
                }                            
            }
        }

        else if(ind >= COLLATED_LAG && ((sample_size - ind) > COLLATED_LAG)){
            for(t = -COLLATED_LAG; t < COLLATED_LAG+1; ++t){
                for(k = -w; k < w+1; ++k){
                    for(j = -w; j < w+1; ++j){                                    
                        if(std::isfinite(sst(y+k,x+j,ind+t))){
                            sst_sum += weights[t+COLLATED_LAG]*sst(y+k,x+j,ind+t);
                            sst_count += weights[t+COLLATED_LAG];
                        }
                    }
                }                            
            }
        }

        else if(sample_size - ind <= COLLATED_LAG){
            for(t = -COLLATED_LAG; t < (sample_size - ind); ++t){
                for(k = -w; k < w+1; ++k){
                    for(j = -w; j < w+1; ++j){                                    
                        if(std::isfinite(sst(y+k,x+j,ind+t))){
                            sst_sum += weights[t+COLLATED_LAG]*sst(y+k,x+j,ind+t);
                            sst_count += weights[t+COLLATED_LAG];
                        }
                    }
                }                            
            }
        }

        else eprintf("collate_smooth: index %d invalid\n", ind);

        
        if(sst_count < 1){
            sst_count = 0;
            sst_sum = 0;
        }
        
        rhs[i] = sst_sum ;
        v(i,i) = sst_count;
        if(v(i,i)!=0){                      
            if(first_valid == -1){
                first_valid = i;
            }
            if(last_valid < i){
                last_valid = i;
            }
        }
    }
}

void
set_rhs_v(const Mat1f &approx,const Mat1f &clear_samples, VectorXf &rhs, MatrixXf &v, const vector<int> &collated_inds, const float g,
            const int x, const int y, int &first_valid, int &last_valid, const int sample_size, const vector<float> weights,vector<float> &sw)
{
    int i, t, ind;
    int half_ind = COLLATED_LAG;
    int collated_size = collated_inds.size();
    float approx_sum, approx_count, clear_sum, clear_count, clear_w, approx_w,clear_p;
    

    //for each hour in collated
    //auto start = std::chrono::system_clock::now();
    for(i = 0; i < collated_size; ++i){
        //for each value in current collated hour
        clear_sum = clear_count = approx_sum = approx_count = approx_w = clear_w = clear_p = 0;
        ind = collated_inds[i]; 
        
        if(ind < COLLATED_LAG){
            for(t = -ind; t < COLLATED_LAG+1; ++t){
                if(std::isfinite(clear_samples(y,x,ind+t))){
                    clear_sum += weights[t+COLLATED_LAG]*clear_samples(y,x,ind+t);
                    clear_w += weights[t+COLLATED_LAG];
                    clear_count++;
                }
                if(std::isfinite(approx(y,x,ind+t))){
                    approx_sum += weights[t+COLLATED_LAG]*approx(y,x,ind+t);
                    approx_w += weights[t+COLLATED_LAG];
                    approx_count++;
                }                       
            }
            half_ind = (ind + COLLATED_LAG) /2;
            clear_p = clear_count / (COLLATED_LAG+1);
        }

        else if(ind >= COLLATED_LAG && ((sample_size - ind) > COLLATED_LAG)){
            for(t = -COLLATED_LAG; t < COLLATED_LAG+1; ++t){
                if(std::isfinite(clear_samples(y,x,ind+t))){
                    clear_sum += weights[t+COLLATED_LAG]*clear_samples(y,x,ind+t);
                    clear_w += weights[t+COLLATED_LAG];
                    clear_count++;
                }
                if(std::isfinite(approx(y,x,ind+t))){
                    approx_sum += weights[t+COLLATED_LAG]*approx(y,x,ind+t);
                    approx_w += weights[t+COLLATED_LAG];
                    approx_count++;
                }                    
            }
            half_ind = COLLATED_LAG;
            clear_p = clear_count / (2 * COLLATED_LAG +1);
        }

        else if(sample_size - ind <= COLLATED_LAG){
            for(t = -COLLATED_LAG; t < (sample_size - ind); ++t){
                
                if(std::isfinite(clear_samples(y,x,ind+t))){
                    clear_sum += weights[t+COLLATED_LAG] *clear_samples(y,x,ind+t);
                    clear_w += weights[t+COLLATED_LAG];
                    clear_count++;
                }
                if(std::isfinite(approx(y,x,ind+t))){
                    approx_sum += weights[t+COLLATED_LAG]*approx(y,x,ind+t);
                    approx_w+= weights[t+COLLATED_LAG];
                    approx_count++;
                }                    
                half_ind = (COLLATED_LAG + sample_size - ind)/2;
                clear_p = clear_count / (COLLATED_LAG+1);
            }
        }

        else eprintf("collate_smooth: index %d invalid\n", ind);

        if(clear_p < 0) clear_p = 0;
        if(clear_p > 1) clear_p = 1;
        sw[i] = (1-clear_p)*g;
        

        //if(clear_count < 3 && approx_count ==0){
        if(approx_count ==0){
            clear_w =0;
            clear_sum = 0;
            sw[i] = g;
        }
        else if(approx_count < half_ind){
            approx_w = 0;
            approx_sum = 0;
            sw[i] = g;
        }

        if(clear_p > .8){
            approx_w = 0;
            approx_sum = 0;
        }
        
        //sw[i]= 0.5;
        rhs[i] = MU_APPROX*approx_sum + MU_CLEAR*clear_sum;
        v(i,i) = MU_APPROX*approx_w + MU_CLEAR*clear_w;
        
        /*
        if((y==4500 && x==1690) || (y==793 && x==2463)){
            printf("point = (%d,%d):\nat %d sw = %f rhs = %f v = %f clear_p = %f\nclear_count = %f approx_count = %f\n",y,x,i,sw[i],rhs[i],v(i,i), clear_p, clear_count, approx_count);
        }
        */
        if(v(i,i)!=0){                      
            if(first_valid == -1){
                first_valid = i;
            }
            if(last_valid < i){
                last_valid = i;
            }
        }
    }
}