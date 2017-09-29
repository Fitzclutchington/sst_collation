void
apply_l2p_flags(const Mat1w &l2p_mask,Mat1f &bt11_clear, int ind, bool slice)
{
    int x,y;
    if(slice){
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                if( l2p_mask(y,x) == 255 ){
                    bt11_clear(y,x) = NAN;
                }
            }
        }
    }
    

    else{
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                if( l2p_mask(y,x) == 255 ){
                    bt11_clear(y,x,ind) = NAN;
                }
            }
        }
    }    
}

/*

Cloud mask determined by the difference of neighboring pixels diagonal

*/
void 
compute_diagmask(const Mat1f &samples, Mat1w &clear_samples, int cur_ind)
{
    int x,y,i, ind;
    int w = DIAG_LAG;
    float DD,prod;
    //Mat1f DD_test(3,dims);
    Mat1f window(1,DIAG_SIZE);
    Mat1f u = Mat1f::ones(1,DIAG_SIZE);
    Mat1f powers;
    Mat1f D(1,DIAG_SIZE);
    u = u/sqrt(sum(u).val[0]);
    //printf("starting diag loop\n");
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(!(clear_samples(y,x,cur_ind)&L2P)){
    
                prod = 0;
                DD=0;
                for(i=0;i<DIAG_SIZE;i++){
                    ind = cur_ind-w+i;
                    window(i) = samples(y,x,ind);
                    prod+=(window(i)*u(i));
                }
                
                for(i=0;i<DIAG_SIZE;i++){
                    D(i) = window(i) - prod*u(i);
                    D(i) = D(i)*D(i);
                    DD+=D(i);
                }
                if(DD > T_DIAG || std::isnan(DD)){
                    clear_samples(y,x,cur_ind) = clear_samples(y,x,cur_ind) | DIAG;
                }

            }
        }
    }

}

/*

Cloud mask determined by the difference of neighboring pixels

*/
void 
compute_nnmask(const Mat1f &samples, Mat1w &clear_samples, int cur_ind)
{
    int x, y, ind_right, ind_left;
    float D_right, D_left;  
    ind_right =  cur_ind + 1; 
    ind_left = cur_ind - 1; 
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(!(clear_samples(y, x, cur_ind) & L2P)){                
                D_right = fabs(samples(y, x, cur_ind) - samples(y, x, ind_right));
                D_left = fabs(samples(y, x, cur_ind) - samples(y, x, ind_left));
                if((std::isfinite(D_right) && D_right > T_NN) || (std::isfinite(D_left) && D_left > T_NN) || (std::isnan(D_right) && std::isnan(D_left))){
                    clear_samples(y, x, cur_ind) = clear_samples(y,x,cur_ind) | NN;
                }
            }
        }
    }
}

void 
compute_bt3_ratio(const Mat1f &bt08, const Mat1f &bt11, const Mat1f &bt12, Mat1w &clear_mask, int cur_ind)
{
    int y,x;
    Mat1f bt_ratio(HEIGHT,WIDTH);
    bt_ratio = bt08 + 0.8*bt11 - 1.8*bt12;
  
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(!(clear_mask(y, x, cur_ind) & L2P) && bt_ratio(y,x) < T_BT3_RATIO){
                clear_mask(y,x,cur_ind) = clear_mask(y,x,cur_ind) | BT3_RATIO;
            }
        }
    }
}

void
compute_bt2_ratio(const Mat1f &bt08, const Mat1f &bt10, Mat1w &clear_mask, int cur_ind)
{
    int y,x;
    Mat1f bt_ratio(HEIGHT,WIDTH);
    bt_ratio = 100* (bt10 - bt08)/(bt10 +bt08);
  
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(!(clear_mask(y, x, cur_ind) & L2P) && (bt_ratio(y,x) < T_BT2_RATIO)){
                clear_mask(y,x,cur_ind) = clear_mask(y,x,cur_ind) | BT2_RATIO;
            }
        }
    }
}

void
compute_bt12_mask(const Mat1f &sst, const Mat1f &bt12, Mat1w &clear_mask, int cur_ind)
{
    int y,x;
  
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(!(clear_mask(y, x, cur_ind) & L2P) && (bt12(y,x) > sst(y,x))){
                clear_mask(y,x,cur_ind) = clear_mask(y,x,cur_ind) | BT12_TEST;
            }
        }
    } 
}

void
compute_median_mask(const Mat1f &sst, const Mat1f &median_sst, Mat1w &clear_mask, int cur_ind)
{
    int y,x;
    Mat1f median_diff = abs(median_sst - sst);

    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(!(clear_mask(y, x, cur_ind) & L2P) && (median_diff(y,x) > T_MEDIAN)){
                clear_mask(y,x,cur_ind) = clear_mask(y,x,cur_ind) | MEDIAN_TEST;
            }
        }
    } 
}

void compute_std_mask(const Mat1f &sst, const Mat1f &median_sst, Mat1w &clear_mask, int cur_ind)
{
    int y,x;
    Mat1f std_median(HEIGHT,WIDTH);
    Mat1f std_sst(HEIGHT,WIDTH);

    stdfilter(median_sst, std_median, 7);
    stdfilter(sst, std_sst, 7);
    Mat1f std_diff = abs(std_median - std_sst);

    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(!(clear_mask(y, x, cur_ind) & L2P) && (std_diff(y,x) > T_STD)){
                clear_mask(y,x,cur_ind) = clear_mask(y,x,cur_ind) | STD_TEST;
            }
        }
    } 
}

void
mask_l2p(Mat1w &bt_mask,const Mat1b &l2p_mask, int cur_ind)
{
    int y,x;
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(l2p_mask(y,x)){
                bt_mask(y,x,cur_ind) = bt_mask(y,x,cur_ind) | L2P;
            }
        }
    }
}

void
apply_mask_slice(const Mat1w &mask, Mat1f &samples, int cur_ind, bool dim)
{
    // if dim is true mask is a 3d array
    // if dim is false mask is a 2d array
    // necessary to deal with 3rd index
    int y,x;
    if(dim){
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                if(mask(y,x,cur_ind) != 0){
                    samples(y,x,cur_ind) = NAN;
                }
            }
        }
    }
    else{
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                if(mask(y,x) != 0){
                    samples(y,x,cur_ind) = NAN;
                }
            }
        }
    }
}

void
generate_mask_3d(const Mat1f &clear_samples, Mat1w &clear_mask, int time_size)
{
    int y,x,t;
    clear_mask.setTo(0);
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            for(t = 0; t < time_size;  ++t){
                if(std::isfinite(clear_samples(y, x,t))){
                    clear_mask(y, x, t) = 255;
                }
            }
        }
    }
}

void
remove_sza(const Mat1f &sza,Mat1f &clear_samples,int cur_ind)
{
    // if dim is true mask is a 3d array
    // if dim is false mask is a 2d array
    // necessary to deal with 3rd index
    int y,x;
  
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(fabs(sza(y,x)) > 60){
                clear_samples(y,x,cur_ind) = NAN;
            }
        }
    }
    
}

void
mask_histogram(const Mat1f &histogram, Mat1b &hist_mask, bool slice)
{
    int y,x,i;
    if(slice){
        for(y = 0; y < HEIGHT_HIST; ++y){
            for(x = 0; x < WIDTH_HIST; ++x){
                if(histogram(y,x) < 1)  hist_mask(y,x) = 0;
                else hist_mask(y,x) = 255;
                
            }
        }
    }
    else{
        for(y = 0; y < HEIGHT_HIST; ++y){
            for(x = 0; x < WIDTH_HIST; ++x){
                for(i = 0; i < DEPTH_HIST; ++i){
                    if(histogram(y,x,i) < 1)  hist_mask(y,x,i) = 0;
                    else hist_mask(y,x,i) = 255;
                }                
            }
        }
    }
}

