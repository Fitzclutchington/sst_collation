// Compute the gradient magnitude of image src into dst.
//
// src -- source image
// dst -- destination image of gradient magnitude (output)
// dX, dY -- derivative in the X and Y directions (output)
//
void
gradientmag(const Mat1f &src, Mat1f &dst)
{
    //printf("about to seg fault\n");
    Mat1f dX(HEIGHT,WIDTH);
    Mat1f dY(HEIGHT,WIDTH);
    int y,x;
    Mat h = (Mat_<float>(5,1) <<
            0.036420, 0.248972, 0.429217, 0.248972, 0.036420);
    Mat hp = (Mat_<float>(5,1) <<
              0.108415, 0.280353, 0, -0.280353, -0.108415);

    sepFilter2D(src, dX, -1, h, hp);
    // We negate h here to fix the sign of Dy
    sepFilter2D(src, dY, -1, hp, -h);
    for(y=0;y<HEIGHT;y++){
      for(x=0;x<WIDTH;x++){
          dst(y,x) = sqrt(dX(y,x)*dX(y,x) + dY(y,x)*dY(y,x));
      }
    }
}

void
compute_gradient(const Mat1f &bt11, Mat1f &mags, int *dims, int *inds)
{
    int a,b,i,x,y;
	  Mat1f tempmags(HEIGHT,WIDTH);
    Mat1f tempmat(HEIGHT,WIDTH);
    #pragma omp for
    for( i = 0; i < dims[2];i++){
    	for(a=0;a<HEIGHT;a++){
            for(b=0;b<WIDTH;b++){
                tempmat(a,b) = bt11(a,b,inds[i]);
            }
        }

    	gradientmag(tempmat,tempmags);
        for( y = 0; y < HEIGHT; y++){
        	for( x =0; x < WIDTH; x++){
        		mags(y,x,inds[i]) = tempmags(y,x);
        	}
        }
    }
}


void
calculate_diffs(const Mat1f &a, const Mat1f &b, Mat1f &diffs,int cur_ind, bool absval){
  //printf("calculating diffs\n");

  //printf("suxxess at time %d\n",t);
  for(int y = 0;y<HEIGHT;y++){
    for(int x = 0;x<WIDTH;x++){
      diffs(y,x) = NAN;
      if(!std::isnan(a(y,x,cur_ind)) && !std::isnan(b(y,x,cur_ind))){
        if(absval){
            diffs(y,x) = fabs(a(y,x,cur_ind) - b(y,x,cur_ind));
        }
        else{
            diffs(y,x) = a(y,x,cur_ind) - b(y,x,cur_ind);
        }
      }
    }
  }
}   


void
calculate_sums(const Mat1f &a, const Mat1f &b, Mat1f &sums,int* dims, int *inds){
  //printf("calculating diffs\n");
  int time = dims[2];
    for(int t = 0; t<time; t++){
      //printf("suxxess at time %d\n",t);
      for(int y = 0;y<HEIGHT;y++){
        for(int x = 0;x<WIDTH;x++){
                sums(y,x,inds[t]) = a(y,x,inds[t]) + b(y,x,inds[t]);            
          }
        }
    }
    
}

void 
pwl_interp_1d (const vector<int> &xd, const vector<float> &yd, const vector<int> &xi, 
                vector<float> &yi, float threshold_y, float threshold_x)


//  Purpose:
//
//    PWL_INTERP_1D evaluates the piecewise linear interpolant.
//
//  Discussion:
//
//    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
//    linear function which interpolates the data (XD(I),YD(I)) for I = 1
//    to ND.
//
//  Parameters:
//
//    Input, vector<int> xd, the data points.
//
//    Input, vector<float> yd, the data values.
//
//    Input, vector<float> yr_r, the data points of reference
//
//    Input, vector<int> xi, the interpolation points.
//
//
{
  int i;
  int k;
  int pos =0;
  float t;
  int nd = xd.size();
  int ni = xi.size();
  //printf("nd = %d ni = %d\n",nd,ni);
  //if there is only one data point, set all interpolation
  //points to the value at that point
  if ( nd == 1 ){
    for ( i = 0; i < ni; i++ ){
      yi.push_back(NAN);
    }
  }

  else{
    //for xi before xd interval
    while(xi[pos]<xd[0] && pos < ni){
        yi.push_back(NAN);
        pos++;
    }

    if(xi[pos] < xd[nd-1] && pos < ni){
      for ( k = 0; k< nd-1; k++ ){

        if(xi[pos] < xd[k+1]){
          if(fabs(yd[k+1] - yd[k]) < threshold_y && fabs(xd[k+1] - xd[k]) < threshold_x){

            while(xi[pos] < xd[k+1]){

              t = ( xi[pos] - xd[k+1] ) / (float)( xd[k] - xd[k+1] );
              yi.push_back(( 1.0 - t ) * yd[k+1] + t * yd[k]);
              pos++;
            }
          }
          else{

            while(xi[pos] < xd[k+1]){
              
              yi.push_back(NAN);
              pos++;
            }
          }
        }
      }
    }

    while(xi[pos] > xd[nd-1] && pos < ni){
        //printf("outside loop pos=%d\n",pos);
        yi.push_back(NAN);
        pos++; 
    }
  }
}


string
convert_int_to_string(int val){
  std::stringstream ss; 
  ss << val;
  string str = ss.str();
  return str;
}


void
windowed_nanmean_3d(const Mat1f &time_count,const Mat1f &time_sum,const Mat1b &l2p_mask,Mat1f &smooth_output,  const int* window,float threshold)
{

  int y,x,i,j;
  int y_dim = window[0];
  int x_dim = window[1];

  float count,sum,left_sum,left_count;
  smooth_output.setTo(NAN);

  //calculate stats in time dimension
  for(y=y_dim;y<HEIGHT-y_dim;y++){
    x=x_dim;
    sum = 0;
    count  = 0;
    left_sum = 0;
    left_count = 0;
    //find first first valid pixel
    while(l2p_mask(y,x) == 255){
    // && ice_masks(y,x,t_dim) != 0){
      x++;
    }
    // calc first window sum and count

    for(i=-y_dim;i<y_dim+1;i++){
      for(j=-x_dim;j<x_dim+1;j++){
        //if(y+i > HEIGHT || x +j > WIDTH){
        //  printf("x = %d , y= %d\n", y+i, x+j);
        //}          
        count+= time_count(y+i,x+j);
        sum += time_sum(y+i,x+j);         
        
      }
    }


    if(count > threshold && x<WIDTH-x_dim){

      smooth_output(y,x) = sum /count;
    }

    // now repeat for all other valid x values
    for(x=x+1;x<WIDTH-x_dim;x++){
      //remove all values of previous left
      sum -= left_sum;
      count -= left_count;

      //calculate new left and new right
      left_sum = 0;
      left_count = 0;
      for(i=-y_dim;i<y_dim+1;i++){  
        left_count+=time_count(y+i,x-x_dim);
        left_sum += time_sum(y+i,x-x_dim);  
        count+= time_count(y+i,x+x_dim);
        sum += time_sum(y+i,x+x_dim);     
      }


      //printf("count = %f sum = %f\n", count, sum);
      if(l2p_mask(y,x) == 0){
        // calculate this pixels value at 3dsmooth(y,x)
        if(count>threshold){
          smooth_output(y,x) = sum /count;
        }
      }
    }
  }
}

void
get_landborders(const Mat1b &land_mask, Mat1b &border_mask,int kernel_size)
{
  int x,y;
  Mat element = getStructuringElement( MORPH_RECT,
                                       Size( kernel_size, kernel_size ) );
  dilate(land_mask, border_mask,element);
  for(y=0;y<HEIGHT;y++){
    for(x=0;x<WIDTH;x++){
      border_mask(y,x) -= land_mask(y,x);
    }
  }
}


void
compute_eigenvals(const Mat1f &bt08,const Mat1f &bt10,const Mat1f &bt11,const Mat1f &bt12,
                  const Mat1b border_mask, Mat1b &clear_mask,int cur_ind)
{
  int y,x,i,j,k;
  int y_delta = 1;
  int x_delta = 1;
  int t_delta = 0;

  int min_num = (2*y_delta +1) *(2*x_delta + 1)*(2*t_delta+1)/2;
  int count_dim = 0;
  float bt08_sum,bt10_sum,bt11_sum,bt12_sum,count,window_sum,row_sum, res_mean;
  float temp_bt08;
  float temp_bt10;
  float temp_bt11;
  float temp_bt12;


  float bt08_mean;
  float bt10_mean;
  float bt11_mean;
  float bt12_mean;

  vector<float> valid_bt08;
  vector<float> valid_bt10;
  vector<float> valid_bt11;
  vector<float> valid_bt12;

  vector<int> left_inds;

  Vector4f ones(1,1,1,1);
  Vector4f e1;
  MatrixXf r;
  Matrix4f A;

  for(y=y_delta;y<HEIGHT-y_delta;y++){
    for(x=x_delta;x<WIDTH-x_delta;x++){
      if(clear_mask(y,x,cur_ind) == 0 && border_mask(y,x) == 0){
        // calc first window
        // we know that first left are nans so we don't calculate left inds     
        bt08_sum=bt10_sum=bt11_sum=bt12_sum=0;
        valid_bt08.clear();
        valid_bt10.clear();
        valid_bt11.clear();
        valid_bt12.clear();
        for(i=-y_delta;i<y_delta+1;i++){
          for(j=-x_delta;j<x_delta+1;j++){              
            for(k=-t_delta;k<t_delta+1;k++){
              //t = ((((cur_ind+k)%FILTER_TIME_SIZE)+FILTER_TIME_SIZE) % FILTER_TIME_SIZE);
              temp_bt08 = bt08(y+i,x+j,cur_ind+k);
              temp_bt10 = bt10(y+i,x+j,cur_ind+k);
              temp_bt11 = bt11(y+i,x+j,cur_ind+k);
              temp_bt12 = bt12(y+i,x+j,cur_ind+k);

              if(!std::isnan(temp_bt08) && !std::isnan(temp_bt10) && !std::isnan(temp_bt11) && !std::isnan(temp_bt12)){
                valid_bt08.push_back(temp_bt08);
                valid_bt10.push_back(temp_bt10);
                valid_bt11.push_back(temp_bt11);
                valid_bt12.push_back(temp_bt12);

                bt08_sum+= temp_bt08;
                bt10_sum+=temp_bt10;
                bt11_sum+=temp_bt11;
                bt12_sum+=temp_bt12;
              }
            }
          }
        }
  
        //if numberof pixels in window is greater tan threshold
        // calculate the mean of the norm of the pixels
        // projected into the second eigenvector
        count = valid_bt08.size();
        count_dim = valid_bt08.size();
        
        if(count > min_num){
          bt08_mean =bt08_sum/count;
          bt10_mean =bt10_sum/count;
          bt11_mean =bt11_sum/count;
          bt12_mean =bt12_sum/count;

          MatrixXf window(count_dim,4);
          for(i = 0; i < count; ++i){
            window(i,0) = valid_bt08[i] - bt08_mean;
            window(i,1) = valid_bt10[i] - bt10_mean;
            window(i,2) = valid_bt11[i] - bt11_mean;
            window(i,3) = valid_bt12[i] - bt12_mean;
          }
          
          A = (window.transpose()*window);
          e1 = A*(A*ones);
          e1/=sqrt(e1.transpose()*e1);
          r = window - (window*e1)*e1.transpose();
          window_sum =0;
          for(i = 0;i < count; ++i){
            row_sum = 0;
            row_sum+=r(i,0)*r(i,0);
            row_sum+=r(i,1)*r(i,1);
            row_sum+=r(i,2)*r(i,2);
            row_sum+=r(i,3)*r(i,3);
            row_sum = sqrt(row_sum);
            window_sum += row_sum;
          }
          
          res_mean = window_sum/(float)valid_bt08.size();
          //vals(y,x) = res_mean;
          if(res_mean > T_EIGEN){
            clear_mask(y,x,cur_ind) = clear_mask(y,x,cur_ind) | EIGEN;
          }
        }
      }
    }
  }
  //SAVENC(vals);
}


void
calculate_bt_ratio(const Mat1f &bt08, const Mat1f &bt11,const Mat1f &bt12, 
                   Mat1f &bt_ratio,int ind)
{
  int y,x;
  for(y = 0; y < HEIGHT; ++y){
    for(x = 0; x < WIDTH; ++x){
      bt_ratio(y,x) = bt08(y,x,ind) + 0.8*bt11(y,x,ind) - (1+0.8)*bt12(y,x,ind);
    }
  }
}

string
generate_filename(const string file_loc)
{
  string div = "/";
  size_t temp_pos = 0;
  size_t pos = temp_pos;
  string filename;
  while(temp_pos != string::npos){
    temp_pos = file_loc.find(div,pos+1);
    if(temp_pos != string::npos){
      pos = temp_pos;
    }
  }
  
  filename = file_loc.substr(pos);
  return filename;
}


string
generate_foldername(const string file_loc)
{
  string div = "/";
  size_t temp_pos = 0;
  size_t pos = temp_pos;
  string filename;
  while(temp_pos != string::npos){
    temp_pos = file_loc.find(div,pos+1);
    if(temp_pos != string::npos){
      pos = temp_pos;
    }
  }

  filename = file_loc.substr(0,pos);
  return filename;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// remove_high_derivatives:
//
// Input:
// Mat1f smooth - matrix that needs derivatives removed
// Mat1b l2p_mask - a HEIGHT x WIDTH matrix where 2555 is a land or invalid pixel and 0 is a valid pixel 
// int time_size - length in time of smooth matrix
//
// This function removes neighboring pixels with high derivatives determined by (fabs(first - second) > T_DERIV))
//
// Note: the final slice of the matrix is ignored since there is no t+1 pixel to 
// examine to determine if there is a high derivative
//
void
remove_high_derivatives(Mat1f &smooth, const Mat1b &l2p_mask, int time_size)
{
  int x,y,t;
  float first, second;
  for(y=0;y<HEIGHT;y++){
    for(x=0;x<WIDTH;x++){
      if(l2p_mask(y,x) == 0){
        // -1 so that we don't over extend boundary
        for(t=0;t<time_size-1;t++){
          first = smooth(y,x,t);
          second = smooth(y,x,t+1);
          if(std::isnan(first) || (std::isfinite(second) && (fabs(first - second) > T_DERIV)) || std::isnan(second)){ 
            smooth(y,x,t) = NAN;
          }
        }
      }
    }
  }
}

void
remove_last_value(Mat1f &collated_interp, int time_size)
{
  int y,x,t;
  for(y = 0; y < HEIGHT; ++y){
    for(x = 0; x < WIDTH; ++x){
      for(t = 0; t < time_size-1; ++t){
        if(std::isfinite(collated_interp(y,x,t)) && std::isnan(collated_interp(y,x,t+1))){
          collated_interp(y,x,t) = NAN;
        }
      }
    }
  }
}

void
combine_l2pmasks(const Mat1b &land_mask,const Mat1b &invalid_mask, const Mat1b &ice_mask, Mat1b &l2p_mask)
{
  int y,x;
  for(y=0;y<HEIGHT;++y){
    for(x=0;x<WIDTH;++x){
      l2p_mask(y,x) = 0;
      if(land_mask(y,x) == 255 || invalid_mask(y,x) == 255 || ice_mask(y,x) == 255){
        l2p_mask(y,x) = 255;
      }
    }
  }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// update_sums:
//
// Input:
// Mat1f time_count - matrix of the number of finite values in the window for each pixel
// Mat1f time_sum - matrix of the sum of finite values in the window for each pixel
// Mat1f masked_data - matrix that contains all the masked original data
// int cur_ind - index of the granule that is being processed
// sting mode - determines finctionallity of function:
//              "right" - only adds right granule
//              "left"  - only removes left granule
//              "both"  - removes left and adds right
//
// update time_count and time_sum by adding the right granule value and removing the left granule
//
void
update_sums(Mat1f &time_count,Mat1f &time_sum, const Mat1f &masked_data, int cur_ind, int lag, string mode)
{
  
  int y,x;
  float val;
  int left = cur_ind - lag -1;
  int right = cur_ind + lag;

  if(mode == "right"){
    for(y = 0; y < HEIGHT; ++y){
      for(x = 0; x < WIDTH; ++x){
        val = masked_data(y,x,right);
        if(std::isfinite(val)){
          time_count(y,x)++;
          time_sum(y,x) += val;
        }
      }
    }
  }

  else if(mode == "left"){
    for(y = 0; y < HEIGHT; ++y){
      for(x = 0; x < WIDTH; ++x){
        val = masked_data(y,x,left);
        if(std::isfinite(val)){
          time_count(y,x)--;
          time_sum(y,x) -= val;
        }
      }
    }
  }

  else if(mode == "both"){
    for(y = 0; y < HEIGHT; ++y){
      for(x = 0; x < WIDTH; ++x){
        val = masked_data(y,x,left);
        if(std::isfinite(val)){
          time_count(y,x)--;
          time_sum(y,x) -= val;
        }

        val = masked_data(y,x,right);
        if(std::isfinite(val)){
          time_count(y,x)++;
          time_sum(y,x) += val;
        }
      }
    }
  }

}

void
remove_long_interpolation(const Mat1f &reinstated_clear, Mat1f &collated, vector<int> collated_inds, int sample_size)
{
  int y,x,t, i, ind;
  int time_size = collated_inds.size();
  int clear_count;
  for(y = 0; y < HEIGHT; ++y){
    for(x = 0; x < WIDTH; ++x){
      for(i = 0; i < time_size; ++i){
        ind = collated_inds[i];
        clear_count = 0;
        if(ind < COLLATED_LAG){
          for(t = -ind; t < COLLATED_LAG+1; ++t){
            if(std::isfinite(reinstated_clear(y,x,ind+t))){
              clear_count++;
            }                  
          }
        }

        else if(ind >= COLLATED_LAG && ((sample_size - ind) > COLLATED_LAG)){
          for(t = -COLLATED_LAG; t < COLLATED_LAG+1; ++t){
            if(std::isfinite(reinstated_clear(y,x,ind+t))){
              clear_count++;
            }            
          }
        }

        else if(sample_size - ind <= COLLATED_LAG){
          for(t = -COLLATED_LAG; t < (sample_size - ind); ++t){
            if(std::isfinite(reinstated_clear(y,x,ind+t))){
              clear_count++;
            }      
          }
        }

        else eprintf("remove_long_interpolation: index %d invalid\n", ind);

        if( clear_count == 0 ) collated(y,x,i) = NAN;
      }
    }
  }

}


void
compute_angle(const Mat1f &above_collated,const Mat1f &collated_interp,Mat1f &angle,int collated_interp_size)
{
  int y,x,t;
  float a,c; 
  double acc_a = 0, acc_c = 0, acc=0;
  for(y = 0; y < HEIGHT; ++y){
    for(x = 0; x < WIDTH; ++x){
      acc_a = acc_c = acc = 0;
      
      for(t = 0; t < collated_interp_size; ++t){
        a = above_collated(y,x,t);
        c = collated_interp(y,x,t);
        if(std::isfinite(a) && std::isfinite(c)){
          acc += a*c;
          acc_a += a*a;
          acc_c += c*c;
        }
      }
      
      if(acc_a !=0 && acc_c != 0){
        angle(y,x) = acc / (sqrt(acc_a) * sqrt(acc_c));
      }

    }
  }
}

void
compute_var(const Mat1f &above_collated,const Mat1f &collated_interp,Mat1f &var,int collated_interp_size)
{
  int y,x,t;
  float a, c, count; 
  double acc=0;
  for(y = 0; y < HEIGHT; ++y){
    for(x = 0; x < WIDTH; ++x){
      acc = count = 0;
      
      for(t = 0; t < collated_interp_size; ++t){
        a = above_collated(y,x,t);
        c = collated_interp(y,x,t);
        if(std::isfinite(a) && std::isfinite(c)){
          acc += (a - c)*(a - c);
          count++;
        }
      }
      
      if(count != 0){
        var(y,x) = acc / (count - 1);
      }

    }
  }
}


// if a is above b 
// above = a
void
compute_above_samples(const Mat1f &a, const Mat1f &b, Mat1f &above, const Mat1b &l2p_mask, const int time_size)
{
    int y,x,i;
    above.setTo(NAN);
    for(i = 0; i < time_size; ++i){
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                if(l2p_mask(y,x) == 0){
                  if(a(y,x,i) > b(y,x,i) - DELTA_SST && std::isfinite(a(y,x,i)) && std::isfinite(b(y,x,i))){
                      above(y,x,i) = a(y,x,i);
                  }
                  else if( std::isfinite(b(y,x,i))){
                      above(y,x,i) = b(y,x,i);
                  }
                }
            }
        }
    }
}

void
compute_diff_mean(const Mat1f &reinstated,const Mat1f &collated_interp,Mat1f &mean,int collated_interp_size)
{
  int y,x,t;
  float a, c, count; 
  double acc=0;
  for(y = 0; y < HEIGHT; ++y){
    for(x = 0; x < WIDTH; ++x){
      acc = count = 0;
      
      for(t = 0; t < collated_interp_size; ++t){
        a = reinstated(y,x,t);
        c = collated_interp(y,x,t);
        if(std::isfinite(a) && std::isfinite(c)){
          acc += (a - c)*(a - c);
          count++;
        }
      }
      
      if(count != 0){
        mean(y,x) = acc / (count );
      }

    }
  }
}

void
set_gamma(MatrixXf &Gamma, const float g)
{
    int i;
    int rows = Gamma.rows();
    int cols = Gamma.cols();

    if( rows != cols ) eprintf("Gamma is a %dx%d, must be an nxn matrix\n", rows, cols);

    Gamma.setZero();
    Gamma(0,0) = g;Gamma(0,1) = -g;
    Gamma(rows-1,rows-1) = g; Gamma(rows-1,rows-2) = -g;
    for(i=1;i<rows-1;i++){
        Gamma(i,i)   = 2*g;
        Gamma(i,i-1) = -g;
        Gamma(i,i+1) = -g;
    }
}

void
set_gamma_var(MatrixXf &Gamma, const vector<float> sw)
{
    int i;
    int rows = Gamma.rows();
    int cols = Gamma.cols();
    int size = sw.size();

    if( rows != cols ) eprintf("Gamma is a %dx%d, must be an nxn matrix\n", rows, cols);
    if( size != rows ) eprintf("sw is size %d rows is size %d\n",size,rows);
    Gamma.setZero();
    Gamma(0,0) = sw[0];Gamma(0,1) = -sw[0];
    Gamma(rows-1,rows-1) = sw[size-2]; Gamma(rows-1,rows-2) = -sw[size-2];
    for(i=1;i<rows-1;i++){
        Gamma(i,i)   = sw[i-1] + sw[i];
        Gamma(i,i-1) = -sw[i-1];
        Gamma(i,i+1) = -sw[i];
    }
}

void
compute_gaussian_weights(vector<float> &weights)
{
  int i;
  int size = 2*COLLATED_LAG + 1;
  float w;
  float denominator = 2*pow(SIGMA_H, 2);
  for(i = 0; i < size; ++i){
    w = exp(-pow((-1 + i/(float) COLLATED_LAG),2)/denominator);
    weights.push_back(w);
  }  
}

void
remove_gaps(Mat1f &smooth_samples, Mat1f &clear_samples, const Mat1b &l2p_mask,int sample_size, int collated_size)
{
  int y,x,t, i, segment_num, first_valid, start, end, count, valid_size, diff;
  vector<int> valid_inds;
  vector<int> segment_start;
  vector<int> segment_end;

  for(y = 0; y < HEIGHT; ++y){
    for(x = 0; x < WIDTH; ++x){
      if(l2p_mask(y,x) == 0){
        valid_inds.clear();
        segment_start.clear();
        segment_end.clear();
        count = 0;
        // first get array of valid pixels
        for(t = 0; t < collated_size; ++t){
          if(std::isfinite(smooth_samples(y,x,t))){
            valid_inds.push_back(t);
          }
        }

        //take sequential diffs
        valid_size = valid_inds.size();
        if(valid_size > 1){
          first_valid = valid_inds[0];
          for(t = 0; t < valid_size -1; ++t){
            diff = valid_inds[t+1] - valid_inds[t];
            if(diff > 1){
              segment_start.push_back(first_valid);
              segment_end.push_back(valid_inds[t]);
              first_valid = valid_inds[t+1];
            }
          }
          if(valid_inds[t] == collated_size-1){
              segment_start.push_back(first_valid);
              segment_end.push_back(valid_inds[t]);
          }


          // count num clear in gap
          segment_num = segment_start.size();

          for(t = 0; t < segment_num; ++t){
            start = segment_start[t];
            end = segment_end[t];
            diff = end - start;
            for(i = start ; i < end+1; ++i){
              if(std::isfinite(clear_samples(y,x,i))){
                count++;
              }
            }

            if(count == 0){
              for(i = start ; i < end+1; ++i){
                smooth_samples(y,x,i) = NAN;
              }

              //handle end of time series, since clear and smooth are 
              // not the same size
              if(i == collated_size - 1){
                for(;i < sample_size; ++i){
                  smooth_samples(y,x,i) = NAN;
                }
              }
            }
            count = 0;
          }
        }
      }
    }
  }
}

void
compute_dt_thresholds(Mat1f &t_cold, Mat1f &t_warm, string ref_file)
{
  Mat1f sza(HEIGHT,WIDTH);
  int y,x;
  get_var(ref_file, sza, "satellite_zenith_angle");
  for(y = 0; y < HEIGHT; ++y){
        for(x = 0;x < WIDTH; ++x){
            if(fabs(sza(y,x)) > MAX_SZA){
                t_cold(y,x) = T_COLD_DT*0.3;
                t_warm(y,x) = T_WARM_DT*0.3;
            }             
            else{
                t_cold(y,x) = T_COLD_DT*pow(cos((90.0/MAX_SZA)*PI*sza(y,x)/180.0),EXP);
                t_warm(y,x) = T_WARM_DT*pow(cos((90.0/MAX_SZA)*PI*sza(y,x)/180.0),EXP);
            } 
        }
    }
    sza.release();
}

void
compute_reinstated(Mat1f &reinstated_clear, Mat1f &original_sst, const Mat1f &smooth, Mat1b mask,
                   Mat1f &t_cold, Mat1f &t_warm, int collated_interp_size, string ref_file)
{
  int i, y, x;
  float DD, D_ref, D_right, D_left;
  bool nn,diag;
  Mat1f reference(HEIGHT,WIDTH);
  get_var(ref_file,reference,"sst_reynolds");

  for(i = 1; i < collated_interp_size-1; ++i){        
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                if(i > 0 && i < collated_interp_size -1){
                    DD = original_sst(y,x,i) - smooth(y,x,i);
                    D_ref = original_sst(y,x,i) - reference(y,x); 
                    nn = !(mask(y,x,i) & NN);
                    diag = !(mask(y,x,i) & DIAG);
                    if((std::isfinite(DD) && DD > -DELTA_SST) && D_ref > t_cold(y,x) && D_ref < t_warm(y,x) && original_sst(y,x,i) > MIN_TEMP && nn && diag){
                        reinstated_clear(y,x,i) = original_sst(y,x,i);
                    }
                }
            }
        }
    }

    t_cold.release();
    t_warm.release();
    original_sst.release();
    reference.release();
}