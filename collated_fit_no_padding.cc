#include "collation_functions.cc"

void
collate_samples(const Mat1f &clear_samples, const Mat1f &approx, Mat1f &collated, 
    const vector<int> collated_inds,const Mat1b &l2p_mask, int sample_size, const string ref_file, bool window, float g)
{
    
    // 1) set up right hand side-> rhs[k] = MU*clear_sum + MU* approx sum -> k = time step
    // 2) set up left hand side -> diagonal = mu*clear_count+mu*approx_count
    // 3) add gamma to left hand side for tridiagonal matrix
    // 4) these should all be premade eigen arrays determined by the number of collated values being produced
    // 5) collated = inv(lhs)*rhs

    // read reference for sst_reynolds and o2ld
    // if abs(collated - reynolds ) > 7 && o2ld <= 50 -> sset to nan
    int y,x,i;
    int collated_size = collated_inds.size();
    int first_valid, last_valid;
    
    vector<float> sw(collated_size);

    vector<float> weights;

    VectorXf collated_vals(collated_size);
    VectorXf rhs(collated_size);
    MatrixXf v(collated_size,collated_size);
    MatrixXf gamma(collated_size,collated_size); 

    Mat1f reference(HEIGHT,WIDTH);
    Mat1f o2ld(HEIGHT,WIDTH);
    Mat1f g_mat,v_mat, rhs_mat;

    rhs.setZero(); v.setZero();

    get_var(ref_file,reference,"sst_reynolds");
    get_var(ref_file,o2ld,"ocean_to_land_dist");

   
    compute_gaussian_weights(weights);
    /*
    for(i = 0; i < 2 * COLLATED_LAG + 1; ++i){
        printf("weights at %d = %f\n", i, weights[i]);
    }
    */
    collated.setTo(NAN);
    //for each pixel
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(l2p_mask(y,x) ==0){
                first_valid = -1; last_valid = -1;
                if(window) set_rhs_v_window(clear_samples, rhs, v, collated_inds, x,  y, first_valid, last_valid, sample_size, weights);
                else set_rhs_v(approx, clear_samples, rhs, v, collated_inds, g, x,  y, first_valid, last_valid, sample_size, weights, sw);
        

                if(last_valid != -1 && first_valid != -1){

                    /*
                    for(i = 0; i < first_valid; ++i){
                        v(i) = v(first_valid);
                        rhs[i] = rhs[first_valid];
                    }
                    for(i = last_valid; i < collated_size; ++i){
                        v(i) = v(last_valid);
                        rhs[i] = rhs[last_valid];
                    }
                    */
                    /*
                    // Handle the case where the first and last valid pixels are not 
                    // at the ends of the v matrix
                    if(first_valid != 0){
                        gamma(first_valid,first_valid) = g;    //Takes value of gamma from Top right
                        gamma(first_valid, first_valid-1) = 0; 
                        gamma(first_valid-1, first_valid-1) = g;
                        gamma(first_valid-1, first_valid) = 0;
                    }
                    if(last_valid != collated_size-1){
                        gamma(last_valid,last_valid) = g;
                        gamma(last_valid,last_valid+1) = 0;
                        gamma(last_valid+1,last_valid+1) = g;
                        gamma(last_valid+1,last_valid) = 0;
                    }
                    */
                    set_gamma_var(gamma, sw);

                    collated_vals = (v + gamma).llt().solve(rhs);


                    if(!window){
                        for(i = first_valid;i <= last_valid; ++i){
                            collated(y,x,i) = collated_vals(i);
                            if(abs(reference(y,x) - collated(y,x,i)) > 7 && o2ld(y,x) <=50){
                                collated(y,x,i) = NAN;
                            }
                        }
                    }
                    else{
                        for(i = first_valid;i <= last_valid; ++i){
                            collated(y,x,i) = collated_vals(i);
                        }
                    }
                    /*
                    // Restore gamma matrix for next pixel
                    if(first_valid != 0){
                        gamma(first_valid,first_valid) = 2*g;    //Takes value of gamma from Top right
                        gamma(first_valid, first_valid-1) = -g; 
                        gamma(first_valid-1, first_valid-1) = 2*g;
                        gamma(first_valid-1, first_valid) = -g;
                    }
                    if(last_valid != collated_size-1){
                        gamma(last_valid,last_valid) = 2*g;
                        gamma(last_valid,last_valid+1) = -g;
                        gamma(last_valid+1,last_valid+1) = 2*g;
                        gamma(last_valid+1,last_valid) = -g;
                    }
                    */
                }
            }
        }
    }
    reference.release();
    o2ld.release();
}

