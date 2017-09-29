// filter clouds
// input:
// paths - paths to  original granules
// l2p_mask - mask for l2p flags
// border_mask - mask around border between land and ocean
// clear_paths - path to masks
//
// This function applies the nn_mask, diag_mask, and eigen_mask 
// and assigns a bit in the mask file to denote whether the mask
// detected clouds or not
void
filter_clouds(const vector<string> &paths, const Mat1b &l2p_mask, const Mat1b &border_mask, 
              vector<string> &clear_paths){

    int j;
    int sample_size = paths.size();  
    int dims[3] = {HEIGHT,WIDTH,sample_size};
    // Cur controls where the mask is read from


    Mat1f bt11;
    Mat1f bt12;
    Mat1f bt08;
    Mat1f bt10;
    Mat1w bt_mask(3,dims);
    //Mat1f sst(3,dims);
    

    //sst.release(); // release sst back to memory
    
    // allocate space for brightness temps
    bt08.create(3,dims);
    bt10.create(3,dims);
    bt11.create(3,dims);
    bt12.create(3,dims);

    // Fill data arrays with inital brightness_temperature/sst data
    for(j=0;j<sample_size;j++){
        readgranule(paths[j], bt11,bt12,bt08,bt10, j);
        read_mask(clear_paths[j], bt_mask, j); // open mask
        apply_l2p_flags( l2p_mask, bt11, j,false); // set l2p flags to NAN
        apply_l2p_flags( l2p_mask, bt12, j,false);
        apply_l2p_flags( l2p_mask, bt10, j,false);
        apply_l2p_flags( l2p_mask, bt08, j,false);
    }

    // compute eigen mask 
    // offset counter to accoutn for FILTER_WINDOW
    for(j = FILTER_WINDOW_LAG; j < sample_size - FILTER_WINDOW_LAG; ++j ){
        compute_eigenvals(bt08,bt10,bt11,bt12,border_mask, bt_mask,j);
        printf("computed eigen mask for file %d/%d\n",j,sample_size - FILTER_WINDOW_LAG);
    }
    

    save_mat(clear_paths, bt_mask, "brightness_temperature_11um2",true,HEIGHT, WIDTH);
    
}
       
