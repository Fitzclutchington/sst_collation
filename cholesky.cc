void
compute_gamma_vector(vector<float> &g, const float gamma, int sample_size)
{
	int i;
	float g_2;
	g_2 = 1 + 2* gamma;

	// first value
	g.push_back(1 + gamma);

	//middle values
	for(i=1; i < sample_size-1; ++i) g.push_back(g_2);

	// last value
	g.push_back(1 + gamma);
}

void
compute_diagonal(const vector<float> g, vector<float> &d, vector<float> &u, float gamma, int sample_size)
{
	int i;
	float ui;

	d.clear(); u.clear();
	d.push_back(sqrt(g[0]));
	u.push_back(-gamma/d[0]);

	for(i = 1; i < sample_size-1; ++i){
		ui = u[i-1];
		d.push_back(sqrt(g[i]- ui*ui));
		u.push_back(-gamma/d[i]);
	}
	ui = u[i-1];
	d.push_back(sqrt(g[i]-ui*ui));
}

void
least_square_decomposition(const Mat1f &sst, Mat1f &ls_approx, const Mat1b &l2p_mask, float gamma, int sample_size, bool nan)
{
	int y,x,i;
	vector<float> d, u, g, g_a, s_a, y_m;

	int final_ind = sample_size - 1;


	ls_approx.setTo(NAN);
	compute_gamma_vector(g,gamma, sample_size);
	compute_diagonal(g, d, u, gamma, sample_size);
	
	if(!nan){
		for(y = 0; y < HEIGHT; ++y){
			for(x = 0; x < WIDTH; ++x){

				if(l2p_mask(y,x) == 0){	
					y_m.clear();
					y_m.push_back(sst(y,x,0)/d[0]);
					for(i = 1; i < sample_size; ++i){
						y_m.push_back((sst(y,x,i)-u[i-1]*y_m[i-1])/d[i]);
					}

					ls_approx(y,x,final_ind) = y_m[final_ind]/d[final_ind];
					for(i = final_ind -1; i >= 0; --i){
						ls_approx(y,x,i) = (y_m[i] - u[i]*ls_approx(y,x,i+1))/d[i];
					}
				}
			}
		}
	}
	else{
		for(y = 0; y < HEIGHT; ++y){
			for(x = 0; x < WIDTH; ++x){

				if(l2p_mask(y,x) == 0){	
					
					g_a.clear();
					s_a.clear();
					y_m.clear();
					for(i =0; i < sample_size; ++i){
						g_a.push_back(g[i]);
						s_a.push_back(sst(y,x,i));
						if(std::isnan(sst(y,x,i))){
							g_a[i]--;
							s_a[i] = 0;
						}
					}

					compute_diagonal(g_a, d, u, gamma, sample_size);

					y_m.push_back(s_a[0]/d[0]);
					for(i = 1; i < sample_size; ++i){
						y_m.push_back((s_a[i]-u[i-1]*y_m[i-1])/d[i]);
					}

					ls_approx(y,x,final_ind) = y_m[final_ind]/d[final_ind];
					for(i = final_ind -1; i >= 0; --i){
						ls_approx(y,x,i) = (y_m[i] - u[i]*ls_approx(y,x,i+1))/d[i];
					}
				}
			}
		}		
	}
}

