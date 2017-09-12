
#define SAVENC(X)	if(DEBUG)savenc(#X ".nc", (X))
#define CHECKMAT(M, T) CV_Assert((M).type() == (T) && (M).isContinuous())
//#define CHECKMAT(M) CV_Assert((M).isContinuous())



void
eprintf(const char *fmt, ...)
{
	va_list args;

	fflush(stdout);
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	if(fmt[0] != '\0' && fmt[strlen(fmt)-1] == ':'){
		fprintf(stderr, " %s", strerror(errno));
	}
	fprintf(stderr, "\n");

	exit(2);
}

// Print out error for NetCDF error number n and exit the program.
void
ncfatal(int n, const char *fmt, ...)
{
	va_list args;

	fflush(stdout);
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	fprintf(stderr, ": %s\n", nc_strerror(n));
	exit(2);
}

char*
fileprefix(const char *path)
{
	const char *b;
	char *p, *s;
	
	b = strrchr(path, '/');
	if(b == nullptr){
		b = path;
	}else{
		b++;
	}
	p = strdup(b);
	s = strrchr(p, '.');
	if(s != nullptr){
		*s = '\0';
	}
	return p;
}

const char*
type2str(int type)
{
	switch(type){
	default:       return "UnknownType";
	case CV_8UC1:  return "CV_8UC1";
	case CV_8SC1:  return "CV_8SC1";
	case CV_16UC1: return "CV_16UC1";
	case CV_16SC1: return "CV_16SC1";
	case CV_32SC1: return "CV_32SC1";
	case CV_32FC1: return "CV_32FC1";
	case CV_64FC1: return "CV_64FC1";
	}
}

void
savenc(const char *path, const Mat &mat, bool compress=false)
{
	int i, n, ncid, dims, varid, xtype;
	int dimids[MAXDIMS] = {};
	char *name;
	const char *dimnames[MAXDIMS] = {
		"dim0",
		"dim1",
		"dim2",
		"dim3",
		"dim4",
	};
	
	dims = mat.dims;
	if(mat.channels() > 1){
		dims++;
	}
	if(dims > MAXDIMS){
		eprintf("savenc: too many dimensions %d\n", dims);
	}
	
	n = nc_create(path, NC_NETCDF4, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "savenc: creating %s failed", path);
	}
	for(i = 0; i < mat.dims; i++){
		n = nc_def_dim(ncid, dimnames[i], mat.size[i], &dimids[i]);
		if(n != NC_NOERR){
			ncfatal(n, "savenc: creating dim %d failed", i);
		}
	}
	if(mat.channels() > 1){
		n = nc_def_dim(ncid, dimnames[i], mat.channels(), &dimids[i]);
		if(n != NC_NOERR){
			ncfatal(n, "savenc: creating dim %d failed", i);
		}
	}
	
	xtype = -1;
	switch(mat.depth()){
	default:
		eprintf("savenc: unsupported type %s\n", type2str(mat.type()));
		break;
	case CV_8U:	xtype = NC_UBYTE; break;
	case CV_8S:	xtype = NC_BYTE; break;
	case CV_16U:	xtype = NC_USHORT; break;
	case CV_16S:	xtype = NC_SHORT; break;
	case CV_32S:	xtype = NC_INT; break;
	case CV_32F:	xtype = NC_FLOAT; break;
	case CV_64F:	xtype = NC_DOUBLE; break;
	}
	
	n = nc_def_var(ncid, "data", xtype, dims, dimids, &varid);
	if(n != NC_NOERR){
		ncfatal(n, "savenc: creating variable failed");
	}
	if(compress){	// enable compression?
		n = nc_def_var_deflate(ncid, varid, 0, 1, 1);
		if(n != NC_NOERR){
			ncfatal(n, "savenc: setting deflate parameters failed");
		}
	}

	n = nc_put_var(ncid, varid, mat.data);
	if(n != NC_NOERR){
		ncfatal(n, "savenc: writing variable failed");
	}
	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "savenc: closing %s failed", path);
	}

	name = fileprefix(path);
	printf("%s = loadnc(\"%s\")\n", name, path);
	free(name);
}

int
readvar(int ncid, const char *name, Mat &img)
{
	int i, varid, n, ndims, dimids[MAXDIMS], ishape[MAXDIMS], cvt;
	size_t shape[MAXDIMS] = {};
	nc_type nct;
	
	n = nc_inq_varid(ncid, name, &varid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_inq_varid failed for variable %s", name);
	}
	n = nc_inq_var(ncid, varid, nullptr, &nct, &ndims, dimids, nullptr);
	if(n != NC_NOERR){
		ncfatal(n, "nc_inq_var failed for variable %s", name);
	}
	if(ndims > MAXDIMS){
		eprintf("number of dimensions %d > MAXDIMS=%d\n", ndims, MAXDIMS);
	}
	
	for(i = 0; i < ndims; i++){
		n = nc_inq_dimlen(ncid, dimids[i], &shape[i]);
		if(n != NC_NOERR){
			ncfatal(n, "nc_inq_dimlen failed for dim %d", dimids[i]);
		}
	}
	
	cvt = -1;
	switch(nct){
	default:
		eprintf("unknown netcdf data type");
		break;
	case NC_BYTE:	cvt = CV_8SC1; break;
	case NC_UBYTE:	cvt = CV_8UC1; break;
	case NC_SHORT:	cvt = CV_16SC1; break;
	case NC_USHORT:	cvt = CV_16UC1; break;
	case NC_INT:	cvt = CV_32SC1; break;
	case NC_FLOAT:	cvt = CV_32FC1; break;
	case NC_DOUBLE:	cvt = CV_64FC1; break;
	}
	
	for(i = 0; i < MAXDIMS; i++){
		ishape[i] = shape[i];
	}
	img = Mat(ndims, ishape, cvt);
	n = nc_get_var(ncid, varid, img.data);
	if(n != NC_NOERR){
		ncfatal(n, "readvar: nc_get_var '%s' failed", name);
	}
	return varid;
}

void
get_file_names(const char *pathsfile, vector<string> &paths)
{
	std::ifstream f(pathsfile);
	std::string line; 
	while (std::getline(f, line)){
		paths.push_back(line);

	}
	f.close();

	if(paths.size() < TIMEINT){
		eprintf("You must have at least %d granules\n",TIMEINT);
	}
	
}

int
get_current_files(const vector<string> &paths, vector<string> &curPaths, int nfiles, int loc=0)
{
	curPaths.clear();
	for(int i=0;i<nfiles;i++){
		curPaths.push_back(paths[loc+i]);
	}
	return loc + 1;
}

void
get_icemask(const string pathsfile, Mat1b &ice_mask)
{
	int ncid,y,x;
	int n = nc_open(pathsfile.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s\n", pathsfile.c_str());
	}

	Mat1b img;
	readvar(ncid, "l2p_flags", img);
	if(img.dims != 2 || img.size[0] != HEIGHT || img.size[1] != HEIGHT){
		eprintf("unpexpected dimensions\n");
	}
	if(img.type() != CV_8UC1){
		eprintf("unpexpected type\n");
	}
    
    ice_mask = img & (1 << 2);

	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}
}

void
open_LUT(const string pathsfile, Mat1b &lut, int *dims)
{
	int ncid, i,j,k,l;
	int n = nc_open(pathsfile.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s\n", pathsfile.c_str());
	}

	Mat1b img;
	readvar(ncid, "LUT", img);
	if(img.dims != 4 || img.size[0] != dims[0] || img.size[1] != dims[1] || img.size[2] != dims[2] || img.size[3] != dims[3]){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_8UC1){
		eprintf("unpexpected type\n");
	}
    
    int d[4] = {0,0,0,0};
	for(i = 0; i < dims[0]; i++){
		d[0] = i;
		for(j = 0; j < dims[1]; j++){
			d[1] =j;
			for(k = 0; k < dims[2]; k++){
				d[2] = k;
				for(l = 0; l < dims[3]; l++){	
					d[3] =l;
					lut(d) = img(d);
				}
			}		 
		}
	}
	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}
}

void
read_mask(const string pathsfile, Mat1b &mask, int cur_ind)
{
	int ncid, y,x;
	int n = nc_open(pathsfile.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s\n", pathsfile.c_str());
	}

	Mat1b img;
	readvar(ncid, "brightness_temperature_11um2", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != WIDTH){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_8UC1){
		eprintf("unpexpected type\n");
	}
   
    if(cur_ind < 0){
		for(y = 0; y < HEIGHT; y++){		
			for(x = 0; x < WIDTH; x++){
			
					//printf("val of img = %d", img(0,y,x));	

				mask(y,x) = img(0,y,x);
			}
		}
	}
	else{
		for(y = 0; y < HEIGHT; y++){		
			for(x = 0; x < WIDTH; x++){
			
					//printf("val of img = %d", img(0,y,x));	

				mask(y,x,cur_ind) = img(0,y,x);
			}
		}
	}		 

	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}
}

void
read_mask_binary(const string pathsfile, Mat1s &mask, int cur_ind)
{
	int ncid, y,x;
	int n = nc_open(pathsfile.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s\n", pathsfile.c_str());
	}

	Mat1s img;
	readvar(ncid, "brightness_temperature_11um2", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != WIDTH){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
   
    if(cur_ind < 0){
		for(y = 0; y < HEIGHT; y++){		
			for(x = 0; x < WIDTH; x++){
			
					//printf("val of img = %d", img(0,y,x));	

				mask(y,x) = img(0,y,x);
			}
		}
	}
	else{
		for(y = 0; y < HEIGHT; y++){		
			for(x = 0; x < WIDTH; x++){
			
					//printf("val of img = %d", img(0,y,x));	

				mask(y,x,cur_ind) = img(0,y,x);
			}
		}
	}		 

	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}
}

void
get_l2pmask(const char *pathsfile, Mat1b &land_mask, Mat1b &l2p_mask)
{
    int ncid;
	int n = nc_open(pathsfile, 0, &ncid);
	short val;

	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s\n", pathsfile);
	}

	Mat1s img;
	readvar(ncid, "l2p_flags", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}

    for(int y = 0; y < HEIGHT; ++y){
		for(int x = 0; x < WIDTH; ++x){
			val = img(0, y, x);
			if((val & 2)){
				land_mask(y, x) = 255;
				l2p_mask(y, x) = 255;
			}
			else{
				land_mask(y, x) = 0;
				l2p_mask(y,x) = 0;
			}
		}
	}

	for(int y = 0; y < HEIGHT; ++y){
		for(int x = 0; x < WIDTH; ++x){
			val = img(0, y, x);
			if((val & 256) || (val & (1 << 2))){
				l2p_mask(y, x) = 255;
			}
			else{
				l2p_mask(y,x) = 0;
			}
		}
	}
	
	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}
}

void
read_acspo(const string pathsfile, Mat1b &mask, int cur_ind)
{
	int ncid, y,x;
	int n = nc_open(pathsfile.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s\n", pathsfile.c_str());
	}

	Mat1s img;
	readvar(ncid, "l2p_flags", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
    
    if(cur_ind < 0){
		for(y = 0; y < HEIGHT; y++){		
			for(x = 0; x < WIDTH; x++){
				if(img(0,y,x) & -16384){
				 	mask(y, x) = mask(y,x) | ACSPO;
				 }
				 else{
				 	mask(y,x) =  0;
				 }
			}
		}
	}
	else{
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 if(img(0,y,x) & -16384){
				 	mask(y, x, cur_ind) = mask(y,x,cur_ind) | ACSPO;
				 }
				 else{
				 	mask(y,x,cur_ind) = 0;
				 }
			}
		}
	}
	

	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}
}

void
read_reference(const string path,Mat1f &data,string variable)
{
	//printf("%s\n", paths[i].c_str());
	float val;
	int ncid;
	int n = nc_open(path.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s", path.c_str());
	}
		
	Mat img;

	readvar(ncid, variable.c_str(), img);
	if(img.dims != 2 || img.size[0] != HEIGHT|| img.size[1] != WIDTH){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_32FC1){
		eprintf("unpexpected type\n");
	}
	
    
	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 val = img.at<float>(y, x);
			 data(y,x) = NAN;
			 if(!std::isnan(val)){
			 	data(y, x) = val;
			 }
		}
	}
	

	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}
}
void
readgranule_oneband(const string path,Mat1f &bt11, int ind,string variable)
{
	// TODO: deal with non-continuous time interval

	//printf("%s\n", paths[i].c_str());
	float val;
	int ncid;
	int n = nc_open(path.c_str(), 0, &ncid);
	float scale,offset;
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s", path.c_str());
	}
		
	Mat img;

	int varid = readvar(ncid, variable.c_str(), img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT|| img.size[2] != WIDTH){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_32FC1 && img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
    if(img.type() == CV_32FC1){
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<float>(0,y, x);
				 bt11(y,x,ind) = NAN;
				 if(val > 0 && val < 400){
				 	bt11(y, x, ind) = val;
				 }
			}
		}
	}
	else{
		n = nc_get_att_float(ncid, varid, "add_offset", &offset);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed add_offset");
		}
		n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed scale");
		}

		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 float val = img.at<short>(0,y, x)*scale + offset;
				 bt11(y,x,ind) = NAN;
				 if(val > 0 && val < 400){
				 	bt11(y, x, ind) = val;
				 }
			}
		}
	}

	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}
}

void
get_var(const string path, Mat1f &mat, const string variable)
{
	// TODO: deal with non-continuous time interval
    
	//printf("%s\n", paths[i].c_str());
	int ncid;
	float val;
	float offset, scale;
	int n = nc_open(path.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s", path.c_str());
	}
	
	Mat img;

	int varid = readvar(ncid, variable.c_str(), img);
	if(img.dims != 2 && img.dims != 3){
		if(img.dims ==2 && (img.size[0] != HEIGHT|| img.size[1] != WIDTH)){
			printf("unpexpected dimensions\n");
		}
		if(img.dims ==3 && (img.size[0] != 1 || img.size[1] != HEIGHT|| img.size[2] != WIDTH)){
			printf("unpexpected dimensions\n");
		}
	}
	if(img.type() != CV_32FC1 && img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	if(img.type() == CV_16SC1){
		n = nc_get_att_float(ncid, varid, "add_offset", &offset);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 if(img.dims == 3){
				 	val = img.at<short>(0, y, x)*scale + offset;
				 }
				 else{
				 	val = img.at<short>(y, x)*scale + offset;
				 }
				 mat(y,x) = NAN;
				 if(val > -200){
				 	mat(y, x) = val;
				 }
			}
		}
	}

	else{
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 if(img.dims == 3){
				 	val = img.at<float>(0, y, x);
				 }
				 else{
				 	val = img.at<float>(y, x);
				 }
				 mat(y,x) = NAN;
				 if(val > -200){
				 	mat(y, x) = val;
				 }
			}
		}
	}

	
	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}

}


void
open_hist(const string path, Mat1f &mat)
{
	// TODO: deal with non-continuous time interval
    
	//printf("%s\n", paths[i].c_str());
	int ncid;
	int n = nc_open(path.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s", path.c_str());
	}
	
	Mat img;
	printf("in function\n");
	readvar(ncid, "histogram", img);

	if(img.dims ==3 && (img.size[0] != 1 || img.size[1] != HEIGHT_HIST|| img.size[2] != WIDTH_HIST)){
		printf("unpexpected dimensions\n");
	}
	
	if(img.type() != CV_32FC1 && img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	printf("starting loop\n");
	for(int y = 0; y < HEIGHT_HIST; y++){
		for(int x = 0; x < WIDTH_HIST; x++){
			 
			 mat(y, x) = img.at<float>(0, y, x);

		}
	}
	
	printf("finished loop\n");
	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}

}

void
open_hist_3d(const string path, Mat1f &mat)
{
	// TODO: deal with non-continuous time interval
    
	//printf("%s\n", paths[i].c_str());
	int ncid;
	int n = nc_open(path.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s", path.c_str());
	}
	
	Mat1f img;
	printf("in function\n");
	readvar(ncid, "data", img);

	if(img.dims ==3 && ( img.size[0] != HEIGHT_HIST|| img.size[1] != WIDTH_HIST || img.size[2] != DEPTH_HIST)){
		printf("unpexpected dimensions\n");
	}
	
	if(img.type() != CV_32FC1){
		eprintf("unpexpected type\n");
	}
	
	printf("starting loop\n");
	for(int y = 0; y < HEIGHT_HIST; y++){
		
		for(int x = 0; x < WIDTH_HIST; x++){
			
			 for(int t = 0; t < DEPTH_HIST; t++){
			 	
			 	mat(y, x,t) = img(y,x,t);
			 }

		}
	}
	
	printf("finished loop\n");
	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}

}

void
save_test_nc_float(const Mat samples, const char *filename){
	string width_name = "ni";
	string height_name = "nj";
    std::vector<int> dimid(3);
    
    int ncid;
	int n = nc_create(filename,NC_NETCDF4,&ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_create failed");
	}

	n = nc_def_dim(ncid, "time", 1, &dimid[0]);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_dim failed");
	}

	n = nc_def_dim(ncid, height_name.c_str(), HEIGHT, &dimid[1]);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_dim failed");
	}

	n = nc_def_dim(ncid, width_name.c_str(), WIDTH, &dimid[2]);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_dim failed");
	}
    
    int varid;
    if(samples.type()== CV_32FC1){
	    n = nc_def_var(ncid, "brightness_temperature_11um2", NC_FLOAT, dimid.size(), dimid.data(), &varid);
		if(n != NC_NOERR){
			ncfatal(n, "nc_def_var failed");
		}

	
		n = nc_put_var_float(ncid, varid, (float*)samples.data);
		if(n != NC_NOERR)
			ncfatal(n, "nc_put_var_float failed");
	}

	else if(samples.type() == CV_8UC1){
	    n = nc_def_var(ncid, "brightness_temperature_11um2", NC_UBYTE, dimid.size(), dimid.data(), &varid);
		if(n != NC_NOERR){
			ncfatal(n, "nc_def_var failed");
		}

		n = nc_put_var_uchar(ncid, varid, (uchar*)samples.data);
		if(n != NC_NOERR)
			ncfatal(n, "nc_put_var_float failed");
	}

	else if(samples.type() == CV_16SC1){
	    n = nc_def_var(ncid, "brightness_temperature_11um2", NC_SHORT, dimid.size(), dimid.data(), &varid);
		if(n != NC_NOERR){
			ncfatal(n, "nc_def_var failed");
		}

		n = nc_put_var_short(ncid, varid, (short*)samples.data);
		if(n != NC_NOERR)
			ncfatal(n, "nc_put_var_float failed");
	}


	else{
		eprintf("Cannot save this data type\n");
	}

    n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "savenc: closing %s failed", filename);
	}
}

void
readgranule(const string path, Mat1f &bt11, Mat1f &bt12, Mat1f &bt08, Mat1f &bt10, int ind)
{
	// TODO: deal with non-continuous time interval 
	int ncid;
	float val;

	int n = nc_open(path.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s", path.c_str());
	}
	
	Mat1s img;

	int varid = readvar(ncid, "brightness_temperature_11um2", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	float scale, offset;
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 bt11(y,x,ind) = NAN;
			 val = img(0, y, x)*scale + offset;
			 if(val > 0){
			 	bt11(y, x,ind) = val;
			 }
		}
	}
	
	varid = readvar(ncid, "brightness_temperature_12um3", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 val = img(0, y, x)*scale + offset;
			 bt12(y,x,ind) = NAN;
			 if(val > 0){
			 	bt12(y, x,ind) = val;
			 }
		}
	}

	varid = readvar(ncid, "brightness_temperature_08um6", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 val = img(0, y, x)*scale + offset;
			 bt08(y,x,ind) = NAN;
			 if(val > 0){
			 	bt08(y, x,ind) = val;
			 }
		}
	}
	
	varid = readvar(ncid, "brightness_temperature_10um4", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 val = img(0, y, x)*scale + offset;
			 bt10(y,x,ind) = NAN;
			 if(val > 0){
			 	bt10(y, x,ind) = val;
			 }
		}
	}
	/*
	varid = readvar(ncid, "sea_surface_temperature", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 val = img(0, y, x)*scale + offset;
			 sst(y,x,ind) = NAN;
			 if(val > 0){
			 	sst(y, x,ind) = val;
			 }
		}
	}
    */
	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}

}

void
readgranule_1d(const string path, Mat1f &bt08, Mat1f &bt10, Mat1f &bt11, Mat1f &bt12)
{
	// TODO: deal with non-continuous time interval 
	int ncid;
	float val;

	int n = nc_open(path.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s", path.c_str());
	}
	
	Mat1s img;

	int varid = readvar(ncid, "brightness_temperature_11um2", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	float scale, offset;
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 bt11(y,x) = NAN;
			 val = img(0, y, x)*scale + offset;
			 if(val > 0){
			 	bt11(y, x) = val;
			 }
		}
	}
	
	varid = readvar(ncid, "brightness_temperature_12um3", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 val = img(0, y, x)*scale + offset;
			 bt12(y,x) = NAN;
			 if(val > 0){
			 	bt12(y, x) = val;
			 }
		}
	}

	varid = readvar(ncid, "brightness_temperature_08um6", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 val = img(0, y, x)*scale + offset;
			 bt08(y,x) = NAN;
			 if(val > 0){
			 	bt08(y, x) = val;
			 }
		}
	}
	
	varid = readvar(ncid, "brightness_temperature_10um4", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1){
		eprintf("unpexpected type\n");
	}
	
	n = nc_get_att_float(ncid, varid, "add_offset", &offset);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			 val = img(0, y, x)*scale + offset;
			 bt10(y,x) = NAN;
			 if(val > 0){
			 	bt10(y, x) = val;
			 }
		}
	}
    
	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}

}

void
readgranule_fullbands(const string path, Mat1f &bt11, Mat1f &bt12, Mat1f &bt08, Mat1f &bt10, Mat1f &sst, int ind)
{
	// TODO: deal with non-continuous time interval 
	
	int ncid;
	float val;
    float scale, offset;

	int n = nc_open(path.c_str(), 0, &ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_open failed for %s", path.c_str());
	}
	
	Mat img;
	
	int varid = readvar(ncid, "brightness_temperature_11um2", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1 && img.type() != CV_32FC1){
		eprintf("unpexpected type\n");
	}

	if(img.type() == CV_16SC1){
		
		n = nc_get_att_float(ncid, varid, "add_offset", &offset);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
	
	
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 bt11(y,x,ind) = NAN;
				 val = img.at<short>(0, y, x)*scale + offset;
				 if(val > 0){
				 	bt11(y, x,ind) = val;
				 }
			}
		}
	}

	else{
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 bt11(y,x,ind) = NAN;
				 val = img.at<float>(0, y, x);
				 if(val > 0){
				 	bt11(y, x,ind) = val;
				 }
			}
		}
	}
		
	varid = readvar(ncid, "brightness_temperature_12um3", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1 && img.type() != CV_32FC1){
		eprintf("unpexpected type\n");
	}
	
	if(img.type() == CV_16SC1){
		n = nc_get_att_float(ncid, varid, "add_offset", &offset);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		

		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<short>(0, y, x)*scale + offset;
				 bt12(y,x,ind) = NAN;
				 if(val > 0){
				 	bt12(y, x,ind) = val;
				 }
			}
		}
	}

	else{
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<float>(0, y, x);
				 bt12(y,x,ind) = NAN;
				 if(val > 0){
				 	bt12(y, x,ind) = val;
				 }
			}
		}
	}

	varid = readvar(ncid, "brightness_temperature_08um6", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1 && img.type() != CV_32FC1){
		eprintf("unpexpected type\n");
	}
	
	if(img.type() == CV_16SC1){
		n = nc_get_att_float(ncid, varid, "add_offset", &offset);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		

		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<short>(0, y, x)*scale + offset;
				 bt08(y,x,ind) = NAN;
				 if(val > 0){
				 	bt08(y, x,ind) = val;
				 }
			}
		}
	}

	else{
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<float>(0, y, x);
				 bt08(y,x,ind) = NAN;
				 if(val > 0){
				 	bt08(y, x,ind) = val;
				 }
			}
		}
	}
	
	varid = readvar(ncid, "brightness_temperature_10um4", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1 && img.type() != CV_32FC1){
		eprintf("unpexpected type\n");
	}
	if(img.type() == CV_16SC1){
		n = nc_get_att_float(ncid, varid, "add_offset", &offset);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		

		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<short>(0, y, x)*scale + offset;
				 bt10(y,x,ind) = NAN;
				 if(val > 0){
				 	bt10(y, x,ind) = val;
				 }
			}
		}
	}

	else{
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<float>(0, y, x);
				 bt10(y,x,ind) = NAN;
				 if(val > 0){
				 	bt10(y, x,ind) = val;
				 }
			}
		}
	}

	varid = readvar(ncid, "sea_surface_temperature", img);
	if(img.dims != 3 || img.size[0] != 1 || img.size[1] != HEIGHT || img.size[2] != HEIGHT){
		printf("unpexpected dimensions\n");
	}
	if(img.type() != CV_16SC1 && img.type() != CV_32FC1){
		eprintf("unpexpected type\n");
	}
	if(img.type() == CV_16SC1){
		n = nc_get_att_float(ncid, varid, "add_offset", &offset);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		n = nc_get_att_float(ncid, varid, "scale_factor", &scale);
		if(n != NC_NOERR){
			ncfatal(n, "nc_get_att_float failed");
		}
		

		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<short>(0, y, x)*scale + offset;
				 sst(y,x,ind) = NAN;
				 if(val > 0){
				 	sst(y, x,ind) = val;
				 }
			}
		}
	}

	else{
		for(int y = 0; y < HEIGHT; y++){
			for(int x = 0; x < WIDTH; x++){
				 val = img.at<float>(0, y, x);
				 sst(y,x,ind) = NAN;
				 if(val > 0){
				 	sst(y, x,ind) = val;
				 }
			}
		}
	}

	n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_close failed");
	}

}

void
save_test_nc_fullbands(const Mat1f &bt08,const Mat1f & bt10,const Mat1f & bt11, const Mat1f &bt12,
	                   const Mat1f &sst,const char* filename)
{
	string width_name = "ni";
	string height_name = "nj";
    std::vector<int> dimid(3);
    
    int ncid;
	int n = nc_create(filename,NC_NETCDF4,&ncid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_create failed");
	}

	n = nc_def_dim(ncid, "time", 1, &dimid[0]);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_dim failed");
	}

	n = nc_def_dim(ncid, height_name.c_str(), HEIGHT, &dimid[1]);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_dim failed");
	}

	n = nc_def_dim(ncid, width_name.c_str(), WIDTH, &dimid[2]);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_dim failed");
	}
    
    int varid;
    n = nc_def_var(ncid, "brightness_temperature_10um4", NC_FLOAT, dimid.size(), dimid.data(), &varid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_var failed");
	}

	n = nc_put_var_float(ncid, varid, (float*)bt10.data);
	if(n != NC_NOERR)
		ncfatal(n, "nc_put_var_float failed");

	n = nc_def_var(ncid, "brightness_temperature_08um6", NC_FLOAT, dimid.size(), dimid.data(), &varid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_var failed");
	}

	n = nc_put_var_float(ncid, varid, (float*)bt08.data);
	if(n != NC_NOERR)
		ncfatal(n, "nc_put_var_float failed");

	n = nc_def_var(ncid, "brightness_temperature_11um2", NC_FLOAT, dimid.size(), dimid.data(), &varid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_var failed");
	}

	n = nc_put_var_float(ncid, varid, (float*)bt11.data);
	if(n != NC_NOERR)
		ncfatal(n, "nc_put_var_float failed");
	
	n = nc_def_var(ncid, "brightness_temperature_12um3", NC_FLOAT, dimid.size(), dimid.data(), &varid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_var failed");
	}

	n = nc_put_var_float(ncid, varid, (float*)bt12.data);
	if(n != NC_NOERR)
		ncfatal(n, "nc_put_var_float failed");

	n = nc_def_var(ncid, "sea_surface_temperature", NC_FLOAT, dimid.size(), dimid.data(), &varid);
	if(n != NC_NOERR){
		ncfatal(n, "nc_def_var failed");
	}

	n = nc_put_var_float(ncid, varid, (float*)sst.data);
	if(n != NC_NOERR)
		ncfatal(n, "nc_put_var_float failed");

    n = nc_close(ncid);
	if(n != NC_NOERR){
		ncfatal(n, "savenc: closing %s failed", filename);
	}
}

void
save_and_update(const string filename, const Mat &samples,string variable, bool mode, int height, int width)
{
	//mode true- create file
	//mode false - append to file
	string width_name = "ni";
	string height_name = "nj";
    std::vector<int> dimid(3);
    
    int ncid,n, varid;
    if(mode==true){
		n = nc_create(filename.c_str(),NC_NETCDF4,&ncid);
		if(n != NC_NOERR){
			ncfatal(n, "nc_create failed");
		}

		n = nc_def_dim(ncid, "time", 1, &dimid[0]);
		if(n != NC_NOERR){
			ncfatal(n, "nc_def_dim failed");
		}

		n = nc_def_dim(ncid, height_name.c_str(), height, &dimid[1]);
		if(n != NC_NOERR){
			ncfatal(n, "nc_def_dim failed");
		}

		n = nc_def_dim(ncid, width_name.c_str(), width, &dimid[2]);
		if(n != NC_NOERR){
			ncfatal(n, "nc_def_dim failed");
		}


		if(samples.type()== CV_32FC1){
		    n = nc_def_var(ncid, variable.c_str(), NC_FLOAT, dimid.size(), dimid.data(), &varid);
			if(n != NC_NOERR){
				ncfatal(n, "nc_def_var failed");
			}

		
			n = nc_put_var_float(ncid, varid, (float*)samples.data);
			if(n != NC_NOERR)
				ncfatal(n, "nc_put_var_float failed");
		}

		else if(samples.type() == CV_8UC1){
		    n = nc_def_var(ncid, variable.c_str(), NC_UBYTE, dimid.size(), dimid.data(), &varid);
			if(n != NC_NOERR){
				ncfatal(n, "nc_def_var failed");
			}

			n = nc_put_var_uchar(ncid, varid, (uchar*)samples.data);
			if(n != NC_NOERR)
				ncfatal(n, "nc_put_var_float failed");
		}

	    n = nc_close(ncid);
		if(n != NC_NOERR){
			ncfatal(n, "savenc: closing %s failed", filename.c_str());
		}

	}
	else{
		n = nc_open(filename.c_str(), NC_WRITE, &ncid);
		if(n != NC_NOERR){
			ncfatal(n, "nc_open failed for %s", filename.c_str());
		}

		dimid[0] = 0; dimid[1] = 1; dimid[2] = 2;
	    n = nc_def_var(ncid, variable.c_str(), NC_FLOAT, dimid.size(), dimid.data() , &varid);
		if(n != NC_NOERR){
			ncfatal(n, "nc_def_var failed");
		}

		n = nc_put_var_float(ncid, varid, (float*)samples.data);
		if(n != NC_NOERR)
			ncfatal(n, "nc_put_var_float failed");

	    n = nc_close(ncid);
		if(n != NC_NOERR){
			ncfatal(n, "savenc: closing %s failed", filename.c_str());
		}
	}


}


void
save_mat(vector<string> paths, Mat &samples, string variable,bool create, int height, int width)
{
	int i,y,x;
	int time_size = paths.size();
	Mat save_slice;
	if(samples.type()== CV_32FC1){
		save_slice.create(height,width,CV_32FC1);
	}

	else if(samples.type() == CV_8UC1){
		save_slice.create(height,width,CV_8UC1);
	}
	else{
		eprintf("WRONG TYPE\n");
	}
	for(i=0;i<time_size;i++){


        if(samples.type()== CV_32FC1){
			// place slice from 3d matrix into 2d matrix in order to save
	        for(y=0;y<height;y++){
	            for(x=0;x<width;x++){
	                save_slice.at<float>(y,x) = samples.at<float>(y,x,i);
	            }
	        }
		}

		else if(samples.type() == CV_8UC1){
		    // place slice from 3d matrix into 2d matrix in order to save
	        for(y=0;y<height;y++){
	            for(x=0;x<width;x++){
	                save_slice.at<uchar>(y,x) = samples.at<uchar>(y,x,i);
	            }
	        }
		}

        save_and_update(paths[i], save_slice,variable, create, height, width);
        printf("generated file %s with variable %s\n", paths[i].c_str(),variable.c_str());

    }
}

