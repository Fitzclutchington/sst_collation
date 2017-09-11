const int WIDTH = 5500;
const int HEIGHT = 5500;
const int WIDTH_HIST = 101;
const int HEIGHT_HIST = 144;
const int DEPTH_HIST = 90;

const int SAMPLING_RATE = 6; //samples per hour: 10 min for himawari
const int TIMEINT = 144;
const int MAXDIMS = 5;
const int TIME_SIZE = 144;
const float PI = 3.14159265358979323;

const bool DEBUG = true;
const float PR_CLEAR = 0.1;
const int   MIN_POINTS = 6;
const float DELTA_SST = 0.1;

const int DIAG_SIZE = 3;
const int DIAG_LAG = DIAG_SIZE/2;

const int FILTER_TIME_SIZE = 5;
const int FILTER_WINDOW_LAG = FILTER_TIME_SIZE/2;
const int COLLATED_SMOOTH_LAG = 7;

const float T_DIAG = 2*DELTA_SST*sqrt(DIAG_SIZE);
const float T_NN = 0.5;//0.35;//2*DELTA_SST;
const float MIN_TEMP = 271.15;
const float T_COLD_DT = -5;
const float T_WARM_DT = 10;
const float T_INTERP = 1.2;//0.5; //change to .3

const int COLLATED_LAG = 6;
const float GAMMA2 = 1;//COLLATED_LAG;
const float GAMMA = 6;

const float MU_CLEAR = 0.9;
const float MU_APPROX = 1 - MU_CLEAR;

const int LAND_KERNEL = 7;

const float T_DERIV = 0.04;
const float T_EIGEN = 0.1;
const float T_LS = 0.5;
const float T_ABOVE = 0.4;

const int PASS_THRESH = 1000;
const int INTERP_DIST = 7;
const float MAX_SZA = 55.0;
const float EXP = 1.0/30;
const float SIGMA_H = 1.0/3;

// Bit mask for each cloud mask
const uchar L2P = 1;
const uchar NN = 2;
const uchar DIAG = 4;
const uchar EIGEN = 8;
const uchar LS = 16;
const uchar ACSPO = 32;