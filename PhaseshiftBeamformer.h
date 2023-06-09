
#ifndef PSBF_H
#define PSBF_H


#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <random>

#include <ap_int.h>     // Xilinx hls的整型库
#include <ap_fixed.h>

#include <hls_math.h>

//#include "constants.h"
#include "Template.h"


using namespace std;



//typedef float data_psb;
typedef ap_fixed<W_IN,IW_IN> data_psb;
typedef float data_tb;

//typedef ap_fixed<W_IN,IW_IN> data_typeB;



void PhaseshiftBeamformer(
		data_psb cov_Mat_re[NUMELEMENTS],
		data_psb cov_Mat_im[NUMELEMENTS],
		data_psb weightsRe[NUMELEMENTS],
		data_psb weightsIm[NUMELEMENTS],
		data_psb *outRe,
		data_psb *outIm);



#endif


