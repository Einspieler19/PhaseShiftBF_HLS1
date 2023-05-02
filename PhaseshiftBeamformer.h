
#ifndef PSBF_H
#define PSBF_H


#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include <ap_int.h>     // Xilinx hls的整型库
#include <ap_fixed.h>

#include <hls_math.h>

//#include "constants.h"
#include "Template.h"


using namespace std;



typedef float data_psb;
//typedef ap_fixed<W_IN,IW_IN> data_psb;
typedef float data_tb;

typedef float data_typeA;
typedef float data_typeB;
//typedef ap_fixed<W_IN,IW_IN> data_typeB;

struct my_Output{
data_psb a;
data_psb b;
};


my_Output PhaseshiftBeamformer(
		data_psb cov_Mat_re[NUMELEMENTS],
		data_psb cov_Mat_im[NUMELEMENTS],
		data_psb weightsRe[NUMELEMENTS],
		data_psb weightsIm[NUMELEMENTS]);



#endif


