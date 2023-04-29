
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


#include "constants.h"


using namespace std;


#define FLOAT_DATA // Used to select error tolerance in test program

#define W_IN    64
#define IW_IN   32
#define N 1001
#define NUM_ELEMENTS 7





typedef ap_fixed<W_IN,IW_IN> data_psb;
typedef float data_tb;

typedef float data_typeA;
typedef ap_fixed<W_IN,IW_IN> data_typeB;



void PhaseshiftBeamformer(
		data_psb (*cov_Mat_re)[NUM_ELEMENTS],
		data_psb (*cov_Mat_im)[NUM_ELEMENTS],
		data_psb steeringAngle,
		data_psb* weightsRe,
		data_psb* weightsIm,
		data_psb* y_re,
		data_psb* y_im);



#endif

