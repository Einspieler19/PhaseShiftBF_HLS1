
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


#include <hls_stream.h>
#include <hls_math.h>


//#include "constants.h"
#include "Template.h"


using namespace std;




typedef ap_fixed<W_IN,IW_IN> data_psb;
typedef float data_tb;

typedef float data_typeA;
typedef ap_fixed<W_IN,IW_IN> data_typeB;




void PhaseshiftBeamformer();



#endif

