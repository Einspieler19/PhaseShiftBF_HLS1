
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

#define W_IN    48
#define IW_IN   24
#define N 1001
#define NUM_ELEMENTS 7





typedef ap_fixed<W_IN,IW_IN> data_t;
typedef float matrix_t;

void PhaseshiftBeamformer(
    matrix_t rx_re[N],
    matrix_t rx_im[N],
    matrix_t incidentAngle,
    matrix_t y_re[N],
    matrix_t y_im[N]);


#endif

