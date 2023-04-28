#include "PhaseshiftBeamformer.h"


void PhaseshiftBeamformer(
    matrix_t rx_re[N],
    matrix_t rx_im[N],
    matrix_t incidentAngle,
    matrix_t y_re[N],
    matrix_t y_im[N])
{
	// 位置初始化
	matrix_t Elementpos[NUM_ELEMENTS] = {-1.5, -1, -0.5, 0.0, 0.5, 1, 1.5};
    
// 权重
    matrix_t w2_re[NUM_ELEMENTS], w2_im[NUM_ELEMENTS];

    matrix_t c = 299792458;
    matrix_t fc = 300000000;

// w是复数，分实部和虚部计算
    for (int i = 0; i < NUM_ELEMENTS; i++) {
    	//float a =hls::cos(2 * M_PI * hls::sinf(M_PI/3) * 0 * c / fc)/NUM_ELEMENTS;
    	//cout<<"a:"<<a<<endl;
    	w2_re[i] = hls::cos(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[i] * c / fc);
        w2_im[i] = hls::sin(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[i] * c / fc);

    }

    // Perform beamforming
    for (int i = 0; i < N; i++) {
        y_re[i] = 0.0;
        y_im[i] = 0.0;
 // 初始化
        for (int j = 0; j < NUM_ELEMENTS; j++) {
            matrix_t rx_w2_re = rx_re[i] * w2_re[j] - rx_im[i] * w2_im[j]; // 实部
            matrix_t rx_w2_im = rx_re[i] * w2_im[j] + rx_im[i] * w2_re[j]; // 虚部
            //cout<<"rx_re:"<<rx_re[i]<<endl;
           // cout<<"rx_im:"<<rx_im[i]<<endl;
            //cout<<"w2_re:"<<w2_re[j]<<endl;
           // cout<<"w2_im:"<<w2_im[j]<<endl;
           // cout<<"rx_w2_im:"<<rx_w2_im<<endl;
 // j从1到7：天线入射的叠加
            y_re[i] += rx_w2_re;
            y_im[i] += rx_w2_im;

        }

 // 除以7：功率归一
        y_re[i] /= NUM_ELEMENTS;
        y_im[i] /= NUM_ELEMENTS;

    }
}
