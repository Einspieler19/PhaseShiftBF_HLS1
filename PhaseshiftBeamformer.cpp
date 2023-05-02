#include "PhaseshiftBeamformer.h"

my_Output PhaseshiftBeamformer(
		data_psb cov_Mat_re[NUMELEMENTS],
		data_psb cov_Mat_im[NUMELEMENTS],
		data_psb weightsRe[NUMELEMENTS],
		data_psb weightsIm[NUMELEMENTS])
{
 // 初始化
	data_psb rx_weightsRe[NUMELEMENTS];
	data_psb rx_weightsIm[NUMELEMENTS];

	data_psb sum_re = 0.0;
	data_psb sum_im = 0.0;;

	for (int j = 0; j < NUMELEMENTS; j++) {
		data_psb a = cov_Mat_re[j] * weightsRe[j];
		data_psb b = cov_Mat_im[j] * weightsIm[j];
		data_psb c = cov_Mat_re[j] * weightsIm[j];
		data_psb d = cov_Mat_im[j] * weightsRe[j];
		rx_weightsRe[j] = a - b; // 实部
		rx_weightsIm[j] = c + d; // 虚部

		//data_psb rx_weightsRe = cov_Mat_re[j] * weightsRe[j] - cov_Mat_im[j] * weightsIm[j]; // 实部
        //data_psb rx_weightsIm = cov_Mat_re[j] * weightsIm[j] + cov_Mat_im[j] * weightsRe[j]; // 虚部
        // j从1到7：天线入射的叠加

        }
	PhaseshiftBeamformer_label0:for (int j = 0; j < NUMELEMENTS; j++) {
		sum_re += rx_weightsRe[j];
		sum_im += rx_weightsIm[j];
	}

        my_Output result;
        result.a = sum_re / NUMELEMENTS;
        result.b = sum_im / NUMELEMENTS;
        return result;
    }

