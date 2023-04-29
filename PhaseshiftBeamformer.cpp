#include "PhaseshiftBeamformer.h"


void PhaseshiftBeamformer(
		data_psb (*cov_Mat_re)[NUMELEMENTS],
		data_psb (*cov_Mat_im)[NUMELEMENTS],
		data_psb steeringAngle,
		data_psb* weightsRe,
		data_psb* weightsIm,
		data_psb* y_re,
		data_psb* y_im)
{

    //std::ofstream outFile;

    // 打开文件
    //outFile.open("HWoutput.dat");
    // Perform beamforming


    for (int i = 0; i < SIGNALLENGTH; i++) {
        y_re[i] = 0.0;
        y_im[i] = 0.0;
 // 初始化
        for (int j = 0; j < NUMELEMENTS; j++) {
            data_psb rx_weightsRe = cov_Mat_re[i][j] * weightsRe[j] - cov_Mat_im[i][j] * weightsIm[j]; // 实部
            data_psb rx_weightsIm = cov_Mat_re[i][j] * weightsIm[j] + cov_Mat_im[i][j] * weightsRe[j]; // 虚部
            //cout << "rx_weightsRe:   " << rx_weightsRe[j] <<endl;
            //cout << cov_Mat_re[i][j] << "cov_Mat_re" <<endl;
 // j从1到7：天线入射的叠加
            y_re[i] += rx_weightsRe;
            y_im[i] += rx_weightsIm;
        }
 // 除以7：功率归一

    y_re[i] /= NUMELEMENTS;
    y_im[i] /= NUMELEMENTS;

        // 写入数据
        //outFile << y_re[i] << endl;
    }

    // 关闭文件
    //outFile.close();
}
