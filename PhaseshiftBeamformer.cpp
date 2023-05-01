#include "PhaseshiftBeamformer.h"


void PhaseshiftBeamformer(){
	data_psb cov_Mat_re[SIGNALLENGTH][NUMELEMENTS], cov_Mat_im[SIGNALLENGTH][NUMELEMENTS];
	data_psb weightsRe[NUMELEMENTS],weightsIm[NUMELEMENTS];
	data_psb y_re[SIGNALLENGTH],y_im[SIGNALLENGTH];

	std::ifstream inFile;
	// 打开文件
	inFile.open("Noised.dat");
	for (int i = 0; i < SIGNALLENGTH; i++) {
	            for (int j = 0; j < NUMELEMENTS; j++) {
	            	inFile >> cov_Mat_re[i][j];
	            	inFile >> cov_Mat_im[i][j];
	            }
	    }
	inFile.close();


	// 打开文件
	inFile.open("Weights.dat");
	for (int j = 0; j < NUMELEMENTS; j++) {
		inFile >> weightsRe[j];
		inFile >> weightsIm[j];
	    }
	inFile.close();

    std::ofstream outFile;
    // 打开文件
    outFile.open("HWout.dat");
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
    outFile << y_re[i] << endl;
    outFile << y_im[i] << endl;
    }

    // 关闭文件
    outFile.close();
}
