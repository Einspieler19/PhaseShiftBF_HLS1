
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "PhaseshiftBeamformer.h"

using namespace std;



#define N 1001
#define NUM_ELEMENTS 7


////////////////////////// 确定对比精度 //////////////////////////
#ifdef data_tb_DATA
#define ABS_ERR_THRESH 0.1
#else
#define ABS_ERR_THRESH 0.001
#endif


////////////////////////// 是否输出对比数据 //////////////////////////
#define WINDOW_FN_DEBUG 0


void GenerateInput(data_tb* re, data_tb* im) {
    // Generate input signal
    data_tb t[N];
    for (int i = 0; i < N; i++) {
        t[i] = i;
    }

    std::ofstream outFile;
    // 打开文件
    outFile.open("Signal.txt");

    data_tb fsignal = 0.01;
    //data_tb re[N];
    //data_tb im[N];
    for (int i = 0; i < N; i++) {
        re[i] = sin(2 * M_PI * fmod(fsignal * t[i], 1.0));
        im[i] = 0;
        // 写入数据
        outFile << re[i] << endl;
    }

    // 关闭文件
    outFile.close();


}


// 添加噪声函数
void addNoise(data_tb (*cov_Mat_re)[NUM_ELEMENTS], data_tb (*cov_Mat_im)[NUM_ELEMENTS], data_tb variance) {
    // 设置随机数种子
    srand((unsigned)time(NULL));
    double randNum_re;
    double randNum_im;
    ofstream FILE;
    //Save the results to a file
    FILE.open ("Noised.dat");
    for (int i = 0; i < N; i++) {
    	for (int j = 0; j < NUM_ELEMENTS; j++) {
    	randNum_re = (double)rand() / RAND_MAX-0.5;
    	randNum_im = (double)rand() / RAND_MAX-0.5;
    	cov_Mat_re[i][j] += sqrt(variance) * randNum_re;
    	cov_Mat_im[i][j] += sqrt(variance) * randNum_im;
    	}
    	FILE << cov_Mat_re[i][1] << endl;

    	}
    FILE.close();
}

void collectWave(data_tb* rx_re, data_tb* rx_im, data_tb (*cov_Mat_re)[NUM_ELEMENTS], data_tb (*cov_Mat_im)[NUM_ELEMENTS], data_tb incidentAngle)
{

	// 位置初始化
	    data_tb Elementpos[NUM_ELEMENTS] = {-1.5, -1, -0.5, 0.0, 0.5, 1, 1.5};

	// 权重
	    data_tb v2_re[NUM_ELEMENTS], v2_im[NUM_ELEMENTS];

	    data_tb c = 3e8;
	    data_tb fc = 300e6;
	    // 打开文件
	    //outFile.open("Weight.txt");
	// Position
	    for (int i = 0; i < NUM_ELEMENTS; i++) {
	        v2_re[i] = hls::cos(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[i] * c / fc);
	        v2_im[i] = hls::sin(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[i] * c / fc);
	        // 写入数据
	        //outFile << weightsRe[i] << weightsIm[i] << "i" << i <<  endl;
	    }

	    for (int i = 0; i < N; i++) {
	        for (int j = 0; j < NUM_ELEMENTS; j++) {
	        	cov_Mat_re[i][j] =  rx_re[i] * v2_re[j] - rx_im[i] * v2_im[j]; // 实部
	        	cov_Mat_im[i][j] = rx_re[i] * v2_im[j] + rx_im[i] * v2_re[j]; // 虚部
	        	}
	    }
}



void computeWeights(
		data_tb steeringAngle,
		data_tb* weightsRe,
		data_tb* weightsIm)
{

    std::ofstream outFile;
    // 打开文件
    outFile.open("Weights.dat");

	// 位置初始化
	data_tb Elementpos[NUM_ELEMENTS] = {-1.5, -1, -0.5, 0.0, 0.5, 1, 1.5};
    // 权重

    data_tb c = 3e8;
    data_tb fc = 300e6;

    // w是复数，分实部和虚部计算

    for (int i = 0; i < NUM_ELEMENTS; i++) {
    	weightsRe[i] = hls::cos(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * c / fc);
        weightsIm[i] = hls::sin(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * c / fc);
        // 写入数据
        outFile << weightsRe[i] << weightsIm[i] << "i" << i <<  endl;
    }

    // 关闭文件
    outFile.close();

}

data_typeB datatypeConverter(data_typeA data_i){
	data_typeB data_o = data_typeB(data_i);
	return data_o;
}

void parametersConverter(
		data_tb (*cov_Mat_re)[NUM_ELEMENTS],
		data_tb (*cov_Mat_im)[NUM_ELEMENTS],
		//data_tb* steeringAngle,
		data_tb* weightsRe,
		data_tb* weightsIm,
		data_psb (*hw_cov_Mat_re)[NUM_ELEMENTS],
		data_psb (*hw_cov_Mat_im)[NUM_ELEMENTS],
		//data_psb* hw_steeringAngle,
		data_psb* hw_weightsRe,
		data_psb* hw_weightsIm,
		data_tb* y_re,
		data_tb* y_im,
		data_psb* hw_y_re,
		data_psb* hw_y_im){

	for (int i = 0; i < N; i++) {
		hw_y_re[i] = datatypeConverter(y_re[i]);
		hw_y_im[i] = datatypeConverter(y_im[i]);
	        for (int j = 0; j < NUM_ELEMENTS; j++) {
	        	hw_cov_Mat_re[i][j] = datatypeConverter(cov_Mat_re[i][j]);
	        	hw_cov_Mat_im[i][j] = datatypeConverter(cov_Mat_im[i][j]);
	        }
	}
	for (int j = 0; j < NUM_ELEMENTS; j++) {
	hw_weightsRe[j] = datatypeConverter(weightsRe[j]);
	hw_weightsIm[j] = datatypeConverter(weightsIm[j]);

	//cout << "hw_weightsRe[i]:   " << hw_weightsRe[i] <<endl;
	//cout << "weightsRe[i]:   " << weightsRe[i] <<endl;
	}
}





void SW_PhaseshiftBeamformer(
		data_tb (*cov_Mat_re)[NUM_ELEMENTS],
		data_tb (*cov_Mat_im)[NUM_ELEMENTS],
		data_tb steeringAngle,
		data_tb* weightsRe,
		data_tb* weightsIm,
		data_tb* y_re,
		data_tb* y_im)
{

    std::ofstream outFile;

    // 打开文件
    outFile.open("SWoutput.dat");
    // Perform beamforming


    for (int i = 0; i < N; i++) {
        y_re[i] = 0.0;
        y_im[i] = 0.0;
 // 初始化
        for (int j = 0; j < NUM_ELEMENTS; j++) {
            data_tb rx_weightsRe = cov_Mat_re[i][j] * weightsRe[j] - cov_Mat_im[i][j] * weightsIm[j]; // 实部
            data_tb rx_weightsIm = cov_Mat_re[i][j] * weightsIm[j] + cov_Mat_im[i][j] * weightsRe[j]; // 虚部
            //cout << rx_weightsRe << "rx_weightsRe" <<endl;
 // j从1到7：天线入射的叠加
            y_re[i] += rx_weightsRe;
            y_im[i] += rx_weightsIm;
        }
 // 除以7：功率归一

    y_re[i] /= NUM_ELEMENTS;
    y_im[i] /= NUM_ELEMENTS;

        // 写入数据
        outFile << y_re[i] << endl;
    }

    // 关闭文件
    outFile.close();
}




int main()
{

    // Initialize input
    data_tb rx_re[N], rx_im[N];
    data_psb hw_rx_re[N], hw_rx_im[N];


    data_tb cov_Mat_re[N][NUM_ELEMENTS], cov_Mat_im[N][NUM_ELEMENTS];
    data_psb hw_cov_Mat_re[N][NUM_ELEMENTS], hw_cov_Mat_im[N][NUM_ELEMENTS];


    // Set incident angle
    data_tb incidentAngle = M_PI/3;
    data_tb steeringAngle = -M_PI/3;
    data_psb hw_steeringAngle = -M_PI/3;
    data_tb variance = 1;
    data_tb weightsRe[NUM_ELEMENTS], weightsIm[NUM_ELEMENTS];
    data_psb hw_weightsRe[NUM_ELEMENTS], hw_weightsIm[NUM_ELEMENTS];

    // Initialize output
    data_psb hw_result_re[N], hw_result_im[N];
    data_tb sw_result_re[N], sw_result_im[N];


	GenerateInput(rx_re, rx_im);

	collectWave(rx_re, rx_im, cov_Mat_re, cov_Mat_im, incidentAngle);

	addNoise(cov_Mat_re, cov_Mat_im, variance);

	computeWeights(steeringAngle, weightsRe, weightsIm);

    // Call beamforming SW function
    SW_PhaseshiftBeamformer(cov_Mat_re, cov_Mat_im, steeringAngle, weightsRe, weightsIm, sw_result_re, sw_result_im);

    parametersConverter(
    		cov_Mat_re,
    		cov_Mat_im,
     		weightsRe,
    		weightsIm,
    		hw_cov_Mat_re,
    		hw_cov_Mat_im,
    		hw_weightsRe,
    		hw_weightsIm,
			rx_re,
			rx_im,
			hw_rx_re,
			hw_rx_im
			);


    //PhaseshiftBeamformer(hw_cov_Mat_re, hw_cov_Mat_im, hw_steeringAngle, hw_weightsRe, hw_weightsIm, hw_sw_result_re, hw_sw_result_im);

    // Call DUT
    PhaseshiftBeamformer(hw_cov_Mat_re, hw_cov_Mat_im, hw_steeringAngle, hw_weightsRe, hw_weightsIm, hw_result_re, hw_result_im);

////////////////////////// 结果对比 //////////////////////////
	unsigned err_cnt = 0;
	
	
//检查精度
   cout << "Checking results against a tolerance of " << ABS_ERR_THRESH << endl;
   cout << fixed << setprecision(5);

//循环求差
   for (unsigned i = 0; i < N; i++) {
      data_tb abs_err_re = (data_tb)hw_result_re[i] - (data_tb)sw_result_re[i];
      data_tb abs_err_im = (data_tb)hw_result_im[i] - (data_tb)sw_result_im[i];


////////////////////////// 输出每组差值 //////////////////////////
#if WINDOW_FN_DEBUG
      cout << "i = " << i << "\thw_result = " << hw_result_re[i] << " + j " << hw_result_im[i];
      cout << "\t sw_result = " << sw_result_re[i] << " + j " << sw_result_im[i] << endl;
#endif

////////////////////////// 超阈值报错 //////////////////////////

   if ((fabs(abs_err_re) > ABS_ERR_THRESH)|(fabs(abs_err_im) > ABS_ERR_THRESH)) {
         cout << "Error threshold exceeded: i = " << i;
         cout << "  Expected: "  << sw_result_re[i] << " + j " << sw_result_im[i];
         cout << "  Got: "  << hw_result_re[i] << " + j " << hw_result_im[i];
         cout << "  Delta: " << abs_err_re << " + j " << abs_err_im << endl;
         err_cnt++;
         cout << endl;
      }

   }

   if (err_cnt) {
      cout << "!!! TEST FAILED - " << err_cnt;
      cout << " results out of tolerance." << endl;
   } else
      cout << "Test Passed" << endl;

   // Only return 0 on success
   return err_cnt;
}

