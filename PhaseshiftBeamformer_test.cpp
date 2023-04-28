#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "PhaseshiftBeamformer.h"

using namespace std;

typedef matrix_t matrix_t;

#define N 1001
#define NUM_ELEMENTS 7


////////////////////////// 确定对比精度 //////////////////////////
#ifdef matrix_t_DATA
#define ABS_ERR_THRESH 0.1
#else
#define ABS_ERR_THRESH 0.001
#endif


////////////////////////// 是否输出对比数据 //////////////////////////
#define WINDOW_FN_DEBUG 1


void GenerateInput(matrix_t* re, matrix_t* im) {
    // Generate input signal
    matrix_t t[N];
    for (int i = 0; i < N; i++) {
        t[i] = i;
    }

    std::ofstream outFile;
    // 打开文件
    outFile.open("Signal.txt");

    matrix_t fsignal = 0.01;
    //matrix_t re[N];
    //matrix_t im[N];
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
void addNoise(matrix_t (*cov_Mat_re)[NUM_ELEMENTS], matrix_t (*cov_Mat_im)[NUM_ELEMENTS], matrix_t variance) {
    // 设置随机数种子
    srand((unsigned)time(NULL));
    double randNum_re;
    double randNum_im;
    ofstream FILE;
    //Save the results to a file
    FILE.open ("result_noised.dat");
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

void Collectwave(matrix_t* rx_re, matrix_t* rx_im, matrix_t (*cov_Mat_re)[NUM_ELEMENTS], matrix_t (*cov_Mat_im)[NUM_ELEMENTS], matrix_t incidentAngle)
{

	// 位置初始化
	    matrix_t Elementpos[NUM_ELEMENTS] = {-1.5, -1, -0.5, 0.0, 0.5, 1, 1.5};

	// 权重
	    matrix_t v2_re[NUM_ELEMENTS], v2_im[NUM_ELEMENTS];

	    matrix_t c = 3e8;
	    matrix_t fc = 300e6;
	    // 打开文件
	    //outFile.open("Weight.txt");
	// Position
	    for (int i = 0; i < NUM_ELEMENTS; i++) {
	        v2_re[i] = hls::cos(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[i] * c / fc);
	        v2_im[i] = hls::sin(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[i] * c / fc);
	        // 写入数据
	        //outFile << w2_re[i] << w2_im[i] << "i" << i <<  endl;
	    }

	    for (int i = 0; i < N; i++) {
	        for (int j = 0; j < NUM_ELEMENTS; j++) {
	        	cov_Mat_re[i][j] =  rx_re[i] * v2_re[j] - rx_im[i] * v2_im[j]; // 实部
	        	cov_Mat_im[i][j] = rx_re[i] * v2_im[j] + rx_im[i] * v2_re[j]; // 虚部
	        	}
	    }
}



void SW_PhaseshiftBeamformer(
		matrix_t (*cov_Mat_re)[NUM_ELEMENTS],
		matrix_t (*cov_Mat_im)[NUM_ELEMENTS],
		matrix_t steeringAngle,
		matrix_t* y_re,
		matrix_t* y_im)
{

    std::ofstream outFile;
    // 打开文件
    outFile.open("Weights.dat");

	// 位置初始化
	matrix_t Elementpos[NUM_ELEMENTS] = {-1.5, -1, -0.5, 0.0, 0.5, 1, 1.5};
    // 权重
    matrix_t w2_re[NUM_ELEMENTS], w2_im[NUM_ELEMENTS];
    matrix_t c = 3e8;
    matrix_t fc = 300e6;


    // w是复数，分实部和虚部计算

    for (int i = 0; i < NUM_ELEMENTS; i++) {
        w2_re[i] = hls::cos(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * c / fc);
        w2_im[i] = hls::sin(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * c / fc);
        // 写入数据
        outFile << w2_re[i] << w2_im[i] << "i" << i <<  endl;
    }

    // 关闭文件
    outFile.close();

    // 打开文件
    outFile.open("output.dat");
    // Perform beamforming


    for (int i = 0; i < N; i++) {
        y_re[i] = 0.0;
        y_im[i] = 0.0;
 // 初始化
        for (int j = 0; j < NUM_ELEMENTS; j++) {
            matrix_t rx_w2_re = cov_Mat_re[i][j] * w2_re[j] - cov_Mat_im[i][j] * w2_im[j]; // 实部
            matrix_t rx_w2_im = cov_Mat_re[i][j] * w2_im[j] + cov_Mat_im[i][j] * w2_re[j]; // 虚部
            //cout << rx_w2_re << "rx_w2_re" <<endl;
 // j从1到7：天线入射的叠加
            y_re[i] += rx_w2_re;
            y_im[i] += rx_w2_im;
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
    matrix_t rx_re[N], rx_im[N];
    //for (int i = 0; i < N; i++) {
    //    matrix_t rand_val = (matrix_t)rand() / RAND_MAX;
    //    rx_re[i] = sin(2 * M_PI * 0.01 * i) + 0.1 * rand_val;
    //   rx_im[i] = 0.1 * rand_val;
    //}
    matrix_t cov_Mat_re[N][NUM_ELEMENTS], cov_Mat_im[N][NUM_ELEMENTS];


    // Set incident angle
    matrix_t incidentAngle = M_PI/3;
    matrix_t steeringAngle = -M_PI/3;
    matrix_t variance = 1;

    // Initialize output
    matrix_t hw_result_re[N], hw_result_im[N];
    matrix_t sw_result_re[N], sw_result_im[N];

	GenerateInput(rx_re, rx_im);

	Collectwave(rx_re, rx_im, cov_Mat_re, cov_Mat_im, incidentAngle);

	addNoise( cov_Mat_re, cov_Mat_im, variance);
    // Call beamforming SW function
    SW_PhaseshiftBeamformer(cov_Mat_re, cov_Mat_im, steeringAngle, sw_result_re, sw_result_im);

    // Call DUT
    PhaseshiftBeamformer(rx_re, rx_im, incidentAngle, hw_result_re, hw_result_im);

////////////////////////// 结果对比 //////////////////////////
	unsigned err_cnt = 0;
	
	
//检查精度
   cout << "Checking results against a tolerance of " << ABS_ERR_THRESH << endl;
   cout << fixed << setprecision(5);

//循环求差
   for (unsigned i = 0; i < N; i++) {
      matrix_t abs_err_re = hw_result_re[i] - sw_result_re[i];
      matrix_t abs_err_im = hw_result_im[i] - sw_result_im[i];


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
      }
   cout << endl;
   }

   if (err_cnt) {
      cout << "!!! TEST FAILED - " << err_cnt;
      cout << " results out of tolerance." << endl;
   } else
      cout << "Test Passed" << endl;

   // Only return 0 on success
   return err_cnt;
}

