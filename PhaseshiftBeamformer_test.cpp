#include "PhaseshiftBeamformer.h"

using namespace std;


////////////////////////// 确定对比精度 //////////////////////////
#ifdef data_tb_DATA
#define ABS_ERR_THRESH 0.1
#else
#define ABS_ERR_THRESH 0.001
#endif



void GenerateInput(data_tb* re, data_tb* im) {
    // Generate input signal
    data_tb t[SIGNALLENGTH];
    for (int i = 0; i < SIGNALLENGTH; i++) {
        t[i] = i;
    }

    std::ofstream outFile;
    // 打开文件
    outFile.open("Signal.txt");

    data_tb fsignal = 0.01;


    for (int i = 0; i < SIGNALLENGTH; i++) {
        re[i] = sin(2 * M_PI * fmod(fsignal * t[i], 1.0));
        im[i] = 0;
        // 写入数据
        outFile << re[i] << endl;
    }

    // 关闭文件
    outFile.close();


}


// 添加噪声函数
void addNoise(data_tb (*cov_Mat_re)[NUMELEMENTS], data_tb (*cov_Mat_im)[NUMELEMENTS], data_tb variance) {
    // 设置随机数种子
    srand((unsigned)time(NULL));
    double randNum_re;
    double randNum_im;
    ofstream FILE;
    //Save the results to a file
    FILE.open ("Noised.dat");
    for (int i = 0; i < SIGNALLENGTH; i++) {
    	for (int j = 0; j < NUMELEMENTS; j++) {
    	randNum_re = (double)rand() / RAND_MAX-0.5;
    	randNum_im = (double)rand() / RAND_MAX-0.5;
    	cov_Mat_re[i][j] += sqrt(variance) * randNum_re;
    	cov_Mat_im[i][j] += sqrt(variance) * randNum_im;
    	}
    	FILE << cov_Mat_re[i][1] << endl;

    	}
    FILE.close();
}

void collectWave(data_tb* rx_re, data_tb* rx_im, data_tb (*cov_Mat_re)[NUMELEMENTS], data_tb (*cov_Mat_im)[NUMELEMENTS], data_tb incidentAngle)
{

	// 位置初始化
	    data_tb Elementpos[NUMELEMENTS];
	    generateElementpos<data_tb>(Elementpos);

	// 权重
	    data_tb v2_re[NUMELEMENTS], v2_im[NUMELEMENTS];

	    // 打开文件
	    //outFile.open("Weight.txt");
	// Position
	    for (int i = 0; i < NUMELEMENTS; i++) {
	        v2_re[i] = hls::cos(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
	        v2_im[i] = hls::sin(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
	        // 写入数据
	        //outFile << weightsRe[i] << weightsIm[i] << "i" << i <<  endl;
	    }

	    for (int i = 0; i < SIGNALLENGTH; i++) {
	        for (int j = 0; j < NUMELEMENTS; j++) {
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
    data_tb Elementpos[NUMELEMENTS];
    generateElementpos<data_tb>(Elementpos);

    // 权重

    // w是复数，分实部和虚部计算

    for (int i = 0; i < NUMELEMENTS; i++) {
    	weightsRe[i] = hls::cos(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
        weightsIm[i] = hls::sin(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
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
		data_tb (*cov_Mat_re)[NUMELEMENTS],
		data_tb (*cov_Mat_im)[NUMELEMENTS],
		//data_tb* steeringAngle,
		data_tb* weightsRe,
		data_tb* weightsIm,
		data_psb (*hw_cov_Mat_re)[NUMELEMENTS],
		data_psb (*hw_cov_Mat_im)[NUMELEMENTS],
		//data_psb* hw_steeringAngle,
		data_psb* hw_weightsRe,
		data_psb* hw_weightsIm,
		data_tb* y_re,
		data_tb* y_im,
		data_psb* hw_y_re,
		data_psb* hw_y_im){

	for (int i = 0; i < SIGNALLENGTH; i++) {
		hw_y_re[i] = datatypeConverter(y_re[i]);
		hw_y_im[i] = datatypeConverter(y_im[i]);
	        for (int j = 0; j < NUMELEMENTS; j++) {
	        	hw_cov_Mat_re[i][j] = datatypeConverter(cov_Mat_re[i][j]);
	        	hw_cov_Mat_im[i][j] = datatypeConverter(cov_Mat_im[i][j]);
	        }
	}
	for (int j = 0; j < NUMELEMENTS; j++) {
	hw_weightsRe[j] = datatypeConverter(weightsRe[j]);
	hw_weightsIm[j] = datatypeConverter(weightsIm[j]);

	//cout << "hw_weightsRe[i]:   " << hw_weightsRe[i] <<endl;
	//cout << "weightsRe[i]:   " << weightsRe[i] <<endl;
	}
}





void SW_PhaseshiftBeamformer(
		data_tb (*cov_Mat_re)[NUMELEMENTS],
		data_tb (*cov_Mat_im)[NUMELEMENTS],
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


    for (int i = 0; i < SIGNALLENGTH; i++) {
        y_re[i] = 0.0;
        y_im[i] = 0.0;
 // 初始化
        for (int j = 0; j < NUMELEMENTS; j++) {
            data_tb rx_weightsRe = cov_Mat_re[i][j] * weightsRe[j] - cov_Mat_im[i][j] * weightsIm[j]; // 实部
            data_tb rx_weightsIm = cov_Mat_re[i][j] * weightsIm[j] + cov_Mat_im[i][j] * weightsRe[j]; // 虚部
            //cout << rx_weightsRe << "rx_weightsRe" <<endl;
 // j从1到7：天线入射的叠加
            y_re[i] += rx_weightsRe;
            y_im[i] += rx_weightsIm;
        }
 // 除以7：功率归一

    y_re[i] /= NUMELEMENTS;
    y_im[i] /= NUMELEMENTS;

        // 写入数据
        outFile << y_re[i] << endl;
    }
    // 关闭文件
    outFile.close();
}



int main()
{

    // Initialize input
    data_tb rx_re[SIGNALLENGTH], rx_im[SIGNALLENGTH];
    data_psb hw_rx_re[SIGNALLENGTH], hw_rx_im[SIGNALLENGTH];


    data_tb cov_Mat_re[SIGNALLENGTH][NUMELEMENTS], cov_Mat_im[SIGNALLENGTH][NUMELEMENTS];
    data_psb hw_cov_Mat_re[SIGNALLENGTH][NUMELEMENTS], hw_cov_Mat_im[SIGNALLENGTH][NUMELEMENTS];


    // Set incident angle
    data_tb incidentAngle = DOANGLE;
    data_tb steeringAngle = STEERINGANGLE;
    data_psb hw_steeringAngle = STEERINGANGLE;
    data_tb variance = NOISEVARIANCE;
    data_tb weightsRe[NUMELEMENTS], weightsIm[NUMELEMENTS];
    data_psb hw_weightsRe[NUMELEMENTS], hw_weightsIm[NUMELEMENTS];

    // Initialize output
    data_psb hw_result_re[SIGNALLENGTH], hw_result_im[SIGNALLENGTH];
    data_tb sw_result_re[SIGNALLENGTH], sw_result_im[SIGNALLENGTH];


	GenerateInput(rx_re, rx_im);

	collectWave(rx_re, rx_im, cov_Mat_re, cov_Mat_im, incidentAngle);

	//addNoise(cov_Mat_re, cov_Mat_im, variance);

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

    data_psb in_cov_Mat_re[NUMELEMENTS];
    data_psb in_cov_Mat_im[NUMELEMENTS];
    data_psb y_re[SIGNALLENGTH];
    data_psb y_im[SIGNALLENGTH];

    std::ofstream outFile;


    // 打开文件
    outFile.open("HWoutput.dat");
    for (int i = 0; i < SIGNALLENGTH; i++) {
        for (int j = 0; j < NUMELEMENTS; j++) {
        	in_cov_Mat_re[j] = hw_cov_Mat_re[i][j];
        	in_cov_Mat_im[j] = hw_cov_Mat_im[i][j];
        }

        // Perform beamforming

        my_Output hw_myResult = PhaseshiftBeamformer(in_cov_Mat_re, in_cov_Mat_im, hw_weightsRe, hw_weightsIm);
        //cout<< hw_myResult.a <<endl;
        hw_result_re[i] = hw_myResult.a;
        hw_result_im[i] = hw_myResult.b;
        // 写入数据
        outFile << hw_result_re[i] << endl;
    }
    // 关闭文件
    outFile.close();

////////////////////////// 结果对比 //////////////////////////
	unsigned err_cnt = 0;


//检查精度
   cout << "Checking results against a tolerance of " << ABS_ERR_THRESH << endl;
   cout << fixed << setprecision(5);
   // 打开文件
       outFile.open("SWoutput2.dat");
       for (int i = 0; i < SIGNALLENGTH; i++) {
       	// 写入数据
       	outFile << sw_result_re[i] << endl;
       }
       // 关闭文件
       outFile.close();
//循环求差
       // 打开文件
           outFile.open("SWoutput3.dat");
   for (unsigned i = 0; i < SIGNALLENGTH; i++) {
      data_tb abs_err_re = (data_tb)hw_result_re[i] - sw_result_re[i];
      data_tb abs_err_im = (data_tb)hw_result_im[i] - sw_result_im[i];

      //outFile << "\t sw_result = " << sw_result_re[i] << " + j " << sw_result_im[i] << endl;

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

   }       // 关闭文件
   outFile.close();

   if (err_cnt) {
      cout << "!!! TEST FAILED - " << err_cnt;
      cout << " results out of tolerance." << endl;
   } else
      cout << "Test Passed" << endl;

   // Only return 0 on success

   return 0;
}

