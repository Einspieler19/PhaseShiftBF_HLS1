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
my_complex_Array<data_tb, NUMELEMENTS> addNoise(
		data_tb collected_re[NUMELEMENTS],
		data_tb collected_im[NUMELEMENTS],
		data_tb variance)
{
	my_complex_Array<data_tb, NUMELEMENTS> result;
	// 定义随机数引擎和分布
	std::random_device rd;  // 用于生成随机种子
	std::mt19937 gen(rd()); // 随机数引擎
	std::uniform_real_distribution<double> dis(0.0, 1.0); // 生成 [0.0, 1.0) 范围内的随机浮点数

    ofstream FILE;
    //Save the results to a file
    FILE.open ("Noised.dat");
    for (int j = 0; j < NUMELEMENTS; j++) {
    	data_tb randNum_re = (data_tb)dis(gen);
    	data_tb randNum_im = (data_tb)dis(gen);

    	result.re[j] = collected_re[j] + sqrt(variance) * randNum_re;
    	result.im[j] = collected_im[j] + sqrt(variance) * randNum_im;

    	FILE<<result.im[j]<<endl;

    }
    FILE.close();
    return result;
}

my_complex_Array<data_tb, NUMELEMENTS> collectWave(
		data_tb rx_re,
		data_tb rx_im,
		data_tb cov_Mat_re[NUMELEMENTS],
		data_tb cov_Mat_im[NUMELEMENTS],
		data_tb incidentAngle)
{

	// 位置初始化
	    data_tb Elementpos[NUMELEMENTS];
	    generateElementpos<data_tb>(Elementpos);

	// 权重
	    data_tb v2_re[NUMELEMENTS], v2_im[NUMELEMENTS];
	    my_complex_Array<data_tb, NUMELEMENTS> result;
	    for (int j = 0; j < NUMELEMENTS; j++) {
	    	v2_re[j] = hls::cos(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[j] * SPEEDOFLIGHT / CENTERFREQ);
	    	v2_im[j] = hls::sin(2 * M_PI * hls::sinf(incidentAngle) * Elementpos[j] * SPEEDOFLIGHT / CENTERFREQ);

	    	cov_Mat_re[j] = rx_re * v2_re[j] - rx_im * v2_im[j]; // 实部
	    	cov_Mat_im[j] = rx_re * v2_im[j] + rx_im * v2_re[j]; // 虚部

	    	result.re[j] = cov_Mat_re[j];
	    	result.im[j] = cov_Mat_im[j];
	    	        	}
	    return result;
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

    // w是复数，分实部和虚部计算

    for (int i = 0; i < NUMELEMENTS; i++) {
    	weightsRe[i] = hls::cos(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
        weightsIm[i] = hls::sin(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
        // 写入数据
        outFile << weightsRe[i] << "+"<< weightsIm[i] << "i" <<  endl;
    }

    // 关闭文件
    outFile.close();

}


my_complex_Value<data_tb> SW_PhaseshiftBeamformer(
		data_tb cov_Mat_re[NUMELEMENTS],
		data_tb cov_Mat_im[NUMELEMENTS],
		data_tb weightsRe[NUMELEMENTS],
		data_tb weightsIm[NUMELEMENTS])
{
	data_tb sum_re = 0.0;
	data_tb sum_im = 0.0;
	data_tb rx_weightsRe = 0.0;
	data_tb rx_weightsIm = 0.0;
    // Perform beamforming
        for (int j = 0; j < NUMELEMENTS; j++) {
            rx_weightsRe = cov_Mat_re[j] * weightsRe[j] - cov_Mat_im[j] * weightsIm[j]; // 实部
            rx_weightsIm = cov_Mat_re[j] * weightsIm[j] + cov_Mat_im[j] * weightsRe[j]; // 虚部
            // j从1到7：天线入射的叠加
            sum_re += rx_weightsRe;
            sum_im += rx_weightsIm;
            }
 // 除以7：功率归一
    sum_re /= NUMELEMENTS;
    sum_im /= NUMELEMENTS;

    my_complex_Value<data_tb> result;
    result.re = sum_re;
    result.im = sum_im;
    return result;
}

int main()
{
    // Set Consts
    data_tb incidentAngle = DOANGLE;
    data_tb steeringAngle = STEERINGANGLE;
    data_tb variance = NOISEVARIANCE;

    // Initialize input
    data_tb rx_re[SIGNALLENGTH], rx_im[SIGNALLENGTH];
    data_tb all_cov_Mat_re[SIGNALLENGTH][NUMELEMENTS], all_cov_Mat_im[SIGNALLENGTH][NUMELEMENTS];

    data_tb weightsRe[NUMELEMENTS], weightsIm[NUMELEMENTS];
    data_psb hw_weightsRe[NUMELEMENTS], hw_weightsIm[NUMELEMENTS];

    // Initialize output
    data_tb sw_result_re[SIGNALLENGTH], sw_result_im[SIGNALLENGTH];
    data_psb hw_result_re[SIGNALLENGTH], hw_result_im[SIGNALLENGTH];

    std::ofstream outFile1;
    std::ofstream outFile2;

	GenerateInput(rx_re, rx_im);

	computeWeights(steeringAngle, weightsRe, weightsIm);

    // Call beamforming SW function

	my_complex_Array<data_psb, NUMELEMENTS> sw_covMats;
	my_complex_Array<data_psb, NUMELEMENTS> hw_collectedWave;
	my_complex_Array<data_psb, NUMELEMENTS> sw_collectedWave;
	my_complex_Array<data_psb, NUMELEMENTS> hw_Weights;

	my_complex_Value<data_psb> sw_myResult;
	my_complex_Value<data_psb> hw_myResult;

    // 打开文件
    outFile1.open("SWoutput.dat");
    outFile2.open("HWoutput.dat");

    for (int i = 0; i < SIGNALLENGTH; i++) {
    	// assign
    	sw_covMats = complexarrayConverter<data_tb, data_tb, NUMELEMENTS>(all_cov_Mat_re[i], all_cov_Mat_im[i]);

    	sw_collectedWave = collectWave(rx_re[i], rx_im[i], sw_covMats.re, sw_covMats.im, incidentAngle);
    	//sw_collectedWave = addNoise(sw_collectedWave.re,sw_collectedWave.im,variance);

        hw_collectedWave = complexarrayConverter<data_tb, data_psb, NUMELEMENTS>(sw_collectedWave.re, sw_collectedWave.im);

        hw_Weights = complexarrayConverter<data_tb, data_psb, NUMELEMENTS>(weightsRe, weightsIm);

        sw_myResult = SW_PhaseshiftBeamformer(sw_collectedWave.re, sw_collectedWave.im, weightsRe, weightsIm);
        hw_myResult = PhaseshiftBeamformer(hw_collectedWave.re, hw_collectedWave.im, hw_Weights.re, hw_Weights.im);


        sw_result_re[i] = sw_myResult.re;
        sw_result_im[i] = sw_myResult.im;

        hw_result_re[i] = hw_myResult.re;
        hw_result_im[i] = hw_myResult.im;

        // 写入数据
        outFile1 << sw_result_re[i] << endl;
        outFile2 << hw_result_re[i] << endl;
    }
    // 关闭文件
    outFile1.close();
    outFile2.close();


////////////////////////// 结果对比 //////////////////////////
	unsigned err_cnt = 0;


//检查精度
   cout << "Checking results against a tolerance of " << ABS_ERR_THRESH << endl;
   cout << fixed << setprecision(5);
   //循环求差
   for (unsigned i = 0; i < SIGNALLENGTH; i++) {
      data_tb abs_err_re = (data_tb)hw_result_re[i] - sw_result_re[i];
      data_tb abs_err_im = (data_tb)hw_result_im[i] - sw_result_im[i];

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

   return 0;
}

