
#include "PhaseshiftBeamformer.h"


using namespace std;

/////////////////////////////////// 检测部分 ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

////////////////////////// 生成信号 //////////////////////////
void generateInput(data_tb* re, data_tb* im) {
    // 时间信号
    data_tb t[SIGNALLENGTH];
    for (int i = 0; i < SIGNALLENGTH; i++) {
        t[i] = i;
    }

    // 信号写进文件
    std::ofstream outFile;
    outFile.open("Signal.txt");
    data_tb fsignal = MYSAMPLERATE_PERSAMPLE;

	// 跟着t[i]逐点生成信号RE, 初始IM全为0
    for (int i = 0; i < SIGNALLENGTH; i++) {
        re[i] = sin(2 * M_PI * fmod(fsignal * t[i], 1.0));
        im[i] = 0.0;

        outFile << re[i] << endl;     //写下
    }
    outFile.close();    // 关闭文件
}


////////////////////////// 将噪声表读到数组中 //////////////////////////
// 实部噪声矩阵 Noise_re[1001][7]
// 虚部噪声矩阵 Noise_im[1001][7]

void getNoisedata(data_tb Noise_re[SIGNALLENGTH][NUMELEMENTS],data_tb Noise_im[SIGNALLENGTH][NUMELEMENTS])
{
    std::ifstream noiseFile;
    noiseFile.open("/home/dian/ModelComposer/MyProject/HLS_Trials/PhaseshiftBeamformer/PhaseshiftBeamformer_prj/Noise.dat"); //按绝对路径读取噪声表

	// 循环读取噪声到两个矩阵
    for (int i = 0; i < SIGNALLENGTH; i++) {
    	for (int j = 0; j < NUMELEMENTS; j++) {
    		noiseFile >> Noise_re[i][j];
    		noiseFile >> Noise_im[i][j];
    	}
    }
    noiseFile.close();
}


////////////////////// 使用噪声矩阵为相关阵添加噪声 //////////////////////////

my_complex_Array<data_tb, NUMELEMENTS> addNoise( // 返回加噪的相关阵 实部 + 虚部
		data_tb collected_re[NUMELEMENTS], // 入射相关阵 实部
		data_tb collected_im[NUMELEMENTS], // 入射相关阵 虚部
		data_tb noise_re[NUMELEMENTS], // 添加噪声 实部
		data_tb noise_im[NUMELEMENTS], // 添加噪声 虚部
		data_tb variance)
{
	my_complex_Array<data_tb, NUMELEMENTS> result; //加噪结果
	//循环加噪
    for (int j = 0; j < NUMELEMENTS; j++) {
    	result.re[j] = collected_re[j] + sqrt(variance) * noise_re[j] ;
    	result.im[j] = collected_im[j] + sqrt(variance) * noise_im[j] ;
    }
    return result;
}

////////////////////// 入射角 * 天线阵 = 阵列流形v //////////////////////////
////////////////////// 信号 * 阵列流形v = 相关阵r //////////////////////////
//每snapshot
my_complex_Array<data_tb, NUMELEMENTS> collectWave(
		data_tb rx_re,
		data_tb rx_im,
		data_tb manifoldRe[NUMELEMENTS],
		data_tb manifoldIm[NUMELEMENTS])
{
		my_complex_Array<data_tb, NUMELEMENTS> result;

	for (int j = 0; j < NUMELEMENTS; j++) {
			// 求相关阵
	    	result.re[j] = rx_re * manifoldRe[j] - rx_im * manifoldIm[j]; // 实部
	    	result.im[j] = rx_re * manifoldIm[j] + rx_im * manifoldRe[j]; // 虚部
	    }
	    return result;
}


/////////////////////// 定向角 * 天线阵 = 权重向量w ////////////////////////
void computeWeights(
		data_tb steeringAngle,
		data_tb* weightsRe,
		data_tb* weightsIm)
{
	// 初始化天线阵向量
    data_tb Elementpos[NUMELEMENTS];
    generateElementpos<data_tb>(Elementpos);

	// 求权重向量
    for (int i = 0; i < NUMELEMENTS; i++) {
    	weightsRe[i] = hls::cos(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
        weightsIm[i] = hls::sin(2 * M_PI * hls::sinf(steeringAngle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
    }
}

/////////////////////// SW 赋型 ////////////////////////

void SW_PhaseshiftBeamformer(
		data_tb cov_Mat_re[NUMELEMENTS], // 入射相关阵 实部
		data_tb cov_Mat_im[NUMELEMENTS], // 入射相关阵 虚部
		data_tb weightsRe[NUMELEMENTS], // 权重 实部
		data_tb weightsIm[NUMELEMENTS], // 权重 虚部
		data_tb *outRe, // 赋型求和 实部
		data_tb *outIm) // 赋型求和 虚部
{
	// 初始化
	data_tb sum_re = 0.0;
	data_tb sum_im = 0.0;
	data_tb rx_weightsRe = 0.0;
	data_tb rx_weightsIm = 0.0;
    // Perform beamforming
    for (int j = 0; j < NUMELEMENTS; j++) {
	    // 收集信号*权重 (s[j]*v[j])*w[j]
        rx_weightsRe = cov_Mat_re[j] * weightsRe[j] - cov_Mat_im[j] * weightsIm[j]; // 求实部
        rx_weightsIm = cov_Mat_re[j] * weightsIm[j] + cov_Mat_im[j] * weightsRe[j]; // 求虚部

        // 所有天线叠加(j从1到7)
        sum_re += rx_weightsRe;
        sum_im += rx_weightsIm;
    }
    // 功率归一(除以天线数量)
    sum_re /= NUMELEMENTS;
    sum_im /= NUMELEMENTS;

	// 第i个赋型信号
    *outRe = sum_re;
    *outIm = sum_im;
}


/////////////////////////////////// 主函数 ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


int main()
{
    // Set Consts
    data_tb incidentAngle = DOANGLE;
    data_tb steeringAngle = STEERINGANGLE;
    data_tb variance = NOISEVARIANCE;

    // 输入信号
    data_tb rx_re[SIGNALLENGTH], rx_im[SIGNALLENGTH];

	// 收集信号：全长入射阵
    data_tb all_cov_Mat_re[SIGNALLENGTH][NUMELEMENTS], all_cov_Mat_im[SIGNALLENGTH][NUMELEMENTS];

    data_tb manifoldRe[NUMELEMENTS], manifoldIm[NUMELEMENTS];

	// 权重
    data_tb weightsRe[NUMELEMENTS], weightsIm[NUMELEMENTS];
    data_psb hw_weightsRe[NUMELEMENTS], hw_weightsIm[NUMELEMENTS];

    // 赋型结果
    data_tb sw_result_re[SIGNALLENGTH], sw_result_im[SIGNALLENGTH];
    data_psb hw_result_re[SIGNALLENGTH], hw_result_im[SIGNALLENGTH];

	// 读文件
    std::ofstream outFile1;
    std::ofstream outFile2;

/////////////////////// 中间变量 ////////////////////////
	// snapshot入射阵
	my_complex_Array<data_tb, NUMELEMENTS> sw_collectedWave;
    my_complex_Array<data_psb, NUMELEMENTS> hw_collectedWave;

    // 权重
    my_complex_Array<data_psb, NUMELEMENTS> hw_Weights;

	// 赋型结果
    my_complex_Value<data_tb> sw_myResult;
    my_complex_Value<data_psb> hw_myResult;

	// 噪声
    data_tb Noise_re[SIGNALLENGTH][NUMELEMENTS];
    data_tb Noise_im[SIGNALLENGTH][NUMELEMENTS];
    my_complex_Array<data_tb, NUMELEMENTS> array_Noise;

    for(int i = 0; i < SIGNALLENGTH;i++){
	    // 入射信号
    	rx_re[i]=0.0;
    	rx_im[i]=0.0;

    	// 赋型结果 sw&hw
    	sw_result_re[i]=0.0;
    	sw_result_im[i]=0.0;
    	hw_result_re[i]=0.0;
    	hw_result_im[i]=0.0;

    	for(int j = 0; j < NUMELEMENTS;j++){
    	//噪声
    				Noise_re[i][j]=0.0;
    				Noise_im[i][j]=0.0;
		// 全长入射阵
    				all_cov_Mat_re[i][j]=0.0;
    				all_cov_Mat_im[i][j]=0.0;
    				}
    		}
    		for(int j = 0; j < NUMELEMENTS;j++){
		// 权重 sw&hw
    			weightsRe[j]=0.0;
    			weightsIm[j]=0.0;
    			manifoldRe[j]=0.0;
    			manifoldIm[j]=0.0;
    			hw_weightsRe[j]=0.0;
    			hw_weightsIm[j]=0.0;
    		}

	///////////////////////开始运行////////////////////////

	// 生成入射信号
	generateInput(rx_re, rx_im);


	computeManifoldvector<data_tb>(incidentAngle, manifoldRe, manifoldIm);

	// 计算权重
	computeManifoldvector<data_tb>(steeringAngle, weightsRe, weightsIm);


	// 读取噪声
	getNoisedata(Noise_re,Noise_im);

    // 打开文件
    outFile1.open("SWoutput.dat");
    outFile2.open("HWoutput.dat");

    for (int i = 0; i < SIGNALLENGTH; i++) {
	/////////////////////// 按snapshot递值 ////////////////////////

	// 使用已有入射信号流
	// 收集snapshot入射信号
	sw_collectedWave = collectWave(rx_re[i], rx_im[i], manifoldRe, manifoldIm);

	// 将噪声矩阵 递为 每snapshot的天线阵噪声分布
	array_Noise = complexarrayConverter<data_tb, data_tb, NUMELEMENTS>(Noise_re[i], Noise_im[i]);

	// 为收集信号加噪
    sw_collectedWave = addNoise(sw_collectedWave.re,sw_collectedWave.im, array_Noise.re, array_Noise.im, variance);

	// 将加噪的收集信号转为hw
    hw_collectedWave = complexarrayConverter<data_tb, data_psb, NUMELEMENTS>(sw_collectedWave.re, sw_collectedWave.im);

	// 将权重转为hw
    hw_Weights = complexarrayConverter<data_tb, data_psb, NUMELEMENTS>(weightsRe, weightsIm);

	// 信号* 阵列流形 * 权重 = 赋型结果
    SW_PhaseshiftBeamformer(sw_collectedWave.re, sw_collectedWave.im, weightsRe, weightsIm, sw_result_re + i, sw_result_im + i);
    PhaseshiftBeamformer(hw_collectedWave.re, hw_collectedWave.im, hw_Weights.re, hw_Weights.im, hw_result_re + i, hw_result_im + i);

	// 写下结果
    outFile1 << sw_result_re[i] << endl;
    outFile2 << hw_result_re[i] << endl;
    }
    // 关闭文件
    outFile1.close();
    outFile2.close();



/////////////////////////////////// 检测部分 ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

////////////////////////// 确定对比精度 //////////////////////////
#ifdef data_tb_DATA
#define ABS_ERR_THRESH 0.1
#else
#define ABS_ERR_THRESH 0.001
#endif

	unsigned err_cnt = 0;


//检查精度
   cout << "Checking results against a tolerance of " << ABS_ERR_THRESH << endl;
   cout << fixed <<
		   setprecision(5);
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

   return err_cnt;
}

