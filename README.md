#include "test_qr_inverse.hpp"
#include "kernel_qr_inverse.hpp"

#include "src/utils.hpp"
#include "hw/utils/x_matrix_utils.hpp"
#include "src/matrix_test_utils.hpp"
#include "src/type_test_utils.hpp"

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <cmath> 

// ---------------------------------------------------------------------------------------------
// Main test program
// ---------------------------------------------------------------------------------------------


int main() {
    
    long unsigned num_tests = 1;
    // 每类矩阵有多少个Test case?
    // long unsigned num_tests = (ROWSCOLSA >= 16 ? 5 : 20); // Small default for HLS
    unsigned int debug = 1;
    // 是否输出矩阵    

    double ratio_threshold = 30.0;
    // ratio多少算是差别太大

    double mat_type = 0; // Specify to only run a single matrix type.
	// 0:全测试
	// n:第n个测试

    unsigned print_precision = 10;

   unsigned allowed_ulp_mismatch = 0;
	//调用matrices equal的参数

	// 变量定义
    // ====================================================================
    int qr_inverse_return = 0; // Return code from hls::qr_inverse
    // 求逆成功了吗


    // Matrix arrays
    MATRIX_IN_T A[ROWSCOLSA][ROWSCOLSA]; // The input array.  Cast from A_generated

    MATRIX_OUT_T Weights[ROWSCOLSA][ROWSCOLSA];          // The inverse result from the DUT
    MATRIX_OUT_T Weights_expected[ROWSCOLSA][ROWSCOLSA]; // The inverse result from LAPACK in target format

    // Test variables
    QR_INV_TYPE A_cast[ROWSCOLSA][ROWSCOLSA]; // Cast A back to LAPACK compatible type for analysis
    // A切换到另一种数据类型，用于求A的norm

    QR_INV_TYPE I_delta[ROWSCOLSA][ROWSCOLSA];
	// 差矩阵1

    MATRIX_OUT_T Weights_delta_vs_lapack[ROWSCOLSA]
                                         [ROWSCOLSA]; // The difference between the DUT's Weights and LAPACK's Weights
	// 差矩阵2


    // Non-complex type
    double A_norm;
    double Weights_norm;


    unsigned int pass_fail = 0; // Pass=0 Fail =1
	// 程序总结果

    bool matched_lapack_Weights;
    // 是否一致
    
    QR_INV_BASE_TYPE I_DUT_ratio, I_LAPACK_ratio;
    // 有多不一致


	// 检测输入数据是什么数据类型并打印
	// 好像不影响程序 没什么卵用
    printf("Running %lu %s tests per matrix type on %d x %d matrices\n", num_tests,
           x_is_float(Weights_expected[0][0])
               ? "single precision"
               : x_is_double(Weights_expected[0][0])
                     ? "double precision"
                     : x_is_fixed(Weights_expected[0][0]) ? "fixed point" : "Unknown type",
           ROWSCOLSA, ROWSCOLSA);


                // 循环读文件
                // ====================================================================
    //---------------- New code post-test review ----------
    for (unsigned int imat = 0; imat < NUM_MAT_TYPES; imat++) {
        // Test which matrix type to run
        if ((mat_type == 0 || imat + 1 == mat_type) && (skip_mat_type[0] - 1 != imat && skip_mat_type[1] - 1 != imat)) {
            for (long unsigned i = 0; i < num_tests; i++) {
                if ((imat == 10) && i > 0) {
                    // Skip the too large one
                    break;
                }

====================================================================
                // Read input matrix and golden matrix from files
                // 找根路径
                std::string data_path = std::string(DATA_PATH);
                std::string base_path;


				// 是哪个数据类型的，从哪读数据
                if (x_is_complex(A[0][0]) == true) {
                    base_path = data_path.append("/complex/");
                } else {
                    base_path = data_path.append("/float/");
                }
            
                std::string file_A =
                    base_path + "A_matType_" + std::to_string(imat + 1) + "_" + std::to_string(i) + ".txt";
                std::string file_Weights =
                    base_path + "TestPoint1_A_matType_" + std::to_string(imat + 1) + "_" + std::to_string(i) + ".txt";


				// 文件size, 读取的参数
                int A_size = ROWSCOLSA * ROWSCOLSA;
                int Weights_size = ROWSCOLSA * ROWSCOLSA;
                
				// 数组指针, 读取的参数
                MATRIX_IN_T* A_ptr = reinterpret_cast<MATRIX_IN_T*>(A);
                MATRIX_OUT_T* Weights_ptr = reinterpret_cast<MATRIX_OUT_T*>(Weights_expected);


				// 读取
                readTxt(file_A, A_ptr, A_size);
                readTxt(file_Weights, Weights_ptr, Weights_size);


                // cast
                for (int r = 0; r < ROWSCOLSA; r++) {
                    for (int c = 0; c < ROWSCOLSA; c++) {
                        A_cast[r][c] = A[r][c];
                    }
                }


                // 进出DUT的流
                hls::stream<MATRIX_IN_T> matrixAStrm;
                hls::stream<MATRIX_OUT_T> matrixMMSEH;


                // 写入流
                for (int r = 0; r < ROWSCOLSA; r++) {
                    for (int c = 0; c < ROWSCOLSA; c++) {
                        matrixAStrm.write(A[r][c]);
                    }
                }

	
                // 噪音
                float var_Noise = 1;

                // 文件size, 读取的参数
                qr_inverse_return = kernel_qr_inverse_0(matrixAStrm, matrixMMSEH, var_Noise);


                // 读出流
                for (int r = 0; r < ROWSCOLSA; r++) {
                    for (int c = 0; c < ROWSCOLSA; c++) {
                    	matrixMMSEH.read(Weights[r][c]);
                    }
                }
                



//                 Check for NaNs in result
                if (anyNaN<ROWSCOLSA, ROWSCOLSA>(Weights) == 1 && !testing_singular_matrix && imat != 10) {
                    printf("ERROR: Caught NaN in Weights\n");
                    xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T, xf::solver::NoTranspose>(
                        Weights, "   ", print_precision, 0);
                    printf("TB:Fail\n");
                    return (32);
                }


                // Test results
                // ====================================================================

                // Basic check cell by cell check based on allowed_ulp_mismatch value.
                // 求差方法1
                // 用一个函数计算与目标结果的差
                msub<ROWSCOLSA, ROWSCOLSA, MATRIX_IN_T, QR_INV_TYPE>(Weights, Weights_expected, I_delta);
                
                // 求差方法2
                // 用一个函数计算与目标结果是否一致，差值多少
                matched_lapack_Weights = are_matrices_equal<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T>(
                    (MATRIX_OUT_T*)Weights, (MATRIX_OUT_T*)Weights_expected, allowed_ulp_mismatch,
                    (MATRIX_OUT_T*)Weights_delta_vs_lapack);


                // 求norm
                A_norm = norm1_dbl<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(A_cast);

                Weights_norm = norm1_dbl<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(Weights_expected);

                I_delta_norm = norm1<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(I_delta);
                
                // 验证norm的异常值
            
                if (isinf(A_norm)) {
                    // Should never be Inf - if it is, we've probably overflowed
                    printf("ERROR: Caught unexpected Inf for A_norm\n");
                    printf("TB:Fail\n");
                    return (4);
                }

                if (isinf(Weights_norm)) {
                    // Should never be Inf - if it is, we've probably overflowed
                    printf("ERROR: Caught unexpected Inf for Weights_norm\n");
                    printf("TB:Fail\n");
                    return (5);
                }


                // norm的比率
                 I_DUT_ratio =(double)I_delta_norm / (double)A_norm;


                // Check that the norm values are OK and we are not comparing two bad ratios
                // 检验I_DUT_ratio异常值

                if (isnan(I_DUT_ratio)) {
                    // Should only be NaN if A_norm was zero, so check that
                    if (A_norm != 0) {
                        printf("ERROR: Caught unexpected NaN for I_DUT_ratio\n");
                        printf("TB:Fail\n");
                        return (6);
                    }
                }

                if (isinf(I_DUT_ratio)) {
                 // Should never be Inf
                	printf("ERROR: Caught unexpected Inf for I_DUT_ratio\n");
                    printf("TB:Fail\n");
                return (8);
                }

                if (I_DUT_ratio < 0) {
                    // Should never be less than zero - if it is, it's either an error code or something went badly
                    // wrong
                    printf("ERROR: Caught unexpected negative I_DUT_ratio\n");
                    printf("TB:Fail\n");
                    return (12);
                }


                // Determine if pass or fail.
                // o Check DUT ratio against test threshold, default taken from LAPACK
                if (I_DUT_ratio > ratio_threshold ) {
                    std::cout << "ERROR: I_DUT_ratio(" << I_DUT_ratio << ") > ratio_threshold(" << ratio_threshold
                              << "). " << std::endl;
                    pass_fail = 1; // Test run fails
                }

                // Print matrices for debug
                // 输出所有矩阵
                if ( debug > 0 || I_DUT_ratio > ratio_threshold) {
                    printf("testing_singular_matrix = %s \n", testing_singular_matrix ? "true" : "false");
                    printf("  Channel Matrix=\n");
                    xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_IN_T, xf::solver::NoTranspose>(
                        A, "   ", print_precision, 1);
                    printf("  Weights=\n");
                    xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T, xf::solver::NoTranspose>(
                        Weights, "   ", print_precision, 1);
                    printf("  Weights_expected (var_Noise=1)=\n");
                    xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T, xf::solver::NoTranspose>(
                        Weights_expected, "   ", print_precision, 1);
                    printf("  Weights_delta=\n");
                    xf::solver::print_matrix<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, xf::solver::NoTranspose>(
                        I_delta, "   ", print_precision, 1);
                    printf("  ratio=\n");
                    std::cout<<I_DUT_ratio<<std::endl;
                    printf("  matched? ");
                    std::cout<< matched_lapack_Weights <<std::endl;
                    }

            } // End of test loop
            printf("\n");
        } // Type test
    }     // Matrix type loop


    if (pass_fail == 1) {
        std::cout << "TB:Fail" << std::endl;
    } else {
        std::cout << "TB:Pass" << std::endl;
    }
    std::cout << "" << std::endl;
    return (pass_fail);
}
