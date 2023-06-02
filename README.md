

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
#include <cmath> // for abs

// ---------------------------------------------------------------------------------------------
// Main test program
// ---------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    // Variables set by command line
    // ====================================================================
    long unsigned num_tests = 1;
    // long unsigned num_tests = (ROWSCOLSA >= 16 ? 5 : 20); // Small default for HLS
    unsigned allowed_ulp_mismatch = 0;
    unsigned int debug = 1;
    double ratio_threshold = 30.0;
    double mat_type = 0; // Specify to only run a single matrix type.
    unsigned int skip_mat_type[2] = {0};
    int skip_mat_type_itr = 0;
    unsigned print_precision = 10;


    // Parse command line options
    // ====================================================================
    std::vector<std::string> args(argv + 1, argv + argc);
    for (std::vector<std::string>::iterator i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            std::cout << "Syntax: main.exe [-num_tests <long unsigned> -mat_type <unsigned> -skip_mat_type <unsigned>]"
                      << std::endl;
            return 0;
        } else if (*i == "-num_tests") {
            num_tests = (long unsigned)atol((*++i).c_str());
            printf("num_tests as long unsigned = %lu\n", num_tests);
        } else if (*i == "-max_ulp_mismatch_in_l") {
            allowed_ulp_mismatch = (unsigned)atol((*++i).c_str());
        } else if (*i == "-prec") {
            print_precision = (unsigned)atol((*++i).c_str());
        } else if (*i == "-ratio_threshold") {
            ratio_threshold = (double)atol((*++i).c_str());
        } else if (*i == "-mat_type") {
            mat_type = (unsigned)atol((*++i).c_str());
        } else if (*i == "-skip_mat_type") {
            skip_mat_type[skip_mat_type_itr] = (unsigned)atol((*++i).c_str());
            skip_mat_type_itr++;
        } else if (*i == "-debug") {
            debug = (unsigned)atol((*++i).c_str());
        } else {
            printf("Unknown command line option %s.  Try again\n", i->c_str());
            exit(1);
        }
    }

    if (skip_mat_type[0] != 0) {
        std::cout << "INFO: Skipping matrix type " << skip_mat_type[0] << std::endl;
    }
    if (skip_mat_type[1] != 0) {
        std::cout << "INFO: Skipping matrix type " << skip_mat_type[1] << std::endl;
    }

    int qr_inverse_return = 0; // Return code from hls::qr_inverse
    bool testing_singular_matrix = false;

    // Matrix arrays
    MATRIX_IN_T A[ROWSCOLSA][ROWSCOLSA]; // The input array.  Cast from A_generated

    MATRIX_OUT_T Weights[ROWSCOLSA][ROWSCOLSA];          // The inverse result from the DUT
    MATRIX_OUT_T Weights_expected[ROWSCOLSA][ROWSCOLSA]; // The inverse result from LAPACK in target format

    MATRIX_OUT_T I[ROWSCOLSA][ROWSCOLSA]; // The identity matrix to compare against

    // Test variables
    QR_INV_TYPE A_cast[ROWSCOLSA][ROWSCOLSA]; // Cast A back to LAPACK compatible type for analysis
    MATRIX_OUT_T Weights_delta_vs_lapack[ROWSCOLSA]
                                         [ROWSCOLSA]; // The difference between the DUT's Weights and LAPACK's Weights
    MATRIX_IN_T I_restored[ROWSCOLSA][ROWSCOLSA];     // I recreated from Weights*A
    MATRIX_IN_T I_restored_lapack[ROWSCOLSA][ROWSCOLSA]; // I recreated from LAPACK's Weights*A

    QR_INV_TYPE I_delta[ROWSCOLSA]
                       [ROWSCOLSA]; // The delta values will be passed to xLANSY so need to use LAPACK compatible types
    QR_INV_TYPE I_delta_lapack[ROWSCOLSA][ROWSCOLSA]; //

    // Non-complex type
    double A_norm;
    double Weights_norm;
    QR_INV_BASE_TYPE I_delta_norm;
    QR_INV_BASE_TYPE I_delta_lapack_norm;
    QR_INV_BASE_TYPE I_DUT_ratio, I_LAPACK_ratio;

    double I_imat_max_ratio_diff[NUM_MAT_TYPES];
    double I_imat_min_ratio_diff[NUM_MAT_TYPES];
    double I_imat_max_ratio[NUM_MAT_TYPES];
    double I_imat_min_ratio[NUM_MAT_TYPES];

    unsigned int I_ratio_same = 0;
    unsigned int I_ratio_better = 0;
    unsigned int I_ratio_worse = 0;
    unsigned int I_imat_ratio_same[NUM_MAT_TYPES];
    unsigned int I_imat_ratio_better[NUM_MAT_TYPES];
    unsigned int I_imat_ratio_worse[NUM_MAT_TYPES];

    // Zero values
    for (int i = 0; i < NUM_MAT_TYPES; i++) {
        I_imat_max_ratio_diff[i] = 0;
        I_imat_min_ratio_diff[i] = 0;
        I_imat_ratio_same[i] = 0;
        I_imat_ratio_better[i] = 0;
        I_imat_ratio_worse[i] = 0;
        I_imat_max_ratio[i] = 0;
        I_imat_min_ratio[i] = ratio_threshold;
    }

    double I_ratio_difference = 0;
    unsigned int pass_fail = 0; // Pass=0 Fail =1

    bool matched_lapack_Weights;

    printf("Running %lu %s tests per matrix type on %d x %d matrices\n", num_tests,
           x_is_float(Weights_expected[0][0])
               ? "single precision"
               : x_is_double(Weights_expected[0][0])
                     ? "double precision"
                     : x_is_fixed(Weights_expected[0][0]) ? "fixed point" : "Unknown type",
           ROWSCOLSA, ROWSCOLSA);


    // Generate results table header
    std::cout << "RESULTS_TABLE,Test,IMAT,Weights Matching,DUT Ratio,LAPACK Ratio,Relative Ratio Difference"
              << std::endl;

    // Create I to compare against later
    for (int r = 0; r < ROWSCOLSA; r++) {
        for (int c = 0; c < ROWSCOLSA; c++) {
            if (r == c) {
                I[r][c] = 1.0;
            } else {
                I[r][c] = 0.0;
            }
        }
    }

    //---------------- New code post-test review ----------
    for (unsigned int imat = 0; imat < NUM_MAT_TYPES; imat++) {
        // Test which matrix type to run
        if ((mat_type == 0 || imat + 1 == mat_type) && (skip_mat_type[0] - 1 != imat && skip_mat_type[1] - 1 != imat)) {
            for (long unsigned i = 0; i < num_tests; i++) {
//                if (imat >= 4) {
//                    testing_singular_matrix = true;
//                } else {
//                    testing_singular_matrix = false;
//                }
                if ((imat == 10) && i > 0) {
                    // Skip the too large one
                    break;
                }

                // Get reference results
                // ====================================================================
                // Read input matrix and golden matrix from files
                std::string data_path = std::string(DATA_PATH);
                std::string base_path;

                if (x_is_complex(A[0][0]) == true) {
                    base_path = data_path.append("/complex/");
                } else {
                    base_path = data_path.append("/float/");
                }
                std::string file_A =
                    base_path + "A_matType_" + std::to_string(imat + 1) + "_" + std::to_string(i) + ".txt";
                std::string file_Weights =
                    base_path + "TestPoint1_A_matType_" + std::to_string(imat + 1) + "_" + std::to_string(i) + ".txt";

                int A_size = ROWSCOLSA * ROWSCOLSA;
                int Weights_size = ROWSCOLSA * ROWSCOLSA;

                MATRIX_IN_T* A_ptr = reinterpret_cast<MATRIX_IN_T*>(A);
                MATRIX_OUT_T* Weights_ptr = reinterpret_cast<MATRIX_OUT_T*>(Weights_expected);

                readTxt(file_A, A_ptr, A_size);
                readTxt(file_Weights, Weights_ptr, Weights_size);

                for (int r = 0; r < ROWSCOLSA; r++) {
                    for (int c = 0; c < ROWSCOLSA; c++) {
                        A_cast[r][c] = A[r][c];
                    }
                }

                hls::stream<MATRIX_IN_T> matrixAStrm;
                hls::stream<MATRIX_OUT_T> matrixMMSEH;


                for (int r = 0; r < ROWSCOLSA; r++) {
                    for (int c = 0; c < ROWSCOLSA; c++) {
                        matrixAStrm.write(A[r][c]);
                    }
                }

                float var_Noise = 1;

                qr_inverse_return = kernel_qr_inverse_0(matrixAStrm, matrixMMSEH, var_Noise);

                for (int r = 0; r < ROWSCOLSA; r++) {
                    for (int c = 0; c < ROWSCOLSA; c++) {
                    	matrixMMSEH.read(Weights[r][c]);
                    }
                }

                if (qr_inverse_return != 0 && !testing_singular_matrix) {
                    printf("ERROR: Input matrix was not singular, but QR Inverse thinks it is!\n");
                    printf("TB:Fail\n");
                    return (1);
//                    return (0);
                } else if (qr_inverse_return != 0 && testing_singular_matrix) {
                    printf("INFO: Singular matrix was correctly detected\n");
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
                matched_lapack_Weights = are_matrices_equal<ROWSCOLSA, ROWSCOLSA, MATRIX_OUT_T>(
                    (MATRIX_OUT_T*)Weights, (MATRIX_OUT_T*)Weights_expected, allowed_ulp_mismatch,
                    (MATRIX_OUT_T*)Weights_delta_vs_lapack);


//                mmult<false, false, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA>(Weights, A_cast,
//                                                                                                      I_restored);
//                mmult<false, false, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA, ROWSCOLSA>(
//                    Weights_expected, A_cast, I_restored_lapack);
//                msub<ROWSCOLSA, ROWSCOLSA, MATRIX_IN_T, QR_INV_TYPE>(I, I_restored, I_delta);
//                msub<ROWSCOLSA, ROWSCOLSA, MATRIX_IN_T, QR_INV_TYPE>(I, I_restored_lapack, I_delta_lapack);




                msub<ROWSCOLSA, ROWSCOLSA, MATRIX_IN_T, QR_INV_TYPE>(Weights, Weights_expected, I_delta);

                // REVISIT: is A_cast the appropriate format to use here?
                // norm1 as used in SPOT01
                A_norm = norm1_dbl<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(A_cast);

                if (isinf(A_norm)) {
                    // Should never be Inf - if it is, we've probably overflowed
                    printf("ERROR: Caught unexpected Inf for A_norm\n");
                    printf("TB:Fail\n");
                    return (4);
                }


                // REVISIT: which version of Weights should we use here?  Probably LAPACK's, since it's more likely to
                // be correct - hopefully!
                // norm1 as used in SPOT01
                Weights_norm = norm1_dbl<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(Weights_expected);


                if (isinf(Weights_norm)) {
                    // Should never be Inf - if it is, we've probably overflowed
                    printf("ERROR: Caught unexpected Inf for Weights_norm\n");
                    printf("TB:Fail\n");
                    return (5);
                }

                // norm1 as used in SPOT01
                I_delta_norm = norm1<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(I_delta);

                // norm1 as used in SPOT01
                I_delta_lapack_norm = norm1<ROWSCOLSA, ROWSCOLSA, QR_INV_TYPE, QR_INV_BASE_TYPE>(I_delta_lapack);

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


//                if (I_DUT_ratio == 0) {
//                    if (!(imat == 0 || imat == 1 || testing_singular_matrix)) {
//                        // Neither diagonal nor upper-triangular, so there should be error in the reconstruction, but we
//                        // didn't detect that, so fail
//                        printf("ERROR: Caught unexpected Zero for I_DUT_ratio in %d\n",imat+1);
//                        std::cout << "RESULTS_TABLE," << i << "," << imat + 1 << "," << matched_lapack_Weights << ","
//                                  << I_DUT_ratio << "," << I_LAPACK_ratio << "," << I_ratio_difference << std::endl;
//                        printf("TB:Fail\n");
//                        return (10);
//                    }
//                }

                if (I_DUT_ratio < 0) {
                    // Should never be less than zero - if it is, it's either an error code or something went badly
                    // wrong
                    printf("ERROR: Caught unexpected negative I_DUT_ratio\n");
                    printf("TB:Fail\n");
                    return (12);
                }

//                std::cout << "RESULTS_TABLE," << i << "," << imat + 1 << "," << matched_lapack_Weights << ","
//                          << I_DUT_ratio << "," << I_LAPACK_ratio << "," << I_ratio_difference << std::endl;
//
                // Determine if pass or fail.
                // o Check DUT ratio against test threshold, default taken from LAPACK
                if (I_DUT_ratio > ratio_threshold && !testing_singular_matrix) {
                    std::cout << "ERROR: I_DUT_ratio(" << I_DUT_ratio << ") > ratio_threshold(" << ratio_threshold
                              << "). I LAPACK ratio = " << I_LAPACK_ratio << std::endl;
                    pass_fail = 1; // Test run fails
                }

                // Print matrices for debug
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
