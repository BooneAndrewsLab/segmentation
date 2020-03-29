#include <ctime>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cmath>
#include "fmath.hpp"

using namespace std;

inline
double fastgamma3(double x) {
    double result, sum, num, denom;
    double one = 1;
    result = 0;
    sum = 0;

    if (x >= 0.01f && x <= one) {
        double const coef1 = 6.69569585833067770821885e+6;
        double const coef2 = 407735.985300921332020398;
        double const coef3 = 1.29142492667105836457693e+6;
        double const coef4 = 1.00000000000000000000000000e+00;
        double const coef5 = 6.69558099277749024219574e+6;
        double const coef6 = 4.27571696102861619139483e+6;
        double const coef7 = -2.89391642413453042503323e+6;
        double const coef8 = 317457.367152592609873458;

        num = coef1 + x * (coef2 + x * (coef3));//MiniMaxApproximation calculated by Mathematica 8
        denom = coef4 +
                x * (coef5 + x * (coef6 + x * (coef7 + x * (coef8))));//MiniMaxApproximation calculated by Mathematica 8
        return num / denom;
    } else if (1. >= one && x <= 171.) {
        double const coef_1 = 0.08333333333333333333333333;
        double const coef_2 = 0.00347222222222222222222222;
        double const coef_3 = -0.00268132716049382716049383;
        double const coef_4 = -0.000229472093621399176954733;
        double const coef_5 = 0.000784039221720066627474035;
        double const coef_6 = 0.0000697281375836585777429399;
        double const coef_7 = -0.000592166437353693882864836;
        double const coef_8 = -0.0000517179090826059219337058;
        double const coef_9 = 0.000839498720672087279993358;
        double const coef_10 = 0.0000720489541602001055908572;
        double ln, power, pi_sqrt, two_pi, arg;

        two_pi = 2 * M_PI;
        double invx = 1 / x;
        ln = exp(-x);
        arg = x - 0.5;

        power = pow(x, arg);
        pi_sqrt = sqrt(two_pi);

        sum = ln * power * pi_sqrt;
        result = one + invx * (coef_1 + invx * (coef_2 + invx * (coef_3 + invx * (coef_4 + invx * (coef_5 + invx *
                                                                                                            (coef_6 +
                                                                                                             invx *
                                                                                                             (coef_7 +
                                                                                                              invx *
                                                                                                              (coef_8 +
                                                                                                               invx *
                                                                                                               (coef_9 +
                                                                                                                invx *
                                                                                                                (coef_10))))))))));
    }

    return sum * result;
}

/**
 * Run mixture model on _im_ and save the result in _segmented_. Input size is m x n
 */
void cmm(const unsigned short *im, unsigned char *segmented, unsigned int M, unsigned int N, int nucNum) {
    const unsigned char components = 3; // Number of components
    const unsigned char component_tdist[components] = {2, 2, 1}; // t-distributions per component
    const double var = 1e8;
    const double learning_rate = 1e-6;
    const int m_nuc_size = 210;
    const int m_cell_size = 1900;
    const std::size_t image_size = M * N;

    // Matrix indexers, 2D, 3D, 4D
#define i2(i, j)       (((i) * N) + (j))
#define i3(i, j, k)    (((k) * M * N) + ((i) * N) + (j))
#define i4(i, j, k, l) (((k) * 2 * M * N) + ((l) * M * N) + ((i) * N) + (j))

    unsigned char max_iterations = 100; // # Maximum iterations
    unsigned char iter = 0;
    std::size_t i, j, k, m, n; // reserved iterators
    std::size_t a, b, c, d; // box blur indexes
    std::size_t imn, imnk, imnkj; // matrix index

    double x, y, z, t1, t2, t3, t4;
    double B = 12;
    double Bold = 0;

    double U[3][2] = {{0.,     500.},
                      {17500., 22500.},
                      {50000., 0.}}; // Initialize means
    double S[3][2] = {{var, var},
                      {var, var},
                      {var, 0.}}; // Initialize variance
    double V[3][2] = {{1.,   1.},
                      {100., 100.},
                      {100., 0.}}; // Initial degree of freedom for t-distribution
    double E[3][2] = {{.5, .5},
                      {.5, .5},
                      {1., 0.}}; // Mixing proportions for t-distributions
    double LogLike[max_iterations];

    // Expected nucleus and cell area
    double nucleus_area = 1. * nucNum * m_nuc_size / image_size;
    double cytoplasm_area = 1. * nucNum * m_cell_size / image_size;
    double background_area = 1. - cytoplasm_area;

    if (background_area < 0.1) {
        background_area = 0.1;
        cytoplasm_area = 1. - background_area - nucleus_area;
    }

    if (nucleus_area < 0.1) {
        nucleus_area = 0.1;
        cytoplasm_area = 0.2;
        background_area = 1. - nucleus_area - cytoplasm_area;
    }

    double area[3] = {background_area, cytoplasm_area, nucleus_area};

    unsigned short im_max = *max_element(im, im + image_size);
    x = pow(2, 16) / im_max;

    std::size_t size_4d = image_size * components * 2 * sizeof(double);
    auto *Y = (double *) malloc(size_4d);
    auto *u = (double *) malloc(size_4d);
    auto *StudPDFVal = (double *) malloc(size_4d);

    std::size_t size_3d = image_size * components * sizeof(double);
    auto *AveLocZ = (double *) malloc(size_3d);
    auto *MP = (double *) malloc(size_3d);
    auto *Z = (double *) malloc(size_3d);

    std::size_t size_2d = image_size * sizeof(double);
    auto *temp = (double *) malloc(size_2d);
    auto *sumYk = (double *) malloc(size_2d);
    auto *sumZ = (double *) malloc(size_2d);
    auto *sumMP = (double *) malloc(size_2d);
    auto *im_ad = (double *) malloc(size_2d);

    for (i = 0; i < image_size; i++) {
        im_ad[i] = round(im[i] * x);
    }

    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < components; k++)
                MP[i3(i, j, k)] = area[k];

    while (iter < max_iterations) {
        int start_s = clock();
        std::fill(sumZ, sumZ + image_size, 0);
        std::fill(sumMP, sumMP + image_size, 0);

        // E-step
        for (k = 0; k < components; k++) {
            // Reset temp accumulator back to 0
            std::fill(temp, temp + image_size, 0);
            std::fill(sumYk, sumYk + image_size, 0);

            for (m = 0; m < component_tdist[k]; m++) {
                double exponent = (V[k][m] + 1) / 2.;
                double dofcovar = V[k][m] * S[k][m];
                double c3 = (fastgamma3((V[k][m] / 2.) + 0.5) * pow(S[k][m], -0.5)) /
                            (sqrt(V[k][m] * M_PI) * fastgamma3(V[k][m] / 2.));

                for (i = 0; i < M; i++) {
                    for (j = 0; j < N; j++) {
                        imn = i2(i, j);
                        imnk = i3(i, j, k);
                        imnkj = i4(i, j, k, m);

                        StudPDFVal[imnkj] = c3 / pow(1 + pow(im_ad[imn] - U[k][m], 2) / dofcovar, exponent);
                        temp[imn] += E[k][m] * StudPDFVal[imnkj];

                        Y[imnkj] = E[k][m] * StudPDFVal[imnkj];
                        u[imnkj] = (V[k][m] + 1) / (V[k][m] + pow(im_ad[imn] - U[k][m], 2) / S[k][m]);
                        sumYk[imn] += Y[imnkj];

                        // run this calculation on last loop when all accumulators are full
                        if (m == component_tdist[k] - 1) {
                            Z[imnk] = MP[imnk] * temp[imn];
                            sumZ[imn] += Z[imnk];

                            for (n = 0; n < component_tdist[k]; n++) {
                                Y[i4(i, j, k, n)] /= sumYk[imn];
                            }

                            if (k == components - 1) {
                                for (n = 0; n < components; n++) {
                                    Z[i3(i, j, n)] /= sumZ[imn];
                                }
                            }
                        }
                    }
                }
            }
        }

        // M-step
        for (k = 0; k < components; k++) {
            for (m = 0; m < component_tdist[k]; m++) {
                x = 0;
                y = 0;
                t1 = 0;
                t2 = 0;

                for (i = 0; i < M; i++) {
                    for (j = 0; j < N; j++) {
                        imn = i2(i, j);
                        imnk = i3(i, j, k);
                        imnkj = i4(i, j, k, m);

                        z = Z[imnk] * Y[imnkj];
                        // for V
                        t1 += z;

                        z *= u[imnkj];
                        // for U
                        x += z * im_ad[imn];
                        y += z;

                        // for V
                        t3 = 0;
                        for (n = 0; n < component_tdist[k]; n++) {
                            t3 += Y[i4(i, j, k, n)];
                        }
                        t2 += Z[imnk] * t3;

                        if (m == 0) {
                            // image blurring submatrix sum
                            if (i == 0) {
                                temp[imn] = Z[imnk]; // Copy first row
                            } else {
                                temp[imn] = Z[imnk] + temp[i2(i - 1, j)]; // sum columns
                            }
                        }
                    }
                }

                U[k][m] = x / y;
                E[k][m] = t1 / t2;
            }

            // image blurring submatrix sum rows
            for (i = 0; i < M; i++) {
                for (j = 0; j < N; j++) {
                    if (j > 0) {
                        temp[i2(i, j)] += temp[i2(i, j - 1)];
                    }
                }
            }

            // Blurring using submatrix sum
            for (i = 0; i < M; i++) {
                for (j = 0; j < N; j++) {
                    imn = i2(i, j);
                    imnk = i3(i, j, k);

                    a = i - 3;
                    b = i + 3;
                    c = j - 3;
                    d = j + 3;

                    a = a > M ? 0 : a; // Remember, it's an unsigned int -- tli
                    b = b >= M ? M - 1 : b; // rbi
                    c = c > N ? 0 : c; // tlj
                    d = d >= N ? N - 1 : d; // rbj

                    AveLocZ[imnk] = temp[i2(b, d)];
                    if (a > 0) AveLocZ[imnk] -= temp[i2(a - 1, d)];
                    if (c > 0) AveLocZ[imnk] -= temp[i2(b, c - 1)];
                    if (a > 0 && c > 0) AveLocZ[imnk] += temp[i2(a - 1, c - 1)];

                    AveLocZ[imnk] /= 49;
                    MP[imnk] = fmath::expd(B * AveLocZ[imnk]);
                    sumMP[imn] += MP[imnk];
                }
            }
        }

        for (k = 0; k < components; k++) {
            for (i = 0; i < M; i++) {
                for (j = 0; j < N; j++) {
                    MP[i3(i, j, k)] /= sumMP[i2(i, j)];
                }
            }

            for (m = 0; m < component_tdist[k]; m++) {
                if (S[k][m] > 500) {
                    t1 = 0;
                    t2 = 0;

                    for (i = 0; i < M; i++) {
                        for (j = 0; j < N; j++) {
                            imn = i2(i, j);
                            imnk = i3(i, j, k);
                            imnkj = i4(i, j, k, m);

                            t3 = im_ad[imn] - U[k][m];
                            t1 += Z[imnk] * Y[imnkj] * u[imnkj] * t3 * t3;
                            t2 += Z[imnk] * Y[imnkj];
                        }
                    }
                    S[k][m] = t1 / t2;
                } else {
                    S[k][m] = 500;
                }
            }
        }

        int maxU1 = 0;
        if (U[1][0] < U[1][1]) {
            maxU1 = 1;
        }

        if (U[2][0] < U[1][maxU1]) {
            double Utmp = U[2][0];
            U[2][0] = U[1][maxU1];
            U[1][maxU1] = Utmp;

            double Stmp = S[2][0];
            S[2][0] = S[1][maxU1];
            S[1][maxU1] = Stmp;

            double Vtmp = V[2][0];
            V[2][0] = V[1][maxU1];
            V[1][maxU1] = Vtmp;
        }

        while (fabs(B - Bold) > 0.05) {
            Bold = B;
            t4 = 0;

            for (i = 0; i < M; i++) {
                for (j = 0; j < N; j++) {
                    t1 = 0;
                    t2 = 0;

                    for (k = 0; k < components; k++) {
                        imnk = i3(i, j, k);
                        t3 = fmath::expd(B * AveLocZ[imnk]);
                        t1 += t3 * AveLocZ[imnk];
                        t2 += t3;
                    }

                    t3 = t1 / t2;
                    t1 = 0;

                    for (k = 0; k < components; k++) {
                        imnk = i3(i, j, k);
                        t1 += (AveLocZ[imnk] - t3) * Z[imnk];
                    }

                    t4 += -t1;
                }
            }

            B = Bold - t4 * learning_rate;
        }

        LogLike[iter] = 0;
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                t1 = 0;
                for (k = 0; k < components; k++) {
                    t2 = 0;
                    for (m = 0; m < component_tdist[k]; m++) {
                        t2 += E[k][m] * StudPDFVal[i4(i, j, k, m)];
                    }

                    t1 += MP[i3(i, j, k)] * t2;
                }

                LogLike[iter] += log(t1);
            }
        }

        int stop_s = clock();

        if (iter > 0) {
            printf("Iterations = %d LogLikelihood = %f DiffLikelihood = %f Time = %fms\n", iter, LogLike[iter],
                   fabs(LogLike[iter - 1] - LogLike[iter]), (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000);

            if (fabs(LogLike[iter - 1] - LogLike[iter]) < 3000) break;
        }

        iter++;
    }

    free(Y);
    free(u);
    free(StudPDFVal);

    free(AveLocZ);
    free(MP);

    free(temp);
    free(sumYk);
    free(sumZ);
    free(sumMP);
    free(im_ad);

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            imn = i2(i, j);
            segmented[imn] = 0;
            t1 = Z[i3(i, j, 0)];
            for (k = 1; k < components; k++) {
                imnk = i3(i, j, k);
                if (Z[imnk] > t1) {
                    segmented[imn] = k;
                    t1 = Z[imnk];
                }
            }
        }
    }

    free(Z);
}
