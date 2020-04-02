#include "fastmm.h"

/**
 * Run mixture model on _im_ and save the result in _segmented_. Input size is m x n
 */
void cmm(const unsigned short *im, unsigned char *segmented, unsigned int dim1, unsigned int dim2, int nuc_num,
         CmmConfig config, cbfunc callback, void *user_data) {
    const size_t image_size = dim1 * dim2;

    // Matrix indexers, 2D, 3D, 4D
#define i2(i, j)       (((i) * dim2) + (j))
#define i3(i, j, k)    (((k) * dim1 * dim2) + ((i) * dim2) + (j))
#define i4(i, j, k, l) (((k) * 2 * dim1 * dim2) + ((l) * dim1 * dim2) + ((i) * dim2) + (j))

    unsigned char iter = 0;
    size_t i, j, k, m, n; // reserved iterators
    size_t a, b, c, d; // box blur indexes
    size_t imn, imnk, imnkj; // matrix index

    double x, y, z, t1, t2, t3, t4;
    double B = 12;
    double Bold = 0;

    const unsigned int components = 3; // Number of components
    const unsigned int component_tdist[components] = {2, 2, 1}; // t-distributions per component

    double U[3][2] = {{0.,     500.},
                      {17500., 22500.},
                      {50000., 0.}}; // Initialize means
    double S[3][2] = {{config.var, config.var},
                      {config.var, config.var},
                      {config.var, 0.}}; // Initialize variance
    double V[3][2] = {{1.,   1.},
                      {100., 100.},
                      {100., 0.}}; // Initial degree of freedom for t-distribution
    double E[3][2] = {{.5, .5},
                      {.5, .5},
                      {1., 0.}}; // Mixing proportions for t-distributions
    double LogLike[config.max_iterations];

    // Expected nucleus and cell area
    double nucleus_area = 1. * nuc_num * config.m_nuc_size / image_size;
    double cytoplasm_area = 1. * nuc_num * config.m_cell_size / image_size;
    double background_area = 1. - cytoplasm_area - nucleus_area;

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

    vector<double> Y(image_size * components * 2);
    vector<double> u(image_size * components * 2);
    vector<double> StudPDFVal(image_size * components * 2);

    vector<double> AveLocZ(image_size * components);
    vector<double> MP(image_size * components);
    vector<double> Z(image_size * components);

    vector<double> temp(image_size);
    vector<double> sumZ(image_size);
    vector<double> sumMP(image_size);
    vector<double> im_ad(image_size);

    for (i = 0; i < image_size; i++) {
        im_ad[i] = round(im[i] * x);
    }

    for (i = 0; i < dim1; i++)
        for (j = 0; j < dim2; j++)
            for (k = 0; k < components; k++)
                MP[i3(i, j, k)] = area[k];

    double exponent, dofcovar, c3;
    int start_s;
    while (iter < config.max_iterations) {
        start_s = clock();
        fill(sumZ.begin(), sumZ.end(), 0);
        fill(sumMP.begin(), sumMP.end(), 0);

        // E-step
        for (k = 0; k < components; k++) {
            // Reset temp accumulator back to 0
            fill(temp.begin(), temp.end(), 0);

            for (m = 0; m < component_tdist[k]; m++) {
                exponent = (V[k][m] + 1) / 2.;
                dofcovar = V[k][m] * S[k][m];
                c3 = (fmath::fastgamma3((V[k][m] / 2.) + 0.5) * pow(S[k][m], -0.5)) /
                     (sqrt(V[k][m] * M_PI) * fmath::fastgamma3(V[k][m] / 2.));

                for (i = 0; i < dim1; i++) {
                    for (j = 0; j < dim2; j++) {
                        imn = i2(i, j);
                        imnk = i3(i, j, k);
                        imnkj = i4(i, j, k, m);

                        t1 = pow(im_ad[imn] - U[k][m], 2);

                        StudPDFVal[imnkj] = c3 / pow(1 + t1 / dofcovar, exponent);
                        Y[imnkj] = E[k][m] * StudPDFVal[imnkj];
                        u[imnkj] = (V[k][m] + 1) / (V[k][m] + t1 / S[k][m]);

                        temp[imn] += Y[imnkj];

                        // run this calculation on last loop when all accumulators are full
                        if (m == component_tdist[k] - 1) {
                            Z[imnk] = MP[imnk] * temp[imn];
                            sumZ[imn] += Z[imnk];

                            for (n = 0; n < component_tdist[k]; n++) {
                                Y[i4(i, j, k, n)] /= temp[imn];
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

        if (callback != NULL) {
            callback(iter, Z.data(), components, dim1, dim2, user_data);
        }

        // M-step
        for (k = 0; k < components; k++) {
            for (m = 0; m < component_tdist[k]; m++) {
                x = 0;
                y = 0;
                t1 = 0;
                t2 = 0;

                for (i = 0; i < dim1; i++) {
                    for (j = 0; j < dim2; j++) {
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
            for (i = 0; i < dim1; i++) {
                for (j = 0; j < dim2; j++) {
                    if (j > 0) {
                        temp[i2(i, j)] += temp[i2(i, j - 1)];
                    }
                }
            }

            // Blurring using submatrix sum
            for (i = 0; i < dim1; i++) {
                for (j = 0; j < dim2; j++) {
                    imn = i2(i, j);
                    imnk = i3(i, j, k);

                    a = i - 3;
                    b = i + 3;
                    c = j - 3;
                    d = j + 3;

                    a = a > dim1 ? 0 : a; // Remember, it's an unsigned int -- tli
                    b = b >= dim1 ? dim1 - 1 : b; // rbi
                    c = c > dim2 ? 0 : c; // tlj
                    d = d >= dim2 ? dim2 - 1 : d; // rbj

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
            for (i = 0; i < dim1; i++) {
                for (j = 0; j < dim2; j++) {
                    MP[i3(i, j, k)] /= sumMP[i2(i, j)];
                }
            }

            for (m = 0; m < component_tdist[k]; m++) {
                if (S[k][m] > 500) {
                    t1 = 0;
                    t2 = 0;

                    for (i = 0; i < dim1; i++) {
                        for (j = 0; j < dim2; j++) {
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

            for (i = 0; i < dim1; i++) {
                for (j = 0; j < dim2; j++) {
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

            B = Bold - t4 * config.learning_rate;
        }

        LogLike[iter] = 0;
        for (i = 0; i < dim1; i++) {
            for (j = 0; j < dim2; j++) {
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

        if (iter > 0) {
            if (config.debug) {
                printf("Iterations = %d LogLikelihood = %f DiffLikelihood = %f Time=%.3fms\n", iter, LogLike[iter],
                       fabs(LogLike[iter - 1] - LogLike[iter]),
                       (clock() - start_s) / double(CLOCKS_PER_SEC) * 1000
                );
            }

            if (fabs(LogLike[iter - 1] - LogLike[iter]) < 3000) break;
        }

        iter++;
    }

    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim2; j++) {
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
}
