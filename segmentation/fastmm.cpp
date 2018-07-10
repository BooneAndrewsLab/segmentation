#include <boost/multi_array.hpp>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <random>
#include <cstdint>
#include <fstream>
#include "fmath.hpp"

using namespace std;

typedef boost::multi_array<double, 4> farray4;
typedef boost::multi_array<double, 3> farray3;
typedef boost::multi_array<double, 2> farray2;

#define BEGIN_ITER_K_M                  \
    for ( k = 0; k < K; k++ ) {         \
        for ( m = 0; m < Kj[k]; m++ ) {

#define BEGIN_ITER_IMAGE            \
    for ( i = 0; i < M; i++ ) {     \
        for ( j = 0; j < N; j++ ) {

#define END_ITER \
        }            \
    }                \

inline
double fastgamma3(double x)
{
    double result, sum, num, denom;
    double one = 1;
    result = 0;
    sum = 0;

    if(x >= 0.01f && x <= one) {
        double const coef1 = 6.69569585833067770821885e+6;
        double const coef2 = 407735.985300921332020398;
        double const coef3 = 1.29142492667105836457693e+6;
        double const coef4 = 1.00000000000000000000000000e+00;
        double const coef5 = 6.69558099277749024219574e+6;
        double const coef6 = 4.27571696102861619139483e+6;
        double const coef7 = -2.89391642413453042503323e+6;
        double const coef8 = 317457.367152592609873458;

        num = coef1+x*(coef2+x*(coef3));//MiniMaxApproximation calculated by Mathematica 8
        denom = coef4+x*(coef5+x*(coef6+x*(coef7+x*(coef8))));//MiniMaxApproximation calculated by Mathematica 8
        return num/denom;
    } else if ( 1. >= one && x <= 171.) {
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
        double ln,power,pi_sqrt,two_pi,arg;

        two_pi = 2*M_PI;
        double invx = 1/x;
        ln = exp(-x);
        arg = x-0.5;

        power = pow(x,arg);
        pi_sqrt = sqrt(two_pi);

        sum = ln*power*pi_sqrt;
        result = one+invx*(coef_1+invx*(coef_2+invx*(coef_3+invx*(coef_4+invx*(coef_5+invx*(coef_6+invx*(coef_7+invx*(coef_8+invx*(coef_9+invx*(coef_10))))))))));
    }

    return sum*result;
}

inline
double exp1(double x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}

inline
double exp2(double x) {
  x = 1.0 + x / 1024;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x;
  return x;
}

inline
double fexp(double x) {
    return x < 2 ? exp1(x) : exp2(x);
}

inline
double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

tuple<vector<unsigned short>, int, int> read() {
    ifstream in("./im_invert.csv");
    vector<unsigned short> fields;
    int width = 0, height = 0;

    if (in) {
        string line;
        while (getline(in, line)) {
            stringstream sep(line);
            string field;
            height++;
            width = 0;

            while (getline(sep, field, ',')) {
                fields.push_back(stoi(field));
                width++;
            }
        }
    }

    return make_tuple(fields, height, width);
}

/**
 * Run mixture model on _im_ and save the result in _segmented_. Input size is m x n
 */
int cmm (const unsigned short* im, unsigned char* segmented, unsigned int M, unsigned int N, int nucNum) {
    const unsigned char K = 3; // Number of components
    const unsigned char Kj[3] = {2,2,1}; // t-distributions per component
//    const unsigned char Nh = 7; // Neighborhood for averaging
    const double var = 1e8;
    const double learning_rate = 1e-6;
    const int m_nuc_size = 210;
    const int m_cell_size = 1900;

    unsigned char max_iterations = 100; // # Maximum iterations
    unsigned char iter = 0;
    std::size_t i, j, k, m, n; // reserved iterators
    std::size_t a, b, c, d; // box blur indexes
    std::size_t p, q; // random iterators

    double x, y, t1, t2, t3, t4;
    double B = 12;
    double Bold = 0;
    double Ptot = 0;

    fmath::PowGenerator fpow2(2);
    fmath::PowGenerator fpowm05(-0.5);

    double U[3][2] = {{ 0., 500.}, {17500., 22500.}, {50000., 0.}}; // Initialize means
    double S[3][2] = {{var,  var}, {   var,    var}, {   var, 0.}}; // Initialize variance
    double V[3][2] = { {1.,   1.}, {  100.,   100.}, {  100., 0.}}; // Initial degree of freedom for t-distribution
    double E[3][2] = { {.5,   .5}, {    .5,     .5}, {    1., 0.}}; // Mixing proportions for t-distributions
    double LogLike[max_iterations];

    // Expected nucleus and cell area
    double NucArea = nucNum * m_nuc_size / (M * N);
    double CytArea = nucNum * m_cell_size / (M * N);
    double BackArea = 1. - CytArea;

    if (BackArea < 0.1) {
        BackArea = 0.1;
        CytArea = 1. - BackArea - NucArea;
    }

    if (NucArea < 0.1) {
        NucArea = 0.1;
        CytArea = 0.2;
        BackArea = 1. - NucArea - CytArea;
    }

    double area[3] = {BackArea, CytArea, NucArea};
    // printf("Percent area - Background %.1f - Cytoplasm %.1f - Nucleus %.1f\n", area[0], area[1], area[2]);

    unsigned short im_max = *max_element(im,im+M*N);
    x = pow(2, 16) / im_max;

    farray4 Y(boost::extents[K][M][N][2]);
    farray4 u(boost::extents[K][M][N][2]);
    farray4 Temp1(boost::extents[K][M][N][2]);
    farray4 StudPDFVal(boost::extents[K][M][N][2]);

    farray3 AveLocZ(boost::extents[M][N][K]);
    farray3 MP(boost::extents[M][N][K]);
    farray3 Z(boost::extents[M][N][K]);

    farray2 temp(boost::extents[M][N]);
    farray2 sumYk(boost::extents[M][N]);
    farray2 sumZ(boost::extents[M][N]);
    farray2 sumMP(boost::extents[M][N]);
    farray2 im_ad(boost::extents[M][N]);

    for ( i = 0; i < M; i++ ) {
        for ( j = 0; j < N; j++ ) {
            im_ad[i][j] = round(im[i * N + j] * x);
        }
    }

    for ( i = 0; i < M; i++)
        for ( j = 0; j < N; j++)
            for ( k = 0; k < K; k++ )
                MP[i][j][k] = area[k];

/* USELESS CODE IN PYTHON IMPLEMENTATION??? KEEP HERE JUST TO AVOID ANOTHER 2 DAYS OF GOOGLING */
//    GiNaC::symbol v("v");
//    GiNaC::ex digamma_x_2;
//
//    digamma_x_2 = log(v)-1/(2*v)-1/(12*pow(v,2))+1/(120*pow(v,4))-1/(252*pow(v,6))+1/(240*pow(v,8))-5/(660*pow(v,10))+691/(32760*pow(v,12))-1/(12*pow(v,14));
//    digamma_x_2 = digamma_x_2.subs(v==v/2);

    while (iter < max_iterations) {
        int start_s = clock();
        std::fill( sumZ.origin(), sumZ.origin() + sumZ.num_elements(), 0 );
        std::fill( sumMP.origin(), sumMP.origin() + sumMP.num_elements(), 0 );

        // E-step
        for ( k = 0; k < K; k++ ) {
            // Reset temp accumulator back to 0
            std::fill( temp.origin(), temp.origin() + temp.num_elements(), 0 );
            std::fill( sumYk.origin(), sumYk.origin() + sumYk.num_elements(), 0 );

            for ( m = 0; m < Kj[k]; m++ ) {
                double exponent = (V[k][m] + 1) / 2.;
//                fmath::PowGenerator fpowm(exponent);
                double dofcovar = V[k][m] * S[k][m];
                double c3 = (fastgamma3((V[k][m] / 2.) + 0.5) * pow(S[k][m], -0.5)) / (sqrt(V[k][m] * M_PI) * fastgamma3(V[k][m] / 2.));
//                double c3 = (fastgamma3((V[k][m] / 2.) + 0.5) * fastPow(S[k][m], -0.5)) / (sqrt(V[k][m] * M_PI) * fastgamma3(V[k][m] / 2.));
//                double c3 = (fastgamma3((V[k][m] / 2.) + 0.5) * fpowm05.get(S[k][m])) / (sqrt(V[k][m] * M_PI) * fastgamma3(V[k][m] / 2.));

                BEGIN_ITER_IMAGE
                    StudPDFVal[k][i][j][m] = c3 / pow(1 + pow(im_ad[i][j] - U[k][m], 2) / dofcovar, exponent);
//                    StudPDFVal[k][i][j][m] = c3 / fastPow(1 + fastPow(im_ad[i][j] - U[k][m], 2) / dofcovar, exponent);
//                    StudPDFVal[k][i][j][m] = c3 / fpowm.get(1 + fpow2.get(im_ad[i][j] - U[k][m]) / dofcovar);
                    temp[i][j] += E[k][m] * StudPDFVal[k][i][j][m];
                    Y[k][i][j][m] = E[k][m] * StudPDFVal[k][i][j][m];
                    u[k][i][j][m] = (V[k][m] + 1) / (V[k][m] + pow(im_ad[i][j] - U[k][m], 2) / S[k][m]);
//                    u[k][i][j][m] = (V[k][m] + 1) / (V[k][m] + fastPow(im_ad[i][j] - U[k][m], 2) / S[k][m]);
//                    u[k][i][j][m] = (V[k][m] + 1) / (V[k][m] + fpow2.get(im_ad[i][j] - U[k][m]) / S[k][m]);
                    sumYk[i][j] += Y[k][i][j][m];

                    // run this calculation on last loop when all accumulators are full
                    if (m == Kj[k] - 1) {
                        Z[i][j][k] = MP[i][j][k] * temp[i][j];
                        sumZ[i][j] += Z[i][j][k];

                        for ( n = 0; n < Kj[k]; n++ ) {
                            Y[k][i][j][n] /= sumYk[i][j];
                        }

                        if (k == K - 1) {
                            for ( n = 0; n < K; n++ ) {
                                Z[i][j][n] /= sumZ[i][j];
                            }
                        }
                    }
                END_ITER
            }
        }

        // M-step
        for ( k = 0; k < K; k++ ) {
            for ( m = 0; m < Kj[k]; m++ ) {
                x = 0;
                y = 0;
                t1 = 0;
                t2 = 0;

                BEGIN_ITER_IMAGE
                    // for U
                    x += Z[i][j][k] * Y[k][i][j][m] * u[k][i][j][m] * im_ad[i][j];
                    y += Z[i][j][k] * Y[k][i][j][m] * u[k][i][j][m];

                    // for V
//                    t1 += Z[i][j][k] * Y[k][i][j][m] * (log(u[k][i][j][m]) - u[k][i][j][m]);
//                    t2 += Z[i][j][k] * Y[k][i][j][m];

                    t1 += Z[i][j][k] * Y[k][i][j][m];
                    t3 = 0;
                    for ( n = 0; n < Kj[k]; n++ ) {
                        t3 += Y[k][i][j][n];
                    }
                    t2 += Z[i][j][k] * t3;

                    if ( m == 0 ) {
                        // image blurring submatrix sum
                        if (i == 0) {
                            temp[i][j] = Z[i][j][k]; // Copy first row
                        } else {
                            temp[i][j] = Z[i][j][k] + temp[i-1][j]; // sum columns
                        }
                    }
                END_ITER

                U[k][m] = x / y;
                E[k][m] = t1 / t2;
            }

            // image blurring submatrix sum rows
            BEGIN_ITER_IMAGE
                if (j > 0) {
                    temp[i][j] += temp[i][j-1];
                }
            END_ITER

            // Blurring using submatrix sum
            BEGIN_ITER_IMAGE
                a = i - 3;
                b = i + 3;
                c = j - 3;
                d = j + 3;

                a = a > M ? 0 : a; // Remember, it's an unsigned int -- tli
                b = b >= M ? M-1 : b; // rbi
                c = c > N ? 0 : c; // tlj
                d = d >= N ? N-1 : d; // rbj

                AveLocZ[i][j][k] = temp[b][d];
                if (a > 0) AveLocZ[i][j][k] -= temp[a-1][d];
                if (c > 0) AveLocZ[i][j][k] -= temp[b][c-1];
                if (a > 0 && c > 0) AveLocZ[i][j][k] += temp[a-1][c-1];

                AveLocZ[i][j][k] /= 49;
//                MP[i][j][k] = fexp(B * AveLocZ[i][j][k]);
                MP[i][j][k] = fmath::expd(B * AveLocZ[i][j][k]);
                sumMP[i][j] += MP[i][j][k];
            END_ITER
        }

        for ( k = 0; k < K; k++ ) {
            BEGIN_ITER_IMAGE
                MP[i][j][k] /= sumMP[i][j];
            END_ITER

            for ( m = 0; m < Kj[k]; m++ ) {
                if ( S[k][m] > 500 ) {
                    t1 = 0;
                    t2 = 0;

                    BEGIN_ITER_IMAGE
                        t3 = im_ad[i][j] - U[k][m];
                        t1 += Z[i][j][k] * Y[k][i][j][m] * u[k][i][j][m] * t3 * t3;
                        t2 += Z[i][j][k] * Y[k][i][j][m];
                    END_ITER
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

        while(fabs(B - Bold) > 0.05) {
            Bold = B;
            t4 = 0;

            BEGIN_ITER_IMAGE
                t1 = 0;
                t2 = 0;

                for ( k = 0; k < K; k++ ) {
//                    t3 = fexp(B * AveLocZ[i][j][k]);
                    t3 = fmath::expd(B * AveLocZ[i][j][k]);
                    t1 += t3 * AveLocZ[i][j][k];
                    t2 += t3;
                }

                t3 = t1 / t2;
                t1 = 0;

                for ( k = 0; k < K; k++ ) {
                    t1 += (AveLocZ[i][j][k] - t3) * Z[i][j][k];
                }

                t4 += -t1;
            END_ITER

            B = Bold - t4 * learning_rate;
        }

        LogLike[iter] = 0;
        BEGIN_ITER_IMAGE
            t1 = 0;
            for ( k = 0; k < K; k++ ) {
                t2 = 0;
                for ( m = 0; m < Kj[k]; m++ ) {
                    t2 += E[k][m] * StudPDFVal[k][i][j][m];;
                }

                t1 += MP[i][j][k] * t2;
            }

//            LogLike[iter] += log(t1);
            LogLike[iter] += fmath::log(t1);
        END_ITER

        Ptot = LogLike[iter];
        int stop_s = clock();

        if (iter > 0) {
            // printf("Iterations = %d LogLikelihood = %f DiffLikelihood = %f Time = %fms\n", iter, LogLike[iter], fabs(LogLike[iter - 1] - LogLike[iter]), (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);

            if (fabs(LogLike[iter - 1] - LogLike[iter]) < 3000) break;
        }

        iter++;
    }

    BEGIN_ITER_IMAGE
        segmented[i * N + j] = 0;
        t1 = Z[i][j][0];
        for ( k = 1; k < K; k++ ) {
            if (Z[i][j][k] > t1) {
                segmented[i * N + j] = k;
                t1 = Z[i][j][k];
                // printf("SEGMENTED %lu %lu %lu\n", i, j, k);
            }
        }
    END_ITER

    return Ptot;
}

int main() {
    auto image = read();
    vector<unsigned char> segmented(get<0>(image).size(), 0);
//    vector<double> debugarr(get<0>(image).size() * 3, 0);

    cmm(get<0>(image).data(), segmented.data(), get<1>(image), get<2>(image), 59);

    return 0;
}
