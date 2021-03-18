#include <ctime>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>

#ifdef __x86_64__
	#include "fmath.hpp"
	#define mm_gamma(x) fmath::fastgamma3(x)
	#define mm_exp(x) fmath::expd(x)
#else
	#define mm_gamma(x) tgamma(x)
	#define mm_exp(x) exp(x)
#endif

using namespace std;

/*
 *  Adapted from https://stackoverflow.com/a/7455490
 *  cmath::isnan has problems with -Ofast
 */
static inline int safeisnan(float f) {
    union { float f; uint32_t x; } u = { f };
    return (u.x << 1) > 0xff000000u;
}

typedef void (*cbfunc)(int step, double *x, int m, int n, int o, void *user_data);

struct CmmConfig {
    double learning_rate = 1e-6;
    double var = 1e8;
    int m_nuc_size = 210;
    int m_cell_size = 1900;
    char max_iterations = 100;
    bool debug = false;
};

void cmm(const unsigned short *im, unsigned char *segmented, unsigned int dim1, unsigned int dim2, int nuc_num,
         CmmConfig config, cbfunc callback, void *user_data);
