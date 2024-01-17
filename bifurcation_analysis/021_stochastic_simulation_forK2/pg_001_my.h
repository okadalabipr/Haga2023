//============================================================================81
// Simulate TGF-beta expression at the single-cell level.
//
// Author: Keita Iida (Haga, Iida, Okada, 2023)
//============================================================================80
#include <stdio.h>
#include "param.h"
//----------------------------------------------------------------------------80
// math.sci.hiroshima-u.ac.jp_m-mat_MT_VERSIONS_C-LANG_mt19937-64.c
//----------------------------------------------------------------------------80
void init_genrand64(unsigned long long seed);
void init_by_array64(unsigned long long init_key[],
                     unsigned long long key_length);
unsigned long long genrand64_int64(void);
long long genrand64_int63(void);
double genrand64_real1(void);
double genrand64_real2(void);
double genrand64_real3(void);
//----------------------------------------------------------------------------80
// pg_001_func.c
//----------------------------------------------------------------------------80
double f_tgfb1( double z, double p );
double f_thbs1( double x, double z );
double f_fmod( double y, double z );
double Rand_Exp( double beta );
double Rand_Norm( double mu, double sigma );
//----------------------------------------------------------------------------80
// pg_001_outg.c
//----------------------------------------------------------------------------80
void output( double p, long int t, double x[], double y[], double z[] );
