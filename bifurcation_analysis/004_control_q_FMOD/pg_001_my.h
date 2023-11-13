//============================================================================80
// Compute steady state solutions for THBS1, FMOD, and TGFb1 expressions, wrt
// parameters for THBS1 and FMOD, i.e., K1 (p) and K2 (q), respectively.
//
// Author: Keita Iida (Haga, Iida, Okada, 2023)
//============================================================================80
#include <stdio.h>
#include "param.h"
//----------------------------------------------------------------------------80
// pg_001_func.c
//----------------------------------------------------------------------------80
double f( double x, double p, double q );
double dfdx( double x, double p, double q );
double g( double x, double p, double q );
double dgdx( double x, double p, double q );
double x1( double x, double p, double q );
double x2( double x, double p, double q );
//----------------------------------------------------------------------------80
// pg_001_outg.c
//----------------------------------------------------------------------------80
void output0( void );
void output( double p, double q, double x );
