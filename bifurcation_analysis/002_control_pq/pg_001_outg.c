//============================================================================80
// Compute steady state solutions for THBS1, FMOD, and TGFb1 expressions, wrt
// parameters for THBS1 and FMOD, i.e., K1 (p) and K2 (q), respectively.
//
// M. Haga, K. Iida, M. Okada (2023)
//============================================================================80
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "param.h"
#include "pg_001_my.h"
//----------------------------------------------------------------------------80
// output0
//----------------------------------------------------------------------------80
void output0( void )
{
  FILE  *fp;

  fp = fopen( FILENAME, "w" );
  fprintf( fp, "#p q x1 x2 x3 df/dx3\n" );
  fclose( fp );
}
//----------------------------------------------------------------------------80
// output
//----------------------------------------------------------------------------80
void output( double p, double q, double x )
{
  FILE  *fp;

  fp = fopen( FILENAME, "a" );
  fprintf( fp, "%lf %lf %lf %lf %lf %lf\n",
           p, q, x1(x, p, q), x2(x, p, q), x, dfdx( x, p, q ) );
  fclose( fp );
}
