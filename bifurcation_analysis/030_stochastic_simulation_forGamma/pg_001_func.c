//============================================================================80
// Simulate TGF-beta expression at the single-cell level.
//
// Author: Keita Iida (Haga, Iida, Okada, 2023)
//============================================================================80
#include <stdio.h>
#include <math.h>
#include "param.h"
#include "pg_001_my.h"
//----------------------------------------------------------------------------80
// f_tgfb1
//----------------------------------------------------------------------------80
double f_tgfb1( double z, double p )
{
  double  zna, znb, K, L;

  zna = pow( z, NA );
  znb = pow( z, NB );
  K = K1 / (A / D1);
  L = K2 / (B / D2);

  return p * zna / (zna + K * (zna + KA))
    * L * (znb + KB) / (KB + L * (znb + KB));
}
//----------------------------------------------------------------------------80
// f_thbs1
//----------------------------------------------------------------------------80
double f_thbs1( double x, double z )
{
  double  zna;

  zna = pow( z, NA );

  return (A / D1) * zna / (zna + KA);
}
//----------------------------------------------------------------------------80
// f_fmod
//----------------------------------------------------------------------------80
double f_fmod( double y, double z )
{
  double  znb;

  znb = pow( z, NB );

  return (B / D2) * KB / (znb + KB);
}
//----------------------------------------------------------------------------80
// Rand_Exp: random sampling from exponential distribution
//----------------------------------------------------------------------------80
double Rand_Exp( double beta )
{
  double  rnd;

  rnd = genrand64_real2();  //[0,1)-real value
  return -beta * log( 1.0 - rnd );
}
//----------------------------------------------------------------------------80
// Rand_Norm: random sampling from normal distribution
//----------------------------------------------------------------------------80
double Rand_Norm( double mu, double sigma )
{
  double  rnd, rnd_1, rnd_2;

  rnd_1 = genrand64_real2();  //[0,1)-real value
  rnd_2 = genrand64_real2();  //[0,1)-real value
  rnd = sqrt( -2.0 * log( rnd_1 ) ) * sin( 2.0 * M_PI * rnd_2 );
  return mu + sigma * rnd;
}
