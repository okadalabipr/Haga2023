//============================================================================80
// Compute steady state solutions for THBS1, FMOD, and TGFb1 expressions, wrt
// parameters for PDL.
//
// Author: Keita Iida (Haga, Iida, Okada, 2023)
//============================================================================80
#include <stdio.h>
#include <math.h>
#include "param.h"
#include "pg_001_my.h"
//----------------------------------------------------------------------------80
// f
//----------------------------------------------------------------------------80
double f( double x, double p, double q, double r )
{
  double  xna, xnb, KA, KB, CC;

  xna = pow( x, na );
  xnb = pow( x, nb );
  KA  = p / A;
  KB  = q / B;
  CC  = r;

  return CC * xna / (xna + KA * (xna + Ka))
    * KB * (xnb + Kb) / (Kb + KB * (xnb + Kb)) - x;
}
//----------------------------------------------------------------------------80
// dfdx
//----------------------------------------------------------------------------80
double dfdx( double x, double p, double q, double r )
{
  double  xna, xnb, KA, KB, CC, term1, term2, term3, term4;

  xna = pow( x, na );
  xnb = pow( x, nb );
  KA  = p / A;
  KB  = q / B;
  CC  = r;

  term1 = na * KA * Ka * pow( x, na - 1 ) / pow( xna + KA * (xna + Ka), 2 );
  term2 = KB * (xnb + Kb) / (Kb + KB * (xnb + Kb));
  term3 = xna / (xna + KA * (xna + Ka));
  term4 = nb * KB * Kb * pow( x, nb - 1 ) / pow( Kb + KB * (xnb + Kb), 2 );

  return CC * (term1 * term2 + term3 * term4) - 1.0;
}
//----------------------------------------------------------------------------80
// x1
//----------------------------------------------------------------------------80
double x1( double x, double p, double q, double r )
{
  double  xna = pow( x, na );

  return A * xna / (xna + Ka);
}
//----------------------------------------------------------------------------80
// x2
//----------------------------------------------------------------------------80
double x2( double x, double p, double q, double r )
{
  double  xnb = pow( x, nb );

  return B * Kb / (xnb + Kb);
}

