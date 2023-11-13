//============================================================================80
// Compute steady state solutions for THBS1, FMOD, and TGFb1 expressions, wrt
// parameters for THBS1 and FMOD, i.e., K1 (p) and K2 (q), respectively.
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
double f( double x, double p, double q )
{
  double  xna, xnb, KA, KB;

  xna = pow( x, na );
  xnb = pow( x, nb );
  KA  = p / A;
  KB  = q / B;

  return C * xna / (xna + KA * (xna + Ka))
    * KB * (xnb + Kb) / (Kb + KB * (xnb + Kb)) - x;
}
//----------------------------------------------------------------------------80
// dfdx
//----------------------------------------------------------------------------80
double dfdx( double x, double p, double q )
{
  double  xna, xnb, KA, KB, term1, term2, term3, term4;

  xna = pow( x, na );
  xnb = pow( x, nb );
  KA  = p / A;
  KB  = q / B;

  term1 = na * KA * Ka * pow( x, na - 1 ) / pow( xna + KA * (xna + Ka), 2 );
  term2 = KB * (xnb + Kb) / (Kb + KB * (xnb + Kb));
  term3 = xna / (xna + KA * (xna + Ka));
  term4 = nb * KB * Kb * pow( x, nb - 1 ) / pow( Kb + KB * (xnb + Kb), 2 );

  return C * (term1 * term2 + term3 * term4) - 1.0;
}
//----------------------------------------------------------------------------80
// g
//----------------------------------------------------------------------------80
double g( double x, double p, double q )
{
  double  xna, xnb, KA, KB;

  xna = pow( x, na );
  xnb = pow( x, nb );
  KA  = p / A;
  KB  = q / B;

  return C * KB * pow( x, na - 1 ) * (xnb + Kb)
    - (xna + KA * (xna + Ka)) * (Kb + KB * (xnb + Kb));
}
//----------------------------------------------------------------------------80
// dgdx
//----------------------------------------------------------------------------80
double dgdx( double x, double p, double q )
{
  double  xna, xnb, KA, KB, term1, term2, term3, term4;

  xna = pow( x, na );
  xnb = pow( x, nb );
  KA  = p / A;
  KB  = q / B;

  term1 = (na - 1) * C * KB * pow( x, na - 2 ) * (xnb + Kb);
  term2 = nb * C * KB * pow( x, na - 1 ) * pow( x, nb - 1 );
  term3 = na * (1 + KA) * pow( x, na - 1 ) * (Kb + KB * (xnb + Kb));
  term4 = (xna + KA * (xna + Ka)) * nb * KB * pow( x, nb - 1 );

  return  term1 + term2 - term3 - term4;
}
//----------------------------------------------------------------------------80
// x1
//----------------------------------------------------------------------------80
double x1( double x, double p, double q )
{
  double  xna = pow( x, na );

  return A * xna / (xna + Ka);
}
//----------------------------------------------------------------------------80
// x2
//----------------------------------------------------------------------------80
double x2( double x, double p, double q )
{
  double  xnb = pow( x, nb );

  return B * Kb / (xnb + Kb);
}

