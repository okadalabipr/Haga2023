//============================================================================80
// Compute steady state solutions for THBS1, FMOD, and TGFb1 expressions, wrt
// parameters for PDL.
//
// Author: Keita Iida (Haga, Iida, Okada, 2023)
//============================================================================80
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "param.h"
#include "pg_001_my.h"
//----------------------------------------------------------------------------80
// main
//----------------------------------------------------------------------------80
int main()
{
  int     ite;                  // Number of iteration of Newton's method
  double  p, q, r, dp, dq, dr;  // Control parameters and the increments
  double  x0, dx0;              // Initial guess and the increment
                                // for Newton's method
  double  x, new_x;             // Numerical solutions
  double  dydx, dx;             // Temporal variables for Newtons' method
  //------------------------------------------------------60
  // Initialize
  //------------------------------------------------------60
  output0();

  //------------------------------------------------------60
  // Loop for p
  //------------------------------------------------------60
  p = PMIN;
  while( p <= PMAX ){
    dp = 0.01;
    //----------------------------------------------------60
    // Loop for q
    //----------------------------------------------------60
    q = QMIN;
    while( q <= QMAX ){
      dq = 0.001;
      //--------------------------------------------------60
      // Loop for r
      //--------------------------------------------------60
      r = RMIN;
      while( r <= RMAX ){
        if( (r >= 0.074) & (r < 0.0745) ){
          dr = 0.00001;
        }else if( (r >= 0.06) & (r < 0.1) ){
          dr = 0.0001;
        }else{
          dr = 0.001;
        }
        //------------------------------------------------60
        // Loop for x0
        //------------------------------------------------60
        x0 = X0MIN;
        while( x0 <= X0MAX ){
          dx0 = 0.01;
          x = x0;
          //----------------------------------------------60
          // Newton's method
          //----------------------------------------------60
          ite = 0;
          while( 1 ){
            dydx = dfdx( x, p, q, r );
            if( fabs( dydx ) < EPS ){
              break;
            }
            dx = f( x, p, q, r ) / dydx;
            new_x = x - dx;
            if( fabs( dx ) < EPS ){
              output( p, q, r, x );
              break;
            }
            x = new_x;
            ite++;
            if( ite > ITEMAX ){
              break;
            }
          }
          x0 += dx0;
        }
        r += dr;
      }
      q += dq;
    }
    p += dp;
  }

  return 0;
}
