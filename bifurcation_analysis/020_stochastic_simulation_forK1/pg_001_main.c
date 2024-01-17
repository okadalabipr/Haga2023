//============================================================================80
// Simulate TGF-beta expression at the single-cell level.
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
  long int  t;
  int       k;
  double    x[NCELL], new_x[NCELL],  // THBS1
            y[NCELL], new_y[NCELL],  // FMOD
            z[NCELL], new_z[NCELL],  // TGFb1
            dN, p;
  //--------------------------------------------------------------------------80
  // Loop for p
  //--------------------------------------------------------------------------80
  p = PMIN;
  while( p <= PMAX ){
    // Initialize x
    for( k=0; k<NCELL; ++k ){
      x[k] = 0.0;
      y[k] = 0.0;
      z[k] = ZMAX * genrand64_real1();  //[0,1]-real value
    }
    //----------------------------------------------------60
    // Loop for t
    //----------------------------------------------------60
    for( t=0; t<=TMAX; ++t ){
      if( t % CUT == 0 ){
        output( p, t, x, y, z );
        printf( "t / CUT = %ld\n", t / CUT );
      }
      //------------------------------40
      // Loop for k
      //------------------------------40
      for( k=0; k<NCELL; ++k ){
        dN = sqrt( DT ) * Rand_Norm( 0, 1 );
        new_x[k] = x[k] + DT * (f_thbs1( x[k], z[k] ) - x[k]) + SIGMA * dN;
        new_y[k] = y[k] + DT * (f_fmod( y[k], z[k] ) - y[k]) + SIGMA * dN;
        new_z[k] = z[k] + DT * (f_tgfb1( z[k], p ) - z[k]);
      }
      // Update
      for( k=0; k<NCELL; ++k ){
        if( new_x[k] < 0.0 ){
          new_x[k] = 0.0;
        }
        if( new_y[k] < 0.0 ){
          new_y[k] = 0.0;
        }
        if( new_z[k] < 0.0 ){
          new_z[k] = 0.0;
        }
        x[k] = new_x[k];
        y[k] = new_y[k];
        z[k] = new_z[k];
      }
    }
    p = p + DP;
  }
}
