//============================================================================80
// Simulate TGF-beta expression at the single-cell level.
//
// Author: Keita Iida (Haga, Iida, Okada, 2023)
//============================================================================80
#include <stdio.h>
#include <stdlib.h>
#include "param.h"
#include "pg_001_my.h"
//----------------------------------------------------------------------------80
// output
//----------------------------------------------------------------------------80
void output( double p, long int t, double x[], double y[], double z[] )
{
  int   k;
  char  filename[100];
  FILE  *fp;

  sprintf( filename, "%s_THBS1_p_%3.2f.txt", FILE_CELLS, p );
  fp = fopen( filename, "a" );
  fprintf( fp, "%ld ", t );
  for( k=0; k<NCELL; ++k )
    fprintf( fp, "%lf ", x[k] );
  fprintf( fp, "\n" );
  fclose( fp );

  sprintf( filename, "%s_FMOD_p_%3.2f.txt", FILE_CELLS, p );
  fp = fopen( filename, "a" );
  fprintf( fp, "%ld ", t );
  for( k=0; k<NCELL; ++k )
    fprintf( fp, "%lf ", y[k] );
  fprintf( fp, "\n" );
  fclose( fp );

  sprintf( filename, "%s_TGFb1_p_%3.2f.txt", FILE_CELLS, p );
  fp = fopen( filename, "a" );
  fprintf( fp, "%ld ", t );
  for( k=0; k<NCELL; ++k )
    fprintf( fp, "%lf ", z[k] );
  fprintf( fp, "\n" );
  fclose( fp );
}
