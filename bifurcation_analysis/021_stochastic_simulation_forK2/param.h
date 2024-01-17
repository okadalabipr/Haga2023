//============================================================================80
// Simulate TGF-beta expression at the single-cell level.
//
// Author: Keita Iida (Haga, Iida, Okada, 2023)
//============================================================================80
#define TMAX            (100 * 10000)
#define NCELL           10000
#define DT              0.01
#define	CUT             10000

#define FILE_CELLS      "01_outgoing/cells"

#define PMIN            0.0
#define PMAX            2.0
#define DP              0.5

#define A               2.36
#define B               0.33
#define C               1.0
#define D1              1.0
#define D2              1.0
#define D3              1.0
#define KA              0.016
#define KB              0.002
#define NA              1.6
#define NB              1.7
#define K1              16.0
#define K2              1.0

#define ZMAX            0.2    // Initial [TGFb1] is sampled from U(0, ZMAX).
#define MU              0.0    // Mean for Gaussian white nose
#define SIGMA           0.2    // Standard deviation for Gaussian white noise

