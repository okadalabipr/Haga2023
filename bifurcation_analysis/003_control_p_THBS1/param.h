//============================================================================80
// Compute steady state solutions for THBS1, FMOD, and TGFb1 expressions, wrt
// parameters for THBS1 and FMOD, i.e., K1 (p) and K2 (q), respectively.
//
// Author: Keita Iida (Haga, Iida, Okada, 2023)
//============================================================================80
//----------------------------------------------------------------------------80
// File name
//----------------------------------------------------------------------------80
#define FILENAME  "01_outgoing/x.txt"

//----------------------------------------------------------------------------80
// Control parameters
//----------------------------------------------------------------------------80
#define PMIN      10.0        // Parameter range for THBS1 (K1)
#define PMAX      20.0        // Parameter range for THBS1 (K1)
#define QMIN      1.0         // Parameter range for FMOD (K2)
#define QMAX      1.0         // Parameter range for FMOD (K2)
//----------------------------------------------------------------------------80
// Parameters for Newton's method
//----------------------------------------------------------------------------80
#define EPS       0.00000001  // Tolerance for convergence
#define ITEMAX    20          // Maximal iteration step
#define X0MIN     0.00001     // Range for initial guess
#define X0MAX     1.0         // Range for initial guess
//----------------------------------------------------------------------------80
// Model parameters (nondimentionalized)
//----------------------------------------------------------------------------80
#define A        2.36
#define B        0.33
#define C        1.0
#define Ka       0.016
#define Kb       0.002
#define na       1.6
#define nb       1.7

