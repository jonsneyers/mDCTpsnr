/************************************************************************************
 **  Copyright (C) 2009 University of Stuttgart, Thomas Richter                    **
 **                                                                                **
 **  This software is provided 'as-is', without any express or implied             **
 **  warranty.  In no event will the authors be held liable for any damages        **
 **  arising from the use of this software.                                        **
 **                                                                                **
 **  Permission is granted to anyone to use this software for any purpose,         **
 **  including commercial applications, and to alter it and redistribute it        **
 **  freely, subject to the following restrictions:                                **
 **                                                                                **
 **  1. The origin of this software must not be misrepresented; you must not       **
 **     claim that you wrote the original software. If you use this software       **
 **     in a product, an acknowledgment in the product documentation would be      **
 **     appreciated but is not required.                                           **
 **  2. Altered source versions must be plainly marked as such, and must not be    **
 **     misrepresented as being the original software.                             **
 **  3. This notice may not be removed or altered from any source distribution.    **
 **                                                                                **
 **     Thomas Richter                                                             **
 **     thomas.richter@iis.fraunhofer.de                                           **
 **                                                                                **
 ************************************************************************************/

/// Includes
#include "ediff.hpp"
#include "std/math.hpp"


/// Gamma
// Implement the nonlinearity of the CIE Lab color space
static double Gamma(double t)
{ 
  if (t > 0.008856)
    return cbrt(t);
  else
    return 7.787037 * t + 0.137931;
}
///

/// OppToCIE
// Convert the opponent color space back to CIELab
static void OppToCIE(double &l,double &as,double &bs,
		     double  a,double  c1,double c2)
{
  double x  =  0.9795961601 * a - 1.534715701  * c1 + 0.4445976432 * c2;
  double y  =  1.188977906  * a + 0.7643549573 * c1 + 0.1351257479 * c2;
  double z  =  1.231833313  * a + 1.163159260  * c1 + 2.078407588  * c2;
  //
  // D65 whitepoint in XYZ
  double yn = 100.0;
  double xn =  95.4;
  double zn = 108.88;
  //
  // Normalize to 100 and to the white point.
  x *= 100.0 / xn;
  y *= 100.0 / yn;
  z *= 100.0 / zn;
  //
  // Convert now to CIELab
  l  = 116.0 * Gamma(y) - 16.0;
  as = 500.0 * (Gamma(x) - Gamma(y));
  bs = 200.0 * (Gamma(y) - Gamma(z));
}
///

/// DeltaE
double DeltaE(double a1ref,double c1ref,double c2ref,
	      double a1dst,double c1dst,double c2dst)
{
  // Step 1: Convert the opponent color space to CIELab
  double a1,b1,l1;
  double a2,b2,l2;
  //
  OppToCIE(l1,a1,b1,a1ref,c1ref,c2ref);
  OppToCIE(l2,a2,b2,a1dst,c1dst,c2dst);
  //
  double c1 = sqrt(a1*a1 + b1*b1);
  double c2 = sqrt(a2*a2 + b2*b2);
  double cb = 0.5 * (c1 + c2);
  double g  = 0.5 * (1.0 - sqrt(pow(cb,7) / (pow(cb,7) + pow(25,7))));
  // Fixup the value of a, b and L remain unaffected.
  a1        = (1 + g) * a1;
  a2        = (1 + g) * a2;
  // Recompute chroma and hue angle.
  c1        = sqrt(a1*a1 + b1*b1);
  c2        = sqrt(a2*a2 + b2*b2);
  double h1 = atan2(b1,a1);
  double h2 = atan2(b2,a2);
  // CIEDE2000 uses the convention that the cut is at the positive real axis. Urgh.
  if (h1 < 0)
    h1     += 2 * M_PI;
  if (h2 < 0)
    h2     += 2 * M_PI;
  double dl = l1 - l2;
  double dc = c1 - c2;
  double dh = h1 - h2;
  // Fixup for difference larger than 180 degrees.
  if (dh > M_PI) {
    dh -= 2.0*M_PI;
  } else if (dh < -M_PI) {
    dh += 2.0*M_PI;
  }
  double dH = 2.0 * sqrt(c1 * c2) * sin(dh * 0.5);
  // Compute mean lightness, chroma and hue.
  double lm = 0.5 * (l1 + l2);
  double cm = 0.5 * (c1 + c2);
  // Compute mean angle (intermediate result)
  double hm = 0.5 * (h1 + h2); 
  // Again, fixup mean for large angle differences.
  dh        = h1 - h2;
  if (dh < -M_PI || dh > M_PI) {
    if (hm < M_PI) {
      hm   += M_PI;
    } else {
      hm   -= M_PI;
    }
  }
  // Compute sl,sc. First normalize lb for convenience.
  double lb = lm - 50.0;
  double sl = 1.0 + ((0.015 * lb * lb) / sqrt(20.0 + lb * lb));
  double sc = 1.0 + 0.045 * cm;
  // Hue weighting. Tough.
  double t  = 1.0 - 0.17 * cos(hm - 30.0 / 180.0 * M_PI) + 0.24 * cos(2 * hm) + 
    0.32*cos(3 * hm + 6.0 / 180.0 * M_PI) - 0.20 * cos(4 * hm - 63.0 / 180.0 * M_PI); 
  double sh = 1.0 + 0.015 * cm * t;
  // Correction for blue nonlinearities.
  double rc = 2.0 * sqrt(pow(cm,7) / (pow(cm,7) + pow(25,7)));
  double dt = (hm * 180.0 / M_PI - 275.0) / 25.0; // difference to blue "eccentricy"
  dt        = 30.0 * exp(-dt * dt) * M_PI / 180.0;
  double rt = -sin(2 * dt) * rc;
  // Now compute delta E, the final error difference output.
  double de = sqrt((dl * dl)/(sl * sl) + (dc * dc)/(sc * sc) + (dH * dH)/(sh * sh) + rt * dc * dH / (sc * sh));
  //
  return de;
}
///


