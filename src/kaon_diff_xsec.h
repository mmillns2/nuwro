#ifndef KAON_DIFF_XSEC
#define KAON_DIFF_XSEC


#include "vect.h"


double single_kaon_integrand(double ma, double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar);
/*
  
*/


double simpson2D(double thetaStar_min, double thetaStar_max, double phiStar_min, double phiStar_max, int thetaStar_n, 
                 int phiStar_n, double ma, double m1, double m2, double m3, double s, double W, double theta);
/*

*/


double single_kaon_diff_xsec(double ma, double m1, double m2, double m3, double s, double W, double theta);
/*

*/


double single_kaon_diff_xsec_2(double s, double N1_mass, double kaon_mass, double lepton_mass, vect vects[6]);
/*

*/


#endif
