#ifndef SCATTER
#define SCATTER


#include "event1.h"


#include <array>


namespace TwoThreeScatter {
  
enum Four_mom_type { 
  p1,   // incoming nucleon
  p2,   // incoming neutrino
  q1,   // outgoing nucleon
  q2,   // outgoing kaon
  q3,   // outgoing lepton
  q,    // 4-momentum transfer = p2 - q3 = q1 - p1
  max_types
};
 
double kallen(double x, double y, double z);

// 4-momentum components
double Ep1(double s, double ma);
double Ep2(double s, double ma); 
double E3(double s, double m3, double W); 
double q3mod(double s, double m3, double W);
double q12CM(double m1, double m2, double W);
double E12(double s, double m3, double W);
double gamma(double s, double m3, double W); 
double beta(double s, double m3, double W); 
double q10Star(double m1, double m2, double W);
double q20Star(double m1, double m2, double W);
double q11Star(double m1, double m2, double W, double thetaStar, double phiStar);
double q21Star(double m1, double m2, double W, double thetaStar, double phiStar);
double q12Star(double m1, double m2, double W, double thetaStar, double phiStar);
double q22Star(double m1, double m2, double W, double thetaStar, double phiStar);
double q13Star(double m1, double m2, double W, double thetaStar);
double q23Star(double m1, double m2, double W, double thetaStar);
double q10(double m1, double m2, double m3, double s, double W, double thetaStar);
double q20(double m1, double m2, double m3, double s, double W, double thetaStar);
double q11(double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar);
double q21(double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar);
double q12(double m1, double m2, double W, double thetaStar, double phiStar); 
double q22(double m1, double m2, double W, double thetaStar, double phiStar);
double q13(double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar);
double q23(double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar);

// mathematica wrapper functions
vect Momentum(Four_mom_type k, double ma, double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar);
/*
  This function translates the Momentum() fuction outputted by Cform to a (4-)vect
*/

double Pair(vect k1, vect k2);
/*
  This function translates the Pair() fuction outputted by Cform to a dot-product
*/

double Power(double x, int n);
/*
  This function translates the Power() fuction outputted by Cform to a std::pow
*/

// 4D Levi-Civita functionality
using LeviCivitaTable = std::array<std::array<std::array<std::array<int, 4>, 4>, 4>, 4>;

constexpr LeviCivitaTable generateTable() { // C++17 & earlier
  LeviCivitaTable lc{}; 

  lc[0][1][2][3] =  1; lc[0][1][3][2] = -1;
  lc[0][2][3][1] =  1; lc[0][2][1][3] = -1;
  lc[0][3][1][2] =  1; lc[0][3][2][1] = -1;

  lc[1][3][2][0] =  1; lc[1][3][0][2] = -1;
  lc[1][2][0][3] =  1; lc[1][2][3][0] = -1;
  lc[1][0][3][2] =  1; lc[1][0][2][3] = -1;

  lc[2][0][1][3] =  1; lc[2][0][3][1] = -1;
  lc[2][1][3][0] =  1; lc[2][1][0][3] = -1;
  lc[2][3][0][1] =  1; lc[2][3][1][0] = -1;

  lc[3][2][1][0] =  1; lc[3][2][0][1] = -1;
  lc[3][0][2][1] =  1; lc[3][0][1][2] = -1;
  lc[3][1][0][2] =  1; lc[3][1][2][0] = -1;

  return lc;
}

constexpr LeviCivitaTable table = generateTable();

constexpr bool hasDuplicateIndices(int i, int j, int k, int l) {
  return (i == j) || (i == k) || (i == l) || (j == k) || (j == l) || (k == l);
}

constexpr int getLeviCivitaValue(int i, int j, int k, int l) {
  return hasDuplicateIndices(i, j, k, l) ? 0 : table[i][j][k][l];
}

double Eps(vect k1, vect k2, vect k3, vect k4);
/*
  This function translates the Eps() fuction outputted by Cform to a contraction of the levi-civita with four (4-)vect
*/

}

#endif
