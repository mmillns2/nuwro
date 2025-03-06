#include "23Scatter.h"



namespace TwoThreeScatter {

double kallen(double x, double y, double z)
{
	double ret{ std::pow(x, 4) + std::pow(y, 4) + std::pow(z, 4) - 2*std::pow(x, 2)*std::pow(y, 2) 
				 - 2*std::pow(x, 2)*std::pow(z, 2) - 2*std::pow(y, 2)*std::pow(z, 2) }; 	
	if(fabs(ret) < 1e-5)
		    ret = 0;
	assert(ret >= 0);
	return std::sqrt(ret);
}

double Ep1(double s, double ma) {
	return (s + ma*ma) / (2*std::sqrt(s));
}

double Ep2(double s, double ma) {
	return std::sqrt(s) - Ep1(s, ma);
}

double E3(double s, double m3, double W) {
  return (s + m3*m3 - W*W) / (2*std::sqrt(s));
}

double q3mod(double s, double m3, double W) {
  return kallen(std::sqrt(s), W, m3) / (2*std::sqrt(s));
}

double q12CM(double m1, double m2, double W) {
	return kallen(W, m1, m2) / (2*W);
}

double E12(double s, double m3, double W) {
  return (s + W*W - m3*m3) / (2*std::sqrt(s));
}

double gamma(double s, double m3, double W) {
  double ret{ E12(s, m3, W) / W };
  if(fabs(1 - ret) < 0.0000001) ret = 1;
  return ret;
}

double beta(double s, double m3, double W) {
  return std::sqrt(1 - 1 / (gamma(s, m3, W) * gamma(s, m3, W)));
}

double q10Star(double m1, double m2, double W) {
  return (W*W + m1*m1 - m2*m2) / (2*W);
}

double q20Star(double m1, double m2, double W) {
  return (W*W + m2*m2 - m1*m1) / (2*W);
}

double q11Star(double m1, double m2, double W, double thetaStar, double phiStar) {
  return q12CM(m1, m2, W) * std::sin(thetaStar) * std::cos(phiStar);
}

double q21Star(double m1, double m2, double W, double thetaStar, double phiStar) {
  return -q12CM(m1, m2, W) * std::sin(thetaStar) * std::cos(phiStar);
}

double q12Star(double m1, double m2, double W, double thetaStar, double phiStar) {
  return q12CM(m1, m2, W) * std::sin(thetaStar) * std::sin(phiStar);
}

double q22Star(double m1, double m2, double W, double thetaStar, double phiStar) {
  return -q12CM(m1, m2, W) * std::sin(thetaStar) * std::sin(phiStar);
}

double q13Star(double m1, double m2, double W, double thetaStar) {
  return q12CM(m1, m2, W) * std::cos(thetaStar);
}

double q23Star(double m1, double m2, double W, double thetaStar) {
  return -q12CM(m1, m2, W) * std::cos(thetaStar);
}

double q10(double m1, double m2, double m3, double s, double W, double thetaStar) {
  return gamma(s, m3, W) * (q10Star(m1, m2, W) + beta(s, m3, W) * q13Star(m1, m2, W, thetaStar));
}

double q20(double m1, double m2, double m3, double s, double W, double thetaStar) {
  return gamma(s, m3, W) * (q20Star(m1, m2, W) + beta(s, m3, W) * q23Star(m1, m2, W, thetaStar));
}

double q11(double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar) {
  return q11Star(m1, m2, W, thetaStar, phiStar) * std::cos(theta + M_PI) + 
         gamma(s, m3, W) * (beta(s, m3, W) * q10Star(m1, m2, W) + q13Star(m1, m2, W, thetaStar)) * std::sin(theta + M_PI);
}

double q21(double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar) {
  return q21Star(m1, m2, W, thetaStar, phiStar) * std::cos(theta + M_PI) + 
         gamma(s, m3, W) * (beta(s, m3, W) * q20Star(m1, m2, W) + q23Star(m1, m2, W, thetaStar)) * std::sin(theta + M_PI);
}

double q12(double m1, double m2, double W, double thetaStar, double phiStar) {
  return q12Star(m1, m2, W, thetaStar, phiStar);
}

double q22(double m1, double m2, double W, double thetaStar, double phiStar) {
  return q22Star(m1, m2, W, thetaStar, phiStar);
}

double q13(double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar) {
  return -q11Star(m1, m2, W, thetaStar, phiStar) * std::sin(theta + M_PI) + 
         gamma(s, m3, W) * (beta(s, m3, W) * q10Star(m1, m2, W) + q13Star(m1, m2, W, thetaStar)) * std::cos(theta + M_PI);
}

double q23(double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar) {
  return -q21Star(m1, m2, W, thetaStar, phiStar) * std::sin(theta + M_PI) + 
         gamma(s, m3, W) * (beta(s, m3, W) * q20Star(m1, m2, W) + q23Star(m1, m2, W, thetaStar)) * std::cos(theta + M_PI);
}

vect Momentum(Four_mom_type k, double ma, double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar) {
  switch(k){  
  case Four_mom_type::p1: {
    double k0{ Ep1(s, ma) };
    double k1{ 0 };
    double k2{ 0 };
    double k3{ std::sqrt(Ep1(s, ma)*Ep1(s, ma) - ma*ma) };
    return vect{ k0, k1, k2, k3 };
  }

  case Four_mom_type::p2: {
    double k0{ Ep2(s, ma) };
    double k1{ 0 };
    double k2{ 0 };
    double k3{ Ep2(s, ma) };
    return vect{ k0, k1, k2, k3 };
  }
  case Four_mom_type::q1: {
    double k0{ q10(m1, m2, m3, s, W, thetaStar) };
    double k1{ q11(m1, m2, m3, s, W, theta, thetaStar, phiStar) };
    double k2{ q12(m1, m2, W, thetaStar, phiStar) };
    double k3{ q13(m1, m2, m3, s, W, theta, thetaStar, phiStar) };
    return vect{ k0, k1, k2, k3 };
  }
  case Four_mom_type::q2: {
    double k0{ q20(m1, m2, m3, s, W, thetaStar) };
    double k1{ q21(m1, m2, m3, s, W, theta, thetaStar, phiStar) };
    double k2{ q22(m1, m2, W, thetaStar, phiStar) };
    double k3{ q23(m1, m2, m3, s, W, theta, thetaStar, phiStar) };
    return vect{ k0, k1, k2, k3 };
  }
  case Four_mom_type::q3: {
    double k0{ E3(s, m3, W) };
    double k1{ q3mod(s, m3, W) * std::sin(theta) };
    double k2{ 0 };
    double k3{ q3mod(s, m3, W) * std::cos(theta) };
    return vect{ k0, k1, k2, k3 };
  }
  case Four_mom_type::q: {
    double k0{ Ep2(s, ma) - E3(s, m3, W) };
    double k1{ -q3mod(s, m3, W) * std::sin(theta) };
    double k2{ 0 };
    double k3{ Ep2(s, ma) - q3mod(s, m3, W) * std::cos(theta) };
    return vect{ k0, k1, k2, k3 };
  }
  default:
    return {};
  }
  return {};
}

double Pair(vect k1, vect k2) {
  return -k1*k2;
}

double Power(double x, int n) {
  return std::pow(x, n);
}

double Eps(vect k1, vect k2, vect k3, vect k4) {
  double ret{ 0 };
  for(int i{ 0 }; i < 4; i++) {
    for(int j{ 0 }; j < 4; j++) {
      for(int k{ 0 }; k < 4; k++) {  
        for(int l{ 0 }; l < 4; l++) {
          if(!hasDuplicateIndices(i, j, k, l)) { 
            ret += k1[i]*k2[j]*k3[k]*k4[l]*table[i][j][k][l];
          }
        }
      }
    }
  }
  return ret;
}

}
