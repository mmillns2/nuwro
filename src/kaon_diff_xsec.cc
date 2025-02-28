#include "23Scatter.h"
#include "kaon_diff_xsec.h"
#include "single_kaon_matrix_elements/matrix_elements.h"



double single_kaon_integrand(double ma, double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar) {

  constexpr double GF{ 1.16639e-5 };

  // this will include the neccessary contributions to the matrix element
  // neutrino or anti-neutrinio?
  // have a switch cas of Params here

  double mat{ 0 };  // add all contributing matrix elements to mat

  mat += TwoThreeScatter::singlekaon::CT(ma, m1, m2, m3, s, W, theta, thetaStar, phiStar);

  return 0.25 * GF*GF * mat * std::sin(thetaStar);
}


double simpson2D(double thetaStar_min, double thetaStar_max, double phiStar_min, double phiStar_max, int thetaStar_n, 
                 int phiStar_n, double ma, double m1, double m2, double m3, double s, double W, double theta) {

  double h_thetaStar{ (thetaStar_max - thetaStar_min) / thetaStar_n };
  double h_phiStar{ (phiStar_max - phiStar_min) / phiStar_n };

  double integral{ 0 };

  for(int i = 0; i <= thetaStar_n; i++) {

    double thetaStar{ thetaStar_min + i*h_thetaStar };
    double w_thetaStar{ static_cast<double>((i == 0 || i == thetaStar_n) ? 1 : (i % 3 == 0 ? 2 : 3)) };

    for(int j = 0; j <= phiStar_n; j++) {

      double phiStar{ phiStar_min + j*h_phiStar };
      double w_phiStar{ static_cast<double>((j == 0 || j == phiStar_n) ? 1 : (j % 3 == 0 ? 2 : 3)) };
    
      integral += w_thetaStar * w_phiStar * single_kaon_integrand(ma, m1, m2, m3, s, W, theta, thetaStar, phiStar);
    }
  }
  integral *= std::pow(3.0/8.0, 2) * h_thetaStar * h_phiStar;
  return integral;
}


double single_kaon_diff_xsec(double ma, double m1, double m2, double m3, double s, double W, double theta) {
  // step numbers
  int thetaStar_n{ 36 };
  int phiStar_n{ 36 };

  // integration limits
  double thetaStar_min{ 0 };
  double thetaStar_max{ M_PI };
  double phiStar_min{ 0 };
  double phiStar_max{ 2*M_PI };

  // integral
  double p1p2{ TwoThreeScatter::Ep1(s, ma)*TwoThreeScatter::Ep2(s, ma) - 
               TwoThreeScatter::Ep2(s, ma)*std::sqrt(TwoThreeScatter::Ep1(s, ma)*TwoThreeScatter::Ep1(s, ma) - ma*ma) };
  
  double numerator{ TwoThreeScatter::kallen(std::sqrt(s), W, m3) * TwoThreeScatter::kallen(W, m1, m2) };
  double denominator{ 128 * std::pow(2*M_PI, 5) * p1p2 * W * s };

  double integral = simpson2D(thetaStar_min, thetaStar_max, phiStar_min, phiStar_max, thetaStar_n, phiStar_n, ma, m1, m2, m3, s, W, theta);

  return (numerator / denominator) * integral;
}
