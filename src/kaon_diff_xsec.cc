#include "23Scatter.h"
#include "kaon_diff_xsec.h"
#include "matrix_elements.h"



double single_kaon_integrand(double ma, double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar) {

  constexpr double GF{ 1.16639e-11 };

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
  int thetaStar_n{ 18 };
  int phiStar_n{ 18 };

  // integration limits
  double thetaStar_min{ 0 };
  double thetaStar_max{ 3.14 };
  double phiStar_min{ 0 };
  double phiStar_max{ 2*3.14 };

  // integral
  double p1p2{ TwoThreeScatter::Ep1(s, ma)*TwoThreeScatter::Ep2(s, ma) - 
               TwoThreeScatter::Ep2(s, ma)*std::sqrt(TwoThreeScatter::Ep1(s, ma)*TwoThreeScatter::Ep1(s, ma) - ma*ma) };
  
	//std::cout << "attemting k1\n";
	//double k1 = TwoThreeScatter::kallen(std::sqrt(s), W, m3);
	//std::cout << "made k1\n";
	//std::cout << "attemting k2\n";
	//double k2 = TwoThreeScatter::kallen(W, m1, m2);
	//std::cout << "made k2\n";
  double numerator{ TwoThreeScatter::kallen(std::sqrt(s), W, m3) * TwoThreeScatter::kallen(W, m1, m2) };
	//double numerator = k1 * k2;
	//std::cout << "numerator made\n";
  double denominator{ 128 * std::pow(2*3.14, 5) * p1p2 * W * s };


  double integral = simpson2D(thetaStar_min, thetaStar_max, phiStar_min, phiStar_max, thetaStar_n, phiStar_n, ma, m1, m2, m3, s, W, theta);

  return (numerator / denominator) * integral;
}


// take 2


// { p1, p2, q1, q2, q3, q }
double single_kaon_sum_square_matrix_element(vect vects[6]) {

  constexpr double GF{ 1.16639e-11 };

  // this will include the neccessary contributions to the matrix element
  // neutrino or anti-neutrinio?
  // have a switch cas of Params here

  double mat{ 0 };  // add all contributing matrix elements to mat

  mat += TwoThreeScatter::singlekaon::CT(vects);

  return 0.25 * GF*GF * mat;
   
}


// { p1, p2, q1, q2, q3, q }
double single_kaon_diff_xsec_2(double s, double N1_mass, double kaon_mass, double lepton_mass, vect vects[6]) {
  
  constexpr int num_decays{ 100 };

  vect q12{ vects[2] + vects[3] };

  particle hadron_blob;
  hadron_blob.p4() = q12;
  hadron_blob.set_mass(std::sqrt(q12*q12));

  vec pair_dir{ q12.v() };
  vec vcms{ (vects[0] + vects[1]).v() };

  double diff_xsec{ 0 };
  
  for(int i{ 0 }; i < num_decays; i++) {

    particle N1_star;
    particle kaon_star;
    N1_star.set_mass(N1_mass);
    kaon_star.set_mass(kaon_mass);
    hadron_blob.decay(N1_star, kaon_star); 

    vect q1_star{ N1_star };
    vect q2_star{ kaon_star };

    q1_star.boost(-vcms);
    q2_star.boost(-vcms);

    diff_xsec += single_kaon_sum_square_matrix_element(vects);
  }

  double W{ std::sqrt(q12*q12) };

  double numerator{ TwoThreeScatter::kallen(std::sqrt(s), W, lepton_mass) * TwoThreeScatter::kallen(W, N1_mass, kaon_mass) };
  double denominator{ 128 * std::pow(2*M_PI, 5) * (vects[0]*vects[1]) * W * s };

  return 4*M_PI * (1.0/num_decays) * (numerator / denominator) * diff_xsec;
}







