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
  
  double numerator{ TwoThreeScatter::kallen(std::sqrt(s), W, m3) * TwoThreeScatter::kallen(W, m1, m2) };
  double denominator{ 128 * std::pow(2*3.14, 5) * p1p2 * W * s };

  double integral = simpson2D(thetaStar_min, thetaStar_max, phiStar_min, phiStar_max, thetaStar_n, phiStar_n, ma, m1, m2, m3, s, W, theta);

  return (numerator / denominator) * integral;
}


// take 2


// { p1, p2, q1, q2, q3, q }
double single_kaon_sum_square_matrix_element(double ma, double m1, double m2, double m3, vect vects[6]) {

  constexpr double GF{ 1.16639e-11 };

  // this will include the neccessary contributions to the matrix element
  // neutrino or anti-neutrinio?
  // have a switch cas of Params here

  double mat{ 0 };  // add all contributing matrix elements to mat

  //mat += TwoThreeScatter::singlekaon::CT(ma, m1, m2, m3, vects);
  //mat += TwoThreeScatter::singlekaon::CrossSigma(ma, m1, m2, m3, vects);
  //mat += TwoThreeScatter::singlekaon::CrossLambda(ma, m1, m2, m3, vects);
  //mat += TwoThreeScatter::singlekaon::PionInFlight(ma, m1, m2, m3, vects);
  //mat += TwoThreeScatter::singlekaon::EtaInFlight(ma, m1, m2, m3, vects);
  //mat += TwoThreeScatter::singlekaon::KaonPole(ma, m1, m2, m3, vects);

	double p1 =  TwoThreeScatter::singlekaon::EtaTest1(ma, m1, m2, m3, vects);
	double p2 = TwoThreeScatter::singlekaon::EtaTest2(ma, m1, m2, m3, vects);
	double p3 = TwoThreeScatter::singlekaon::EtaTest3(ma, m1, m2, m3, vects);
	double p4 = TwoThreeScatter::singlekaon::EtaTest4(ma, m1, m2, m3, vects);
 
	if((p1+p2)<(p3+p4)){
		std::cout<< "failiure\n";
	} 

  return GF*GF * mat;
   
}


// { p1, p2, q1, q2, q3, q }
double single_kaon_diff_xsec_2(double s, double N1_mass, double kaon_mass, double lepton_mass, vect vects[6]) {
  
  constexpr int num_decays{ 1 };

  vect q12{ vects[2] + vects[3] };

	// integrate over W
	double W_min{ kaon_mass + N1_mass };
	double W_max{ std::sqrt(s) - lepton_mass };
	//double W{ W_min + W_max*frandom() };
	double dW{ W_max - W_min };

  particle hadron_blob;
  hadron_blob.p4() = q12;
  hadron_blob.set_mass(std::sqrt(q12*q12));
  //hadron_blob.set_mass(W);

  vec pair_dir{ q12.v() };
  vec vcms{ (vects[0] + vects[1]).v() };
	
  double diff_xsec{ 0 };
  
  for(int i{ 0 }; i < num_decays; i++) {

    particle N1_star;
    particle kaon_star;
    N1_star.set_mass(N1_mass);
    kaon_star.set_mass(kaon_mass);
    if(!hadron_blob.decay(N1_star, kaon_star)) { return 0; } // decay forbidden 

    	vect q1_star{ N1_star };
    	vect q2_star{ kaon_star };

//    q1_star.boost(-vcms);
  //  q2_star.boost(-vcms);

    	diff_xsec += single_kaon_sum_square_matrix_element(N1_mass, N1_mass, kaon_mass, lepton_mass, vects);
  }

	vect q{ vects[5] };
	double Q2{ -q*q };
	double dipoleff{ ( 1/(1 + (Q2/1e6)) ) * ( 1/(1 + (Q2/1e6)) ) }; 
	//double dipoleff{ 1 }; 

  double W{ std::sqrt(q12*q12) };

  double numerator{ TwoThreeScatter::kallen(std::sqrt(s), W, lepton_mass) * TwoThreeScatter::kallen(W, N1_mass, kaon_mass) };
  double denominator{ 128 * std::pow(2*M_PI, 5) * (vects[0]*vects[1]) * W * s };

  return dipoleff * 4*M_PI * (1.0/num_decays) * dW * (numerator / denominator) * diff_xsec;
}







