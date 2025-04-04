#include "23phaseSpace.h"


#include <cmath>
#include <cassert>


namespace phase23 {

// ------ CONSTRUCTOR(s) - private ------

TwoThreePhaseSpace::TwoThreePhaseSpace(double ma_, double m1_, double m2_, double m3_, double EbeamMin_, double EbeamMax_) 
  : EbeamMin{ std::max(0.0, EbeamMin_ - 1000) }, EbeamMax{ EbeamMax_ + 10000 }, ma{ ma_ }, m1{ m1_ }, m2{ m2_ }, m3{ m3_ } {
		
	//s{ sCalc(ma_, Ebeam) }
	//double physical{ s + m3*m3 - 2*m3*(Ep1() + Ep2()) };
	//assert(physical > 0);
	//WMax = std::sqrt(physical);	
	//assert(WMax > WMin);
}


// ------ PUBLIC MEMBER FUNCTIONS ------

double TwoThreePhaseSpace::xsec(double Ebeam, vect p1_v, vect& q1_v, vect& q2_v, vect& q3_v) {

  double theta{ std::asin(frandom()) };

	std::shared_ptr<TH3D> curr{ getHistogram(Ebeam, theta) }; 
	if(!curr) return 0.0;	// nullptr => forbidden kinematics

  // initialize variables
  double W{ };
  double thetaStar{ };
  double phiStar{ };

  // set variables based on histogram
  if(curr) curr->GetRandom3(W, theta, thetaStar);

  // set 4-momentum in CMS frame 
  q1_v = Momentum(s, q1, W, theta, thetaStar, phiStar);
  q2_v = Momentum(s, q2, W, theta, thetaStar, phiStar);
  q3_v = Momentum(s, q3, W, theta, thetaStar, phiStar);

  // boost to lab frame
  vec vlab{ p1_v.v() };
  q1_v.boost(-vlab);
  q2_v.boost(-vlab);
  q3_v.boost(-vlab);

  // maybe do some conservation checks here?

  // extract cross section
  double xsec{ curr ? curr->Integral("WIDTH") : 0 };
  assert(xsec >= 0);

  return xsec;
}

double TwoThreePhaseSpace::kallen(double x, double y, double z) const {
	double ret{ std::pow(x, 4) + std::pow(y, 4) + std::pow(z, 4) - 2*std::pow(x, 2)*std::pow(y, 2) 
				 - 2*std::pow(x, 2)*std::pow(z, 2) - 2*std::pow(y, 2)*std::pow(z, 2) }; 	
	if(fabs(ret) < 1e-10)
		    ret = 0;
	assert(ret >= 0);
	return std::sqrt(ret);
}


// ------ PRIVATE MEMBER FUNCTIONS ------

void TwoThreePhaseSpace::xsecPrecompute() {

  m_phaseSpace.resize(nEbeam);
	for(auto& thetaVec : m_phaseSpace) thetaVec.resize(ntheta);

  double EbeamStep{ (EbeamMax - EbeamMin) / nEbeam };
  double thetaStep{ (thetaMax - thetaMin) / ntheta };

	for(int iEbeam{ 0 }; iEbeam < nEbeam; ++iEbeam) {

		double Ebeam{ EbeamMin + iEbeam*EbeamStep };
		assert(Ebeam >= EbeamMin && Ebeam <= EbeamMax);

		// check physical
		double s{ sCalc(ma, Ebeam) };
		double sMin{ (m1 + m2 + m3)*(m1 + m2 + m3) };
		bool forbidden{ s < sMin };
		bool unphysical{ false };
		double WMax{ };
		double WMin{ };
		if(!forbidden) {
			WMin = m1 + m2;
			double temp{ s + m3*m3 - 2*m3*(Ep1(s) + Ep2(s)) };
			bool unphysical = (temp < 0);
			WMax = unphysical ? -1 : std::sqrt(temp);	
			unphysical = (WMax <= WMin);	// outside phase space => 0 cross section
		}
		
		for(int itheta{ 0 }; itheta < ntheta; ++itheta) {

			if(unphysical || forbidden) {	// set all TH3Ds to nullptr => 0 cross section
				m_phaseSpace[iEbeam][itheta] = nullptr;
				continue;
			}

			assert(WMax > WMin);

			double theta{ thetaMin + itheta*thetaStep };
			assert(theta >= thetaMin && theta <= thetaMax);

			// create TH3
			std::string name{ "hist_theta_" + std::to_string(itheta) + "_" + [this]() -> std::string {
					HISTID++;
					return std::to_string(HISTID);
				}()
			};

			m_phaseSpace[iEbeam][itheta] = std::make_shared<TH3D>(TH3D(name.c_str(), name.c_str(), nW, WMin, WMax, 
																			nthetaStar, thetaStarMin, thetaStarMax,
																			nphiStar, phiStarMin, phiStarMax));

			// compute cross secions
			for(int iW{ 1 }; iW < m_phaseSpace[iEbeam][itheta]->GetNbinsX() + 1; ++iW) {

				double W{ m_phaseSpace[iEbeam][itheta]->GetXaxis()->GetBinCenter(iW) }; 
				assert(W >= WMin && W <= WMax);

				for(int ithetaStar{ 1 }; ithetaStar < m_phaseSpace[iEbeam][itheta]->GetNbinsY() + 1; ++ithetaStar) {

				double thetaStar{ m_phaseSpace[iEbeam][itheta]->GetYaxis()->GetBinCenter(ithetaStar) };
				assert(thetaStar >= thetaStarMin && thetaStar <= thetaStarMax);

					for(int iphiStar{ 1 }; iphiStar < m_phaseSpace[iEbeam][itheta]->GetNbinsZ() + 1; ++iphiStar) {
						
						double phiStar{ m_phaseSpace[iEbeam][itheta]->GetZaxis()->GetBinCenter(iphiStar) };
						assert(phiStar >= phiStarMin && phiStar <= phiStarMax);

						vect p1_v{ Momentum(s, p1, W, theta, thetaStar, phiStar) };
						vect p2_v{ Momentum(s, p2, W, theta, thetaStar, phiStar) };
						vect q1_v{ Momentum(s, q1, W, theta, thetaStar, phiStar) };
						vect q2_v{ Momentum(s, q2, W, theta, thetaStar, phiStar) };
						vect q3_v{ Momentum(s, q3, W, theta, thetaStar, phiStar) };

						double s_new1{ (p1_v + p2_v)*(p1_v + p2_v) };
						double s_new2{ (q1_v + q2_v + q3_v)*(q1_v + q2_v + q3_v) };
						assert(fabs(s - s_new1) < 1e-5); 
						assert(fabs(s_new2 - s_new1) < 1e-5); 

						vect q_new1{ p2_v - q3_v };
						vect q_new2{ q1_v + q2_v - p1_v };
						assert(fabs(q_new1[0] - q_new2[0]) < 1e-5); 
						assert(fabs(q_new1[1] - q_new2[1]) < 1e-5); 
						assert(fabs(q_new1[2] - q_new2[2]) < 1e-5); 
						assert(fabs(q_new1[3] - q_new2[3]) < 1e-5); 

						double xsec{ xsecCalc(s, W, theta, thetaStar, phiStar) };
						assert(xsec >= 0);

						// bin cross section
						m_phaseSpace[iEbeam][itheta]->SetBinContent(iW, ithetaStar, iphiStar, xsec);
					}
				}
			}
		}
	}
}

int TwoThreePhaseSpace::getThetaIndex(double theta) const {

  double thetaStep{ (thetaMax - thetaMin) / ntheta };
  int i{ static_cast<int>((theta - thetaMin) / thetaStep) };

  assert(i >= 0);
  assert(i < ntheta);

  //if(i < 0) return 0;
  //if(i >= ntheta) return ntheta - 1;

  return i;
}

int TwoThreePhaseSpace::getEbeamIndex(double Ebeam) const {

  double EbeamStep{ (EbeamMax - EbeamMin) / nEbeam };
  int i{ static_cast<int>((Ebeam - EbeamMin) / EbeamStep) };

	//std::cout << Ebeam << '\n';
  assert(i >= 0);
	
	if(i >= nEbeam) i = nEbeam - 1; 

  assert(i < nEbeam);

  return i;
}

double TwoThreePhaseSpace::sCalc(double ma_, double Ebeam) const {
  double ret{ (ma_ + Ebeam)*(ma_ + Ebeam) - Ebeam*Ebeam }; 
	//std::cout << "Ebeam=" << Ebeam << '\n';
	//std::cout << "ma=" << ma_ << '\n';
	//std::cout << "s=" << ret << '\n';
  assert(ret > 0); // add better physical check
  return ret;
}

double TwoThreePhaseSpace::dipoleffCalc(double s, double W, double theta, double thetaStar, double phiStar) const {
  vect q_{ Momentum(s, q, W, theta, thetaStar, phiStar) }; 
  double Q2{ -q_ * q_ };
  double M2{ 1000 * 1000 };
  assert(fabs(Q2 + M2) > 1e-5); // Q2 != M2
  return 1 / ( ( 1 + (Q2/M2) ) * ( 1 + (Q2/M2) ) );
}

double TwoThreePhaseSpace::xsecCalc(double s, double W, double theta, double thetaStar, double phiStar) const {
  double mat{ matrixElement(s, W, theta, thetaStar, phiStar) };  
  double k1{ kallen(std::sqrt(s), W, m3) };
  double k2{ kallen(W, m1, m2) };
  double p1p2{ Pair(Momentum(s, p1, W, theta, thetaStar, phiStar), Momentum(s, p2, W, theta, thetaStar, phiStar)) };
  double dipoleff{ dipoleffCalc(s, W, theta, thetaStar, phiStar) };
	dipoleff = 1;
  assert(fabs(p1p2) > 1e-5);  // p1p2 != 0
  assert(fabs(W) > 1e-5);  // W != 0
  assert(fabs(s) > 1e-5);  // s != 0
  return (mat * k1 * k2 * dipoleff) / (128 * std::pow(2*M_PI, 5) * p1p2 * W * s);
}

double TwoThreePhaseSpace::Ep1(double s) const {
  assert(s > 0);
	return (s + ma*ma) / (2*std::sqrt(s));
}

double TwoThreePhaseSpace::Ep2(double s) const {
	double numerator{ s - ma*ma };
	double denominator{ 2*(Ep1(s) + std::sqrt(Ep1(s)*Ep1(s) - ma*ma)) };
	double ret{ numerator / denominator };
  assert(ret > 0);
  return ret;
}

double TwoThreePhaseSpace::E3(double s, double W) const {
  assert(s > 0);
  return (s + m3*m3 - W*W) / (2*std::sqrt(s));
}

double TwoThreePhaseSpace::q3mod(double s, double W) const {
  assert(s > 0);
  return kallen(std::sqrt(s), W, m3) / (2*std::sqrt(s));
}

double TwoThreePhaseSpace::q12CM(double W) const {
  assert(W > 0);
	return kallen(W, m1, m2) / (2*W);
}

double TwoThreePhaseSpace::E12(double s, double W) const {
  assert(s > 0);
  double ret{ (s + W*W - m3*m3) / (2*std::sqrt(s)) };
  assert(ret > 0);
  return ret;
}

double TwoThreePhaseSpace::gamma(double s, double W) const {
  double ret{ E12(s, W) / W };
  if(fabs(1 - ret) < 0.0000001) ret = 1;
  assert(ret >= 1);
  return ret;
}

double TwoThreePhaseSpace::beta(double s, double W) const {
  double ret{ 1 - 1 / (gamma(s, W) * gamma(s, W)) };
  if(fabs(ret) < 1e-5) ret = 0;
  assert(ret >= 0);
  return std::sqrt(ret);
}

double TwoThreePhaseSpace::q10Star(double W) const {
  assert(W > 0);
  return (W*W + m1*m1 - m2*m2) / (2*W);
}

double TwoThreePhaseSpace::q20Star(double W) const {
  assert(W > 0);
  return (W*W + m2*m2 - m1*m1) / (2*W);
}

double TwoThreePhaseSpace::q11Star(double W, double thetaStar, double phiStar) const {
  return q12CM(W) * std::sin(thetaStar) * std::cos(phiStar);
}

double TwoThreePhaseSpace::q21Star(double W, double thetaStar, double phiStar) const {
  return -q12CM(W) * std::sin(thetaStar) * std::cos(phiStar);
}

double TwoThreePhaseSpace::q12Star(double W, double thetaStar, double phiStar) const {
  return q12CM(W) * std::sin(thetaStar) * std::sin(phiStar);
}

double TwoThreePhaseSpace::q22Star(double W, double thetaStar, double phiStar) const {
  return -q12CM(W) * std::sin(thetaStar) * std::sin(phiStar);
}

double TwoThreePhaseSpace::q13Star(double W, double thetaStar) const {
  return q12CM(W) * std::cos(thetaStar);
}

double TwoThreePhaseSpace::q23Star(double W, double thetaStar) const {
  return -q12CM(W) * std::cos(thetaStar);
}

double TwoThreePhaseSpace::q10(double s, double W, double thetaStar) const {
  return gamma(s, W) * (q10Star(W) + beta(s, W) * q13Star(W, thetaStar));
}

double TwoThreePhaseSpace::q20(double s, double W, double thetaStar) const {
  return gamma(s, W) * (q20Star(W) + beta(s, W) * q23Star(W, thetaStar));
}

double TwoThreePhaseSpace::q11(double s, double W, double theta, double thetaStar, double phiStar) const {
  return q11Star(W, thetaStar, phiStar) * std::cos(theta + M_PI) + 
         gamma(s, W) * (beta(s, W) * q10Star(W) + q13Star(W, thetaStar)) * std::sin(theta + M_PI);
}

double TwoThreePhaseSpace::q21(double s, double W, double theta, double thetaStar, double phiStar) const {
  return q21Star(W, thetaStar, phiStar) * std::cos(theta + M_PI) + 
         gamma(s, W) * (beta(s, W) * q20Star(W) + q23Star(W, thetaStar)) * std::sin(theta + M_PI);
}

double TwoThreePhaseSpace::q12(double W, double thetaStar, double phiStar) const {
  return q12Star(W, thetaStar, phiStar);
}

double TwoThreePhaseSpace::q22(double W, double thetaStar, double phiStar) const {
  return q22Star(W, thetaStar, phiStar);
}

double TwoThreePhaseSpace::q13(double s, double W, double theta, double thetaStar, double phiStar) const {
  return -q11Star(W, thetaStar, phiStar) * std::sin(theta + M_PI) + 
         gamma(s, W) * (beta(s, W) * q10Star(W) + q13Star(W, thetaStar)) * std::cos(theta + M_PI);
}

double TwoThreePhaseSpace::q23(double s, double W, double theta, double thetaStar, double phiStar) const {
  return -q21Star(W, thetaStar, phiStar) * std::sin(theta + M_PI) + 
         gamma(s, W) * (beta(s, W) * q20Star(W) + q23Star(W, thetaStar)) * std::cos(theta + M_PI);
}

vect TwoThreePhaseSpace::Momentum(double s, Four_mom_type k,double W, double theta, double thetaStar, double phiStar) const {
  switch(k){  
  case Four_mom_type::p1: {
    double k0{ Ep1(s) };
    double k1{ 0 };
    double k2{ 0 };
    double k3{ -std::sqrt(Ep1(s)*Ep1(s) - ma*ma) };
    return vect{ k0, k1, k2, k3 };
  }

  case Four_mom_type::p2: {
    double k0{ Ep2(s) };
    double k1{ 0 };
    double k2{ 0 };
    double k3{ Ep2(s) };
    return vect{ k0, k1, k2, k3 };
  }
  case Four_mom_type::q1: {
    double k0{ q10(s, W, thetaStar) };
    double k1{ q11(s, W, theta, thetaStar, phiStar) };
    double k2{ q12(W, thetaStar, phiStar) };
    double k3{ q13(s, W, theta, thetaStar, phiStar) };
    return vect{ k0, k1, k2, k3 };
  }
  case Four_mom_type::q2: {
    double k0{ q20(s, W, thetaStar) };
    double k1{ q21(s, W, theta, thetaStar, phiStar) };
    double k2{ q22(W, thetaStar, phiStar) };
    double k3{ q23(s, W, theta, thetaStar, phiStar) };
    return vect{ k0, k1, k2, k3 };
  }
  case Four_mom_type::q3: {
    double k0{ E3(s, W) };
    double k1{ q3mod(s, W) * std::sin(theta) };
    double k2{ 0 };
    double k3{ q3mod(s, W) * std::cos(theta) };
    return vect{ k0, k1, k2, k3 };
  }
  case Four_mom_type::q: {  // p2 - q3
    double k0{ Ep2(s) - E3(s, W) };
    double k1{ -q3mod(s, W) * std::sin(theta) };
    double k2{ 0 };
    double k3{ Ep2(s) - q3mod(s, W) * std::cos(theta) };
    return vect{ k0, k1, k2, k3 };
  }
  default:
    return {};
  }
  return {};
}

double TwoThreePhaseSpace::Pair(vect k1, vect k2) const {
  return k1*k2;
}

double TwoThreePhaseSpace::Power(double x, int n) const {
  return std::pow(x, n);
}

double TwoThreePhaseSpace::Power(vect x, int n) const {
  assert(n == 2);
  return x*x;
}

double TwoThreePhaseSpace::Eps(vect k1, vect k2, vect k3, vect k4) const {
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
 
constexpr bool TwoThreePhaseSpace::hasDuplicateIndices(int i, int j, int k, int l) const {
  return (i == j) || (i == k) || (i == l) || (j == k) || (j == l) || (k == l);
}

constexpr int TwoThreePhaseSpace::getLeviCivitaValue(int i, int j, int k, int l) const {
  return hasDuplicateIndices(i, j, k, l) ? 0 : table[i][j][k][l];
}



//****************** Single Kaon *********************

namespace singlek {


//------ CONSTRUCTOR(s) - PRIVATE ------



SingleKaonPhaseSpace::SingleKaonPhaseSpace(Channel channel, double m3_, double EbeamMin_, double EbeamMax_)
  : TwoThreePhaseSpace{ 
    [channel]() -> double {   
      switch (channel) {    
      case Channel::pp: return mProton; 
      case Channel::np: return mNeutron; 
      case Channel::nn: return mNeutron; 
			default: return mProton;					// change later to -1?
      }
    }(), 
    [channel]() -> double {   
      switch (channel) {    
      case Channel::pp: return mProton; 
      case Channel::np: return mProton; 
      case Channel::nn: return mNeutron; 
			default: return mProton;
      }
    }(),
    [channel]() -> double {   
      switch (channel) {    
      case Channel::pp: return mKaonCharged; 
      case Channel::np: return mKaon0; 
      case Channel::nn: return mKaonCharged; 
			default: return mProton;
      }
    }(), m3_, EbeamMin_, EbeamMax_ }, m_channel{ channel } { 

  // set form factors
  switch(channel) {

  case Channel::pp: 
    ACT = 2;
    BCT = -F;
    ACRSigma = -(D - F)/2;
    ACRLambda = D + 3*F;
    AKP = 2;
    APi = -1;
    AEta = 1;
    break;

  case Channel::np: 
    ACT = 1;
    BCT = -(D + F);
    ACRSigma = (D - F)/2;
    ACRLambda = D + 3*F;
    AKP = 1;
    APi = -2;
    AEta = 0;
    break;

  case Channel::nn: 
    ACT = 1;
    BCT = D - F;
    ACRSigma = -(D - F);
    ACRLambda = 0;
    AKP = 1;
    APi = 1;
    AEta = 1;
    break;
	
	default:
		std::cout << "Error in SingleKaonPhaseSpace constructor\n";
		break;
  }

	//compute phase space
  xsecPrecompute(); 
}


//------ PRIVATE MEMBER FUNCTIONS ------


double SingleKaonPhaseSpace::matrixElement(double s, double W, double theta, double thetaStar, double phiStar) const {

  double ret{ 0 };

  ret += CT(s, W, theta, thetaStar, phiStar);
	//double temp{ CrossLambda(s, W, theta, thetaStar, phiStar) };
  //ret += (temp < 0 ? 0 : temp);
  //ret += CrossSigma(s, W, theta, thetaStar, phiStar);
  //ret += PionInFlight(s, W, theta, thetaStar, phiStar);
  //ret += EtaInFlight(s, W, theta, thetaStar, phiStar);
  //ret += KaonPole(s, W, theta, thetaStar, phiStar);
	
	//std::cout << "ret=" << ret << '\n';

	//if(ret < 0) return 0;

  assert(ret > 0);

  return GF * GF * ret;
}


//------ MATRIX ELEMENTS ------


double SingleKaonPhaseSpace::CT(double s, double W, double theta, double thetaStar, double phiStar) const {
  // Not included GF squared factor here; have included the missing 1/4 factor
  return (16*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
  (32*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
  (16*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))/Power(fPi,2) - 
  (16*Power(ACT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
  (16*Power(ACT,2)*Power(BCT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
  (16*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))/Power(fPi,2) - 
  (32*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
  (16*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))/Power(fPi,2);
}

double SingleKaonPhaseSpace::CrossLambda(double s, double W, double theta, double thetaStar, double phiStar) const {
  // Not included GF squared factor here; have included the missing 1/4 factor
  return (Power(ACRLambda,2)*Power(Vus,2)*((32*Power(D + 3*F,2)*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2))/9. + 
  32*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) + 
  32*Power(mLambda,2)*(-((Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  (32*Power(D + 3*F,2)*Power(mLambda,2)*(-((Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
         m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
       Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/9. + 
  32*(4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
     m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  (32*Power(D + 3*F,2)*(4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/9. + 
  (16*Power(D + 3*F,2)*((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2)) + 
  (16*Power(D + 3*F,2)*Power(mLambda,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (4*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    ((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2)) - 
  (16*Power(D + 3*F,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
     (4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     (-4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        (4*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) - 
       m1*ma*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2)) + 
  (2*Power(kappaP,2)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
  (2*Power(kappaP,2)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
  (2*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
       Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
  (2*Power(kappaP,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
       Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) + (2*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) + (2*Power(kappaP,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
    4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2)))/(32*Power(fPi,2)*Power(-Power(mLambda,2) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar)),2)) + 
  (Power(ACRLambda,2)*Power(Vus,2)*((-64*(D + 3*F)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2))/3. - 
  64*mLambda*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  (ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
     m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
  (64*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    (ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))))/9. + 
  (64*Power(D + 3*F,2)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/9. + 
  64*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
  (64*(D + 3*F)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/3. - 
  (64*(D + 3*F)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    (-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))) + ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/3. + 
  (64*(D + 3*F)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/3. + 
  (64*(D + 3*F)*mLambda*(ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/3. + 
  (64*(D + 3*F)*mLambda*(ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/3. + 
  (32*(D + 3*F)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*m1*ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/3. + 
  (64*(D + 3*F)*Power(mLambda,2)*((Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/3. - 
  (64*(D + 3*F)*(4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/3. - 
  (64*Power(D + 3*F,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/9. - 
  64*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  (16*kappaP*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*mLambda*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))) - 
    2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     ((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*(-4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/ma + 
  (16*kappaP*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/ma + 
  (16*kappaP*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma + 
  (16*kappaP*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma - 
  (32*Power(D + 3*F,2)*mLambda*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2)) - 
  (32*Power(D + 3*F,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
       4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2)) + 
  (32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) - 
    (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2)) + 
  (32*Power(D + 3*F,2)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    ((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/
  (9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) + (16*kappaP*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma + 
  (16*kappaP*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma - 
  (4*(D + 3*F)*kappaP*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/(3.*ma) - 
  (4*Power(kappaP,2)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(ma,2) - 
  (4*Power(kappaP,2)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
  (16*(D + 3*F)*kappaP*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/(3.*ma) - 
  (4*Power(kappaP,2)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((-2*m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*(m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       (-2*m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
  (32*Power(D + 3*F,2)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/(9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) - 
  (32*Power(D + 3*F,2)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/(9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) + 
  (4*Power(kappaP,2)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
  (32*Power(D + 3*F,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
          2*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/
  (9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) - (32*Power(D + 3*F,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
          2*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/
  (9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) + (32*Power(D + 3*F,2)*Power(mLambda,2)*
  (Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        ((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/
  (9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) + (4*Power(kappaP,2)*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/
  Power(ma,2) - (32*(D + 3*F)*kappaP*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
  (16*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*(4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (32*(D + 3*F)*kappaP*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/
  (3.*ma) + (32*(D + 3*F)*kappaP*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/
  (3.*ma) + (32*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
  (32*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
  (32*(D + 3*F)*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/(3.*ma)
  + (4*Power(kappaP,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (-4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
          4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((-(m1*ma) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (-(m1*ma) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
  (4*Power(kappaP,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (-(m1*ma) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (-(m1*ma) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
  (32*Power(D + 3*F,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       (-((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/(9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) + 
  (4*Power(kappaP,2)*mLambda*(-2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
  (4*Power(kappaP,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(4*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
       Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
  (16*kappaP*mLambda*(2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*kappaP*mLambda*(2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (4*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((2*m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*(m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (2*m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
  (4*Power(kappaP,2)*mLambda*(-2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) - (4*Power(kappaP,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) + (32*(D + 3*F)*kappaP*mLambda*(-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     ((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/(3.*ma) + 
  (32*(D + 3*F)*kappaP*mLambda*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     ((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/(3.*ma) - 
  (16*kappaP*mLambda*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     ((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/ma - 
  (32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
        (-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*
        (-(m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/(9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) + 
  (32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
        (-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*
        (-(m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/(9.*(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)))) + 
  (4*Power(kappaP,2)*mLambda*(-2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
              (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
  (4*Power(kappaP,2)*mLambda*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
              (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
  (4*Power(kappaP,2)*(-(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
       (4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
       (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
         Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*(4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - m1*ma*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
              (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
                2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
                2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
                Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))))/Power(ma,2)))/
  (32*Power(fPi,2)*Power(-Power(mLambda,2) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar)),2));
}

double SingleKaonPhaseSpace::CrossSigma(double s, double W, double theta, double thetaStar, double phiStar) const {
  // Not included GF squared factor here; have included the missing 1/4 factor
  return (Power(ACRSigma,2)*Power(Vus,2)*(32*Power(D - F,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
  Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) + 32*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) + 
  32*Power(mSigma,2)*(-((Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  32*Power(D - F,2)*Power(mSigma,2)*(-((Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  32*(4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
     m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  32*Power(D - F,2)*(4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  (16*Power(D - F,2)*((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2) + 
  (16*Power(D - F,2)*Power(mSigma,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (4*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    ((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2) - 
  (16*Power(D - F,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
     (4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     (-4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        (4*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) - 
       m1*ma*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2) + 
  (2*Power(2*kappaN + kappaP,2)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
  (2*Power(2*kappaN + kappaP,2)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) - 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
  (2*Power(2*kappaN + kappaP,2)*Power(mSigma,2)*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
       Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
  (2*Power(2*kappaN + kappaP,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
       Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) + (2*Power(2*kappaN + kappaP,2)*Power(mSigma,2)*
  (Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) - (2*Power(2*kappaN + kappaP,2)*(-(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
       (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
         2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
         Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
            Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
         2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
         m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
         Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
         Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
            Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
            Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))) + 
    4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2)))/(8.*Power(fPi,2)*Power(-Power(mSigma,2) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar)),2)) + 
  (Power(ACRSigma,2)*Power(Vus,2)*(64*(D - F)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - 
  64*mSigma*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  (ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
     m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
  64*Power(D - F,2)*mSigma*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
  (ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
     m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
  64*Power(D - F,2)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
  64*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  64*(D - F)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
  64*(D - F)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
  (-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))) + ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  64*(D - F)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  64*(D - F)*mSigma*(ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
  2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  32*(D - F)*mSigma*(2*ma*(-(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
  4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  32*(D - F)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
  2*m1*ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  64*(D - F)*(4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
  64*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  64*Power(D - F,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  32*(D - F)*Power(mSigma,2)*(2*(-(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-4*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
     4*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
  (16*(2*kappaN + kappaP)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*mSigma*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))) - 
    2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     ((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*(-4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*Power(mSigma,2)*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/ma + 
  (16*(2*kappaN + kappaP)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/ma - 
  (16*(2*kappaN + kappaP)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma + 
  (16*(2*kappaN + kappaP)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma + 
  (32*(D - F)*(2*kappaN + kappaP)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma - 
  (32*Power(D - F,2)*mSigma*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2) - 
  (32*Power(D - F,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
       4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2) + 
  (32*Power(D - F,2)*mSigma*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) - 
    (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2) + 
  (32*Power(D - F,2)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    ((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/
  (Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) + (16*(2*kappaN + kappaP)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma + 
  (16*(2*kappaN + kappaP)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma - 
  (16*(2*kappaN + kappaP)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma - 
  (4*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    4*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/ma + 
  (4*Power(2*kappaN + kappaP,2)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(ma,2) - 
  (4*Power(2*kappaN + kappaP,2)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/Power(ma,2) + 
  (16*(D - F)*(2*kappaN + kappaP)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/ma - 
  (4*Power(2*kappaN + kappaP,2)*Power(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
  (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((-2*m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*(m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       (-2*m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
  (32*Power(D - F,2)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) - 
  (32*Power(D - F,2)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) + 
  (4*Power(2*kappaN + kappaP,2)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
  (32*Power(D - F,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
          2*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/
  (Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) - (32*Power(D - F,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       ((m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
          2*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/
  (Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) + (32*Power(D - F,2)*Power(mSigma,2)*
  (Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        ((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))/
  (Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) + (4*Power(2*kappaN + kappaP,2)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/
  Power(ma,2) + (32*(D - F)*(2*kappaN + kappaP)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*Power(mSigma,2)*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*(4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma - 
  (16*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma - 
  (32*(D - F)*(2*kappaN + kappaP)*mSigma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (32*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma - 
  (32*(D - F)*(2*kappaN + kappaP)*Power(mSigma,2)*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     ((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (4*Power(2*kappaN + kappaP,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (-4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
          4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((-(m1*ma) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (m1*ma - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (-(m1*ma) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
  (4*Power(2*kappaN + kappaP,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (-(m1*ma) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (-(m1*ma) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
  (32*Power(D - F,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (-(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       (-((m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) - 
  (16*(2*kappaN + kappaP)*mSigma*(2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (16*(2*kappaN + kappaP)*mSigma*(2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/ma + 
  (4*Power(2*kappaN + kappaP,2)*Power(mSigma,2)*(Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*((2*m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*(m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          (2*m1*ma + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) - 
    2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
  (4*Power(2*kappaN + kappaP,2)*mSigma*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) - (4*Power(2*kappaN + kappaP,2)*mSigma*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-(m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
          (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) - (4*Power(2*kappaN + kappaP,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) + (4*Power(2*kappaN + kappaP,2)*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(-4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*m1*ma*Power(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)),2) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       4*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/
  Power(ma,2) - (32*(D - F)*(2*kappaN + kappaP)*mSigma*(-2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     ((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/ma - 
  (32*(D - F)*(2*kappaN + kappaP)*mSigma*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     ((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/ma - 
  (16*(2*kappaN + kappaP)*mSigma*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     ((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       m1*ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/ma - 
  (32*Power(D - F,2)*mSigma*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
        (-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*
        (-(m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) + 
  (32*Power(D - F,2)*mSigma*(2*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*
        (-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(ma*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*
        (-(m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          m1*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
       ma*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))) + 
  (4*Power(2*kappaN + kappaP,2)*mSigma*(-2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
              (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
  (4*Power(2*kappaN + kappaP,2)*mSigma*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    4*m1*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*m1*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       ma*(2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*
              (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) - 
  (4*Power(2*kappaN + kappaP,2)*(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
     (4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) - 
    2*(4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2)*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*
        (4*Power(Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)),2) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
        (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar))*((Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - m1*ma*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*
              (-Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
                2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
                2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
             Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-2*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + 
                Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))))/Power(ma,2)))/
  (8.*Power(fPi,2)*Power(-Power(mSigma,2) + Pair(Momentum(s, p1,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, p1,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar)),2));
}

double SingleKaonPhaseSpace::PionInFlight(double s, double W, double theta, double thetaStar, double phiStar) const {
  // Not included GF squared factor here; have included the missing 1/4 factor
  return (-2*Power(APi,2)*Power(D + F,2)*Power(ma,2)*Power(Vus,2)*(m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-4*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    4*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*(2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/(Power(fPi,2)*Power(Power(mPi,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar)),2));
}

double SingleKaonPhaseSpace::EtaInFlight(double s, double W, double theta, double thetaStar, double phiStar) const {
  // Not included GF squared factor here; have included the missing 1/4 factor
  return (-2*Power(AEta,2)*Power(D - 3*F,2)*Power(ma,2)*Power(Vus,2)*(m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*(-4*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    4*Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 
       2*(2*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))))/(Power(fPi,2)*Power(Power(mEta,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)-Momentum(s, q2,W,theta,thetaStar,phiStar)),2));
}

double SingleKaonPhaseSpace::KaonPole(double s, double W, double theta, double thetaStar, double phiStar) const {
  // Not included GF squared factor here; have included the missing 1/4 factor
  return (Power(AKP,2)*Power(Vus,2)*((m1*ma - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + 2*m1*ma*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar)) + Pair(Momentum(s, q1,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar))) + 
  m1*ma*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)) - Pair(Momentum(s, p1,W,theta,thetaStar,phiStar),Momentum(s, q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q2,W,theta,thetaStar,phiStar),Momentum(s, q2,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, p2,W,theta,thetaStar,phiStar))*(Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar)) - 2*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s, p2,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s, q3,W,theta,thetaStar,phiStar),Momentum(s, q3,W,theta,thetaStar,phiStar))))/(2*Power(fPi,2)*Power(Power(m2,2) - Pair(Momentum(s, q,W,theta,thetaStar,phiStar),Momentum(s, q,W,theta,thetaStar,phiStar)),2));
}


}
}


