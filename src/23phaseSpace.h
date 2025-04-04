#pragma once


#include <memory>
#include <vector>
#include <map>
#include <tuple>
#include <string>


#include "event1.h"
#include "TH3.h"
#include "TH3D.h"


namespace phase23 {

enum Four_mom_type { 
  p1,   // incoming nucleon
  p2,   // incoming neutrino
  q1,   // outgoing nucleon
  q2,   // outgoing kaon
  q3,   // outgoing lepton
  q,    // 4-momentum transfer = p2 - q3 = q1 - p1
  max_types
};


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

static int HISTID{ 1000 };


class TwoThreePhaseSpace {

protected:
  virtual double matrixElement(double s, double W, double theta, double thetaStar, double phiStar) const = 0;
  
public:
  std::shared_ptr<TH3D> getHistogram(double Ebeam, double theta) const { return m_phaseSpace[getEbeamIndex(Ebeam)][getThetaIndex(theta)]; }
  double xsec(double Ebeam, vect p1_v, vect& q1_v, vect& q2_v, vect& q3_v);
  double kallen(double x, double y, double z) const;

public:
  // prevent copying
  TwoThreePhaseSpace(const TwoThreePhaseSpace&) = delete;
  TwoThreePhaseSpace& operator=(const TwoThreePhaseSpace&) = delete;
	TwoThreePhaseSpace(TwoThreePhaseSpace&&) noexcept = default;
  TwoThreePhaseSpace& operator=(TwoThreePhaseSpace&&) noexcept = default;

protected:
  TwoThreePhaseSpace(double ma_, double m1_, double m2_, double m3_, double EbeamMin_, double EbeamMax_);

public:
  virtual ~TwoThreePhaseSpace() { }

protected:
  void xsecPrecompute();
  int getThetaIndex(double theta) const;
  int getEbeamIndex(double Ebeam) const;
  double sCalc(double ma_, double Ebeam) const;
  double dipoleffCalc(double s, double W, double theta, double thetaStar, double phiStar) const;
  double xsecCalc(double s, double W, double theta, double thetaStar, double phiStar) const;

protected:
  std::vector<std::vector<std::shared_ptr<TH3D>>> m_phaseSpace;

  static constexpr int ntheta{ 12 };
  static constexpr int nthetaStar{ 12 };
  static constexpr int nphiStar{ 12 };
  static constexpr int nW{ 12 };
  static constexpr int nEbeam{ 10 };

  static constexpr double thetaMin{ 0.0 };
  static constexpr double thetaMax{ 3.14159 };
  static constexpr double thetaStarMin{ 0.0 };
  static constexpr double thetaStarMax{ 3.14159 };
  static constexpr double phiStarMin{ 0.0 };
  static constexpr double phiStarMax{ 2*3.14159 };

	double EbeamMin;
	double EbeamMax;

  double ma;
  double m1;
  double m2;
  double m3;

private:
  // 4-momentum components
  double Ep1(double s) const;
  double Ep2(double s) const; 
  double E3(double s, double W) const; 
  double q3mod(double s, double W) const;
  double q12CM(double W) const;
  double E12(double s, double W) const;
  double gamma(double s, double W) const; 
  double beta(double s, double W) const; 
  double q10Star(double W) const;
  double q20Star(double W) const;
  double q11Star(double W, double thetaStar, double phiStar) const;
  double q21Star(double W, double thetaStar, double phiStar) const;
  double q12Star(double W, double thetaStar, double phiStar) const;
  double q22Star(double W, double thetaStar, double phiStar) const;
  double q13Star(double W, double thetaStar) const;
  double q23Star(double W, double thetaStar) const;
  double q10(double s, double W, double thetaStar) const;
  double q20(double s, double W, double thetaStar) const;
  double q11(double s, double W, double theta, double thetaStar, double phiStar) const;
  double q21(double s, double W, double theta, double thetaStar, double phiStar) const;
  double q12(double W, double thetaStar, double phiStar) const; 
  double q22(double W, double thetaStar, double phiStar) const;
  double q13(double s, double W, double theta, double thetaStar, double phiStar) const;
  double q23(double s, double W, double theta, double thetaStar, double phiStar) const;

protected:
  // mathematica wrapper functions
  vect Momentum(double s, Four_mom_type k, double W, double theta, double thetaStar, double phiStar) const;
  double Pair(vect k1, vect k2) const;
  double Power(double x, int n) const;
  double Power(vect x, int n) const;
  double Eps(vect k1, vect k2, vect k3, vect k4) const;

protected:
  // 4D Levi-Civita functionality
  constexpr bool hasDuplicateIndices(int i, int j, int k, int l) const;
  constexpr int getLeviCivitaValue(int i, int j, int k, int l) const;
};


namespace singlek {

class SingleKaonPhaseSpace : public TwoThreePhaseSpace {
  
public:
  enum class Channel {
    pp,
    np,
    nn,
  };

public:
  // prevent copying
  SingleKaonPhaseSpace(const SingleKaonPhaseSpace&) = delete;
  SingleKaonPhaseSpace& operator=(const SingleKaonPhaseSpace&) = delete;
  SingleKaonPhaseSpace(SingleKaonPhaseSpace&&) noexcept = default;
  SingleKaonPhaseSpace& operator=(SingleKaonPhaseSpace&&) noexcept = default;
  SingleKaonPhaseSpace(Channel channel, double m3_, double EbeamMin_, double EbeamMax_);
  ~SingleKaonPhaseSpace() { }

private:
  // constants
  static constexpr double GF{ 1.16639e-11 };
  static constexpr double D{ 0.804 };
  static constexpr double F{ 0.463 };
  static constexpr double Vus{ 0.22 };
  static constexpr double fPi{ 92.4 };
  static constexpr double mSigma{ 1192.46 };
  static constexpr double mLambda{ 1115.68 };
  static constexpr double mPi{ 134.98 };
  static constexpr double mEta{ 547.86 };
  static constexpr double kappaP{ 1.7928 };
  static constexpr double kappaN{ -1.9130 };
  static constexpr double mProton{ 938.272 };
  static constexpr double mNeutron{ 939.565 };
  static constexpr double mKaon0{ 497.648 };
  static constexpr double mKaonCharged{ 493.677 };

private:
  // members
  Channel m_channel;

  // form factors
  double ACT;
  double BCT;
  double ACRSigma;
  double ACRLambda;
  double AKP;
  double APi;
  double AEta;

private:
  // matrix elements
  double CT(double s, double W, double theta, double thetaStar, double phiStar) const;
  double CrossLambda(double s, double W, double theta, double thetaStar, double phiStar) const;
  double CrossSigma(double s, double W, double theta, double thetaStar, double phiStar) const;
  double PionInFlight(double s, double W, double theta, double thetaStar, double phiStar) const;
  double EtaInFlight(double s, double W, double theta, double thetaStar, double phiStar) const;
  double KaonPole(double s, double W, double theta, double thetaStar, double phiStar) const;

private:
  double matrixElement(double s, double W, double theta, double thetaStar, double phiStar) const override;
};


class SingleKaonPP : public SingleKaonPhaseSpace {

public:
  // prevent copying
  SingleKaonPP(const SingleKaonPP&) = delete;
  SingleKaonPP& operator=(const SingleKaonPP&) = delete;
  SingleKaonPP(SingleKaonPP&&) noexcept = default;
  SingleKaonPP& operator=(SingleKaonPP&&) noexcept = default;
  SingleKaonPP(double m3_, double EbeamMin_, double EbeamMax_) : SingleKaonPhaseSpace{ Channel::pp, m3_, EbeamMin_, EbeamMax_ } {}
  ~SingleKaonPP() { }

  static SingleKaonPP& construct(double m3_, double EbeamMin_, double EbeamMax_) {
    static SingleKaonPP instance{ m3_, EbeamMin_, EbeamMax_ };
    return instance;
  }

  static SingleKaonPP& get() {
    return construct(0.0, 0.0, 0.0);
  }
  
};

class SingleKaonNP : public SingleKaonPhaseSpace {

public:
  // prevent copying
  SingleKaonNP(const SingleKaonNP&) = delete;
  SingleKaonNP& operator=(const SingleKaonNP&) = delete;
  SingleKaonNP(SingleKaonNP&&) noexcept = default;
  SingleKaonNP& operator=(SingleKaonNP&&) noexcept = default;
  SingleKaonNP(double m3_, double EbeamMin_, double EbeamMax_) : SingleKaonPhaseSpace{ Channel::np, m3_, EbeamMin_, EbeamMax_ } {}
  ~SingleKaonNP() { }

  static SingleKaonNP& construct(double m3_, double EbeamMin_, double EbeamMax_) {
    static SingleKaonNP instance{ m3_, EbeamMin_, EbeamMax_ };
    return instance;
  }

  static SingleKaonNP& get() {
    return construct(0.0, 0.0, 0.0);
  }
  
};

class SingleKaonNN : public SingleKaonPhaseSpace {

public:
  // prevent copying
  SingleKaonNN(const SingleKaonNN&) = delete;
  SingleKaonNN& operator=(const SingleKaonNN&) = delete;
  SingleKaonNN(SingleKaonNN&&) noexcept = default;
  SingleKaonNN& operator=(SingleKaonNN&&) noexcept = default;
  SingleKaonNN(double m3_, double EbeamMin_, double EbeamMax_) : SingleKaonPhaseSpace{ Channel::nn, m3_, EbeamMin_, EbeamMax_ } {}
  ~SingleKaonNN() { }

  static SingleKaonNN& construct(double m3_, double EbeamMin_, double EbeamMax_) {
    static SingleKaonNN instance{ m3_, EbeamMin_, EbeamMax_ };
    return instance;
  }

  static SingleKaonNN& get() {
    return construct(0.0, 0.0, 0.0);
  }
  
};

}
}





