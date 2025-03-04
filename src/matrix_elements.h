#ifndef SINGLE_KAON_MATRIX_ELEMENTS
#define SINGLE_KAON_MATRIX_ELEMENTS


#include "23Scatter.h"


namespace TwoThreeScatter {
namespace singlekaon { 

// constants
constexpr double D{ 0.804 };
constexpr double F{ 0.463 };
constexpr double Vus{ 0.22 };
constexpr double fPi{ 92.4e-3 };
constexpr double mSigma{  };
constexpr double mLambda{  };
constexpr double mPi{  };
constexpr double mEta{  };
constexpr double kappaP{ 1.7928 };
constexpr double kappaN{ -1.9130 };

// form factors
inline double ACT;
inline double BCT;
inline double ACRSigma;
inline double ACRLambda;
inline double AKP;
inline double APi;
inline double AEta;


// matrix elements
double CT(double ma, double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar);

}
}


#endif
