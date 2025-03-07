#include "matrix_elements.h"



namespace TwoThreeScatter {
namespace singlekaon {


double CT(vect vects[6]) {
  return 0;
}

  
double CT(double ma, double m1, double m2, double m3, double s, double W, double theta, double thetaStar, double phiStar) {
  return (64*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(p1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),
       Momentum(q3,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar))*
      Pair(Momentum(p2,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),Momentum(q1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
   (128*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(p1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),
       Momentum(q3,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar))*
      Pair(Momentum(p2,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),Momentum(q1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(p1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),
       Momentum(q3,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar))*
      Pair(Momentum(p2,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),Momentum(q1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar)))/Power(fPi,2) - 
   (64*Power(ACT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(p2,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),
       Momentum(q3,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(p2,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),
       Momentum(q3,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(p1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),
       Momentum(p2,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar))*
      Pair(Momentum(q1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),Momentum(q3,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar)))/Power(fPi,2) - 
   (128*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(p1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),
       Momentum(p2,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar))*
      Pair(Momentum(q1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),Momentum(q3,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(p1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),
       Momentum(p2,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar))*
      Pair(Momentum(q1,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar),Momentum(q3,ma,m1,m2,m3,s,W,theta,thetaStar,phiStar)))/Power(fPi,2);
}

}
}
