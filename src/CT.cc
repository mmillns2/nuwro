#include "matrix_elements.h"



namespace TwoThreeScatter {
namespace singlekaon {


double CT(double ma, double m1, double m2, double m3, vect vects[6]) {
  return (64*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(p1,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects))*Pair(Momentum(p2,m1,m2,m3,vects),Momentum(q1,m1,m2,m3,vects)))/
    Power(fPi,2) + (128*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(p1,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects))*
      Pair(Momentum(p2,m1,m2,m3,vects),Momentum(q1,m1,m2,m3,vects)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(p1,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects))*
      Pair(Momentum(p2,m1,m2,m3,vects),Momentum(q1,m1,m2,m3,vects)))/Power(fPi,2) - 
   (64*Power(ACT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(p2,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(p2,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(p1,m1,m2,m3,vects),Momentum(p2,m1,m2,m3,vects))*Pair(Momentum(q1,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects)))/
    Power(fPi,2) - (128*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(p1,m1,m2,m3,vects),Momentum(p2,m1,m2,m3,vects))*
      Pair(Momentum(q1,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(p1,m1,m2,m3,vects),Momentum(p2,m1,m2,m3,vects))*
      Pair(Momentum(q1,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects)))/Power(fPi,2); 
}


double CT_sem1(double ma, double m1, double m2, double m3, vect vects[6]) {	

	double a = Pair(Momentum(p1,m1,m2,m3,vects),Momentum(p2,m1,m2,m3,vects));
	double b = Pair(Momentum(p1,m1,m2,m3,vects),Momentum(q1,m1,m2,m3,vects));
	double c = Pair(Momentum(p1,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects));
	double d = Pair(Momentum(q1,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects));
	double e = Pair(Momentum(p2,m1,m2,m3,vects),Momentum(q3,m1,m2,m3,vects));
	double f = Pair(Momentum(q1,m1,m2,m3,vects),Momentum(p2,m1,m2,m3,vects));

  double front{ ACT*ACT * (Vus)*(Vus) * (1/((fPi)*(fPi))) }; 
	return front * ( 64*(a*d + f*c) + 64*BCT*(a*d + c*f + b*e - (3*(e*(b + m1*ma)))) );
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
