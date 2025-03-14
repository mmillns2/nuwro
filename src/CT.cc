#include "matrix_elements.h"



namespace TwoThreeScatter {
namespace singlekaon {


double CT(double ma, double m1, double m2, double m3, vect vects[6]) {
  return (64*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)))/Power(fPi,2) + 
   (128*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)))/Power(fPi,2) - 
   (64*Power(ACT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(p2,vects),Momentum(q3,vects)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(p2,vects),Momentum(q3,vects)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))/Power(fPi,2) - 
   (128*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))/Power(fPi,2) + 
   (64*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))/Power(fPi,2); 
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


double CrossLambda(double ma, double m1, double m2, double m3, vect vects[6]) {
  return (Power(ACRLambda,2)*Power(Vus,2)*((32*(D + 3*F)*(-((2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
                 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
               (-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) + 4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                 (2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                  Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 2*
             ((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p1,vects),Momentum(q2,vects)) + 
               (ma*mLambda + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
             (2*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/(3.*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2))
         - (32*Power(D + 3*F,2)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + m1*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(2*ma*mLambda*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Power(mLambda,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*((ma + 2*mLambda)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  ma*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
               2*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               m1*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/
        (9.*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2)) - 
       (32*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - m1*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(-2*ma*mLambda*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Power(mLambda,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               2*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*((ma + 2*mLambda)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  ma*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
               2*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/
        (Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2)) - (16*Power(D + 3*F,2)*
          (-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             ((m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(q,vects),Momentum(q1,vects))*
                ((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               Pair(Momentum(q,vects),Momentum(q,vects))*((-Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) + (-(m1*ma) + m1*mLambda + 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(-2*Pair(Momentum(q,vects),Momentum(q1,vects))*
                (-2*(ma*mLambda + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*(2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                     Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
               Pair(Momentum(q,vects),Momentum(q,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                   (m1*ma + 2*m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p1,vects),Momentum(q1,vects))*(mLambda*(2*ma + mLambda) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  ma*(-2*mLambda*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))))*
          (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/
        (9.*Power(Power(m2,2) - Power(q,2),2)*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2)) + 
       (8*(D + 3*F)*kappaP*(m1*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)))*
             (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - 4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(mLambda*(2*ma + mLambda)) - Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  Pair(Momentum(q2,vects),Momentum(q2,vects)))) + (2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
             (-4*mLambda*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) + 
               2*(-ma + mLambda)*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Pair(Momentum(q2,vects),Momentum(q2,vects))*((ma + 2*mLambda)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  ma*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))) - 
            2*m1*((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p1,vects),Momentum(q2,vects)) + 
               (ma*mLambda + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
             (2*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
            m1*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - 4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(mLambda*(2*ma + mLambda)) - Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  Pair(Momentum(q2,vects),Momentum(q2,vects))))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            (-4*mLambda*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) + 
               2*(-ma + mLambda)*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Pair(Momentum(q2,vects),Momentum(q2,vects))*((ma + 2*mLambda)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  ma*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))*
             (2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
            2*m1*((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p1,vects),Momentum(q2,vects)) + 
               (ma*mLambda + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
             (2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
        (3.*ma*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2)) + 
       (32*Power(D + 3*F,2)*(-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
               ((m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
                ((-Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects))*
                   Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(q,vects),Momentum(q3,vects))*
                   ((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                     (m1*ma - m1*mLambda - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
                  Power(mLambda,2)*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
                  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               ((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                   Pair(Momentum(q,vects),Momentum(q1,vects)) - Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q,vects))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) - m1*ma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  m1*mLambda*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p2,vects),Momentum(q1,vects))*((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                      Pair(Momentum(q,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))))*
                Pair(Momentum(q3,vects),Momentum(q3,vects))) + Pair(Momentum(q2,vects),Momentum(q2,vects))*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
                (mLambda*(2*ma + mLambda) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(-(Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                     (mLambda*(2*ma + mLambda) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
                  Pair(Momentum(q,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (m1*ma + 2*m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mLambda*(2*ma + mLambda) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mLambda*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))) + 
                  2*(ma*mLambda + Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   (Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
                     Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
               (-2*ma*mLambda*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
                  Power(mLambda,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                   (mLambda*(2*ma + mLambda) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p2,vects),Momentum(q,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (m1*ma + 2*m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mLambda*(2*ma + mLambda) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mLambda*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)))))/(9.*(Power(m2,2) - Power(q,2))*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2)) + 
       (16*kappaP*(Pair(Momentum(q2,vects),Momentum(q2,vects))*(-2*ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) - 4*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) - 4*m1*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) - 4*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) - ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) - ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) - 2*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) + 4*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) + m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)))*
                (2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
               2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*m1*ma*mLambda*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               4*m1*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               4*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                  Pair(Momentum(q2,vects),Momentum(q2,vects)))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               2*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               2*m1*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(-2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               4*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               mLambda*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*m1*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               m1*(Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                (2*Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
            4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects))) - 
               mLambda*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
        (ma*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2)) + 
       (16*(D + 3*F)*kappaP*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
             (ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               2*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               2*m1*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               2*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
                (2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*m1*ma*mLambda*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*m1*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                  Pair(Momentum(q2,vects),Momentum(q2,vects)))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
               ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               2*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               2*m1*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               mLambda*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               m1*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               m1*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               m1*(Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
               mLambda*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
        (3.*ma*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2)) - 
       (4*Power(kappaP,2)*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
             (4*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*m1*ma*Power(mLambda,2)*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*m1*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*ma*mLambda*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Power(mLambda,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2)*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + 4*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   ((ma + 2*mLambda)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + ma*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
                  2*Pair(Momentum(q1,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*
                      (2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
                     2*(ma*mLambda + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mLambda + Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
               4*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects))) - 4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*m1*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + m1*mLambda*Pair(Momentum(q3,vects),Momentum(q3,vects))))
              + 2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
                (m1*(-ma + mLambda)*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                     (-Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
               2*Pair(Momentum(p2,vects),Momentum(q1,vects))*((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(2*m1*(ma - mLambda)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) - 4*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Power(mLambda,2)*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - m1*ma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  m1*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
        (Power(ma,2)*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2)) - 
       (4*Power(kappaP,2)*(-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*((2*m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)))*
                   Pair(Momentum(p2,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
                  2*(m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*m1*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
                  (2*m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(4*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(mLambda,2)*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*m1*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*ma*mLambda*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Power(mLambda,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2)*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + 4*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (mLambda*(2*ma + mLambda) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
               4*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*ma*Power(mLambda,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 4*m1*mLambda*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*ma*mLambda*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + Power(mLambda,2)*Pair(Momentum(p1,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) - 2*ma*mLambda*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               2*m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
                (2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   (mLambda*(2*ma + mLambda) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (2*m1*ma + 4*m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mLambda*(2*ma + mLambda) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     2*ma*(-(mLambda*Pair(Momentum(q1,vects),Momentum(q2,vects))) + m1*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))) - 
                  4*(ma*mLambda + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*
                   (Pair(Momentum(p1,vects),Momentum(p1,vects))*(m1*ma + 2*m1*mLambda + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 
                        2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mLambda*(2*ma + mLambda) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mLambda*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mLambda,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                ((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(2*m1*(ma - mLambda)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) - 4*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Power(mLambda,2)*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - Power(mLambda,2)*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)) - 2*m1*ma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  2*m1*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   ((-Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                     2*(m1*(-ma + mLambda) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  2*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                     (-Power(mLambda,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                      Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
                     ((Power(mLambda,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                        (m1*ma - m1*mLambda - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                      Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(Power(ma,2)*Power(mLambda + p1 - q2,2)*Power(mLambda - p1 + q2,2))))/(8.*Power(fPi,2));
}


double PionInFlight(double ma, double m1, double m2, double m3, vect vects[6]) {
  return (8*Power(APi,2)*Power(D + F,2)*Power(ma,2)*Power(Vus,2)*(m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*
     (Pair(Momentum(p2,vects),Momentum(q3,vects))*(-4*Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
       4*Pair(Momentum(p2,vects),Momentum(q2,vects))*(2*Pair(Momentum(q2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
       Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
          2*(2*Pair(Momentum(q2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
   (Power(fPi,2)*Power(mPi + q - q2,2)*Power(mPi - q + q2,2));
}

double EtaInFlight(double ma, double m1, double m2, double m3, vect vects[6]) {
  return (8*Power(AEta,2)*Power(D - 3*F,2)*Power(ma,2)*Power(Vus,2)*(m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*
     (Pair(Momentum(p2,vects),Momentum(q3,vects))*(-4*Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
       4*Pair(Momentum(p2,vects),Momentum(q2,vects))*(2*Pair(Momentum(q2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
       Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
          2*(2*Pair(Momentum(q2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
   (Power(fPi,2)*Power(mEta + q - q2,2)*Power(mEta - q + q2,2));
}

double CrossSigma(double ma, double m1, double m2, double m3, vect vects[6]) {
  return (Power(ACRSigma,2)*Power(Vus,2)*((32*(D - F)*((2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
             (-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) + 4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                Pair(Momentum(q2,vects),Momentum(q2,vects))) + 2*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                Pair(Momentum(p1,vects),Momentum(q2,vects)) + (ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
             (2*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/((Power (p1 - q2, 2) - Power (mSigma, 2))) - 
       (32*Power(D - F,2)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(2*ma*mSigma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/
        ((Power (p1 - q2, 2) - Power (mSigma, 2))) - (32*
          (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
               m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(-2*ma*mSigma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/
        ((Power (p1 - q2, 2) - Power (mSigma, 2))) - (16*Power(D - F,2)*
          (-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             ((m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(q,vects),Momentum(q1,vects))*
                ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               Pair(Momentum(q,vects),Momentum(q,vects))*((-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) + (-(m1*ma) + m1*mSigma + 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(-2*Pair(Momentum(q,vects),Momentum(q1,vects))*
                (-2*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*(2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                     Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
               Pair(Momentum(q,vects),Momentum(q,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                   (m1*ma + 2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  ma*(-2*mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))))*
          (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/
        (Power(Power(m2,2) - Power(q,2),2)*(Power (p1 - q2, 2) - Power (mSigma, 2))) + 
       (8*(2*kappaN + kappaP)*(-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects)) + 
               mSigma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               (m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
                ((-ma + mSigma)*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  m1*(-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               (m1*(Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
                  (-2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + (ma - mSigma)*Pair(Momentum(p2,vects),Momentum(q1,vects)))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   ((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
                  2*m1*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
               (-2*m1*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
                  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                     Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p2,vects),Momentum(q1,vects))*((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                     ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
        (ma*(Power (p1 - q2, 2) - Power (mSigma, 2))) + 
       (8*(2*kappaN + kappaP)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects)) + 
               mSigma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               (m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
                ((-ma + mSigma)*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  m1*(-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               (m1*(Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
                  (-2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + (ma - mSigma)*Pair(Momentum(p2,vects),Momentum(q1,vects)))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                  (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   ((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
                  2*m1*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               (-2*m1*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
                  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                     Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p2,vects),Momentum(q1,vects))*((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                     ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
        (ma*(Power (p1 - q2, 2) - Power (mSigma, 2))) + 
       (32*Power(D - F,2)*(-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
               ((m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
                ((-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects))*
                   Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(q,vects),Momentum(q3,vects))*
                   ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                     (m1*ma - m1*mSigma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
                  Power(mSigma,2)*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
                  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                   Pair(Momentum(q,vects),Momentum(q1,vects)) - Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q,vects))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) - m1*ma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  m1*mSigma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p2,vects),Momentum(q1,vects))*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                      Pair(Momentum(q,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))))*
                Pair(Momentum(q3,vects),Momentum(q3,vects))) + Pair(Momentum(q2,vects),Momentum(q2,vects))*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
                (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(-(Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                     (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
                  Pair(Momentum(q,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (m1*ma + 2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))) + 
                  2*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   (Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
                     Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
               (-2*ma*mSigma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
                  Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                   (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p2,vects),Momentum(q,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (m1*ma + 2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)))))/((Power(m2,2) - Power(q,2))*(Power (p1 - q2, 2) - Power (mSigma, 2))) + 
       (4*Power(2*kappaN + kappaP,2)*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
             (4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*m1*ma*Power(mSigma,2)*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*m1*mSigma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*ma*mSigma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2)*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + 4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   ((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
                  2*Pair(Momentum(q1,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*
                      (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
                     2*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
               4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*mSigma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects))) - 4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + m1*mSigma*Pair(Momentum(q3,vects),Momentum(q3,vects))))\
             + 2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
                (m1*(-ma + mSigma)*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                     (-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
               2*Pair(Momentum(p2,vects),Momentum(q1,vects))*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(2*m1*(ma - mSigma)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) - 4*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Power(mSigma,2)*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - m1*ma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  m1*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
        (Power(ma,2)*(Power (p1 - q2, 2) - Power (mSigma, 2))) + 
       (4*Power(2*kappaN + kappaP,2)*(-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*((2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
                  2*(m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
                  (2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(mSigma,2)*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*m1*mSigma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*ma*mSigma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2)*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + 4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
               4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*ma*Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 4*m1*mSigma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*ma*mSigma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (2*m1*ma + 4*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     2*ma*(-(mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects))) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))) - 
                  4*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*
                   (Pair(Momentum(p1,vects),Momentum(p1,vects))*(m1*ma + 2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 
                        2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(2*m1*(ma - mSigma)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) - 4*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Power(mSigma,2)*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - Power(mSigma,2)*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)) - 2*m1*ma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  2*m1*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   ((-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                     2*(m1*(-ma + mSigma) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  2*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                     (-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                      Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
                     ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                        (m1*ma - m1*mSigma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                      Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(Power(ma,2)*(Power (p1 - q2, 2) - Power (mSigma, 2)))))/(8.*Power(fPi,2));
}

double KaonPole(double ma, double m1, double m2, double m3, vect vects[6]) {
  return (2*Power(AKP,2)*Power(Vus,2)*((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
       2*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 2*m1*ma*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
       2*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
       2*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
       2*Pair(Momentum(p1,vects),Momentum(q,vects))*(Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
       m1*ma*Pair(Momentum(q2,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
     (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
       Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(Power(fPi,2)*Power(Power(m2,2) - Power(q,2),2));
}

}
}
