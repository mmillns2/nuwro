#include "matrix_elements.h"

namespace TwoThreeScatter
{
namespace singlekaon
{
  double CT(double ma, double m1, double m2, double m3, vect vects[6])
  {
      // Not included GF squared factor here; have included the missing 1/4 factor
return (16*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)))/Power(fPi,2) + 
(32*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)))/Power(fPi,2) + 
(16*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)))/Power(fPi,2) - 
(16*Power(ACT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(p2,vects),Momentum(q3,vects)))/Power(fPi,2) + 
(16*Power(ACT,2)*Power(BCT,2)*m1*ma*Power(Vus,2)*Pair(Momentum(p2,vects),Momentum(q3,vects)))/Power(fPi,2) + 
(16*Power(ACT,2)*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))/Power(fPi,2) - 
(32*Power(ACT,2)*BCT*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))/Power(fPi,2) + 
(16*Power(ACT,2)*Power(BCT,2)*Power(Vus,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))/Power(fPi,2);
  }

  double CrossSigma(double ma, double m1, double m2, double m3, vect vects[6])
   {
      // Not included GF squared factor here; have included the missing 1/4 factor
return (Power(ACRSigma,2)*Power(Vus,2)*(32*Power(D - F,2)*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2) + 32*
(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2) + 
32*Power(mSigma,2)*(-((Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   Pair(Momentum(q2,vects),Momentum(q2,vects))) + 2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
32*Power(D - F,2)*Power(mSigma,2)*(-((Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   Pair(Momentum(q2,vects),Momentum(q2,vects))) + 2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
32*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
 (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
   m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
 Pair(Momentum(q2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
32*Power(D - F,2)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
 (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
   m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
 Pair(Momentum(q2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
(16*Power(D - F,2)*((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2) + 
(16*Power(D - F,2)*Power(mSigma,2)*(Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (4*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
     2*Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
  ((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) - 
     2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2) - 
(16*Power(D - F,2)*(Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q,vects))*
   (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
   (-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p1,vects),Momentum(q2,vects))*
      (4*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
        2*Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))) - 
     m1*ma*Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2) + 
(2*Power(2*kappaN + kappaP,2)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
      (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) - 
(2*Power(2*kappaN + kappaP,2)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) + 
(2*Power(2*kappaN + kappaP,2)*Power(mSigma,2)*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q3,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*
           (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects)))) - 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(q2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
     Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) + 
(2*Power(2*kappaN + kappaP,2)*(Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q3,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*
           (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects)))) - 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) - 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(q2,vects),Momentum(q3,vects))*
         (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
     Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) + (2*Power(2*kappaN + kappaP,2)*Power(mSigma,2)*
(Pair(Momentum(q2,vects),Momentum(q2,vects))*(-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) - (2*Power(2*kappaN + kappaP,2)*(-(Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
     (-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
       2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
       Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
          Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
       2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
       m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
       Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
       Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
          Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
          Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))) + 
  4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(q1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2)))/(8.*Power(fPi,2)*Power(-Power(mSigma,2) + Pair(Momentum(p1,vects)-Momentum(q2,vects),Momentum(p1,vects)-Momentum(q2,vects)),2)) + 
(Power(ACRSigma,2)*Power(Vus,2)*(64*(D - F)*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2) - 
64*mSigma*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
(ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
   m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
   ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
64*Power(D - F,2)*mSigma*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
(ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
   m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
   ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
64*Power(D - F,2)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
64*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
64*(D - F)*Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
64*(D - F)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
(-(m1*Pair(Momentum(p1,vects),Momentum(p2,vects))) + ma*Pair(Momentum(p2,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
64*(D - F)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
(m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + ma*Pair(Momentum(p2,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
64*(D - F)*mSigma*(ma*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
32*(D - F)*mSigma*(2*ma*(-(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))) + 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
4*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
32*(D - F)*Pair(Momentum(q2,vects),Momentum(q2,vects))*(-4*Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
2*m1*ma*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
2*Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
64*(D - F)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
 (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
Pair(Momentum(p1,vects),Momentum(p1,vects))*(-(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))) + 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) - 
64*Pair(Momentum(q2,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
   2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
64*Power(D - F,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p1,vects),Momentum(q2,vects))*(-2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
   m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
32*(D - F)*Power(mSigma,2)*(2*(-(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))) + 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
Pair(Momentum(p1,vects),Momentum(q2,vects))*(-4*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   4*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
(16*(2*kappaN + kappaP)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*mSigma*(2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))) - 
  2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   ((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*(-4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) + m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q2,vects)) + Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*Power(mSigma,2)*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q2,vects)) - Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))))/ma + 
(16*(2*kappaN + kappaP)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))))/ma - 
(16*(2*kappaN + kappaP)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     ma*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
  ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma + 
(16*(2*kappaN + kappaP)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma + 
(32*(D - F)*(2*kappaN + kappaP)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma - 
(32*Power(D - F,2)*mSigma*(m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
  2*ma*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
  ma*Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2) - 
(32*Power(D - F,2)*(Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
     4*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(-2*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
     Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2) + 
(32*Power(D - F,2)*mSigma*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(q,vects),Momentum(q,vects)) - 
  (m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
     ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q,vects)) - 
     2*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2) + 
(32*Power(D - F,2)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  ((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/
(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) + (16*(2*kappaN + kappaP)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma + 
(16*(2*kappaN + kappaP)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma - 
(16*(2*kappaN + kappaP)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma - 
(4*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(4*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects)) + 
  2*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
  4*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  4*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  4*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma + 
(4*Power(2*kappaN + kappaP,2)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(ma,2) - 
(4*Power(2*kappaN + kappaP,2)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(ma,2) + 
(16*(D - F)*(2*kappaN + kappaP)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  ma*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))/ma - 
(4*Power(2*kappaN + kappaP,2)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*((-2*m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     2*(m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*(2*m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     (-2*m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) - 
(32*Power(D - F,2)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) - 
(32*Power(D - F,2)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(-(m1*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))) + 
     m1*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) + 
(4*Power(2*kappaN + kappaP,2)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(q2,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) - 
(32*Power(D - F,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
   Pair(Momentum(q2,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
      ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) + 
        2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) + 
        2*(Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
           Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) - (32*Power(D - F,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
   Pair(Momentum(q2,vects),Momentum(q3,vects)) - m1*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
      ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) + 
        2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) + 
        2*(Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
           Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) + (32*Power(D - F,2)*Power(mSigma,2)*
(Pair(Momentum(q2,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) - Pair(Momentum(p2,vects),Momentum(p2,vects))*
      ((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     ((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) - 
        Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
      (Pair(Momentum(q,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
        Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
     (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) + (4*Power(2*kappaN + kappaP,2)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(q1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
      (Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
Power(ma,2) + (32*(D - F)*(2*kappaN + kappaP)*(Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*Power(mSigma,2)*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*(4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma - 
(16*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma - 
(32*(D - F)*(2*kappaN + kappaP)*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(32*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma - 
(32*(D - F)*(2*kappaN + kappaP)*Power(mSigma,2)*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(4*Power(2*kappaN + kappaP,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
   Pair(Momentum(q2,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q3,vects))*
   (-4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
        4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*((-(m1*ma) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        (m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        (-(m1*ma) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) - 
(4*Power(2*kappaN + kappaP,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q3,vects))*
   (4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
        4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        (-(m1*ma) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        (-(m1*ma) + 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) + 
(32*Power(D - F,2)*(Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (-(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     (-((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(q,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
        Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
     (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q1,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) - 
(16*(2*kappaN + kappaP)*mSigma*(2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*(2*kappaN + kappaP)*mSigma*(2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(4*Power(2*kappaN + kappaP,2)*Power(mSigma,2)*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*((2*m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        2*(m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
        2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        (2*m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) - 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
        2*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        2*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) - 
(4*Power(2*kappaN + kappaP,2)*mSigma*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) - (4*Power(2*kappaN + kappaP,2)*mSigma*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
        (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
     ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) - (4*Power(2*kappaN + kappaP,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p1,vects),Momentum(q2,vects))*(4*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
        2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     4*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     2*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
      (-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) + (4*Power(2*kappaN + kappaP,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p1,vects),Momentum(q2,vects))*(-4*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      Pair(Momentum(p2,vects),Momentum(q3,vects)) + 2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     4*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     2*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) - (32*(D - F)*(2*kappaN + kappaP)*mSigma*(-2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   ((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/ma - 
(32*(D - F)*(2*kappaN + kappaP)*mSigma*(2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   ((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/ma - 
(16*(2*kappaN + kappaP)*mSigma*(2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   ((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/ma - 
(32*Power(D - F,2)*mSigma*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*
      Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q,vects))*
      (-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*
      (-(m1*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))) + 
        m1*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     ma*(Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
           Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) + 
(32*Power(D - F,2)*mSigma*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*
      Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q,vects))*
      (-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*
      (-(m1*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))) + 
        m1*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     ma*(Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
           Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects))) + 
(4*Power(2*kappaN + kappaP,2)*mSigma*(-2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
            (Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/Power(ma,2) + 
(4*Power(2*kappaN + kappaP,2)*mSigma*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
            (Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/Power(ma,2) - 
(4*Power(2*kappaN + kappaP,2)*(Pair(Momentum(p1,vects),Momentum(q1,vects))*
   (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  2*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
           Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - m1*ma*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
            (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
              2*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
              2*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
              Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))))/Power(ma,2)))/
(8.*Power(fPi,2)*Power(-Power(mSigma,2) + Pair(Momentum(p1,vects)-Momentum(q2,vects),Momentum(p1,vects)-Momentum(q2,vects)),2));
  }

  double CrossLambda(double ma, double m1, double m2, double m3, vect vects[6])
  {
      // Not included GF squared factor here; have included the missing 1/4 factor
return (Power(ACRLambda,2)*Power(Vus,2)*((32*Power(D + 3*F,2)*
(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2))/9. + 
32*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2) + 
32*Power(mLambda,2)*(-((Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   Pair(Momentum(q2,vects),Momentum(q2,vects))) + 2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
(32*Power(D + 3*F,2)*Power(mLambda,2)*(-((Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
       m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
     Pair(Momentum(q2,vects),Momentum(q2,vects))) + 2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/9. + 
32*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
 (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
   Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
   m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
 Pair(Momentum(q2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
(32*Power(D + 3*F,2)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   Pair(Momentum(q2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/9. + 
(16*Power(D + 3*F,2)*((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(9.*Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2)) + 
(16*Power(D + 3*F,2)*Power(mLambda,2)*(Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (4*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
     2*Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
  ((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) - 
     2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(9.*Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2)) - 
(16*Power(D + 3*F,2)*(Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q,vects))*
   (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
   (-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p1,vects),Momentum(q2,vects))*
      (4*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
        2*Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))) - 
     m1*ma*Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(9.*Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2)) + 
(2*Power(kappaP,2)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
      (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) + 
(2*Power(kappaP,2)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) + 
(2*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q3,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*
           (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects)))) - 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(q2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
     Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) + 
(2*Power(kappaP,2)*(Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q3,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*
           (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects)))) - 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) - 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(q2,vects),Momentum(q3,vects))*
         (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
     Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) + (2*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) + (2*Power(kappaP,2)*(Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (-2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) - 
  4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))) - 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(q1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2)))/(32*Power(fPi,2)*Power(-Power(mLambda,2) + Pair(Momentum(p1,vects)-Momentum(q2,vects),Momentum(p1,vects)-Momentum(q2,vects)),2)) + 
(Power(ACRLambda,2)*Power(Vus,2)*((-64*(D + 3*F)*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2))/3. - 
64*mLambda*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
(ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
   m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
   ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
(64*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  (ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
     m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))))/9. + 
(64*Power(D + 3*F,2)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))))/9. + 
64*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
(64*(D + 3*F)*Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))))/3. - 
(64*(D + 3*F)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  (-(m1*Pair(Momentum(p1,vects),Momentum(p2,vects))) + ma*Pair(Momentum(p2,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))))/3. + 
(64*(D + 3*F)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  (m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + ma*Pair(Momentum(p2,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))))/3. + 
(64*(D + 3*F)*mLambda*(ma*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/3. + 
(64*(D + 3*F)*mLambda*(ma*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(-2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/3. + 
(32*(D + 3*F)*Pair(Momentum(q2,vects),Momentum(q2,vects))*(4*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  2*m1*ma*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/3. + 
(64*(D + 3*F)*Power(mLambda,2)*((Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/3. - 
(64*(D + 3*F)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(-(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/3. - 
(64*Power(D + 3*F,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/9. - 
64*Pair(Momentum(q2,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q2,vects))*
 (2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
   2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
   Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
(16*kappaP*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*mLambda*(2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))) - 
  2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   ((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*(-4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) + m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q2,vects)) + Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q2,vects)) - Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))))/ma + 
(16*kappaP*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))))/ma + 
(16*kappaP*Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
     Pair(Momentum(p2,vects),Momentum(q3,vects))) + m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) - ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q1,vects),Momentum(q3,vects)) + ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(q3,vects))*
   (ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*(-2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma + 
(16*kappaP*Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
     Pair(Momentum(p2,vects),Momentum(q3,vects))) + m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) - ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q1,vects),Momentum(q3,vects)) + ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q3,vects))*
   (ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     ma*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma - 
(32*Power(D + 3*F,2)*mLambda*(m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
  2*ma*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
  ma*Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(9.*Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2)) - 
(32*Power(D + 3*F,2)*(Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
     4*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(-2*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
     Pair(Momentum(q,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(9.*Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2)) + 
(32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(q,vects),Momentum(q,vects)) - 
  (m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
     ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q,vects)) - 
     2*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
(Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(9.*Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2)) + 
(32*Power(D + 3*F,2)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  ((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/
(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) + (16*kappaP*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma + 
(16*kappaP*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/ma - 
(4*(D + 3*F)*kappaP*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-4*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects)) - 
  2*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
  4*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  4*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  4*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(3.*ma) - 
(4*Power(kappaP,2)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/Power(ma,2) - 
(4*Power(kappaP,2)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  2*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) + 
(16*(D + 3*F)*kappaP*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
  ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))/(3.*ma) - 
(4*Power(kappaP,2)*Power(Pair(Momentum(q2,vects),Momentum(q2,vects)),2)*
(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*((-2*m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     2*(m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(q3,vects))*(2*m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     (-2*m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) - 
(32*Power(D + 3*F,2)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) - 
(32*Power(D + 3*F,2)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(-(m1*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))) + 
     m1*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) + 
(4*Power(kappaP,2)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(q2,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/Power(ma,2) - 
(32*Power(D + 3*F,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
   Pair(Momentum(q2,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
      ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) + 
        2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) + 
        2*(Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
           Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) - (32*Power(D + 3*F,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
   Pair(Momentum(q2,vects),Momentum(q3,vects)) - m1*ma*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
      ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) + 
        2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     ((m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) + 
        2*(Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
           Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) + (32*Power(D + 3*F,2)*Power(mLambda,2)*
(Pair(Momentum(q2,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) - Pair(Momentum(p2,vects),Momentum(p2,vects))*
      ((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     ((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) - 
        Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
      (Pair(Momentum(q,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
        Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
     (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) + (4*Power(kappaP,2)*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(q1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
      (Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
     2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
Power(ma,2) - (32*(D + 3*F)*kappaP*(Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/(3.*ma) + 
(16*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*(4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) - ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) + ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(32*(D + 3*F)*kappaP*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
(3.*ma) + (32*(D + 3*F)*kappaP*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
(3.*ma) + (32*(D + 3*F)*kappaP*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/(3.*ma) + 
(32*(D + 3*F)*kappaP*mLambda*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(q2,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/(3.*ma) + 
(32*(D + 3*F)*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   ((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     ma*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p2,vects),Momentum(q1,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/(3.*ma)
+ (4*Power(kappaP,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
   Pair(Momentum(q2,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q3,vects))*
   (-4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
        4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*((-(m1*ma) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        (m1*ma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        (-(m1*ma) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) - 
(4*Power(kappaP,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(-2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
  m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
   Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))*
   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q3,vects))*
   (4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
        4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        (-(m1*ma) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        (-(m1*ma) + 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) + 
(32*Power(D + 3*F,2)*(Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (-(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     (-((m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + 
        Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(q,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) - 
        Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
     (Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q1,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) + 
(4*Power(kappaP,2)*mLambda*(-2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*
      (Pair(Momentum(q1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
        2*Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
     m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*
           (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects)))) - 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) - 
(4*Power(kappaP,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p1,vects),Momentum(q2,vects))*(4*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     4*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(q3,vects))*(4*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        2*Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (Pair(Momentum(q2,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
        2*Pair(Momentum(p2,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
     Pair(Momentum(q1,vects),Momentum(q3,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*
           (Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(q2,vects),Momentum(q3,vects)))) - 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) - 
(16*kappaP*mLambda*(2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(16*kappaP*mLambda*(2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (-Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/ma + 
(4*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
   (2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*((2*m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
        2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        2*(m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
        2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        (2*m1*ma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) - 
  2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
      (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
        2*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
        2*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/Power(ma,2) + 
(4*Power(kappaP,2)*mLambda*(-2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - 2*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q1,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) - (4*Power(kappaP,2)*Pair(Momentum(q2,vects),Momentum(q2,vects))*
(Pair(Momentum(p1,vects),Momentum(q2,vects))*(4*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) + 
     Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
        2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
     4*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
     2*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
        2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
        Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(p1,vects),Momentum(p1,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
      (-2*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(-2*Pair(Momentum(p2,vects),Momentum(q2,vects))*
         (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(q2,vects),Momentum(q3,vects))*(-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/
Power(ma,2) + (32*(D + 3*F)*kappaP*mLambda*(-2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   ((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(3.*ma) + 
(32*(D + 3*F)*kappaP*mLambda*(2*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p2,vects),Momentum(q1,vects)) - Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   ((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects))) + 
     m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*(-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(3.*ma) - 
(16*kappaP*mLambda*(2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
  2*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   ((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p2,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     m1*ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/ma - 
(32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*
      Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q,vects))*
      (-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*
      (-(m1*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))) + 
        m1*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     ma*(Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
           Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) + 
(32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
   (Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
     Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q2,vects))*
      Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(q,vects))*
      (-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(ma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
      Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*
      (-(m1*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects))) + 
        m1*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
     ma*(Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
        Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
           Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q,vects))*(-(Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(9.*(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)))) + 
(4*Power(kappaP,2)*mLambda*(-2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
            (Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/Power(ma,2) + 
(4*Power(kappaP,2)*mLambda*(2*m1*Pair(Momentum(p1,vects),Momentum(q2,vects))*
   (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
  4*m1*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
   (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(q2,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(p1,vects))*
      (Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     ma*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
         (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q1,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
            (Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/Power(ma,2) + 
(4*Power(kappaP,2)*(-(Pair(Momentum(p1,vects),Momentum(q1,vects))*
     (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
     (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
       Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
  2*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*
      (-Pair(Momentum(p2,vects),Momentum(p2,vects)) + Pair(Momentum(p2,vects),Momentum(q3,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
     Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
      (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) - Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
      (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
     Pair(Momentum(p1,vects),Momentum(p1,vects))*((Pair(Momentum(p2,vects),Momentum(p2,vects)) - Pair(Momentum(p2,vects),Momentum(q3,vects)))*
         Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
           Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - m1*ma*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
        Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
            (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
              2*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
              2*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
           Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
              Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))))/Power(ma,2)))/
(32*Power(fPi,2)*Power(-Power(mLambda,2) + Pair(Momentum(p1,vects)-Momentum(q2,vects),Momentum(p1,vects)-Momentum(q2,vects)),2));
  }

  double KaonPole(double ma, double m1, double m2, double m3, vect vects[6])
  {
    // Not included GF squared factor here; have included the missing 1/4 factor
     return (Power(AKP,2)*Power(Vus,2)*((m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) + 
       2*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 2*m1*ma*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
       2*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
       2*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
       2*Pair(Momentum(p1,vects),Momentum(q,vects))*(Pair(Momentum(q,vects),Momentum(q1,vects)) + Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
       m1*ma*Pair(Momentum(q2,vects),Momentum(q2,vects)) - Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
     (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
       Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/(2*Power(fPi,2)*Power(Power(m2,2) - Pair(Momentum(q,vects),Momentum(q,vects)),2));
  }

  double PionInFlight(double ma, double m1, double m2, double m3, vect vects[6])
  {
    // Not included GF squared factor here; have included the missing 1/4 factor
return (-2*Power(APi,2)*Power(D + F,2)*Power(ma,2)*Power(Vus,2)*(m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*
(Pair(Momentum(p2,vects),Momentum(q3,vects))*(-4*Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  4*Pair(Momentum(p2,vects),Momentum(q2,vects))*(2*Pair(Momentum(q2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*(2*Pair(Momentum(q2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/(Power(fPi,2)*Power(Power(mPi,2) - Pair(Momentum(q,vects)-Momentum(q2,vects),Momentum(q,vects)-Momentum(q2,vects)),2));
  }

  double EtaInFlight(double ma, double m1, double m2, double m3, vect vects[6])
  {
     // Not included GF squared factor here; have included the missing 1/4 factor
return (-2*Power(AEta,2)*Power(D - 3*F,2)*Power(ma,2)*Power(Vus,2)*(m1*ma - Pair(Momentum(p1,vects),Momentum(q1,vects)))*
(Pair(Momentum(p2,vects),Momentum(q3,vects))*(-4*Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  4*Pair(Momentum(p2,vects),Momentum(q2,vects))*(2*Pair(Momentum(q2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
  Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
     2*(2*Pair(Momentum(q2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))))))/(Power(fPi,2)*Power(Power(mEta,2) - Pair(Momentum(q,vects)-Momentum(q2,vects),Momentum(q,vects)-Momentum(q2,vects)),2));
  }
  
}
}
