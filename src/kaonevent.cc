#include "jednostki.h"
#include <cassert>          
#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "params.h"
#include "event1.h"
#include "hypevent.h"
#include "kinematics.h"
#include "pdg.h"
#include "nucleus.h"
#include <cstdlib>
#include "hiperon_sigma.h"
#include "hyperon_interaction.h"
#include "scatter.h"
#include "dis/LeptonMass.h"

#include "kaon_diff_xsec.h"
#include "23Scatter.h"
#include "single_kaon_matrix_elements/matrix_elements.h"

#define LOCALKF localkf_O


double kaonevent(params& p, event& e, nucleus& t)
{
  /*

  add some description here

  */

  // flags
  e.flag.kaon = true;// need to add kaon flag to event1
  e.weight = 0;

  // particles
  particle nu = e.in[0];  // initial antineutrino
  particle N0 = e.in[1];  // initial nucleon
  particle lepton;        // outgoing antilepton
  particle kaon;          // created kaon
  particle N1;            // outgoing nucleon

  // same coordinates for initial/final particles
  lepton.r = N0.r;
  kaon.r = N0.r;
  N1.r = N0.r;

  // adjust the lepton flavour
  lepton.pdg = nu.pdg - (nu.pdg > 0 ? 1 : -1);

  // set final particle masses
  kaon.set_mass(PDG::mass(kaon.pdg));
  lepton.set_mass(PDG::mass(lepton.pdg));
  N1.set_mass(PDG::mass(N1.pdg));

  particle N0_Eb = N0; // nucleon with 4 momentum adjusted for binding energy
  //N0_Eb.t -= _E_bind; will implement later
 
  // scatter
  constexpr int num = 3; // number of produced particles
  particle particles[num] = { lepton, kaon, N1 };
  double q2 = scatter_n(num, nu, N0_Eb, particles);
  
  if(q2==0) return 0; //indicates interaction is forbidden by kinematics

  vect nu4 = nu;

  // not sure whats going on here
  if(p.kaon_effmass) nu4.boost (-N0_Eb.v());  // go to target frame
  else nu4.boost (-N0.v());  

  double Enu0=nu4.t;    // neutrino energy in target frame   

  // Four momenta of particles involved
  vect v1(nu);
  vect v2(N0_Eb);
  vect v3(lepton);
  vect v4(kaon);
  vect v5(N1);

  //boost these to CMS for calculations
  vec vcms = (vect(nu) + vect(N0_Eb)).v();

  v1.boost(-vcms);
  v2.boost(-vcms);
  v3.boost(-vcms);
  v4.boost(-vcms);
  v5.boost(-vcms);

  // Generate vector with direction of Kaon in CMS
  vec cms_dir = vec(v3)/vec(v3).length();

  /////////////////////////////////////////////////////////
  // Generate Cross Sections
  /////////////////////////////////////////////////////////

  double kin = v1.length(); //incoming neutrino momentum
  double kout = v3.length(); //outgoing lepton momentum  

  // kinematic variables for cross-section calculation
  double Nuc0_mass;
  double Nuc1_mass;
  double Kaon_mass;
  double Lepton_mass{ lepton.mass() };

  double s{ e.s() };
  double W{ e.W() };
  double theta{ std::acos(e.costheta()) };

  // Switch to use effective masses of particles in cross section calculation
  if(p.kaon_effmass){    // change to kaon effective mass
     Nuc0_mass = sqrt(N0_Eb*N0_Eb);
     Nuc1_mass = sqrt(N1*N1);
     Kaon_mass = sqrt(kaon*kaon);
  }
  else {
     Nuc0_mass = N0.mass();
     Nuc1_mass = N1.mass();
     Kaon_mass = kaon.mass();
  }

  // set form factors
  if(N0_Eb.pdg == 2212 && N1.pdg == 2212) {  // proton-proton
    TwoThreeScatter::singlekaon::ACT = 2; 
    TwoThreeScatter::singlekaon::BCT = -TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F) / 2;
    TwoThreeScatter::singlekaon::ACRLambda = TwoThreeScatter::singlekaon::D + 3*TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::AKP = 2;
    TwoThreeScatter::singlekaon::APi = -1;
    TwoThreeScatter::singlekaon::AEta = 1;
  }
  else if(N0_Eb.pdg == 2112 && N1.pdg == 2112) {  // neutron-neutron
    TwoThreeScatter::singlekaon::ACT = 1; 
    TwoThreeScatter::singlekaon::BCT = TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F);
    TwoThreeScatter::singlekaon::ACRLambda = 0;
    TwoThreeScatter::singlekaon::AKP = 1;
    TwoThreeScatter::singlekaon::APi = 1;
    TwoThreeScatter::singlekaon::AEta = 1;
  }
  else {  // neutron-proton
    TwoThreeScatter::singlekaon::ACT = 1; 
    TwoThreeScatter::singlekaon::BCT = -TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F);
    TwoThreeScatter::singlekaon::ACRLambda = TwoThreeScatter::singlekaon::D + 3*TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::AKP = 1;
    TwoThreeScatter::singlekaon::APi = -2;
    TwoThreeScatter::singlekaon::AEta = 0;
}

  // differential cross section
  double xsec = single_kaon_diff_xsec(Nuc0_mass, Nuc1_mass, Kaon_mass, Lepton_mass, s, W, theta); 

  
  e.temp.push_back(lepton);
  e.temp.push_back(kaon);
  e.temp.push_back(N1);
  e.out.push_back(lepton);
  e.out.push_back(kaon);
  e.out.push_back(N1);

  e.weight = xsec / cm2;

  return e.weight * cm2;

}
