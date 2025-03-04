#include "generatormt.h"
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
#include "matrix_elements.h"

#define LOCALKF localkf_O


double kaonevent(params& p, event& e, nucleus& t)
{
  /*

  add some description here

  */

  // flags
  e.flag.kaon = true;// need to add kaon flag to event1
  e.weight = 0;

	// disable for anti neutrinos for now
	if(e.in[0].pdg < 0) return 0;

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
  lepton.set_mass(PDG::mass(lepton.pdg));

  particle N0_Eb = N0; // nucleon with 4 momentum adjusted for binding energy
  //N0_Eb.t -= _E_bind; will implement later
 
  

  // cross section calculation

  // kinematic variables for cross-section calculation
  double Nuc0_mass{ N0.mass() };
  double Nuc1_mass{ };
  double Kaon_mass{ };
  double Lepton_mass{ lepton.mass() }; 

  double s{ e.s() };
  double W{ };
  double theta{ };

  double xsec_proton{ };    // cross section for final state proton
  double xsec_neutron{ };   // cross section for final state neutron

  if(N0_Eb.pdg == 2212) {   // intital proton -> final state proton
    // final proton
    // set final kaon
    kaon.pdg = 321;
    kaon.set_mass(PDG::mass(kaon.pdg));
    Kaon_mass = kaon.mass();

    // set final nucleon
    N1.pdg = 2212;          // proton
    N1.set_mass(PDG::mass(N1.pdg));
    Nuc1_mass = N1.mass();

    // set form factors
    TwoThreeScatter::singlekaon::ACT = 2; 
    TwoThreeScatter::singlekaon::BCT = -TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F) / 2;
    TwoThreeScatter::singlekaon::ACRLambda = TwoThreeScatter::singlekaon::D + 3*TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::AKP = 2;
    TwoThreeScatter::singlekaon::APi = -1;
    TwoThreeScatter::singlekaon::AEta = 1;

    // set kinematic variables
    W = e.W();
    theta = std::acos(e.costheta());

    // differential cross section
    xsec_proton = single_kaon_diff_xsec(Nuc0_mass, Nuc1_mass, Kaon_mass, Lepton_mass, s, W, theta);
    
    // neutron
    xsec_neutron = 0;
  }

  else if(N0_Eb.pdg == 2112) {    // initial neutron -> final state neutron or proton
    // final proton
    // set final kaon
    kaon.pdg = 311;
    kaon.set_mass(PDG::mass(kaon.pdg));
    Kaon_mass = kaon.mass();

    // set final nucleon
    N1.pdg = 2212;          // proton
    N1.set_mass(PDG::mass(N1.pdg));
    Nuc1_mass = N1.mass();

    // set form factors
    TwoThreeScatter::singlekaon::ACT = 1; 
    TwoThreeScatter::singlekaon::BCT = -TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F);
    TwoThreeScatter::singlekaon::ACRLambda = TwoThreeScatter::singlekaon::D + 3*TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::AKP = 1;
    TwoThreeScatter::singlekaon::APi = -2;
    TwoThreeScatter::singlekaon::AEta = 0;
    
    // set kinematic variables
    W = e.W();
    theta = std::acos(e.costheta());

    // differential cross section
    xsec_proton = single_kaon_diff_xsec(Nuc0_mass, Nuc1_mass, Kaon_mass, Lepton_mass, s, W, theta);
    
    // neutron
    // set final kaon
    kaon.pdg = 321;
    kaon.set_mass(PDG::mass(kaon.pdg));
    Kaon_mass = kaon.mass();

    // set final nucleon
    N1.pdg = 2112;          // neutron
    N1.set_mass(PDG::mass(N1.pdg));
    Nuc1_mass = N1.mass();

    // set form factors
    TwoThreeScatter::singlekaon::ACT = 1; 
    TwoThreeScatter::singlekaon::BCT = TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F;
    TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F);
    TwoThreeScatter::singlekaon::ACRLambda = 0;
    TwoThreeScatter::singlekaon::AKP = 1;
    TwoThreeScatter::singlekaon::APi = 1;
    TwoThreeScatter::singlekaon::AEta = 1;

    // differential cross section
    xsec_neutron = single_kaon_diff_xsec(Nuc0_mass, Nuc1_mass, Kaon_mass, Lepton_mass, s, W, theta);
  }

  // differential cross section
  double xsec = (frandom() > xsec_proton / (xsec_proton + xsec_neutron)) ? xsec_neutron : xsec_proton;

  // scatter
  constexpr int num = 3; // number of produced particles
  particle particles[num] = { lepton, kaon, N1 };
  double q2 = scatter_n(num, nu, N0_Eb, particles);
  
  if(q2==0) return 0; //indicates interaction is forbidden by kinematics

  e.temp.push_back(lepton);
  e.temp.push_back(kaon);
  e.temp.push_back(N1);
  e.out.push_back(lepton);
  e.out.push_back(kaon);
  e.out.push_back(N1);

  e.weight = xsec / cm2;

  return e.weight * cm2;

}
