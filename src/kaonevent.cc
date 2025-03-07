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
  lepton.set_mass(PDG::mass(lepton.pdg));

	// assume first a proton in produced
	N1.pdg = 2212;
	N1.set_mass(PDG::mass(N1.pdg));

	switch(N0.pdg) {
		case PDG::pdg_proton:		// proton -> proton
			kaon.pdg = 321;

			// set form factors
    	TwoThreeScatter::singlekaon::ACT = 2; 
    	TwoThreeScatter::singlekaon::BCT = -TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F) / 2;
    	TwoThreeScatter::singlekaon::ACRLambda = TwoThreeScatter::singlekaon::D + 3*TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::AKP = 2;
    	TwoThreeScatter::singlekaon::APi = -1;
    	TwoThreeScatter::singlekaon::AEta = 1;
			break;

		case PDG::pdg_neutron:	// neutron -> proton
			kaon.pdg = 311;

			// set form factors
    	TwoThreeScatter::singlekaon::ACT = 1; 
    	TwoThreeScatter::singlekaon::BCT = -TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F);
    	TwoThreeScatter::singlekaon::ACRLambda = TwoThreeScatter::singlekaon::D + 3*TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::AKP = 1;
    	TwoThreeScatter::singlekaon::APi = -2;
    	TwoThreeScatter::singlekaon::AEta = 0;
			break;
	}
	
	kaon.set_mass(PDG::mass(kaon.pdg));

  particle N0_Eb = N0; // nucleon with 4 momentum adjusted for binding energy
  //N0_Eb.t -= _E_bind; will implement later
 
  // scatter
 	constexpr int num = 3; // number of produced particles
 	particle particles[num] = { lepton, kaon, N1 };
 	int q2 = scatter_n(num, nu, N0_Eb, particles);
 	if(q2==0){
		//std::cout << "scatter error\n";
	return 0; //indicates interaction is forbidden by kinematics
	}
	//std::cout << "scatter worked\n";

	lepton = particles[0];
	kaon = particles[1];
	N1 = particles[2];

	vect va{ N0_Eb };
	vect vb{ nu };
	vect v1{ N1 };
	vect v2{ kaon };
	vect v3{ lepton };

	//std::cout << "va: " << va << '\n';
	//std::cout << "vb: " << vb << '\n';
	//std::cout << "v1: " << v1 << '\n';
	//std::cout << "v2: " << v2 << '\n';
	//std::cout << "v3: " << v3 << '\n';

	vec vcms = (va + vb).v();

	va.boost(-vcms);
	vb.boost(-vcms);
	v1.boost(-vcms);
	v2.boost(-vcms);
	v3.boost(-vcms);

	vect v12{ v1 + v2 };

	// generate vector with direction of W in CMS
	vec cms_dir = vec(v3) / vec(v3).length();

  // cross section calculation

  // kinematic variables for cross-section calculation
  double Nuc0_mass{ N0.mass() };
  double Nuc1_mass{ N1.mass() };
  double Kaon_mass{ kaon.mass() };
  double Lepton_mass{ lepton.mass() }; 

  double s{ e.s() };
  double W{ std::sqrt(v12*v12) };
	//std::cout << "W: " << W << '\n';
	
  double theta{ std::acos(v2.v()*vb.v() / (v2.v().length()*vb.v().length())) };

  // differential cross section
  double xsec = single_kaon_diff_xsec(Nuc0_mass, Nuc1_mass, Kaon_mass, Lepton_mass, s, W, theta);


  if(N0_Eb.pdg == 2112) {    // neutron -> neutron
		
		double xsec_2{ 0 };

		// new particles
		particle lepton_2{ lepton };

		particle kaon_2{ kaon };
		kaon_2.pdg = 321;
		kaon_2.set_mass(PDG::mass(kaon_2.pdg));

		particle N1_2{ N1 };
		N1_2.pdg = 2112;
		N1_2.set_mass(PDG::mass(N1_2.pdg));

		particle hadron_blob;
		vect v12_2{ vect{ N1_2 } + vect{ kaon_2 } };
		hadron_blob.set_mass(std::sqrt( v12_2*v12_2 ));
		hadron_blob.p4() = v12_2;
		//std::cout << "v12: " << vect{ hadron_blob } << '\n';

		vect cms_Eb = vect{ N0_Eb } + vect{ nu };

		if(rescale_momenta(cms_Eb, cms_dir, lepton_2, hadron_blob)) {

			v3 = vect{ lepton_2 };
			v12 = vect{ hadron_blob };
			//std::cout << "Rescaled v12: " << v12 << '\n';

			v3.boost(-vcms);
			v12.boost(-vcms);

			W = std::sqrt(v12 * v12);
			//std::cout << "W: " << W << '\n';

			theta = std::acos(v2.v()*vb.v() / (v2.v().length()*vb.v().length()));

			// set form factors
    	TwoThreeScatter::singlekaon::ACT = 1; 
    	TwoThreeScatter::singlekaon::BCT = TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F);
    	TwoThreeScatter::singlekaon::ACRLambda = 0;
    	TwoThreeScatter::singlekaon::AKP = 1;
    	TwoThreeScatter::singlekaon::APi = 1;
    	TwoThreeScatter::singlekaon::AEta = 1;

			Nuc1_mass = N1_2.mass();
			Kaon_mass = kaon_2.mass();

			xsec_2 = single_kaon_diff_xsec(Nuc0_mass, Nuc1_mass, Kaon_mass, Lepton_mass, s, W, theta);

			if(frandom() < xsec_2 / (xsec + xsec_2)) {
				lepton = lepton_2;
				hadron_blob.decay(N1_2, kaon_2);
				N1 = N1_2;
				kaon = kaon_2;
			}

			xsec += xsec_2;
		}	
	}

  e.temp.push_back(lepton);
  e.temp.push_back(kaon);
  e.temp.push_back(N1);
  e.out.push_back(lepton);
  e.out.push_back(kaon);
  e.out.push_back(N1);

  e.weight = xsec / cm2;

  return e.weight * cm2;

//	return 0;
}


double kaonevent_2(params& p, event& e, nucleus& t)
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
  lepton.set_mass(PDG::mass(lepton.pdg));

	// assume first a proton in produced
	N1.pdg = 2212;
	N1.set_mass(PDG::mass(N1.pdg));

	switch(N0.pdg) {
		case PDG::pdg_proton:		// proton -> proton
			kaon.pdg = 321;

			// set form factors
    	TwoThreeScatter::singlekaon::ACT = 2; 
    	TwoThreeScatter::singlekaon::BCT = -TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F) / 2;
    	TwoThreeScatter::singlekaon::ACRLambda = TwoThreeScatter::singlekaon::D + 3*TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::AKP = 2;
    	TwoThreeScatter::singlekaon::APi = -1;
    	TwoThreeScatter::singlekaon::AEta = 1;
			break;

		case PDG::pdg_neutron:	// neutron -> proton
			kaon.pdg = 311;

			// set form factors
    	TwoThreeScatter::singlekaon::ACT = 1; 
    	TwoThreeScatter::singlekaon::BCT = -TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F);
    	TwoThreeScatter::singlekaon::ACRLambda = TwoThreeScatter::singlekaon::D + 3*TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::AKP = 1;
    	TwoThreeScatter::singlekaon::APi = -2;
    	TwoThreeScatter::singlekaon::AEta = 0;
			break;
	}
	
	kaon.set_mass(PDG::mass(kaon.pdg));

  particle N0_Eb = N0; // nucleon with 4 momentum adjusted for binding energy
  //N0_Eb.t -= _E_bind; will implement later
 
  // scatter
 	constexpr int num = 3; // number of produced particles
 	particle particles[num] = { lepton, kaon, N1 };
 	int scattered = scatter_n(num, nu, N0_Eb, particles);
 	if(scattered == 0){	return 0; }

	lepton = particles[0];
	kaon = particles[1];
	N1 = particles[2];

	vect p1{ N0_Eb };
	vect p2{ nu };
	vect q1{ N1 };
	vect q2{ kaon };
	vect q3{ lepton };

	vec vcms = (p1 + p2).v();

	p1.boost(-vcms);
	p2.boost(-vcms);
	q1.boost(-vcms);
	q2.boost(-vcms);
	q3.boost(-vcms);

  vect q{ p2 - q3 };

	// generate vector with direction of W in CMS
	vec cms_dir = vec{ q3 } / vec{ q3 }.length();

  // cross section calculation

  // kinematic variables for cross-section calculation
  double Nuc0_mass{ N0.mass() };
  double Nuc1_mass{ N1.mass() };
  double Kaon_mass{ kaon.mass() };
  double Lepton_mass{ lepton.mass() }; 

  double s{ e.s() };
	
  // differential cross section
  vect cms_vects[6] = { p1, p2, q1, q2, q3, q };
  double xsec = single_kaon_diff_xsec_2(s, Nuc1_mass, Kaon_mass, Lepton_mass, cms_vects);


  if(N0_Eb.pdg == 2112) {    // neutron -> neutron
		
		double xsec_2{ 0 };

		// new particles
		particle lepton_2{ lepton };

		particle kaon_2{ kaon };
		kaon_2.pdg = 321;
		kaon_2.set_mass(PDG::mass(kaon_2.pdg));

		particle N1_2{ N1 };
		N1_2.pdg = 2112;
		N1_2.set_mass(PDG::mass(N1_2.pdg));

		particle hadron_blob;
		vect q12{ vect{ N1_2 } + vect{ kaon_2 } };
		hadron_blob.set_mass(std::sqrt( q12*q12 ));
		hadron_blob.p4() = q12;

		vect cms_Eb = vect{ N0_Eb } + vect{ nu };

		if(rescale_momenta(cms_Eb, cms_dir, lepton_2, hadron_blob)) {

			hadron_blob.decay(N1_2, kaon_2);

      q1 = vect{ N1_2 };
      q2 = vect{ kaon_2 };
			q3 = vect{ lepton_2 };

      q1.boost(-vcms);
      q2.boost(-vcms);
			q3.boost(-vcms);

      q = p2 - q3;

			// set form factors
    	TwoThreeScatter::singlekaon::ACT = 1; 
    	TwoThreeScatter::singlekaon::BCT = TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F;
    	TwoThreeScatter::singlekaon::ACRSigma = -(TwoThreeScatter::singlekaon::D - TwoThreeScatter::singlekaon::F);
    	TwoThreeScatter::singlekaon::ACRLambda = 0;
    	TwoThreeScatter::singlekaon::AKP = 1;
    	TwoThreeScatter::singlekaon::APi = 1;
    	TwoThreeScatter::singlekaon::AEta = 1;

			Nuc1_mass = N1_2.mass();
			Kaon_mass = kaon_2.mass();
      
      cms_vects[2] = q1;
      cms_vects[3] = q2;
      cms_vects[4] = q3;
      cms_vects[5] = q;

			xsec_2 = single_kaon_diff_xsec_2(s, Nuc1_mass, Kaon_mass, Lepton_mass, cms_vects);

			if(frandom() < xsec_2 / (xsec + xsec_2)) {
				lepton = lepton_2;
				N1 = N1_2;
				kaon = kaon_2;
			}

			xsec += xsec_2;
		}	
	}

  e.temp.push_back(lepton);
  e.temp.push_back(kaon);
  e.temp.push_back(N1);
  e.out.push_back(lepton);
  e.out.push_back(kaon);
  e.out.push_back(N1);

  e.weight = xsec / cm2;

  return e.weight * cm2;

//	return 0;
}

