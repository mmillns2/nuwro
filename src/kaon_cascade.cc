#include "kaon_cascade.h"


double get_xsec(TGraph* g, double p) {
	  if(!g) return 0;
	  double p_min = g->GetX()[0];
	  double p_max = g->GetX()[g->GetN() - 1];
	  double slope = g->Eval(p_min) / p_min;
	  double const_max = g->Eval(p_max);
	  if(p < p_min) return p*slope;
	  if(p > p_max) return const_max;
	  return g->Eval(p);
}

// E = CMS energy
// Plab = hyperon momentum in nucleon rest frame
// xsecs = { PK+ -> PK+, NK+ -> NK+, PK- -> PK-, NK- -> NK- }
// kaon_state = { Kplus, Kminus, Kzero, KzeroBar } (initial kaon state)

void kaon_exp_xsec(double E, double Plab, double xsecs[], kaon_state state) {
  
  constexpr int n_PK{ 47 };
  double pVals_PK[n_PK]{ 145, 175, 205, 235, 265, 295, 325, 355, 385, 432, 479, 500, 525, 565, 603, 613, 646, 689, 726, 731, 772, 813, 857, 865, 899, 910, 939, 970, 1098, 1170, 1207, 1255, 1310, 1370, 1400, 1450, 1495, 1540, 1600, 1640, 1705, 1740, 1790, 1880, 1965, 2050, 2125 };
  double xsecVals_PK[n_PK]{ 11.8, 12.1, 11.6, 11.8, 11.6, 12.3, 11.1, 11.4, 11.4, 11.345, 11.92, 12.6, 11.34, 12.085, 11.86, 12.1, 11.675, 11.88, 12.2, 11.38, 11.6, 11.295, 11.255, 11.98, 11.315, 11.75, 11.38, 11.62, 10.75, 10.55, 10.69, 10.32, 10.4, 9.83, 9.51, 9.17, 9.35, 8.79, 8.55, 8.35, 8.46, 8.04, 7.56, 7.36, 6.72, 6.19, 5.96 };
  TGraph* g_PK = new TGraph(n_PK, pVals_PK, xsecVals_PK);

  constexpr int n_NK{ 13 };
  double pVals_NK[n_NK]{ 640, 720, 780, 850, 900, 980, 1060, 1130, 1210, 1290, 1350, 1420, 1510 };
  double xsecVals_NK[n_NK]{ 5.62973, 5.94389, 6.53765, 6.56593, 5.71142, 4.8035, 4.75009, 4.24115, 3.51544, 3.18243, 3.11018, 2.61066, 2.11429 };
  TGraph* g_NK = new TGraph(n_NK, pVals_NK, xsecVals_NK);

  //std::array<double> pVals_PK_anti{  };
  //std::array<double> xsecVals_PK_anti{  };

  //std::array<double> pVals_NK_anti{  };
  //std::array<double> xsecVals_NK_anti{  };

  
  double xsec_PK{ get_xsec(g_PK, Plab) };
  double xsec_NK{ get_xsec(g_NK, Plab) };
  double xsec_PK_anti{  };
  double xsec_NK_anti{  };

  if(state == Kplus) {           // initial K+
    xsecs[0] = xsec_PK;          // P K+ -> P K+
    xsecs[1] = xsec_NK;          // N K+ -> N K+
  }

  else if(state == Kminus) {
    xsecs[0] = xsec_PK_anti;     // P K- -> P K-
    xsecs[1] = xsec_NK_anti;     // N K- -> N K-
  }
  
	delete g_PK;
	delete g_NK;
}

void get_kaon_state(kaon_state state, double xsecs[], int& ij, particle p[]) {

  double rand_num{ frandom() };
  double xsec_total{ xsecs[0] + xsecs[1] };
  double xsec_P{ xsecs[0] };  // P K -> P K
  double xsec_N{ xsecs[1] };  // N K -> N K

  if(state == Kplus) {  // X K+ -> X K+
    
    if(rand_num < xsec_P/xsec_total) {  // P K+ -> P K+ 
      ij = 0;
      p[0].set_pdg_and_mass(PDG::pdg_proton, PDG::mass_proton);
		} 

    else {  // N K+ -> N K+
      ij = 1;
      p[0].set_pdg_and_mass(PDG::pdg_neutron, PDG::mass_neutron);
    }

    p[1].set_pdg_and_mass(PDG::pdg_KP, PDG::mass_KP);
  }

  else if(state == Kminus) {  // X K- -> X K-

    if(rand_num < xsec_P/xsec_total) {  // P K- -> P K- 
      ij = 2;
      p[0].set_pdg_and_mass(PDG::pdg_proton, PDG::mass_proton);
		} 

    else {  // N K- -> N K-
      ij = 1;
      p[0].set_pdg_and_mass(PDG::pdg_neutron, PDG::mass_neutron);
    }

    p[1].set_pdg_and_mass(-PDG::pdg_KP, PDG::mass_KP);
  }
  
  else { }
}





