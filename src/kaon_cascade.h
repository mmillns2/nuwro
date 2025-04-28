#pragma once


#include "particle.h"
#include "TGraph.h"


enum kaon_state { Kplus, Kminus, Kzero, Kzerobar };

double get_xsec(TGraph* g, double p);
double get_xsec_0(TGraph* g, double p);

void kaon_exp_xsec(double E, double Plab, double sigma[], kaon_state state);

void get_kaon_state(kaon_state state, double xsecs[], int& ij, particle p[]);
