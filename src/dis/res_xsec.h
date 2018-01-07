#ifndef _RES_XSEC_H
#define _RES_XSEC_H

#include "res_kinematics.h"

//! store xsec / spp parameters and xsec contributions from different channels
struct res_xsec {
  int j, k, l;       //!< indices for SPP table (see singlepion.cc)
  int final_charge;  //!< total electric charge of the pion-nucleon system

  double dis_pip, dis_pi0, dis_pim;        //!< contribution to RES from DIS
  double delta_pip, delta_pi0, delta_pim;  //!< contribution to RES from Delta

  double dis_total, delta_total;  //!< total contribution from RES and DIS

  double from_dis;  //!< the strength of DIS background

  bool is_cc;  //!< true for CC events

  res_xsec(res_kinematics &kin, const bool cc);  //!< initilize common parameters

  //! return (normalized) total cross section (Delta + DIS or just Pythia)
  inline double get_total(const double &jacobian) {
    return 1e-38 * jacobian * (delta_total ? (dis_total + delta_total) : from_dis);
  }

  //! Pythia part is skipped if True
  inline bool is_no_dis() { return from_dis == 0; }

  //! get SPP contribution for the background for given pion
  double get_dis_spp(const int pion_code, const res_kinematics &kin, const params &p);

  //! get Delta contribution for given pion
  double get_delta_spp(const int pion_pdg, const int nucleon_pdg, const res_kinematics &kin, const params &p);

  //! calculate cross sections for the no-Pythia scenario
  void set_xsec_nopythia(const res_kinematics &kin, const params &p);

  //! calcaulte cross sections for SPP with Pythia contribution
  void set_xsec(res_kinematics &kin, const params &p, const int pion_pdg, const int nucleon_pdg,
                const double neutrino_energy);

  //! return random (xsec based) final pion pdg
  int get_pion_pdg();
};

//! cross section reduction to remove contribution from pion-less delta decay
double pdd_red(double energy);

#endif