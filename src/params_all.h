#define PARAMS_ALL()\
PARAM(int,random_seed,0)\
PARAM(int,number_of_events,100000)\
PARAM(int,number_of_test_events,1000000)\
PARAM(int,save_test_events,0)\
PARAM(int,user_events,0)\
PARAM(line,user_params,"")\
PARAM(int,beam_type,0)\
PARAM(line,beam_energy,"1000")\
PARAM(int,beam_particle,14)\
PARAM(vec,beam_direction,"0 0 1")\
PARAM(line,beam_content,"")\
PARAM(string,beam_folder,"flux")\
PARAM(int,beam_file_first,1)\
PARAM(int,beam_file_limit,0)\
PARAM(double,beam_pot_per_file,-1)\
PARAM(string,geom_length_units,"cm")\
PARAM(double,geom_density_convert,1.0)\
PARAM(string,beam_length_units,"cm")\
PARAM(bool,beam_weighted,0)\
PARAM(vec,beam_offset,"0 0 0")\
PARAM(int,beam_placement,0)\
PARAM(string,beam_inputroot,"")\
PARAM(string,beam_inputroot_flux,"")\
PARAM(string,beam_inputroot_nue,"")\
PARAM(string,beam_inputroot_nueb,"")\
PARAM(string,beam_inputroot_numu,"")\
PARAM(string,beam_inputroot_numub,"")\
PARAM(string,beam_inputroot_nutau,"")\
PARAM(string,beam_inputroot_nutaub,"")\
PARAM(line,beam_atmo_files,"")\
PARAM(int,beam_test_only,0)\
PARAM(int,target_type,0)\
PARAM(int,nucleus_p,6)\
PARAM(int,nucleus_n,6)\
PARAM(double,nucleus_E_b,34)\
PARAM(double,nucleus_kf,220)\
PARAM(line,target_content,"")\
PARAM(string,geo_file,"target/ND280_v9r7p5.root")\
PARAM(string,geo_name,"ND280Geometry_v9r7p5")\
PARAM(string,geo_volume,"")\
PARAM(vec,geo_o,"0 0 0")\
PARAM(vec,geo_d,"2000 2000 5000")\
PARAM(int,nucleus_target,2)\
PARAM(int,nucleus_model,1)\
PARAM(bool,dyn_qel_cc,1)\
PARAM(bool,dyn_qel_nc,0)\
PARAM(bool,dyn_res_cc,1)\
PARAM(bool,dyn_res_nc,0)\
PARAM(bool,dyn_dis_cc,1)\
PARAM(bool,dyn_dis_nc,0)\
PARAM(bool,dyn_coh_cc,1)\
PARAM(bool,dyn_coh_nc,0)\
PARAM(bool,dyn_mec_cc,1)\
PARAM(bool,dyn_mec_nc,0)\
PARAM(bool,dyn_hyp_cc,1)\
PARAM(bool,dyn_kaon,1)\
PARAM(bool,dyn_lep,1)\
PARAM(bool,dyn_qel_el,0)\
PARAM(bool,dyn_res_el,0)\
PARAM(double,el_costh_lab,90)\
PARAM(double,el_costh_del,2)\
PARAM(string,eel_alg,"old")\
PARAM(int,qel_kinematics,0)\
PARAM(int,qel_vector_ff_set,2)\
PARAM(int,qel_axial_ff_set,1)\
PARAM(int,qel_rpa,1)\
PARAM(int,qel_strange,1)\
PARAM(int,qel_strangeEM,0)\
PARAM(double,delta_s,0)\
PARAM(double,qel_cc_vector_mass,840)\
PARAM(double,qel_cc_axial_mass,1030)\
PARAM(double,qel_minerva_ff_scale,0)\
PARAM(double,qel_deuterium_ff_scale,0)\
PARAM(double,qel_nc_axial_mass,1030)\
PARAM(double,qel_s_axial_mass,1030)\
PARAM(bool,flux_correction,1)\
PARAM(int,sf_method,0)\
PARAM(bool,sf_fsi,1)\
PARAM(bool,sf_coulomb,1)\
PARAM(int,sf_pb,1)\
PARAM(bool,cc_smoothing,1)\
PARAM(int,delta_FF_set,1)\
PARAM(int,e_spp_ff_set,4)\
PARAM(int,delta_selfenergy,0)\
PARAM(double,pion_axial_mass,0.94)\
PARAM(double,pion_C5A,1.19)\
PARAM(int,delta_angular,2)\
PARAM(int,spp_precision,500)\
PARAM(double,res_dis_cut,1600)\
PARAM(double,bkgrscaling,0.0)\
PARAM(bool,coh_mass_correction,1)\
PARAM(bool,coh_new,1)\
PARAM(int,coh_kind,2)\
PARAM(int,mec_kind,3)\
PARAM(double,mec_ratio_pp,0.85)\
PARAM(double,mec_ratio_ppp,0.8)\
PARAM(double,mec_central_motion,0.0)\
PARAM(double,mec_back_to_back_smearing,0.0)\
PARAM(int,mec_pb_trials,30)\
PARAM(bool,MEC_pauli_blocking,1)\
PARAM(double,MEC_cm_direction,0.0)\
PARAM(int,mec_scaling,0)\
PARAM(bool,hyp_lambda,1)\
PARAM(bool,hyp_sigma_zero,1)\
PARAM(bool,hyp_sigma_minus,1)\
PARAM(double,hyp_g2_Re_part,0)\
PARAM(double,hyp_g2_Im_part,0)\
PARAM(bool,hyp_su3_sym_breaking,0)\
PARAM(double,hyp_axial_mass,1030)\
PARAM(bool,hyp_effmass,true)\
PARAM(double,hyp_Lambda_Eb,27)\
PARAM(double,hyp_Sigma_Eb,-70)\
PARAM(bool, kaon_Kplus, 1)\
PARAM(bool, kaon_Kminus, 1)\
PARAM(bool, kaon_Kzero, 1)\
PARAM(bool, kaon_KzeroBar, 1)\
PARAM(bool, kaon_CT, 1)\
PARAM(bool, kaon_CRSigma, 1)\
PARAM(bool, kaon_CRLambda, 1)\
PARAM(bool, kaon_KaonPole, 1)\
PARAM(bool, kaon_PionInFlight, 1)\
PARAM(bool, kaon_EtaInFlight, 1)\
PARAM(bool, kaon_CT_anti, 1)\
PARAM(bool, kaon_Sigma_anti, 1)\
PARAM(bool, kaon_SigmaStar_anti, 1)\
PARAM(bool, kaon_Lambda_anti, 1)\
PARAM(bool, kaon_KaonPole_anti, 1)\
PARAM(bool, kaon_PionInFlight_anti, 1)\
PARAM(bool, kaon_EtaInFlight_anti, 1)\
PARAM(bool, kaon_dipoleff, 1)\
PARAM(double, kaon_MF, 1000)\
PARAM(bool, kaon_fullBeam, 1)\
PARAM(int, kaon_beam_bins, 10)\
PARAM(int, kaon_W_bins, 10)\
PARAM(int, kaon_theta_bins, 10)\
PARAM(int, kaon_thetaStar_bins, 10)\
PARAM(int, kaon_phiStar_bins, 10)\
PARAM(bool,kaskada_on,1)\
PARAM(double,kaskada_w,7)\
PARAM(bool,kaskada_redo,0)\
PARAM(bool,kaskada_writeall,0)\
PARAM(string,formation_zone,"fz-new")\
PARAM(double,tau,8)\
PARAM(double,formation_length,1)\
PARAM(bool,first_step,1)\
PARAM(double,step,0.2)\
PARAM(double,kaskada_NN_mfp_scale,1)\
PARAM(int,kaskada_NN_xsec,2)\
PARAM(int,kaskada_NN_inel,2)\
PARAM(int,kaskada_NN_angle,3)\
PARAM(int,kaskada_NN_corr,1)\
PARAM(int,kaskada_piN_xsec,1)\
PARAM(bool,pauli_blocking,1)\
PARAM(bool,mixed_order,1)\
PARAM(bool,kaskada_events,0)\
PARAM(string,kaskada_events_file,"")\

