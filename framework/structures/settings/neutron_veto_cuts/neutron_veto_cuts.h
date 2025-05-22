//
// Created by Alon Sportes on 21/01/2025.
//

#ifndef NEUTRON_VETO_CUTS_H
#define NEUTRON_VETO_CUTS_H

#include <iostream>

namespace neutron_veto_cuts {

// Cut settings -------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool Apply_CTOF_veto = false;

// NOTE: Step 0 cuts - disabled following RG-M meeting (29-01-25) and Adi meeting (02-02-25)
bool Apply_Step0_Cuts = true;

// Cut values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

int Num_of_e_cut = 1;
int Num_of_p_cut = 1;

double dVz_pCD_cut = 4;
double P_pCD_lcut = 0.3;
double P_pCD_ucut = 1.5;
double pCD_chi2_lcut = -2.0;
double pCD_chi2_ucut = 2.5;

double dVz_pFD_cut = 5;
double P_pFD_lcut = 0.4;
double P_pFD_ucut = 3.0;
double pFD_chi2_lcut = -2.0;
double pFD_chi2_ucut = 2.5;

double P_miss_lcut = 0.2;
double P_miss_ucut = 1.5;
double Theta_miss_lcut = 40.;
double Theta_miss_ucut = 135.;
// double M_miss_lcut = 0.0;
// double M_miss_ucut = 9999.;
double M_miss_lcut = 0.85;
double M_miss_ucut = 1.05;
// double M_miss_lcut = 0.7;
// double M_miss_ucut = 1.2;

double Beta_n_lcut = 0.15;
double Beta_n_ucut = 0.80;
double Theta_n_lcut = Theta_miss_lcut;
double Theta_n_ucut = Theta_miss_ucut;
int Status_n_cut = 0;

/* Good and bad neutron definitions */
double GN_theta_n_miss_ucut = 20.;
// double GN_theta_n_miss_lcut = 0.;
double GN_dpp_ucut = 0.4;
double GN_dpp_lcut = -0.3;
// double BN_theta_n_miss_ucut = 9999.;
double BN_theta_n_miss_lcut = 40.;
double BN_dpp_ucut = -1;
// double BN_dpp_lcut = -9999;

/* Step 0 cuts */
double dBeta_n_cut = 0.01;
double Vz_n_lcut = -40.;
double Vz_n_ucut = 45.;
double ToF_n_lcut = 0.;
double ToF_n_ucut = 20.;

/* Step 1 cuts */
double Edep_CND_lcut = 5.;

/* Step 2 cuts (PARTIAL!) */
int Cluster_size_cut = 1;
int CND1_LayerMult_cut = 1;
int CND2andCND3_LayerMult_ucut = 2;

};  // namespace neutron_veto_functions

#endif  // NEUTRON_VETO_CUTS_H
