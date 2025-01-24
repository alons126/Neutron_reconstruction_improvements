//
// Created by Alon Sportes on 21/01/2025.
//

#ifndef VETOCUTS_H
#define VETOCUTS_H

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
double M_miss_lcut = 0.85;
double M_miss_ucut = 1.05;
// double M_miss_lcut = 0.7;
// double M_miss_ucut = 1.2;

double Beta_n_lcut = 0.15;
double Beta_n_ucut = 0.80;
double Theta_n_lcut = Theta_miss_lcut;
double Theta_n_ucut = Theta_miss_ucut;
int Status_n_cut = 0;

#endif //VETOCUTS_H
