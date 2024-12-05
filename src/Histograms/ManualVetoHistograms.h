#ifndef MANUALVETOHISTOGRAMS_H
#define MANUALVETOHISTOGRAMS_H

#include <cstdlib>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

// ======================================================================================================================================================================
// Andrew's histograms
// ======================================================================================================================================================================

#pragma region /* Andrew's histograms - start */

/////////////////////////////////////
// Prepare histograms
/////////////////////////////////////

vector<TH1 *> hist_list_1_A;
vector<TH2 *> hist_list_2_A;

gStyle->SetTitleXSize(0.05);
gStyle->SetTitleYSize(0.05);

gStyle->SetTitleXOffset(0.8);
gStyle->SetTitleYOffset(0.8);

char temp_name_A[100];
char temp_title_A[100];

// (e,e'p) plots
// ======================================================================================================================================================================

/* Proton histograms (from Erin) */
TH1D *h_p_multiplicity_epCD = new TH1D("p_multiplicity_epCD", "Number of CD Protons in Event", 10, 0, 10);
hist_list_1_A.push_back(h_p_multiplicity_epCD);
TH2D *h_p_angles_epCD = new TH2D("p_angles_epCD", "CD Proton Angular Distribution;#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
hist_list_2.push_back(h_p_angles_epCD);

TH2D *h_dbeta_p_epCD = new TH2D("dbeta_p_epCD", "#Delta #beta vs CD proton momentum;P_{p} [GeV/c];#Delta#beta", 50, 0, 3, 50, -0.2, 0.2);
hist_list_2.push_back(h_dbeta_p_epCD);
TH1D *h_dVz_p_epCD = new TH1D("dVz_p_epCD", "Vertex correlation CD between proton and electron;dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
hist_list_1_A.push_back(h_dVz_p_epCD);
TH1D *h_Chi2pid_p_epCD = new TH1D("Chi2pid_p_epCD", "CD Proton #chi^{2}_{p};#chi^{2}_{p};Counts", 50, -6, 6);
hist_list_1_A.push_back(h_Chi2pid_p_epCD);

TH2D *h_dbeta_p_epFD = new TH2D("dbeta_p_epFD", "#Delta #beta vs FD proton momentum;P_{p} [GeV/c];#Delta#beta", 50, 0, 3, 50, -0.2, 0.2);
hist_list_2.push_back(h_dbeta_p_epFD);
TH1D *h_dVz_p_epFD = new TH1D("dVz_p_epFD", "Vertex correlation FD between proton and electron;dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
hist_list_1_A.push_back(h_dVz_p_epFD);
TH1D *h_Chi2pid_p_epFD = new TH1D("Chi2pid_p_epFD", "FD Proton #chi^{2}_{p};#chi^{2}_{p};Counts", 50, -6, 6);
hist_list_1_A.push_back(h_Chi2pid_p_epFD);

/* Neutron histograms (from Erin) */
TH1D *h_n_multiplicity_allN_epCD = new TH1D("n_multiplicity_allN_epCD", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_allN_epCD);
TH1D *h_n_multiplicity_goodN_epCD = new TH1D("n_multiplicity_goodN_epCD", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_goodN_epCD);
TH1D *h_n_multiplicity_badN_epCD = new TH1D("n_multiplicity_badN_epCD", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_badN_epCD);

TH1D *h_n_multiplicity_allN_epFD = new TH1D("n_multiplicity_allN_epFD", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_allN_epFD);
TH1D *h_n_multiplicity_goodN_epFD = new TH1D("n_multiplicity_goodN_epFD", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_goodN_epFD);
TH1D *h_n_multiplicity_badN_epFD = new TH1D("n_multiplicity_badN_epFD", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_badN_epFD);

/* Kinematical variables */
TH1D *h_theta_n_epCD = new TH1D("theta_n_epCD", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_epCD);
TH1D *h_phi_n_epCD = new TH1D("phi_n_epCD", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_epCD);
TH2D *h_theta_n_VS_phi_n_epCD = new TH2D("theta_n_VS_phi_n_epCD", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_epCD);
TH2D *h_theta_n_VS_beta_n_epCD = new TH2D("theta_VS_beta_epCD", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_epCD);

TH1D *h_P_n_epCD = new TH1D("P_n_epCD", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_epCD);
TH2D *h_P_n_VS_theta_n_epCD = new TH2D("P_n_VS_theta_n_epCD", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_epCD);

TH1D *h_P_miss_epCD = new TH1D("P_miss_epCD", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_epCD);
TH2D *h_P_miss_VS_theta_miss_epCD = new TH2D("P_miss_VS_theta_miss_epCD", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_epCD);

TH1D *h_P_n_minus_P_miss_epCD = new TH1D("P_n_minus_P_miss_epCD", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_epCD);
TH1D *h_P_n_x_minus_P_miss_x_epCD = new TH1D("P_n_x_minus_P_miss_x_epCD", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_epCD);
TH1D *h_P_n_y_minus_P_miss_y_epCD = new TH1D("P_n_y_minus_P_miss_y_epCD", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_epCD);
TH1D *h_P_n_z_minus_P_miss_z_epCD = new TH1D("P_n_z_minus_P_miss_z_epCD", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_epCD);

TH2D *h_P_n_VS_P_miss_epCD = new TH2D("P_n_VS_P_miss_epCD", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_epCD);
TH2D *h_P_n_x_VS_P_miss_x_epCD = new TH2D("P_n_x_VS_P_miss_x_epCD", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_epCD);
TH2D *h_P_n_y_VS_P_miss_y_epCD = new TH2D("P_n_y_VS_P_miss_y_epCD", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_epCD);
TH2D *h_P_n_z_VS_P_miss_z_epCD = new TH2D("P_n_z_VS_P_miss_z_epCD", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_epCD);

TH1D *h_E_p_CD_epCD = new TH1D("E_p_CD_epCD", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_epCD);
TH1D *h_E_p_FD_epCD = new TH1D("E_p_FD_epCD", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_epCD);
TH1D *h_E_miss_epCD = new TH1D("E_miss_epCD", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_epCD);
TH1D *h_M_miss_epCD = new TH1D("M_miss_epCD", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_epCD);
TH2D *h_M_miss_VS_P_n_epCD = new TH2D("M_miss_VS_P_n_epCD", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_epCD);
TH2D *h_M_miss_VS_P_miss_epCD = new TH2D("M_miss_P_miss_epCD", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_epCD);

TH2D *h_theta_P_n_P_p_epCD = new TH2D("theta_P_n_P_p_epCD", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_epCD);

TH1D *h_xB_epCD = new TH1D("xB_epCD", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_epCD);

/* Detector responses */
TH1D *h_Edep_epCD = new TH1D("n_Edep_epCD", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_epCD);
TH2D *h_P_n_VS_Edep_epCD = new TH2D("P_n_VS_Edep_epCD", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_epCD);
TH2D *h_P_miss_VS_Edep_epCD = new TH2D("P_miss_VS_Edep_epCD", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_epCD);
TH2D *h_dpp_VS_Edep_epCD = new TH2D("dpp_VS_Edep_epCD", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_epCD);

// Checks on which events have neutrons (Andrew)
// ======================================================================================================================================================================
TH2D *h_xB_mmiss_epFD = new TH2D("xB_mmiss_epFD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
hist_list_2_A.push_back(h_xB_mmiss_epFD);
TH2D *h_xB_mmiss_epnFD = new TH2D("xB_mmiss_epnFD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
hist_list_2_A.push_back(h_xB_mmiss_epnFD);
TH2D *h_xB_mmiss_epn_goodN_pFD = new TH2D("xB_mmiss_epn_goodFD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
hist_list_2_A.push_back(h_xB_mmiss_epn_goodN_pFD);
TH2D *h_xB_mmiss_epn_badN_pFD = new TH2D("xB_mmiss_epn_badFD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
hist_list_2_A.push_back(h_xB_mmiss_epn_badN_pFD);
TH2D *h_xB_mmiss_epCD = new TH2D("xB_mmiss_epCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
hist_list_2_A.push_back(h_xB_mmiss_epCD);
TH2D *h_xB_mmiss_epnCD = new TH2D("xB_mmiss_epnCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
hist_list_2_A.push_back(h_xB_mmiss_epnCD);
TH2D *h_xB_mmiss_epn_goodN_pCD = new TH2D("xB_mmiss_epn_goodCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
hist_list_2_A.push_back(h_xB_mmiss_epn_goodN_pCD);
TH2D *h_xB_mmiss_epn_badN_pCD = new TH2D("xB_mmiss_epn_badCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
hist_list_2_A.push_back(h_xB_mmiss_epn_badN_pCD);

// Step Zero (Andrew)
// ======================================================================================================================================================================

/* Neutron histograms (from Erin) */
TH1D *h_n_multiplicity_allN_epCD_Step0 = new TH1D("n_multiplicity_allN_epCD_Step0", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_allN_epCD_Step0);
TH1D *h_n_multiplicity_goodN_epCD_Step0 = new TH1D("n_multiplicity_goodN_epCD_Step0", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_goodN_epCD_Step0);
TH1D *h_n_multiplicity_badN_epCD_Step0 = new TH1D("n_multiplicity_badN_epCD_Step0", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_badN_epCD_Step0);

/* Kinematical variables */
TH1D *h_theta_n_goodN_Step0 = new TH1D("theta_n_goodN_Step0", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_goodN_Step0);
TH1D *h_phi_n_goodN_Step0 = new TH1D("phi_n_goodN_Step0", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_goodN_Step0);
TH2D *h_theta_n_VS_phi_n_goodN_Step0 = new TH2D("theta_n_VS_phi_n_goodN_Step0", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_goodN_Step0);
TH2D *h_theta_n_VS_beta_n_goodN_Step0 = new TH2D("theta_VS_beta_goodN_Step0", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_goodN_Step0);
TH1D *h_theta_n_badN_Step0 = new TH1D("theta_n_badN_Step0", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_badN_Step0);
TH1D *h_phi_n_badN_Step0 = new TH1D("phi_n_badN_Step0", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_badN_Step0);
TH2D *h_theta_n_VS_phi_n_badN_Step0 = new TH2D("theta_n_VS_phi_n_badN_Step0", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_badN_Step0);
TH2D *h_theta_n_VS_beta_n_badN_Step0 = new TH2D("theta_VS_beta_badN_Step0", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_badN_Step0);

TH1D *h_P_n_goodN_Step0 = new TH1D("P_n_goodN_Step0", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_goodN_Step0);
TH2D *h_P_n_VS_theta_n_goodN_Step0 = new TH2D("P_n_VS_theta_n_goodN_Step0", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_goodN_Step0);
TH1D *h_P_n_badN_Step0 = new TH1D("P_n_badN_Step0", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_badN_Step0);
TH2D *h_P_n_VS_theta_n_badN_Step0 = new TH2D("P_n_VS_theta_n_badN_Step0", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_badN_Step0);

TH1D *h_P_miss_goodN_Step0 = new TH1D("P_miss_goodN_Step0", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_goodN_Step0);
TH2D *h_P_miss_VS_theta_miss_goodN_Step0 = new TH2D("P_miss_VS_theta_miss_goodN_Step0", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_goodN_Step0);
TH1D *h_P_miss_badN_Step0 = new TH1D("P_miss_badN_Step0", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_badN_Step0);
TH2D *h_P_miss_VS_theta_miss_badN_Step0 = new TH2D("P_miss_VS_theta_miss_badN_Step0", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_badN_Step0);

TH1D *h_P_n_minus_P_miss_goodN_Step0 = new TH1D("P_n_minus_P_miss_goodN_Step0", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_goodN_Step0);
TH1D *h_P_n_x_minus_P_miss_x_goodN_Step0 = new TH1D("P_n_x_minus_P_miss_x_goodN_Step0", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_goodN_Step0);
TH1D *h_P_n_y_minus_P_miss_y_goodN_Step0 = new TH1D("P_n_y_minus_P_miss_y_goodN_Step0", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_goodN_Step0);
TH1D *h_P_n_z_minus_P_miss_z_goodN_Step0 = new TH1D("P_n_z_minus_P_miss_z_goodN_Step0", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_goodN_Step0);
TH1D *h_P_n_minus_P_miss_badN_Step0 = new TH1D("P_n_minus_P_miss_badN_Step0", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_badN_Step0);
TH1D *h_P_n_x_minus_P_miss_x_badN_Step0 = new TH1D("P_n_x_minus_P_miss_x_badN_Step0", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_badN_Step0);
TH1D *h_P_n_y_minus_P_miss_y_badN_Step0 = new TH1D("P_n_y_minus_P_miss_y_badN_Step0", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_badN_Step0);
TH1D *h_P_n_z_minus_P_miss_z_badN_Step0 = new TH1D("P_n_z_minus_P_miss_z_badN_Step0", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_badN_Step0);

TH2D *h_P_n_VS_P_miss_goodN_Step0 = new TH2D("P_n_VS_P_miss_goodN_Step0", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_goodN_Step0);
TH2D *h_P_n_x_VS_P_miss_x_goodN_Step0 = new TH2D("P_n_x_VS_P_miss_x_goodN_Step0", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_goodN_Step0);
TH2D *h_P_n_y_VS_P_miss_y_goodN_Step0 = new TH2D("P_n_y_VS_P_miss_y_goodN_Step0", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_goodN_Step0);
TH2D *h_P_n_z_VS_P_miss_z_goodN_Step0 = new TH2D("P_n_z_VS_P_miss_z_goodN_Step0", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_goodN_Step0);
TH2D *h_P_n_VS_P_miss_badN_Step0 = new TH2D("P_n_VS_P_miss_badN_Step0", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_badN_Step0);
TH2D *h_P_n_x_VS_P_miss_x_badN_Step0 = new TH2D("P_n_x_VS_P_miss_x_badN_Step0", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_badN_Step0);
TH2D *h_P_n_y_VS_P_miss_y_badN_Step0 = new TH2D("P_n_y_VS_P_miss_y_badN_Step0", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_badN_Step0);
TH2D *h_P_n_z_VS_P_miss_z_badN_Step0 = new TH2D("P_n_z_VS_P_miss_z_badN_Step0", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_badN_Step0);

TH1D *h_E_p_CD_goodN_Step0 = new TH1D("E_p_CD_goodN_Step0", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_goodN_Step0);
TH1D *h_E_p_FD_goodN_Step0 = new TH1D("E_p_FD_goodN_Step0", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_goodN_Step0);
TH1D *h_E_miss_goodN_Step0 = new TH1D("E_miss_goodN_Step0", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_goodN_Step0);
TH1D *h_M_miss_goodN_Step0 = new TH1D("M_miss_goodN_Step0", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_goodN_Step0);
TH2D *h_M_miss_VS_P_n_goodN_Step0 = new TH2D("M_miss_VS_P_n_goodN_Step0", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_goodN_Step0);
TH2D *h_M_miss_VS_P_miss_goodN_Step0 = new TH2D("M_miss_P_miss_goodN_Step0", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_goodN_Step0);
TH1D *h_E_p_CD_badN_Step0 = new TH1D("E_p_CD_badN_Step0", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_badN_Step0);
TH1D *h_E_p_FD_badN_Step0 = new TH1D("E_p_FD_badN_Step0", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_badN_Step0);
TH1D *h_E_miss_badN_Step0 = new TH1D("E_miss_badN_Step0", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_badN_Step0);
TH1D *h_M_miss_badN_Step0 = new TH1D("M_miss_badN_Step0", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_badN_Step0);
TH2D *h_M_miss_VS_P_n_badN_Step0 = new TH2D("M_miss_VS_P_n_badN_Step0", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_badN_Step0);
TH2D *h_M_miss_VS_P_miss_badN_Step0 = new TH2D("M_miss_P_miss_badN_Step0", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_badN_Step0);

TH2D *h_theta_P_n_P_p_goodN_Step0 = new TH2D("theta_P_n_P_p_goodN_Step0", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_goodN_Step0);
TH2D *h_theta_P_n_P_p_badN_Step0 = new TH2D("theta_P_n_P_p_badN_Step0", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_badN_Step0);

TH1D *h_xB_goodN_Step0 = new TH1D("xB_goodN_Step0", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_goodN_Step0);
TH1D *h_xB_badN_Step0 = new TH1D("xB_badN_Step0", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_badN_Step0);

/* Detector responses */
TH1D *h_Edep_goodN_Step0 = new TH1D("n_Edep_goodN_Step0", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_goodN_Step0);
TH2D *h_P_n_VS_Edep_goodN_Step0 = new TH2D("P_n_VS_Edep_goodN_Step0", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_goodN_Step0);
TH2D *h_P_miss_VS_Edep_goodN_Step0 = new TH2D("P_miss_VS_Edep_goodN_Step0", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_goodN_Step0);
TH2D *h_dpp_VS_Edep_goodN_Step0 = new TH2D("dpp_VS_Edep_goodN_Step0", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_goodN_Step0);
TH1D *h_Edep_badN_Step0 = new TH1D("n_Edep_badN_Step0", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_badN_Step0);
TH2D *h_P_n_VS_Edep_badN_Step0 = new TH2D("P_n_VS_Edep_badN_Step0", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_badN_Step0);
TH2D *h_P_miss_VS_Edep_badN_Step0 = new TH2D("P_miss_VS_Edep_badN_Step0", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_badN_Step0);
TH2D *h_dpp_VS_Edep_badN_Step0 = new TH2D("dpp_VS_Edep_badN_Step0", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_badN_Step0);

TH2D *h_pnRes_theta_nmiss_Step0 = new TH2D("pnRes_theta_nmiss_Step0", "(P_{miss}-P_{n})/P_{miss} vs. #theta_{n,miss};(P_{miss}-P_{n})/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step0);

TH1D *h_ToF_goodN_Step0 = new TH1D("ToF_goodN_Step0", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_goodN_Step0);
TH1D *h_ToF_badN_Step0 = new TH1D("ToF_badN_Step0", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_badN_Step0);

TH1D *h_beta_goodN_Step0 = new TH1D("beta_goodN_Step0", "#beta of CND Neutrons;#beta;Counts", 50, 0, 1.1);
hist_list_1_A.push_back(h_beta_goodN_Step0);
TH1D *h_beta_badN_Step0 = new TH1D("beta_badN_Step0", "#beta of CND Neutrons;#beta;Counts", 50, 0, 1.1);
hist_list_1_A.push_back(h_beta_badN_Step0);

// TH1D *h_Edep_goodN_Step0 = new TH1D("Edep_goodN_Step0", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 100);
// hist_list_1_A.push_back(h_Edep_goodN_Step0);
// TH1D *h_Edep_badN_Step0 = new TH1D("Edep_badN_Step0", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 100);
// hist_list_1_A.push_back(h_Edep_badN_Step0);

TH2D *h_beta_Edep_goodN_Step0 = new TH2D("Edep_beta_goodN_Step0", "#beta vs. E_{dep} of CND Neutrons;#beta;E_{dep} [MeV]", 50, 0, 1.1, 50, 0, 100);
hist_list_2_A.push_back(h_beta_Edep_goodN_Step0);
TH2D *h_beta_Edep_badN_Step0 = new TH2D("Edep_beta_badN_Step0", "#beta vs. E_{dep} of CND Neutrons;#beta;E_{dep} [MeV]", 50, 0, 1.1, 50, 0, 100);
hist_list_2_A.push_back(h_beta_Edep_badN_Step0);

TH2D *h_Edep_ToF_goodN_Step0 = new TH2D("Edep_ToF_goodN_Step0", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_goodN_Step0);
TH2D *h_Edep_ToF_badN_Step0 = new TH2D("Edep_ToF_badN_Step0", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_badN_Step0);

// Step One (After Beta Cut) (Andrew)
// ======================================================================================================================================================================

/* Neutron histograms (from Erin) */
TH1D *h_n_multiplicity_allN_epCD_Step1 = new TH1D("n_multiplicity_allN_epCD_Step1", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_allN_epCD_Step1);
TH1D *h_n_multiplicity_goodN_epCD_Step1 = new TH1D("n_multiplicity_goodN_epCD_Step1", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_goodN_epCD_Step1);
TH1D *h_n_multiplicity_badN_epCD_Step1 = new TH1D("n_multiplicity_badN_epCD_Step1", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_badN_epCD_Step1);

/* Kinematical variables */
TH1D *h_theta_n_goodN_Step1 = new TH1D("theta_n_goodN_Step1", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_goodN_Step1);
TH1D *h_phi_n_goodN_Step1 = new TH1D("phi_n_goodN_Step1", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_goodN_Step1);
TH2D *h_theta_n_VS_phi_n_goodN_Step1 = new TH2D("theta_n_VS_phi_n_goodN_Step1", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_goodN_Step1);
TH2D *h_theta_n_VS_beta_n_goodN_Step1 = new TH2D("theta_VS_beta_goodN_Step1", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_goodN_Step1);
TH1D *h_theta_n_badN_Step1 = new TH1D("theta_n_badN_Step1", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_badN_Step1);
TH1D *h_phi_n_badN_Step1 = new TH1D("phi_n_badN_Step1", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_badN_Step1);
TH2D *h_theta_n_VS_phi_n_badN_Step1 = new TH2D("theta_n_VS_phi_n_badN_Step1", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_badN_Step1);
TH2D *h_theta_n_VS_beta_n_badN_Step1 = new TH2D("theta_VS_beta_badN_Step1", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_badN_Step1);

TH1D *h_P_n_goodN_Step1 = new TH1D("P_n_goodN_Step1", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_goodN_Step1);
TH2D *h_P_n_VS_theta_n_goodN_Step1 = new TH2D("P_n_VS_theta_n_goodN_Step1", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_goodN_Step1);
TH1D *h_P_n_badN_Step1 = new TH1D("P_n_badN_Step1", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_badN_Step1);
TH2D *h_P_n_VS_theta_n_badN_Step1 = new TH2D("P_n_VS_theta_n_badN_Step1", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_badN_Step1);

TH1D *h_P_miss_goodN_Step1 = new TH1D("P_miss_goodN_Step1", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_goodN_Step1);
TH2D *h_P_miss_VS_theta_miss_goodN_Step1 = new TH2D("P_miss_VS_theta_miss_goodN_Step1", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_goodN_Step1);
TH1D *h_P_miss_badN_Step1 = new TH1D("P_miss_badN_Step1", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_badN_Step1);
TH2D *h_P_miss_VS_theta_miss_badN_Step1 = new TH2D("P_miss_VS_theta_miss_badN_Step1", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_badN_Step1);

TH1D *h_P_n_minus_P_miss_goodN_Step1 = new TH1D("P_n_minus_P_miss_goodN_Step1", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_goodN_Step1);
TH1D *h_P_n_x_minus_P_miss_x_goodN_Step1 = new TH1D("P_n_x_minus_P_miss_x_goodN_Step1", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_goodN_Step1);
TH1D *h_P_n_y_minus_P_miss_y_goodN_Step1 = new TH1D("P_n_y_minus_P_miss_y_goodN_Step1", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_goodN_Step1);
TH1D *h_P_n_z_minus_P_miss_z_goodN_Step1 = new TH1D("P_n_z_minus_P_miss_z_goodN_Step1", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_goodN_Step1);
TH1D *h_P_n_minus_P_miss_badN_Step1 = new TH1D("P_n_minus_P_miss_badN_Step1", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_badN_Step1);
TH1D *h_P_n_x_minus_P_miss_x_badN_Step1 = new TH1D("P_n_x_minus_P_miss_x_badN_Step1", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_badN_Step1);
TH1D *h_P_n_y_minus_P_miss_y_badN_Step1 = new TH1D("P_n_y_minus_P_miss_y_badN_Step1", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_badN_Step1);
TH1D *h_P_n_z_minus_P_miss_z_badN_Step1 = new TH1D("P_n_z_minus_P_miss_z_badN_Step1", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_badN_Step1);

TH2D *h_P_n_VS_P_miss_goodN_Step1 = new TH2D("P_n_VS_P_miss_goodN_Step1", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_goodN_Step1);
TH2D *h_P_n_x_VS_P_miss_x_goodN_Step1 = new TH2D("P_n_x_VS_P_miss_x_goodN_Step1", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_goodN_Step1);
TH2D *h_P_n_y_VS_P_miss_y_goodN_Step1 = new TH2D("P_n_y_VS_P_miss_y_goodN_Step1", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_goodN_Step1);
TH2D *h_P_n_z_VS_P_miss_z_goodN_Step1 = new TH2D("P_n_z_VS_P_miss_z_goodN_Step1", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_goodN_Step1);
TH2D *h_P_n_VS_P_miss_badN_Step1 = new TH2D("P_n_VS_P_miss_badN_Step1", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_badN_Step1);
TH2D *h_P_n_x_VS_P_miss_x_badN_Step1 = new TH2D("P_n_x_VS_P_miss_x_badN_Step1", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_badN_Step1);
TH2D *h_P_n_y_VS_P_miss_y_badN_Step1 = new TH2D("P_n_y_VS_P_miss_y_badN_Step1", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_badN_Step1);
TH2D *h_P_n_z_VS_P_miss_z_badN_Step1 = new TH2D("P_n_z_VS_P_miss_z_badN_Step1", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_badN_Step1);

TH1D *h_E_p_CD_goodN_Step1 = new TH1D("E_p_CD_goodN_Step1", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_goodN_Step1);
TH1D *h_E_p_FD_goodN_Step1 = new TH1D("E_p_FD_goodN_Step1", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_goodN_Step1);
TH1D *h_E_miss_goodN_Step1 = new TH1D("E_miss_goodN_Step1", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_goodN_Step1);
TH1D *h_M_miss_goodN_Step1 = new TH1D("M_miss_goodN_Step1", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_goodN_Step1);
TH2D *h_M_miss_VS_P_n_goodN_Step1 = new TH2D("M_miss_VS_P_n_goodN_Step1", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_goodN_Step1);
TH2D *h_M_miss_VS_P_miss_goodN_Step1 = new TH2D("M_miss_P_miss_goodN_Step1", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_goodN_Step1);
TH1D *h_E_p_CD_badN_Step1 = new TH1D("E_p_CD_badN_Step1", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_badN_Step1);
TH1D *h_E_p_FD_badN_Step1 = new TH1D("E_p_FD_badN_Step1", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_badN_Step1);
TH1D *h_E_miss_badN_Step1 = new TH1D("E_miss_badN_Step1", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_badN_Step1);
TH1D *h_M_miss_badN_Step1 = new TH1D("M_miss_badN_Step1", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_badN_Step1);
TH2D *h_M_miss_VS_P_n_badN_Step1 = new TH2D("M_miss_VS_P_n_badN_Step1", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_badN_Step1);
TH2D *h_M_miss_VS_P_miss_badN_Step1 = new TH2D("M_miss_P_miss_badN_Step1", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_badN_Step1);

TH2D *h_theta_P_n_P_p_goodN_Step1 = new TH2D("theta_P_n_P_p_goodN_Step1", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_goodN_Step1);
TH2D *h_theta_P_n_P_p_badN_Step1 = new TH2D("theta_P_n_P_p_badN_Step1", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_badN_Step1);

TH1D *h_xB_goodN_Step1 = new TH1D("xB_goodN_Step1", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_goodN_Step1);
TH1D *h_xB_badN_Step1 = new TH1D("xB_badN_Step1", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_badN_Step1);

/* Detector responses */
TH1D *h_Edep_goodN_Step1 = new TH1D("n_Edep_goodN_Step1", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_goodN_Step1);
TH2D *h_P_n_VS_Edep_goodN_Step1 = new TH2D("P_n_VS_Edep_goodN_Step1", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_goodN_Step1);
TH2D *h_P_miss_VS_Edep_goodN_Step1 = new TH2D("P_miss_VS_Edep_goodN_Step1", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_goodN_Step1);
TH2D *h_dpp_VS_Edep_goodN_Step1 = new TH2D("dpp_VS_Edep_goodN_Step1", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_goodN_Step1);
TH1D *h_Edep_badN_Step1 = new TH1D("n_Edep_badN_Step1", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_badN_Step1);
TH2D *h_P_n_VS_Edep_badN_Step1 = new TH2D("P_n_VS_Edep_badN_Step1", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_badN_Step1);
TH2D *h_P_miss_VS_Edep_badN_Step1 = new TH2D("P_miss_VS_Edep_badN_Step1", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_badN_Step1);
TH2D *h_dpp_VS_Edep_badN_Step1 = new TH2D("dpp_VS_Edep_badN_Step1", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_badN_Step1);

TH2D *h_pnRes_theta_nmiss_Step1 = new TH2D("pnRes_theta_nmiss_Step1", "(P_{miss}-P_{n})/P_{miss} vs. #theta_{n,miss};(P_{miss}-P_{n})/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step1);

TH1D *h_pmiss_goodN_Step1 = new TH1D("pmiss_goodN_Step1", "P_{miss} Step1;P_{miss} [GeV/c];Counts", 25, 0.25, 1.0);
hist_list_1_A.push_back(h_pmiss_goodN_Step1);
TH1D *h_pmiss_badN_Step1 = new TH1D("pmiss_badN_Step1", "P_{miss} Step1;P_{miss} [GeV/c];Counts", 25, 0.25, 1.0);
hist_list_1_A.push_back(h_pmiss_badN_Step1);

TH1D *h_ToF_goodN_Step1 = new TH1D("ToF_goodN_Step1", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_goodN_Step1);
TH1D *h_ToF_badN_Step1 = new TH1D("ToF_badN_Step1", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_badN_Step1);

TH1D *h_edep_goodN_Step1 = new TH1D("edep_goodN_Step1", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 100, 0, 50);
hist_list_1_A.push_back(h_edep_goodN_Step1);
TH1D *h_edep_badN_Step1 = new TH1D("edep_badN_Step1", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 100, 0, 50);
hist_list_1_A.push_back(h_edep_badN_Step1);

TH2D *h_Edep_ToF_goodN_Step1 = new TH2D("Edep_ToF_goodN_Step1", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_goodN_Step1);
TH2D *h_Edep_ToF_badN_Step1 = new TH2D("Edep_ToF_badN_Step1", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_badN_Step1);

TH1D *h_edep_over_edepCTOT_goodN_Step1 = new TH1D("edep_over_edepCTOT_goodN_Step1", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 100, 0, 5);
hist_list_1_A.push_back(h_edep_over_edepCTOT_goodN_Step1);
TH1D *h_edep_over_edepCTOT_badN_Step1 = new TH1D("edep_over_edepCTOT_badN_Step1", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 100, 0, 5);
hist_list_1_A.push_back(h_edep_over_edepCTOT_badN_Step1);

TH1D *h_edep_goodN_withNearbyPos_Step1 = new TH1D("edep_goodN_withNearbyPos_Step1", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 100, 0, 50);
hist_list_1_A.push_back(h_edep_goodN_withNearbyPos_Step1);
TH1D *h_edep_badN_withNearbyPos_Step1 = new TH1D("edep_badN_withNearbyPos_Step1", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 100, 0, 50);
hist_list_1_A.push_back(h_edep_badN_withNearbyPos_Step1);

TH1D *h_sdiff_pos_goodN_Step1_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_goodN_Step1_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
    h_sdiff_pos_goodN_Step1_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
    hist_list_1_A.push_back(h_sdiff_pos_goodN_Step1_layer[k]);
}

TH1D *h_sdiff_pos_badN_Step1_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_badN_Step1_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
    h_sdiff_pos_badN_Step1_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
    hist_list_1_A.push_back(h_sdiff_pos_badN_Step1_layer[k]);
}

TH2D *h_sdiff_pos_mom_goodN_Step1_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_mom_goodN_Step1_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
    h_sdiff_pos_mom_goodN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
    hist_list_2_A.push_back(h_sdiff_pos_mom_goodN_Step1_layer[k]);
}

TH2D *h_sdiff_pos_mom_badN_Step1_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_mom_badN_Step1_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
    h_sdiff_pos_mom_badN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
    hist_list_2_A.push_back(h_sdiff_pos_mom_badN_Step1_layer[k]);
}

TH2D *h_sdiff_pos_z_goodN_Step1_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_z_goodN_Step1_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
    h_sdiff_pos_z_goodN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
    hist_list_2_A.push_back(h_sdiff_pos_z_goodN_Step1_layer[k]);
}

TH2D *h_sdiff_pos_z_badN_Step1_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_z_badN_Step1_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
    h_sdiff_pos_z_badN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
    hist_list_2_A.push_back(h_sdiff_pos_z_badN_Step1_layer[k]);
}

TH2D *h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_goodN_Step1_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
    h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
    hist_list_2_A.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[k]);
}

TH2D *h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_badN_Step1_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
    h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
    hist_list_2_A.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[k]);
}

TH2D *h_diff_ToFc_z_Edep_noNear_goodN_Step1 = new TH2D("diff_ToFc_z_Edep_noNear_goodN_Step1", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
hist_list_2_A.push_back(h_diff_ToFc_z_Edep_noNear_goodN_Step1);
TH2D *h_diff_ToFc_z_Edep_noNear_badN_Step1 = new TH2D("diff_ToFc_z_Edep_noNear_badN_Step1", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
hist_list_2_A.push_back(h_diff_ToFc_z_Edep_noNear_badN_Step1);

TH2D *h_diff_ToFc_z_Edep_yesNear_goodN_Step1 = new TH2D("diff_ToFc_z_Edep_yesNear_goodN_Step1", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
hist_list_2_A.push_back(h_diff_ToFc_z_Edep_yesNear_goodN_Step1);
TH2D *h_diff_ToFc_z_Edep_yesNear_badN_Step1 = new TH2D("diff_ToFc_z_Edep_yesNear_badN_Step1", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
hist_list_2_A.push_back(h_diff_ToFc_z_Edep_yesNear_badN_Step1);

// Step Two (After applying Phi Diff Charge Track cut) (Andrew)
// ======================================================================================================================================================================

/* Neutron histograms (from Erin) */
TH1D *h_n_multiplicity_allN_epCD_Step2 = new TH1D("n_multiplicity_allN_epCD_Step2", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_allN_epCD_Step2);
TH1D *h_n_multiplicity_goodN_epCD_Step2 = new TH1D("n_multiplicity_goodN_epCD_Step2", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_goodN_epCD_Step2);
TH1D *h_n_multiplicity_badN_epCD_Step2 = new TH1D("n_multiplicity_badN_epCD_Step2", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_badN_epCD_Step2);

/* Kinematical variables */
TH1D *h_theta_n_goodN_Step2 = new TH1D("theta_n_goodN_Step2", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_goodN_Step2);
TH1D *h_phi_n_goodN_Step2 = new TH1D("phi_n_goodN_Step2", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_goodN_Step2);
TH2D *h_theta_n_VS_phi_n_goodN_Step2 = new TH2D("theta_n_VS_phi_n_goodN_Step2", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_goodN_Step2);
TH2D *h_theta_n_VS_beta_n_goodN_Step2 = new TH2D("theta_VS_beta_goodN_Step2", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_goodN_Step2);
TH1D *h_theta_n_badN_Step2 = new TH1D("theta_n_badN_Step2", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_badN_Step2);
TH1D *h_phi_n_badN_Step2 = new TH1D("phi_n_badN_Step2", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_badN_Step2);
TH2D *h_theta_n_VS_phi_n_badN_Step2 = new TH2D("theta_n_VS_phi_n_badN_Step2", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_badN_Step2);
TH2D *h_theta_n_VS_beta_n_badN_Step2 = new TH2D("theta_VS_beta_badN_Step2", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_badN_Step2);

TH1D *h_P_n_goodN_Step2 = new TH1D("P_n_goodN_Step2", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_goodN_Step2);
TH2D *h_P_n_VS_theta_n_goodN_Step2 = new TH2D("P_n_VS_theta_n_goodN_Step2", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_goodN_Step2);
TH1D *h_P_n_badN_Step2 = new TH1D("P_n_badN_Step2", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_badN_Step2);
TH2D *h_P_n_VS_theta_n_badN_Step2 = new TH2D("P_n_VS_theta_n_badN_Step2", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_badN_Step2);

TH1D *h_P_miss_goodN_Step2 = new TH1D("P_miss_goodN_Step2", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_goodN_Step2);
TH2D *h_P_miss_VS_theta_miss_goodN_Step2 = new TH2D("P_miss_VS_theta_miss_goodN_Step2", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_goodN_Step2);
TH1D *h_P_miss_badN_Step2 = new TH1D("P_miss_badN_Step2", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_badN_Step2);
TH2D *h_P_miss_VS_theta_miss_badN_Step2 = new TH2D("P_miss_VS_theta_miss_badN_Step2", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_badN_Step2);

TH1D *h_P_n_minus_P_miss_goodN_Step2 = new TH1D("P_n_minus_P_miss_goodN_Step2", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_goodN_Step2);
TH1D *h_P_n_x_minus_P_miss_x_goodN_Step2 = new TH1D("P_n_x_minus_P_miss_x_goodN_Step2", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_goodN_Step2);
TH1D *h_P_n_y_minus_P_miss_y_goodN_Step2 = new TH1D("P_n_y_minus_P_miss_y_goodN_Step2", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_goodN_Step2);
TH1D *h_P_n_z_minus_P_miss_z_goodN_Step2 = new TH1D("P_n_z_minus_P_miss_z_goodN_Step2", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_goodN_Step2);
TH1D *h_P_n_minus_P_miss_badN_Step2 = new TH1D("P_n_minus_P_miss_badN_Step2", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_badN_Step2);
TH1D *h_P_n_x_minus_P_miss_x_badN_Step2 = new TH1D("P_n_x_minus_P_miss_x_badN_Step2", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_badN_Step2);
TH1D *h_P_n_y_minus_P_miss_y_badN_Step2 = new TH1D("P_n_y_minus_P_miss_y_badN_Step2", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_badN_Step2);
TH1D *h_P_n_z_minus_P_miss_z_badN_Step2 = new TH1D("P_n_z_minus_P_miss_z_badN_Step2", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_badN_Step2);

TH2D *h_P_n_VS_P_miss_goodN_Step2 = new TH2D("P_n_VS_P_miss_goodN_Step2", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_goodN_Step2);
TH2D *h_P_n_x_VS_P_miss_x_goodN_Step2 = new TH2D("P_n_x_VS_P_miss_x_goodN_Step2", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_goodN_Step2);
TH2D *h_P_n_y_VS_P_miss_y_goodN_Step2 = new TH2D("P_n_y_VS_P_miss_y_goodN_Step2", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_goodN_Step2);
TH2D *h_P_n_z_VS_P_miss_z_goodN_Step2 = new TH2D("P_n_z_VS_P_miss_z_goodN_Step2", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_goodN_Step2);
TH2D *h_P_n_VS_P_miss_badN_Step2 = new TH2D("P_n_VS_P_miss_badN_Step2", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_badN_Step2);
TH2D *h_P_n_x_VS_P_miss_x_badN_Step2 = new TH2D("P_n_x_VS_P_miss_x_badN_Step2", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_badN_Step2);
TH2D *h_P_n_y_VS_P_miss_y_badN_Step2 = new TH2D("P_n_y_VS_P_miss_y_badN_Step2", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_badN_Step2);
TH2D *h_P_n_z_VS_P_miss_z_badN_Step2 = new TH2D("P_n_z_VS_P_miss_z_badN_Step2", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_badN_Step2);

TH1D *h_E_p_CD_goodN_Step2 = new TH1D("E_p_CD_goodN_Step2", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_goodN_Step2);
TH1D *h_E_p_FD_goodN_Step2 = new TH1D("E_p_FD_goodN_Step2", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_goodN_Step2);
TH1D *h_E_miss_goodN_Step2 = new TH1D("E_miss_goodN_Step2", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_goodN_Step2);
TH1D *h_M_miss_goodN_Step2 = new TH1D("M_miss_goodN_Step2", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_goodN_Step2);
TH2D *h_M_miss_VS_P_n_goodN_Step2 = new TH2D("M_miss_VS_P_n_goodN_Step2", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_goodN_Step2);
TH2D *h_M_miss_VS_P_miss_goodN_Step2 = new TH2D("M_miss_P_miss_goodN_Step2", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_goodN_Step2);
TH1D *h_E_p_CD_badN_Step2 = new TH1D("E_p_CD_badN_Step2", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_badN_Step2);
TH1D *h_E_p_FD_badN_Step2 = new TH1D("E_p_FD_badN_Step2", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_badN_Step2);
TH1D *h_E_miss_badN_Step2 = new TH1D("E_miss_badN_Step2", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_badN_Step2);
TH1D *h_M_miss_badN_Step2 = new TH1D("M_miss_badN_Step2", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_badN_Step2);
TH2D *h_M_miss_VS_P_n_badN_Step2 = new TH2D("M_miss_VS_P_n_badN_Step2", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_badN_Step2);
TH2D *h_M_miss_VS_P_miss_badN_Step2 = new TH2D("M_miss_P_miss_badN_Step2", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_badN_Step2);

TH2D *h_theta_P_n_P_p_goodN_Step2 = new TH2D("theta_P_n_P_p_goodN_Step2", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_goodN_Step2);
TH2D *h_theta_P_n_P_p_badN_Step2 = new TH2D("theta_P_n_P_p_badN_Step2", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_badN_Step2);

TH1D *h_xB_goodN_Step2 = new TH1D("xB_goodN_Step2", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_goodN_Step2);
TH1D *h_xB_badN_Step2 = new TH1D("xB_badN_Step2", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_badN_Step2);

/* Detector responses */
TH1D *h_Edep_goodN_Step2 = new TH1D("n_Edep_goodN_Step2", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_goodN_Step2);
TH2D *h_P_n_VS_Edep_goodN_Step2 = new TH2D("P_n_VS_Edep_goodN_Step2", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_goodN_Step2);
TH2D *h_P_miss_VS_Edep_goodN_Step2 = new TH2D("P_miss_VS_Edep_goodN_Step2", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_goodN_Step2);
TH2D *h_dpp_VS_Edep_goodN_Step2 = new TH2D("dpp_VS_Edep_goodN_Step2", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_goodN_Step2);
TH1D *h_Edep_badN_Step2 = new TH1D("n_Edep_badN_Step2", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_badN_Step2);
TH2D *h_P_n_VS_Edep_badN_Step2 = new TH2D("P_n_VS_Edep_badN_Step2", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_badN_Step2);
TH2D *h_P_miss_VS_Edep_badN_Step2 = new TH2D("P_miss_VS_Edep_badN_Step2", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_badN_Step2);
TH2D *h_dpp_VS_Edep_badN_Step2 = new TH2D("dpp_VS_Edep_badN_Step2", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_badN_Step2);

TH2D *h_pnRes_theta_nmiss_Step2 = new TH2D("pnRes_theta_nmiss_Step2", "(P_{miss}-P_{n})/P_{miss} vs. #theta_{n,miss};(P_{miss}-P_{n})/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step2);

TH1D *h_ToF_goodN_Step2 = new TH1D("ToF_goodN_Step2", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_goodN_Step2);
TH1D *h_ToF_badN_Step2 = new TH1D("ToF_badN_Step2", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_badN_Step2);

TH2D *h_Edep_ToF_goodN_Step2 = new TH2D("Edep_ToF_goodN_Step2", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_goodN_Step2);
TH2D *h_Edep_ToF_badN_Step2 = new TH2D("Edep_ToF_badN_Step2", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_badN_Step2);

TH1D *h_sdiff_pos_goodN_Step2_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_goodN_Step2_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
    h_sdiff_pos_goodN_Step2_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
    hist_list_1_A.push_back(h_sdiff_pos_goodN_Step2_layer[k]);
}

TH1D *h_sdiff_pos_badN_Step2_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_pos_badN_Step2_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
    h_sdiff_pos_badN_Step2_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
    hist_list_1_A.push_back(h_sdiff_pos_badN_Step2_layer[k]);
}

TH1D *h_sdiff_allhit_goodN_Step2_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_allhit_goodN_Step2_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus Random Hit Sector (Layer Difference = %d)", k - 3);
    h_sdiff_allhit_goodN_Step2_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
    hist_list_1_A.push_back(h_sdiff_allhit_goodN_Step2_layer[k]);
}

TH2D *h_sdiff_ldiff_allhit_goodN_Step2 = new TH2D("sdiff_ldiff_allhit_goodN_Step2", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
hist_list_2_A.push_back(h_sdiff_ldiff_allhit_goodN_Step2);

TH1D *h_sdiff_allhit_badN_Step2_layer[7];

for (int k = 0; k < 7; k++)
{
    sprintf(temp_name_A, "sdiff_allhit_badN_Step2_layer_%d", k - 3);
    sprintf(temp_title_A, "Nuetral Sector minus Random Hit Sector (Layer Difference = %d)", k - 3);
    h_sdiff_allhit_badN_Step2_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
    hist_list_1_A.push_back(h_sdiff_allhit_badN_Step2_layer[k]);
}

TH2D *h_sdiff_ldiff_allhit_badN_Step2 = new TH2D("sdiff_ldiff_allhit_badN_Step2", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
hist_list_2_A.push_back(h_sdiff_ldiff_allhit_badN_Step2);

TH1D *h_numberNearby_goodN_Step2 = new TH1D("numberNearby_goodN_Step2", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
hist_list_1_A.push_back(h_numberNearby_goodN_Step2);
TH1D *h_numberNearby_badN_Step2 = new TH1D("numberNearby_badN_Step2", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
hist_list_1_A.push_back(h_numberNearby_badN_Step2);

TH2D *h_numberNearby_momN_goodN_Step2 = new TH2D("numberNearby_momN_goodN_Step2", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
hist_list_2_A.push_back(h_numberNearby_momN_goodN_Step2);
TH2D *h_numberNearby_momN_badN_Step2 = new TH2D("numberNearby_momN_badN_Step2", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
hist_list_2_A.push_back(h_numberNearby_momN_badN_Step2);

TH1D *h_NearbyEdep_goodN_Step2 = new TH1D("NearbyEdep_goodN_Step2", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_NearbyEdep_goodN_Step2);
TH1D *h_NearbyEdep_badN_Step2 = new TH1D("NearbyEdep_badN_Step2", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_NearbyEdep_badN_Step2);

TH1D *h_nsector_goodN_Step2 = new TH1D("nsector_goodN_Step2", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
hist_list_1_A.push_back(h_nsector_goodN_Step2);
TH1D *h_nsector_badN_Step2 = new TH1D("nsector_badN_Step2", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
hist_list_1_A.push_back(h_nsector_badN_Step2);

// Step Three (After applying Phi Diff Charge Track cut) (Andrew)
// ======================================================================================================================================================================

/* Neutron histograms (from Erin) */
TH1D *h_n_multiplicity_allN_epCD_Step3 = new TH1D("n_multiplicity_allN_epCD_Step3", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_allN_epCD_Step3);
TH1D *h_n_multiplicity_goodN_epCD_Step3 = new TH1D("n_multiplicity_goodN_epCD_Step3", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_goodN_epCD_Step3);
TH1D *h_n_multiplicity_badN_epCD_Step3 = new TH1D("n_multiplicity_badN_epCD_Step3", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_badN_epCD_Step3);

/* Kinematical variables */
TH1D *h_theta_n_goodN_Step3 = new TH1D("theta_n_goodN_Step3", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_goodN_Step3);
TH1D *h_phi_n_goodN_Step3 = new TH1D("phi_n_goodN_Step3", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_goodN_Step3);
TH2D *h_theta_n_VS_phi_n_goodN_Step3 = new TH2D("theta_n_VS_phi_n_goodN_Step3", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_goodN_Step3);
TH2D *h_theta_n_VS_beta_n_goodN_Step3 = new TH2D("theta_VS_beta_goodN_Step3", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_goodN_Step3);
TH1D *h_theta_n_badN_Step3 = new TH1D("theta_n_badN_Step3", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_badN_Step3);
TH1D *h_phi_n_badN_Step3 = new TH1D("phi_n_badN_Step3", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_badN_Step3);
TH2D *h_theta_n_VS_phi_n_badN_Step3 = new TH2D("theta_n_VS_phi_n_badN_Step3", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_badN_Step3);
TH2D *h_theta_n_VS_beta_n_badN_Step3 = new TH2D("theta_VS_beta_badN_Step3", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_badN_Step3);

TH1D *h_P_n_goodN_Step3 = new TH1D("P_n_goodN_Step3", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_goodN_Step3);
TH2D *h_P_n_VS_theta_n_goodN_Step3 = new TH2D("P_n_VS_theta_n_goodN_Step3", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_goodN_Step3);
TH1D *h_P_n_badN_Step3 = new TH1D("P_n_badN_Step3", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_badN_Step3);
TH2D *h_P_n_VS_theta_n_badN_Step3 = new TH2D("P_n_VS_theta_n_badN_Step3", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_badN_Step3);

TH1D *h_P_miss_goodN_Step3 = new TH1D("P_miss_goodN_Step3", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_goodN_Step3);
TH2D *h_P_miss_VS_theta_miss_goodN_Step3 = new TH2D("P_miss_VS_theta_miss_goodN_Step3", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_goodN_Step3);
TH1D *h_P_miss_badN_Step3 = new TH1D("P_miss_badN_Step3", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_badN_Step3);
TH2D *h_P_miss_VS_theta_miss_badN_Step3 = new TH2D("P_miss_VS_theta_miss_badN_Step3", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_badN_Step3);

TH1D *h_P_n_minus_P_miss_goodN_Step3 = new TH1D("P_n_minus_P_miss_goodN_Step3", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_goodN_Step3);
TH1D *h_P_n_x_minus_P_miss_x_goodN_Step3 = new TH1D("P_n_x_minus_P_miss_x_goodN_Step3", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_goodN_Step3);
TH1D *h_P_n_y_minus_P_miss_y_goodN_Step3 = new TH1D("P_n_y_minus_P_miss_y_goodN_Step3", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_goodN_Step3);
TH1D *h_P_n_z_minus_P_miss_z_goodN_Step3 = new TH1D("P_n_z_minus_P_miss_z_goodN_Step3", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_goodN_Step3);
TH1D *h_P_n_minus_P_miss_badN_Step3 = new TH1D("P_n_minus_P_miss_badN_Step3", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_badN_Step3);
TH1D *h_P_n_x_minus_P_miss_x_badN_Step3 = new TH1D("P_n_x_minus_P_miss_x_badN_Step3", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_badN_Step3);
TH1D *h_P_n_y_minus_P_miss_y_badN_Step3 = new TH1D("P_n_y_minus_P_miss_y_badN_Step3", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_badN_Step3);
TH1D *h_P_n_z_minus_P_miss_z_badN_Step3 = new TH1D("P_n_z_minus_P_miss_z_badN_Step3", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_badN_Step3);

TH2D *h_P_n_VS_P_miss_goodN_Step3 = new TH2D("P_n_VS_P_miss_goodN_Step3", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_goodN_Step3);
TH2D *h_P_n_x_VS_P_miss_x_goodN_Step3 = new TH2D("P_n_x_VS_P_miss_x_goodN_Step3", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_goodN_Step3);
TH2D *h_P_n_y_VS_P_miss_y_goodN_Step3 = new TH2D("P_n_y_VS_P_miss_y_goodN_Step3", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_goodN_Step3);
TH2D *h_P_n_z_VS_P_miss_z_goodN_Step3 = new TH2D("P_n_z_VS_P_miss_z_goodN_Step3", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_goodN_Step3);
TH2D *h_P_n_VS_P_miss_badN_Step3 = new TH2D("P_n_VS_P_miss_badN_Step3", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_badN_Step3);
TH2D *h_P_n_x_VS_P_miss_x_badN_Step3 = new TH2D("P_n_x_VS_P_miss_x_badN_Step3", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_badN_Step3);
TH2D *h_P_n_y_VS_P_miss_y_badN_Step3 = new TH2D("P_n_y_VS_P_miss_y_badN_Step3", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_badN_Step3);
TH2D *h_P_n_z_VS_P_miss_z_badN_Step3 = new TH2D("P_n_z_VS_P_miss_z_badN_Step3", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_badN_Step3);

TH1D *h_E_p_CD_goodN_Step3 = new TH1D("E_p_CD_goodN_Step3", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_goodN_Step3);
TH1D *h_E_p_FD_goodN_Step3 = new TH1D("E_p_FD_goodN_Step3", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_goodN_Step3);
TH1D *h_E_miss_goodN_Step3 = new TH1D("E_miss_goodN_Step3", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_goodN_Step3);
TH1D *h_M_miss_goodN_Step3 = new TH1D("M_miss_goodN_Step3", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_goodN_Step3);
TH2D *h_M_miss_VS_P_n_goodN_Step3 = new TH2D("M_miss_VS_P_n_goodN_Step3", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_goodN_Step3);
TH2D *h_M_miss_VS_P_miss_goodN_Step3 = new TH2D("M_miss_P_miss_goodN_Step3", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_goodN_Step3);
TH1D *h_E_p_CD_badN_Step3 = new TH1D("E_p_CD_badN_Step3", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_badN_Step3);
TH1D *h_E_p_FD_badN_Step3 = new TH1D("E_p_FD_badN_Step3", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_badN_Step3);
TH1D *h_E_miss_badN_Step3 = new TH1D("E_miss_badN_Step3", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_badN_Step3);
TH1D *h_M_miss_badN_Step3 = new TH1D("M_miss_badN_Step3", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_badN_Step3);
TH2D *h_M_miss_VS_P_n_badN_Step3 = new TH2D("M_miss_VS_P_n_badN_Step3", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_badN_Step3);
TH2D *h_M_miss_VS_P_miss_badN_Step3 = new TH2D("M_miss_P_miss_badN_Step3", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_badN_Step3);

TH2D *h_theta_P_n_P_p_goodN_Step3 = new TH2D("theta_P_n_P_p_goodN_Step3", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_goodN_Step3);
TH2D *h_theta_P_n_P_p_badN_Step3 = new TH2D("theta_P_n_P_p_badN_Step3", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_badN_Step3);

TH1D *h_xB_goodN_Step3 = new TH1D("xB_goodN_Step3", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_goodN_Step3);
TH1D *h_xB_badN_Step3 = new TH1D("xB_badN_Step3", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_badN_Step3);

/* Detector responses */
TH1D *h_Edep_goodN_Step3 = new TH1D("n_Edep_goodN_Step3", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_goodN_Step3);
TH2D *h_P_n_VS_Edep_goodN_Step3 = new TH2D("P_n_VS_Edep_goodN_Step3", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_goodN_Step3);
TH2D *h_P_miss_VS_Edep_goodN_Step3 = new TH2D("P_miss_VS_Edep_goodN_Step3", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_goodN_Step3);
TH2D *h_dpp_VS_Edep_goodN_Step3 = new TH2D("dpp_VS_Edep_goodN_Step3", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_goodN_Step3);
TH1D *h_Edep_badN_Step3 = new TH1D("n_Edep_badN_Step3", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_badN_Step3);
TH2D *h_P_n_VS_Edep_badN_Step3 = new TH2D("P_n_VS_Edep_badN_Step3", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_badN_Step3);
TH2D *h_P_miss_VS_Edep_badN_Step3 = new TH2D("P_miss_VS_Edep_badN_Step3", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_badN_Step3);
TH2D *h_dpp_VS_Edep_badN_Step3 = new TH2D("dpp_VS_Edep_badN_Step3", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_badN_Step3);

TH2D *h_pnRes_theta_nmiss_Step3 = new TH2D("pnRes_theta_nmiss_Step3", "(P_{miss}-P_{n})/P_{miss} vs. #theta_{n,miss};(P_{miss}-P_{n})/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step3);

TH1D *h_ToF_goodN_Step3 = new TH1D("ToF_goodN_Step3", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_goodN_Step3);
TH1D *h_ToF_badN_Step3 = new TH1D("ToF_badN_Step3", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_badN_Step3);

TH2D *h_Edep_ToF_goodN_Step3 = new TH2D("Edep_ToF_goodN_Step3", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_goodN_Step3);
TH2D *h_Edep_ToF_badN_Step3 = new TH2D("Edep_ToF_badN_Step3", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_badN_Step3);

TH2D *h_sdiff_ldiff_allhit_goodN_Step3 = new TH2D("sdiff_ldiff_allhit_goodN_Step3", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
hist_list_2_A.push_back(h_sdiff_ldiff_allhit_goodN_Step3);
TH2D *h_sdiff_ldiff_allhit_badN_Step3 = new TH2D("sdiff_ldiff_allhit_badN_Step3", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
hist_list_2_A.push_back(h_sdiff_ldiff_allhit_badN_Step3);

TH2D *h_sdiff_ldiff_CTOFhit_goodN_Step3 = new TH2D("sdiff_ldiff_CTOFhit_goodN_Step3", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 3, 0.5, 3.5);
hist_list_2_A.push_back(h_sdiff_ldiff_CTOFhit_goodN_Step3);
TH2D *h_sdiff_ldiff_CTOFhit_badN_Step3 = new TH2D("sdiff_ldiff_CTOFhit_badN_Step3", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 3, 0.5, 3.5);
hist_list_2_A.push_back(h_sdiff_ldiff_CTOFhit_badN_Step3);

TH1D *h_numberCTOF_goodN_Step3 = new TH1D("numberCTOF_goodN_Step3", "Number of CTOF Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
hist_list_1_A.push_back(h_numberCTOF_goodN_Step3);
TH1D *h_numberCTOF_badN_Step3 = new TH1D("numberCTOF_badN_Step3", "Number of CTOF Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
hist_list_1_A.push_back(h_numberCTOF_badN_Step3);

TH2D *h_numberCTOF_momN_goodN_Step3 = new TH2D("numberCTOF_momN_goodN_Step3", "Number of CTOF Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
hist_list_2_A.push_back(h_numberCTOF_momN_goodN_Step3);
TH2D *h_numberCTOF_momN_badN_Step3 = new TH2D("numberCTOF_momN_badN_Step3", "Number of CTOF Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
hist_list_2_A.push_back(h_numberCTOF_momN_badN_Step3);

// Step Four (After applying Phi Diff CND hit cut) (Andrew)
// ======================================================================================================================================================================

/* Neutron histograms (from Erin) */
TH1D *h_n_multiplicity_allN_epCD_Step4 = new TH1D("n_multiplicity_allN_epCD_Step4", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_allN_epCD_Step4);
TH1D *h_n_multiplicity_goodN_epCD_Step4 = new TH1D("n_multiplicity_goodN_epCD_Step4", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_goodN_epCD_Step4);
TH1D *h_n_multiplicity_badN_epCD_Step4 = new TH1D("n_multiplicity_badN_epCD_Step4", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_badN_epCD_Step4);

/* Kinematical variables */
TH1D *h_theta_n_goodN_Step4 = new TH1D("theta_n_goodN_Step4", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_goodN_Step4);
TH1D *h_phi_n_goodN_Step4 = new TH1D("phi_n_goodN_Step4", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_goodN_Step4);
TH2D *h_theta_n_VS_phi_n_goodN_Step4 = new TH2D("theta_n_VS_phi_n_goodN_Step4", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_goodN_Step4);
TH2D *h_theta_n_VS_beta_n_goodN_Step4 = new TH2D("theta_VS_beta_goodN_Step4", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_goodN_Step4);
TH1D *h_theta_n_badN_Step4 = new TH1D("theta_n_badN_Step4", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_badN_Step4);
TH1D *h_phi_n_badN_Step4 = new TH1D("phi_n_badN_Step4", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_badN_Step4);
TH2D *h_theta_n_VS_phi_n_badN_Step4 = new TH2D("theta_n_VS_phi_n_badN_Step4", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_badN_Step4);
TH2D *h_theta_n_VS_beta_n_badN_Step4 = new TH2D("theta_VS_beta_badN_Step4", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_badN_Step4);

TH1D *h_P_n_goodN_Step4 = new TH1D("P_n_goodN_Step4", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_goodN_Step4);
TH2D *h_P_n_VS_theta_n_goodN_Step4 = new TH2D("P_n_VS_theta_n_goodN_Step4", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_goodN_Step4);
TH1D *h_P_n_badN_Step4 = new TH1D("P_n_badN_Step4", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_badN_Step4);
TH2D *h_P_n_VS_theta_n_badN_Step4 = new TH2D("P_n_VS_theta_n_badN_Step4", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_badN_Step4);

TH1D *h_P_miss_goodN_Step4 = new TH1D("P_miss_goodN_Step4", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_goodN_Step4);
TH2D *h_P_miss_VS_theta_miss_goodN_Step4 = new TH2D("P_miss_VS_theta_miss_goodN_Step4", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_goodN_Step4);
TH1D *h_P_miss_badN_Step4 = new TH1D("P_miss_badN_Step4", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_badN_Step4);
TH2D *h_P_miss_VS_theta_miss_badN_Step4 = new TH2D("P_miss_VS_theta_miss_badN_Step4", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_badN_Step4);

TH1D *h_P_n_minus_P_miss_goodN_Step4 = new TH1D("P_n_minus_P_miss_goodN_Step4", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_goodN_Step4);
TH1D *h_P_n_x_minus_P_miss_x_goodN_Step4 = new TH1D("P_n_x_minus_P_miss_x_goodN_Step4", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_goodN_Step4);
TH1D *h_P_n_y_minus_P_miss_y_goodN_Step4 = new TH1D("P_n_y_minus_P_miss_y_goodN_Step4", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_goodN_Step4);
TH1D *h_P_n_z_minus_P_miss_z_goodN_Step4 = new TH1D("P_n_z_minus_P_miss_z_goodN_Step4", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_goodN_Step4);
TH1D *h_P_n_minus_P_miss_badN_Step4 = new TH1D("P_n_minus_P_miss_badN_Step4", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_badN_Step4);
TH1D *h_P_n_x_minus_P_miss_x_badN_Step4 = new TH1D("P_n_x_minus_P_miss_x_badN_Step4", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_badN_Step4);
TH1D *h_P_n_y_minus_P_miss_y_badN_Step4 = new TH1D("P_n_y_minus_P_miss_y_badN_Step4", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_badN_Step4);
TH1D *h_P_n_z_minus_P_miss_z_badN_Step4 = new TH1D("P_n_z_minus_P_miss_z_badN_Step4", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_badN_Step4);

TH2D *h_P_n_VS_P_miss_goodN_Step4 = new TH2D("P_n_VS_P_miss_goodN_Step4", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_goodN_Step4);
TH2D *h_P_n_x_VS_P_miss_x_goodN_Step4 = new TH2D("P_n_x_VS_P_miss_x_goodN_Step4", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_goodN_Step4);
TH2D *h_P_n_y_VS_P_miss_y_goodN_Step4 = new TH2D("P_n_y_VS_P_miss_y_goodN_Step4", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_goodN_Step4);
TH2D *h_P_n_z_VS_P_miss_z_goodN_Step4 = new TH2D("P_n_z_VS_P_miss_z_goodN_Step4", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_goodN_Step4);
TH2D *h_P_n_VS_P_miss_badN_Step4 = new TH2D("P_n_VS_P_miss_badN_Step4", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_badN_Step4);
TH2D *h_P_n_x_VS_P_miss_x_badN_Step4 = new TH2D("P_n_x_VS_P_miss_x_badN_Step4", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_badN_Step4);
TH2D *h_P_n_y_VS_P_miss_y_badN_Step4 = new TH2D("P_n_y_VS_P_miss_y_badN_Step4", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_badN_Step4);
TH2D *h_P_n_z_VS_P_miss_z_badN_Step4 = new TH2D("P_n_z_VS_P_miss_z_badN_Step4", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_badN_Step4);

TH1D *h_E_p_CD_goodN_Step4 = new TH1D("E_p_CD_goodN_Step4", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_goodN_Step4);
TH1D *h_E_p_FD_goodN_Step4 = new TH1D("E_p_FD_goodN_Step4", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_goodN_Step4);
TH1D *h_E_miss_goodN_Step4 = new TH1D("E_miss_goodN_Step4", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_goodN_Step4);
TH1D *h_M_miss_goodN_Step4 = new TH1D("M_miss_goodN_Step4", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_goodN_Step4);
TH2D *h_M_miss_VS_P_n_goodN_Step4 = new TH2D("M_miss_VS_P_n_goodN_Step4", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_goodN_Step4);
TH2D *h_M_miss_VS_P_miss_goodN_Step4 = new TH2D("M_miss_P_miss_goodN_Step4", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_goodN_Step4);
TH1D *h_E_p_CD_badN_Step4 = new TH1D("E_p_CD_badN_Step4", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_badN_Step4);
TH1D *h_E_p_FD_badN_Step4 = new TH1D("E_p_FD_badN_Step4", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_badN_Step4);
TH1D *h_E_miss_badN_Step4 = new TH1D("E_miss_badN_Step4", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_badN_Step4);
TH1D *h_M_miss_badN_Step4 = new TH1D("M_miss_badN_Step4", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_badN_Step4);
TH2D *h_M_miss_VS_P_n_badN_Step4 = new TH2D("M_miss_VS_P_n_badN_Step4", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_badN_Step4);
TH2D *h_M_miss_VS_P_miss_badN_Step4 = new TH2D("M_miss_P_miss_badN_Step4", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_badN_Step4);

TH2D *h_theta_P_n_P_p_goodN_Step4 = new TH2D("theta_P_n_P_p_goodN_Step4", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_goodN_Step4);
TH2D *h_theta_P_n_P_p_badN_Step4 = new TH2D("theta_P_n_P_p_badN_Step4", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_badN_Step4);

TH1D *h_xB_goodN_Step4 = new TH1D("xB_goodN_Step4", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_goodN_Step4);
TH1D *h_xB_badN_Step4 = new TH1D("xB_badN_Step4", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_badN_Step4);

/* Detector responses */
TH1D *h_Edep_goodN_Step4 = new TH1D("n_Edep_goodN_Step4", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_goodN_Step4);
TH2D *h_P_n_VS_Edep_goodN_Step4 = new TH2D("P_n_VS_Edep_goodN_Step4", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_goodN_Step4);
TH2D *h_P_miss_VS_Edep_goodN_Step4 = new TH2D("P_miss_VS_Edep_goodN_Step4", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_goodN_Step4);
TH2D *h_dpp_VS_Edep_goodN_Step4 = new TH2D("dpp_VS_Edep_goodN_Step4", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_goodN_Step4);
TH1D *h_Edep_badN_Step4 = new TH1D("n_Edep_badN_Step4", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_badN_Step4);
TH2D *h_P_n_VS_Edep_badN_Step4 = new TH2D("P_n_VS_Edep_badN_Step4", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_badN_Step4);
TH2D *h_P_miss_VS_Edep_badN_Step4 = new TH2D("P_miss_VS_Edep_badN_Step4", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_badN_Step4);
TH2D *h_dpp_VS_Edep_badN_Step4 = new TH2D("dpp_VS_Edep_badN_Step4", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_badN_Step4);

TH2D *h_pnRes_theta_nmiss_Step4 = new TH2D("pnRes_theta_nmiss_Step4", "(P_{miss}-P_{n})/P_{miss} vs. #theta_{n,miss};(P_{miss}-P_{n})/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step4);

TH1D *h_ToF_goodN_Step4 = new TH1D("ToF_goodN_Step4", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_goodN_Step4);
TH1D *h_ToF_badN_Step4 = new TH1D("ToF_badN_Step4", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_badN_Step4);

TH2D *h_Edep_ToF_goodN_Step4 = new TH2D("Edep_ToF_goodN_Step4", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_goodN_Step4);
TH2D *h_Edep_ToF_badN_Step4 = new TH2D("Edep_ToF_badN_Step4", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_badN_Step4);

// Step Five (After event selection cuts) (Andrew)
// ======================================================================================================================================================================

/* Neutron histograms (from Erin) */
TH1D *h_n_multiplicity_allN_epCD_Step5 = new TH1D("n_multiplicity_allN_epCD_Step5", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_allN_epCD_Step5);
TH1D *h_n_multiplicity_goodN_epCD_Step5 = new TH1D("n_multiplicity_goodN_epCD_Step5", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_goodN_epCD_Step5);
TH1D *h_n_multiplicity_badN_epCD_Step5 = new TH1D("n_multiplicity_badN_epCD_Step5", "Number of Neutrons in Event", 10, 0, 10);
hist_list_1.push_back(h_n_multiplicity_badN_epCD_Step5);

/* Kinematical variables */
TH1D *h_theta_n_goodN_Step5 = new TH1D("theta_n_goodN_Step5", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_goodN_Step5);
TH1D *h_phi_n_goodN_Step5 = new TH1D("phi_n_goodN_Step5", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_goodN_Step5);
TH2D *h_theta_n_VS_phi_n_goodN_Step5 = new TH2D("theta_n_VS_phi_n_goodN_Step5", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_goodN_Step5);
TH2D *h_theta_n_VS_beta_n_goodN_Step5 = new TH2D("theta_VS_beta_goodN_Step5", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_goodN_Step5);
TH1D *h_theta_n_badN_Step5 = new TH1D("theta_n_badN_Step5", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
hist_list_1_A.push_back(h_theta_n_badN_Step5);
TH1D *h_phi_n_badN_Step5 = new TH1D("phi_n_badN_Step5", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 50, -180, 180);
hist_list_1_A.push_back(h_phi_n_badN_Step5);
TH2D *h_theta_n_VS_phi_n_badN_Step5 = new TH2D("theta_n_VS_phi_n_badN_Step5", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 50, -180, 180, 50, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_phi_n_badN_Step5);
TH2D *h_theta_n_VS_beta_n_badN_Step5 = new TH2D("theta_VS_beta_badN_Step5", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
hist_list_2_A.push_back(h_theta_n_VS_beta_n_badN_Step5);

TH1D *h_P_n_goodN_Step5 = new TH1D("P_n_goodN_Step5", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_goodN_Step5);
TH2D *h_P_n_VS_theta_n_goodN_Step5 = new TH2D("P_n_VS_theta_n_goodN_Step5", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_goodN_Step5);
TH1D *h_P_n_badN_Step5 = new TH1D("P_n_badN_Step5", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_badN_Step5);
TH2D *h_P_n_VS_theta_n_badN_Step5 = new TH2D("P_n_VS_theta_n_badN_Step5", "Neutron Momentum vs Theta;#theta [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
hist_list_2_A.push_back(h_P_n_VS_theta_n_badN_Step5);

TH1D *h_P_miss_goodN_Step5 = new TH1D("P_miss_goodN_Step5", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_goodN_Step5);
TH2D *h_P_miss_VS_theta_miss_goodN_Step5 = new TH2D("P_miss_VS_theta_miss_goodN_Step5", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_goodN_Step5);
TH1D *h_P_miss_badN_Step5 = new TH1D("P_miss_badN_Step5", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
hist_list_1_A.push_back(h_P_miss_badN_Step5);
TH2D *h_P_miss_VS_theta_miss_badN_Step5 = new TH2D("P_miss_VS_theta_miss_badN_Step5", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_theta_miss_badN_Step5);

TH1D *h_P_n_minus_P_miss_goodN_Step5 = new TH1D("P_n_minus_P_miss_goodN_Step5", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_goodN_Step5);
TH1D *h_P_n_x_minus_P_miss_x_goodN_Step5 = new TH1D("P_n_x_minus_P_miss_x_goodN_Step5", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_goodN_Step5);
TH1D *h_P_n_y_minus_P_miss_y_goodN_Step5 = new TH1D("P_n_y_minus_P_miss_y_goodN_Step5", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_goodN_Step5);
TH1D *h_P_n_z_minus_P_miss_z_goodN_Step5 = new TH1D("P_n_z_minus_P_miss_z_goodN_Step5", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_goodN_Step5);
TH1D *h_P_n_minus_P_miss_badN_Step5 = new TH1D("P_n_minus_P_miss_badN_Step5", "P_{n}-P_{miss} [GeV/c];P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_minus_P_miss_badN_Step5);
TH1D *h_P_n_x_minus_P_miss_x_badN_Step5 = new TH1D("P_n_x_minus_P_miss_x_badN_Step5", "P_{x,n}-P_{x,miss} Distribution; P_{x,n}-P_{x,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_x_minus_P_miss_x_badN_Step5);
TH1D *h_P_n_y_minus_P_miss_y_badN_Step5 = new TH1D("P_n_y_minus_P_miss_y_badN_Step5", "P_{y,n}-P_{y,miss} Distribution; P_{y,n}-P_{y,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_y_minus_P_miss_y_badN_Step5);
TH1D *h_P_n_z_minus_P_miss_z_badN_Step5 = new TH1D("P_n_z_minus_P_miss_z_badN_Step5", "P_{z,n}-P_{z,miss} Distribution; P_{z,n}-P_{z,miss} [GeV/c];Counts", 50, -1.5, 1.5);
hist_list_1_A.push_back(h_P_n_z_minus_P_miss_z_badN_Step5);

TH2D *h_P_n_VS_P_miss_goodN_Step5 = new TH2D("P_n_VS_P_miss_goodN_Step5", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_goodN_Step5);
TH2D *h_P_n_x_VS_P_miss_x_goodN_Step5 = new TH2D("P_n_x_VS_P_miss_x_goodN_Step5", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_goodN_Step5);
TH2D *h_P_n_y_VS_P_miss_y_goodN_Step5 = new TH2D("P_n_y_VS_P_miss_y_goodN_Step5", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_goodN_Step5);
TH2D *h_P_n_z_VS_P_miss_z_goodN_Step5 = new TH2D("P_n_z_VS_P_miss_z_goodN_Step5", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_goodN_Step5);
TH2D *h_P_n_VS_P_miss_badN_Step5 = new TH2D("P_n_VS_P_miss_badN_Step5", "P_{n} vs P_{miss} [GeV/c];P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_VS_P_miss_badN_Step5);
TH2D *h_P_n_x_VS_P_miss_x_badN_Step5 = new TH2D("P_n_x_VS_P_miss_x_badN_Step5", "P_{x,n} vs P_{x,miss} Distribution;P_{x,n} [GeV/c];P_{x,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_x_VS_P_miss_x_badN_Step5);
TH2D *h_P_n_y_VS_P_miss_y_badN_Step5 = new TH2D("P_n_y_VS_P_miss_y_badN_Step5", "P_{y,n} vs P_{y,miss} Distribution;P_{y,n} [GeV/c];P_{y,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_y_VS_P_miss_y_badN_Step5);
TH2D *h_P_n_z_VS_P_miss_z_badN_Step5 = new TH2D("P_n_z_VS_P_miss_z_badN_Step5", "P_{z,n} vs P_{z,miss} Distribution;P_{z,n} [GeV/c];P_{z,miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
hist_list_1_A.push_back(h_P_n_z_VS_P_miss_z_badN_Step5);

TH1D *h_E_p_CD_goodN_Step5 = new TH1D("E_p_CD_goodN_Step5", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_goodN_Step5);
TH1D *h_E_p_FD_goodN_Step5 = new TH1D("E_p_FD_goodN_Step5", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_goodN_Step5);
TH1D *h_E_miss_goodN_Step5 = new TH1D("E_miss_goodN_Step5", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_goodN_Step5);
TH1D *h_M_miss_goodN_Step5 = new TH1D("M_miss_goodN_Step5", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_goodN_Step5);
TH2D *h_M_miss_VS_P_n_goodN_Step5 = new TH2D("M_miss_VS_P_n_goodN_Step5", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_goodN_Step5);
TH2D *h_M_miss_VS_P_miss_goodN_Step5 = new TH2D("M_miss_P_miss_goodN_Step5", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_goodN_Step5);
TH1D *h_E_p_CD_badN_Step5 = new TH1D("E_p_CD_badN_Step5", "CD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_CD_badN_Step5);
TH1D *h_E_p_FD_badN_Step5 = new TH1D("E_p_FD_badN_Step5", "FD Proton Energy;E_{p} [GeV]", 50, 0, 1.5);
hist_list_1_A.push_back(h_E_p_FD_badN_Step5);
TH1D *h_E_miss_badN_Step5 = new TH1D("E_miss_badN_Step5", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_E_miss_badN_Step5);
TH1D *h_M_miss_badN_Step5 = new TH1D("M_miss_badN_Step5", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
hist_list_1_A.push_back(h_M_miss_badN_Step5);
TH2D *h_M_miss_VS_P_n_badN_Step5 = new TH2D("M_miss_VS_P_n_badN_Step5", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_n_badN_Step5);
TH2D *h_M_miss_VS_P_miss_badN_Step5 = new TH2D("M_miss_P_miss_badN_Step5", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
hist_list_2_A.push_back(h_M_miss_VS_P_miss_badN_Step5);

TH2D *h_theta_P_n_P_p_goodN_Step5 = new TH2D("theta_P_n_P_p_goodN_Step5", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_goodN_Step5);
TH2D *h_theta_P_n_P_p_badN_Step5 = new TH2D("theta_P_n_P_p_badN_Step5", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
hist_list_2_A.push_back(h_theta_P_n_P_p_badN_Step5);

TH1D *h_xB_goodN_Step5 = new TH1D("xB_goodN_Step5", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_goodN_Step5);
TH1D *h_xB_badN_Step5 = new TH1D("xB_badN_Step5", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
hist_list_1_A.push_back(h_xB_badN_Step5);

/* Detector responses */
TH1D *h_Edep_goodN_Step5 = new TH1D("n_Edep_goodN_Step5", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_goodN_Step5);
TH2D *h_P_n_VS_Edep_goodN_Step5 = new TH2D("P_n_VS_Edep_goodN_Step5", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_goodN_Step5);
TH2D *h_P_miss_VS_Edep_goodN_Step5 = new TH2D("P_miss_VS_Edep_goodN_Step5", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_goodN_Step5);
TH2D *h_dpp_VS_Edep_goodN_Step5 = new TH2D("dpp_VS_Edep_goodN_Step5", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_goodN_Step5);
TH1D *h_Edep_badN_Step5 = new TH1D("n_Edep_badN_Step5", "Neutron Energy Deposition;Energy [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_badN_Step5);
TH2D *h_P_n_VS_Edep_badN_Step5 = new TH2D("P_n_VS_Edep_badN_Step5", "Neutron Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_n_VS_Edep_badN_Step5);
TH2D *h_P_miss_VS_Edep_badN_Step5 = new TH2D("P_miss_VS_Edep_badN_Step5", "Missing Momentum vs Neutron Energy Deposition;E_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
hist_list_2_A.push_back(h_P_miss_VS_Edep_badN_Step5);
TH2D *h_dpp_VS_Edep_badN_Step5 = new TH2D("dpp_VS_Edep_badN_Step5", "(P_{miss}-P_{n})/P_{miss} vs Neutron Energy Deposition;E_{dep} [MeV];(P_{miss}-P_{n})/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
hist_list_2_A.push_back(h_dpp_VS_Edep_badN_Step5);

TH2D *h_pnRes_theta_nmiss_Step5 = new TH2D("pnRes_theta_nmiss_Step5", "(P_{miss}-P_{n})/P_{miss} vs. #theta_{n,miss};(P_{miss}-P_{n})/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step5);

TH1D *h_ToF_goodN_Step5 = new TH1D("ToF_goodN_Step5", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_goodN_Step5);
TH1D *h_ToF_badN_Step5 = new TH1D("ToF_badN_Step5", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
hist_list_1_A.push_back(h_ToF_badN_Step5);

TH2D *h_Edep_ToF_goodN_Step5 = new TH2D("Edep_ToF_goodN_Step5", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_goodN_Step5);
TH2D *h_Edep_ToF_badN_Step5 = new TH2D("Edep_ToF_badN_Step5", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_Edep_ToF_badN_Step5);

TH1D *h_pmiss_goodN_Step5 = new TH1D("pmiss_goodN_Step5", "P_{miss} good N Step5;P_{miss} [GeV/c];Counts", 25, 0.25, 1.0);
hist_list_1_A.push_back(h_pmiss_goodN_Step5);
TH1D *h_pmiss_allN_Step5 = new TH1D("pmiss_allN_Step5", "P_{miss} all N Step5;P_{miss} [GeV/c];Counts", 25, 0.25, 1.0);
hist_list_1_A.push_back(h_pmiss_allN_Step5);

TH1D *h_Edep_infront_goodN_Step5 = new TH1D("Edep_infront_goodN_Step5", "E_{dep} of Hit infront CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_infront_goodN_Step5);
TH1D *h_Edep_infront_badN_Step5 = new TH1D("Edep_infront_badN_Step5", "E_{dep} of Hit infront CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_infront_badN_Step5);

TH1D *h_Edep_behind_goodN_Step5 = new TH1D("Edep_behind_goodN_Step5", "E_{dep} of Hit behind CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_behind_goodN_Step5);
TH1D *h_Edep_behind_badN_Step5 = new TH1D("Edep_behind_badN_Step5", "E_{dep} of Hit behind CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
hist_list_1_A.push_back(h_Edep_behind_badN_Step5);

TH2D *h_diff_ToFc_z_Edep_goodN_Step5 = new TH2D("diff_ToFc_z_Edep_goodN_Step5", "ToF*c - z vs. E_{dep} of CND Neutrons;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
hist_list_2_A.push_back(h_diff_ToFc_z_Edep_goodN_Step5);
TH2D *h_diff_ToFc_z_Edep_badN_Step5 = new TH2D("diff_ToFc_z_Edep_badN_Step5", "ToF*c - z vs. E_{dep} of CND Neutrons;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
hist_list_2_A.push_back(h_diff_ToFc_z_Edep_badN_Step5);

TH2D *h_diff_ToFc_z_Edep_goodN_Step5_layer[3];

for (int k = 0; k < 3; k++)
{
    sprintf(temp_name_A, "diff_ToFc_z_goodN_Step5_layer_%d", k + 1);
    sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
    h_diff_ToFc_z_Edep_goodN_Step5_layer[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_goodN_Step5_layer[k]);
}

TH2D *h_diff_ToFc_z_Edep_badN_Step5_layer[3];

for (int k = 0; k < 3; k++)
{
    sprintf(temp_name_A, "diff_ToFc_z_badN_Step5_layer_%d", k + 1);
    sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
    h_diff_ToFc_z_Edep_badN_Step5_layer[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_badN_Step5_layer[k]);
}

TH1D *h_phidiff_en_goodN_Step5 = new TH1D("phidiff_en_goodN_Step5", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
hist_list_1_A.push_back(h_phidiff_en_goodN_Step5);
TH1D *h_phidiff_en_badN_Step5 = new TH1D("phidiff_en_badN_Step5", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
hist_list_1_A.push_back(h_phidiff_en_badN_Step5);

TH1D *h_TP_goodN_Step5 = new TH1D("TP_goodN_Step5", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
hist_list_1_A.push_back(h_TP_goodN_Step5);
TH1D *h_TP_badN_Step5 = new TH1D("TP_badN_Step5", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
hist_list_1_A.push_back(h_TP_badN_Step5);

TH1D *h_Z_goodN_Step5 = new TH1D("Z_goodN_Step5", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
hist_list_1_A.push_back(h_Z_goodN_Step5);
TH1D *h_Z_badN_Step5 = new TH1D("Z_badN_Step5", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
hist_list_1_A.push_back(h_Z_badN_Step5);

TH2D *h_beta_Edep_goodN_Step5 = new TH2D("Edep_beta_goodN_Step5", "#beta vs. E_{dep} of CND Neutrons;#beta;E_{dep} [MeV]", 50, 0, 1.1, 50, 0, 100);
hist_list_2_A.push_back(h_beta_Edep_goodN_Step5);
TH2D *h_beta_Edep_badN_Step5 = new TH2D("Edep_beta_badN_Step5", "#beta vs. E_{dep} of CND Neutrons;#beta;E_{dep} [MeV]", 50, 0, 1.1, 50, 0, 100);
hist_list_2_A.push_back(h_beta_Edep_badN_Step5);

TH2D *h_ToF_Edep_goodN_Step5 = new TH2D("ToF_Edep_goodN_Step5", "ToF vs. E_{dep} of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_ToF_Edep_goodN_Step5);
TH2D *h_ToF_Edep_badN_Step5 = new TH2D("ToF_Edep_badN_Step5", "ToF vs. E_{dep} of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
hist_list_2_A.push_back(h_ToF_Edep_badN_Step5);

TH2D *h_TP_Edep_goodN_Step5 = new TH2D("TP_Edep_goodN_Step5", "TP vs. E_{dep} of CND Neutrons;TP [ns/m];E_{dep} [MeV]", 150, 0, 50, 50, 0, 100);
hist_list_2_A.push_back(h_TP_Edep_goodN_Step5);
TH2D *h_TP_Edep_badN_Step5 = new TH2D("TP_Edep_badN_Step5", "TP vs. E_{dep} of CND Neutrons;TP [ns/m];E_{dep} [MeV]", 150, 0, 50, 50, 0, 100);
hist_list_2_A.push_back(h_TP_Edep_badN_Step5);

for (int i = 0; i < hist_list_1_A.size(); i++)
{
    hist_list_1_A[i]->Sumw2();
    hist_list_1_A[i]->GetXaxis()->CenterTitle();
    hist_list_1_A[i]->GetYaxis()->CenterTitle();
}

for (int i = 0; i < hist_list_2_A.size(); i++)
{
    hist_list_2_A[i]->Sumw2();
    hist_list_2_A[i]->GetXaxis()->CenterTitle();
    hist_list_2_A[i]->GetYaxis()->CenterTitle();
}

#pragma endregion /* Andrew's histograms - end */

#endif // MANUALVETOHISTOGRAMS_H
