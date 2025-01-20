//
// Created by Alon Sportes on 19/01/2025.
//

#ifndef VETOHISTOGRAMS_H
#define VETOHISTOGRAMS_H

#include <cstdlib>
#include <iostream>
#include <vector>

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

#include "../../functions/GeneralFunctions.h"
#include "../../constants.h"

using namespace std;

// ======================================================================================================================================================================
// Manual neutron veto histograms
// ======================================================================================================================================================================

class VetoHistograms {
public:
#pragma region /* Veto histograms - start */

    /////////////////////////////////////
    // Prepare histograms
    /////////////////////////////////////

    vector<TH1 *> HistoList;

    char temp_name[300];
    char temp_title[300];

    // (e,e'p) plots
    // ======================================================================================================================================================================

#pragma region /* (e,e'p) plots - start */

    /* Proton histograms (from Erin) */
    TH1D *h_p_multiplicity_BPID_epCD;
    TH1D *h_P_p_BPID_epCD;
    TH2D *h_theta_p_VS_phi_p_BPID_epCD;
    TH1D *h_p_multiplicity_APID_epCD;
    TH1D *h_P_p_APID_epCD;
    TH2D *h_theta_p_VS_phi_p_APID_epCD;

    TH1D *h_p_multiplicity_BPID_epFD;
    TH1D *h_P_p_BPID_epFD;
    TH2D *h_theta_p_VS_phi_p_BPID_epFD;
    TH1D *h_p_multiplicity_APID_epFD;
    TH1D *h_P_p_APID_epFD;
    TH2D *h_theta_p_VS_phi_p_APID_epFD;

    TH2D *h_dbeta_p_VS_P_p_BPID_epCD;
    TH1D *h_dVz_p_BPID_epCD;
    TH1D *h_Chi2pid_p_BPID_epCD;
    TH2D *h_dbeta_p_VS_P_p_APID_epCD;
    TH1D *h_dVz_p_APID_epCD;
    TH1D *h_Chi2pid_p_APID_epCD;

    TH2D *h_dbeta_p_VS_P_p_BPID_epFD;
    TH1D *h_dVz_p_BPID_epFD;
    TH1D *h_Chi2pid_p_BPID_epFD;
    TH2D *h_dbeta_p_VS_P_p_APID_epFD;
    TH1D *h_dVz_p_APID_epFD;
    TH1D *h_Chi2pid_p_APID_epFD;

    /* Missing variabels */
    TH1D *h_P_miss_BmissC_epCD;
    TH1D *h_theta_miss_BmissC_epCD;
    TH2D *h_P_miss_VS_theta_miss_BmissC_epCD;
    TH1D *h_P_miss_AmissC_epCD;
    TH1D *h_theta_miss_AmissC_epCD;
    TH2D *h_P_miss_VS_theta_miss_AmissC_epCD;

    TH1D *h_P_miss_BmissC_epFD;
    TH1D *h_theta_miss_BmissC_epFD;
    TH2D *h_P_miss_VS_theta_miss_BmissC_epFD;
    TH1D *h_P_miss_AmissC_epFD;
    TH1D *h_theta_miss_AmissC_epFD;
    TH2D *h_P_miss_VS_theta_miss_AmissC_epFD;

    TH1D *h_E_p_BmissC_epCD;
    TH1D *h_E_miss_BmissC_epCD;
    TH1D *h_M_miss_BmissC_epCD;
    TH1D *h_E_p_AmissC_epCD;
    TH1D *h_E_miss_AmissC_epCD;
    TH1D *h_M_miss_AmissC_epCD;

    TH1D *h_E_p_BmissC_epFD;
    TH1D *h_E_miss_BmissC_epFD;
    TH1D *h_M_miss_BmissC_epFD;
    TH1D *h_E_p_AmissC_epFD;
    TH1D *h_E_miss_AmissC_epFD;
    TH1D *h_M_miss_AmissC_epFD;

    /* Checks on which events have neutrons (Andrew) */
    TH1D *h_xB_BmissC_epCD;
    TH2D *h_xB_VS_M_miss_BmissC_epCD;
    TH1D *h_xB_AmissC_epCD;
    TH2D *h_xB_VS_M_miss_AmissC_epCD;

    TH1D *h_xB_BmissC_epFD;
    TH2D *h_xB_VS_M_miss_BmissC_epFD;
    TH1D *h_xB_AmissC_epFD;
    TH2D *h_xB_VS_M_miss_AmissC_epFD;

    TH2D *h_xB_VS_M_miss_epCDn;
    TH2D *h_xB_VS_M_miss_epFDn;

    TH2D *h_xB_VS_M_miss_goodN_epCDn;
    TH2D *h_xB_VS_M_miss_badN_epCDn;

    TH2D *h_xB_VS_M_miss_goodN_epFDn;
    TH2D *h_xB_VS_M_miss_badN_epFDn;

    /* Kinematical variables */
    TH1D *h_theta_n_epCDn;
    TH1D *h_phi_n_epCDn;
    TH2D *h_theta_n_VS_phi_n_epCDn;
    TH2D *h_theta_n_VS_beta_n_epCDn;

    TH1D *h_theta_n_epFDn;
    TH1D *h_phi_n_epFDn;
    TH2D *h_theta_n_VS_phi_n_epFDn;
    TH2D *h_theta_n_VS_beta_n_epFDn;

    TH1D *h_P_n_epCDn;
    TH2D *h_P_n_VS_theta_n_epCDn;

    TH1D *h_P_n_epFDn;
    TH2D *h_P_n_VS_theta_n_epFDn;

    TH1D *h_P_miss_epCDn;
    TH2D *h_P_miss_VS_theta_miss_epCDn;

    TH1D *h_P_miss_epFDn;
    TH2D *h_P_miss_VS_theta_miss_epFDn;

    TH1D *h_dpp_allN_epCDn;
    TH1D *h_dpp_goodN_epCDn;
    TH1D *h_dpp_badN_epCDn;
    TH1D *h_dpp_allN_for_theta_n_miss_less_than_25_epCDn;

    TH1D *h_dpp_allN_epFDn;
    TH1D *h_dpp_goodN_epFDn;
    TH1D *h_dpp_badN_epFDn;
    TH1D *h_dpp_allN_for_theta_n_miss_less_than_25_epFDn;

    TH1D *h_theta_n_miss_allN_epCDn;
    TH1D *h_theta_n_miss_goodN_epCDn;
    TH1D *h_theta_n_miss_badN_epCDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_05_epCDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_03_epCDn;

    TH1D *h_theta_n_miss_allN_epFDn;
    TH1D *h_theta_n_miss_goodN_epFDn;
    TH1D *h_theta_n_miss_badN_epFDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_05_epFDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_03_epFDn;

    TH2D *h_dpp_VS_theta_n_miss_epCDn;

    TH2D *h_dpp_VS_theta_n_miss_epFDn;

    TH1D *h_E_p_epCDn;
    TH1D *h_E_miss_epCDn;
    TH1D *h_M_miss_epCDn;
    TH2D *h_M_miss_VS_P_n_epCDn;
    TH2D *h_M_miss_VS_theta_n_epCDn;
    TH2D *h_M_miss_VS_phi_n_epCDn;
    TH2D *h_M_miss_VS_P_miss_epCDn;
    TH2D *h_M_miss_VS_theta_miss_epCDn;
    TH2D *h_M_miss_VS_phi_miss_epCDn;

    TH1D *h_E_p_epFDn;
    TH1D *h_E_miss_epFDn;
    TH1D *h_M_miss_epFDn;
    TH2D *h_M_miss_VS_P_n_epFDn;
    TH2D *h_M_miss_VS_theta_n_epFDn;
    TH2D *h_M_miss_VS_phi_n_epFDn;
    TH2D *h_M_miss_VS_P_miss_epFDn;
    TH2D *h_M_miss_VS_theta_miss_epFDn;
    TH2D *h_M_miss_VS_phi_miss_epFDn;

    TH1D *h_P_n_minus_P_miss_epCDn;
    TH1D *h_P_n_x_minus_P_miss_x_epCDn;
    TH1D *h_P_n_y_minus_P_miss_y_epCDn;
    TH1D *h_P_n_z_minus_P_miss_z_epCDn;

    TH1D *h_P_n_minus_P_miss_epFDn;
    TH1D *h_P_n_x_minus_P_miss_x_epFDn;
    TH1D *h_P_n_y_minus_P_miss_y_epFDn;
    TH1D *h_P_n_z_minus_P_miss_z_epFDn;

    TH2D *h_P_n_VS_P_miss_epCDn;
    TH2D *h_P_n_x_VS_P_miss_x_epCDn;
    TH2D *h_P_n_y_VS_P_miss_y_epCDn;
    TH2D *h_P_n_z_VS_P_miss_z_epCDn;

    TH2D *h_P_n_VS_P_miss_epFDn;
    TH2D *h_P_n_x_VS_P_miss_x_epFDn;
    TH2D *h_P_n_y_VS_P_miss_y_epFDn;
    TH2D *h_P_n_z_VS_P_miss_z_epFDn;

    TH1D *h_theta_n_p_epCDn;
    TH2D *h_theta_p_n_VS_P_p_epCDn;

    TH1D *h_theta_n_p_epFDn;
    TH2D *h_theta_p_n_VS_P_p_epFDn;

    TH1D *h_xB_epCDn;

    TH1D *h_xB_epFDn;

    /* Detector responses */
    TH1D *h_Edep_CND_epCDn;
    TH2D *h_P_n_VS_Edep_CND_epCDn;
    TH2D *h_theta_n_VS_Edep_CND_epCDn;
    TH2D *h_phi_n_VS_Edep_CND_epCDn;
    TH2D *h_P_miss_VS_Edep_CND_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND_epCDn;
    TH2D *h_dpp_VS_Edep_CND_epCDn;
    TH2D *h_beta_n_VS_Edep_CND_epCDn;
    TH2D *h_E_miss_VS_Edep_CND_epCDn;
    TH2D *h_M_miss_VS_Edep_CND_epCDn;
    TH2D *h_path_VS_Edep_CND_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_epCDn;
    TH2D *h_ToF_VS_Edep_CND_epCDn;
    TH2D *h_nSector_VS_Edep_CND_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_epCDn;

    TH1D *h_Edep_CND_epFDn;
    TH2D *h_P_n_VS_Edep_CND_epFDn;
    TH2D *h_theta_n_VS_Edep_CND_epFDn;
    TH2D *h_phi_n_VS_Edep_CND_epFDn;
    TH2D *h_P_miss_VS_Edep_CND_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND_epFDn;
    TH2D *h_dpp_VS_Edep_CND_epFDn;
    TH2D *h_beta_n_VS_Edep_CND_epFDn;
    TH2D *h_E_miss_VS_Edep_CND_epFDn;
    TH2D *h_M_miss_VS_Edep_CND_epFDn;
    TH2D *h_path_VS_Edep_CND_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_epFDn;
    TH2D *h_ToF_VS_Edep_CND_epFDn;
    TH2D *h_nSector_VS_Edep_CND_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_epFDn;

    TH1D *h_Edep_CTOF_epCDn;
    TH2D *h_P_n_VS_Edep_CTOF_epCDn;
    TH2D *h_theta_n_VS_Edep_CTOF_epCDn;
    TH2D *h_phi_n_VS_Edep_CTOF_epCDn;
    TH2D *h_P_miss_VS_Edep_CTOF_epCDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_epCDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_epCDn;
    TH2D *h_dpp_VS_Edep_CTOF_epCDn;
    TH2D *h_beta_n_VS_Edep_CTOF_epCDn;
    TH2D *h_E_miss_VS_Edep_CTOF_epCDn;
    TH2D *h_M_miss_VS_Edep_CTOF_epCDn;
    TH2D *h_path_VS_Edep_CTOF_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_epCDn;
    TH2D *h_ToF_VS_Edep_CTOF_epCDn;
    TH2D *h_nSector_VS_Edep_CTOF_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_epCDn;

    TH1D *h_Edep_CTOF_epFDn;
    TH2D *h_P_n_VS_Edep_CTOF_epFDn;
    TH2D *h_theta_n_VS_Edep_CTOF_epFDn;
    TH2D *h_phi_n_VS_Edep_CTOF_epFDn;
    TH2D *h_P_miss_VS_Edep_CTOF_epFDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_epFDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_epFDn;
    TH2D *h_dpp_VS_Edep_CTOF_epFDn;
    TH2D *h_beta_n_VS_Edep_CTOF_epFDn;
    TH2D *h_E_miss_VS_Edep_CTOF_epFDn;
    TH2D *h_M_miss_VS_Edep_CTOF_epFDn;
    TH2D *h_path_VS_Edep_CTOF_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_epFDn;
    TH2D *h_ToF_VS_Edep_CTOF_epFDn;
    TH2D *h_nSector_VS_Edep_CTOF_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_epFDn;

    // TH1D *h_Edep_single_epCDn;
    // TH2D *h_P_n_VS_Edep_single_epCDn;
    // TH2D *h_theta_n_VS_Edep_single_epCDn;
    // TH2D *h_phi_n_VS_Edep_single_epCDn;
    // TH2D *h_P_miss_VS_Edep_single_epCDn;
    // TH2D *h_theta_miss_VS_Edep_single_epCDn;
    // TH2D *h_phi_miss_VS_Edep_single_epCDn;
    // TH2D *h_dpp_VS_Edep_single_epCDn;
    // TH2D *h_beta_n_VS_Edep_single_epCDn;
    // TH2D *h_E_miss_VS_Edep_single_epCDn;
    // TH2D *h_M_miss_VS_Edep_single_epCDn;
    // TH2D *h_path_VS_Edep_single_epCDn;
    // TH2D *h_theta_n_miss_VS_Edep_single_epCDn;
    // TH2D *h_ToF_VS_Edep_single_epCDn;
    // TH2D *h_nSector_VS_Edep_single_epCDn;

    // TH1D *h_Edep_single_epFDn;
    // TH2D *h_P_n_VS_Edep_single_epFDn;
    // TH2D *h_theta_n_VS_Edep_single_epFDn;
    // TH2D *h_phi_n_VS_Edep_single_epFDn;
    // TH2D *h_P_miss_VS_Edep_single_epFDn;
    // TH2D *h_theta_miss_VS_Edep_single_epFDn;
    // TH2D *h_phi_miss_VS_Edep_single_epFDn;
    // TH2D *h_dpp_VS_Edep_single_epFDn;
    // TH2D *h_beta_n_VS_Edep_single_epFDn;
    // TH2D *h_E_miss_VS_Edep_single_epFDn;
    // TH2D *h_M_miss_VS_Edep_single_epFDn;
    // TH2D *h_path_VS_Edep_single_epFDn;
    // TH2D *h_theta_n_miss_VS_Edep_single_epFDn;
    // TH2D *h_ToF_VS_Edep_single_epFDn;
    // TH2D *h_nSector_VS_Edep_single_epFDn;

    TH1D *h_Edep_CND1_epCDn;
    TH2D *h_P_n_VS_Edep_CND1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND1_epCDn;
    TH2D *h_dpp_VS_Edep_CND1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND1_epCDn;
    TH2D *h_path_VS_Edep_CND1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_epCDn;
    TH2D *h_ToF_VS_Edep_CND1_epCDn;
    TH2D *h_nSector_VS_Edep_CND1_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_epCDn;

    TH1D *h_Edep_CND1_epFDn;
    TH2D *h_P_n_VS_Edep_CND1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND1_epFDn;
    TH2D *h_dpp_VS_Edep_CND1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND1_epFDn;
    TH2D *h_path_VS_Edep_CND1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_epFDn;
    TH2D *h_ToF_VS_Edep_CND1_epFDn;
    TH2D *h_nSector_VS_Edep_CND1_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_epFDn;

    TH1D *h_Edep_CND2_epCDn;
    TH2D *h_P_n_VS_Edep_CND2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND2_epCDn;
    TH2D *h_dpp_VS_Edep_CND2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND2_epCDn;
    TH2D *h_path_VS_Edep_CND2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_epCDn;
    TH2D *h_ToF_VS_Edep_CND2_epCDn;
    TH2D *h_nSector_VS_Edep_CND2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_epCDn;

    TH1D *h_Edep_CND2_epFDn;
    TH2D *h_P_n_VS_Edep_CND2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND2_epFDn;
    TH2D *h_dpp_VS_Edep_CND2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND2_epFDn;
    TH2D *h_path_VS_Edep_CND2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_epFDn;
    TH2D *h_ToF_VS_Edep_CND2_epFDn;
    TH2D *h_nSector_VS_Edep_CND2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_epFDn;

    TH1D *h_Edep_CND3_epCDn;
    TH2D *h_P_n_VS_Edep_CND3_epCDn;
    TH2D *h_theta_n_VS_Edep_CND3_epCDn;
    TH2D *h_phi_n_VS_Edep_CND3_epCDn;
    TH2D *h_P_miss_VS_Edep_CND3_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND3_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND3_epCDn;
    TH2D *h_dpp_VS_Edep_CND3_epCDn;
    TH2D *h_beta_n_VS_Edep_CND3_epCDn;
    TH2D *h_E_miss_VS_Edep_CND3_epCDn;
    TH2D *h_M_miss_VS_Edep_CND3_epCDn;
    TH2D *h_path_VS_Edep_CND3_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_epCDn;
    TH2D *h_ToF_VS_Edep_CND3_epCDn;
    TH2D *h_nSector_VS_Edep_CND3_epCDn;

    TH1D *h_Edep_CND3_epFDn;
    TH2D *h_P_n_VS_Edep_CND3_epFDn;
    TH2D *h_theta_n_VS_Edep_CND3_epFDn;
    TH2D *h_phi_n_VS_Edep_CND3_epFDn;
    TH2D *h_P_miss_VS_Edep_CND3_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND3_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND3_epFDn;
    TH2D *h_dpp_VS_Edep_CND3_epFDn;
    TH2D *h_beta_n_VS_Edep_CND3_epFDn;
    TH2D *h_E_miss_VS_Edep_CND3_epFDn;
    TH2D *h_M_miss_VS_Edep_CND3_epFDn;
    TH2D *h_path_VS_Edep_CND3_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_epFDn;
    TH2D *h_ToF_VS_Edep_CND3_epFDn;
    TH2D *h_nSector_VS_Edep_CND3_epFDn;

    TH1D *h_Size_CND1_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND1_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_epCDn;
    TH2D *h_P_n_VS_Size_CND1_epCDn;
    TH2D *h_theta_n_VS_Size_CND1_epCDn;
    TH2D *h_phi_n_VS_Size_CND1_epCDn;
    TH2D *h_P_miss_VS_Size_CND1_epCDn;
    TH2D *h_theta_miss_VS_Size_CND1_epCDn;
    TH2D *h_phi_miss_VS_Size_CND1_epCDn;
    TH2D *h_dpp_VS_Size_CND1_epCDn;
    TH2D *h_beta_n_VS_Size_CND1_epCDn;
    TH2D *h_E_miss_VS_Size_CND1_epCDn;
    TH2D *h_M_miss_VS_Size_CND1_epCDn;
    TH2D *h_path_VS_Size_CND1_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_epCDn;
    TH2D *h_ToF_VS_Size_CND1_epCDn;
    TH2D *h_nSector_VS_Size_CND1_epCDn;

    TH1D *h_Size_CND1_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND1_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_epFDn;
    TH2D *h_P_n_VS_Size_CND1_epFDn;
    TH2D *h_theta_n_VS_Size_CND1_epFDn;
    TH2D *h_phi_n_VS_Size_CND1_epFDn;
    TH2D *h_P_miss_VS_Size_CND1_epFDn;
    TH2D *h_theta_miss_VS_Size_CND1_epFDn;
    TH2D *h_phi_miss_VS_Size_CND1_epFDn;
    TH2D *h_dpp_VS_Size_CND1_epFDn;
    TH2D *h_beta_n_VS_Size_CND1_epFDn;
    TH2D *h_E_miss_VS_Size_CND1_epFDn;
    TH2D *h_M_miss_VS_Size_CND1_epFDn;
    TH2D *h_path_VS_Size_CND1_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_epFDn;
    TH2D *h_ToF_VS_Size_CND1_epFDn;
    TH2D *h_nSector_VS_Size_CND1_epFDn;

    TH1D *h_Size_CND2_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND2_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_epCDn;
    TH2D *h_P_n_VS_Size_CND2_epCDn;
    TH2D *h_theta_n_VS_Size_CND2_epCDn;
    TH2D *h_phi_n_VS_Size_CND2_epCDn;
    TH2D *h_P_miss_VS_Size_CND2_epCDn;
    TH2D *h_theta_miss_VS_Size_CND2_epCDn;
    TH2D *h_phi_miss_VS_Size_CND2_epCDn;
    TH2D *h_dpp_VS_Size_CND2_epCDn;
    TH2D *h_beta_n_VS_Size_CND2_epCDn;
    TH2D *h_E_miss_VS_Size_CND2_epCDn;
    TH2D *h_M_miss_VS_Size_CND2_epCDn;
    TH2D *h_path_VS_Size_CND2_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_epCDn;
    TH2D *h_ToF_VS_Size_CND2_epCDn;
    TH2D *h_nSector_VS_Size_CND2_epCDn;

    TH1D *h_Size_CND2_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND2_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_epFDn;
    TH2D *h_P_n_VS_Size_CND2_epFDn;
    TH2D *h_theta_n_VS_Size_CND2_epFDn;
    TH2D *h_phi_n_VS_Size_CND2_epFDn;
    TH2D *h_P_miss_VS_Size_CND2_epFDn;
    TH2D *h_theta_miss_VS_Size_CND2_epFDn;
    TH2D *h_phi_miss_VS_Size_CND2_epFDn;
    TH2D *h_dpp_VS_Size_CND2_epFDn;
    TH2D *h_beta_n_VS_Size_CND2_epFDn;
    TH2D *h_E_miss_VS_Size_CND2_epFDn;
    TH2D *h_M_miss_VS_Size_CND2_epFDn;
    TH2D *h_path_VS_Size_CND2_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_epFDn;
    TH2D *h_ToF_VS_Size_CND2_epFDn;
    TH2D *h_nSector_VS_Size_CND2_epFDn;

    TH1D *h_Size_CND3_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND3_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_epCDn;
    TH2D *h_P_n_VS_Size_CND3_epCDn;
    TH2D *h_theta_n_VS_Size_CND3_epCDn;
    TH2D *h_phi_n_VS_Size_CND3_epCDn;
    TH2D *h_P_miss_VS_Size_CND3_epCDn;
    TH2D *h_theta_miss_VS_Size_CND3_epCDn;
    TH2D *h_phi_miss_VS_Size_CND3_epCDn;
    TH2D *h_dpp_VS_Size_CND3_epCDn;
    TH2D *h_beta_n_VS_Size_CND3_epCDn;
    TH2D *h_E_miss_VS_Size_CND3_epCDn;
    TH2D *h_M_miss_VS_Size_CND3_epCDn;
    TH2D *h_path_VS_Size_CND3_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_epCDn;
    TH2D *h_ToF_VS_Size_CND3_epCDn;
    TH2D *h_nSector_VS_Size_CND3_epCDn;

    TH1D *h_Size_CND3_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND3_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_epFDn;
    TH2D *h_P_n_VS_Size_CND3_epFDn;
    TH2D *h_theta_n_VS_Size_CND3_epFDn;
    TH2D *h_phi_n_VS_Size_CND3_epFDn;
    TH2D *h_P_miss_VS_Size_CND3_epFDn;
    TH2D *h_theta_miss_VS_Size_CND3_epFDn;
    TH2D *h_phi_miss_VS_Size_CND3_epFDn;
    TH2D *h_dpp_VS_Size_CND3_epFDn;
    TH2D *h_beta_n_VS_Size_CND3_epFDn;
    TH2D *h_E_miss_VS_Size_CND3_epFDn;
    TH2D *h_M_miss_VS_Size_CND3_epFDn;
    TH2D *h_path_VS_Size_CND3_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_epFDn;
    TH2D *h_ToF_VS_Size_CND3_epFDn;
    TH2D *h_nSector_VS_Size_CND3_epFDn;

    TH2D *h_Size_CND1_VS_Size_CND2_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_epCDn;

    TH2D *h_Size_CND1_VS_Size_CND2_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_epFDn;

    TH1D *h_LayerMult_CND1_epCDn;
    TH1D *h_LayerMult_CND2_epCDn;
    TH1D *h_LayerMult_CND3_epCDn;

    TH1D *h_LayerMult_CND1_epFDn;
    TH1D *h_LayerMult_CND2_epFDn;
    TH1D *h_LayerMult_CND3_epFDn;

    TH2D *h_LayerMult_CND1_VS_LayerMult_CND2_epCDn;
    TH2D *h_LayerMult_CND1_VS_LayerMult_CND3_epCDn;
    TH2D *h_LayerMult_CND2_VS_LayerMult_CND3_epCDn;

    TH2D *h_LayerMult_CND1_VS_LayerMult_CND2_epFDn;
    TH2D *h_LayerMult_CND1_VS_LayerMult_CND3_epFDn;
    TH2D *h_LayerMult_CND2_VS_LayerMult_CND3_epFDn;

    TH1D *h_ToF_n_epCDn;
    TH1D *h_ToF_zoomout_epCDn;
    TH2D *h_P_n_VS_ToF_n_epCDn;
    TH2D *h_theta_n_VS_ToF_n_epCDn;
    TH2D *h_phi_n_VS_ToF_n_epCDn;
    TH2D *h_P_miss_VS_ToF_n_epCDn;
    TH2D *h_theta_miss_VS_ToF_n_epCDn;
    TH2D *h_phi_miss_VS_ToF_n_epCDn;
    TH2D *h_dpp_VS_ToF_n_epCDn;
    TH2D *h_beta_n_VS_ToF_n_epCDn;
    TH2D *h_E_miss_VS_ToF_n_epCDn;
    TH2D *h_M_miss_VS_ToF_n_epCDn;
    TH2D *h_path_VS_ToF_n_epCDn;
    TH2D *h_theta_n_miss_VS_ToF_n_epCDn;
    TH2D *h_nSector_VS_ToF_n_epCDn;

    TH1D *h_ToF_n_epFDn;
    TH1D *h_ToF_zoomout_epFDn;
    TH2D *h_P_n_VS_ToF_n_epFDn;
    TH2D *h_theta_n_VS_ToF_n_epFDn;
    TH2D *h_phi_n_VS_ToF_n_epFDn;
    TH2D *h_P_miss_VS_ToF_n_epFDn;
    TH2D *h_theta_miss_VS_ToF_n_epFDn;
    TH2D *h_phi_miss_VS_ToF_n_epFDn;
    TH2D *h_dpp_VS_ToF_n_epFDn;
    TH2D *h_beta_n_VS_ToF_n_epFDn;
    TH2D *h_E_miss_VS_ToF_n_epFDn;
    TH2D *h_M_miss_VS_ToF_n_epFDn;
    TH2D *h_path_VS_ToF_n_epFDn;
    TH2D *h_theta_n_miss_VS_ToF_n_epFDn;
    TH2D *h_nSector_VS_ToF_n_epFDn;

#pragma endregion /* (e,e'p) plots - end */

    // Step Zero (Andrew)
    // ======================================================================================================================================================================

#pragma region /* Step Zero (Andrew) - start */

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step0;
    TH1D *h_n_multiplicity_goodN_epCDn_Step0;
    TH1D *h_n_multiplicity_badN_epCDn_Step0;

    TH1D *h_n_multiplicity_allN_epFDn_Step0;
    TH1D *h_n_multiplicity_goodN_epFDn_Step0;
    TH1D *h_n_multiplicity_badN_epFDn_Step0;

    /* Step0 cuts */
    TH1D *h_dbeta_n_BS0C_Step0_epCDn;
    TH2D *h_dbeta_n_VS_P_n_BS0C_Step0_epCDn;
    TH2D *h_dbeta_n_VS_ToF_BS0C_Step0_epCDn;
    TH1D *h_dbeta_n_AS0C_Step0_epCDn;
    TH2D *h_dbeta_n_VS_P_n_AS0C_Step0_epCDn;
    TH2D *h_dbeta_n_VS_ToF_AS0C_Step0_epCDn;

    TH1D *h_dbeta_n_BS0C_Step0_epFDn;
    TH2D *h_dbeta_n_VS_P_n_BS0C_Step0_epFDn;
    TH2D *h_dbeta_n_VS_ToF_BS0C_Step0_epFDn;
    TH1D *h_dbeta_n_AS0C_Step0_epFDn;
    TH2D *h_dbeta_n_VS_P_n_AS0C_Step0_epFDn;
    TH2D *h_dbeta_n_VS_ToF_AS0C_Step0_epFDn;

    TH1D *h_Vhit_z_n_BS0C_Step0_epCDn;
    TH1D *h_Vhit_z_n_AS0C_Step0_epCDn;

    TH1D *h_Vhit_z_n_BS0C_Step0_epFDn;
    TH1D *h_Vhit_z_n_AS0C_Step0_epFDn;

    TH1D *h_ToF_n_BS0C_Step0_epCDn;
    TH1D *h_ToF_n_AS0C_Step0_epCDn;

    TH1D *h_ToF_n_BS0C_Step0_epFDn;
    TH1D *h_ToF_n_AS0C_Step0_epFDn;

    /* Kinematical variables */
    TH1D *h_theta_n_goodN_Step0_epCDn;
    TH1D *h_theta_n_badN_Step0_epCDn;
    TH1D *h_phi_n_goodN_Step0_epCDn;
    TH1D *h_phi_n_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_phi_n_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_phi_n_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_beta_n_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_beta_n_badN_Step0_epCDn;

    TH1D *h_theta_n_goodN_Step0_epFDn;
    TH1D *h_theta_n_badN_Step0_epFDn;
    TH1D *h_phi_n_goodN_Step0_epFDn;
    TH1D *h_phi_n_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_phi_n_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_phi_n_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_beta_n_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_beta_n_badN_Step0_epFDn;

    TH1D *h_P_n_goodN_Step0_epCDn;
    TH1D *h_P_n_badN_Step0_epCDn;
    TH2D *h_P_n_VS_theta_n_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_theta_n_badN_Step0_epCDn;

    TH1D *h_P_n_goodN_Step0_epFDn;
    TH1D *h_P_n_badN_Step0_epFDn;
    TH2D *h_P_n_VS_theta_n_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_theta_n_badN_Step0_epFDn;

    TH1D *h_P_miss_goodN_Step0_epCDn;
    TH1D *h_P_miss_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_theta_miss_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_theta_miss_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_phi_miss_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_phi_miss_badN_Step0_epCDn;

    TH1D *h_P_miss_goodN_Step0_epFDn;
    TH1D *h_P_miss_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_theta_miss_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_theta_miss_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_phi_miss_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_phi_miss_badN_Step0_epFDn;

    TH1D *h_dpp_allN_Step0_epCDn;
    TH1D *h_dpp_goodN_Step0_epCDn;
    TH1D *h_dpp_badN_Step0_epCDn;
    TH1D *h_dpp_allN_for_theta_n_miss_less_than_25_Step0_epCDn;

    TH1D *h_dpp_allN_Step0_epFDn;
    TH1D *h_dpp_goodN_Step0_epFDn;
    TH1D *h_dpp_badN_Step0_epFDn;
    TH1D *h_dpp_allN_for_theta_n_miss_less_than_25_Step0_epFDn;

    TH1D *h_theta_n_miss_allN_Step0_epCDn;
    TH1D *h_theta_n_miss_goodN_Step0_epCDn;
    TH1D *h_theta_n_miss_badN_Step0_epCDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_05_Step0_epCDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_03_Step0_epCDn;

    TH1D *h_theta_n_miss_allN_Step0_epFDn;
    TH1D *h_theta_n_miss_goodN_Step0_epFDn;
    TH1D *h_theta_n_miss_badN_Step0_epFDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_05_Step0_epFDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_03_Step0_epFDn;

    TH2D *h_dpp_VS_theta_n_miss_allN_Step0_epCDn;

    TH2D *h_dpp_VS_theta_n_miss_allN_Step0_epFDn;

    TH1D *h_E_p_goodN_Step0_epCDn;
    TH1D *h_E_p_badN_Step0_epCDn;
    TH1D *h_E_miss_goodN_Step0_epCDn;
    TH1D *h_E_miss_badN_Step0_epCDn;
    TH1D *h_M_miss_goodN_Step0_epCDn;
    TH1D *h_M_miss_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_P_n_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_P_n_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_theta_n_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_theta_n_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_phi_n_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_phi_n_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_P_miss_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_P_miss_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_theta_miss_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_theta_miss_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_phi_miss_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_phi_miss_badN_Step0_epCDn;

    TH1D *h_E_p_goodN_Step0_epFDn;
    TH1D *h_E_p_badN_Step0_epFDn;
    TH1D *h_E_miss_goodN_Step0_epFDn;
    TH1D *h_E_miss_badN_Step0_epFDn;
    TH1D *h_M_miss_goodN_Step0_epFDn;
    TH1D *h_M_miss_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_P_n_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_P_n_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_theta_n_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_theta_n_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_phi_n_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_phi_n_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_P_miss_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_P_miss_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_theta_miss_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_theta_miss_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_phi_miss_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_phi_miss_badN_Step0_epFDn;

    TH1D *h_P_n_minus_P_miss_goodN_Step0_epCDn;
    TH1D *h_P_n_minus_P_miss_badN_Step0_epCDn;
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step0_epCDn;
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step0_epCDn;
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step0_epCDn;
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step0_epCDn;
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step0_epCDn;
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step0_epCDn;

    TH1D *h_P_n_minus_P_miss_goodN_Step0_epFDn;
    TH1D *h_P_n_minus_P_miss_badN_Step0_epFDn;
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step0_epFDn;
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step0_epFDn;
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step0_epFDn;
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step0_epFDn;
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step0_epFDn;
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step0_epFDn;

    TH2D *h_P_n_VS_P_miss_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_P_miss_badN_Step0_epCDn;
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step0_epCDn;
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step0_epCDn;
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step0_epCDn;
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step0_epCDn;
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step0_epCDn;
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step0_epCDn;

    TH2D *h_P_n_VS_P_miss_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_P_miss_badN_Step0_epFDn;
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step0_epFDn;
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step0_epFDn;
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step0_epFDn;
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step0_epFDn;
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step0_epFDn;
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step0_epFDn;

    TH1D *h_theta_n_p_goodN_Step0_epCDn;
    TH1D *h_theta_n_p_badN_Step0_epCDn;
    TH2D *h_theta_n_p_VS_P_p_goodN_Step0_epCDn;
    TH2D *h_theta_n_p_VS_P_p_badN_Step0_epCDn;

    TH1D *h_theta_n_p_goodN_Step0_epFDn;
    TH1D *h_theta_n_p_badN_Step0_epFDn;
    TH2D *h_theta_n_p_VS_P_p_goodN_Step0_epFDn;
    TH2D *h_theta_n_p_VS_P_p_badN_Step0_epFDn;

    TH1D *h_xB_goodN_Step0_epCDn;
    TH1D *h_xB_badN_Step0_epCDn;

    TH1D *h_xB_goodN_Step0_epFDn;
    TH1D *h_xB_badN_Step0_epFDn;

    /* Detector responses */
    TH1D *h_Edep_CND_goodN_Step0_epCDn;
    TH1D *h_Edep_CND_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step0_epCDn;

    TH1D *h_Edep_CND_goodN_Step0_epFDn;
    TH1D *h_Edep_CND_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step0_epFDn;

    TH1D *h_Edep_CTOF_goodN_Step0_epCDn;
    TH1D *h_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step0_epCDn;

    TH1D *h_Edep_CTOF_goodN_Step0_epFDn;
    TH1D *h_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step0_epFDn;

    TH1D *h_Edep_single_goodN_Step0_epCDn;
    TH1D *h_Edep_single_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_path_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_path_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_single_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_single_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_single_badN_Step0_epCDn;

    TH1D *h_Edep_single_goodN_Step0_epFDn;
    TH1D *h_Edep_single_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_path_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_path_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_single_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_single_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_single_badN_Step0_epFDn;

    TH1D *h_Edep_CND1_goodN_Step0_epCDn;
    TH1D *h_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step0_epCDn;

    TH1D *h_Edep_CND1_goodN_Step0_epFDn;
    TH1D *h_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step0_epFDn;

    TH1D *h_Edep_CND2_goodN_Step0_epCDn;
    TH1D *h_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CND2_badN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step0_epCDn;

    TH1D *h_Edep_CND2_goodN_Step0_epFDn;
    TH1D *h_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CND2_badN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step0_epFDn;

    TH1D *h_Edep_CND3_goodN_Step0_epCDn;
    TH1D *h_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_path_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Edep_CND3_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Edep_CND3_badN_Step0_epCDn;

    TH1D *h_Edep_CND3_goodN_Step0_epFDn;
    TH1D *h_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_path_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Edep_CND3_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Edep_CND3_badN_Step0_epFDn;

    TH1D *h_Size_CND1_goodN_Step0_epCDn;
    TH1D *h_Size_CND1_badN_Step0_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_path_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_path_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Size_CND1_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Size_CND1_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Size_CND1_badN_Step0_epCDn;

    TH1D *h_Size_CND1_goodN_Step0_epFDn;
    TH1D *h_Size_CND1_badN_Step0_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_path_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_path_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Size_CND1_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Size_CND1_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Size_CND1_badN_Step0_epFDn;

    TH1D *h_Size_CND2_goodN_Step0_epCDn;
    TH1D *h_Size_CND2_badN_Step0_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_path_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_path_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Size_CND2_badN_Step0_epCDn;

    TH1D *h_Size_CND2_goodN_Step0_epFDn;
    TH1D *h_Size_CND2_badN_Step0_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_path_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_path_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Size_CND2_badN_Step0_epFDn;

    TH1D *h_Size_CND3_goodN_Step0_epCDn;
    TH1D *h_Size_CND3_badN_Step0_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_P_n_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_dpp_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_E_p_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_path_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_path_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_ToF_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_ToF_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_nSector_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_Size_CND3_badN_Step0_epCDn;

    TH1D *h_Size_CND3_goodN_Step0_epFDn;
    TH1D *h_Size_CND3_badN_Step0_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_P_n_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_dpp_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_E_p_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_path_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_path_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_ToF_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_ToF_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_nSector_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_Size_CND3_badN_Step0_epFDn;

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step0_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step0_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step0_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step0_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step0_epCDn;

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step0_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step0_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step0_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step0_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step0_epFDn;

    TH1D *h_ToF_goodN_Step0_epCDn;
    TH1D *h_ToF_badN_Step0_epCDn;
    TH2D *h_P_n_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_P_n_VS_ToF_badN_Step0_epCDn;
    TH2D *h_theta_n_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_theta_n_VS_ToF_badN_Step0_epCDn;
    TH2D *h_phi_n_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_phi_n_VS_ToF_badN_Step0_epCDn;
    TH2D *h_P_miss_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_P_miss_VS_ToF_badN_Step0_epCDn;
    TH2D *h_theta_miss_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_theta_miss_VS_ToF_badN_Step0_epCDn;
    TH2D *h_phi_miss_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_phi_miss_VS_ToF_badN_Step0_epCDn;
    TH2D *h_dpp_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_dpp_VS_ToF_badN_Step0_epCDn;
    TH2D *h_beta_n_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_beta_n_VS_ToF_badN_Step0_epCDn;
    TH2D *h_E_p_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_E_p_VS_ToF_badN_Step0_epCDn;
    TH2D *h_E_miss_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_E_miss_VS_ToF_badN_Step0_epCDn;
    TH2D *h_M_miss_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_M_miss_VS_ToF_badN_Step0_epCDn;
    TH2D *h_path_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_path_VS_ToF_badN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_theta_n_miss_VS_ToF_badN_Step0_epCDn;
    TH2D *h_nSector_VS_ToF_goodN_Step0_epCDn;
    TH2D *h_nSector_VS_ToF_badN_Step0_epCDn;

    TH1D *h_ToF_goodN_Step0_epFDn;
    TH1D *h_ToF_badN_Step0_epFDn;
    TH2D *h_P_n_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_P_n_VS_ToF_badN_Step0_epFDn;
    TH2D *h_theta_n_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_theta_n_VS_ToF_badN_Step0_epFDn;
    TH2D *h_phi_n_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_phi_n_VS_ToF_badN_Step0_epFDn;
    TH2D *h_P_miss_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_P_miss_VS_ToF_badN_Step0_epFDn;
    TH2D *h_theta_miss_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_theta_miss_VS_ToF_badN_Step0_epFDn;
    TH2D *h_phi_miss_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_phi_miss_VS_ToF_badN_Step0_epFDn;
    TH2D *h_dpp_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_dpp_VS_ToF_badN_Step0_epFDn;
    TH2D *h_beta_n_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_beta_n_VS_ToF_badN_Step0_epFDn;
    TH2D *h_E_p_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_E_p_VS_ToF_badN_Step0_epFDn;
    TH2D *h_E_miss_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_E_miss_VS_ToF_badN_Step0_epFDn;
    TH2D *h_M_miss_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_M_miss_VS_ToF_badN_Step0_epFDn;
    TH2D *h_path_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_path_VS_ToF_badN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_theta_n_miss_VS_ToF_badN_Step0_epFDn;
    TH2D *h_nSector_VS_ToF_goodN_Step0_epFDn;
    TH2D *h_nSector_VS_ToF_badN_Step0_epFDn;

    TH1D *h_beta_n_goodN_Step0_epCDn;
    TH1D *h_beta_n_badN_Step0_epCDn;

    TH1D *h_beta_n_goodN_Step0_epFDn;
    TH1D *h_beta_n_badN_Step0_epFDn;

#pragma endregion /* Step Zero (Andrew) - end */

    // Step One (After Edep_CND Cut) (Andrew)
    // ======================================================================================================================================================================

#pragma region /* Step One (After Edep_CND Cut) (Andrew) - start */

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step1;
    TH1D *h_n_multiplicity_goodN_epCDn_Step1;
    TH1D *h_n_multiplicity_badN_epCDn_Step1;

    TH1D *h_n_multiplicity_allN_epFDn_Step1;
    TH1D *h_n_multiplicity_goodN_epFDn_Step1;
    TH1D *h_n_multiplicity_badN_epFDn_Step1;

    /* Step1 cuts */
    /*
    TH2D *h_dbeta_n_VS_P_n_BS1C_Step1_epCDn;
    TH2D *h_dbeta_n_VS_ToF_BS1C_Step1_epCDn;
    TH2D *h_dbeta_n_VS_P_n_AS1C_Step1_epCDn;
    TH2D *h_dbeta_n_VS_ToF_AS1C_Step1_epCDn;

    TH2D *h_dbeta_n_VS_P_n_BS1C_Step1_epFDn;
    TH2D *h_dbeta_n_VS_ToF_BS1C_Step1_epFDn;
    TH2D *h_dbeta_n_VS_P_n_AS1C_Step1_epFDn;
    TH2D *h_dbeta_n_VS_ToF_AS1C_Step1_epFDn;

    TH1D *h_Vhit_z_n_BS1C_Step1_epCDn;
    TH1D *h_Vhit_z_n_AS1C_Step1_epCDn;

    TH1D *h_Vhit_z_n_BS1C_Step1_epFDn;
    TH1D *h_Vhit_z_n_AS1C_Step1_epFDn;

    TH1D *h_ToF_n_BS1C_Step1_epCDn;
    TH1D *h_ToF_n_AS1C_Step1_epCDn;

    TH1D *h_ToF_n_BS1C_Step1_epFDn;
    TH1D *h_ToF_n_AS1C_Step1_epFDn;

    TH1D *h_beta_n_BS1C_Step1_epFDn;
    TH1D *h_beta_n_AS1C_Step1_epFDn;
 */

    /* Kinematical variables */
    TH1D *h_theta_n_goodN_Step1_epCDn;
    TH1D *h_theta_n_badN_Step1_epCDn;
    TH1D *h_phi_n_goodN_Step1_epCDn;
    TH1D *h_phi_n_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_phi_n_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_phi_n_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_beta_n_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_beta_n_badN_Step1_epCDn;

    TH1D *h_theta_n_goodN_Step1_epFDn;
    TH1D *h_theta_n_badN_Step1_epFDn;
    TH1D *h_phi_n_goodN_Step1_epFDn;
    TH1D *h_phi_n_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_phi_n_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_phi_n_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_beta_n_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_beta_n_badN_Step1_epFDn;

    TH1D *h_P_n_goodN_Step1_epCDn;
    TH1D *h_P_n_badN_Step1_epCDn;
    TH2D *h_P_n_VS_theta_n_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_theta_n_badN_Step1_epCDn;

    TH1D *h_P_n_goodN_Step1_epFDn;
    TH1D *h_P_n_badN_Step1_epFDn;
    TH2D *h_P_n_VS_theta_n_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_theta_n_badN_Step1_epFDn;

    TH1D *h_P_miss_goodN_Step1_epCDn;
    TH1D *h_P_miss_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_theta_miss_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_theta_miss_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_phi_miss_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_phi_miss_badN_Step1_epCDn;

    TH1D *h_P_miss_goodN_Step1_epFDn;
    TH1D *h_P_miss_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_theta_miss_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_theta_miss_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_phi_miss_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_phi_miss_badN_Step1_epFDn;

    TH1D *h_dpp_allN_Step1_epCDn;
    TH1D *h_dpp_goodN_Step1_epCDn;
    TH1D *h_dpp_badN_Step1_epCDn;
    TH1D *h_dpp_allN_for_theta_n_miss_less_than_25_Step1_epCDn;

    TH1D *h_dpp_allN_Step1_epFDn;
    TH1D *h_dpp_goodN_Step1_epFDn;
    TH1D *h_dpp_badN_Step1_epFDn;
    TH1D *h_dpp_allN_for_theta_n_miss_less_than_25_Step1_epFDn;

    TH1D *h_theta_n_miss_allN_Step1_epCDn;
    TH1D *h_theta_n_miss_goodN_Step1_epCDn;
    TH1D *h_theta_n_miss_badN_Step1_epCDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_05_Step1_epCDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_03_Step1_epCDn;

    TH1D *h_theta_n_miss_allN_Step1_epFDn;
    TH1D *h_theta_n_miss_goodN_Step1_epFDn;
    TH1D *h_theta_n_miss_badN_Step1_epFDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_05_Step1_epFDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_03_Step1_epFDn;

    TH2D *h_dpp_VS_theta_n_miss_allN_Step1_epCDn;

    TH2D *h_dpp_VS_theta_n_miss_allN_Step1_epFDn;

    TH1D *h_E_p_goodN_Step1_epCDn;
    TH1D *h_E_p_badN_Step1_epCDn;
    TH1D *h_E_miss_goodN_Step1_epCDn;
    TH1D *h_E_miss_badN_Step1_epCDn;
    TH1D *h_M_miss_goodN_Step1_epCDn;
    TH1D *h_M_miss_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_P_n_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_P_n_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_theta_n_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_theta_n_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_phi_n_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_phi_n_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_P_miss_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_P_miss_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_theta_miss_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_theta_miss_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_phi_miss_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_phi_miss_badN_Step1_epCDn;

    TH1D *h_E_p_goodN_Step1_epFDn;
    TH1D *h_E_p_badN_Step1_epFDn;
    TH1D *h_E_miss_goodN_Step1_epFDn;
    TH1D *h_E_miss_badN_Step1_epFDn;
    TH1D *h_M_miss_goodN_Step1_epFDn;
    TH1D *h_M_miss_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_P_n_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_P_n_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_theta_n_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_theta_n_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_phi_n_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_phi_n_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_P_miss_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_P_miss_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_theta_miss_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_theta_miss_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_phi_miss_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_phi_miss_badN_Step1_epFDn;

    TH1D *h_P_n_minus_P_miss_goodN_Step1_epCDn;
    TH1D *h_P_n_minus_P_miss_badN_Step1_epCDn;
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step1_epCDn;
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step1_epCDn;
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step1_epCDn;
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step1_epCDn;
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step1_epCDn;
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step1_epCDn;

    TH1D *h_P_n_minus_P_miss_goodN_Step1_epFDn;
    TH1D *h_P_n_minus_P_miss_badN_Step1_epFDn;
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step1_epFDn;
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step1_epFDn;
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step1_epFDn;
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step1_epFDn;
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step1_epFDn;
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step1_epFDn;

    TH2D *h_P_n_VS_P_miss_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_P_miss_badN_Step1_epCDn;
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step1_epCDn;
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step1_epCDn;
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step1_epCDn;
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step1_epCDn;
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step1_epCDn;
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step1_epCDn;

    TH2D *h_P_n_VS_P_miss_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_P_miss_badN_Step1_epFDn;
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step1_epFDn;
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step1_epFDn;
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step1_epFDn;
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step1_epFDn;
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step1_epFDn;
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step1_epFDn;

    TH1D *h_theta_n_p_goodN_Step1_epCDn;
    TH1D *h_theta_n_p_badN_Step1_epCDn;
    TH2D *h_theta_n_p_VS_P_p_goodN_Step1_epCDn;
    TH2D *h_theta_n_p_VS_P_p_badN_Step1_epCDn;

    TH1D *h_theta_n_p_goodN_Step1_epFDn;
    TH1D *h_theta_n_p_badN_Step1_epFDn;
    TH2D *h_theta_n_p_VS_P_p_goodN_Step1_epFDn;
    TH2D *h_theta_n_p_VS_P_p_badN_Step1_epFDn;

    TH1D *h_xB_goodN_Step1_epCDn;
    TH1D *h_xB_badN_Step1_epCDn;

    TH1D *h_xB_goodN_Step1_epFDn;
    TH1D *h_xB_badN_Step1_epFDn;

    /* Detector responses */
    TH1D *h_Edep_CND_goodN_Step1_epCDn;
    TH1D *h_Edep_CND_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step1_epCDn;

    TH1D *h_Edep_CND_goodN_Step1_epFDn;
    TH1D *h_Edep_CND_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step1_epFDn;

    TH1D *h_Edep_CTOF_goodN_Step1_epCDn;
    TH1D *h_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step1_epCDn;

    TH1D *h_Edep_CTOF_goodN_Step1_epFDn;
    TH1D *h_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step1_epFDn;

    TH1D *h_Edep_single_goodN_Step1_epCDn;
    TH1D *h_Edep_single_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_path_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_path_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_single_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_single_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_single_badN_Step1_epCDn;

    TH1D *h_Edep_single_goodN_Step1_epFDn;
    TH1D *h_Edep_single_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_path_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_path_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_single_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_single_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_single_badN_Step1_epFDn;

    TH1D *h_Edep_CND1_goodN_Step1_epCDn;
    TH1D *h_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step1_epCDn;

    TH1D *h_Edep_CND1_goodN_Step1_epFDn;
    TH1D *h_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step1_epFDn;

    TH1D *h_Edep_CND2_goodN_Step1_epCDn;
    TH1D *h_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CND2_badN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step1_epCDn;

    TH1D *h_Edep_CND2_goodN_Step1_epFDn;
    TH1D *h_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CND2_badN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step1_epFDn;

    TH1D *h_Edep_CND3_goodN_Step1_epCDn;
    TH1D *h_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_path_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Edep_CND3_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Edep_CND3_badN_Step1_epCDn;

    TH1D *h_Edep_CND3_goodN_Step1_epFDn;
    TH1D *h_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_path_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Edep_CND3_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Edep_CND3_badN_Step1_epFDn;

    TH1D *h_Size_CND1_goodN_Step1_epCDn;
    TH1D *h_Size_CND1_badN_Step1_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_path_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_path_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Size_CND1_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Size_CND1_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Size_CND1_badN_Step1_epCDn;

    TH1D *h_Size_CND1_goodN_Step1_epFDn;
    TH1D *h_Size_CND1_badN_Step1_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_path_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_path_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Size_CND1_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Size_CND1_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Size_CND1_badN_Step1_epFDn;

    TH1D *h_Size_CND2_goodN_Step1_epCDn;
    TH1D *h_Size_CND2_badN_Step1_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_path_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_path_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Size_CND2_badN_Step1_epCDn;

    TH1D *h_Size_CND2_goodN_Step1_epFDn;
    TH1D *h_Size_CND2_badN_Step1_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_path_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_path_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Size_CND2_badN_Step1_epFDn;

    TH1D *h_Size_CND3_goodN_Step1_epCDn;
    TH1D *h_Size_CND3_badN_Step1_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_P_n_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_dpp_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_E_p_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_path_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_path_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_ToF_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_ToF_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_nSector_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_Size_CND3_badN_Step1_epCDn;

    TH1D *h_Size_CND3_goodN_Step1_epFDn;
    TH1D *h_Size_CND3_badN_Step1_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_P_n_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_dpp_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_E_p_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_path_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_path_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_ToF_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_ToF_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_nSector_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_Size_CND3_badN_Step1_epFDn;

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step1_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step1_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step1_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step1_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step1_epCDn;

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step1_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step1_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step1_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step1_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step1_epFDn;

    TH1D *h_ToF_goodN_Step1_epCDn;
    TH1D *h_ToF_badN_Step1_epCDn;
    TH2D *h_P_n_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_P_n_VS_ToF_badN_Step1_epCDn;
    TH2D *h_theta_n_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_theta_n_VS_ToF_badN_Step1_epCDn;
    TH2D *h_phi_n_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_phi_n_VS_ToF_badN_Step1_epCDn;
    TH2D *h_P_miss_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_P_miss_VS_ToF_badN_Step1_epCDn;
    TH2D *h_theta_miss_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_theta_miss_VS_ToF_badN_Step1_epCDn;
    TH2D *h_phi_miss_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_phi_miss_VS_ToF_badN_Step1_epCDn;
    TH2D *h_dpp_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_dpp_VS_ToF_badN_Step1_epCDn;
    TH2D *h_beta_n_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_beta_n_VS_ToF_badN_Step1_epCDn;
    TH2D *h_E_p_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_E_p_VS_ToF_badN_Step1_epCDn;
    TH2D *h_E_miss_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_E_miss_VS_ToF_badN_Step1_epCDn;
    TH2D *h_M_miss_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_M_miss_VS_ToF_badN_Step1_epCDn;
    TH2D *h_path_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_path_VS_ToF_badN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_theta_n_miss_VS_ToF_badN_Step1_epCDn;
    TH2D *h_nSector_VS_ToF_goodN_Step1_epCDn;
    TH2D *h_nSector_VS_ToF_badN_Step1_epCDn;

    TH1D *h_ToF_goodN_Step1_epFDn;
    TH1D *h_ToF_badN_Step1_epFDn;
    TH2D *h_P_n_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_P_n_VS_ToF_badN_Step1_epFDn;
    TH2D *h_theta_n_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_theta_n_VS_ToF_badN_Step1_epFDn;
    TH2D *h_phi_n_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_phi_n_VS_ToF_badN_Step1_epFDn;
    TH2D *h_P_miss_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_P_miss_VS_ToF_badN_Step1_epFDn;
    TH2D *h_theta_miss_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_theta_miss_VS_ToF_badN_Step1_epFDn;
    TH2D *h_phi_miss_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_phi_miss_VS_ToF_badN_Step1_epFDn;
    TH2D *h_dpp_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_dpp_VS_ToF_badN_Step1_epFDn;
    TH2D *h_beta_n_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_beta_n_VS_ToF_badN_Step1_epFDn;
    TH2D *h_E_p_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_E_p_VS_ToF_badN_Step1_epFDn;
    TH2D *h_E_miss_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_E_miss_VS_ToF_badN_Step1_epFDn;
    TH2D *h_M_miss_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_M_miss_VS_ToF_badN_Step1_epFDn;
    TH2D *h_path_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_path_VS_ToF_badN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_theta_n_miss_VS_ToF_badN_Step1_epFDn;
    TH2D *h_nSector_VS_ToF_goodN_Step1_epFDn;
    TH2D *h_nSector_VS_ToF_badN_Step1_epFDn;

    TH1D *h_beta_n_goodN_Step1_epCDn;
    TH1D *h_beta_n_badN_Step1_epCDn;

    TH1D *h_beta_n_goodN_Step1_epFDn;
    TH1D *h_beta_n_badN_Step1_epFDn;

#pragma endregion /* Step One (After Edep_CND Cut) (Andrew) - end */

    // Step Two (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================

#pragma region /* Step Two (After applying Phi Diff Charge Track cut) (Andrew) - start */

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step2;
    TH1D *h_n_multiplicity_goodN_epCDn_Step2;
    TH1D *h_n_multiplicity_badN_epCDn_Step2;

    TH1D *h_n_multiplicity_allN_epFDn_Step2;
    TH1D *h_n_multiplicity_goodN_epFDn_Step2;
    TH1D *h_n_multiplicity_badN_epFDn_Step2;

    /* Step2 prep plots */
    /* ToF * c - v_hit_3v.Z() plots */
    TH1D *h_ToF_c_minus_VhitZ_BC_allN_Step2prep_epCDn;
    TH1D *h_ToF_c_minus_VhitZ_BC_goodN_Step2prep_epCDn;
    TH1D *h_ToF_c_minus_VhitZ_BC_badN_Step2prep_epCDn;
    // TH1D *h_ToF_c_minus_VhitZ_AC_allN_Step2prep_epCDn;
    // TH1D *h_ToF_c_minus_VhitZ_AC_goodN_Step2prep_epCDn;
    // TH1D *h_ToF_c_minus_VhitZ_AC_badN_Step2prep_epCDn;

    TH1D *h_ToF_c_minus_VhitZ_BC_allN_Step2prep_epFDn;
    TH1D *h_ToF_c_minus_VhitZ_BC_goodN_Step2prep_epFDn;
    TH1D *h_ToF_c_minus_VhitZ_BC_badN_Step2prep_epFDn;
    // TH1D *h_ToF_c_minus_VhitZ_AC_allN_Step2prep_epFDn;
    // TH1D *h_ToF_c_minus_VhitZ_AC_goodN_Step2prep_epFDn;
    // TH1D *h_ToF_c_minus_VhitZ_AC_badN_Step2prep_epFDn;

    TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_BC_allN_Step2prep_epCDn;
    TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_BC_goodN_Step2prep_epCDn;
    TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_BC_badN_Step2prep_epCDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_AC_allN_Step2prep_epCDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_AC_goodN_Step2prep_epCDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_AC_badN_Step2prep_epCDn;

    TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_BC_allN_Step2prep_epFDn;
    TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_BC_goodN_Step2prep_epFDn;
    TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_BC_badN_Step2prep_epFDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_AC_allN_Step2prep_epFDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_AC_goodN_Step2prep_epFDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_VhitZ_AC_badN_Step2prep_epFDn;

    TH2D *h_ToF_c_minus_VhitZ_VS_ToF_BC_allN_Step2prep_epCDn;
    TH2D *h_ToF_c_minus_VhitZ_VS_ToF_BC_goodN_Step2prep_epCDn;
    TH2D *h_ToF_c_minus_VhitZ_VS_ToF_BC_badN_Step2prep_epCDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_ToF_AC_allN_Step2prep_epCDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_ToF_AC_goodN_Step2prep_epCDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_ToF_AC_badN_Step2prep_epCDn;

    TH2D *h_ToF_c_minus_VhitZ_VS_ToF_BC_allN_Step2prep_epFDn;
    TH2D *h_ToF_c_minus_VhitZ_VS_ToF_BC_goodN_Step2prep_epFDn;
    TH2D *h_ToF_c_minus_VhitZ_VS_ToF_BC_badN_Step2prep_epFDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_ToF_AC_allN_Step2prep_epFDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_ToF_AC_goodN_Step2prep_epFDn;
    // TH2D *h_ToF_c_minus_VhitZ_VS_ToF_AC_badN_Step2prep_epFDn;

    // TH1D *h_Edep_CND_goodN_Step2prep_epCDn;
    // TH1D *h_Edep_CND_badN_Step2prep_epCDn;

    // TH1D *h_Edep_CND_goodN_Step2prep_epFDn;
    // TH1D *h_Edep_CND_badN_Step2prep_epFDn;

    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step2prep_epCDn;
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step2prep_epCDn;

    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step2prep_epFDn;
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step2prep_epFDn;

    TH1D *h_Edep_CND_goodN_withNearbyPos_Step2prep_epCDn;
    TH1D *h_Edep_CND_badN_withNearbyPos_Step2prep_epCDn;

    TH1D *h_Edep_CND_goodN_withNearbyPos_Step2prep_epFDn;
    TH1D *h_Edep_CND_badN_withNearbyPos_Step2prep_epFDn;

    TH1D *h_sdiff_pos_goodN_Step2prep_layer_epCDn[7];
    TH1D *h_sdiff_pos_badN_Step2prep_layer_epCDn[7];

    TH1D *h_sdiff_pos_goodN_Step2prep_layer_epFDn[7];
    TH1D *h_sdiff_pos_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_mom_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_mom_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_mom_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_mom_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_VhitZ_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_VhitZ_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_epFDn[7];

    TH1D *h_theta_n_goodN_Step2prep_layer_epCDn[7];
    TH1D *h_theta_n_badN_Step2prep_layer_epCDn[7];

    TH1D *h_theta_n_goodN_Step2prep_layer_epFDn[7];
    TH1D *h_theta_n_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epFDn[7];

    TH1D *h_phi_n_goodN_Step2prep_layer_epCDn[7];
    TH1D *h_phi_n_badN_Step2prep_layer_epCDn[7];

    TH1D *h_phi_n_goodN_Step2prep_layer_epFDn[7];
    TH1D *h_phi_n_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_phi_n_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_phi_n_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_phi_n_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_phi_n_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_ToF_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_ToF_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_ToF_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_ToF_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_path_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_path_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_path_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_path_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_beta_n_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_beta_n_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_beta_n_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_beta_n_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_miss_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_miss_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_miss_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_miss_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_dpp_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_dpp_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_dpp_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_dpp_badN_Step2prep_layer_epFDn[7];

    TH1D *h_dToF_goodN_Step2prep_layer_epCDn[7];
    TH1D *h_dToF_badN_Step2prep_layer_epCDn[7];

    TH1D *h_dToF_goodN_Step2prep_layer_epFDn[7];
    TH1D *h_dToF_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_dToF_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_dToF_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_dToF_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_dToF_badN_Step2prep_layer_epFDn[7];

    TH1D *h_dToF_rel_pos_goodN_Step2prep_layer_epCDn[7];
    TH1D *h_dToF_rel_pos_badN_Step2prep_layer_epCDn[7];

    TH1D *h_dToF_rel_pos_goodN_Step2prep_layer_epFDn[7];
    TH1D *h_dToF_rel_pos_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_dToF_rel_pos_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_dToF_rel_pos_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_dToF_rel_pos_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_dToF_rel_pos_badN_Step2prep_layer_epFDn[7];

    TH1D *h_dToF_rel_n_goodN_Step2prep_layer_epCDn[7];
    TH1D *h_dToF_rel_n_badN_Step2prep_layer_epCDn[7];

    TH1D *h_dToF_rel_n_goodN_Step2prep_layer_epFDn[7];
    TH1D *h_dToF_rel_n_badN_Step2prep_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_dToF_rel_n_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_dToF_rel_n_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_dToF_rel_n_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_dToF_rel_n_badN_Step2prep_layer_epFDn[7];

    TH2D *h_diff_ToFc_z_VS_Edep_noNear_goodN_Step2prep_epCDn;
    TH2D *h_diff_ToFc_z_VS_Edep_noNear_badN_Step2prep_epCDn;

    TH2D *h_diff_ToFc_z_VS_Edep_noNear_goodN_Step2prep_epFDn;
    TH2D *h_diff_ToFc_z_VS_Edep_noNear_badN_Step2prep_epFDn;

    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step2prep_epCDn;
    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_badN_Step2prep_epCDn;

    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step2prep_epFDn;
    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_badN_Step2prep_epFDn;

    /* Step2 cuts */
    TH1D *h_Size_CND1_BS2C_Step2_epCDn;
    TH1D *h_Size_CND1_AS2C_Step2_epCDn;
    TH1D *h_Size_CND2_BS2C_Step2_epCDn;
    TH1D *h_Size_CND2_AS2C_Step2_epCDn;
    TH1D *h_Size_CND3_BS2C_Step2_epCDn;
    TH1D *h_Size_CND3_AS2C_Step2_epCDn;

    TH1D *h_Size_CND1_BS2C_Step2_epFDn;
    TH1D *h_Size_CND1_AS2C_Step2_epFDn;
    TH1D *h_Size_CND2_BS2C_Step2_epFDn;
    TH1D *h_Size_CND2_AS2C_Step2_epFDn;
    TH1D *h_Size_CND3_BS2C_Step2_epFDn;
    TH1D *h_Size_CND3_AS2C_Step2_epFDn;

    TH2D *h_Size_CND1_VS_Size_CND2_BS2C_Step2_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND2_AS2C_Step2_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_BS2C_Step2_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_AS2C_Step2_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_BS2C_Step2_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_AS2C_Step2_epCDn;

    TH2D *h_Size_CND1_VS_Size_CND2_BS2C_Step2_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND2_AS2C_Step2_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_BS2C_Step2_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_AS2C_Step2_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_BS2C_Step2_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_AS2C_Step2_epFDn;

    TH1D *h_LayerMult_CND1_BS2C_Step2_epCDn;
    TH1D *h_LayerMult_CND1_AS2C_Step2_epCDn;
    TH1D *h_LayerMult_CND2_BS2C_Step2_epCDn;
    TH1D *h_LayerMult_CND2_AS2C_Step2_epCDn;
    TH1D *h_LayerMult_CND3_BS2C_Step2_epCDn;
    TH1D *h_LayerMult_CND3_AS2C_Step2_epCDn;

    TH1D *h_LayerMult_CND1_BS2C_Step2_epFDn;
    TH1D *h_LayerMult_CND1_AS2C_Step2_epFDn;
    TH1D *h_LayerMult_CND2_BS2C_Step2_epFDn;
    TH1D *h_LayerMult_CND2_AS2C_Step2_epFDn;
    TH1D *h_LayerMult_CND3_BS2C_Step2_epFDn;
    TH1D *h_LayerMult_CND3_AS2C_Step2_epFDn;

    TH2D *h_LayerMult_CND1_VS_LayerMult_CND2_BS2C_Step2_epCDn;
    TH2D *h_LayerMult_CND1_VS_LayerMult_CND2_AS2C_Step2_epCDn;
    TH2D *h_LayerMult_CND1_VS_LayerMult_CND3_BS2C_Step2_epCDn;
    TH2D *h_LayerMult_CND1_VS_LayerMult_CND3_AS2C_Step2_epCDn;
    TH2D *h_LayerMult_CND2_VS_LayerMult_CND3_BS2C_Step2_epCDn;
    TH2D *h_LayerMult_CND2_VS_LayerMult_CND3_AS2C_Step2_epCDn;

    TH2D *h_LayerMult_CND1_VS_LayerMult_CND2_BS2C_Step2_epFDn;
    TH2D *h_LayerMult_CND1_VS_LayerMult_CND2_AS2C_Step2_epFDn;
    TH2D *h_LayerMult_CND1_VS_LayerMult_CND3_BS2C_Step2_epFDn;
    TH2D *h_LayerMult_CND1_VS_LayerMult_CND3_AS2C_Step2_epFDn;
    TH2D *h_LayerMult_CND2_VS_LayerMult_CND3_BS2C_Step2_epFDn;
    TH2D *h_LayerMult_CND2_VS_LayerMult_CND3_AS2C_Step2_epFDn;
    /*
    TH2D *h_dbeta_n_VS_P_n_BS1C_Step2_epCDn;
    TH2D *h_dbeta_n_VS_ToF_BS1C_Step2_epCDn;
    TH2D *h_dbeta_n_VS_P_n_AS1C_Step2_epCDn;
    TH2D *h_dbeta_n_VS_ToF_AS1C_Step2_epCDn;

    TH2D *h_dbeta_n_VS_P_n_BS1C_Step2_epFDn;
    TH2D *h_dbeta_n_VS_ToF_BS1C_Step2_epFDn;
    TH2D *h_dbeta_n_VS_P_n_AS1C_Step2_epFDn;
    TH2D *h_dbeta_n_VS_ToF_AS1C_Step2_epFDn;

    TH1D *h_Vhit_z_n_BS1C_Step2_epCDn;
    TH1D *h_Vhit_z_n_AS1C_Step2_epCDn;

    TH1D *h_Vhit_z_n_BS1C_Step2_epFDn;
    TH1D *h_Vhit_z_n_AS1C_Step2_epFDn;

    TH1D *h_ToF_n_BS1C_Step2_epCDn;
    TH1D *h_ToF_n_AS1C_Step2_epCDn;

    TH1D *h_ToF_n_BS1C_Step2_epFDn;
    TH1D *h_ToF_n_AS1C_Step2_epFDn;

    TH1D *h_beta_n_BS1C_Step2_epFDn;
    TH1D *h_beta_n_AS1C_Step2_epFDn;
 */

    /* ToF * c - v_hit_3v.Z() plots */
    // TODO: move from here!
    // TH1D *h_ToF_c_minus_VhitZ_BC_allN_Step2_epCDn;
    // TH1D *h_ToF_c_minus_VhitZ_BC_goodN_Step2_epCDn;
    // TH1D *h_ToF_c_minus_VhitZ_BC_badN_Step2_epCDn;
    TH1D *h_ToF_c_minus_VhitZ_AC_allN_Step2_epCDn;
    TH1D *h_ToF_c_minus_VhitZ_AC_goodN_Step2_epCDn;
    TH1D *h_ToF_c_minus_VhitZ_AC_badN_Step2_epCDn;

    // TH1D *h_ToF_c_minus_VhitZ_BC_allN_Step2_epFDn;
    // TH1D *h_ToF_c_minus_VhitZ_BC_goodN_Step2_epFDn;
    // TH1D *h_ToF_c_minus_VhitZ_BC_badN_Step2_epFDn;
    TH1D *h_ToF_c_minus_VhitZ_AC_allN_Step2_epFDn;
    TH1D *h_ToF_c_minus_VhitZ_AC_goodN_Step2_epFDn;
    TH1D *h_ToF_c_minus_VhitZ_AC_badN_Step2_epFDn;

    // TH1D *h_Edep_CND_goodN_Step2_test_epCDn;
    // TH1D *h_Edep_CND_badN_Step2_test_epCDn;

    // TH1D *h_Edep_CND_goodN_Step2_test_epFDn;
    // TH1D *h_Edep_CND_badN_Step2_test_epFDn;

    /* Kinematical variables */
    TH1D *h_theta_n_goodN_Step2_epCDn;
    TH1D *h_theta_n_badN_Step2_epCDn;
    TH1D *h_phi_n_goodN_Step2_epCDn;
    TH1D *h_phi_n_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_phi_n_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_phi_n_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_beta_n_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_beta_n_badN_Step2_epCDn;

    TH1D *h_theta_n_goodN_Step2_epFDn;
    TH1D *h_theta_n_badN_Step2_epFDn;
    TH1D *h_phi_n_goodN_Step2_epFDn;
    TH1D *h_phi_n_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_phi_n_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_phi_n_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_beta_n_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_beta_n_badN_Step2_epFDn;

    TH1D *h_P_n_goodN_Step2_epCDn;
    TH1D *h_P_n_badN_Step2_epCDn;
    TH2D *h_P_n_VS_theta_n_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_theta_n_badN_Step2_epCDn;

    TH1D *h_P_n_goodN_Step2_epFDn;
    TH1D *h_P_n_badN_Step2_epFDn;
    TH2D *h_P_n_VS_theta_n_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_theta_n_badN_Step2_epFDn;

    TH1D *h_P_miss_goodN_Step2_epCDn;
    TH1D *h_P_miss_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_theta_miss_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_theta_miss_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_phi_miss_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_phi_miss_badN_Step2_epCDn;

    TH1D *h_P_miss_goodN_Step2_epFDn;
    TH1D *h_P_miss_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_theta_miss_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_theta_miss_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_phi_miss_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_phi_miss_badN_Step2_epFDn;

    TH1D *h_dpp_allN_Step2_epCDn;
    TH1D *h_dpp_goodN_Step2_epCDn;
    TH1D *h_dpp_badN_Step2_epCDn;
    TH1D *h_dpp_allN_for_theta_n_miss_less_than_25_Step2_epCDn;

    TH1D *h_dpp_allN_Step2_epFDn;
    TH1D *h_dpp_goodN_Step2_epFDn;
    TH1D *h_dpp_badN_Step2_epFDn;
    TH1D *h_dpp_allN_for_theta_n_miss_less_than_25_Step2_epFDn;

    TH1D *h_theta_n_miss_allN_Step2_epCDn;
    TH1D *h_theta_n_miss_goodN_Step2_epCDn;
    TH1D *h_theta_n_miss_badN_Step2_epCDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_05_Step2_epCDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_03_Step2_epCDn;

    TH1D *h_theta_n_miss_allN_Step2_epFDn;
    TH1D *h_theta_n_miss_goodN_Step2_epFDn;
    TH1D *h_theta_n_miss_badN_Step2_epFDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_05_Step2_epFDn;
    TH1D *h_theta_n_miss_allN_for_dpp_less_than_03_Step2_epFDn;

    TH2D *h_dpp_VS_theta_n_miss_allN_Step2_epCDn;

    TH2D *h_dpp_VS_theta_n_miss_allN_Step2_epFDn;

    TH1D *h_E_p_goodN_Step2_epCDn;
    TH1D *h_E_p_badN_Step2_epCDn;
    TH1D *h_E_miss_goodN_Step2_epCDn;
    TH1D *h_E_miss_badN_Step2_epCDn;
    TH1D *h_M_miss_goodN_Step2_epCDn;
    TH1D *h_M_miss_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_P_n_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_P_n_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_theta_n_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_theta_n_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_phi_n_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_phi_n_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_P_miss_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_P_miss_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_theta_miss_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_theta_miss_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_phi_miss_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_phi_miss_badN_Step2_epCDn;

    TH1D *h_E_p_goodN_Step2_epFDn;
    TH1D *h_E_p_badN_Step2_epFDn;
    TH1D *h_E_miss_goodN_Step2_epFDn;
    TH1D *h_E_miss_badN_Step2_epFDn;
    TH1D *h_M_miss_goodN_Step2_epFDn;
    TH1D *h_M_miss_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_P_n_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_P_n_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_theta_n_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_theta_n_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_phi_n_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_phi_n_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_P_miss_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_P_miss_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_theta_miss_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_theta_miss_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_phi_miss_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_phi_miss_badN_Step2_epFDn;

    TH1D *h_P_n_minus_P_miss_goodN_Step2_epCDn;
    TH1D *h_P_n_minus_P_miss_badN_Step2_epCDn;
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step2_epCDn;
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step2_epCDn;
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step2_epCDn;
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step2_epCDn;
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step2_epCDn;
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step2_epCDn;

    TH1D *h_P_n_minus_P_miss_goodN_Step2_epFDn;
    TH1D *h_P_n_minus_P_miss_badN_Step2_epFDn;
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step2_epFDn;
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step2_epFDn;
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step2_epFDn;
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step2_epFDn;
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step2_epFDn;
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step2_epFDn;

    TH2D *h_P_n_VS_P_miss_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_P_miss_badN_Step2_epCDn;
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step2_epCDn;
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step2_epCDn;
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step2_epCDn;
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step2_epCDn;
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step2_epCDn;
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step2_epCDn;

    TH2D *h_P_n_VS_P_miss_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_P_miss_badN_Step2_epFDn;
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step2_epFDn;
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step2_epFDn;
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step2_epFDn;
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step2_epFDn;
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step2_epFDn;
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step2_epFDn;

    TH1D *h_theta_n_p_goodN_Step2_epCDn;
    TH1D *h_theta_n_p_badN_Step2_epCDn;
    TH2D *h_theta_n_p_VS_P_p_goodN_Step2_epCDn;
    TH2D *h_theta_n_p_VS_P_p_badN_Step2_epCDn;

    TH1D *h_theta_n_p_goodN_Step2_epFDn;
    TH1D *h_theta_n_p_badN_Step2_epFDn;
    TH2D *h_theta_n_p_VS_P_p_goodN_Step2_epFDn;
    TH2D *h_theta_n_p_VS_P_p_badN_Step2_epFDn;

    TH1D *h_xB_goodN_Step2_epCDn;
    TH1D *h_xB_badN_Step2_epCDn;

    TH1D *h_xB_goodN_Step2_epFDn;
    TH1D *h_xB_badN_Step2_epFDn;

    /* Detector responses */
    TH1D *h_Edep_CND_goodN_Step2_epCDn;
    TH1D *h_Edep_CND_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step2_epCDn;

    TH1D *h_Edep_CND_goodN_Step2_epFDn;
    TH1D *h_Edep_CND_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step2_epFDn;

    TH1D *h_Edep_CTOF_goodN_Step2_epCDn;
    TH1D *h_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step2_epCDn;

    TH1D *h_Edep_CTOF_goodN_Step2_epFDn;
    TH1D *h_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step2_epFDn;

    TH1D *h_Edep_single_goodN_Step2_epCDn;
    TH1D *h_Edep_single_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_path_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_path_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_single_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_single_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_single_badN_Step2_epCDn;

    TH1D *h_Edep_single_goodN_Step2_epFDn;
    TH1D *h_Edep_single_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_path_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_path_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_single_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_single_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_single_badN_Step2_epFDn;

    TH1D *h_Edep_CND1_goodN_Step2_epCDn;
    TH1D *h_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step2_epCDn;

    TH1D *h_Edep_CND1_goodN_Step2_epFDn;
    TH1D *h_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step2_epFDn;

    TH1D *h_Edep_CND2_goodN_Step2_epCDn;
    TH1D *h_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CND2_badN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step2_epCDn;

    TH1D *h_Edep_CND2_goodN_Step2_epFDn;
    TH1D *h_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CND2_badN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step2_epFDn;

    TH1D *h_Edep_CND3_goodN_Step2_epCDn;
    TH1D *h_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_path_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Edep_CND3_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Edep_CND3_badN_Step2_epCDn;

    TH1D *h_Edep_CND3_goodN_Step2_epFDn;
    TH1D *h_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_path_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Edep_CND3_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Edep_CND3_badN_Step2_epFDn;

    TH1D *h_Size_CND1_goodN_Step2_epCDn;
    TH1D *h_Size_CND1_badN_Step2_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_path_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_path_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Size_CND1_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Size_CND1_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Size_CND1_badN_Step2_epCDn;

    TH1D *h_Size_CND1_goodN_Step2_epFDn;
    TH1D *h_Size_CND1_badN_Step2_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_path_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_path_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Size_CND1_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Size_CND1_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Size_CND1_badN_Step2_epFDn;

    TH1D *h_Size_CND2_goodN_Step2_epCDn;
    TH1D *h_Size_CND2_badN_Step2_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_path_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_path_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Size_CND2_badN_Step2_epCDn;

    TH1D *h_Size_CND2_goodN_Step2_epFDn;
    TH1D *h_Size_CND2_badN_Step2_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_path_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_path_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Size_CND2_badN_Step2_epFDn;

    TH1D *h_Size_CND3_goodN_Step2_epCDn;
    TH1D *h_Size_CND3_badN_Step2_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_Edep_CND_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_P_n_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_dpp_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_E_p_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_path_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_path_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_ToF_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_ToF_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_nSector_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_Size_CND3_badN_Step2_epCDn;

    TH1D *h_Size_CND3_goodN_Step2_epFDn;
    TH1D *h_Size_CND3_badN_Step2_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_Edep_CND_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_Edep_CND1_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_Edep_CND2_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_Edep_CND3_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_P_n_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_dpp_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_E_p_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_path_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_path_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_ToF_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_ToF_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_nSector_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_Size_CND3_badN_Step2_epFDn;

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step2_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step2_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step2_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step2_epCDn;
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step2_epCDn;

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step2_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step2_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step2_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step2_epFDn;
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step2_epFDn;

    TH1D *h_ToF_goodN_Step2_epCDn;
    TH1D *h_ToF_badN_Step2_epCDn;
    TH2D *h_P_n_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_P_n_VS_ToF_badN_Step2_epCDn;
    TH2D *h_theta_n_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_theta_n_VS_ToF_badN_Step2_epCDn;
    TH2D *h_phi_n_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_phi_n_VS_ToF_badN_Step2_epCDn;
    TH2D *h_P_miss_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_P_miss_VS_ToF_badN_Step2_epCDn;
    TH2D *h_theta_miss_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_theta_miss_VS_ToF_badN_Step2_epCDn;
    TH2D *h_phi_miss_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_phi_miss_VS_ToF_badN_Step2_epCDn;
    TH2D *h_dpp_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_dpp_VS_ToF_badN_Step2_epCDn;
    TH2D *h_beta_n_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_beta_n_VS_ToF_badN_Step2_epCDn;
    TH2D *h_E_p_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_E_p_VS_ToF_badN_Step2_epCDn;
    TH2D *h_E_miss_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_E_miss_VS_ToF_badN_Step2_epCDn;
    TH2D *h_M_miss_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_M_miss_VS_ToF_badN_Step2_epCDn;
    TH2D *h_path_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_path_VS_ToF_badN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_theta_n_miss_VS_ToF_badN_Step2_epCDn;
    TH2D *h_nSector_VS_ToF_goodN_Step2_epCDn;
    TH2D *h_nSector_VS_ToF_badN_Step2_epCDn;

    TH1D *h_ToF_goodN_Step2_epFDn;
    TH1D *h_ToF_badN_Step2_epFDn;
    TH2D *h_P_n_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_P_n_VS_ToF_badN_Step2_epFDn;
    TH2D *h_theta_n_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_theta_n_VS_ToF_badN_Step2_epFDn;
    TH2D *h_phi_n_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_phi_n_VS_ToF_badN_Step2_epFDn;
    TH2D *h_P_miss_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_P_miss_VS_ToF_badN_Step2_epFDn;
    TH2D *h_theta_miss_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_theta_miss_VS_ToF_badN_Step2_epFDn;
    TH2D *h_phi_miss_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_phi_miss_VS_ToF_badN_Step2_epFDn;
    TH2D *h_dpp_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_dpp_VS_ToF_badN_Step2_epFDn;
    TH2D *h_beta_n_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_beta_n_VS_ToF_badN_Step2_epFDn;
    TH2D *h_E_p_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_E_p_VS_ToF_badN_Step2_epFDn;
    TH2D *h_E_miss_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_E_miss_VS_ToF_badN_Step2_epFDn;
    TH2D *h_M_miss_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_M_miss_VS_ToF_badN_Step2_epFDn;
    TH2D *h_path_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_path_VS_ToF_badN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_theta_n_miss_VS_ToF_badN_Step2_epFDn;
    TH2D *h_nSector_VS_ToF_goodN_Step2_epFDn;
    TH2D *h_nSector_VS_ToF_badN_Step2_epFDn;

    TH1D *h_beta_n_goodN_Step2_epCDn;
    TH1D *h_beta_n_badN_Step2_epCDn;

    TH1D *h_beta_n_goodN_Step2_epFDn;
    TH1D *h_beta_n_badN_Step2_epFDn;

    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step2_epCDn;
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step2_epCDn;

    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step2_epFDn;
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step2_epFDn;

    TH1D *h_Edep_CND_goodN_withNearbyPos_Step2_epCDn;
    TH1D *h_Edep_CND_badN_withNearbyPos_Step2_epCDn;

    TH1D *h_Edep_CND_goodN_withNearbyPos_Step2_epFDn;
    TH1D *h_Edep_CND_badN_withNearbyPos_Step2_epFDn;

    TH1D *h_sdiff_pos_goodN_Step2_layer_epCDn[7];
    TH1D *h_sdiff_pos_badN_Step2_layer_epCDn[7];

    TH1D *h_sdiff_pos_goodN_Step2_layer_epFDn[7];
    TH1D *h_sdiff_pos_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_mom_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_mom_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_mom_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_mom_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_VhitZ_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_VhitZ_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_VhitZ_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_VhitZ_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2_layer_epFDn[7];

    TH1D *h_theta_n_goodN_Step2_layer_epCDn[7];
    TH1D *h_theta_n_badN_Step2_layer_epCDn[7];

    TH1D *h_theta_n_goodN_Step2_layer_epFDn[7];
    TH1D *h_theta_n_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_badN_Step2_layer_epFDn[7];

    TH1D *h_phi_n_goodN_Step2_layer_epCDn[7];
    TH1D *h_phi_n_badN_Step2_layer_epCDn[7];

    TH1D *h_phi_n_goodN_Step2_layer_epFDn[7];
    TH1D *h_phi_n_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_phi_n_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_phi_n_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_phi_n_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_phi_n_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_ToF_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_ToF_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_ToF_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_ToF_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_path_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_path_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_path_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_path_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_beta_n_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_beta_n_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_beta_n_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_beta_n_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_Edep_CND_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_Edep_CND_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_Edep_CND_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_Edep_CND_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_miss_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_miss_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_miss_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_miss_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_dpp_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_dpp_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_dpp_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_dpp_badN_Step2_layer_epFDn[7];

    TH1D *h_dToF_goodN_Step2_layer_epCDn[7];
    TH1D *h_dToF_badN_Step2_layer_epCDn[7];

    TH1D *h_dToF_goodN_Step2_layer_epFDn[7];
    TH1D *h_dToF_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_dToF_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_dToF_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_dToF_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_dToF_badN_Step2_layer_epFDn[7];

    TH1D *h_dToF_rel_pos_goodN_Step2_layer_epCDn[7];
    TH1D *h_dToF_rel_pos_badN_Step2_layer_epCDn[7];

    TH1D *h_dToF_rel_pos_goodN_Step2_layer_epFDn[7];
    TH1D *h_dToF_rel_pos_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_dToF_rel_pos_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_dToF_rel_pos_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_dToF_rel_pos_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_dToF_rel_pos_badN_Step2_layer_epFDn[7];

    TH1D *h_dToF_rel_n_goodN_Step2_layer_epCDn[7];
    TH1D *h_dToF_rel_n_badN_Step2_layer_epCDn[7];

    TH1D *h_dToF_rel_n_goodN_Step2_layer_epFDn[7];
    TH1D *h_dToF_rel_n_badN_Step2_layer_epFDn[7];

    TH2D *h_sdiff_pos_VS_dToF_rel_n_goodN_Step2_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_dToF_rel_n_badN_Step2_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_dToF_rel_n_goodN_Step2_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_dToF_rel_n_badN_Step2_layer_epFDn[7];

    TH2D *h_diff_ToFc_z_VS_Edep_noNear_goodN_Step2_epCDn;
    TH2D *h_diff_ToFc_z_VS_Edep_noNear_badN_Step2_epCDn;

    TH2D *h_diff_ToFc_z_VS_Edep_noNear_goodN_Step2_epFDn;
    TH2D *h_diff_ToFc_z_VS_Edep_noNear_badN_Step2_epFDn;

    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step2_epCDn;
    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_badN_Step2_epCDn;

    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step2_epFDn;
    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_badN_Step2_epFDn;

    TH1D *h_numberNearby_goodN_Step2_epCDn;
    TH1D *h_numberNearby_badN_Step2_epCDn;

    TH1D *h_numberNearby_goodN_Step2_epFDn;
    TH1D *h_numberNearby_badN_Step2_epFDn;

    TH2D *h_numberNearby_momN_goodN_Step2_epCDn;
    TH2D *h_numberNearby_momN_badN_Step2_epCDn;

    TH2D *h_numberNearby_momN_goodN_Step2_epFDn;
    TH2D *h_numberNearby_momN_badN_Step2_epFDn;

    TH1D *h_NearbyEdep_goodN_Step2_epCDn;
    TH1D *h_NearbyEdep_badN_Step2_epCDn;

    TH1D *h_NearbyEdep_goodN_Step2_epFDn;
    TH1D *h_NearbyEdep_badN_Step2_epFDn;

    TH1D *h_nsector_goodN_Step2_epCDn;
    TH1D *h_nsector_badN_Step2_epCDn;

    TH1D *h_nsector_goodN_Step2_epFDn;
    TH1D *h_nsector_badN_Step2_epFDn;

#pragma endregion /* Step Two (After applying Phi Diff Charge Track cut) (Andrew) - end */

    // Step Three (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step3;
    TH1D *h_n_multiplicity_goodN_epCDn_Step3;
    TH1D *h_n_multiplicity_badN_epCDn_Step3;

    TH1D *h_n_multiplicity_allN_epFDn_Step3;
    TH1D *h_n_multiplicity_goodN_epFDn_Step3;
    TH1D *h_n_multiplicity_badN_epFDn_Step3;

    // Step Four (After applying Phi Diff CND hit cut) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step4;
    TH1D *h_n_multiplicity_goodN_epCDn_Step4;
    TH1D *h_n_multiplicity_badN_epCDn_Step4;

    TH1D *h_n_multiplicity_allN_epFDn_Step4;
    TH1D *h_n_multiplicity_goodN_epFDn_Step4;
    TH1D *h_n_multiplicity_badN_epFDn_Step4;

    // Step Five (After event selection cuts) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step5;
    TH1D *h_n_multiplicity_goodN_epCDn_Step5;
    TH1D *h_n_multiplicity_badN_epCDn_Step5;

    TH1D *h_n_multiplicity_allN_epFDn_Step5;
    TH1D *h_n_multiplicity_goodN_epFDn_Step5;
    TH1D *h_n_multiplicity_badN_epFDn_Step5;

    // Constructor
    // ======================================================================================================================================================================

    VetoHistograms();

    // InitHistograms function
    // ======================================================================================================================================================================

    void InitHistograms();

    // UpdateBPIDpCDHistograms function
    // ======================================================================================================================================================================

    void UpdateBPIDpCDHistograms(bool pInCD, bool pInFD, TVector3 P_miss_3v, double E_p, double E_miss, double M_miss, double xB, double weight);

    // UpdateAPIDpCDHistograms function
    // ======================================================================================================================================================================

    void UpdateAPIDpCDHistograms(bool pInCD, bool pInFD, TVector3 P_miss_3v, double E_p, double E_miss, double M_miss, double xB, double weight);

    // UpdateBPIDpFDHistograms function
    // ======================================================================================================================================================================

    void UpdateBPIDpFDHistograms(bool pInFD, bool pInFD, TVector3 P_miss_3v, double E_p, double E_miss, double M_miss, double xB, double weight);

    // UpdateAPIDpFDHistograms function
    // ======================================================================================================================================================================

    void UpdateAPIDpFDHistograms(bool pInFD, bool pInFD, TVector3 P_miss_3v, double E_p, double E_miss, double M_miss, double xB, double weight);

    // UpdateProtonMultiBCHistograms function
    // ======================================================================================================================================================================

    void UpdateProtonMultiBCHistograms(int counter_pCD_multiplicity_BPID, int counter_pFD_multiplicity_BPID, double weight);

    // UpdateProtonMultiACHistograms function
    // ======================================================================================================================================================================

    void UpdateProtonMultiACHistograms(int counter_pCD_multiplicity_APID, int counter_pFD_multiplicity_APID, double weight);

    // UpdateBmissCHistograms function
    // ======================================================================================================================================================================

    void UpdateBmissCHistograms(bool pInCD, bool pInFD, TVector3 P_miss_3v, double E_p, double E_miss, double M_miss, double xB, double weight);

    // UpdateAmissCHistograms function
    // ======================================================================================================================================================================

    void UpdateAmissCHistograms(bool pInCD, bool pInFD, TVector3 P_miss_3v, double E_p, double E_miss, double M_miss, double xB, double weight);

    // UpdatePreStepHistograms function
    // ======================================================================================================================================================================

    void UpdatePreStepHistograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 P_p_3v, TVector3 P_miss_3v, TVector3 P_n_3v,
                                 double E_p, double E_miss, double M_miss, double xB, double dpp, double theta_n_miss, double Edep_CND,
                                 double Edep_CND1, double Edep_CND2, double Edep_CND3, double Edep_CTOF, double nSector, double Size_CND1,
                                 double Size_CND2, double Size_CND3, double LayerMult_CND1, double LayerMult_CND2, double LayerMult_CND3, double beta,
                                 double path, double ToF, double weight);

    // UpdateBS0CHistograms function
    // ======================================================================================================================================================================

    void UpdateBS0CHistograms(bool pInCD, bool pInFD, TVector3 P_n_3v, TVector3 v_hit_3v, double beta, double path, double ToF,
                              double weight);

    // UpdateAS0CHistograms function
    // ======================================================================================================================================================================

    void UpdateAS0CHistograms(bool pInCD, bool pInFD, TVector3 P_n_3v, TVector3 v_hit_3v, double beta, double path, double ToF,
                              double weight);

    // UpdateStep0Histograms function
    // ======================================================================================================================================================================

    void UpdateStep0Histograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 P_p_3v, TVector3 P_miss_3v, TVector3 P_n_3v,
                               double E_p, double E_miss, double M_miss, double xB, double dpp, double theta_n_miss, double Edep_CND,
                               double Edep_CND1, double Edep_CND2, double Edep_CND3, double Edep_CTOF, double nSector, double Size_CND1,
                               double Size_CND2, double Size_CND3, double LayerMult_CND1, double LayerMult_CND2, double LayerMult_CND3, double beta,
                               double path, double ToF, double weight);

    // UpdateStep1Histograms function
    // ======================================================================================================================================================================

    void UpdateStep1Histograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 P_p_3v, TVector3 P_miss_3v, TVector3 P_n_3v,
                               double E_p, double E_miss, double M_miss, double xB, double dpp, double theta_n_miss, double Edep_CND,
                               double Edep_CND1, double Edep_CND2, double Edep_CND3, double Edep_CTOF, double nSector, double Size_CND1,
                               double Size_CND2, double Size_CND3, double LayerMult_CND1, double LayerMult_CND2, double LayerMult_CND3, double beta,
                               double path, double ToF, double weight);

    // UpdateStep2prepBCHistograms function
    // ======================================================================================================================================================================

    void UpdateStep2prepBCHistograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 v_hit_3v, double ToF, double weight);

    // UpdateStep2prepHistograms function
    // ======================================================================================================================================================================

    void UpdateStep2prepHistograms(bool pInCD, bool pInFD, bool isGN, bool isBN, int ldiff, int sdiff,
                                   TVector3 p_C_3v, TVector3 v_hit_3v, TVector3 P_n_3v, double dToF, double dToF_rel_pos,
                                   double dToF_rel_n, double dpp, double theta_n_miss, double Edep_CND, double beta, double path,
                                   double ToF, double weight);

    // UpdateMonitorStep2prepHistograms1 function
    // ======================================================================================================================================================================

    void UpdateMonitorStep2prepHistograms1(bool Nearby_clusters_from_cPart_tracks, bool pInCD, bool pInFD, bool isGN, bool isBN, double Edep_CND,
                                           double Edep_CTOF_pos, double weight);

    // UpdateMonitorStep2prepHistograms2 function
    // ======================================================================================================================================================================

    void UpdateMonitorStep2prepHistograms2(bool Nearby_clusters_from_cPart_tracks, bool pInCD, bool pInFD, bool isGN, bool isBN, double Edep_CND,
                                           double ToF, TVector3 v_hit_3v, double weight);

    // UpdateBS0CHistograms function
    // ======================================================================================================================================================================

    void UpdateBS2CHistograms(bool pInCD, bool pInFD, double Size_CND1, double Size_CND2, double Size_CND3, double LayerMult_CND1,
                              double LayerMult_CND2, double LayerMult_CND3, double weight);

    // UpdateAS0CHistograms function
    // ======================================================================================================================================================================

    void UpdateAS2CHistograms(bool pInCD, bool pInFD, double Size_CND1, double Size_CND2, double Size_CND3, double LayerMult_CND1,
                              double LayerMult_CND2, double LayerMult_CND3, double weight);

    // UpdateStep2Histograms function
    // ======================================================================================================================================================================

    void UpdateStep2Histograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 P_p_3v, TVector3 P_miss_3v, TVector3 P_n_3v,
                               double E_p, double E_miss, double M_miss, double xB, double dpp, double theta_n_miss, double Edep_CND,
                               double Edep_CND1, double Edep_CND2, double Edep_CND3, double Edep_CTOF, double nSector, double Size_CND1,
                               double Size_CND2, double Size_CND3, double LayerMult_CND1, double LayerMult_CND2, double LayerMult_CND3, double beta,
                               double path, double ToF, double weight);

    // UpdateStep2Histograms2 function
    // ======================================================================================================================================================================

    void UpdateStep2Histograms2(bool pInCD, bool pInFD, bool isGN, bool isBN, int ldiff, int sdiff,
                                TVector3 p_C_3v, TVector3 v_hit_3v, TVector3 P_n_3v, double dToF, double dToF_rel_pos,
                                double dToF_rel_n, double dpp, double theta_n_miss, double Edep_CND, double beta, double path,
                                double ToF, double weight);

    // UpdateMultiplicityHistograms function
    // ======================================================================================================================================================================

    void UpdateMultiplicityHistograms(bool pInCD, bool pInFD, bool isGN, bool isBN,
                                      int counter_n_multiplicity_allN_epCDn_Step0, int counter_n_multiplicity_goodN_epCDn_Step0,
                                      int counter_n_multiplicity_badN_epCDn_Step0,
                                      int counter_n_multiplicity_allN_epCDn_Step1, int counter_n_multiplicity_goodN_epCDn_Step1,
                                      int counter_n_multiplicity_badN_epCDn_Step1,
                                      int counter_n_multiplicity_allN_epCDn_Step2, int counter_n_multiplicity_goodN_epCDn_Step2,
                                      int counter_n_multiplicity_badN_epCDn_Step2,
                                      int counter_n_multiplicity_allN_epCDn_Step3, int counter_n_multiplicity_goodN_epCDn_Step3,
                                      int counter_n_multiplicity_badN_epCDn_Step3,
                                      int counter_n_multiplicity_allN_epCDn_Step4, int counter_n_multiplicity_goodN_epCDn_Step4,
                                      int counter_n_multiplicity_badN_epCDn_Step4,
                                      int counter_n_multiplicity_allN_epCDn_Step5, int counter_n_multiplicity_goodN_epCDn_Step5,
                                      int counter_n_multiplicity_badN_epCDn_Step5,
                                      int counter_n_multiplicity_allN_epFDn_Step0, int counter_n_multiplicity_goodN_epFDn_Step0,
                                      int counter_n_multiplicity_badN_epFDn_Step0,
                                      int counter_n_multiplicity_allN_epFDn_Step1, int counter_n_multiplicity_goodN_epFDn_Step1,
                                      int counter_n_multiplicity_badN_epFDn_Step1,
                                      int counter_n_multiplicity_allN_epFDn_Step2, int counter_n_multiplicity_goodN_epFDn_Step2,
                                      int counter_n_multiplicity_badN_epFDn_Step2,
                                      int counter_n_multiplicity_allN_epFDn_Step3, int counter_n_multiplicity_goodN_epFDn_Step3,
                                      int counter_n_multiplicity_badN_epFDn_Step3,
                                      int counter_n_multiplicity_allN_epFDn_Step4, int counter_n_multiplicity_goodN_epFDn_Step4,
                                      int counter_n_multiplicity_badN_epFDn_Step4,
                                      int counter_n_multiplicity_allN_epFDn_Step5, int counter_n_multiplicity_goodN_epFDn_Step5,
                                      int counter_n_multiplicity_badN_epFDn_Step5,
                                      double weight);

    // GetHistogramEntries function
    // ======================================================================================================================================================================

    double GetHistogramEntries(const std::vector<TH1 *> &HistoList, const std::string &histName);

    // extractStep function
    // ======================================================================================================================================================================

    std::string extractStep(const std::string &input);

    // SkippingCondition function
    // ======================================================================================================================================================================

    bool SkippingCondition(string HistoName, int canvas_ind);

    // replaceSubstring function
    // ======================================================================================================================================================================

    std::string replaceSubstring(const std::string &input, const std::string &toReplace, const std::string &replaceWith);

    // SectionPlotter function
    // ======================================================================================================================================================================

    void SectionPlotter(int n_col, int n_row, TCanvas *myCanvas, TCanvas *myText, TCanvas *myTable, vector<TH1 *> HistoList,
                        string PDFFile, string Constraint1 = "", string Constraint2 = "", bool LogScale2D = false);

    // HistPrinter function
    // ======================================================================================================================================================================

    void HistPrinter(vector<TH1 *> HistoList, string PDFFile, bool LogScale2D = false);

    void PlotHistograms(string PDFFile) {
        HistPrinter(HistoList, PDFFile);
    };
};


#endif //VETOHISTOGRAMS_H
