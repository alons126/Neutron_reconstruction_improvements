//
// Created by Alon Sportes on 13/03/2025.
//

#ifndef VARIABLEHISTOGRAMS_H
#define VARIABLEHISTOGRAMS_H

#include <cstdlib>
#include <iostream>
#include <vector>
//
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TTree.h"
//
#include "../../namespaces/utilities/utilities.cpp"
//
using namespace std;
using namespace utilities;

class VariableHistograms {
   private:
    bool OnlySeparatedN = false;

    TH1D *h_Var_1D;
    TH2D *h_P_n_VS_Var_2D;
    TH2D *h_theta_n_VS_Var_2D;
    TH2D *h_phi_n_VS_Var_2D;
    TH2D *h_P_miss_VS_Var_2D;
    TH2D *h_theta_miss_VS_Var_2D;
    TH2D *h_phi_miss_VS_Var_2D;
    TH2D *h_dpp_VS_Var_2D;
    TH2D *h_theta_n_miss_VS_Var_2D;
    TH2D *h_beta_n_VS_Var_2D;
    TH2D *h_E_p_VS_Var_2D;
    TH2D *h_E_miss_VS_Var_2D;
    TH2D *h_path_VS_Var_2D;
    TH2D *h_ToF_VS_Var_2D;
    TH2D *h_nSector_VS_Var_2D;
    TH2D *h_Edep_CND_VS_Var_2D;
    TH2D *h_Edep_CND1_VS_Var_2D;
    TH2D *h_Edep_CND2_VS_Var_2D;
    TH2D *h_Edep_CND3_VS_Var_2D;
    TH2D *h_Size_CND1_VS_Var_2D;
    TH2D *h_Size_CND2_VS_Var_2D;
    TH2D *h_Size_CND3_VS_Var_2D;
    TH2D *h_LayerMulti_CND1_VS_Var_2D;
    TH2D *h_LayerMulti_CND2_VS_Var_2D;
    TH2D *h_LayerMulti_CND3_VS_Var_2D;

    TH1D *h_Var_goodN_1D;
    TH1D *h_Var_badN_1D;
    TH2D *h_P_n_VS_Var_goodN_2D;
    TH2D *h_P_n_VS_Var_badN_2D;
    TH2D *h_theta_n_VS_Var_goodN_2D;
    TH2D *h_theta_n_VS_Var_badN_2D;
    TH2D *h_phi_n_VS_Var_goodN_2D;
    TH2D *h_phi_n_VS_Var_badN_2D;
    TH2D *h_P_miss_VS_Var_goodN_2D;
    TH2D *h_P_miss_VS_Var_badN_2D;
    TH2D *h_theta_miss_VS_Var_goodN_2D;
    TH2D *h_theta_miss_VS_Var_badN_2D;
    TH2D *h_phi_miss_VS_Var_goodN_2D;
    TH2D *h_phi_miss_VS_Var_badN_2D;
    TH2D *h_dpp_VS_Var_goodN_2D;
    TH2D *h_dpp_VS_Var_badN_2D;
    TH2D *h_theta_n_miss_VS_Var_goodN_2D;
    TH2D *h_theta_n_miss_VS_Var_badN_2D;
    TH2D *h_beta_n_VS_Var_goodN_2D;
    TH2D *h_beta_n_VS_Var_badN_2D;
    TH2D *h_E_p_VS_Var_goodN_2D;
    TH2D *h_E_p_VS_Var_badN_2D;
    TH2D *h_E_miss_VS_Var_goodN_2D;
    TH2D *h_E_miss_VS_Var_badN_2D;
    TH2D *h_path_VS_Var_goodN_2D;
    TH2D *h_path_VS_Var_badN_2D;
    TH2D *h_ToF_VS_Var_goodN_2D;
    TH2D *h_ToF_VS_Var_badN_2D;
    TH2D *h_nSector_VS_Var_goodN_2D;
    TH2D *h_nSector_VS_Var_badN_2D;
    TH2D *h_Edep_CND_VS_Var_goodN_2D;
    TH2D *h_Edep_CND_VS_Var_badN_2D;
    TH2D *h_Edep_CND1_VS_Var_goodN_2D;
    TH2D *h_Edep_CND1_VS_Var_badN_2D;
    TH2D *h_Edep_CND2_VS_Var_goodN_2D;
    TH2D *h_Edep_CND2_VS_Var_badN_2D;
    TH2D *h_Edep_CND3_VS_Var_goodN_2D;
    TH2D *h_Edep_CND3_VS_Var_badN_2D;
    TH2D *h_Size_CND1_VS_Var_goodN_2D;
    TH2D *h_Size_CND1_VS_Var_badN_2D;
    TH2D *h_Size_CND2_VS_Var_goodN_2D;
    TH2D *h_Size_CND2_VS_Var_badN_2D;
    TH2D *h_Size_CND3_VS_Var_goodN_2D;
    TH2D *h_Size_CND3_VS_Var_badN_2D;
    TH2D *h_LayerMulti_CND1_VS_Var_goodN_2D;
    TH2D *h_LayerMulti_CND1_VS_Var_badN_2D;
    TH2D *h_LayerMulti_CND2_VS_Var_goodN_2D;
    TH2D *h_LayerMulti_CND2_VS_Var_badN_2D;
    TH2D *h_LayerMulti_CND3_VS_Var_goodN_2D;
    TH2D *h_LayerMulti_CND3_VS_Var_badN_2D;

   public:
    // Constructor -----------------------------------------------------------------------------------------------------------------------------------------------------------

    VariableHistograms() = default;

    // Constructor -----------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitHistograms(vector<vector<TH1 *>> &HistoList, const string VarName, const string Step, const string FinalState, const string TitleAdditions = "", const bool separateGoodBadN = false);
};

#endif  // VARIABLEHISTOGRAMS_H
