//
// Created by Alon Sportes on 20/01/2025.
//

#include "VetoHistograms.cpp"

// UpdateBmissCHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateBmissCHistograms(bool pInCD, bool pInFD, TVector3 P_miss_3v, double E_p, double E_miss, double M_miss, double xB,
                                            double weight) {
    if (pInCD) {
        h_P_miss_BmissC_epCD->Fill(P_miss_3v.Mag(), weight);
        h_theta_miss_BmissC_epCD->Fill(P_miss_3v.Theta() * 180 / M_PI, weight);
        h_P_miss_VS_theta_miss_BmissC_epCD->Fill(P_miss_3v.Theta() * 180 / M_PI, P_miss_3v.Mag(), weight);

        // h_beta_n_BmissC_epCD->Fill(AllParticles[itr1]->par()->getBeta(), weight);

        h_E_p_BmissC_epCD->Fill(E_p, weight);
        h_E_miss_BmissC_epCD->Fill(E_miss, weight);
        h_M_miss_BmissC_epCD->Fill(M_miss, weight);

        h_xB_BmissC_epCD->Fill(xB, weight);
        h_xB_VS_M_miss_BmissC_epCD->Fill(xB, M_miss, weight);
    } else if (pInFD) {
        h_P_miss_BmissC_epFD->Fill(P_miss_3v.Mag(), weight);
        h_theta_miss_BmissC_epFD->Fill(P_miss_3v.Theta() * 180 / M_PI, weight);
        h_P_miss_VS_theta_miss_BmissC_epFD->Fill(P_miss_3v.Theta() * 180 / M_PI, P_miss_3v.Mag(), weight);

        // h_beta_n_BmissC_epFD->Fill(AllParticles[itr1]->par()->getBeta(), weight);

        h_E_p_BmissC_epFD->Fill(E_p, weight);
        h_E_miss_BmissC_epFD->Fill(E_miss, weight);
        h_M_miss_BmissC_epFD->Fill(M_miss, weight);

        h_xB_BmissC_epFD->Fill(xB, weight);
        h_xB_VS_M_miss_BmissC_epFD->Fill(xB, M_miss, weight);
    }
}

// UpdateAmissCHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateAmissCHistograms(bool pInCD, bool pInFD, TVector3 P_miss_3v, double E_p, double E_miss, double M_miss, double xB,
                                            double weight) {
    if (pInCD) {
        h_P_miss_AmissC_epCD->Fill(P_miss_3v.Mag(), weight);
        h_theta_miss_AmissC_epCD->Fill(P_miss_3v.Theta() * 180 / M_PI, weight);
        h_P_miss_VS_theta_miss_AmissC_epCD->Fill(P_miss_3v.Theta() * 180 / M_PI, P_miss_3v.Mag(), weight);

        // h_beta_n_AmissC_epCD->Fill(AllParticles[itr1]->par()->getBeta(), weight);

        h_E_p_AmissC_epCD->Fill(E_p, weight);
        h_E_miss_AmissC_epCD->Fill(E_miss, weight);
        h_M_miss_AmissC_epCD->Fill(M_miss, weight);

        h_xB_AmissC_epCD->Fill(xB, weight);
        h_xB_VS_M_miss_AmissC_epCD->Fill(xB, M_miss, weight);
    } else if (pInFD) {
        h_P_miss_AmissC_epFD->Fill(P_miss_3v.Mag(), weight);
        h_theta_miss_AmissC_epFD->Fill(P_miss_3v.Theta() * 180 / M_PI, weight);
        h_P_miss_VS_theta_miss_AmissC_epFD->Fill(P_miss_3v.Theta() * 180 / M_PI, P_miss_3v.Mag(), weight);

        // h_beta_n_AmissC_epFD->Fill(AllParticles[itr1]->par()->getBeta(), weight);

        h_E_p_AmissC_epFD->Fill(E_p, weight);
        h_E_miss_AmissC_epFD->Fill(E_miss, weight);
        h_M_miss_AmissC_epFD->Fill(M_miss, weight);

        h_xB_AmissC_epFD->Fill(xB, weight);
        h_xB_VS_M_miss_AmissC_epFD->Fill(xB, M_miss, weight);
    }
}

// UpdateAmissCHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdatePreStepHistograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 P_p_3v, TVector3 P_miss_3v, TVector3 P_n_3v,
                                             double E_p, double E_miss, double M_miss, double xB, double dpp, double theta_n_miss,
                                             double Edep_CND, double Edep_CND1, double Edep_CND2, double Edep_CND3, double Edep_CTOF,
                                             double nSector, double Size_CND1, double Size_CND2, double Size_CND3,
                                             double LayerMult_CND1, double LayerMult_CND2, double LayerMult_CND3,
                                             double beta, double path, double ToF, double weight) {
    if (pInCD) {
        h_xB_VS_M_miss_epCDn->Fill(xB, M_miss, weight);

        h_theta_n_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
        h_theta_n_VS_phi_n_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
        h_theta_n_VS_beta_n_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

        h_P_n_epCDn->Fill(P_n_3v.Mag(), weight);
        h_P_n_VS_theta_n_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

        h_P_miss_epCDn->Fill(P_miss_3v.Mag(), weight);
        h_P_miss_VS_theta_miss_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);

        h_dpp_allN_epCDn->Fill(dpp, weight);
        h_theta_n_miss_allN_epCDn->Fill(theta_n_miss, weight);
        h_dpp_VS_theta_n_miss_epCDn->Fill(dpp, theta_n_miss, weight);

        if (isGN) {
            h_dpp_goodN_epCDn->Fill(dpp, weight);
            h_theta_n_miss_goodN_epCDn->Fill(theta_n_miss, weight);
        } else if (isBN) {
            h_dpp_badN_epCDn->Fill(dpp, weight);
            h_theta_n_miss_badN_epCDn->Fill(theta_n_miss, weight);
        }

        if (theta_n_miss < 25.) {
            h_dpp_allN_for_theta_n_miss_less_than_25_epCDn->Fill(dpp, weight);
        }

        if (dpp < 0.5) {
            h_theta_n_miss_allN_for_dpp_less_than_05_epCDn->Fill(theta_n_miss, weight);

            if (dpp < 0.3) {
                h_theta_n_miss_allN_for_dpp_less_than_03_epCDn->Fill(theta_n_miss, weight);
            }
        }

        h_E_p_epCDn->Fill(E_p, weight);
        h_E_miss_epCDn->Fill(E_miss, weight);
        h_M_miss_epCDn->Fill(M_miss, weight);
        h_M_miss_VS_P_n_epCDn->Fill(P_n_3v.Mag(), M_miss, weight);
        h_M_miss_VS_theta_n_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
        h_M_miss_VS_phi_n_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
        h_M_miss_VS_P_miss_epCDn->Fill(P_miss_3v.Mag(), M_miss, weight);
        h_M_miss_VS_theta_miss_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
        h_M_miss_VS_phi_miss_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

        h_P_n_minus_P_miss_epCDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
        h_P_n_x_minus_P_miss_x_epCDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
        h_P_n_y_minus_P_miss_y_epCDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
        h_P_n_z_minus_P_miss_z_epCDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

        h_P_n_VS_P_miss_epCDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
        h_P_n_x_VS_P_miss_x_epCDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
        h_P_n_y_VS_P_miss_y_epCDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
        h_P_n_z_VS_P_miss_z_epCDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

        h_theta_n_p_epCDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
        h_theta_p_n_VS_P_p_epCDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

        h_xB_epCDn->Fill(xB, weight);

        h_Edep_CND_epCDn->Fill(Edep_CND, weight);
        h_P_n_VS_Edep_CND_epCDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CND_epCDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CND_epCDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CND_epCDn->Fill(Edep_CND, dpp, weight);
        h_beta_n_VS_Edep_CND_epCDn->Fill(Edep_CND, beta, weight);
        h_E_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, E_miss, weight);
        h_M_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, M_miss, weight);
        h_path_VS_Edep_CND_epCDn->Fill(Edep_CND, path * 100, weight);
        h_theta_n_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, theta_n_miss, weight);
        h_ToF_VS_Edep_CND_epCDn->Fill(Edep_CND, ToF, weight);
        h_nSector_VS_Edep_CND_epCDn->Fill(Edep_CND, nSector, weight);
        h_Edep_CND1_VS_Edep_CND_epCDn->Fill(Edep_CND, Edep_CND1, weight);
        h_Edep_CND2_VS_Edep_CND_epCDn->Fill(Edep_CND, Edep_CND2, weight);
        h_Edep_CND3_VS_Edep_CND_epCDn->Fill(Edep_CND, Edep_CND3, weight);

        h_Edep_CTOF_epCDn->Fill(Edep_CTOF, weight);
        h_P_n_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, dpp, weight);
        h_beta_n_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, beta, weight);
        h_E_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, E_miss, weight);
        h_M_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, M_miss, weight);
        h_path_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, path * 100, weight);
        h_theta_n_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
        h_ToF_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, ToF, weight);
        h_nSector_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, nSector, weight);
        h_Edep_CND1_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, Edep_CND1, weight);
        h_Edep_CND2_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, Edep_CND2, weight);
        h_Edep_CND3_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, Edep_CND3, weight);

        // h_Edep_single_epCDn->Fill(Edep_single, weight);
        // h_P_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n_3v.Mag(), weight);
        // h_theta_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n_3v.Theta() * 180. / M_PI, weight);
        // h_phi_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n_3v.Phi() * 180. / M_PI, weight);
        // h_P_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss_3v.Mag(), weight);
        // h_theta_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss_3v.Theta() * 180. / M_PI, weight);
        // h_phi_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss_3v.Phi() * 180. / M_PI, weight);
        // h_dpp_VS_Edep_single_epCDn->Fill(Edep_single, dpp, weight);
        // h_beta_n_VS_Edep_single_epCDn->Fill(Edep_single, beta, weight);
        // h_E_miss_VS_Edep_single_epCDn->Fill(Edep_single, E_miss, weight);
        // h_M_miss_VS_Edep_single_epCDn->Fill(Edep_single, M_miss, weight);
        // h_path_VS_Edep_single_epCDn->Fill(Edep_single, path * 100, weight);
        // h_theta_n_miss_VS_Edep_single_epCDn->Fill(Edep_single, theta_n_miss, weight);
        // h_ToF_VS_Edep_single_epCDn->Fill(Edep_single, ToF, weight);
        // h_nSector_VS_Edep_single_epCDn->Fill(Edep_single, nSector, weight);

        h_Edep_CND1_epCDn->Fill(Edep_CND1, weight);
        h_P_n_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CND1_epCDn->Fill(Edep_CND1, dpp, weight);
        h_beta_n_VS_Edep_CND1_epCDn->Fill(Edep_CND1, beta, weight);
        h_E_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, E_miss, weight);
        h_M_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, M_miss, weight);
        h_path_VS_Edep_CND1_epCDn->Fill(Edep_CND1, path * 100, weight);
        h_theta_n_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
        h_ToF_VS_Edep_CND1_epCDn->Fill(Edep_CND1, ToF, weight);
        h_nSector_VS_Edep_CND1_epCDn->Fill(Edep_CND1, nSector, weight);
        h_Edep_CND2_VS_Edep_CND1_epCDn->Fill(Edep_CND1, Edep_CND2, weight);
        h_Edep_CND3_VS_Edep_CND1_epCDn->Fill(Edep_CND1, Edep_CND3, weight);

        h_Edep_CND2_epCDn->Fill(Edep_CND2, weight);
        h_P_n_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CND2_epCDn->Fill(Edep_CND2, dpp, weight);
        h_beta_n_VS_Edep_CND2_epCDn->Fill(Edep_CND2, beta, weight);
        h_E_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, E_miss, weight);
        h_M_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, M_miss, weight);
        h_path_VS_Edep_CND2_epCDn->Fill(Edep_CND2, path * 100, weight);
        h_theta_n_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
        h_ToF_VS_Edep_CND2_epCDn->Fill(Edep_CND2, ToF, weight);
        h_nSector_VS_Edep_CND2_epCDn->Fill(Edep_CND2, nSector, weight);
        h_Edep_CND3_VS_Edep_CND2_epCDn->Fill(Edep_CND2, Edep_CND3, weight);

        h_Edep_CND3_epCDn->Fill(Edep_CND3, weight);
        h_P_n_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CND3_epCDn->Fill(Edep_CND3, dpp, weight);
        h_beta_n_VS_Edep_CND3_epCDn->Fill(Edep_CND3, beta, weight);
        h_E_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, E_miss, weight);
        h_M_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, M_miss, weight);
        h_path_VS_Edep_CND3_epCDn->Fill(Edep_CND3, path * 100, weight);
        h_theta_n_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
        h_ToF_VS_Edep_CND3_epCDn->Fill(Edep_CND3, ToF, weight);
        h_nSector_VS_Edep_CND3_epCDn->Fill(Edep_CND3, nSector, weight);

        h_Size_CND1_epCDn->Fill(Size_CND1, weight);
        h_Edep_CND_VS_Size_CND1_epCDn->Fill(Size_CND1, Edep_CND, weight);
        h_Edep_CND1_VS_Size_CND1_epCDn->Fill(Size_CND1, Edep_CND1, weight);
        h_Edep_CND2_VS_Size_CND1_epCDn->Fill(Size_CND1, Edep_CND2, weight);
        h_Edep_CND3_VS_Size_CND1_epCDn->Fill(Size_CND1, Edep_CND3, weight);
        h_P_n_VS_Size_CND1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
        h_theta_n_VS_Size_CND1_epCDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Size_CND1_epCDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Size_CND1_epCDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Size_CND1_epCDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Size_CND1_epCDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Size_CND1_epCDn->Fill(Size_CND1, dpp, weight);
        h_beta_n_VS_Size_CND1_epCDn->Fill(Size_CND1, beta, weight);
        h_E_miss_VS_Size_CND1_epCDn->Fill(Size_CND1, E_miss, weight);
        h_M_miss_VS_Size_CND1_epCDn->Fill(Size_CND1, M_miss, weight);
        h_path_VS_Size_CND1_epCDn->Fill(Size_CND1, path * 100, weight);
        h_theta_n_miss_VS_Size_CND1_epCDn->Fill(Size_CND1, theta_n_miss, weight);
        h_ToF_VS_Size_CND1_epCDn->Fill(Size_CND1, ToF, weight);
        h_nSector_VS_Size_CND1_epCDn->Fill(Size_CND1, nSector, weight);

        h_Size_CND2_epCDn->Fill(Size_CND2, weight);
        h_Edep_CND_VS_Size_CND2_epCDn->Fill(Size_CND2, Edep_CND, weight);
        h_Edep_CND1_VS_Size_CND2_epCDn->Fill(Size_CND2, Edep_CND1, weight);
        h_Edep_CND2_VS_Size_CND2_epCDn->Fill(Size_CND2, Edep_CND2, weight);
        h_Edep_CND3_VS_Size_CND2_epCDn->Fill(Size_CND2, Edep_CND3, weight);
        h_P_n_VS_Size_CND2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
        h_theta_n_VS_Size_CND2_epCDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Size_CND2_epCDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Size_CND2_epCDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Size_CND2_epCDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Size_CND2_epCDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Size_CND2_epCDn->Fill(Size_CND2, dpp, weight);
        h_beta_n_VS_Size_CND2_epCDn->Fill(Size_CND2, beta, weight);
        h_E_miss_VS_Size_CND2_epCDn->Fill(Size_CND2, E_miss, weight);
        h_M_miss_VS_Size_CND2_epCDn->Fill(Size_CND2, M_miss, weight);
        h_path_VS_Size_CND2_epCDn->Fill(Size_CND2, path * 100, weight);
        h_theta_n_miss_VS_Size_CND2_epCDn->Fill(Size_CND2, theta_n_miss, weight);
        h_ToF_VS_Size_CND2_epCDn->Fill(Size_CND2, ToF, weight);
        h_nSector_VS_Size_CND2_epCDn->Fill(Size_CND2, nSector, weight);

        h_Size_CND3_epCDn->Fill(Size_CND3, weight);
        h_Edep_CND_VS_Size_CND3_epCDn->Fill(Size_CND3, Edep_CND, weight);
        h_Edep_CND1_VS_Size_CND3_epCDn->Fill(Size_CND3, Edep_CND1, weight);
        h_Edep_CND2_VS_Size_CND3_epCDn->Fill(Size_CND3, Edep_CND2, weight);
        h_Edep_CND3_VS_Size_CND3_epCDn->Fill(Size_CND3, Edep_CND3, weight);
        h_P_n_VS_Size_CND3_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND3_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND3_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
        h_theta_n_VS_Size_CND3_epCDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Size_CND3_epCDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Size_CND3_epCDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Size_CND3_epCDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Size_CND3_epCDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Size_CND3_epCDn->Fill(Size_CND3, dpp, weight);
        h_beta_n_VS_Size_CND3_epCDn->Fill(Size_CND3, beta, weight);
        h_E_miss_VS_Size_CND3_epCDn->Fill(Size_CND3, E_miss, weight);
        h_M_miss_VS_Size_CND3_epCDn->Fill(Size_CND3, M_miss, weight);
        h_path_VS_Size_CND3_epCDn->Fill(Size_CND3, path * 100, weight);
        h_theta_n_miss_VS_Size_CND3_epCDn->Fill(Size_CND3, theta_n_miss, weight);
        h_ToF_VS_Size_CND3_epCDn->Fill(Size_CND3, ToF, weight);
        h_nSector_VS_Size_CND3_epCDn->Fill(Size_CND3, nSector, weight);

        h_Size_CND1_VS_Size_CND2_epCDn->Fill(Size_CND1, Size_CND2, weight);
        h_Size_CND1_VS_Size_CND3_epCDn->Fill(Size_CND1, Size_CND3, weight);
        h_Size_CND2_VS_Size_CND3_epCDn->Fill(Size_CND2, Size_CND3, weight);

        h_LayerMult_CND1_epCDn->Fill(LayerMult_CND1, weight);
        h_LayerMult_CND2_epCDn->Fill(LayerMult_CND2, weight);
        h_LayerMult_CND3_epCDn->Fill(LayerMult_CND3, weight);

        h_LayerMult_CND1_VS_LayerMult_CND2_epCDn->Fill(LayerMult_CND1, LayerMult_CND2, weight);
        h_LayerMult_CND1_VS_LayerMult_CND3_epCDn->Fill(LayerMult_CND1, LayerMult_CND3, weight);
        h_LayerMult_CND2_VS_LayerMult_CND3_epCDn->Fill(LayerMult_CND2, LayerMult_CND3, weight);

        h_ToF_n_epCDn->Fill(ToF, weight);
        h_ToF_zoomout_epCDn->Fill(ToF, weight);
        h_P_n_VS_ToF_n_epCDn->Fill(ToF, P_n_3v.Mag(), weight);
        h_theta_n_VS_ToF_n_epCDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_ToF_n_epCDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_ToF_n_epCDn->Fill(ToF, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_ToF_n_epCDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_ToF_n_epCDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_ToF_n_epCDn->Fill(ToF, dpp, weight);
        h_beta_n_VS_ToF_n_epCDn->Fill(ToF, beta, weight);
        h_E_miss_VS_ToF_n_epCDn->Fill(ToF, E_miss, weight);
        h_M_miss_VS_ToF_n_epCDn->Fill(ToF, M_miss, weight);
        h_path_VS_ToF_n_epCDn->Fill(ToF, path * 100, weight);
        h_theta_n_miss_VS_ToF_n_epCDn->Fill(ToF, theta_n_miss, weight);
        h_nSector_VS_ToF_n_epCDn->Fill(ToF, nSector, weight);
    } else if (pInFD) {
        h_xB_VS_M_miss_epFDn->Fill(xB, M_miss, weight);

        h_theta_n_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
        h_theta_n_VS_phi_n_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
        h_theta_n_VS_beta_n_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

        h_P_n_epFDn->Fill(P_n_3v.Mag(), weight);
        h_P_n_VS_theta_n_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

        h_P_miss_epFDn->Fill(P_miss_3v.Mag(), weight);
        h_P_miss_VS_theta_miss_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);

        h_dpp_allN_epFDn->Fill(dpp, weight);
        h_theta_n_miss_allN_epFDn->Fill(theta_n_miss, weight);
        h_dpp_VS_theta_n_miss_epFDn->Fill(dpp, theta_n_miss, weight);

        if (isGN) {
            h_dpp_goodN_epFDn->Fill(dpp, weight);
            h_theta_n_miss_goodN_epFDn->Fill(theta_n_miss, weight);
        } else if (isBN) {
            h_dpp_badN_epFDn->Fill(dpp, weight);
            h_theta_n_miss_badN_epFDn->Fill(theta_n_miss, weight);
        }

        if (theta_n_miss < 25.) {
            h_dpp_allN_for_theta_n_miss_less_than_25_epFDn->Fill(dpp, weight);
        }

        if (dpp < 0.5) {
            h_theta_n_miss_allN_for_dpp_less_than_05_epFDn->Fill(theta_n_miss, weight);

            if (dpp < 0.3) {
                h_theta_n_miss_allN_for_dpp_less_than_03_epFDn->Fill(theta_n_miss, weight);
            }
        }

        h_E_p_epFDn->Fill(E_p, weight);
        h_E_miss_epFDn->Fill(E_miss, weight);
        h_M_miss_epFDn->Fill(M_miss, weight);
        h_M_miss_VS_P_n_epFDn->Fill(P_n_3v.Mag(), M_miss, weight);
        h_M_miss_VS_theta_n_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
        h_M_miss_VS_phi_n_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
        h_M_miss_VS_P_miss_epFDn->Fill(P_miss_3v.Mag(), M_miss, weight);
        h_M_miss_VS_theta_miss_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
        h_M_miss_VS_phi_miss_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

        h_P_n_minus_P_miss_epFDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
        h_P_n_x_minus_P_miss_x_epFDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
        h_P_n_y_minus_P_miss_y_epFDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
        h_P_n_z_minus_P_miss_z_epFDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

        h_P_n_VS_P_miss_epFDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
        h_P_n_x_VS_P_miss_x_epFDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
        h_P_n_y_VS_P_miss_y_epFDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
        h_P_n_z_VS_P_miss_z_epFDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

        h_theta_n_p_epFDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
        h_theta_p_n_VS_P_p_epFDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

        h_xB_epFDn->Fill(xB, weight);

        h_Edep_CND_epFDn->Fill(Edep_CND, weight);
        h_P_n_VS_Edep_CND_epFDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CND_epFDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CND_epFDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CND_epFDn->Fill(Edep_CND, dpp, weight);
        h_beta_n_VS_Edep_CND_epFDn->Fill(Edep_CND, beta, weight);
        h_E_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, E_miss, weight);
        h_M_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, M_miss, weight);
        h_path_VS_Edep_CND_epFDn->Fill(Edep_CND, path * 100, weight);
        h_theta_n_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, theta_n_miss, weight);
        h_ToF_VS_Edep_CND_epFDn->Fill(Edep_CND, ToF, weight);
        h_nSector_VS_Edep_CND_epFDn->Fill(Edep_CND, nSector, weight);
        h_Edep_CND1_VS_Edep_CND_epFDn->Fill(Edep_CND, Edep_CND1, weight);
        h_Edep_CND2_VS_Edep_CND_epFDn->Fill(Edep_CND, Edep_CND2, weight);
        h_Edep_CND3_VS_Edep_CND_epFDn->Fill(Edep_CND, Edep_CND3, weight);

        h_Edep_CTOF_epFDn->Fill(Edep_CTOF, weight);
        h_P_n_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, dpp, weight);
        h_beta_n_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, beta, weight);
        h_E_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, E_miss, weight);
        h_M_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, M_miss, weight);
        h_path_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, path * 100, weight);
        h_theta_n_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
        h_ToF_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, ToF, weight);
        h_nSector_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, nSector, weight);
        h_Edep_CND1_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, Edep_CND1, weight);
        h_Edep_CND2_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, Edep_CND2, weight);
        h_Edep_CND3_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, Edep_CND3, weight);

        // h_Edep_single_epFDn->Fill(Edep_single, weight);
        // h_P_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n_3v.Mag(), weight);
        // h_theta_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n_3v.Theta() * 180. / M_PI, weight);
        // h_phi_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n_3v.Phi() * 180. / M_PI, weight);
        // h_P_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss_3v.Mag(), weight);
        // h_theta_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss_3v.Theta() * 180. / M_PI, weight);
        // h_phi_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss_3v.Phi() * 180. / M_PI, weight);
        // h_dpp_VS_Edep_single_epFDn->Fill(Edep_single, dpp, weight);
        // h_beta_n_VS_Edep_single_epFDn->Fill(Edep_single, beta, weight);
        // h_E_miss_VS_Edep_single_epFDn->Fill(Edep_single, E_miss, weight);
        // h_M_miss_VS_Edep_single_epFDn->Fill(Edep_single, M_miss, weight);
        // h_path_VS_Edep_single_epFDn->Fill(Edep_single, path * 100, weight);
        // h_theta_n_miss_VS_Edep_single_epFDn->Fill(Edep_single, theta_n_miss, weight);
        // h_ToF_VS_Edep_single_epFDn->Fill(Edep_single, ToF, weight);
        // h_nSector_VS_Edep_single_epFDn->Fill(Edep_single, nSector, weight);

        h_Edep_CND1_epFDn->Fill(Edep_CND1, weight);
        h_P_n_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CND1_epFDn->Fill(Edep_CND1, dpp, weight);
        h_beta_n_VS_Edep_CND1_epFDn->Fill(Edep_CND1, beta, weight);
        h_E_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, E_miss, weight);
        h_M_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, M_miss, weight);
        h_path_VS_Edep_CND1_epFDn->Fill(Edep_CND1, path * 100, weight);
        h_theta_n_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
        h_ToF_VS_Edep_CND1_epFDn->Fill(Edep_CND1, ToF, weight);
        h_nSector_VS_Edep_CND1_epFDn->Fill(Edep_CND1, nSector, weight);
        h_Edep_CND2_VS_Edep_CND1_epFDn->Fill(Edep_CND1, Edep_CND2, weight);
        h_Edep_CND3_VS_Edep_CND1_epFDn->Fill(Edep_CND1, Edep_CND3, weight);

        h_Edep_CND2_epFDn->Fill(Edep_CND2, weight);
        h_P_n_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CND2_epFDn->Fill(Edep_CND2, dpp, weight);
        h_beta_n_VS_Edep_CND2_epFDn->Fill(Edep_CND2, beta, weight);
        h_E_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, E_miss, weight);
        h_M_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, M_miss, weight);
        h_path_VS_Edep_CND2_epFDn->Fill(Edep_CND2, path * 100, weight);
        h_theta_n_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
        h_ToF_VS_Edep_CND2_epFDn->Fill(Edep_CND2, ToF, weight);
        h_nSector_VS_Edep_CND2_epFDn->Fill(Edep_CND2, nSector, weight);
        h_Edep_CND3_VS_Edep_CND2_epFDn->Fill(Edep_CND2, Edep_CND3, weight);

        h_Edep_CND3_epFDn->Fill(Edep_CND3, weight);
        h_P_n_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
        h_theta_n_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Edep_CND3_epFDn->Fill(Edep_CND3, dpp, weight);
        h_beta_n_VS_Edep_CND3_epFDn->Fill(Edep_CND3, beta, weight);
        h_E_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, E_miss, weight);
        h_M_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, M_miss, weight);
        h_path_VS_Edep_CND3_epFDn->Fill(Edep_CND3, path * 100, weight);
        h_theta_n_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
        h_ToF_VS_Edep_CND3_epFDn->Fill(Edep_CND3, ToF, weight);
        h_nSector_VS_Edep_CND3_epFDn->Fill(Edep_CND3, nSector, weight);

        h_Size_CND1_epFDn->Fill(Size_CND1, weight);
        h_Edep_CND_VS_Size_CND1_epFDn->Fill(Size_CND1, Edep_CND, weight);
        h_Edep_CND1_VS_Size_CND1_epFDn->Fill(Size_CND1, Edep_CND1, weight);
        h_Edep_CND2_VS_Size_CND1_epFDn->Fill(Size_CND1, Edep_CND2, weight);
        h_Edep_CND3_VS_Size_CND1_epFDn->Fill(Size_CND1, Edep_CND3, weight);
        h_P_n_VS_Size_CND1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
        h_theta_n_VS_Size_CND1_epFDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Size_CND1_epFDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Size_CND1_epFDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Size_CND1_epFDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Size_CND1_epFDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Size_CND1_epFDn->Fill(Size_CND1, dpp, weight);
        h_beta_n_VS_Size_CND1_epFDn->Fill(Size_CND1, beta, weight);
        h_E_miss_VS_Size_CND1_epFDn->Fill(Size_CND1, E_miss, weight);
        h_M_miss_VS_Size_CND1_epFDn->Fill(Size_CND1, M_miss, weight);
        h_path_VS_Size_CND1_epFDn->Fill(Size_CND1, path * 100, weight);
        h_theta_n_miss_VS_Size_CND1_epFDn->Fill(Size_CND1, theta_n_miss, weight);
        h_ToF_VS_Size_CND1_epFDn->Fill(Size_CND1, ToF, weight);
        h_nSector_VS_Size_CND1_epFDn->Fill(Size_CND1, nSector, weight);

        h_Size_CND2_epFDn->Fill(Size_CND2, weight);
        h_Edep_CND_VS_Size_CND2_epFDn->Fill(Size_CND2, Edep_CND, weight);
        h_Edep_CND1_VS_Size_CND2_epFDn->Fill(Size_CND2, Edep_CND1, weight);
        h_Edep_CND2_VS_Size_CND2_epFDn->Fill(Size_CND2, Edep_CND2, weight);
        h_Edep_CND3_VS_Size_CND2_epFDn->Fill(Size_CND2, Edep_CND3, weight);
        h_P_n_VS_Size_CND2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
        h_theta_n_VS_Size_CND2_epFDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Size_CND2_epFDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Size_CND2_epFDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Size_CND2_epFDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Size_CND2_epFDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Size_CND2_epFDn->Fill(Size_CND2, dpp, weight);
        h_beta_n_VS_Size_CND2_epFDn->Fill(Size_CND2, beta, weight);
        h_E_miss_VS_Size_CND2_epFDn->Fill(Size_CND2, E_miss, weight);
        h_M_miss_VS_Size_CND2_epFDn->Fill(Size_CND2, M_miss, weight);
        h_path_VS_Size_CND2_epFDn->Fill(Size_CND2, path * 100, weight);
        h_theta_n_miss_VS_Size_CND2_epFDn->Fill(Size_CND2, theta_n_miss, weight);
        h_ToF_VS_Size_CND2_epFDn->Fill(Size_CND2, ToF, weight);
        h_nSector_VS_Size_CND2_epFDn->Fill(Size_CND2, nSector, weight);

        h_Size_CND3_epFDn->Fill(Size_CND3, weight);
        h_Edep_CND_VS_Size_CND3_epFDn->Fill(Size_CND3, Edep_CND, weight);
        h_Edep_CND1_VS_Size_CND3_epFDn->Fill(Size_CND3, Edep_CND1, weight);
        h_Edep_CND2_VS_Size_CND3_epFDn->Fill(Size_CND3, Edep_CND2, weight);
        h_Edep_CND3_VS_Size_CND3_epFDn->Fill(Size_CND3, Edep_CND3, weight);
        h_P_n_VS_Size_CND3_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND3_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
        h_P_n_VS_Size_CND3_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
        h_theta_n_VS_Size_CND3_epFDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_Size_CND3_epFDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_Size_CND3_epFDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_Size_CND3_epFDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_Size_CND3_epFDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_Size_CND3_epFDn->Fill(Size_CND3, dpp, weight);
        h_beta_n_VS_Size_CND3_epFDn->Fill(Size_CND3, beta, weight);
        h_E_miss_VS_Size_CND3_epFDn->Fill(Size_CND3, E_miss, weight);
        h_M_miss_VS_Size_CND3_epFDn->Fill(Size_CND3, M_miss, weight);
        h_path_VS_Size_CND3_epFDn->Fill(Size_CND3, path * 100, weight);
        h_theta_n_miss_VS_Size_CND3_epFDn->Fill(Size_CND3, theta_n_miss, weight);
        h_ToF_VS_Size_CND3_epFDn->Fill(Size_CND3, ToF, weight);
        h_nSector_VS_Size_CND3_epFDn->Fill(Size_CND3, nSector, weight);

        h_Size_CND1_VS_Size_CND2_epFDn->Fill(Size_CND1, Size_CND2, weight);
        h_Size_CND1_VS_Size_CND3_epFDn->Fill(Size_CND1, Size_CND3, weight);
        h_Size_CND2_VS_Size_CND3_epFDn->Fill(Size_CND2, Size_CND3, weight);

        h_LayerMult_CND1_epFDn->Fill(LayerMult_CND1, weight);
        h_LayerMult_CND2_epFDn->Fill(LayerMult_CND2, weight);
        h_LayerMult_CND3_epFDn->Fill(LayerMult_CND3, weight);

        h_LayerMult_CND1_VS_LayerMult_CND2_epFDn->Fill(LayerMult_CND1, LayerMult_CND2, weight);
        h_LayerMult_CND1_VS_LayerMult_CND3_epFDn->Fill(LayerMult_CND1, LayerMult_CND3, weight);
        h_LayerMult_CND2_VS_LayerMult_CND3_epFDn->Fill(LayerMult_CND2, LayerMult_CND3, weight);

        h_ToF_n_epFDn->Fill(ToF, weight);
        h_P_n_VS_ToF_n_epFDn->Fill(ToF, P_n_3v.Mag(), weight);
        h_theta_n_VS_ToF_n_epFDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
        h_phi_n_VS_ToF_n_epFDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
        h_P_miss_VS_ToF_n_epFDn->Fill(ToF, P_miss_3v.Mag(), weight);
        h_theta_miss_VS_ToF_n_epFDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
        h_phi_miss_VS_ToF_n_epFDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
        h_dpp_VS_ToF_n_epFDn->Fill(ToF, dpp, weight);
        h_beta_n_VS_ToF_n_epFDn->Fill(ToF, beta, weight);
        h_E_miss_VS_ToF_n_epFDn->Fill(ToF, E_miss, weight);
        h_M_miss_VS_ToF_n_epFDn->Fill(ToF, M_miss, weight);
        h_path_VS_ToF_n_epFDn->Fill(ToF, path * 100, weight);
        h_theta_n_miss_VS_ToF_n_epFDn->Fill(ToF, theta_n_miss, weight);
        h_nSector_VS_ToF_n_epFDn->Fill(ToF, nSector, weight);
    }
}

// UpdateBS0CHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateBS0CHistograms(bool pInCD, bool pInFD, TVector3 P_n_3v, TVector3 v_hit_3v, double beta, double path, double ToF,
                                          double weight) {
    if (pInCD) {
        h_dbeta_n_BS0C_Step0_epCDn->Fill(beta - (path * 100) / (ToF * c), weight);
        h_dbeta_n_VS_P_n_BS0C_Step0_epCDn->Fill(P_n_3v.Mag(), beta - (path * 100) / (ToF * c), weight);
        h_dbeta_n_VS_ToF_BS0C_Step0_epCDn->Fill(ToF, beta - (path * 100) / (ToF * c), weight);

        h_Vhit_z_n_BS0C_Step0_epCDn->Fill(v_hit_3v.Z(), weight);

        h_ToF_n_BS0C_Step0_epCDn->Fill(ToF, weight);
    } else if (pInFD) {
        h_dbeta_n_BS0C_Step0_epFDn->Fill(beta - (path * 100) / (ToF * c), weight);
        h_dbeta_n_VS_P_n_BS0C_Step0_epFDn->Fill(P_n_3v.Mag(), beta - (path * 100) / (ToF * c), weight);
        h_dbeta_n_VS_ToF_BS0C_Step0_epFDn->Fill(ToF, beta - (path * 100) / (ToF * c), weight);

        h_Vhit_z_n_BS0C_Step0_epFDn->Fill(v_hit_3v.Z(), weight);

        h_ToF_n_BS0C_Step0_epFDn->Fill(ToF, weight);
    }
}

// UpdateAS0CHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateAS0CHistograms(bool pInCD, bool pInFD, TVector3 P_n_3v, TVector3 v_hit_3v, double beta, double path, double ToF,
                                          double weight) {
    if (pInCD) {
        h_dbeta_n_AS0C_Step0_epCDn->Fill(beta - (path * 100) / (ToF * c), weight);
        h_dbeta_n_VS_P_n_AS0C_Step0_epCDn->Fill(P_n_3v.Mag(), beta - (path * 100) / (ToF * c), weight);
        h_dbeta_n_VS_ToF_AS0C_Step0_epCDn->Fill(ToF, beta - (path * 100) / (ToF * c), weight);

        h_Vhit_z_n_AS0C_Step0_epCDn->Fill(v_hit_3v.Z(), weight);

        h_ToF_n_AS0C_Step0_epCDn->Fill(ToF, weight);
    } else if (pInFD) {
        h_dbeta_n_AS0C_Step0_epFDn->Fill(beta - (path * 100) / (ToF * c), weight);
        h_dbeta_n_VS_P_n_AS0C_Step0_epFDn->Fill(P_n_3v.Mag(), beta - (path * 100) / (ToF * c), weight);
        h_dbeta_n_VS_ToF_AS0C_Step0_epFDn->Fill(ToF, beta - (path * 100) / (ToF * c), weight);

        h_Vhit_z_n_AS0C_Step0_epFDn->Fill(v_hit_3v.Z(), weight);

        h_ToF_n_AS0C_Step0_epFDn->Fill(ToF, weight);
    }
}

// UpdateStep0Histograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateStep0Histograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 P_p_3v, TVector3 P_miss_3v, TVector3 P_n_3v,
                                           double E_p, double E_miss, double M_miss, double xB, double dpp, double theta_n_miss,
                                           double Edep_CND, double Edep_CND1, double Edep_CND2, double Edep_CND3, double Edep_CTOF,
                                           double nSector, double Size_CND1, double Size_CND2, double Size_CND3,
                                           double LayerMult_CND1, double LayerMult_CND2, double LayerMult_CND3,
                                           double beta, double path, double ToF, double weight) {
    if (pInCD) {
        h_dpp_allN_Step0_epCDn->Fill(dpp, weight);
        h_theta_n_miss_allN_Step0_epCDn->Fill(theta_n_miss, weight);
        h_dpp_VS_theta_n_miss_allN_Step0_epCDn->Fill(dpp, theta_n_miss, weight);

        if (theta_n_miss < 25.) {
            h_dpp_allN_for_theta_n_miss_less_than_25_Step0_epCDn->Fill(dpp, weight);
        }

        if (dpp < 0.5) {
            h_theta_n_miss_allN_for_dpp_less_than_05_Step0_epCDn->Fill(theta_n_miss, weight);

            if (dpp < 0.3) {
                h_theta_n_miss_allN_for_dpp_less_than_03_Step0_epCDn->Fill(theta_n_miss, weight);
            }
        }

        if (isGN) {
            h_theta_n_goodN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_goodN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI,
                                                       weight);
            h_theta_n_VS_beta_n_goodN_Step0_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_goodN_Step0_epCDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_goodN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_goodN_Step0_epCDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_goodN_Step0_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(),
                                                           weight);
            h_P_miss_VS_phi_miss_goodN_Step0_epCDn->
                    Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_goodN_Step0_epCDn->Fill(dpp, weight);
            h_theta_n_miss_goodN_Step0_epCDn->Fill(theta_n_miss, weight);

            h_E_p_goodN_Step0_epCDn->Fill(E_p, weight);
            h_E_miss_goodN_Step0_epCDn->Fill(E_miss, weight);
            h_M_miss_goodN_Step0_epCDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_goodN_Step0_epCDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_goodN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_goodN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_goodN_Step0_epCDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_goodN_Step0_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_goodN_Step0_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_goodN_Step0_epCDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_goodN_Step0_epCDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_goodN_Step0_epCDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_goodN_Step0_epCDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_goodN_Step0_epCDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_goodN_Step0_epCDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_goodN_Step0_epCDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_goodN_Step0_epCDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_goodN_Step0_epCDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_goodN_Step0_epCDn->
                    Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_goodN_Step0_epCDn->Fill(xB, weight);

            h_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_goodN_Step0_epCDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_goodN_Step0_epCDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_goodN_Step0_epCDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_goodN_Step0_epCDn->Fill(ToF, weight);
            h_P_n_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_goodN_Step0_epCDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_goodN_Step0_epCDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_goodN_Step0_epCDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_goodN_Step0_epCDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_goodN_Step0_epCDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_goodN_epCDn->Fill(xB, M_miss, weight);

            h_beta_n_goodN_Step0_epCDn->Fill(beta, weight);
        } else if (isBN) {
            h_theta_n_badN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_badN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI,
                                                      weight);
            h_theta_n_VS_beta_n_badN_Step0_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_badN_Step0_epCDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_badN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_badN_Step0_epCDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_badN_Step0_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(),
                                                          weight);
            h_P_miss_VS_phi_miss_badN_Step0_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_badN_Step0_epCDn->Fill(dpp, weight);
            h_theta_n_miss_badN_Step0_epCDn->Fill(theta_n_miss, weight);

            h_E_p_badN_Step0_epCDn->Fill(E_p, weight);
            h_E_miss_badN_Step0_epCDn->Fill(E_miss, weight);
            h_M_miss_badN_Step0_epCDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_badN_Step0_epCDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_badN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_badN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_badN_Step0_epCDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_badN_Step0_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_badN_Step0_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_badN_Step0_epCDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_badN_Step0_epCDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_badN_Step0_epCDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_badN_Step0_epCDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_badN_Step0_epCDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_badN_Step0_epCDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_badN_Step0_epCDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_badN_Step0_epCDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_badN_Step0_epCDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_badN_Step0_epCDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_badN_Step0_epCDn->Fill(xB, weight);

            h_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_badN_Step0_epCDn->
                    Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_badN_Step0_epCDn->
                    Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_badN_Step0_epCDn->
                    Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_badN_Step0_epCDn->
                    Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_badN_Step0_epCDn->
                    Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_badN_Step0_epCDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_badN_Step0_epCDn->
                    Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_badN_Step0_epCDn->
                    Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_badN_Step0_epCDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_badN_Step0_epCDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_badN_Step0_epCDn->Fill(ToF, weight);
            h_P_n_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_badN_Step0_epCDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_badN_Step0_epCDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_badN_Step0_epCDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_badN_Step0_epCDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_badN_Step0_epCDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_badN_epCDn->Fill(xB, M_miss, weight);

            h_beta_n_badN_Step0_epCDn->Fill(beta, weight);
        }
    } else if (pInFD) {
        h_dpp_allN_Step0_epFDn->Fill(dpp, weight);
        h_theta_n_miss_allN_Step0_epFDn->Fill(theta_n_miss, weight);
        h_dpp_VS_theta_n_miss_allN_Step0_epFDn->Fill(dpp, theta_n_miss, weight);

        if (theta_n_miss < 25.) {
            h_dpp_allN_for_theta_n_miss_less_than_25_Step0_epFDn->Fill(dpp, weight);
        }

        if (dpp < 0.5) {
            h_theta_n_miss_allN_for_dpp_less_than_05_Step0_epFDn->Fill(theta_n_miss, weight);

            if (dpp < 0.3) {
                h_theta_n_miss_allN_for_dpp_less_than_03_Step0_epFDn->Fill(theta_n_miss, weight);
            }
        }

        if (isGN) {
            h_theta_n_goodN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_goodN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI,
                                                       weight);
            h_theta_n_VS_beta_n_goodN_Step0_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_goodN_Step0_epFDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_goodN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_goodN_Step0_epFDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_goodN_Step0_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(),
                                                           weight);
            h_P_miss_VS_phi_miss_goodN_Step0_epFDn->
                    Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_goodN_Step0_epFDn->Fill(dpp, weight);
            h_theta_n_miss_goodN_Step0_epFDn->Fill(theta_n_miss, weight);

            h_E_p_goodN_Step0_epFDn->Fill(E_p, weight);
            h_E_miss_goodN_Step0_epFDn->Fill(E_miss, weight);
            h_M_miss_goodN_Step0_epFDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_goodN_Step0_epFDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_goodN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_goodN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_goodN_Step0_epFDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_goodN_Step0_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_goodN_Step0_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_goodN_Step0_epFDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_goodN_Step0_epFDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_goodN_Step0_epFDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_goodN_Step0_epFDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_goodN_Step0_epFDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_goodN_Step0_epFDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_goodN_Step0_epFDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_goodN_Step0_epFDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_goodN_Step0_epFDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_goodN_Step0_epFDn->
                    Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_goodN_Step0_epFDn->Fill(xB, weight);

            h_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_goodN_Step0_epFDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_goodN_Step0_epFDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_goodN_Step0_epFDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_goodN_Step0_epFDn->Fill(ToF, weight);
            h_P_n_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_goodN_Step0_epFDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_goodN_Step0_epFDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_goodN_Step0_epFDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_goodN_Step0_epFDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_goodN_Step0_epFDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_goodN_epFDn->Fill(xB, M_miss, weight);

            h_beta_n_goodN_Step0_epFDn->Fill(beta, weight);
        } else if (isBN) {
            h_theta_n_badN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_badN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI,
                                                      weight);
            h_theta_n_VS_beta_n_badN_Step0_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_badN_Step0_epFDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_badN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_badN_Step0_epFDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_badN_Step0_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(),
                                                          weight);
            h_P_miss_VS_phi_miss_badN_Step0_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_badN_Step0_epFDn->Fill(dpp, weight);
            h_theta_n_miss_badN_Step0_epFDn->Fill(theta_n_miss, weight);

            h_E_p_badN_Step0_epFDn->Fill(E_p, weight);
            h_E_miss_badN_Step0_epFDn->Fill(E_miss, weight);
            h_M_miss_badN_Step0_epFDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_badN_Step0_epFDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_badN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_badN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_badN_Step0_epFDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_badN_Step0_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_badN_Step0_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_badN_Step0_epFDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_badN_Step0_epFDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_badN_Step0_epFDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_badN_Step0_epFDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_badN_Step0_epFDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_badN_Step0_epFDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_badN_Step0_epFDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_badN_Step0_epFDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_badN_Step0_epFDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_badN_Step0_epFDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_badN_Step0_epFDn->Fill(xB, weight);

            h_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_badN_Step0_epFDn->
                    Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_badN_Step0_epFDn->
                    Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_badN_Step0_epFDn->
                    Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_badN_Step0_epFDn->
                    Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_badN_Step0_epFDn->
                    Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_badN_Step0_epFDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_badN_Step0_epFDn->
                    Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_badN_Step0_epFDn->
                    Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_badN_Step0_epFDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_badN_Step0_epFDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_badN_Step0_epFDn->Fill(ToF, weight);
            h_P_n_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_badN_Step0_epFDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_badN_Step0_epFDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_badN_Step0_epFDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_badN_Step0_epFDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_badN_Step0_epFDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_badN_epFDn->Fill(xB, M_miss, weight);

            h_beta_n_badN_Step0_epFDn->Fill(beta, weight);
        }
    }
}

// UpdateStep1Histograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateStep1Histograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 P_p_3v, TVector3 P_miss_3v, TVector3 P_n_3v,
                                           double E_p, double E_miss, double M_miss, double xB, double dpp, double theta_n_miss,
                                           double Edep_CND, double Edep_CND1, double Edep_CND2, double Edep_CND3, double Edep_CTOF,
                                           double nSector, double Size_CND1, double Size_CND2, double Size_CND3,
                                           double LayerMult_CND1, double LayerMult_CND2, double LayerMult_CND3,
                                           double beta, double path, double ToF, double weight) {
    if (pInCD) {
        h_dpp_allN_Step1_epCDn->Fill(dpp, weight);
        h_theta_n_miss_allN_Step1_epCDn->Fill(theta_n_miss, weight);
        h_dpp_VS_theta_n_miss_allN_Step1_epCDn->Fill(dpp, theta_n_miss, weight);

        if (theta_n_miss < 25.) {
            h_dpp_allN_for_theta_n_miss_less_than_25_Step1_epCDn->Fill(dpp, weight);
        }

        if (dpp < 0.5) {
            h_theta_n_miss_allN_for_dpp_less_than_05_Step1_epCDn->Fill(theta_n_miss, weight);

            if (dpp < 0.3) {
                h_theta_n_miss_allN_for_dpp_less_than_03_Step1_epCDn->Fill(theta_n_miss, weight);
            }
        }

        if (isGN) {
            h_theta_n_goodN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_goodN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI,
                                                       weight);
            h_theta_n_VS_beta_n_goodN_Step1_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_goodN_Step1_epCDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_goodN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_goodN_Step1_epCDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_goodN_Step1_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(),
                                                           weight);
            h_P_miss_VS_phi_miss_goodN_Step1_epCDn->
                    Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_goodN_Step1_epCDn->Fill(dpp, weight);
            h_theta_n_miss_goodN_Step1_epCDn->Fill(theta_n_miss, weight);

            h_E_p_goodN_Step1_epCDn->Fill(E_p, weight);
            h_E_miss_goodN_Step1_epCDn->Fill(E_miss, weight);
            h_M_miss_goodN_Step1_epCDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_goodN_Step1_epCDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_goodN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_goodN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_goodN_Step1_epCDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_goodN_Step1_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_goodN_Step1_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_goodN_Step1_epCDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_goodN_Step1_epCDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_goodN_Step1_epCDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_goodN_Step1_epCDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_goodN_Step1_epCDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_goodN_Step1_epCDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_goodN_Step1_epCDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_goodN_Step1_epCDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_goodN_Step1_epCDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_goodN_Step1_epCDn->
                    Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_goodN_Step1_epCDn->Fill(xB, weight);

            h_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_goodN_Step1_epCDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_goodN_Step1_epCDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_goodN_Step1_epCDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_goodN_Step1_epCDn->Fill(ToF, weight);
            h_P_n_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_goodN_Step1_epCDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_goodN_Step1_epCDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_goodN_Step1_epCDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_goodN_Step1_epCDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_goodN_Step1_epCDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_goodN_epCDn->Fill(xB, M_miss, weight);

            h_beta_n_goodN_Step1_epCDn->Fill(beta, weight);
        } else if (isBN) {
            h_theta_n_badN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_badN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI,
                                                      weight);
            h_theta_n_VS_beta_n_badN_Step1_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_badN_Step1_epCDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_badN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_badN_Step1_epCDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_badN_Step1_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(),
                                                          weight);
            h_P_miss_VS_phi_miss_badN_Step1_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_badN_Step1_epCDn->Fill(dpp, weight);
            h_theta_n_miss_badN_Step1_epCDn->Fill(theta_n_miss, weight);

            h_E_p_badN_Step1_epCDn->Fill(E_p, weight);
            h_E_miss_badN_Step1_epCDn->Fill(E_miss, weight);
            h_M_miss_badN_Step1_epCDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_badN_Step1_epCDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_badN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_badN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_badN_Step1_epCDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_badN_Step1_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_badN_Step1_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_badN_Step1_epCDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_badN_Step1_epCDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_badN_Step1_epCDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_badN_Step1_epCDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_badN_Step1_epCDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_badN_Step1_epCDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_badN_Step1_epCDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_badN_Step1_epCDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_badN_Step1_epCDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_badN_Step1_epCDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_badN_Step1_epCDn->Fill(xB, weight);

            h_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_badN_Step1_epCDn->
                    Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_badN_Step1_epCDn->
                    Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_badN_Step1_epCDn->
                    Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_badN_Step1_epCDn->
                    Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_badN_Step1_epCDn->
                    Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_badN_Step1_epCDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_badN_Step1_epCDn->
                    Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_badN_Step1_epCDn->
                    Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_badN_Step1_epCDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_badN_Step1_epCDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_badN_Step1_epCDn->Fill(ToF, weight);
            h_P_n_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_badN_Step1_epCDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_badN_Step1_epCDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_badN_Step1_epCDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_badN_Step1_epCDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_badN_Step1_epCDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_badN_epCDn->Fill(xB, M_miss, weight);

            h_beta_n_badN_Step1_epCDn->Fill(beta, weight);
        }
    } else if (pInFD) {
        h_dpp_allN_Step1_epFDn->Fill(dpp, weight);
        h_theta_n_miss_allN_Step1_epFDn->Fill(theta_n_miss, weight);
        h_dpp_VS_theta_n_miss_allN_Step1_epFDn->Fill(dpp, theta_n_miss, weight);

        if (theta_n_miss < 25.) {
            h_dpp_allN_for_theta_n_miss_less_than_25_Step1_epFDn->Fill(dpp, weight);
        }

        if (dpp < 0.5) {
            h_theta_n_miss_allN_for_dpp_less_than_05_Step1_epFDn->Fill(theta_n_miss, weight);

            if (dpp < 0.3) {
                h_theta_n_miss_allN_for_dpp_less_than_03_Step1_epFDn->Fill(theta_n_miss, weight);
            }
        }

        if (isGN) {
            h_theta_n_goodN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_goodN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI,
                                                       weight);
            h_theta_n_VS_beta_n_goodN_Step1_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_goodN_Step1_epFDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_goodN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_goodN_Step1_epFDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_goodN_Step1_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(),
                                                           weight);
            h_P_miss_VS_phi_miss_goodN_Step1_epFDn->
                    Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_goodN_Step1_epFDn->Fill(dpp, weight);
            h_theta_n_miss_goodN_Step1_epFDn->Fill(theta_n_miss, weight);

            h_E_p_goodN_Step1_epFDn->Fill(E_p, weight);
            h_E_miss_goodN_Step1_epFDn->Fill(E_miss, weight);
            h_M_miss_goodN_Step1_epFDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_goodN_Step1_epFDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_goodN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_goodN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_goodN_Step1_epFDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_goodN_Step1_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_goodN_Step1_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_goodN_Step1_epFDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_goodN_Step1_epFDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_goodN_Step1_epFDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_goodN_Step1_epFDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_goodN_Step1_epFDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_goodN_Step1_epFDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_goodN_Step1_epFDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_goodN_Step1_epFDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_goodN_Step1_epFDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_goodN_Step1_epFDn->
                    Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_goodN_Step1_epFDn->Fill(xB, weight);

            h_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_goodN_Step1_epFDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI,
                                                              weight);
            h_phi_miss_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_goodN_Step1_epFDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_goodN_Step1_epFDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_goodN_Step1_epFDn->Fill(ToF, weight);
            h_P_n_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_goodN_Step1_epFDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_goodN_Step1_epFDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_goodN_Step1_epFDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_goodN_Step1_epFDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_goodN_Step1_epFDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_goodN_epFDn->Fill(xB, M_miss, weight);

            h_beta_n_goodN_Step1_epFDn->Fill(beta, weight);
        } else if (isBN) {
            h_theta_n_badN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_badN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI,
                                                      weight);
            h_theta_n_VS_beta_n_badN_Step1_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_badN_Step1_epFDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_badN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_badN_Step1_epFDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_badN_Step1_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(),
                                                          weight);
            h_P_miss_VS_phi_miss_badN_Step1_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_badN_Step1_epFDn->Fill(dpp, weight);
            h_theta_n_miss_badN_Step1_epFDn->Fill(theta_n_miss, weight);

            h_E_p_badN_Step1_epFDn->Fill(E_p, weight);
            h_E_miss_badN_Step1_epFDn->Fill(E_miss, weight);
            h_M_miss_badN_Step1_epFDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_badN_Step1_epFDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_badN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_badN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_badN_Step1_epFDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_badN_Step1_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_badN_Step1_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_badN_Step1_epFDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_badN_Step1_epFDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_badN_Step1_epFDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_badN_Step1_epFDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_badN_Step1_epFDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_badN_Step1_epFDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_badN_Step1_epFDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_badN_Step1_epFDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_badN_Step1_epFDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_badN_Step1_epFDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_badN_Step1_epFDn->Fill(xB, weight);

            h_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_badN_Step1_epFDn->
                    Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_badN_Step1_epFDn->
                    Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_badN_Step1_epFDn->
                    Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_badN_Step1_epFDn->
                    Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_badN_Step1_epFDn->
                    Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_badN_Step1_epFDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_badN_Step1_epFDn->
                    Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_badN_Step1_epFDn->
                    Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_badN_Step1_epFDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_badN_Step1_epFDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_badN_Step1_epFDn->Fill(ToF, weight);
            h_P_n_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_badN_Step1_epFDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_badN_Step1_epFDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_badN_Step1_epFDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_badN_Step1_epFDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_badN_Step1_epFDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_badN_epFDn->Fill(xB, M_miss, weight);

            h_beta_n_badN_Step1_epFDn->Fill(beta, weight);
        }
    }
}

// UpdateStep2prepBCHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateStep2prepBCHistograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 v_hit_3v, double ToF, double weight) {
    if (pInCD) {
        h_ToF_c_minus_VhitZ_BC_allN_Step2prep_epCDn->Fill(ToF * c - v_hit_3v.Z(), weight);
        h_ToF_c_minus_VhitZ_VS_VhitZ_BC_allN_Step2prep_epCDn->
                Fill(ToF * c - v_hit_3v.Z(), v_hit_3v.Z(), weight);
        h_ToF_c_minus_VhitZ_VS_ToF_BC_allN_Step2prep_epCDn->Fill(ToF * c - v_hit_3v.Z(), ToF, weight);

        if (isGN) {
            h_ToF_c_minus_VhitZ_BC_goodN_Step2prep_epCDn->Fill(ToF * c - v_hit_3v.Z(), weight);
            h_ToF_c_minus_VhitZ_VS_VhitZ_BC_goodN_Step2prep_epCDn->Fill(
                ToF * c - v_hit_3v.Z(), v_hit_3v.Z(), weight);
            h_ToF_c_minus_VhitZ_VS_ToF_BC_goodN_Step2prep_epCDn->Fill(ToF * c - v_hit_3v.Z(), ToF, weight);
        } else if (isBN) {
            h_ToF_c_minus_VhitZ_BC_badN_Step2prep_epCDn->Fill(ToF * c - v_hit_3v.Z(), weight);
            h_ToF_c_minus_VhitZ_VS_VhitZ_BC_badN_Step2prep_epCDn->Fill(
                ToF * c - v_hit_3v.Z(), v_hit_3v.Z(), weight);
            h_ToF_c_minus_VhitZ_VS_ToF_BC_badN_Step2prep_epCDn->Fill(ToF * c - v_hit_3v.Z(), ToF, weight);
        }
    } else if (pInFD) {
        h_ToF_c_minus_VhitZ_BC_allN_Step2prep_epFDn->Fill(ToF * c - v_hit_3v.Z(), weight);
        h_ToF_c_minus_VhitZ_VS_VhitZ_BC_allN_Step2prep_epFDn->
                Fill(ToF * c - v_hit_3v.Z(), v_hit_3v.Z(), weight);
        h_ToF_c_minus_VhitZ_VS_ToF_BC_allN_Step2prep_epFDn->Fill(ToF * c - v_hit_3v.Z(), ToF, weight);

        if (isGN) {
            h_ToF_c_minus_VhitZ_BC_goodN_Step2prep_epFDn->Fill(ToF * c - v_hit_3v.Z(), weight);
            h_ToF_c_minus_VhitZ_VS_VhitZ_BC_goodN_Step2prep_epFDn->Fill(
                ToF * c - v_hit_3v.Z(), v_hit_3v.Z(), weight);
            h_ToF_c_minus_VhitZ_VS_ToF_BC_goodN_Step2prep_epFDn->Fill(ToF * c - v_hit_3v.Z(), ToF, weight);
        } else if (isBN) {
            h_ToF_c_minus_VhitZ_BC_badN_Step2prep_epFDn->Fill(ToF * c - v_hit_3v.Z(), weight);
            h_ToF_c_minus_VhitZ_VS_VhitZ_BC_badN_Step2prep_epFDn->Fill(
                ToF * c - v_hit_3v.Z(), v_hit_3v.Z(), weight);
            h_ToF_c_minus_VhitZ_VS_ToF_BC_badN_Step2prep_epFDn->Fill(ToF * c - v_hit_3v.Z(), ToF, weight);
        }
    }
}

// UpdateStep2prepHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateStep2prepHistograms(bool pInCD, bool pInFD, bool isGN, bool isBN, int ldiff, int sdiff,
                                               TVector3 p_C_3v, TVector3 v_hit_3v, TVector3 P_n_3v, double dToF, double dToF_rel_pos,
                                               double dToF_rel_n, double dpp, double theta_n_miss, double Edep_CND, double beta, double path,
                                               double ToF, double weight) {
    if (pInCD) {
        // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step1_layer
        if (isGN) {
            h_sdiff_pos_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, weight);
            h_sdiff_pos_mom_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
            h_sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, v_hit_3v.Z(), weight);
            h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, ToF * c - v_hit_3v.Z(), weight);
            h_theta_n_goodN_Step2prep_layer_epCDn[ldiff + 3]->
                    Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_phi_n_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_ToF_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, ToF, weight);
            h_sdiff_pos_VS_path_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, path * 100, weight);
            h_sdiff_pos_VS_beta_n_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, beta, weight);
            h_sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, Edep_CND, weight);
            h_sdiff_pos_VS_theta_n_miss_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, theta_n_miss, weight);
            h_sdiff_pos_VS_dpp_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, dpp, weight);

            h_dToF_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(dToF, weight);
            h_sdiff_pos_VS_dToF_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, dToF, weight);
            h_dToF_rel_pos_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(dToF_rel_pos, weight);
            h_sdiff_pos_VS_dToF_rel_pos_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_pos, weight);
            h_dToF_rel_n_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(dToF_rel_n, weight);
            h_sdiff_pos_VS_dToF_rel_n_goodN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_n, weight);
        } else if (isBN) {
            h_sdiff_pos_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, weight);
            h_sdiff_pos_mom_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
            h_sdiff_pos_VS_VhitZ_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, v_hit_3v.Z(), weight);
            h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, ToF * c - v_hit_3v.Z(), weight);
            h_theta_n_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_phi_n_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_ToF_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, ToF, weight);
            h_sdiff_pos_VS_path_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, path * 100, weight);
            h_sdiff_pos_VS_beta_n_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, beta, weight);
            h_sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_epCDn[ldiff + 3]->
                    Fill(sdiff, Edep_CND, weight);
            h_sdiff_pos_VS_theta_n_miss_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, theta_n_miss, weight);
            h_sdiff_pos_VS_dpp_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, dpp, weight);

            h_dToF_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(dToF, weight);
            h_sdiff_pos_VS_dToF_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(sdiff, dToF, weight);
            h_dToF_rel_pos_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(dToF_rel_pos, weight);
            h_sdiff_pos_VS_dToF_rel_pos_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_pos, weight);
            h_dToF_rel_n_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(dToF_rel_n, weight);
            h_sdiff_pos_VS_dToF_rel_n_badN_Step2prep_layer_epCDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_n, weight);
        }
    } else if (pInFD) {
        if (isGN) // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step1_layer
        {
            h_sdiff_pos_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, weight);
            h_sdiff_pos_mom_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
            h_sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, v_hit_3v.Z(), weight);
            h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, ToF * c - v_hit_3v.Z(), weight);
            h_theta_n_goodN_Step2prep_layer_epFDn[ldiff + 3]->
                    Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_phi_n_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_ToF_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, ToF, weight);
            h_sdiff_pos_VS_path_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, path * 100, weight);
            h_sdiff_pos_VS_beta_n_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, beta, weight);
            h_sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, Edep_CND, weight);
            h_sdiff_pos_VS_theta_n_miss_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, theta_n_miss, weight);
            h_sdiff_pos_VS_dpp_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, dpp, weight);

            h_dToF_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(dToF, weight);
            h_sdiff_pos_VS_dToF_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, dToF, weight);
            h_dToF_rel_pos_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(dToF_rel_pos, weight);
            h_sdiff_pos_VS_dToF_rel_pos_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_pos, weight);
            h_dToF_rel_n_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(dToF_rel_n, weight);
            h_sdiff_pos_VS_dToF_rel_n_goodN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_n, weight);
        } else if (isBN) {
            h_sdiff_pos_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, weight);
            h_sdiff_pos_mom_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
            h_sdiff_pos_VS_VhitZ_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, v_hit_3v.Z(), weight);
            h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, ToF * c - v_hit_3v.Z(), weight);
            h_theta_n_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_phi_n_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_ToF_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, ToF, weight);
            h_sdiff_pos_VS_path_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, path * 100, weight);
            h_sdiff_pos_VS_beta_n_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, beta, weight);
            h_sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_epFDn[ldiff + 3]->
                    Fill(sdiff, Edep_CND, weight);
            h_sdiff_pos_VS_theta_n_miss_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, theta_n_miss, weight);
            h_sdiff_pos_VS_dpp_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, dpp, weight);

            h_dToF_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(dToF, weight);
            h_sdiff_pos_VS_dToF_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(sdiff, dToF, weight);
            h_dToF_rel_pos_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(dToF_rel_pos, weight);
            h_sdiff_pos_VS_dToF_rel_pos_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_pos, weight);
            h_dToF_rel_n_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(dToF_rel_n, weight);
            h_sdiff_pos_VS_dToF_rel_n_badN_Step2prep_layer_epFDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_n, weight);
        }
    }
}

// UpdateMonitorStep2prepHistograms1 function
// ======================================================================================================================================================================

void VetoHistograms::UpdateMonitorStep2prepHistograms1(bool Nearby_clusters_from_cPart_tracks, bool pInCD, bool pInFD, bool isGN, bool isBN,
                                                       double Edep_CND, double Edep_CTOF_pos, double weight) {
    if (Nearby_clusters_from_cPart_tracks) {
        if (pInCD) {
            if (isGN) {
                h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step2prep_epCDn->Fill(Edep_CND / Edep_CTOF_pos, weight);
            } else if (isBN) {
                h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step2prep_epCDn->Fill(Edep_CND / Edep_CTOF_pos, weight);
            }
        } else if (pInFD) {
            if (isGN) {
                h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step2prep_epFDn->Fill(Edep_CND / Edep_CTOF_pos, weight);
            } else if (isBN) {
                h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step2prep_epFDn->Fill(Edep_CND / Edep_CTOF_pos, weight);
            }
        }
    }
}

// UpdateMonitorStep2prepHistograms2 function
// ======================================================================================================================================================================

void VetoHistograms::UpdateMonitorStep2prepHistograms2(bool Nearby_clusters_from_cPart_tracks, bool pInCD, bool pInFD, bool isGN, bool isBN,
                                                       double Edep_CND, double ToF, TVector3 v_hit_3v, double weight) {
    if (Nearby_clusters_from_cPart_tracks) {
        if (pInCD) {
            if (isGN) {
                h_Edep_CND_goodN_withNearbyPos_Step2prep_epCDn->Fill(Edep_CND, weight);
            } else if (isBN) {
                h_Edep_CND_badN_withNearbyPos_Step2prep_epCDn->Fill(Edep_CND, weight);
            }
        } else if (pInFD) {
            if (isGN) {
                h_Edep_CND_goodN_withNearbyPos_Step2prep_epFDn->Fill(Edep_CND, weight);
            } else if (isBN) {
                h_Edep_CND_badN_withNearbyPos_Step2prep_epFDn->Fill(Edep_CND, weight);
            }
        }
    }

    if (pInCD) {
        if (isGN) {
            if (!Nearby_clusters_from_cPart_tracks)
                h_diff_ToFc_z_VS_Edep_noNear_goodN_Step2prep_epCDn->Fill(
                    ToF * c - v_hit_3v.Z(), Edep_CND, weight);
            else {
                h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step2prep_epCDn->Fill(
                    ToF * c - v_hit_3v.Z(), Edep_CND, weight);
            }
        } else if (isBN) {
            if (!Nearby_clusters_from_cPart_tracks)
                h_diff_ToFc_z_VS_Edep_noNear_badN_Step2prep_epCDn->Fill(
                    ToF * c - v_hit_3v.Z(), Edep_CND, weight);
            else {
                h_diff_ToFc_z_VS_Edep_yesNear_badN_Step2prep_epCDn->Fill(
                    ToF * c - v_hit_3v.Z(), Edep_CND, weight);
            }
        }
    } else if (pInFD) {
        if (isGN) {
            if (!Nearby_clusters_from_cPart_tracks)
                h_diff_ToFc_z_VS_Edep_noNear_goodN_Step2prep_epFDn->Fill(
                    ToF * c - v_hit_3v.Z(), Edep_CND, weight);
            else {
                h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step2prep_epFDn->Fill(
                    ToF * c - v_hit_3v.Z(), Edep_CND, weight);
            }
        } else if (isBN) {
            if (!Nearby_clusters_from_cPart_tracks)
                h_diff_ToFc_z_VS_Edep_noNear_badN_Step2prep_epFDn->Fill(
                    ToF * c - v_hit_3v.Z(), Edep_CND, weight);
            else {
                h_diff_ToFc_z_VS_Edep_yesNear_badN_Step2prep_epFDn->Fill(
                    ToF * c - v_hit_3v.Z(), Edep_CND, weight);
            }
        }
    }
}

// UpdateBS2CHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateBS2CHistograms(bool pInCD, bool pInFD, double Size_CND1, double Size_CND2, double Size_CND3, double LayerMult_CND1,
                                          double LayerMult_CND2, double LayerMult_CND3, double weight) {
    if (pInCD) {
        h_Size_CND1_BS2C_Step2_epCDn->Fill(Size_CND1, weight);
        h_Size_CND2_BS2C_Step2_epCDn->Fill(Size_CND2, weight);
        h_Size_CND3_BS2C_Step2_epCDn->Fill(Size_CND3, weight);

        h_Size_CND1_VS_Size_CND2_BS2C_Step2_epCDn->Fill(Size_CND1, Size_CND2, weight);
        h_Size_CND1_VS_Size_CND3_BS2C_Step2_epCDn->Fill(Size_CND1, Size_CND3, weight);
        h_Size_CND2_VS_Size_CND3_BS2C_Step2_epCDn->Fill(Size_CND2, Size_CND3, weight);

        h_LayerMult_CND1_BS2C_Step2_epCDn->Fill(LayerMult_CND1, weight);
        h_LayerMult_CND2_BS2C_Step2_epCDn->Fill(LayerMult_CND2, weight);
        h_LayerMult_CND3_BS2C_Step2_epCDn->Fill(LayerMult_CND3, weight);

        h_LayerMult_CND1_VS_LayerMult_CND2_BS2C_Step2_epCDn->Fill(LayerMult_CND1, LayerMult_CND2, weight);
        h_LayerMult_CND1_VS_LayerMult_CND3_BS2C_Step2_epCDn->Fill(LayerMult_CND1, LayerMult_CND3, weight);
        h_LayerMult_CND2_VS_LayerMult_CND3_BS2C_Step2_epCDn->Fill(LayerMult_CND2, LayerMult_CND3, weight);
    } else if (pInFD) {
        h_Size_CND1_BS2C_Step2_epFDn->Fill(Size_CND1, weight);
        h_Size_CND2_BS2C_Step2_epFDn->Fill(Size_CND2, weight);
        h_Size_CND3_BS2C_Step2_epFDn->Fill(Size_CND3, weight);

        h_Size_CND1_VS_Size_CND2_BS2C_Step2_epFDn->Fill(Size_CND1, Size_CND2, weight);
        h_Size_CND1_VS_Size_CND3_BS2C_Step2_epFDn->Fill(Size_CND1, Size_CND3, weight);
        h_Size_CND2_VS_Size_CND3_BS2C_Step2_epFDn->Fill(Size_CND2, Size_CND3, weight);

        h_LayerMult_CND1_BS2C_Step2_epFDn->Fill(LayerMult_CND1, weight);
        h_LayerMult_CND2_BS2C_Step2_epFDn->Fill(LayerMult_CND2, weight);
        h_LayerMult_CND3_BS2C_Step2_epFDn->Fill(LayerMult_CND3, weight);

        h_LayerMult_CND1_VS_LayerMult_CND2_BS2C_Step2_epFDn->Fill(LayerMult_CND1, LayerMult_CND2, weight);
        h_LayerMult_CND1_VS_LayerMult_CND3_BS2C_Step2_epFDn->Fill(LayerMult_CND1, LayerMult_CND3, weight);
        h_LayerMult_CND2_VS_LayerMult_CND3_BS2C_Step2_epFDn->Fill(LayerMult_CND2, LayerMult_CND3, weight);
    }
}

// UpdateAS2CHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateAS2CHistograms(bool pInCD, bool pInFD, double Size_CND1, double Size_CND2, double Size_CND3, double LayerMult_CND1,
                                          double LayerMult_CND2, double LayerMult_CND3, double weight) {
    if (pInCD) {
        h_Size_CND1_AS2C_Step2_epCDn->Fill(Size_CND1, weight);
        h_Size_CND2_AS2C_Step2_epCDn->Fill(Size_CND2, weight);
        h_Size_CND3_AS2C_Step2_epCDn->Fill(Size_CND3, weight);

        h_Size_CND1_VS_Size_CND2_AS2C_Step2_epCDn->Fill(Size_CND1, Size_CND2, weight);
        h_Size_CND1_VS_Size_CND3_AS2C_Step2_epCDn->Fill(Size_CND1, Size_CND3, weight);
        h_Size_CND2_VS_Size_CND3_AS2C_Step2_epCDn->Fill(Size_CND2, Size_CND3, weight);

        h_LayerMult_CND1_AS2C_Step2_epCDn->Fill(LayerMult_CND1, weight);
        h_LayerMult_CND2_AS2C_Step2_epCDn->Fill(LayerMult_CND2, weight);
        h_LayerMult_CND3_AS2C_Step2_epCDn->Fill(LayerMult_CND3, weight);

        h_LayerMult_CND1_VS_LayerMult_CND2_AS2C_Step2_epCDn->Fill(LayerMult_CND1, LayerMult_CND2, weight);
        h_LayerMult_CND1_VS_LayerMult_CND3_AS2C_Step2_epCDn->Fill(LayerMult_CND1, LayerMult_CND3, weight);
        h_LayerMult_CND2_VS_LayerMult_CND3_AS2C_Step2_epCDn->Fill(LayerMult_CND2, LayerMult_CND3, weight);
    } else if (pInFD) {
        h_Size_CND1_AS2C_Step2_epFDn->Fill(Size_CND1, weight);
        h_Size_CND2_AS2C_Step2_epFDn->Fill(Size_CND2, weight);
        h_Size_CND3_AS2C_Step2_epFDn->Fill(Size_CND3, weight);

        h_Size_CND1_VS_Size_CND2_AS2C_Step2_epFDn->Fill(Size_CND1, Size_CND2, weight);
        h_Size_CND1_VS_Size_CND3_AS2C_Step2_epFDn->Fill(Size_CND1, Size_CND3, weight);
        h_Size_CND2_VS_Size_CND3_AS2C_Step2_epFDn->Fill(Size_CND2, Size_CND3, weight);

        h_LayerMult_CND1_AS2C_Step2_epFDn->Fill(LayerMult_CND1, weight);
        h_LayerMult_CND2_AS2C_Step2_epFDn->Fill(LayerMult_CND2, weight);
        h_LayerMult_CND3_AS2C_Step2_epFDn->Fill(LayerMult_CND3, weight);

        h_LayerMult_CND1_VS_LayerMult_CND2_AS2C_Step2_epFDn->Fill(LayerMult_CND1, LayerMult_CND2, weight);
        h_LayerMult_CND1_VS_LayerMult_CND3_AS2C_Step2_epFDn->Fill(LayerMult_CND1, LayerMult_CND3, weight);
        h_LayerMult_CND2_VS_LayerMult_CND3_AS2C_Step2_epFDn->Fill(LayerMult_CND2, LayerMult_CND3, weight);
    }
}

// UpdateStep2Histograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateStep2Histograms(bool pInCD, bool pInFD, bool isGN, bool isBN, TVector3 P_p_3v, TVector3 P_miss_3v, TVector3 P_n_3v,
                                           double E_p, double E_miss, double M_miss, double xB, double dpp, double theta_n_miss, double Edep_CND,
                                           double Edep_CND1, double Edep_CND2, double Edep_CND3, double Edep_CTOF, double nSector, double Size_CND1,
                                           double Size_CND2, double Size_CND3, double LayerMult_CND1, double LayerMult_CND2, double LayerMult_CND3,
                                           double beta, double path, double ToF, double weight) {
    if (pInCD) {
        h_dpp_allN_Step2_epCDn->Fill(dpp, weight);
        h_theta_n_miss_allN_Step2_epCDn->Fill(theta_n_miss, weight);
        h_dpp_VS_theta_n_miss_allN_Step2_epCDn->Fill(dpp, theta_n_miss, weight);

        if (theta_n_miss < 25.) {
            h_dpp_allN_for_theta_n_miss_less_than_25_Step2_epCDn->Fill(dpp, weight);
        }

        if (dpp < 0.5) {
            h_theta_n_miss_allN_for_dpp_less_than_05_Step2_epCDn->Fill(theta_n_miss, weight);

            if (dpp < 0.3) {
                h_theta_n_miss_allN_for_dpp_less_than_03_Step2_epCDn->Fill(theta_n_miss, weight);
            }
        }

        if (isGN) {
            h_theta_n_goodN_Step2_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step2_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_goodN_Step2_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
            h_theta_n_VS_beta_n_goodN_Step2_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_goodN_Step2_epCDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_goodN_Step2_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_goodN_Step2_epCDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_goodN_Step2_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
            h_P_miss_VS_phi_miss_goodN_Step2_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_goodN_Step2_epCDn->Fill(dpp, weight);
            h_theta_n_miss_goodN_Step2_epCDn->Fill(theta_n_miss, weight);

            h_E_p_goodN_Step2_epCDn->Fill(E_p, weight);
            h_E_miss_goodN_Step2_epCDn->Fill(E_miss, weight);
            h_M_miss_goodN_Step2_epCDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_goodN_Step2_epCDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_goodN_Step2_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_goodN_Step2_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_goodN_Step2_epCDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_goodN_Step2_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_goodN_Step2_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_goodN_Step2_epCDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_goodN_Step2_epCDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_goodN_Step2_epCDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_goodN_Step2_epCDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_goodN_Step2_epCDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_goodN_Step2_epCDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_goodN_Step2_epCDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_goodN_Step2_epCDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_goodN_Step2_epCDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_goodN_Step2_epCDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_goodN_Step2_epCDn->Fill(xB, weight);

            h_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_goodN_Step2_epCDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_goodN_Step2_epCDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_goodN_Step2_epCDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_goodN_Step2_epCDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_goodN_Step2_epCDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_goodN_Step2_epCDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_goodN_Step2_epCDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_goodN_Step2_epCDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_goodN_Step2_epCDn->Fill(ToF, weight);
            h_P_n_VS_ToF_goodN_Step2_epCDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_goodN_Step2_epCDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_goodN_Step2_epCDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_goodN_Step2_epCDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_goodN_Step2_epCDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_goodN_Step2_epCDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_goodN_Step2_epCDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_goodN_Step2_epCDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_goodN_Step2_epCDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_goodN_Step2_epCDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_goodN_Step2_epCDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_goodN_Step2_epCDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_goodN_Step2_epCDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_goodN_Step2_epCDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_goodN_epCDn->Fill(xB, M_miss, weight);

            h_beta_n_goodN_Step2_epCDn->Fill(beta, weight);
        } else if (isBN) {
            h_theta_n_badN_Step2_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step2_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_badN_Step2_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
            h_theta_n_VS_beta_n_badN_Step2_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_badN_Step2_epCDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_badN_Step2_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_badN_Step2_epCDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_badN_Step2_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
            h_P_miss_VS_phi_miss_badN_Step2_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_badN_Step2_epCDn->Fill(dpp, weight);
            h_theta_n_miss_badN_Step2_epCDn->Fill(theta_n_miss, weight);

            h_E_p_badN_Step2_epCDn->Fill(E_p, weight);
            h_E_miss_badN_Step2_epCDn->Fill(E_miss, weight);
            h_M_miss_badN_Step2_epCDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_badN_Step2_epCDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_badN_Step2_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_badN_Step2_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_badN_Step2_epCDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_badN_Step2_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_badN_Step2_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_badN_Step2_epCDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_badN_Step2_epCDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_badN_Step2_epCDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_badN_Step2_epCDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_badN_Step2_epCDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_badN_Step2_epCDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_badN_Step2_epCDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_badN_Step2_epCDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_badN_Step2_epCDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_badN_Step2_epCDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_badN_Step2_epCDn->Fill(xB, weight);

            h_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_badN_Step2_epCDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_badN_Step2_epCDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_badN_Step2_epCDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_badN_Step2_epCDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_badN_Step2_epCDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_badN_Step2_epCDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_badN_Step2_epCDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_badN_Step2_epCDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_badN_Step2_epCDn->Fill(ToF, weight);
            h_P_n_VS_ToF_badN_Step2_epCDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_badN_Step2_epCDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_badN_Step2_epCDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_badN_Step2_epCDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_badN_Step2_epCDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_badN_Step2_epCDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_badN_Step2_epCDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_badN_Step2_epCDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_badN_Step2_epCDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_badN_Step2_epCDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_badN_Step2_epCDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_badN_Step2_epCDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_badN_Step2_epCDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_badN_Step2_epCDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_badN_epCDn->Fill(xB, M_miss, weight);

            h_beta_n_badN_Step2_epCDn->Fill(beta, weight);
        }
    } else if (pInFD) {
        h_dpp_allN_Step2_epFDn->Fill(dpp, weight);
        h_theta_n_miss_allN_Step2_epFDn->Fill(theta_n_miss, weight);
        h_dpp_VS_theta_n_miss_allN_Step2_epFDn->Fill(dpp, theta_n_miss, weight);

        if (theta_n_miss < 25.) {
            h_dpp_allN_for_theta_n_miss_less_than_25_Step2_epFDn->Fill(dpp, weight);
        }

        if (dpp < 0.5) {
            h_theta_n_miss_allN_for_dpp_less_than_05_Step2_epFDn->Fill(theta_n_miss, weight);

            if (dpp < 0.3) {
                h_theta_n_miss_allN_for_dpp_less_than_03_Step2_epFDn->Fill(theta_n_miss, weight);
            }
        }

        if (isGN) {
            h_theta_n_goodN_Step2_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step2_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_goodN_Step2_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
            h_theta_n_VS_beta_n_goodN_Step2_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_goodN_Step2_epFDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_goodN_Step2_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_goodN_Step2_epFDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_goodN_Step2_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
            h_P_miss_VS_phi_miss_goodN_Step2_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_goodN_Step2_epFDn->Fill(dpp, weight);
            h_theta_n_miss_goodN_Step2_epFDn->Fill(theta_n_miss, weight);

            h_E_p_goodN_Step2_epFDn->Fill(E_p, weight);
            h_E_miss_goodN_Step2_epFDn->Fill(E_miss, weight);
            h_M_miss_goodN_Step2_epFDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_goodN_Step2_epFDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_goodN_Step2_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_goodN_Step2_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_goodN_Step2_epFDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_goodN_Step2_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_goodN_Step2_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_goodN_Step2_epFDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_goodN_Step2_epFDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_goodN_Step2_epFDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_goodN_Step2_epFDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_goodN_Step2_epFDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_goodN_Step2_epFDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_goodN_Step2_epFDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_goodN_Step2_epFDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_goodN_Step2_epFDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_goodN_Step2_epFDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_goodN_Step2_epFDn->Fill(xB, weight);

            h_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_goodN_Step2_epFDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_goodN_Step2_epFDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_goodN_Step2_epFDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_goodN_Step2_epFDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_goodN_Step2_epFDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_goodN_Step2_epFDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_goodN_Step2_epFDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_goodN_Step2_epFDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_goodN_Step2_epFDn->Fill(ToF, weight);
            h_P_n_VS_ToF_goodN_Step2_epFDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_goodN_Step2_epFDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_goodN_Step2_epFDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_goodN_Step2_epFDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_goodN_Step2_epFDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_goodN_Step2_epFDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_goodN_Step2_epFDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_goodN_Step2_epFDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_goodN_Step2_epFDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_goodN_Step2_epFDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_goodN_Step2_epFDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_goodN_Step2_epFDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_goodN_Step2_epFDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_goodN_Step2_epFDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_goodN_epFDn->Fill(xB, M_miss, weight);

            h_beta_n_goodN_Step2_epFDn->Fill(beta, weight);
        } else if (isBN) {
            h_theta_n_badN_Step2_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step2_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_theta_n_VS_phi_n_badN_Step2_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
            h_theta_n_VS_beta_n_badN_Step2_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

            h_P_n_badN_Step2_epFDn->Fill(P_n_3v.Mag(), weight);
            h_P_n_VS_theta_n_badN_Step2_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

            h_P_miss_badN_Step2_epFDn->Fill(P_miss_3v.Mag(), weight);
            h_P_miss_VS_theta_miss_badN_Step2_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
            h_P_miss_VS_phi_miss_badN_Step2_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

            h_dpp_badN_Step2_epFDn->Fill(dpp, weight);
            h_theta_n_miss_badN_Step2_epFDn->Fill(theta_n_miss, weight);

            h_E_p_badN_Step2_epFDn->Fill(E_p, weight);
            h_E_miss_badN_Step2_epFDn->Fill(E_miss, weight);
            h_M_miss_badN_Step2_epFDn->Fill(M_miss, weight);
            h_M_miss_VS_P_n_badN_Step2_epFDn->Fill(P_n_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_n_badN_Step2_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_n_badN_Step2_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_P_miss_badN_Step2_epFDn->Fill(P_miss_3v.Mag(), M_miss, weight);
            h_M_miss_VS_theta_miss_badN_Step2_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, M_miss, weight);
            h_M_miss_VS_phi_miss_badN_Step2_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, M_miss, weight);

            h_P_n_minus_P_miss_badN_Step2_epFDn->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
            h_P_n_x_minus_P_miss_x_badN_Step2_epFDn->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
            h_P_n_y_minus_P_miss_y_badN_Step2_epFDn->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
            h_P_n_z_minus_P_miss_z_badN_Step2_epFDn->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

            h_P_n_VS_P_miss_badN_Step2_epFDn->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
            h_P_n_x_VS_P_miss_x_badN_Step2_epFDn->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
            h_P_n_y_VS_P_miss_y_badN_Step2_epFDn->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
            h_P_n_z_VS_P_miss_z_badN_Step2_epFDn->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

            h_theta_n_p_badN_Step2_epFDn->Fill(P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);
            h_theta_n_p_VS_P_p_badN_Step2_epFDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

            h_xB_badN_Step2_epFDn->Fill(xB, weight);

            h_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, weight);
            h_P_n_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, dpp, weight);
            h_beta_n_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, beta, weight);
            h_E_p_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, E_p, weight);
            h_E_miss_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, E_miss, weight);
            h_M_miss_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, M_miss, weight);
            h_path_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, theta_n_miss, weight);
            h_ToF_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, ToF, weight);
            h_nSector_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, nSector, weight);
            h_Edep_CND1_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND_badN_Step2_epFDn->Fill(Edep_CND, Edep_CND3, weight);

            h_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, weight);
            h_P_n_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, dpp, weight);
            h_beta_n_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, beta, weight);
            h_E_p_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, E_p, weight);
            h_E_miss_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, E_miss, weight);
            h_M_miss_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, M_miss, weight);
            h_path_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, path * 100, weight);
            h_theta_n_miss_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
            h_ToF_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, ToF, weight);
            h_nSector_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, nSector, weight);
            h_Edep_CND1_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, Edep_CND1, weight);
            h_Edep_CND2_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CTOF_badN_Step2_epFDn->Fill(Edep_CTOF, Edep_CND3, weight);

            h_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, weight);
            h_P_n_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, dpp, weight);
            h_beta_n_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, beta, weight);
            h_E_p_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, E_p, weight);
            h_E_miss_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, E_miss, weight);
            h_M_miss_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, M_miss, weight);
            h_path_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
            h_ToF_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, ToF, weight);
            h_nSector_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, nSector, weight);
            h_Edep_CND2_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Edep_CND1_badN_Step2_epFDn->Fill(Edep_CND1, Edep_CND3, weight);

            h_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, weight);
            h_P_n_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, dpp, weight);
            h_beta_n_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, beta, weight);
            h_E_p_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, E_p, weight);
            h_E_miss_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, E_miss, weight);
            h_M_miss_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, M_miss, weight);
            h_path_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
            h_ToF_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, ToF, weight);
            h_nSector_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, nSector, weight);
            h_Edep_CND3_VS_Edep_CND2_badN_Step2_epFDn->Fill(Edep_CND2, Edep_CND3, weight);

            h_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, weight);
            h_P_n_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, dpp, weight);
            h_beta_n_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, beta, weight);
            h_E_p_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, E_p, weight);
            h_E_miss_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, E_miss, weight);
            h_M_miss_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, M_miss, weight);
            h_path_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, path * 100, weight);
            h_theta_n_miss_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
            h_ToF_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, ToF, weight);
            h_nSector_VS_Edep_CND3_badN_Step2_epFDn->Fill(Edep_CND3, nSector, weight);

            h_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, weight);
            h_Edep_CND_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, Edep_CND3, weight);
            h_P_n_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, dpp, weight);
            h_beta_n_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, beta, weight);
            h_E_miss_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, E_miss, weight);
            h_M_miss_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, M_miss, weight);
            h_path_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, path * 100, weight);
            h_theta_n_miss_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, theta_n_miss, weight);
            h_ToF_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, ToF, weight);
            h_nSector_VS_Size_CND1_badN_Step2_epFDn->Fill(Size_CND1, nSector, weight);

            h_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, weight);
            h_Edep_CND_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, Edep_CND3, weight);
            h_P_n_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, dpp, weight);
            h_beta_n_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, beta, weight);
            h_E_miss_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, E_miss, weight);
            h_M_miss_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, M_miss, weight);
            h_path_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, path * 100, weight);
            h_theta_n_miss_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, theta_n_miss, weight);
            h_ToF_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, ToF, weight);
            h_nSector_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND2, nSector, weight);

            h_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, weight);
            h_Edep_CND_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, Edep_CND, weight);
            h_Edep_CND1_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, Edep_CND1, weight);
            h_Edep_CND2_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, Edep_CND2, weight);
            h_Edep_CND3_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, Edep_CND3, weight);
            h_P_n_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_P_n_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Mag(), weight);
            h_theta_n_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, dpp, weight);
            h_beta_n_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, beta, weight);
            h_E_miss_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, E_miss, weight);
            h_M_miss_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, M_miss, weight);
            h_path_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, path * 100, weight);
            h_theta_n_miss_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, theta_n_miss, weight);
            h_ToF_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, ToF, weight);
            h_nSector_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND3, nSector, weight);

            h_Size_CND1_VS_Size_CND2_badN_Step2_epFDn->Fill(Size_CND1, Size_CND2, weight);
            h_Size_CND1_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND1, Size_CND3, weight);
            h_Size_CND2_VS_Size_CND3_badN_Step2_epFDn->Fill(Size_CND2, Size_CND3, weight);

            h_ToF_badN_Step2_epFDn->Fill(ToF, weight);
            h_P_n_VS_ToF_badN_Step2_epFDn->Fill(ToF, P_n_3v.Mag(), weight);
            h_theta_n_VS_ToF_badN_Step2_epFDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_VS_ToF_badN_Step2_epFDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
            h_P_miss_VS_ToF_badN_Step2_epFDn->Fill(ToF, P_miss_3v.Mag(), weight);
            h_theta_miss_VS_ToF_badN_Step2_epFDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
            h_phi_miss_VS_ToF_badN_Step2_epFDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
            h_dpp_VS_ToF_badN_Step2_epFDn->Fill(ToF, dpp, weight);
            h_beta_n_VS_ToF_badN_Step2_epFDn->Fill(ToF, beta, weight);
            h_E_p_VS_ToF_badN_Step2_epFDn->Fill(ToF, E_p, weight);
            h_E_miss_VS_ToF_badN_Step2_epFDn->Fill(ToF, E_miss, weight);
            h_M_miss_VS_ToF_badN_Step2_epFDn->Fill(ToF, M_miss, weight);
            h_path_VS_ToF_badN_Step2_epFDn->Fill(ToF, path * 100, weight);
            h_theta_n_miss_VS_ToF_badN_Step2_epFDn->Fill(ToF, theta_n_miss, weight);
            h_nSector_VS_ToF_badN_Step2_epFDn->Fill(ToF, nSector, weight);

            h_xB_VS_M_miss_badN_epFDn->Fill(xB, M_miss, weight);

            h_beta_n_badN_Step2_epFDn->Fill(beta, weight);
        }
    }
}

// UpdateStep2Histograms2 function
// ======================================================================================================================================================================

void VetoHistograms::UpdateStep2Histograms2(bool pInCD, bool pInFD, bool isGN, bool isBN, int ldiff, int sdiff, TVector3 p_C_3v, TVector3 v_hit_3v,
                                            TVector3 P_n_3v, double dToF, double dToF_rel_pos, double dToF_rel_n, double dpp, double theta_n_miss,
                                            double Edep_CND, double beta, double path, double ToF, double weight) {
    if (pInCD) {
        if (isGN) // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step1_layer
        {
            h_sdiff_pos_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, weight);
            h_sdiff_pos_mom_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
            h_sdiff_pos_VS_VhitZ_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, v_hit_3v.Z(), weight);
            h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, ToF * c - v_hit_3v.Z(), weight);
            h_theta_n_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_sdiff_pos_VS_theta_n_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_phi_n_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_ToF_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, ToF, weight);
            h_sdiff_pos_VS_path_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, path * 100, weight);
            h_sdiff_pos_VS_beta_n_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, beta, weight);
            h_sdiff_pos_VS_Edep_CND_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, Edep_CND, weight);
            h_sdiff_pos_VS_theta_n_miss_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, theta_n_miss, weight);
            h_sdiff_pos_VS_dpp_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, dpp, weight);

            h_dToF_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(dToF, weight);
            h_sdiff_pos_VS_dToF_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, dToF, weight);
            h_dToF_rel_pos_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(dToF_rel_pos, weight);
            h_sdiff_pos_VS_dToF_rel_pos_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_pos, weight);
            h_dToF_rel_n_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(dToF_rel_n, weight);
            h_sdiff_pos_VS_dToF_rel_n_goodN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_n, weight);
        } else if (isBN) {
            h_sdiff_pos_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, weight);
            h_sdiff_pos_mom_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
            h_sdiff_pos_VS_VhitZ_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, v_hit_3v.Z(), weight);
            h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, ToF * c - v_hit_3v.Z(), weight);
            h_theta_n_badN_Step2_layer_epCDn[ldiff + 3]->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_sdiff_pos_VS_theta_n_badN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step2_layer_epCDn[ldiff + 3]->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_phi_n_badN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_ToF_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, ToF, weight);
            h_sdiff_pos_VS_path_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, path * 100, weight);
            h_sdiff_pos_VS_beta_n_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, beta, weight);
            h_sdiff_pos_VS_Edep_CND_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, Edep_CND, weight);
            h_sdiff_pos_VS_theta_n_miss_badN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, theta_n_miss, weight);
            h_sdiff_pos_VS_dpp_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, dpp, weight);

            h_dToF_badN_Step2_layer_epCDn[ldiff + 3]->Fill(dToF, weight);
            h_sdiff_pos_VS_dToF_badN_Step2_layer_epCDn[ldiff + 3]->Fill(sdiff, dToF, weight);
            h_dToF_rel_pos_badN_Step2_layer_epCDn[ldiff + 3]->Fill(dToF_rel_pos, weight);
            h_sdiff_pos_VS_dToF_rel_pos_badN_Step2_layer_epCDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_pos, weight);
            h_dToF_rel_n_badN_Step2_layer_epCDn[ldiff + 3]->Fill(dToF_rel_n, weight);
            h_sdiff_pos_VS_dToF_rel_n_badN_Step2_layer_epCDn[ldiff + 3]->
                    Fill(sdiff, dToF_rel_n, weight);
        }
    } else if (pInFD) {
        if (isGN) // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step1_layer
        {
            h_sdiff_pos_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, weight);
            h_sdiff_pos_mom_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
            h_sdiff_pos_VS_VhitZ_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, v_hit_3v.Z(), weight);
            h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, ToF * c - v_hit_3v.Z(), weight);
            h_theta_n_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_sdiff_pos_VS_theta_n_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_phi_n_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_ToF_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, ToF, weight);
            h_sdiff_pos_VS_path_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, path * 100, weight);
            h_sdiff_pos_VS_beta_n_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, beta, weight);
            h_sdiff_pos_VS_Edep_CND_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, Edep_CND, weight);
            h_sdiff_pos_VS_theta_n_miss_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, theta_n_miss, weight);
            h_sdiff_pos_VS_dpp_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, dpp, weight);

            h_dToF_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(dToF, weight);
            h_sdiff_pos_VS_dToF_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, dToF, weight);
            h_dToF_rel_pos_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(dToF_rel_pos, weight);
            h_sdiff_pos_VS_dToF_rel_pos_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_pos, weight);
            h_dToF_rel_n_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(dToF_rel_n, weight);
            h_sdiff_pos_VS_dToF_rel_n_goodN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_n, weight);
        } else if (isBN) {
            h_sdiff_pos_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, weight);
            h_sdiff_pos_mom_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
            h_sdiff_pos_VS_VhitZ_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, v_hit_3v.Z(), weight);
            h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, ToF * c - v_hit_3v.Z(), weight);
            h_theta_n_badN_Step2_layer_epFDn[ldiff + 3]->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
            h_sdiff_pos_VS_theta_n_badN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Theta() * 180. / M_PI, weight);
            h_phi_n_badN_Step2_layer_epFDn[ldiff + 3]->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_phi_n_badN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, P_n_3v.Phi() * 180. / M_PI, weight);
            h_sdiff_pos_VS_ToF_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, ToF, weight);
            h_sdiff_pos_VS_path_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, path * 100, weight);
            h_sdiff_pos_VS_beta_n_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, beta, weight);
            h_sdiff_pos_VS_Edep_CND_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, Edep_CND, weight);
            h_sdiff_pos_VS_theta_n_miss_badN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, theta_n_miss, weight);
            h_sdiff_pos_VS_dpp_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, dpp, weight);

            h_dToF_badN_Step2_layer_epFDn[ldiff + 3]->Fill(dToF, weight);
            h_sdiff_pos_VS_dToF_badN_Step2_layer_epFDn[ldiff + 3]->Fill(sdiff, dToF, weight);
            h_dToF_rel_pos_badN_Step2_layer_epFDn[ldiff + 3]->Fill(dToF_rel_pos, weight);
            h_sdiff_pos_VS_dToF_rel_pos_badN_Step2_layer_epFDn[ldiff + 3]->Fill(
                sdiff, dToF_rel_pos, weight);
            h_dToF_rel_n_badN_Step2_layer_epFDn[ldiff + 3]->Fill(dToF_rel_n, weight);
            h_sdiff_pos_VS_dToF_rel_n_badN_Step2_layer_epFDn[ldiff + 3]->
                    Fill(sdiff, dToF_rel_n, weight);
        }
    }
}

// UpdateMultiplicityHistograms function
// ======================================================================================================================================================================

void VetoHistograms::UpdateMultiplicityHistograms(bool pInCD, bool pInFD, bool isGN, bool isBN,
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
                                                  double weight) {
    if (pInCD) {
        h_n_multiplicity_allN_epCDn_Step0->Fill(counter_n_multiplicity_allN_epCDn_Step0, weight);
        h_n_multiplicity_goodN_epCDn_Step0->Fill(counter_n_multiplicity_goodN_epCDn_Step0, weight);
        h_n_multiplicity_badN_epCDn_Step0->Fill(counter_n_multiplicity_badN_epCDn_Step0, weight);

        h_n_multiplicity_allN_epCDn_Step1->Fill(counter_n_multiplicity_allN_epCDn_Step1, weight);
        h_n_multiplicity_goodN_epCDn_Step1->Fill(counter_n_multiplicity_goodN_epCDn_Step1, weight);
        h_n_multiplicity_badN_epCDn_Step1->Fill(counter_n_multiplicity_badN_epCDn_Step1, weight);

        h_n_multiplicity_allN_epCDn_Step2->Fill(counter_n_multiplicity_allN_epCDn_Step2, weight);
        h_n_multiplicity_goodN_epCDn_Step2->Fill(counter_n_multiplicity_goodN_epCDn_Step2, weight);
        h_n_multiplicity_badN_epCDn_Step2->Fill(counter_n_multiplicity_badN_epCDn_Step2, weight);

        h_n_multiplicity_allN_epCDn_Step3->Fill(counter_n_multiplicity_allN_epCDn_Step3, weight);
        h_n_multiplicity_goodN_epCDn_Step3->Fill(counter_n_multiplicity_goodN_epCDn_Step3, weight);
        h_n_multiplicity_badN_epCDn_Step3->Fill(counter_n_multiplicity_badN_epCDn_Step3, weight);

        h_n_multiplicity_allN_epCDn_Step4->Fill(counter_n_multiplicity_allN_epCDn_Step4, weight);
        h_n_multiplicity_goodN_epCDn_Step4->Fill(counter_n_multiplicity_goodN_epCDn_Step4, weight);
        h_n_multiplicity_badN_epCDn_Step4->Fill(counter_n_multiplicity_badN_epCDn_Step4, weight);

        h_n_multiplicity_allN_epCDn_Step5->Fill(counter_n_multiplicity_allN_epCDn_Step5, weight);
        h_n_multiplicity_goodN_epCDn_Step5->Fill(counter_n_multiplicity_goodN_epCDn_Step5, weight);
        h_n_multiplicity_badN_epCDn_Step5->Fill(counter_n_multiplicity_badN_epCDn_Step5, weight);
    } else if (pInFD) {
        h_n_multiplicity_allN_epFDn_Step0->Fill(counter_n_multiplicity_allN_epFDn_Step0, weight);
        h_n_multiplicity_goodN_epFDn_Step0->Fill(counter_n_multiplicity_goodN_epFDn_Step0, weight);
        h_n_multiplicity_badN_epFDn_Step0->Fill(counter_n_multiplicity_badN_epFDn_Step0, weight);

        h_n_multiplicity_allN_epFDn_Step1->Fill(counter_n_multiplicity_allN_epFDn_Step1, weight);
        h_n_multiplicity_goodN_epFDn_Step1->Fill(counter_n_multiplicity_goodN_epFDn_Step1, weight);
        h_n_multiplicity_badN_epFDn_Step1->Fill(counter_n_multiplicity_badN_epFDn_Step1, weight);

        h_n_multiplicity_allN_epFDn_Step2->Fill(counter_n_multiplicity_allN_epFDn_Step2, weight);
        h_n_multiplicity_goodN_epFDn_Step2->Fill(counter_n_multiplicity_goodN_epFDn_Step2, weight);
        h_n_multiplicity_badN_epFDn_Step2->Fill(counter_n_multiplicity_badN_epFDn_Step2, weight);

        h_n_multiplicity_allN_epFDn_Step3->Fill(counter_n_multiplicity_allN_epFDn_Step3, weight);
        h_n_multiplicity_goodN_epFDn_Step3->Fill(counter_n_multiplicity_goodN_epFDn_Step3, weight);
        h_n_multiplicity_badN_epFDn_Step3->Fill(counter_n_multiplicity_badN_epFDn_Step3, weight);

        h_n_multiplicity_allN_epFDn_Step4->Fill(counter_n_multiplicity_allN_epFDn_Step4, weight);
        h_n_multiplicity_goodN_epFDn_Step4->Fill(counter_n_multiplicity_goodN_epFDn_Step4, weight);
        h_n_multiplicity_badN_epFDn_Step4->Fill(counter_n_multiplicity_badN_epFDn_Step4, weight);

        h_n_multiplicity_allN_epFDn_Step5->Fill(counter_n_multiplicity_allN_epFDn_Step5, weight);
        h_n_multiplicity_goodN_epFDn_Step5->Fill(counter_n_multiplicity_goodN_epFDn_Step5, weight);
        h_n_multiplicity_badN_epFDn_Step5->Fill(counter_n_multiplicity_badN_epFDn_Step5, weight);
    }
}
