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

    TH2D *h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epCDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epCDn[7];

    TH2D *h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epFDn[7];
    TH2D *h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epFDn[7];

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

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name, "sdiff_pos_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Counts", k - 3);
        h_sdiff_pos_goodN_Step2prep_layer_epCDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Counts", k - 3);
        h_sdiff_pos_badN_Step2prep_layer_epCDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Counts", k - 3);
        h_sdiff_pos_goodN_Step2prep_layer_epFDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Counts", k - 3);
        h_sdiff_pos_badN_Step2prep_layer_epFDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_mom_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Momentum Proton [GeV/c]", k - 3);
        h_sdiff_pos_mom_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0.3, 1.5);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_mom_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Momentum Proton [GeV/c]", k - 3);
        h_sdiff_pos_mom_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0.3, 1.5);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_mom_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Momentum Proton [GeV/c]", k - 3);
        h_sdiff_pos_mom_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0.4, 3.);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_mom_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Momentum Proton [GeV/c]", k - 3);
        h_sdiff_pos_mom_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0.4, 3.);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. V_{hit,z} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};V_{hit,z} [cm]", k - 3);
        h_sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_VhitZ_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. V_{hit,z} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};V_{hit,z} [cm]", k - 3);
        h_sdiff_pos_VS_VhitZ_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_VS_VhitZ_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. V_{hit,z} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};V_{hit,z} [cm]", k - 3);
        h_sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_VS_VhitZ_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_VhitZ_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. V_{hit,z} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};V_{hit,z} [cm]", k - 3);
        h_sdiff_pos_VS_VhitZ_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_VS_VhitZ_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-V_{hit,z} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};ToF*c-V_{hit,z} [cm]", k - 3);
        h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-V_{hit,z} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};ToF*c-V_{hit,z} [cm]", k - 3);
        h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-V_{hit,z} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};ToF*c-V_{hit,z} [cm]", k - 3);
        h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_VS_ToF_c_minus_VhitZ_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-V_{hit,z} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};ToF*c-V_{hit,z} [cm]", k - 3);
        h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_VS_ToF_c_minus_VhitZ_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_VS_theta_n_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #theta_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#theta_{n} [#circ]", k - 3);
        h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 180);
        HistoList.push_back(h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_theta_n_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #theta_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#theta_{n} [#circ]", k - 3);
        h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 180);
        HistoList.push_back(h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_VS_theta_n_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #theta_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#theta_{n} [#circ]", k - 3);
        h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 180);
        HistoList.push_back(h_sdiff_pos_VS_theta_n_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_theta_n_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #theta_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#theta_{n} [#circ]", k - 3);
        h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 180);
        HistoList.push_back(h_sdiff_pos_VS_theta_n_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_VS_phi_n_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #phi_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#phi_{n} [#circ]", k - 3);
        h_sdiff_pos_VS_phi_n_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -180, 180);
        HistoList.push_back(h_sdiff_pos_VS_phi_n_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_phi_n_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #phi_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#phi_{n} [#circ]", k - 3);
        h_sdiff_pos_VS_phi_n_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -180, 180);
        HistoList.push_back(h_sdiff_pos_VS_phi_n_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_VS_phi_n_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #phi_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#phi_{n} [#circ]", k - 3);
        h_sdiff_pos_VS_phi_n_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -180, 180);
        HistoList.push_back(h_sdiff_pos_VS_phi_n_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_phi_n_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #phi_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#phi_{n} [#circ]", k - 3);
        h_sdiff_pos_VS_phi_n_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -180, 180);
        HistoList.push_back(h_sdiff_pos_VS_phi_n_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_VS_ToF_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Neutron ToF (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};t_{ToF}^{n} [ns]", k - 3);
        h_sdiff_pos_VS_ToF_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0., 20.);
        HistoList.push_back(h_sdiff_pos_VS_ToF_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_ToF_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Neutron ToF (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};t_{ToF}^{n} [ns]", k - 3);
        h_sdiff_pos_VS_ToF_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0., 20.);
        HistoList.push_back(h_sdiff_pos_VS_ToF_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_VS_ToF_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Neutron ToF (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};t_{ToF}^{n} [ns]", k - 3);
        h_sdiff_pos_VS_ToF_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0., 20.);
        HistoList.push_back(h_sdiff_pos_VS_ToF_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_ToF_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Neutron ToF (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};t_{ToF}^{n} [ns]", k - 3);
        h_sdiff_pos_VS_ToF_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0., 20.);
        HistoList.push_back(h_sdiff_pos_VS_ToF_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_VS_path_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Neutron path length (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Path length [cm]", k - 3);
        h_sdiff_pos_VS_path_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 20., 60.);
        HistoList.push_back(h_sdiff_pos_VS_path_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_path_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Neutron path length (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Path length [cm]", k - 3);
        h_sdiff_pos_VS_path_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 20., 60.);
        HistoList.push_back(h_sdiff_pos_VS_path_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_VS_path_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Neutron path length (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Path length [cm]", k - 3);
        h_sdiff_pos_VS_path_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 20., 60.);
        HistoList.push_back(h_sdiff_pos_VS_path_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_path_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Neutron path length (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};Path length [cm]", k - 3);
        h_sdiff_pos_VS_path_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 20., 60.);
        HistoList.push_back(h_sdiff_pos_VS_path_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_VS_beta_n_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #beta_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#beta_{n}", k - 3);
        h_sdiff_pos_VS_beta_n_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -0.1, 1.1);
        HistoList.push_back(h_sdiff_pos_VS_beta_n_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_beta_n_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #beta_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#beta_{n}", k - 3);
        h_sdiff_pos_VS_beta_n_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -0.1, 1.1);
        HistoList.push_back(h_sdiff_pos_VS_beta_n_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_VS_beta_n_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #beta_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#beta_{n}", k - 3);
        h_sdiff_pos_VS_beta_n_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -0.1, 1.1);
        HistoList.push_back(h_sdiff_pos_VS_beta_n_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_beta_n_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. #beta_{n} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};#beta_{n}", k - 3);
        h_sdiff_pos_VS_beta_n_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -0.1, 1.1);
        HistoList.push_back(h_sdiff_pos_VS_beta_n_badN_Step2prep_layer_epFDn[k]);


        sprintf(temp_name, "sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. E^{CND}_{dep} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};E^{CND}_{dep} [MeV]", k - 3);
        h_sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -0.1, 1.1);
        HistoList.push_back(h_sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. E^{CND}_{dep} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};E^{CND}_{dep} [MeV]", k - 3);
        h_sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -0.1, 1.1);
        HistoList.push_back(h_sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. E^{CND}_{dep} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};E^{CND}_{dep} [MeV]", k - 3);
        h_sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -0.1, 1.1);
        HistoList.push_back(h_sdiff_pos_VS_Edep_CND_goodN_Step2prep_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. E^{CND}_{dep} (#DeltaL_{n,+} = %d);#DeltaS_{n,+} = S_{n} - S_{+};E^{CND}_{dep} [MeV]", k - 3);
        h_sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -0.1, 1.1);
        HistoList.push_back(h_sdiff_pos_VS_Edep_CND_badN_Step2prep_layer_epFDn[k]);
    }