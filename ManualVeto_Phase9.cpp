#include <cstdlib>
#include <iostream>

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
#include "HipoChain.h"
#include "clas12reader.h"
//
#include "src/classes/VetoHistograms/UpdateHistograms.cpp"
#include "src/classes/clas12ana/clas12ana.cpp"
#include "src/constants.h"
#include "src/functions/Andrews_functions/Andrews_functions.cpp"
#include "src/functions/GeneralFunctions.h"
#include "src/functions/HipoChain_config.cpp"
#include "src/functions/NeutronFunctions.h"
#include "src/functions/neutron-veto/veto_functions.cpp"

using namespace std;
using namespace clas12;

#pragma region /* ManualVeto_Phase9 - start */

int ManualVeto_Phase9(                            //
    const string OutDir, string output_pdf_Erin,  // My arguments
    double Ebeam, bool keep_good, string output_root_Erin, string output_txt_Erin, string input_hipo,
    // Erin's arguments
    string PDFFile, int isMC = 0  // Andrew's arguments
) {
    auto Code_start_time = std::chrono::system_clock::now();  // Start counting running time

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=
    // Printouts
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=

#pragma region /* Printouts 1 - start */

    cout << "\033[33m\n\033[0m";
    cout << "\033[33minput_hipo:\033[0m\t\t" << input_hipo << "\n";
    cout << "\033[33m\n\033[0m";
    cout << "\033[33mOUTDIR:\033[0m\t\t\t" << gSystem->Getenv("OUTDIR") << "\n";
    cout << "\033[33mOutDir:\033[0m\t\t\t" << OutDir << "\n";
    cout << "\033[33moutput_pdf_Erin:\033[0m\t" << output_pdf_Erin << "\n";
    cout << "\033[33moutput_root_Erin:\033[0m\t" << output_root_Erin << "\n";
    cout << "\033[33moutput_txt_Erin:\033[0m\t" << output_txt_Erin << "\n";
    cout << "\033[33mPDFFile:\033[0m\t\t" << PDFFile << "\n\n";

#pragma endregion /* Printouts 1 - end */

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=
    // Initial setup
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=

#pragma region /* Initial setup - start */

    // Delete old output folder
    cout << "\033[33m\nClearing\033[0m '" << OutDir << "'\n";
    system(("rm -r " + OutDir).c_str());
    cout << "\n";

    // Remake old output folder
    cout << "\033[33m\nRemaking\033[0m '" << OutDir << "'\n";
    system(("mkdir -p " + OutDir).c_str());
    cout << "\n\n";

    // Erin's output file names
    TFile *f = new TFile(output_root_Erin.c_str(), "RECREATE");
    TTree *ntree = new TTree("T", "NeutronTree");
    std::ofstream outtxt(output_txt_Erin);

    // Input hipo file
    clas12root::HipoChain chain;
    HipoChain_config(chain, input_hipo);

    auto config_c12 = chain.GetC12Reader();
    chain.SetReaderTags({0});
    const std::unique_ptr<clas12::clas12reader> &c12 = chain.C12ref();
    chain.db()->turnOffQADB();

    int numevent = 0;

    // Set up root tree for TMVA
    Int_t nhits;
    double px, py, pz, momentum;
    Int_t sec[100] = {-1};
    Int_t lay[100] = {-1};
    int event;
    double energy, cnd_energy, ctof_energy, angle_diff;
    int layermult, size, cnd_hits, ctof_hits;
    bool is_CTOF, is_CND1, is_CND2, is_CND3;

    int counter = 0;

    // set up instance of clas12ana
    clas12ana *clasAna = new clas12ana();

    clasAna->readEcalSFPar("src/cuts/paramsSF_LD2_x2.dat");  // TODO: check if applied
    clasAna->readEcalPPar("src/cuts/paramsPI_LD2_x2.dat");   // TODO: check if applied

    clasAna->setProtonPidCuts(true);

#pragma endregion /* Initial setup - end */

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=
    // Veto histograms
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=

    VetoHistograms histograms = VetoHistograms();

#pragma endregion /* Veto histograms - end */

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=
    // Chain loop
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=

#pragma region /* Chain loop - start */

    int EventCounter = 0; /* From Andrew */

    int counter_epXn = 0;
    int counter_pass_step0_cuts = 0, counter_pass_step1_cuts = 0, counter_pass_step2_cuts = 0, counter_pass_step3_cuts = 0, counter_pass_step4_cuts = 0, counter_pass_step5_cuts = 0;
    int counter_n_multiplicity_allN_epCDn = 0, counter_n_multiplicity_goodN_epCDn = 0, counter_n_multiplicity_badN_epCDn = 0;
    int counter_n_multiplicity_allN_epCDn_Step0 = 0, counter_n_multiplicity_goodN_epCDn_Step0 = 0, counter_n_multiplicity_badN_epCDn_Step0 = 0;
    int counter_n_multiplicity_allN_epCDn_Step1 = 0, counter_n_multiplicity_goodN_epCDn_Step1 = 0, counter_n_multiplicity_badN_epCDn_Step1 = 0;
    int counter_n_multiplicity_allN_epCDn_Step2 = 0, counter_n_multiplicity_goodN_epCDn_Step2 = 0, counter_n_multiplicity_badN_epCDn_Step2 = 0;
    int counter_n_multiplicity_allN_epCDn_Step3 = 0, counter_n_multiplicity_goodN_epCDn_Step3 = 0, counter_n_multiplicity_badN_epCDn_Step3 = 0;
    int counter_n_multiplicity_allN_epCDn_Step4 = 0, counter_n_multiplicity_goodN_epCDn_Step4 = 0, counter_n_multiplicity_badN_epCDn_Step4 = 0;
    int counter_n_multiplicity_allN_epCDn_Step5 = 0, counter_n_multiplicity_goodN_epCDn_Step5 = 0, counter_n_multiplicity_badN_epCDn_Step5 = 0;
    int counter_n_multiplicity_allN_epFDn = 0, counter_n_multiplicity_goodN_epFDn = 0, counter_n_multiplicity_badN_epFDn = 0;
    int counter_n_multiplicity_allN_epFDn_Step0 = 0, counter_n_multiplicity_goodN_epFDn_Step0 = 0, counter_n_multiplicity_badN_epFDn_Step0 = 0;
    int counter_n_multiplicity_allN_epFDn_Step1 = 0, counter_n_multiplicity_goodN_epFDn_Step1 = 0, counter_n_multiplicity_badN_epFDn_Step1 = 0;
    int counter_n_multiplicity_allN_epFDn_Step2 = 0, counter_n_multiplicity_goodN_epFDn_Step2 = 0, counter_n_multiplicity_badN_epFDn_Step2 = 0;
    int counter_n_multiplicity_allN_epFDn_Step3 = 0, counter_n_multiplicity_goodN_epFDn_Step3 = 0, counter_n_multiplicity_badN_epFDn_Step3 = 0;
    int counter_n_multiplicity_allN_epFDn_Step4 = 0, counter_n_multiplicity_goodN_epFDn_Step4 = 0, counter_n_multiplicity_badN_epFDn_Step4 = 0;
    int counter_n_multiplicity_allN_epFDn_Step5 = 0, counter_n_multiplicity_goodN_epFDn_Step5 = 0, counter_n_multiplicity_badN_epFDn_Step5 = 0;

    while (chain.Next()) {
        // Display completed (from Andrew)
        EventCounter++;

        if ((EventCounter % 1000000) == 0) {
            cerr << "\n\n";
            cerr << "\033[33m" << EventCounter / 1000000 << " million completed\033[0m\n\n";
        }

#pragma region /* PID & variable definitions - start */

        // PID
        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------=

        // PID (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        clasAna->Run(c12);

        auto Electrons = clasAna->getByPid(11);
        auto Protons = clasAna->getByPid(2212);
        auto Neutrons = clasAna->getByPid(2112);

        auto AllParticles = c12->getDetParticles();

        // Event selection
        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------=

        // Event selection (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // One electron in event:
        if (Electrons.size() != Num_of_e_cut) { continue; }

        // One proton in event:
        if (Protons.size() != Num_of_p_cut) { continue; }

        // At least one neutron in event:
        if (Neutrons.size() < 1) { continue; }

        // Reject particles with the wrong PID
        bool trash = 0;

        for (int i = 0; i < AllParticles.size(); i++) {
            int pid = AllParticles[i]->par()->getPid();

            if (pid != 2112 && pid != 11 && pid != 2212 && pid != 0 && pid != 22) { trash = 1; }
        }

        if (trash == 1) { continue; }

        ++counter_epXn;
        numevent = numevent + 1;

        // Variable definitions
        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------=

        // Variable definitions (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        event = c12->runconfig()->getEvent() << '\n';

        double starttime = c12->event()->getStartTime();

        TVector3 P_b_3v(0, 0, Ebeam);

        // Variable definitions (from Andrew)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        double weight = 1;

        if (isMC) { weight = c12->mcevent()->getWeight(); }

#pragma endregion /* PID & variable definitions - end */

#pragma region /* Electrons - start */

        // Electrons (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        TVector3 P_e_3v(0., 0., 0.);

        double P_e_x = Electrons[0]->par()->getPx();
        double P_e_y = Electrons[0]->par()->getPy();
        double P_e_z = Electrons[0]->par()->getPz();

        P_e_3v.SetXYZ(P_e_x, P_e_y, P_e_z);

        double Vz_e = Electrons[0]->par()->getVz();

        TVector3 P_q_3v = P_b_3v - P_e_3v;      // 3-momentum transfer
        double nu = Ebeam - P_e_3v.Mag();       // Energy transfer
        double Q2 = P_q_3v.Mag2() - (nu * nu);  // 4-momentum transfer squared
        double xB = Q2 / (2 * mN * nu);         // x Bjorken

        // Electrons (from Andrew)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        double EoP_e = (Electrons[0]->cal(PCAL)->getEnergy() + Electrons[0]->cal(ECIN)->getEnergy() + Electrons[0]->cal(ECOUT)->getEnergy()) / P_e_3v.Mag();
        int nphe = Electrons[0]->che(HTCC)->getNphe();

        int e_sector = Electrons[0]->getSector();

        double theta_q = P_q_3v.Theta() * 180 / M_PI;
        double WSq = (mN * mN) - Q2 + (2 * nu * mN);  // Hadronic mass
        double theta_e = P_e_3v.Theta() * 180 / M_PI;

#pragma endregion /* Electrons - end */

#pragma region /* Protons - start */

        // Protons (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        int counter_pCD_multiplicity_BPID = 0, counter_pFD_multiplicity_BPID = 0;
        int counter_pCD_multiplicity_APID = 0, counter_pFD_multiplicity_APID = 0;

        int p_index = -1;

        TVector3 P_p_3v(0., 0., 0.);

        // Technically not optimized - this doesn't address what happens if there are two protons passing cuts
        // TODO: recheck this!
        for (int i = 0; i < Protons.size(); i++) {
            // define quantities
            P_p_3v.SetMagThetaPhi(Protons[i]->getP(), Protons[i]->getTheta(), Protons[i]->getPhi());
            double dbeta = Protons[i]->par()->getBeta() - P_p_3v.Mag() / sqrt(P_p_3v.Mag2() + mP * mP);
            double p_theta = P_p_3v.Theta() * 180. / M_PI;
            double Vz_p = Protons[i]->par()->getVz();
            double chipid = Protons[i]->par()->getChi2Pid();

            if (Protons[i]->getRegion() == CD) {
                ++counter_pCD_multiplicity_BPID;

                histograms.UpdateBPIDpCDHistograms(P_p_3v, p_theta, dbeta, Vz_p, Vz_e, chipid, weight);

                if (fabs(Vz_p - Vz_e) > dVz_pCD_cut) { continue; }

                if (P_p_3v.Mag() < P_pCD_lcut || P_p_3v.Mag() > P_pCD_ucut) { continue; }

                // Moving to chi2 instead of beta cuts
                // if (fabs(dbeta) > 0.05) {
                if (chipid < pCD_chi2_lcut || chipid > pFD_chi2_ucut) { continue; }

                ++counter_pCD_multiplicity_APID;

                histograms.UpdateAPIDpCDHistograms(P_p_3v, p_theta, dbeta, Vz_p, Vz_e, chipid, weight);
            } else if (Protons[i]->getRegion() == FD) {
                ++counter_pFD_multiplicity_BPID;

                histograms.UpdateBPIDpFDHistograms(P_p_3v, p_theta, dbeta, Vz_p, Vz_e, chipid, weight);

                if (fabs(Vz_p - Vz_e) > dVz_pFD_cut) { continue; }

                if (P_p_3v.Mag() < P_pFD_lcut || P_p_3v.Mag() > P_pFD_ucut) { continue; }

                // Moving to chi2 instead of beta cuts
                // if (fabs(dbeta) > 0.03) {
                if (chipid < pFD_chi2_lcut || chipid > pFD_chi2_ucut) { continue; }

                ++counter_pFD_multiplicity_APID;

                histograms.UpdateAPIDpFDHistograms(P_p_3v, p_theta, dbeta, Vz_p, Vz_e, chipid, weight);
            }

            p_index = i;
        }

        histograms.UpdateProtonMultiBCHistograms(counter_pCD_multiplicity_BPID, counter_pFD_multiplicity_BPID, weight);

        if (p_index < 0) { continue; }

        histograms.UpdateProtonMultiACHistograms(counter_pCD_multiplicity_APID, counter_pFD_multiplicity_APID, weight);

        P_p_3v.SetMagThetaPhi(Protons[p_index]->getP(), Protons[p_index]->getTheta(), Protons[p_index]->getPhi());

        double beta_p = Protons[p_index]->par()->getBeta();

        // Determin where is the proton. Moved from angle cuts to getRegion() by the advice of Andrew.
        bool pInFD = (Protons[p_index]->getRegion() == FD);  // My addition
        bool pInCD = (Protons[p_index]->getRegion() == CD);  // My addition

#pragma endregion /* Protons - end */

#pragma region /* Missing momentum - start */

        // Missing momentum (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Missing momentum, energy, mass
        TVector3 P_miss_3v = P_q_3v - P_p_3v;

        momentum = P_miss_3v.Mag();

        double E_p = sqrt(mN * mN + P_p_3v.Mag2());
        double E_miss = Ebeam + mD - P_e_3v.Mag() - E_p;
        double M_miss = sqrt((E_miss * E_miss) - P_miss_3v.Mag2());

#pragma endregion /* Missing momentum - end */

        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // Andrew's manual work
        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------

#pragma region /* Andrew's manual work - start */

#pragma region /* Missing momentum cuts (Andrew) - start */

        //////////////////////////////////////////////
        // Missing momentum cuts
        //////////////////////////////////////////////

        histograms.UpdateBmissCHistograms(pInCD, pInFD, P_miss_3v, E_p, E_miss, M_miss, xB, weight);

        if (P_miss_3v.Mag() < P_miss_lcut || P_miss_3v.Mag() > P_miss_ucut) { continue; }

        if (P_miss_3v.Theta() * 180 / M_PI < Theta_miss_lcut || P_miss_3v.Theta() * 180 / M_PI > Theta_miss_ucut) { continue; }

        if (M_miss < M_miss_lcut || M_miss > M_miss_ucut) { continue; }

        histograms.UpdateAmissCHistograms(pInCD, pInFD, P_miss_3v, E_p, E_miss, M_miss, xB, weight);

#pragma endregion /* Missing momentum cuts (Andrew) - end */

#pragma region /* Neutrons (Andrew) */

        bool pass_step0_cuts = false, pass_step1_cuts = false, pass_step2_cuts = false, pass_step3_cuts = false, pass_step4_cuts = false, pass_step5_cuts = false;

        for (int itr1 = 0; itr1 < AllParticles.size(); itr1++) {
            // Cut out charged particles:
            if (AllParticles[itr1]->par()->getCharge() != 0) { continue; }

            // Why this cut? reco code bug. Neutrons in this angle range are in the BAND and appear in the CND.
            // This bug is probably fixed, yet the cut is still applied to mak sure.
            if (AllParticles[itr1]->getTheta() * 180 / M_PI > Theta_n_ucut) { continue; }

            // Andrew's response checks:
            bool CT = (AllParticles[itr1]->sci(CTOF)->getDetector() == 4);
            bool C1 = (AllParticles[itr1]->sci(CND1)->getDetector() == 3);
            bool C2 = (AllParticles[itr1]->sci(CND2)->getDetector() == 3);
            bool C3 = (AllParticles[itr1]->sci(CND3)->getDetector() == 3);

            // // Erin's response checks:
            // bool CT = (AllParticles[itr1]->sci(CTOF)->getDetector() == 4);
            // bool C1 = (AllParticles[itr1]->sci(CND1)->getLayer() == 1);
            // bool C2 = (AllParticles[itr1]->sci(CND2)->getLayer() == 2);
            // bool C3 = (AllParticles[itr1]->sci(CND3)->getLayer() == 3);

            // Cut out neutrons without a CND hit in one of its layers:
            if (!(C1 || C2 || C3)) { continue; }

            // Use CTOF as a veto for charged particles:
            if (CT) { continue; }

            // Explicit calculation of the neutron's momentum (to bypass cases where P_n is E_dep)
            double theta = AllParticles[itr1]->getTheta() * 180 / M_PI;
            double beta = AllParticles[itr1]->par()->getBeta();
            double gamma = 1 / sqrt(1 - (beta * beta));
            double mom = gamma * beta * mN;
            double ToF = AllParticles[itr1]->getTime() - starttime;

            int detINTlayer = C1 ? 1 : C2 ? 2 : 3;
            auto detlayer = C1 ? CND1 : C2 ? CND2 : CND3;  // CND layer with hit

            int Size_CND1 = AllParticles[itr1]->sci(CND1)->getSize();
            int Size_CND2 = AllParticles[itr1]->sci(CND2)->getSize();
            int Size_CND3 = AllParticles[itr1]->sci(CND3)->getSize();
            int Size_CND = Size_CND1 + Size_CND2 + Size_CND3;

            int LayerMult_CND1 = AllParticles[itr1]->sci(CND1)->getLayermulti();
            int LayerMult_CND2 = AllParticles[itr1]->sci(CND2)->getLayermulti();
            int LayerMult_CND3 = AllParticles[itr1]->sci(CND3)->getLayermulti();
            int LayerMult_CND = LayerMult_CND1 + LayerMult_CND2 + LayerMult_CND3;

            double Edep_CND1 = AllParticles[itr1]->sci(CND1)->getEnergy();
            double Edep_CND2 = AllParticles[itr1]->sci(CND2)->getEnergy();
            double Edep_CND3 = AllParticles[itr1]->sci(CND3)->getEnergy();
            double Edep_CND = Edep_CND1 + Edep_CND2 + Edep_CND3;

            double Edep_single = AllParticles[itr1]->sci(detlayer)->getEnergy();

            double Edep_CTOF = AllParticles[itr1]->sci(CTOF)->getEnergy();

            double nvtx_x = AllParticles[itr1]->par()->getVx();
            double nvtx_y = AllParticles[itr1]->par()->getVy();
            double nvtx_z = AllParticles[itr1]->par()->getVz();
            TVector3 v_nvtx_3v(nvtx_x, nvtx_y, nvtx_z);  // Neutron's vertex location

            TVector3 v_hit_3v;  // Neutron's hit location in CND
            v_hit_3v.SetXYZ(AllParticles[itr1]->sci(detlayer)->getX(), AllParticles[itr1]->sci(detlayer)->getY(), AllParticles[itr1]->sci(detlayer)->getZ());

            TVector3 v_path_3v = v_hit_3v - v_nvtx_3v;  // Direct calculation of neutron's path (in vector form)
            TVector3 P_n_3v;
            P_n_3v.SetMagThetaPhi(mom, v_path_3v.Theta(), v_path_3v.Phi());

            // Why "v_path_3v.Mag() / 100"? unit conversion.
            // TODO: check if this unit conversion is needed!
            double path = v_path_3v.Mag() / 100;
            // double path = v_path_3v.Mag();
            double theta_n_miss = P_n_3v.Angle(P_miss_3v) * 180 / M_PI;
            // Opening angle between calculated neutron's momentum and predicted neutron momentum (= missing momentum)
            double dpp = (P_miss_3v.Mag() - P_n_3v.Mag()) / P_miss_3v.Mag();
            int nSector = AllParticles[itr1]->sci(detlayer)->getSector();
            // Number of CND sector with a neutron hit in the layer detlayer

#pragma region /* Neutron PID cuts - start */

            // Beta cut:
            // Upper: beta > 0.8 -> cut out photons
            // Lower: beta < 0.15 ->
            if (beta < Beta_n_lcut || beta > Beta_n_ucut) { continue; }

            // Why this cut? reco code bug. Neutrons in this angle range are in the BAND and appear in the CND.
            // This bug is probobly fixed, yet the cut is still applied to mak sure.
            // Updated: now is forcing the neutron to be inside the acceptance of the CD
            if ((P_n_3v.Theta() * 180. / M_PI < Theta_n_lcut) || (P_n_3v.Theta() * 180. / M_PI > Theta_n_ucut)) { continue; }

            // Status cut for double-hits (based on Erin's code)
            int status = 0;
            if (C1) { status = status + AllParticles[itr1]->sci(CND1)->getStatus(); }
            if (C2) { status = status + AllParticles[itr1]->sci(CND3)->getStatus(); }
            if (C3) { status = status + AllParticles[itr1]->sci(CND2)->getStatus(); }
            if (status != 0) { continue; }
            // if ((AllParticles[itr1]->sci(CND1)->getStatus() + AllParticles[itr1]->sci(CND2)->getStatus() + AllParticles[itr1]->sci(CND3)->getStatus()) != Status_n_cut) {
            //     continue;
            // }

#pragma endregion /* Neutron PID cuts - end */

            // Check to see if there is a good neutron
            bool isGN = false;
            bool isBN = false;

            // Good neutron definition:
            bool GN_theta_n_miss = (theta_n_miss <= GN_theta_n_miss_ucut);
            bool GN_dpp = ((dpp >= GN_dpp_lcut) && (dpp <= GN_dpp_ucut));
            if (GN_theta_n_miss && GN_dpp) { isGN = true; }

            // Bad neutron definition:
            bool BN_theta_n_miss = (theta_n_miss >= BN_theta_n_miss_lcut);
            bool BN_dpp = (dpp <= BN_dpp_ucut);
            if (BN_theta_n_miss || BN_dpp) { isBN = true; }
            // if (BN_theta_n_miss && BN_dpp) { isBN = true; }
            // if (!((theta_n_miss < 25.) && ((dpp > -0.3) && (dpp < 0.3)))) { isBN = true; }

            if (isGN && isBN) { cout << "\nERROR! good and bad neutrons are overlapping! Aborting...\n", exit(0); }

            // if (!(isGN || isBN)) { continue; }

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn, counter_n_multiplicity_goodN_epCDn, counter_n_multiplicity_badN_epCDn, counter_n_multiplicity_allN_epFDn,
                               counter_n_multiplicity_goodN_epFDn, counter_n_multiplicity_badN_epFDn);

            // FILL HISTOS FOR NEUTRON CANDIDATES
            histograms.UpdatePreStepHistograms(pInCD, pInFD, isGN, isBN, P_q_3v, Q2, P_p_3v, P_miss_3v, P_n_3v, E_p, E_miss, M_miss, xB, dpp, theta_n_miss, Edep_CND, Edep_CND1, Edep_CND2, Edep_CND3,
                                               Edep_CTOF, nSector, Size_CND1, Size_CND2, Size_CND3, LayerMult_CND1, LayerMult_CND2, LayerMult_CND3, beta, beta_p, path, ToF, weight);

            if (!(isGN || isBN)) { continue; }

            //////////////////////////////////////////////
            // Step Zero
            //////////////////////////////////////////////

#pragma region /* Step Zero - start */

            /* Fill BS0C plots */
            histograms.UpdateBS0CHistograms(pInCD, pInFD, P_n_3v, v_hit_3v, beta, path, ToF, weight);

            // Why "path * 100"? unit conversion. Path is in cm; tof is in ns.
            // TODO: check if this unit conversion is needed!
            // A cut on delta beta:
            bool Bad_dBeta_n_CutCondition = (fabs(beta - (path * 100) / (ToF * c)) > dBeta_n_cut);
            if (pInCD) { histograms.Test_dBeta_n_Step0_epCDn.FillTestHistograms(isGN, isBN, beta - (path * 100) / (ToF * c), weight, !Bad_dBeta_n_CutCondition); }
            if (pInFD) { histograms.Test_dBeta_n_Step0_epFDn.FillTestHistograms(isGN, isBN, beta - (path * 100) / (ToF * c), weight, !Bad_dBeta_n_CutCondition); }

            // A cut on the z-component of the CND hit
            // This is a fiducial cut on the range that the CND can reach on the z-axis
            bool Bad_Vz_n_CutCondition = (v_hit_3v.Z() < Vz_n_lcut || v_hit_3v.Z() > Vz_n_ucut);
            if (pInCD) { histograms.Test_Vz_n_Step0_epCDn.FillTestHistograms(isGN, isBN, v_hit_3v.Z(), weight, !Bad_Vz_n_CutCondition); }
            if (pInFD) { histograms.Test_Vz_n_Step0_epFDn.FillTestHistograms(isGN, isBN, v_hit_3v.Z(), weight, !Bad_Vz_n_CutCondition); }

            bool Bad_ToF_n_CutCondition = (ToF < ToF_n_lcut || ToF > ToF_n_ucut);
            if (pInCD) { histograms.Test_ToF_n_Step0_epCDn.FillTestHistograms(isGN, isBN, ToF, weight, !Bad_ToF_n_CutCondition); }
            if (pInFD) { histograms.Test_ToF_n_Step0_epFDn.FillTestHistograms(isGN, isBN, ToF, weight, !Bad_ToF_n_CutCondition); }

            if (Bad_dBeta_n_CutCondition) { continue; }

            if (Bad_Vz_n_CutCondition) { continue; }

            if (Bad_ToF_n_CutCondition) { continue; }

            pass_step0_cuts = true;

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step0, counter_n_multiplicity_goodN_epCDn_Step0, counter_n_multiplicity_badN_epCDn_Step0,
                               counter_n_multiplicity_allN_epFDn_Step0, counter_n_multiplicity_goodN_epFDn_Step0, counter_n_multiplicity_badN_epFDn_Step0);

            /* Fill AS0C plots */
            histograms.UpdateAS0CHistograms(pInCD, pInFD, P_n_3v, v_hit_3v, beta, path, ToF, weight);

            /* Fill other Step0 plots */
            histograms.UpdateStep0Histograms(pInCD, pInFD, isGN, isBN, P_q_3v, Q2, P_p_3v, P_miss_3v, P_n_3v, E_p, E_miss, M_miss, xB, dpp, theta_n_miss, Edep_CND, Edep_CND1, Edep_CND2, Edep_CND3,
                                             Edep_CTOF, nSector, Size_CND1, Size_CND2, Size_CND3, LayerMult_CND1, LayerMult_CND2, LayerMult_CND3, beta, beta_p, path, ToF, weight);

#pragma endregion /* Step Zero - end */

            //////////////////////////////////////////////
            // Step One
            //////////////////////////////////////////////

#pragma region /* Step One - start */

            // Step One = Dep. energy cut

            // Total deposited energy in CND cut:
            // Upper: Edep_CND > (gamma - 1) * mN * 1000 -> the neutron's deposited energy should not exceed its relativistic kinematic energy. Factor 1000 -> convert GeV to MeV!
            // Lower: Edep_CND < 5 ->
            // TODO: add lower Edep_CND cut?
            bool Bad_Edep_CND_CutCondition = ((Edep_CND < Edep_CND_lcut) || (Edep_CND > (gamma - 1) * mN * 1000));
            if (pInCD) { histograms.Test_Edep_CND_Step1_epCDn.FillTestHistograms(isGN, isBN, Edep_CND, weight, !Bad_Edep_CND_CutCondition); }
            if (pInFD) { histograms.Test_Edep_CND_Step1_epFDn.FillTestHistograms(isGN, isBN, Edep_CND, weight, !Bad_Edep_CND_CutCondition); }

            if (Bad_Edep_CND_CutCondition) { continue; }

            pass_step1_cuts = true;

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step1, counter_n_multiplicity_goodN_epCDn_Step1, counter_n_multiplicity_badN_epCDn_Step1,
                               counter_n_multiplicity_allN_epFDn_Step1, counter_n_multiplicity_goodN_epFDn_Step1, counter_n_multiplicity_badN_epFDn_Step1);

            /* Fill other Step1 plots */
            histograms.UpdateStep1Histograms(pInCD, pInFD, isGN, isBN, P_q_3v, Q2, P_p_3v, P_miss_3v, P_n_3v, E_p, E_miss, M_miss, xB, dpp, theta_n_miss, Edep_CND, Edep_CND1, Edep_CND2, Edep_CND3,
                                             Edep_CTOF, nSector, Size_CND1, Size_CND2, Size_CND3, LayerMult_CND1, LayerMult_CND2, LayerMult_CND3, beta, beta_p, path, ToF, weight);

#pragma endregion /* Step One - end */

            //////////////////////////////////////////////
            // Step Two
            //////////////////////////////////////////////

            // TODO: try to veto Nearby_clusters_from_posPart_tracks by looking at ldiff and sdiff vs TOF difference between the neutron and the cPart

#pragma region /* Step Two - start */

            // Step two = cut/veto out neutrons with charged particles close by

#pragma region /* Step 2 preparations - start */
            bool Nearby_clusters_from_posPart_tracks = false;
            bool Nearby_clusters_from_neutPart_tracks = false;

            bool Proper_layer_multi = false;

            /* Filling ToF * c - v_hit_3v.Z() before cut */
            histograms.UpdateStep2prepBCHistograms(pInCD, pInFD, isGN, isBN, v_hit_3v, ToF, weight);

            for (int itr2_pos = 0; itr2_pos < AllParticles.size(); itr2_pos++) {
                // Why skip itr2_pos == 0? it is the electron:
                if (itr2_pos == 0) { continue; }

                if (itr2_pos == itr1) { continue; }

                // Cut negatively charged particles
                // TODO: Maybe it is good to keep the negatively charged particles in the future.
                if (AllParticles[itr2_pos]->par()->getCharge() <= 0) { continue; }

                // Why this cut? because the background (protons) have high probability of hitting the CTOF? all charged particles supposed to have a CTOF hit at the time of
                // writing the code Cut out particles WITHOUT a CTOF hit:
                if (AllParticles[itr2_pos]->sci(CTOF)->getDetector() == 0) { continue; }

                // TODO: what is this? check for sectors with proton hits in any of the layers of the CND and CTOF?
                int vetoSectorbyLayer[4] = {(AllParticles[itr2_pos]->sci(CTOF)->getComponent() + 1) / 2,  // Normalizes CTOF components to CND sectors (since vetoSectorbyLayer is an array if integers)
                                            AllParticles[itr2_pos]->sci(CND1)->getSector(), AllParticles[itr2_pos]->sci(CND2)->getSector(), AllParticles[itr2_pos]->sci(CND3)->getSector()};

                TVector3 p_C_3v;  // Momentum of the charged particle in the itr2_pos-th entry of AllParticles
                p_C_3v.SetMagThetaPhi(AllParticles[itr2_pos]->getP(), AllParticles[itr2_pos]->getTheta(), AllParticles[itr2_pos]->getPhi());

                double Edep_CTOF_pos = AllParticles[itr2_pos]->sci(clas12::CTOF)->getEnergy();
                // E_dep of positively charged particle

                for (int itr3_pos = 0; itr3_pos < 4; itr3_pos++)  //
                {
                    // TODO: why this cut? no hit in the itr3_pos-th layer?
                    if (vetoSectorbyLayer[itr3_pos] == 0) { continue; }

                    int sdiff = nSector - vetoSectorbyLayer[itr3_pos];

                    // sdiff normalization
                    if (sdiff <= -12) {
                        sdiff += 24;
                    } else if (sdiff > 12) {
                        sdiff -= 24;
                    }

                    int ldiff = detINTlayer - itr3_pos;

                    double ToF_n = ToF;                                                                                   // Neutron ToF
                    double ToF_pos = AllParticles[itr2_pos]->getPath() / (AllParticles[itr2_pos]->par()->getBeta() * c);  // Measured pos particle ToF

                    double dToF = ToF_n - ToF_pos;
                    double dToF_rel_pos = dToF / ToF_pos;
                    double dToF_rel_n = dToF / ToF_n;

                    histograms.UpdateStep2prepPosHistograms(pInCD, pInFD, isGN, isBN, ldiff, sdiff, p_C_3v, v_hit_3v, P_n_3v, dToF, dToF_rel_pos, dToF_rel_n, dpp, theta_n_miss, Edep_CND, beta, path,
                                                            ToF, weight);

                    bool Bad_sdiff_of1_CutCondition = (abs(sdiff) <= 1);
                    bool Bad_sdiff_of2_CutCondition = (abs(sdiff) <= 2);
                    if (pInCD) { histograms.Test_sdiff_of1_pos_Step2_layer_epCDn[ldiff + 3].FillTestHistograms(isGN, isBN, sdiff, weight, !Bad_sdiff_of1_CutCondition); }
                    if (pInFD) { histograms.Test_sdiff_of1_pos_Step2_layer_epFDn[ldiff + 3].FillTestHistograms(isGN, isBN, sdiff, weight, !Bad_sdiff_of1_CutCondition); }
                    if (pInCD) { histograms.Test_sdiff_of2_pos_Step2_layer_epCDn[ldiff + 3].FillTestHistograms(isGN, isBN, sdiff, weight, !Bad_sdiff_of2_CutCondition); }
                    if (pInFD) { histograms.Test_sdiff_of2_pos_Step2_layer_epFDn[ldiff + 3].FillTestHistograms(isGN, isBN, sdiff, weight, !Bad_sdiff_of2_CutCondition); }

                    if (                               // Set the cut on neutrons with nearby clusters from charged particle tracks:
                                                       // Bad_sdiff_of1_CutCondition ||                                  // Minimal sdiff is 2
                        Bad_sdiff_of1_CutCondition 
                        // ||  // Minimal sdiff is 3
                        // // isPosNear_PhiCut(sdiff, ldiff, P_n_3v.Phi() * 180. / M_PI) ||  // Phi_n cut
                        isPosNear_dToF(sdiff, ldiff, dToF)  // ToF difference cut
                    ) {
                        Nearby_clusters_from_posPart_tracks = true;
                    }
                }  // End of loop over vetoSectorbyLayer

                histograms.UpdateMonitorStep2prepHistograms1(Nearby_clusters_from_posPart_tracks, pInCD, pInFD, isGN, isBN, Edep_CND, Edep_CTOF_pos, weight);
            }  // End of second loop over AllParticles (step 1)

            histograms.UpdateMonitorStep2prepPosHistograms2(Nearby_clusters_from_posPart_tracks, pInCD, pInFD, isGN, isBN, Edep_CND, ToF, v_hit_3v, weight);

            for (int itr2_neut = itr1 + 1; itr2_neut < AllParticles.size(); itr2_neut++) {
                // Cut charged particles:
                if (AllParticles[itr2_neut]->par()->getCharge() != 0) { continue; }

                bool CT_neut = (AllParticles[itr2_neut]->sci(clas12::CTOF)->getDetector() == 4);
                bool C1_neut = (AllParticles[itr2_neut]->sci(clas12::CND1)->getDetector() == 3);
                bool C2_neut = (AllParticles[itr2_neut]->sci(clas12::CND2)->getDetector() == 3);
                bool C3_neut = (AllParticles[itr2_neut]->sci(clas12::CND3)->getDetector() == 3);

                // Cut out neutrons without a CND hit in one of it's layers:
                if (!(C1_neut || C2_neut || C3_neut)) { continue; }

                // Use CTOF as a veto for charged particles:
                if (CT_neut) { continue; }

                double theta_neut = AllParticles[itr2_neut]->getTheta() * 180 / M_PI;
                double beta_neut = AllParticles[itr2_neut]->par()->getBeta();
                double gamma_neut = 1 / sqrt(1 - (beta_neut * beta_neut));
                double mom_neut = gamma_neut * beta_neut * mN;
                double ToF_neut = AllParticles[itr2_neut]->getTime() - starttime;

                int detINTlayer_neut = C1_neut ? 1 : C2_neut ? 2 : 3;
                auto detlayer_neut = C1_neut ? CND1 : C2_neut ? CND2 : CND3;  // CND layer with hit

                double nvtx_x_neut = AllParticles[itr2_neut]->par()->getVx();
                double nvtx_y_neut = AllParticles[itr2_neut]->par()->getVy();
                double nvtx_z_neut = AllParticles[itr2_neut]->par()->getVz();
                TVector3 v_nvtx_3v_neut(nvtx_x_neut, nvtx_y_neut, nvtx_z_neut);  // Neutron's vertex location

                TVector3 v_hit_3v_neut;  // Neutron's hit location in CND
                v_hit_3v_neut.SetXYZ(AllParticles[itr2_neut]->sci(detlayer_neut)->getX(), AllParticles[itr2_neut]->sci(detlayer_neut)->getY(), AllParticles[itr2_neut]->sci(detlayer_neut)->getZ());

                TVector3 v_path_3v_neut = v_hit_3v_neut - v_nvtx_3v_neut;  // Direct calculation of neutron's path (in vector form)

                TVector3 P_neut_3v;  // Momentum of the charged particle in the itr2_neut-th entry of AllParticles
                P_neut_3v.SetMagThetaPhi(mom_neut, v_path_3v_neut.Theta(), v_path_3v_neut.Phi());

                // Why "v_path_3v.Mag() / 100"? unit conversion.
                // TODO: check if this unit conversion is needed!
                double path_neut = v_path_3v_neut.Mag() / 100;
                // double path = v_path_3v.Mag();
                double theta_n_miss_neut = P_neut_3v.Angle(P_miss_3v) * 180 / M_PI;
                // Opening angle between calculated neutron's momentum and predicted neutron momentum (= missing momentum)
                double dpp_neut = (P_miss_3v.Mag() - P_neut_3v.Mag()) / P_miss_3v.Mag();

                // Beta cut:
                // Upper: beta > 0.8 -> cut out photons
                // Lower: beta < 0.15 ->
                if (beta_neut < Beta_n_lcut || beta_neut > Beta_n_ucut) { continue; }

                // Why this cut? reco code bug. Neutrons in this angle range are in the BAND and appear in the CND.
                // This bug is probobly fixed, yet the cut is still applied to mak sure.
                // Updated: now is forcing the neutron to be inside the acceptance of the CD
                if ((P_neut_3v.Theta() * 180. / M_PI < Theta_n_lcut) || (P_neut_3v.Theta() * 180. / M_PI > Theta_n_ucut)) { continue; }

                // Status cut for double-hits (based on Erin's code)
                int status_neut = 0;
                if (C1_neut) { status_neut = status_neut + AllParticles[itr2_neut]->sci(CND1)->getStatus(); }
                if (C2_neut) { status_neut = status_neut + AllParticles[itr2_neut]->sci(CND3)->getStatus(); }
                if (C3_neut) { status_neut = status_neut + AllParticles[itr2_neut]->sci(CND2)->getStatus(); }
                if (status_neut != 0) { continue; }

                // Check to see if there is a good neutron
                bool isGN_neut = false;
                bool isBN_neut = false;

                // Good neutron definition:
                bool GN_theta_n_miss_neut = (theta_n_miss_neut <= GN_theta_n_miss_ucut);
                bool GN_dpp_neut = ((dpp_neut >= GN_dpp_lcut) && (dpp_neut <= GN_dpp_ucut));
                if (GN_theta_n_miss_neut && GN_dpp_neut) { isGN_neut = true; }

                // Bad neutron definition:
                bool BN_theta_n_miss_neut = (theta_n_miss_neut >= BN_theta_n_miss_lcut);
                bool BN_dpp_neut = (dpp_neut <= BN_dpp_ucut);
                if (BN_theta_n_miss_neut || BN_dpp_neut) { isBN_neut = true; }
                // if (BN_theta_n_miss && BN_dpp) { isBN = true; }
                // if (!((theta_n_miss < 25.) && ((dpp > -0.3) && (dpp < 0.3)))) { isBN = true; }

                if (isGN_neut && isBN_neut) { cout << "\nERROR! good and bad neutrons are overlapping (Step2)! Aborting...\n", exit(0); }

                if (!(isGN_neut || isBN_neut)) { continue; }

                // Why "path * 100"? unit conversion. Path is in cm; tof is in ns.
                // TODO: check if this unit conversion is needed!
                // A cut on delta beta:
                bool Bad_dBeta_n_CutCondition_neut = (fabs(beta_neut - (path_neut * 100) / (ToF_neut * c)) > dBeta_n_cut);
                if (Bad_dBeta_n_CutCondition_neut) { continue; }

                // A cut on the z-component of the CND hit
                // This is a fiducial cut on the range that the CND can reach on the z-axis
                bool Bad_Vz_n_CutCondition_neut = (v_hit_3v_neut.Z() < Vz_n_lcut || v_hit_3v_neut.Z() > Vz_n_ucut);
                if (Bad_Vz_n_CutCondition_neut) { continue; }

                bool Bad_ToF_n_CutCondition_neut = (ToF_neut < ToF_n_lcut || ToF_neut > ToF_n_ucut);
                if (Bad_ToF_n_CutCondition_neut) { continue; }

                // Total deposited energy in CND cut:
                // Upper: Edep_CND > (gamma - 1) * mN * 1000 -> the neutron's deposited energy should not exceed its relativistic kinematic energy. Factor 1000 -> convert GeV to MeV!
                // Lower: Edep_CND < 5 ->
                // TODO: add lower Edep_CND cut?
                double Edep_CND1_neut = AllParticles[itr2_neut]->sci(CND1)->getEnergy();
                double Edep_CND2_neut = AllParticles[itr2_neut]->sci(CND2)->getEnergy();
                double Edep_CND3_neut = AllParticles[itr2_neut]->sci(CND3)->getEnergy();
                double Edep_CND_neut = Edep_CND1_neut + Edep_CND2_neut + Edep_CND3_neut;
                bool Bad_Edep_CND_CutCondition_neut = ((Edep_CND_neut < Edep_CND_lcut) || (Edep_CND_neut > (gamma_neut - 1) * mN * 1000));
                if (Bad_Edep_CND_CutCondition_neut) { continue; }

                // TODO: what is this? check for sectors with proton hits in any of the layers of the CND and CTOF?
                int vetoSectorbyLayer[3] = {AllParticles[itr2_neut]->sci(CND1)->getSector(), AllParticles[itr2_neut]->sci(CND2)->getSector(), AllParticles[itr2_neut]->sci(CND3)->getSector()};
                int vetoSectorbyLayerComponent_neut[3] = {AllParticles[itr2_neut]->sci(CND1)->getComponent(), AllParticles[itr2_neut]->sci(CND2)->getComponent(),
                                                          AllParticles[itr2_neut]->sci(CND3)->getComponent()};

                int vetoSectorbyLayerComponent[3] = {AllParticles[itr1]->sci(CND1)->getComponent(), AllParticles[itr1]->sci(CND2)->getComponent(), AllParticles[itr1]->sci(CND3)->getComponent()};

                for (int itr3_neut = 0; itr3_neut < 3; itr3_neut++)  //
                {
                    int sdiff = nSector - vetoSectorbyLayer[itr3_neut];

                    // sdiff normalization
                    if (sdiff <= -12) {
                        sdiff += 24;
                    } else if (sdiff > 12) {
                        sdiff -= 24;
                    }

                    int ldiff = detINTlayer - (itr3_neut + 1);

                    double ToF_n = ToF;  // Neutron ToF
                    // double ToF_neut = AllParticles[itr2_neut]->getPath() / (AllParticles[itr2_neut]->par()->getBeta() * c);
                    // Measured neut particle ToF

                    double dToF = ToF_n - ToF_neut;
                    double dToF_rel_neut = dToF / ToF_neut;
                    double dToF_rel_n = dToF / ToF_n;

                    histograms.UpdateStep2prepNeutHistograms(pInCD, pInFD, isGN, isBN, ldiff, sdiff, P_neut_3v, v_hit_3v, P_n_3v, dToF, dToF_rel_neut, dToF_rel_n, dpp, theta_n_miss, Edep_CND, beta,
                                                             path, ToF, weight);

                    bool SameSector = (sdiff == 0);
                    bool NeutBeforeN = (ldiff > 0);
                    bool SameComponent = (vetoSectorbyLayerComponent_neut[itr3_neut] == vetoSectorbyLayerComponent[itr3_neut]);
                    bool NeutEarlierN = (dToF > 0);

                    if (  // Set the cut on neutrons with nearby clusters from other neutal particles:
                        SameSector && NeutBeforeN && SameComponent && NeutEarlierN
                        // (sdiff == 0) && (ldiff > 0)
                        // (sdiff == 0) && (ldiff > 0) && (dToF <= 0)
                    ) {
                        Nearby_clusters_from_neutPart_tracks = true;
                    }
                }  // End of loop over vetoSectorbyLayer

                // histograms.UpdateMonitorStep2prepHistograms1(Nearby_clusters_from_posPart_tracks, pInCD, pInFD, isGN, isBN, Edep_CND, Edep_CTOF_neut,
                //                                              weight);
            }  // End of second loop over AllParticles (step 1)

#pragma endregion /* Step 2 preparations - end */

            /* Fill BS2C plots */
            histograms.UpdateBS2CHistograms(pInCD, pInFD, Size_CND1, Size_CND2, Size_CND3, LayerMult_CND1, LayerMult_CND2, LayerMult_CND3, weight);

            // Cutting out neutrons with nearby hits from charged particle tracks
            bool Bad_posTrack_prox_CutCondition = Nearby_clusters_from_posPart_tracks;

            // Cutting out neutrons with nearby hits from other neutrals
            bool Bad_neutTrack_prox_CutCondition = Nearby_clusters_from_neutPart_tracks;

            // Cutting out neutrons cluster width greater than 1
            // Neutrons are neutral (i.e., no curved tracks), and so the can only hit one scintillator paddle (i.e., width = 1)
            bool Bad_Size_CND1_CutCondition = (C1 && (Size_CND1 != Cluster_size_cut));
            bool Bad_Size_CND2_CutCondition = (C2 && (Size_CND2 != Cluster_size_cut));
            bool Bad_Size_CND3_CutCondition = (C3 && (Size_CND3 != Cluster_size_cut));
            if (pInCD) { histograms.Test_Size_CND1_Step2_epCDn.FillTestHistograms(isGN, isBN, Size_CND1, weight, !Bad_Size_CND1_CutCondition); }
            if (pInFD) { histograms.Test_Size_CND1_Step2_epFDn.FillTestHistograms(isGN, isBN, Size_CND1, weight, !Bad_Size_CND1_CutCondition); }
            if (pInCD) { histograms.Test_Size_CND2_Step2_epCDn.FillTestHistograms(isGN, isBN, Size_CND2, weight, !Bad_Size_CND2_CutCondition); }
            if (pInFD) { histograms.Test_Size_CND2_Step2_epFDn.FillTestHistograms(isGN, isBN, Size_CND2, weight, !Bad_Size_CND2_CutCondition); }
            if (pInCD) { histograms.Test_Size_CND3_Step2_epCDn.FillTestHistograms(isGN, isBN, Size_CND3, weight, !Bad_Size_CND3_CutCondition); }
            if (pInFD) { histograms.Test_Size_CND3_Step2_epFDn.FillTestHistograms(isGN, isBN, Size_CND3, weight, !Bad_Size_CND3_CutCondition); }

            // // Cutting out neutrons without:
            // // 1. A hit in CND1 with layer multiplicity of one
            // // 2. A hit in CND2 or CND3 with layer multiplicity of three
            int LayerMult_CND2andCND3 = LayerMult_CND1 + LayerMult_CND2;
            bool Bad_LayerMult_CND1_CutCondition = (C1 && LayerMult_CND != CND1_LayerMult_cut);                        // Condition 1
            bool Bad_LayerMult_CND2andCND3_CutCondition = ((C2 || C3) && LayerMult_CND > CND2andCND3_LayerMult_ucut);  // Condition 2
            if (pInCD) { histograms.Test_LayerMult_CND1_Step2_epCDn.FillTestHistograms(isGN, isBN, LayerMult_CND1, weight, !Bad_LayerMult_CND1_CutCondition); }
            if (pInFD) { histograms.Test_LayerMult_CND1_Step2_epFDn.FillTestHistograms(isGN, isBN, LayerMult_CND1, weight, !Bad_LayerMult_CND1_CutCondition); }
            if (pInCD) { histograms.Test_LayerMult_CND2andCND3_Step2_epCDn.FillTestHistograms(isGN, isBN, LayerMult_CND2andCND3, weight, !Bad_LayerMult_CND2andCND3_CutCondition); }
            if (pInFD) { histograms.Test_LayerMult_CND2andCND3_Step2_epFDn.FillTestHistograms(isGN, isBN, LayerMult_CND2andCND3, weight, !Bad_LayerMult_CND2andCND3_CutCondition); }

            if (Bad_posTrack_prox_CutCondition) { continue; }

            // if (Bad_neutTrack_prox_CutCondition) { continue; }

            // if (Bad_Size_CND1_CutCondition || Bad_Size_CND2_CutCondition || Bad_Size_CND3_CutCondition) { continue; }

            // if (Bad_LayerMult_CND1_CutCondition || Bad_LayerMult_CND2andCND3_CutCondition) { continue; }

            pass_step2_cuts = true;

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step2, counter_n_multiplicity_goodN_epCDn_Step2, counter_n_multiplicity_badN_epCDn_Step2,
                               counter_n_multiplicity_allN_epFDn_Step2, counter_n_multiplicity_goodN_epFDn_Step2, counter_n_multiplicity_badN_epFDn_Step2);

            histograms.UpdateAS2CHistograms(pInCD, pInFD, Size_CND1, Size_CND2, Size_CND3, LayerMult_CND1, LayerMult_CND2, LayerMult_CND3, weight);

            /* Fill other Step2 plots */
            histograms.UpdateStep2Histograms(pInCD, pInFD, isGN, isBN, P_q_3v, Q2, P_p_3v, P_miss_3v, P_n_3v, E_p, E_miss, M_miss, xB, dpp, theta_n_miss, Edep_CND, Edep_CND1, Edep_CND2, Edep_CND3,
                                             Edep_CTOF, nSector, Size_CND1, Size_CND2, Size_CND3, LayerMult_CND1, LayerMult_CND2, LayerMult_CND3, beta, beta_p, path, ToF, weight);

            for (int itr4_pos = 0; itr4_pos < AllParticles.size(); itr4_pos++) {
                // Why skip itr4_pos == 0? it is the electron
                if (itr4_pos == 0) { continue; }

                if (itr4_pos == itr1) { continue; }

                // Cut negatively charged particles
                // TODO: Maybe it is good to keep the nagativly charged particles in the future.
                if (AllParticles[itr4_pos]->par()->getCharge() <= 0) { continue; }

                // Why this cut? because the background (protons) have high probability of hitting the CTOF? all charged particles supposed to have a CTOF hit at the time of
                // writing the code Cut out particles WITHOUT a CTOF hit
                if (AllParticles[itr4_pos]->sci(CTOF)->getDetector() == 0) { continue; }

                int vetoSectorbyLayer[4] = {(AllParticles[itr4_pos]->sci(CTOF)->getComponent() + 1) / 2, AllParticles[itr4_pos]->sci(CND1)->getSector(), AllParticles[itr4_pos]->sci(CND2)->getSector(),
                                            AllParticles[itr4_pos]->sci(CND3)->getSector()};

                TVector3 p_C_3v;  // Momentum of the charged particle in the itr4_pos-th entry of AllParticles
                p_C_3v.SetMagThetaPhi(AllParticles[itr4_pos]->getP(), AllParticles[itr4_pos]->getTheta(), AllParticles[itr4_pos]->getPhi());

                for (int itr5_pos = 0; itr5_pos < 4; itr5_pos++) {
                    // TODO: why this cut? no hit in the itr5_pos-th layer?
                    if (vetoSectorbyLayer[itr5_pos] == 0) { continue; }

                    int sdiff = nSector - vetoSectorbyLayer[itr5_pos];

                    // sdiff normalization
                    if (sdiff <= -12) {
                        sdiff += 24;
                    } else if (sdiff > 12) {
                        sdiff -= 24;
                    }

                    int ldiff = detINTlayer - itr5_pos;

                    double ToF_n = ToF;  // Neutron ToF
                    double ToF_pos = AllParticles[itr4_pos]->getPath() / (AllParticles[itr4_pos]->par()->getBeta() * c);
                    // Measured pos particle ToF

                    double dToF = ToF_n - ToF_pos;
                    double dToF_rel_pos = dToF / ToF_pos;
                    double dToF_rel_n = dToF / ToF_n;

                    histograms.UpdateStep2PosHistograms2(pInCD, pInFD, isGN, isBN, ldiff, sdiff, p_C_3v, v_hit_3v, P_n_3v, dToF, dToF_rel_pos, dToF_rel_n, dpp, theta_n_miss, Edep_CND, beta, path, ToF,
                                                         weight);
                }
            }  // End of third loop over AllParticles (step 2)

            for (int itr4_neut = itr1 + 1; itr4_neut < AllParticles.size(); itr4_neut++) {
                // Cut charged particles:
                if (AllParticles[itr4_neut]->par()->getCharge() != 0) { continue; }

                bool CT_neut = (AllParticles[itr4_neut]->sci(clas12::CTOF)->getDetector() == 4);
                bool C1_neut = (AllParticles[itr4_neut]->sci(clas12::CND1)->getDetector() == 3);
                bool C2_neut = (AllParticles[itr4_neut]->sci(clas12::CND2)->getDetector() == 3);
                bool C3_neut = (AllParticles[itr4_neut]->sci(clas12::CND3)->getDetector() == 3);

                // Cut out neutrons without a CND hit in one of it's layers:
                if (!(C1_neut || C2_neut || C3_neut)) { continue; }

                // Use CTOF as a veto for charged particles:
                if (CT_neut) { continue; }

                double theta_neut = AllParticles[itr4_neut]->getTheta() * 180 / M_PI;
                double beta_neut = AllParticles[itr4_neut]->par()->getBeta();
                double gamma_neut = 1 / sqrt(1 - (beta_neut * beta_neut));
                double mom_neut = gamma_neut * beta_neut * mN;
                double ToF_neut = AllParticles[itr4_neut]->getTime() - starttime;

                int detINTlayer_neut = C1_neut ? 1 : C2_neut ? 2 : 3;
                auto detlayer_neut = C1_neut ? CND1 : C2_neut ? CND2 : CND3;  // CND layer with hit

                double nvtx_x_neut = AllParticles[itr4_neut]->par()->getVx();
                double nvtx_y_neut = AllParticles[itr4_neut]->par()->getVy();
                double nvtx_z_neut = AllParticles[itr4_neut]->par()->getVz();
                TVector3 v_nvtx_3v_neut(nvtx_x_neut, nvtx_y_neut, nvtx_z_neut);  // Neutron's vertex location

                TVector3 v_hit_3v_neut;  // Neutron's hit location in CND
                v_hit_3v_neut.SetXYZ(AllParticles[itr4_neut]->sci(detlayer_neut)->getX(), AllParticles[itr4_neut]->sci(detlayer_neut)->getY(), AllParticles[itr4_neut]->sci(detlayer_neut)->getZ());

                TVector3 v_path_3v_neut = v_hit_3v_neut - v_nvtx_3v_neut;  // Direct calculation of neutron's path (in vector form)

                TVector3 P_neut_3v;  // Momentum of the charged particle in the itr4_neut-th entry of AllParticles
                P_neut_3v.SetMagThetaPhi(mom_neut, v_path_3v_neut.Theta(), v_path_3v_neut.Phi());

                // Why "v_path_3v.Mag() / 100"? unit conversion.
                // TODO: check if this unit conversion is needed!
                double path_neut = v_path_3v_neut.Mag() / 100;
                // double path = v_path_3v.Mag();
                double theta_n_miss_neut = P_neut_3v.Angle(P_miss_3v) * 180 / M_PI;
                // Opening angle between calculated neutron's momentum and predicted neutron momentum (= missing momentum)
                double dpp_neut = (P_miss_3v.Mag() - P_neut_3v.Mag()) / P_miss_3v.Mag();

                // Beta cut:
                // Upper: beta > 0.8 -> cut out photons
                // Lower: beta < 0.15 ->
                if (beta_neut < Beta_n_lcut || beta_neut > Beta_n_ucut) { continue; }

                // Why this cut? reco code bug. Neutrons in this angle range are in the BAND and appear in the CND.
                // This bug is probobly fixed, yet the cut is still applied to mak sure.
                // Updated: now is forcing the neutron to be inside the acceptance of the CD
                if ((P_neut_3v.Theta() * 180. / M_PI < Theta_n_lcut) || (P_neut_3v.Theta() * 180. / M_PI > Theta_n_ucut)) { continue; }

                // Status cut for double-hits (based on Erin's code)
                int status_neut = 0;
                if (C1_neut) { status_neut = status_neut + AllParticles[itr4_neut]->sci(CND1)->getStatus(); }
                if (C2_neut) { status_neut = status_neut + AllParticles[itr4_neut]->sci(CND3)->getStatus(); }
                if (C3_neut) { status_neut = status_neut + AllParticles[itr4_neut]->sci(CND2)->getStatus(); }
                if (status_neut != 0) { continue; }

                // Check to see if there is a good neutron
                bool isGN_neut = false;
                bool isBN_neut = false;

                // Good neutron definition:
                bool GN_theta_n_miss_neut = (theta_n_miss_neut <= GN_theta_n_miss_ucut);
                bool GN_dpp_neut = ((dpp_neut >= GN_dpp_lcut) && (dpp_neut <= GN_dpp_ucut));
                if (GN_theta_n_miss_neut && GN_dpp_neut) { isGN_neut = true; }

                // Bad neutron definition:
                bool BN_theta_n_miss_neut = (theta_n_miss_neut >= BN_theta_n_miss_lcut);
                bool BN_dpp_neut = (dpp_neut <= BN_dpp_ucut);
                if (BN_theta_n_miss_neut || BN_dpp_neut) { isBN_neut = true; }
                // if (BN_theta_n_miss && BN_dpp) { isBN = true; }
                // if (!((theta_n_miss < 25.) && ((dpp > -0.3) && (dpp < 0.3)))) { isBN = true; }

                if (isGN_neut && isBN_neut) { cout << "\nERROR! good and bad neutrons are overlapping (Step2)! Aborting...\n", exit(0); }

                if (!(isGN_neut || isBN_neut)) { continue; }

                // Why "path * 100"? unit conversion. Path is in cm; tof is in ns.
                // TODO: check if this unit conversion is needed!
                // A cut on delta beta:
                bool Bad_dBeta_n_CutCondition_neut = (fabs(beta_neut - (path_neut * 100) / (ToF_neut * c)) > dBeta_n_cut);
                if (Bad_dBeta_n_CutCondition_neut) { continue; }

                // A cut on the z-component of the CND hit
                // This is a fiducial cut on the range that the CND can reach on the z-axis
                bool Bad_Vz_n_CutCondition_neut = (v_hit_3v_neut.Z() < Vz_n_lcut || v_hit_3v_neut.Z() > Vz_n_ucut);
                if (Bad_Vz_n_CutCondition_neut) { continue; }

                bool Bad_ToF_n_CutCondition_neut = (ToF_neut < ToF_n_lcut || ToF_neut > ToF_n_ucut);
                if (Bad_ToF_n_CutCondition_neut) { continue; }

                // Total deposited energy in CND cut:
                // Upper: Edep_CND > (gamma - 1) * mN * 1000 -> the neutron's deposited energy should not exceed its relativistic kinematic energy. Factor 1000 -> convert GeV to MeV!
                // Lower: Edep_CND < 5 ->
                // TODO: add lower Edep_CND cut?
                double Edep_CND1_neut = AllParticles[itr4_neut]->sci(CND1)->getEnergy();
                double Edep_CND2_neut = AllParticles[itr4_neut]->sci(CND2)->getEnergy();
                double Edep_CND3_neut = AllParticles[itr4_neut]->sci(CND3)->getEnergy();
                double Edep_CND_neut = Edep_CND1_neut + Edep_CND2_neut + Edep_CND3_neut;
                bool Bad_Edep_CND_CutCondition_neut = ((Edep_CND_neut < Edep_CND_lcut) || (Edep_CND_neut > (gamma_neut - 1) * mN * 1000));
                if (Bad_Edep_CND_CutCondition_neut) { continue; }

                // TODO: what is this? check for sectors with proton hits in any of the layers of the CND and CTOF?
                int vetoSectorbyLayer[3] = {AllParticles[itr4_neut]->sci(CND1)->getSector(), AllParticles[itr4_neut]->sci(CND2)->getSector(), AllParticles[itr4_neut]->sci(CND3)->getSector()};
                int vetoSectorbyLayerComponent_neut[3] = {AllParticles[itr4_neut]->sci(CND1)->getComponent(), AllParticles[itr4_neut]->sci(CND2)->getComponent(),
                                                          AllParticles[itr4_neut]->sci(CND3)->getComponent()};

                int vetoSectorbyLayerComponent[3] = {AllParticles[itr1]->sci(CND1)->getComponent(), AllParticles[itr1]->sci(CND2)->getComponent(), AllParticles[itr1]->sci(CND3)->getComponent()};

                for (int itr5_neut = 0; itr5_neut < 3; itr5_neut++)  //
                {
                    int sdiff = nSector - vetoSectorbyLayer[itr5_neut];

                    // sdiff normalization
                    if (sdiff <= -12) {
                        sdiff += 24;
                    } else if (sdiff > 12) {
                        sdiff -= 24;
                    }

                    int ldiff = detINTlayer - (itr5_neut + 1);

                    double ToF_n = ToF;  // Neutron ToF
                    // double ToF_neut = AllParticles[itr4_neut]->getPath() / (AllParticles[itr4_neut]->par()->getBeta() * c);
                    // Measured neut particle ToF

                    double dToF = ToF_n - ToF_neut;
                    double dToF_rel_neut = dToF / ToF_neut;
                    double dToF_rel_n = dToF / ToF_n;

                    histograms.UpdateStep2NeutHistograms2(pInCD, pInFD, isGN, isBN, ldiff, sdiff, P_neut_3v, v_hit_3v, P_n_3v, dToF, dToF_rel_neut, dToF_rel_n, dpp, theta_n_miss, Edep_CND, beta, path,
                                                          ToF, weight);

                }  // End of loop over vetoSectorbyLayer

                // histograms.UpdateMonitorStep2prepHistograms1(Nearby_clusters_from_posPart_tracks, pInCD, pInFD, isGN, isBN, Edep_CND, Edep_CTOF_neut,
                //                                              weight);
            }  // End of second loop over AllParticles (step 1)

            /*
            bool AllHitVeto = false;
            int hitsNear = 0;

            for (int row = 0; row < c12->getBank(cnd_hits)->getRows(); row++)
            {
                int hit_sector = c12->getBank(cnd_hits)->getInt(cnd_hit_sector, row);
                int hit_layer = c12->getBank(cnd_hits)->getInt(cnd_hit_layer, row);
                double hit_energy = c12->getBank(cnd_hits)->getFloat(cnd_hit_energy, row);

                int sdiff = nSector - hit_sector;
                // sdiff normalization
                if (sdiff <= -12)
                {
                    sdiff += 24;
                }
                else if (sdiff > 12)
                {
                    sdiff -= 24;
                }
                int ldiff = detINTlayer - hit_layer;

                if ((ldiff == 0) && (sdiff == 0))
                {
                    continue;
                }
                if (isNear(sdiff, ldiff))
                {
                    if (hit_energy > 20)
                    {
                        if (isGN)
                        {
                            h_NearbyEdep_goodN_Step2->Fill(hit_energy, weight);
                        }
                        else if (isBN)
                        {
                            h_NearbyEdep_badN_Step2->Fill(hit_energy, weight);
                        }
                        hitsNear++;
                    }
                }

                // if(ToF*c-v_hit_3v.Z() < 70){
                if (ToF < 6)
                {
                    if (isGN)
                    {
                        h_sdiff_allhit_goodN_Step2_layer[ldiff + 3]->Fill(sdiff, weight);
                        h_sdiff_ldiff_allhit_goodN_Step2->Fill(sdiff, ldiff, weight);
                    }
                    else if (isBN)
                    {
                        h_sdiff_allhit_badN_Step2_layer[ldiff + 3]->Fill(sdiff, weight);
                        h_sdiff_ldiff_allhit_badN_Step2->Fill(sdiff, ldiff, weight);
                    }
                }
                //}
            }

            // Now that we have nearby hits, look at how many we have that are not off time
            // if(ToF*c-v_hit_3v.Z() < 70){
            if (ToF < 6)
            {
                if (pInCD)
                {
                    if (isGN)
                    {
                        h_numberNearby_goodN_Step2_epCDn->Fill(hitsNear, weight);
                        h_numberNearby_momN_goodN_Step2_epCDn->Fill(hitsNear, mom, weight);
                        h_nsector_goodN_Step2_epCDn->Fill(nSector, weight);
                    }
                    else if (isBN)
                    {
                        h_numberNearby_badN_Step2_epCDn->Fill(hitsNear, weight);
                        h_numberNearby_momN_badN_Step2_epCDn->Fill(hitsNear, mom, weight);
                        h_nsector_badN_Step2_epCDn->Fill(nSector, weight);
                    }
                }
                else if (pInFD)
                {
                    if (isGN)
                    {
                        h_numberNearby_goodN_Step2_epFDn->Fill(hitsNear, weight);
                        h_numberNearby_momN_goodN_Step2_epFDn->Fill(hitsNear, mom, weight);
                        h_nsector_goodN_Step2_epFDn->Fill(nSector, weight);
                    }
                    else if (isBN)
                    {
                        h_numberNearby_badN_Step2_epFDn->Fill(hitsNear, weight);
                        h_numberNearby_momN_badN_Step2_epFDn->Fill(hitsNear, mom, weight);
                        h_nsector_badN_Step2_epFDn->Fill(nSector, weight);
                    }
                }
            }
            //}
            if (hitsNear >= 1)
            {
                AllHitVeto = true;
            }
            */

#pragma endregion /* Step Two - end */

            //////////////////////////////////////////////
            // Step Three
            //////////////////////////////////////////////

#pragma region /* Step Three - start */

            /*
            if (AllHitVeto)
            {
                continue;
            }

            bool CTOFHitVeto = false;

            int hitsCTOF = 0;

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step3, counter_n_multiplicity_goodN_epCDn_Step3, counter_n_multiplicity_badN_epCDn_Step3,
                               counter_n_multiplicity_allN_epFDn_Step3, counter_n_multiplicity_goodN_epFDn_Step3, counter_n_multiplicity_badN_epFDn_Step3);
            SetNeutronCounters(isGN, counter_n_multiplicity_allN_Step3, counter_n_multiplicity_goodN_Step3, counter_n_multiplicity_badN_Step3);

            h_pnRes_theta_n_miss_Step3->Fill(dpp, theta_n_miss, weight);

            if (isGN)
            {
                h_ToF_goodN_Step3->Fill(ToF, weight);
                h_Edep_ToF_goodN_Step3->Fill(ToF, edep, weight);
            }
            else if (isBN)
            {
                h_ToF_badN_Step3->Fill(ToF, weight);
                h_Edep_ToF_badN_Step3->Fill(ToF, edep, weight);
            }

            for (int row = 0; row < c12->getBank(cnd_hits)->getRows(); row++)
            {
                int hit_sector = c12->getBank(cnd_hits)->getInt(cnd_hit_sector, row);
                int hit_layer = c12->getBank(cnd_hits)->getInt(cnd_hit_layer, row);
                int sdiff = nSector - hit_sector;
                // sdiff normalization
                if (sdiff <= -12)
                {
                    sdiff += 24;
                }
                else if (sdiff > 12)
                {
                    sdiff -= 24;
                }
                int ldiff = detINTlayer - hit_layer;
                if ((ldiff == 0) && (sdiff == 0))
                {
                    continue;
                }
                if (isGN)
                {
                    h_sdiff_ldiff_allhit_goodN_Step3->Fill(sdiff, ldiff, weight);
                }
                else if (isBN)
                {
                    h_sdiff_ldiff_allhit_badN_Step3->Fill(sdiff, ldiff, weight);
                }
            }

            for (int row = 0; row < c12->getBank(ctof_hits)->getRows(); row++)
            {
                int hit_sector = (c12->getBank(ctof_hits)->getInt(ctof_hit_component, row) + 1) / 2;
                double hit_energy = c12->getBank(ctof_hits)->getFloat(ctof_hit_energy, row);

                int sdiff = nSector - hit_sector;
                // sdiff normalization
                if (sdiff <= -12)
                {
                    sdiff += 24;
                }
                else if (sdiff > 12)
                {
                    sdiff -= 24;
                }
                int ldiff = detINTlayer;

                if ((ldiff == 0) && (sdiff == 0))
                {
                    continue;
                }
                if (isNearCTOF(sdiff, ldiff))
                {
                    // if(isGN){h_NearbyEdep_goodN_Step3->Fill(hit_energy,weight);}
                    // else{h_NearbyEdep_badN_Step3->Fill(hit_energy,weight);}
                    if (hit_energy > 5)
                    {
                        hitsCTOF++;
                    }
                }
                if (isGN)
                {
                    h_sdiff_ldiff_CTOFhit_goodN_Step3->Fill(sdiff, ldiff, weight);
                }
                else if (isBN)
                {
                    h_sdiff_ldiff_CTOFhit_badN_Step3->Fill(sdiff, ldiff, weight);
                }
            }

            if (isGN)
            {
                h_numberCTOF_goodN_Step3->Fill(hitsCTOF, weight);
                h_numberCTOF_momN_goodN_Step3->Fill(hitsCTOF, mom, weight);
            }
            else if (isBN)
            {
                h_numberCTOF_badN_Step3->Fill(hitsCTOF, weight);
                h_numberCTOF_momN_badN_Step3->Fill(hitsCTOF, mom, weight);
            }
            if (hitsCTOF >= 1)
            {
                CTOFHitVeto = true;
            }
            */

#pragma endregion /* Step Three - end */

            //////////////////////////////////////////////
            // Step Four
            //////////////////////////////////////////////

#pragma region /* Step Four - start */

            /*
            if (CTOFHitVeto)
            {
                continue;
            }

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step4, counter_n_multiplicity_goodN_epCDn_Step4, counter_n_multiplicity_badN_epCDn_Step4,
                               counter_n_multiplicity_allN_epFDn_Step4, counter_n_multiplicity_goodN_epFDn_Step4, counter_n_multiplicity_badN_epFDn_Step4);
            SetNeutronCounters(isGN, counter_n_multiplicity_allN_Step4, counter_n_multiplicity_goodN_Step4, counter_n_multiplicity_badN_Step4);

            h_pnRes_theta_n_miss_Step4->Fill(dpp, theta_n_miss, weight);

            if (isGN)
            {
                h_ToF_goodN_Step4->Fill(ToF, weight);
                h_edep_ToF_goodN_Step4->Fill(ToF, edep, weight);
            }
            else if (isBN)
            {
                h_ToF_badN_Step4->Fill(ToF, weight);
                h_Edep_ToF_badN_Step4->Fill(ToF, edep, weight);
            }
            */

#pragma endregion /* Step Four - end */

            //////////////////////////////////////////////
            // Step Five
            //////////////////////////////////////////////

#pragma region /* Step Five - start */

            /*
            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step5, counter_n_multiplicity_goodN_epCDn_Step5, counter_n_multiplicity_badN_epCDn_Step5,
                counter_n_multiplicity_allN_epFDn_Step5, counter_n_multiplicity_goodN_epFDn_Step5, counter_n_multiplicity_badN_epFDn_Step5);
            SetNeutronCounters(isGN, counter_n_multiplicity_allN_Step5, counter_n_multiplicity_goodN_Step5, counter_n_multiplicity_badN_Step5);

            _pnRes_theta_n_miss_Step5->Fill(dpp, theta_n_miss, weight);
            h_pmiss_allN_Step5->Fill(P_miss_3v.Mag(), weight);
            if (isGN)
            {
                h_ToF_goodN_Step5->Fill(ToF, weight);
                h_edep_ToF_goodN_Step5->Fill(ToF, edep, weight);
                h_pmiss_goodN_Step5->Fill(P_miss_3v.Mag(), weight);
                h_diff_ToFc_z_Edep_goodN_Step5->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                h_diff_ToFc_z_Edep_goodN_Step5_layer_epCDn[detINTlayer - 1]->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                h_phidiff_en_goodN_Step5->Fill(get_phi_diff(P_e_3v, P_n_3v), weight);
                h_TP_goodN_Step5->Fill(ToF / path * 100, weight);
                h_Z_goodN_Step5->Fill(v_hit_3v.Z(), weight);
                h_beta_Edep_goodN_Step5->Fill(beta, edep, weight);

                h_ToF_Edep_goodN_Step5->Fill(ToF, edep, weight);
                h_TP_Edep_goodN_Step5->Fill(ToF / path * 100, edep, weight);
            }
            else if (isBN)
            {
                h_ToF_badN_Step5->Fill(ToF, weight);
                h_edep_ToF_badN_Step5->Fill(ToF, edep, weight);
                // if(ToF<8){
                h_diff_ToFc_z_Edep_badN_Step5->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                h_diff_ToFc_z_Edep_badN_Step5_layer_epCDn[detINTlayer - 1]->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                h_phidiff_en_badN_Step5->Fill(get_phi_diff(P_e_3v, P_n_3v), weight);
                h_TP_badN_Step5->Fill(ToF / path * 100, weight);
                h_Z_badN_Step5->Fill(v_hit_3v.Z(), weight);
                h_beta_Edep_badN_Step5->Fill(beta, edep, weight);

                h_ToF_Edep_badN_Step5->Fill(ToF, edep, weight);
                h_TP_Edep_badN_Step5->Fill(ToF / path * 100, edep, weight);
                //}
                if (ToF < 5)
                {
                    if (edep > 20)
                    {
                        // cerr<<"Event="<<c12->runconfig()->getEvent()<<endl;
                        // cerr<<"Neutron Sector = "<<nSector<<endl;
                        // cerr<<"Neutron Z Hit = "<<v_hit_3v.Z()<<endl;
                    }
                }
            }

            for (int row = 0; row < c12->getBank(cnd_hits)->getRows(); row++)
            {
                int hit_sector = c12->getBank(cnd_hits)->getInt(cnd_hit_sector, row);
                int hit_layer = c12->getBank(cnd_hits)->getInt(cnd_hit_layer, row);
                double hit_energy = c12->getBank(cnd_hits)->getFloat(cnd_hit_energy, row);

                int sdiff = nSector - hit_sector;
                // sdiff normalization
                if (sdiff <= -12)
                {
                    sdiff += 24;
                }
                else if (sdiff > 12)
                {
                    sdiff -= 24;
                }
                int ldiff = detINTlayer - hit_layer;

                if ((ldiff == 1) && (sdiff == 0))
                {
                    if (isGN)
                    {
                        h_Edep_infront_goodN_Step5->Fill(hit_energy, weight);
                    }
                    else if (isBN)
                    {
                        h_Edep_infront_badN_Step5->Fill(hit_energy, weight);
                    }
                }

                else if ((ldiff == -1) && (sdiff == 0))
                {
                    if (isGN)
                    {
                        h_Edep_behind_goodN_Step5->Fill(hit_energy, weight);
                    }
                    else if (isBN)
                    {
                        h_Edep_behind_badN_Step5->Fill(hit_energy, weight);
                    }
                }
            }
            for (int row = 0; row < c12->getBank(ctof_hits)->getRows(); row++)
            {
                int hit_sector = (c12->getBank(ctof_hits)->getInt(ctof_hit_component, row) + 1) / 2;
                double hit_energy = c12->getBank(ctof_hits)->getFloat(ctof_hit_energy, row);

                int sdiff = nSector - hit_sector;
                // sdiff normalization
                if (sdiff <= -12)
                {
                    sdiff += 24;
                }
                else if (sdiff > 12)
                {
                    sdiff -= 24;
                }
                int ldiff = detINTlayer;

                if ((ldiff == 1) && (sdiff == 0))
                {
                    if (isGN)
                    {
                        h_Edep_infront_goodN_Step5->Fill(hit_energy, weight);
                    }
                    else if (isBN)
                    {
                        h_Edep_infront_badN_Step5->Fill(hit_energy, weight);
                    }
                }

                else if ((ldiff == -1) && (sdiff == 0))
                {
                    if (isGN)
                    {
                        h_Edep_behind_goodN_Step5->Fill(hit_energy, weight);
                    }
                    else if (isBN)
                    {
                        h_Edep_behind_badN_Step5->Fill(hit_energy, weight);
                    }
                }
            }
            */

#pragma endregion /* Step Five - end */

            //////////////////////////////////////////////
            // Step Six?
            //////////////////////////////////////////////

#pragma region /* Step Six? */

            /*
            if (!isGN)
            {
                h_ToF_badN->Fill(ToF, weight);
                h_ToF_z_badN->Fill(ToF, v_hit_3v.Z(), weight);
                if (ToF < 10)
                {
                    h_TM_badN->Fill(ToF / path * 100, weight);
                    h_beta_badN->Fill(beta, weight);
                    h_mom_badN->Fill(P_n_3v.Mag(), weight);
                    h_Edep_z_badN->Fill(Edep_single, v_hit_3v.Z(), weight);
                    h_Edep_ToF_badN->Fill(Edep_single, ToF, weight);
                    h_beta_z_badN->Fill(beta, v_hit_3v.Z(), weight);
                    if (C1)
                    {
                        h_ToFc_z_1_badN->Fill(ToF * c, v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_1_badN->Fill(ToF * c - v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_Edep_1_badN->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                        h_Edep_z_1_badN->Fill(Edep_single, v_hit_3v.Z(), weight);
                    }
                    if (C2)
                    {
                        h_ToFc_z_2_badN->Fill(ToF * c, v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_2_badN->Fill(ToF * c - v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_Edep_2_badN->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                        h_Edep_z_2_badN->Fill(Edep_single, v_hit_3v.Z(), weight);
                    }
                    if (C3)
                    {
                        h_ToFc_z_3_badN->Fill(ToF * c, v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_3_badN->Fill(ToF * c - v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_Edep_3_badN->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                        h_Edep_z_3_badN->Fill(Edep_single, v_hit_3v.Z(), weight);
                    }

                    h_Edep_mom_badN->Fill(edep, mom, weight);

                    if (v_hit_3v.Z() > 10)
                    {
                        cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
                        cerr << "Run=" << c12->runconfig()->getRun() << endl;
                        cerr << "Event=" << c12->runconfig()->getEvent() << endl;
                        cerr << "Neutron Sector = " << nSector << endl;
                        cerr << "Neutron Z Hit = " << v_hit_3v.Z() << endl;
                        cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
                             << endl
                             << endl;
                    }
                }
            }
            else
            {
                h_ToF_goodN->Fill(ToF, weight);
                h_ToF_z_goodN->Fill(ToF, v_hit_3v.Z(), weight);
                if (ToF < 10)
                {
                    h_TM_goodN->Fill(ToF / path, weight);
                    h_beta_goodN->Fill(beta, weight);
                    h_mom_goodN->Fill(P_n_3v.Mag(), weight);
                    h_Edep_z_goodN->Fill(Edep_single, v_hit_3v.Z(), weight);
                    h_Edep_ToF_goodN->Fill(Edep_single, ToF, weight);
                    h_beta_z_goodN->Fill(beta, v_hit_3v.Z(), weight);
                    if (C1)
                    {
                        h_ToFc_z_1_goodN->Fill(ToF * c, v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_1_goodN->Fill(ToF * c - v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_Edep_1_goodN->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                        h_Edep_z_1_goodN->Fill(Edep_single, v_hit_3v.Z(), weight);
                    }
                    if (C2)
                    {
                        h_ToFc_z_2_goodN->Fill(ToF * c, v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_2_goodN->Fill(ToF * c - v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_Edep_2_goodN->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                        h_Edep_z_2_goodN->Fill(Edep_single, v_hit_3v.Z(), weight);
                    }
                    if (C3)
                    {
                        h_ToFc_z_3_goodN->Fill(ToF * c, v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_3_goodN->Fill(ToF * c - v_hit_3v.Z(), weight);
                        h_diff_ToFc_z_Edep_3_goodN->Fill(ToF * c - v_hit_3v.Z(), edep, weight);
                        h_Edep_z_3_goodN->Fill(Edep_single, v_hit_3v.Z(), weight);
                    }
                }
                h_Edep_mom_goodN->Fill(edep, mom, weight);
            }

            if (edep < 12.5)
            {
                continue;
            }
            if ((edep < -(40.0 / 110.0) * ((ToF * c - v_hit_3v.Z()) - 110)) && C1)
            {
                continue;
            }
            if ((edep < -(32.0 / 110.0) * ((ToF * c - v_hit_3v.Z()) - 110)) && C2)
            {
                continue;
            }
            if ((edep < -(26.0 / 110.0) * ((ToF * c - v_hit_3v.Z()) - 110)) && C3)
            {
                continue;
            }
            */

            // if (C3 && (v_hit_3v.Z() > 25))
            // {
            //     continue;
            // }
            // else if (C2 && (v_hit_3v.Z() > 20))
            // {
            //     continue;
            // }
            // else if (C1 && (v_hit_3v.Z() > 10))
            // {
            //     continue;
            // }

#pragma endregion /* Step Six? */
        }  // End of Andrew's loop over all particles

#pragma endregion /* Neutrons (Andrew) */

#pragma region /* Counters - start */

        histograms.UpdateMultiplicityHistograms(
            pInCD, pInFD, counter_n_multiplicity_allN_epCDn_Step0, counter_n_multiplicity_goodN_epCDn_Step0, counter_n_multiplicity_badN_epCDn_Step0, counter_n_multiplicity_allN_epCDn_Step1,
            counter_n_multiplicity_goodN_epCDn_Step1, counter_n_multiplicity_badN_epCDn_Step1, counter_n_multiplicity_allN_epCDn_Step2, counter_n_multiplicity_goodN_epCDn_Step2,
            counter_n_multiplicity_badN_epCDn_Step2, counter_n_multiplicity_allN_epCDn_Step3, counter_n_multiplicity_goodN_epCDn_Step3, counter_n_multiplicity_badN_epCDn_Step3,
            counter_n_multiplicity_allN_epCDn_Step4, counter_n_multiplicity_goodN_epCDn_Step4, counter_n_multiplicity_badN_epCDn_Step4, counter_n_multiplicity_allN_epCDn_Step5,
            counter_n_multiplicity_goodN_epCDn_Step5, counter_n_multiplicity_badN_epCDn_Step5, counter_n_multiplicity_allN_epFDn_Step0, counter_n_multiplicity_goodN_epFDn_Step0,
            counter_n_multiplicity_badN_epFDn_Step0, counter_n_multiplicity_allN_epFDn_Step1, counter_n_multiplicity_goodN_epFDn_Step1, counter_n_multiplicity_badN_epFDn_Step1,
            counter_n_multiplicity_allN_epFDn_Step2, counter_n_multiplicity_goodN_epFDn_Step2, counter_n_multiplicity_badN_epFDn_Step2, counter_n_multiplicity_allN_epFDn_Step3,
            counter_n_multiplicity_goodN_epFDn_Step3, counter_n_multiplicity_badN_epFDn_Step3, counter_n_multiplicity_allN_epFDn_Step4, counter_n_multiplicity_goodN_epFDn_Step4,
            counter_n_multiplicity_badN_epFDn_Step4, counter_n_multiplicity_allN_epFDn_Step5, counter_n_multiplicity_goodN_epFDn_Step5, counter_n_multiplicity_badN_epFDn_Step5, weight);

        // Count events passing steps
        if (pass_step0_cuts) { ++counter_pass_step0_cuts; }

        if (pass_step1_cuts) { ++counter_pass_step1_cuts; }

        if (pass_step2_cuts) { ++counter_pass_step2_cuts; }

        if (pass_step3_cuts) { ++counter_pass_step3_cuts; }

        if (pass_step4_cuts) { ++counter_pass_step4_cuts; }

        if (pass_step5_cuts) { ++counter_pass_step5_cuts; }

#pragma endregion /* Counters - end */

#pragma endregion /* Andrew's manual work - end */
    }  // closes event loop

#pragma endregion /* Chain loop - end */

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=
    // Wrap up
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=

#pragma region /* Wrap up - start */

    // HistPrinter(HistoList, PDFFile);
    histograms.PlotAndSaveHstograms(PDFFile);

#pragma endregion /* Wrap up - end */

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=
    // Save log file
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=

#pragma region /* Save log file - start */

    // TODO: doesn't work - fix this!

    // Saving setup to log file ------------------------------------------------------------------------------------------------------------------------------------------

    // Saving setup to log file
    ofstream myLogFile;
    myLogFile.open(("./" + OutDir + "/Log_file.txt").c_str());

    myLogFile << "///////////////////////////////////////////////////////////////////////////\n";
    myLogFile << "// Input file was " << input_hipo << "\n";
    myLogFile << "// Beam energy was" << Ebeam << "\n";
    myLogFile << "///////////////////////////////////////////////////////////////////////////\n\n";

    // myLogFile << "Run_Erins_features:\t" << Run_Erins_features << "\n";
    // myLogFile << "Run_Andrews_work:\t" << Run_Andrews_work << "\n\n";

    myLogFile << "Total #(events) in sample:\t" << EventCounter << "\n\n";

    myLogFile << "Total #((e,e'pXn) events):\t" << counter_epXn << "\n\n";

    myLogFile << "Total #((e,e'pXn) events) pass_step0_cuts:\t" << counter_pass_step0_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step1_cuts:\t" << counter_pass_step1_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step2_cuts:\t" << counter_pass_step2_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step3_cuts:\t" << counter_pass_step3_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step4_cuts:\t" << counter_pass_step4_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step5_cuts:\t" << counter_pass_step5_cuts << "\n\n";

    myLogFile.close();

#pragma endregion /* Save log file - end */

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=
    // Printouts
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------=

#pragma region /* Printouts 2 - start */

    cout << "\033[33m\n\033[0m";
    cout << "\033[33minput_hipo:\033[0m\t\t" << input_hipo << "\n";
    cout << "\033[33m\n\033[0m";
    cout << "\033[33mOUTDIR:\033[0m\t\t\t" << gSystem->Getenv("OUTDIR") << "\n";
    cout << "\033[33mOutDir:\033[0m\t\t\t" << OutDir << "\n";
    cout << "\033[33moutput_pdf_Erin:\033[0m\t" << output_pdf_Erin << "\n";
    cout << "\033[33moutput_root_Erin:\033[0m\t" << output_root_Erin << "\n";
    cout << "\033[33moutput_txt_Erin:\033[0m\t" << output_txt_Erin << "\n";
    cout << "\033[33mPDFFile:\033[0m\t\t" << PDFFile << "\n\n\n";

    /* Timing output */
    auto Code_end_time = std::chrono::system_clock::now();
    auto Elapsed_time_seconds = std::chrono::duration_cast<std::chrono::seconds>(Code_end_time - Code_start_time);
    double Elapsed_time_minutes = Elapsed_time_seconds.count() / 60;

    if (Elapsed_time_seconds.count() < 60) {
        std::cout << "\033[33mRunning time:\033[0m\t\t" << Elapsed_time_seconds.count() << " seconds\n\n";
    } else {
        std::cout << "\033[33mRunning time:\033[0m\t\t" << to_string_with_precision(Elapsed_time_minutes, 3) << " minutes\n\n";
    }

#pragma endregion /* Printouts 2 - end */

    return 0;
}  // closes main function

#pragma endregion /* ManualVeto_ToF_n_ep - end */
