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

#include "clas12reader.h"
#include "HipoChain.h"

#include "src/constants.h"
// #include "src/Histograms/ManualVetoHistograms.h"
#include "src/Histograms/HistPrinter.h"
#include "src/functions/GeneralFunctions.h"
#include "src/functions/NeutronFunctions.h"
#include "src/functions/neutron-veto/veto_functions.cpp"
#include "src/functions/Andrews_functions/Andrews_functions.cpp"
#include "src/functions/HipoChain_config.cpp"
#include "src/classes/clas12ana/clas12ana.cpp"

using namespace std;
using namespace clas12;

int ManualVeto_Phase6(                                                                                //
    const string OutDir, string output_pdf_Erin,                                                      // My arguments
    double Ebeam, bool keep_good, string output_root_Erin, string output_txt_Erin, string input_hipo, // Erin's arguments
    string PDFFile, int isMC = 0                                                                      // Andrew's arguments
)
// int main(int argc, char **argv)
{
    auto Code_start_time = std::chrono::system_clock::now(); // Start counting running time

    // ======================================================================================================================================================================
    // Printouts
    // ======================================================================================================================================================================

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

    // ======================================================================================================================================================================
    // Initial setup
    // ======================================================================================================================================================================

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

    clasAna->readEcalSFPar("src/cuts/paramsSF_LD2_x2.dat"); // TODO: check if applied
    clasAna->readEcalPPar("src/cuts/paramsPI_LD2_x2.dat");  // TODO: check if applied

    clasAna->setProtonPidCuts(true);

#pragma endregion /* Initial setup - end */

    // ======================================================================================================================================================================
    // Veto histograms
    // ======================================================================================================================================================================

#pragma region /* Veto histograms - start */

    /////////////////////////////////////
    // Prepare histograms
    /////////////////////////////////////

    vector<TH1 *> HistoList;

    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);

    gStyle->SetTitleXOffset(0.8);
    gStyle->SetTitleYOffset(0.8);

    char temp_name[100];
    char temp_title[100];

    // (e,e'p) plots
    // ======================================================================================================================================================================

    /* Proton histograms (from Erin) */
    TH1D *h_p_multiplicity_BPID_epCD = new TH1D("p_multiplicity_BPID_epCD", "Number of CD Protons in Event (Before PID);CD Proton multiplicity;Counts", 10, 0, 10);
    HistoList.push_back(h_p_multiplicity_BPID_epCD);
    TH1D *h_P_p_BPID_epCD = new TH1D("P_p_BPID_epCD", "CD Proton momentum (Before PID);P_{p} [GeV/C];Counts", 50, 0., 3.);
    HistoList.push_back(h_P_p_BPID_epCD);
    TH2D *h_theta_p_VS_phi_p_BPID_epCD = new TH2D("theta_p_VS_phi_p_BPID_epCD", "#theta_{p} vs #phi_{p} of CD proton (Before PID);#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
    HistoList.push_back(h_theta_p_VS_phi_p_BPID_epCD);
    TH1D *h_p_multiplicity_APID_epCD = new TH1D("p_multiplicity_APID_epCD", "Number of CD Protons in Event (After PID);CD Proton multiplicity;Counts", 10, 0, 10);
    HistoList.push_back(h_p_multiplicity_APID_epCD);
    TH1D *h_P_p_APID_epCD = new TH1D("P_p_APID_epCD", "CD Proton momentum (After PID);P_{p} [GeV/C];Counts", 50, 0., 3.);
    HistoList.push_back(h_P_p_APID_epCD);
    TH2D *h_theta_p_VS_phi_p_APID_epCD = new TH2D("theta_p_VS_phi_p_APID_epCD", "#theta_{p} vs #phi_{p} of CD proton (After PID);#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
    HistoList.push_back(h_theta_p_VS_phi_p_APID_epCD);

    TH1D *h_p_multiplicity_BPID_epFD = new TH1D("p_multiplicity_BPID_epFD", "Number of FD Protons in Event (Before PID);FD Proton multiplicity;Counts", 10, 0, 10);
    HistoList.push_back(h_p_multiplicity_BPID_epFD);
    TH1D *h_P_p_BPID_epFD = new TH1D("P_p_BPID_epFD", "FD Proton momentum (Before PID);P_{p} [GeV/C];Counts", 50, 0., 3.);
    HistoList.push_back(h_P_p_BPID_epFD);
    TH2D *h_theta_p_VS_phi_p_BPID_epFD = new TH2D("theta_p_VS_phi_p_BPID_epFD", "#theta_{p} vs #phi_{p} of FD proton (Before PID);#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
    HistoList.push_back(h_theta_p_VS_phi_p_BPID_epFD);
    TH1D *h_p_multiplicity_APID_epFD = new TH1D("p_multiplicity_APID_epFD", "Number of FD Protons in Even (After PID);FD Proton multiplicity;Counts", 10, 0, 10);
    HistoList.push_back(h_p_multiplicity_APID_epFD);
    TH1D *h_P_p_APID_epFD = new TH1D("P_p_APID_epFD", "FD Proton momentum (After PID);P_{p} [GeV/C];Counts", 50, 0., 3.);
    HistoList.push_back(h_P_p_APID_epFD);
    TH2D *h_theta_p_VS_phi_p_APID_epFD = new TH2D("theta_p_VS_phi_p_APID_epFD", "#theta_{p} vs #phi_{p} of FD proton (After PID);#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
    HistoList.push_back(h_theta_p_VS_phi_p_APID_epFD);

    TH2D *h_dbeta_p_VS_P_p_BPID_epCD = new TH2D("dbeta_p_VS_P_p_BPID_epCD", "#Delta#beta_{p} vs CD proton momentum (Before PID);P_{p} [GeV/c];#Delta#beta_{p}", 50, 0, 3, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_p_VS_P_p_BPID_epCD);
    TH1D *h_dVz_p_BPID_epCD = new TH1D("dVz_p_BPID_epCD", "Vertex correlation between CD proton and electron (Before PID);dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
    HistoList.push_back(h_dVz_p_BPID_epCD);
    TH1D *h_Chi2pid_p_BPID_epCD = new TH1D("Chi2pid_p_BPID_epCD", "CD Proton #chi^{2}_{p} (Before PID);#chi^{2}_{p};Counts", 50, -6, 6);
    HistoList.push_back(h_Chi2pid_p_BPID_epCD);
    TH2D *h_dbeta_p_VS_P_p_APID_epCD = new TH2D("dbeta_p_VS_P_p_APID_epCD", "#Delta#beta_{p} vs CD proton momentum (After PID);P_{p} [GeV/c];#Delta#beta_{p}", 50, 0, 3, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_p_VS_P_p_APID_epCD);
    TH1D *h_dVz_p_APID_epCD = new TH1D("dVz_p_APID_epCD", "Vertex correlation between CD proton and electron (After PID);dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
    HistoList.push_back(h_dVz_p_APID_epCD);
    TH1D *h_Chi2pid_p_APID_epCD = new TH1D("Chi2pid_p_APID_epCD", "CD Proton #chi^{2}_{p} (After PID);#chi^{2}_{p};Counts", 50, -6, 6);
    HistoList.push_back(h_Chi2pid_p_APID_epCD);

    TH2D *h_dbeta_p_VS_P_p_BPID_epFD = new TH2D("dbeta_p_VS_P_p_BPID_epFD", "#Delta#beta_{p} vs FD proton momentum (Before PID);P_{p} [GeV/c];#Delta#beta_{p}", 50, 0, 3, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_p_VS_P_p_BPID_epFD);
    TH1D *h_dVz_p_BPID_epFD = new TH1D("dVz_p_BPID_epFD", "Vertex correlation between FD proton and electron (Before PID);dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
    HistoList.push_back(h_dVz_p_BPID_epFD);
    TH1D *h_Chi2pid_p_BPID_epFD = new TH1D("Chi2pid_p_BPID_epFD", "FD Proton #chi^{2}_{p} (Before PID);#chi^{2}_{p};Counts", 50, -6, 6);
    HistoList.push_back(h_Chi2pid_p_BPID_epFD);
    TH2D *h_dbeta_p_VS_P_p_APID_epFD = new TH2D("dbeta_p_VS_P_p_APID_epFD", "#Delta#beta_{p} vs FD proton momentum (After PID);P_{p} [GeV/c];#Delta#beta_{p}", 50, 0, 3, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_p_VS_P_p_APID_epFD);
    TH1D *h_dVz_p_APID_epFD = new TH1D("dVz_p_APID_epFD", "Vertex correlation between FD proton and electron (After PID);dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
    HistoList.push_back(h_dVz_p_APID_epFD);
    TH1D *h_Chi2pid_p_APID_epFD = new TH1D("Chi2pid_p_APID_epFD", "FD Proton #chi^{2}_{p} (After PID);#chi^{2}_{p};Counts", 50, -6, 6);
    HistoList.push_back(h_Chi2pid_p_APID_epFD);

    /* Missing variabels */
    TH1D *h_P_miss_BmissC_epCD = new TH1D("P_miss_BmissC_epCD", "Missing Momentum (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);P_{miss} [GeV/c]", 50, 0, 1.7);
    HistoList.push_back(h_P_miss_BmissC_epCD);
    TH1D *h_theta_miss_BmissC_epCD = new TH1D("theta_miss_BmissC_epCD", "Missing Momentum (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);#theta_{miss} [#circ]", 50, 0., 180.);
    HistoList.push_back(h_theta_miss_BmissC_epCD);
    TH2D *h_P_miss_VS_theta_miss_BmissC_epCD = new TH2D("P_miss_VS_theta_miss_BmissC_epCD", "Missing Momentum vs #theta_{miss} (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.7);
    HistoList.push_back(h_P_miss_VS_theta_miss_BmissC_epCD);
    TH1D *h_P_miss_AmissC_epCD = new TH1D("P_miss_AmissC_epCD", "Missing Momentum (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);P_{miss} [GeV/c]", 50, 0, 1.7);
    HistoList.push_back(h_P_miss_AmissC_epCD);
    TH1D *h_theta_miss_AmissC_epCD = new TH1D("theta_miss_AmissC_epCD", "Missing Momentum (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);#theta_{miss} [#circ]", 50, 0., 180.);
    HistoList.push_back(h_theta_miss_AmissC_epCD);
    TH2D *h_P_miss_VS_theta_miss_AmissC_epCD = new TH2D("P_miss_VS_theta_miss_AmissC_epCD", "Missing Momentum vs #theta_{miss} (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.7);
    HistoList.push_back(h_P_miss_VS_theta_miss_AmissC_epCD);

    TH1D *h_P_miss_BmissC_epFD = new TH1D("P_miss_BmissC_epFD", "Missing Momentum (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);P_{miss} [GeV/c]", 50, 0, 1.7);
    HistoList.push_back(h_P_miss_BmissC_epFD);
    TH1D *h_theta_miss_BmissC_epFD = new TH1D("theta_miss_BmissC_epFD", "Missing Momentum (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);#theta_{miss} [#circ]", 50, 0., 180.);
    HistoList.push_back(h_theta_miss_BmissC_epFD);
    TH2D *h_P_miss_VS_theta_miss_BmissC_epFD = new TH2D("P_miss_VS_theta_miss_BmissC_epFD", "Missing Momentum vs #theta_{miss} (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.7);
    HistoList.push_back(h_P_miss_VS_theta_miss_BmissC_epFD);
    TH1D *h_P_miss_AmissC_epFD = new TH1D("P_miss_AmissC_epFD", "Missing Momentum (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);P_{miss} [GeV/c]", 50, 0, 1.7);
    HistoList.push_back(h_P_miss_AmissC_epFD);
    TH1D *h_theta_miss_AmissC_epFD = new TH1D("theta_miss_AmissC_epFD", "Missing Momentum (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);#theta_{miss} [#circ]", 50, 0., 180.);
    HistoList.push_back(h_theta_miss_AmissC_epFD);
    TH2D *h_P_miss_VS_theta_miss_AmissC_epFD = new TH2D("P_miss_VS_theta_miss_AmissC_epFD", "Missing Momentum vs #theta_{miss} (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.7);
    HistoList.push_back(h_P_miss_VS_theta_miss_AmissC_epFD);

    TH1D *h_E_p_BmissC_epCD = new TH1D("E_p_BmissC_epCD", "CD Proton Energy (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_BmissC_epCD);
    TH1D *h_E_miss_BmissC_epCD = new TH1D("E_miss_BmissC_epCD", "Missing Energy (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_BmissC_epCD);
    TH1D *h_M_miss_BmissC_epCD = new TH1D("M_miss_BmissC_epCD", "Missing Mass (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_BmissC_epCD);
    TH1D *h_E_p_AmissC_epCD = new TH1D("E_p_AmissC_epCD", "CD Proton Energy (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_AmissC_epCD);
    TH1D *h_E_miss_AmissC_epCD = new TH1D("E_miss_AmissC_epCD", "Missing Energy (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_AmissC_epCD);
    TH1D *h_M_miss_AmissC_epCD = new TH1D("M_miss_AmissC_epCD", "Missing Mass (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_AmissC_epCD);

    TH1D *h_E_p_BmissC_epFD = new TH1D("E_p_BmissC_epFD", "FD Proton Energy (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_BmissC_epFD);
    TH1D *h_E_miss_BmissC_epFD = new TH1D("E_miss_BmissC_epFD", "Missing Energy (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_BmissC_epFD);
    TH1D *h_M_miss_BmissC_epFD = new TH1D("M_miss_BmissC_epFD", "Missing Mass (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_BmissC_epFD);
    TH1D *h_E_p_AmissC_epFD = new TH1D("E_p_AmissC_epFD", "FD Proton Energy (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_AmissC_epFD);
    TH1D *h_E_miss_AmissC_epFD = new TH1D("E_miss_AmissC_epFD", "Missing Energy (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_AmissC_epFD);
    TH1D *h_M_miss_AmissC_epFD = new TH1D("M_miss_AmissC_epFD", "Missing Mass (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_AmissC_epFD);

    /* Checks on which events have neutrons (Andrew) */
    TH1D *h_xB_BmissC_epCD = new TH1D("xB_BmissC_epCD", "x_{B} Distribution (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);x_{B}", 100, 0.0, 2.0);
    HistoList.push_back(h_xB_BmissC_epCD);
    TH2D *h_xB_VS_M_miss_BmissC_epCD = new TH2D("xB_VS_M_miss_BmissC_epCD", "x_{B} vs. M_{miss} (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_BmissC_epCD);
    TH1D *h_xB_AmissC_epCD = new TH1D("xB_AmissC_epCD", "x_{B} Distribution (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);x_{B}", 100, 0.0, 2.0);
    HistoList.push_back(h_xB_AmissC_epCD);
    TH2D *h_xB_VS_M_miss_AmissC_epCD = new TH2D("xB_VS_M_miss_AmissC_epCD", "x_{B} vs. M_{miss} (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_AmissC_epCD);

    TH1D *h_xB_BmissC_epFD = new TH1D("xB_BmissC_epFD", "x_{B} Distribution (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);x_{B}", 100, 0.0, 2.0);
    HistoList.push_back(h_xB_BmissC_epFD);
    TH2D *h_xB_VS_M_miss_BmissC_epFD = new TH2D("xB_VS_M_miss_BmissC_epFD", "x_{B} vs. M_{miss} (Before P_{miss}, #theta_{miss}, and M_{miss} Cuts);x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_BmissC_epFD);
    TH1D *h_xB_AmissC_epFD = new TH1D("xB_AmissC_epFD", "x_{B} Distribution (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);x_{B}", 100, 0.0, 2.0);
    HistoList.push_back(h_xB_AmissC_epFD);
    TH2D *h_xB_VS_M_miss_AmissC_epFD = new TH2D("xB_VS_M_miss_AmissC_epFD", "x_{B} vs. M_{miss} (After P_{miss}, #theta_{miss}, and M_{miss} Cuts);x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_AmissC_epFD);

    TH2D *h_xB_VS_M_miss_epCDn = new TH2D("xB_VS_M_miss_epCDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 50, 0.0, 2.0, 50, 0.5, 1.7);
    HistoList.push_back(h_xB_VS_M_miss_epCDn);
    TH2D *h_xB_VS_M_miss_epFDn = new TH2D("xB_VS_M_miss_epFDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 50, 0.0, 2.0, 50, 0.5, 1.7);
    HistoList.push_back(h_xB_VS_M_miss_epFDn);

    TH2D *h_xB_VS_M_miss_goodN_epCDn = new TH2D("xB_VS_M_miss_goodN_epCDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 50, 0.0, 2.0, 50, 0.5, 1.7);
    HistoList.push_back(h_xB_VS_M_miss_goodN_epCDn);
    TH2D *h_xB_VS_M_miss_badN_epCDn = new TH2D("xB_VS_M_miss_badN_epCDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 50, 0.0, 2.0, 50, 0.5, 1.7);
    HistoList.push_back(h_xB_VS_M_miss_badN_epCDn);

    TH2D *h_xB_VS_M_miss_goodN_epFDn = new TH2D("xB_VS_M_miss_goodN_epFDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 50, 0.0, 2.0, 50, 0.5, 1.7);
    HistoList.push_back(h_xB_VS_M_miss_goodN_epFDn);
    TH2D *h_xB_VS_M_miss_badN_epFDn = new TH2D("xB_VS_M_miss_badN_epFDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 50, 0.0, 2.0, 50, 0.5, 1.7);
    HistoList.push_back(h_xB_VS_M_miss_badN_epFDn);

    /* Kinematical variables */
    TH1D *h_theta_n_epCDn = new TH1D("theta_n_epCDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_epCDn);
    TH1D *h_phi_n_epCDn = new TH1D("phi_n_epCDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_epCDn);
    TH2D *h_theta_n_VS_phi_n_epCDn = new TH2D("theta_n_VS_phi_n_epCDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_epCDn);
    TH2D *h_theta_n_VS_beta_n_epCDn = new TH2D("theta_VS_beta_epCDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_epCDn);

    TH1D *h_theta_n_epFDn = new TH1D("theta_n_epFDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_epFDn);
    TH1D *h_phi_n_epFDn = new TH1D("phi_n_epFDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_epFDn);
    TH2D *h_theta_n_VS_phi_n_epFDn = new TH2D("theta_n_VS_phi_n_epFDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_epFDn);
    TH2D *h_theta_n_VS_beta_n_epFDn = new TH2D("theta_VS_beta_epFDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_epFDn);

    TH1D *h_P_n_epCDn = new TH1D("P_n_epCDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_epCDn);
    TH2D *h_P_n_VS_theta_n_epCDn = new TH2D("P_n_VS_theta_n_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_epCDn);

    TH1D *h_P_n_epFDn = new TH1D("P_n_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_epFDn);
    TH2D *h_P_n_VS_theta_n_epFDn = new TH2D("P_n_VS_theta_n_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_epFDn);

    TH1D *h_P_miss_epCDn = new TH1D("P_miss_epCDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.7);
    HistoList.push_back(h_P_miss_epCDn);
    TH2D *h_P_miss_VS_theta_miss_epCDn = new TH2D("P_miss_VS_theta_miss_epCDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.7);
    HistoList.push_back(h_P_miss_VS_theta_miss_epCDn);

    TH1D *h_P_miss_epFDn = new TH1D("P_miss_epFDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.7);
    HistoList.push_back(h_P_miss_epFDn);
    TH2D *h_P_miss_VS_theta_miss_epFDn = new TH2D("P_miss_VS_theta_miss_epFDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.7);
    HistoList.push_back(h_P_miss_VS_theta_miss_epFDn);

    TH1D *h_dpp_allN_epCDn = new TH1D("dpp_allN_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_allN_epCDn);
    TH1D *h_dpp_goodN_epCDn = new TH1D("dpp_goodN_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_goodN_epCDn);
    TH1D *h_dpp_badN_epCDn = new TH1D("dpp_badN_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_badN_epCDn);

    TH1D *h_dpp_allN_epFDn = new TH1D("dpp_allN_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_allN_epFDn);
    TH1D *h_dpp_goodN_epFDn = new TH1D("dpp_goodN_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_goodN_epFDn);
    TH1D *h_dpp_badN_epFDn = new TH1D("dpp_badN_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_badN_epFDn);

    TH1D *h_theta_n_miss_allN_epCDn = new TH1D("theta_n_miss_allN_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_allN_epCDn);
    TH1D *h_theta_n_miss_goodN_epCDn = new TH1D("theta_n_miss_goodN_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_goodN_epCDn);
    TH1D *h_theta_n_miss_badN_epCDn = new TH1D("theta_n_miss_badN_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_badN_epCDn);

    TH1D *h_theta_n_miss_allN_epFDn = new TH1D("theta_n_miss_allN_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_allN_epFDn);
    TH1D *h_theta_n_miss_goodN_epFDn = new TH1D("theta_n_miss_goodN_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_goodN_epFDn);
    TH1D *h_theta_n_miss_badN_epFDn = new TH1D("theta_n_miss_badN_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_badN_epFDn);

    TH2D *h_dpp_VS_theta_n_miss_epCDn = new TH2D("dpp_VS_theta_n_miss_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -1.5, 1.5, 50, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_epCDn);

    TH2D *h_dpp_VS_theta_n_miss_epFDn = new TH2D("dpp_VS_theta_n_miss_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -1.5, 1.5, 50, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_epFDn);

    TH1D *h_E_p_epCDn = new TH1D("E_p_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_epCDn);
    TH1D *h_E_miss_epCDn = new TH1D("E_miss_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_epCDn);
    TH1D *h_M_miss_epCDn = new TH1D("M_miss_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_epCDn);
    TH2D *h_M_miss_VS_P_n_epCDn = new TH2D("M_miss_VS_P_n_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_epCDn);
    TH2D *h_M_miss_VS_theta_n_epCDn = new TH2D("M_miss_VS_theta_n_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_epCDn);
    TH2D *h_M_miss_VS_phi_n_epCDn = new TH2D("M_miss_VS_phi_n_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_epCDn);
    TH2D *h_M_miss_VS_P_miss_epCDn = new TH2D("M_miss_VS_P_miss_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_epCDn);
    TH2D *h_M_miss_VS_theta_miss_epCDn = new TH2D("M_miss_VS_theta_miss_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_epCDn);
    TH2D *h_M_miss_VS_phi_miss_epCDn = new TH2D("M_miss_VS_phi_miss_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_epCDn);

    TH1D *h_E_p_epFDn = new TH1D("E_p_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_epFDn);
    TH1D *h_E_miss_epFDn = new TH1D("E_miss_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_epFDn);
    TH1D *h_M_miss_epFDn = new TH1D("M_miss_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_epFDn);
    TH2D *h_M_miss_VS_P_n_epFDn = new TH2D("M_miss_VS_P_n_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_epFDn);
    TH2D *h_M_miss_VS_theta_n_epFDn = new TH2D("M_miss_VS_theta_n_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_epFDn);
    TH2D *h_M_miss_VS_phi_n_epFDn = new TH2D("M_miss_VS_phi_n_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_epFDn);
    TH2D *h_M_miss_VS_P_miss_epFDn = new TH2D("M_miss_VS_P_miss_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_epFDn);
    TH2D *h_M_miss_VS_theta_miss_epFDn = new TH2D("M_miss_VS_theta_miss_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_epFDn);
    TH2D *h_M_miss_VS_phi_miss_epFDn = new TH2D("M_miss_VS_phi_miss_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_epFDn);

    TH1D *h_P_n_minus_P_miss_epCDn = new TH1D("P_n_minus_P_miss_epCDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_epCDn);
    TH1D *h_P_n_x_minus_P_miss_x_epCDn = new TH1D("P_n_x_minus_P_miss_x_epCDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_epCDn);
    TH1D *h_P_n_y_minus_P_miss_y_epCDn = new TH1D("P_n_y_minus_P_miss_y_epCDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_epCDn);
    TH1D *h_P_n_z_minus_P_miss_z_epCDn = new TH1D("P_n_z_minus_P_miss_z_epCDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_epCDn);

    TH1D *h_P_n_minus_P_miss_epFDn = new TH1D("P_n_minus_P_miss_epFDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_epFDn);
    TH1D *h_P_n_x_minus_P_miss_x_epFDn = new TH1D("P_n_x_minus_P_miss_x_epFDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_epFDn);
    TH1D *h_P_n_y_minus_P_miss_y_epFDn = new TH1D("P_n_y_minus_P_miss_y_epFDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_epFDn);
    TH1D *h_P_n_z_minus_P_miss_z_epFDn = new TH1D("P_n_z_minus_P_miss_z_epFDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_epFDn);

    TH2D *h_P_n_VS_P_miss_epCDn = new TH2D("P_n_VS_P_miss_epCDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_epCDn);
    TH2D *h_P_n_x_VS_P_miss_x_epCDn = new TH2D("P_n_x_VS_P_miss_x_epCDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_epCDn);
    TH2D *h_P_n_y_VS_P_miss_y_epCDn = new TH2D("P_n_y_VS_P_miss_y_epCDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_epCDn);
    TH2D *h_P_n_z_VS_P_miss_z_epCDn = new TH2D("P_n_z_VS_P_miss_z_epCDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_epCDn);

    TH2D *h_P_n_VS_P_miss_epFDn = new TH2D("P_n_VS_P_miss_epFDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_epFDn);
    TH2D *h_P_n_x_VS_P_miss_x_epFDn = new TH2D("P_n_x_VS_P_miss_x_epFDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_epFDn);
    TH2D *h_P_n_y_VS_P_miss_y_epFDn = new TH2D("P_n_y_VS_P_miss_y_epFDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_epFDn);
    TH2D *h_P_n_z_VS_P_miss_z_epFDn = new TH2D("P_n_z_VS_P_miss_z_epFDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_epFDn);

    TH1D *h_theta_n_p_epCDn = new TH1D("theta_n_p_epCDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_epCDn);
    TH2D *h_theta_n_p_VS_P_p_epCDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_epCDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_epCDn);

    TH1D *h_theta_n_p_epFDn = new TH1D("theta_n_p_epFDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_epFDn);
    TH2D *h_theta_n_p_VS_P_p_epFDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_epFDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_epFDn);

    TH1D *h_xB_epCDn = new TH1D("xB_epCDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_epCDn);

    TH1D *h_xB_epFDn = new TH1D("xB_epFDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_epFDn);

    /* Detector responses */
    TH1D *h_Edep_CND_epCDn = new TH1D("Edep_CND_epCDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_epCDn);
    TH2D *h_P_n_VS_Edep_CND_epCDn = new TH2D("P_n_VS_Edep_CND_epCDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_epCDn);
    TH2D *h_theta_n_VS_Edep_CND_epCDn = new TH2D("theta_n_VS_Edep_CND_epCDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_epCDn);
    TH2D *h_phi_n_VS_Edep_CND_epCDn = new TH2D("phi_n_VS_Edep_CND_epCDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_epCDn);
    TH2D *h_P_miss_VS_Edep_CND_epCDn = new TH2D("P_miss_VS_Edep_CND_epCDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND_epCDn = new TH2D("theta_miss_VS_Edep_CND_epCDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND_epCDn = new TH2D("phi_miss_VS_Edep_CND_epCDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_epCDn);
    TH2D *h_dpp_VS_Edep_CND_epCDn = new TH2D("dpp_VS_Edep_CND_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_epCDn);
    TH2D *h_beta_n_VS_Edep_CND_epCDn = new TH2D("beta_n_VS_Edep_CND_epCDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_epCDn);
    TH2D *h_E_miss_VS_Edep_CND_epCDn = new TH2D("E_miss_VS_Edep_CND_epCDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_epCDn);
    TH2D *h_M_miss_VS_Edep_CND_epCDn = new TH2D("M_miss_VS_Edep_CND_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_epCDn);
    TH2D *h_path_VS_Edep_CND_epCDn = new TH2D("path_VS_Edep_CND_epCDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_epCDn = new TH2D("theta_n_miss_VS_Edep_CND_epCDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_epCDn);
    TH2D *h_ToF_VS_Edep_CND_epCDn = new TH2D("ToF_VS_Edep_CND_epCDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_epCDn);
    TH2D *h_nSector_VS_Edep_CND_epCDn = new TH2D("nSector_VS_Edep_CND_epCDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_epCDn = new TH2D("Edep_CND1_VS_Edep_CND_epCDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_epCDn = new TH2D("Edep_CND2_VS_Edep_CND_epCDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_epCDn = new TH2D("Edep_CND3_VS_Edep_CND_epCDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_epCDn);

    TH1D *h_Edep_CND_epFDn = new TH1D("Edep_CND_epFDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_epFDn);
    TH2D *h_P_n_VS_Edep_CND_epFDn = new TH2D("P_n_VS_Edep_CND_epFDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_epFDn);
    TH2D *h_theta_n_VS_Edep_CND_epFDn = new TH2D("theta_n_VS_Edep_CND_epFDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_epFDn);
    TH2D *h_phi_n_VS_Edep_CND_epFDn = new TH2D("phi_n_VS_Edep_CND_epFDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_epFDn);
    TH2D *h_P_miss_VS_Edep_CND_epFDn = new TH2D("P_miss_VS_Edep_CND_epFDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND_epFDn = new TH2D("theta_miss_VS_Edep_CND_epFDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND_epFDn = new TH2D("phi_miss_VS_Edep_CND_epFDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_epFDn);
    TH2D *h_dpp_VS_Edep_CND_epFDn = new TH2D("dpp_VS_Edep_CND_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_epFDn);
    TH2D *h_beta_n_VS_Edep_CND_epFDn = new TH2D("beta_n_VS_Edep_CND_epFDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_epFDn);
    TH2D *h_E_miss_VS_Edep_CND_epFDn = new TH2D("E_miss_VS_Edep_CND_epFDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_epFDn);
    TH2D *h_M_miss_VS_Edep_CND_epFDn = new TH2D("M_miss_VS_Edep_CND_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_epFDn);
    TH2D *h_path_VS_Edep_CND_epFDn = new TH2D("path_VS_Edep_CND_epFDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_epFDn = new TH2D("theta_n_miss_VS_Edep_CND_epFDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_epFDn);
    TH2D *h_ToF_VS_Edep_CND_epFDn = new TH2D("ToF_VS_Edep_CND_epFDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_epFDn);
    TH2D *h_nSector_VS_Edep_CND_epFDn = new TH2D("nSector_VS_Edep_CND_epFDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_epFDn = new TH2D("Edep_CND1_VS_Edep_CND_epFDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_epFDn = new TH2D("Edep_CND2_VS_Edep_CND_epFDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_epFDn = new TH2D("Edep_CND3_VS_Edep_CND_epFDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_epFDn);

    TH1D *h_Edep_CTOF_epCDn = new TH1D("Edep_CTOF_epCDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_epCDn);
    TH2D *h_P_n_VS_Edep_CTOF_epCDn = new TH2D("P_n_VS_Edep_CTOF_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_epCDn);
    TH2D *h_theta_n_VS_Edep_CTOF_epCDn = new TH2D("theta_n_VS_Edep_CTOF_epCDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_epCDn);
    TH2D *h_phi_n_VS_Edep_CTOF_epCDn = new TH2D("phi_n_VS_Edep_CTOF_epCDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_epCDn);
    TH2D *h_P_miss_VS_Edep_CTOF_epCDn = new TH2D("P_miss_VS_Edep_CTOF_epCDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_epCDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_epCDn = new TH2D("theta_miss_VS_Edep_CTOF_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_epCDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_epCDn = new TH2D("phi_miss_VS_Edep_CTOF_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_epCDn);
    TH2D *h_dpp_VS_Edep_CTOF_epCDn = new TH2D("dpp_VS_Edep_CTOF_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_epCDn);
    TH2D *h_beta_n_VS_Edep_CTOF_epCDn = new TH2D("beta_n_VS_Edep_CTOF_epCDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_epCDn);
    TH2D *h_E_miss_VS_Edep_CTOF_epCDn = new TH2D("E_miss_VS_Edep_CTOF_epCDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_epCDn);
    TH2D *h_M_miss_VS_Edep_CTOF_epCDn = new TH2D("M_miss_VS_Edep_CTOF_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_epCDn);
    TH2D *h_path_VS_Edep_CTOF_epCDn = new TH2D("path_VS_Edep_CTOF_epCDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_epCDn = new TH2D("theta_n_miss_VS_Edep_CTOF_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_epCDn);
    TH2D *h_ToF_VS_Edep_CTOF_epCDn = new TH2D("ToF_VS_Edep_CTOF_epCDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_epCDn);
    TH2D *h_nSector_VS_Edep_CTOF_epCDn = new TH2D("nSector_VS_Edep_CTOF_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_epCDn = new TH2D("Edep_CND1_VS_Edep_CTOF_epCDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_epCDn = new TH2D("Edep_CND2_VS_Edep_CTOF_epCDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_epCDn = new TH2D("Edep_CND3_VS_Edep_CTOF_epCDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_epCDn);

    TH1D *h_Edep_CTOF_epFDn = new TH1D("Edep_CTOF_epFDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_epFDn);
    TH2D *h_P_n_VS_Edep_CTOF_epFDn = new TH2D("P_n_VS_Edep_CTOF_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_epFDn);
    TH2D *h_theta_n_VS_Edep_CTOF_epFDn = new TH2D("theta_n_VS_Edep_CTOF_epFDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_epFDn);
    TH2D *h_phi_n_VS_Edep_CTOF_epFDn = new TH2D("phi_n_VS_Edep_CTOF_epFDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_epFDn);
    TH2D *h_P_miss_VS_Edep_CTOF_epFDn = new TH2D("P_miss_VS_Edep_CTOF_epFDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_epFDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_epFDn = new TH2D("theta_miss_VS_Edep_CTOF_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_epFDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_epFDn = new TH2D("phi_miss_VS_Edep_CTOF_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_epFDn);
    TH2D *h_dpp_VS_Edep_CTOF_epFDn = new TH2D("dpp_VS_Edep_CTOF_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_epFDn);
    TH2D *h_beta_n_VS_Edep_CTOF_epFDn = new TH2D("beta_n_VS_Edep_CTOF_epFDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_epFDn);
    TH2D *h_E_miss_VS_Edep_CTOF_epFDn = new TH2D("E_miss_VS_Edep_CTOF_epFDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_epFDn);
    TH2D *h_M_miss_VS_Edep_CTOF_epFDn = new TH2D("M_miss_VS_Edep_CTOF_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_epFDn);
    TH2D *h_path_VS_Edep_CTOF_epFDn = new TH2D("path_VS_Edep_CTOF_epFDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_epFDn = new TH2D("theta_n_miss_VS_Edep_CTOF_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_epFDn);
    TH2D *h_ToF_VS_Edep_CTOF_epFDn = new TH2D("ToF_VS_Edep_CTOF_epFDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_epFDn);
    TH2D *h_nSector_VS_Edep_CTOF_epFDn = new TH2D("nSector_VS_Edep_CTOF_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_epFDn = new TH2D("Edep_CND1_VS_Edep_CTOF_epFDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_epFDn = new TH2D("Edep_CND2_VS_Edep_CTOF_epFDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_epFDn = new TH2D("Edep_CND3_VS_Edep_CTOF_epFDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_epFDn);

    TH1D *h_Edep_single_epCDn = new TH1D("Edep_single_epCDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_epCDn);
    TH2D *h_P_n_VS_Edep_single_epCDn = new TH2D("P_n_VS_Edep_single_epCDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_epCDn);
    TH2D *h_theta_n_VS_Edep_single_epCDn = new TH2D("theta_n_VS_Edep_single_epCDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_epCDn);
    TH2D *h_phi_n_VS_Edep_single_epCDn = new TH2D("phi_n_VS_Edep_single_epCDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_epCDn);
    TH2D *h_P_miss_VS_Edep_single_epCDn = new TH2D("P_miss_VS_Edep_single_epCDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_epCDn);
    TH2D *h_theta_miss_VS_Edep_single_epCDn = new TH2D("theta_miss_VS_Edep_single_epCDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_epCDn);
    TH2D *h_phi_miss_VS_Edep_single_epCDn = new TH2D("phi_miss_VS_Edep_single_epCDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_epCDn);
    TH2D *h_dpp_VS_Edep_single_epCDn = new TH2D("dpp_VS_Edep_single_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_epCDn);
    TH2D *h_beta_n_VS_Edep_single_epCDn = new TH2D("beta_n_VS_Edep_single_epCDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_epCDn);
    TH2D *h_E_miss_VS_Edep_single_epCDn = new TH2D("E_miss_VS_Edep_single_epCDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_epCDn);
    TH2D *h_M_miss_VS_Edep_single_epCDn = new TH2D("M_miss_VS_Edep_single_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_epCDn);
    TH2D *h_path_VS_Edep_single_epCDn = new TH2D("path_VS_Edep_single_epCDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_single_epCDn = new TH2D("theta_n_miss_VS_Edep_single_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_epCDn);
    TH2D *h_ToF_VS_Edep_single_epCDn = new TH2D("ToF_VS_Edep_single_epCDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_epCDn);
    TH2D *h_nSector_VS_Edep_single_epCDn = new TH2D("nSector_VS_Edep_single_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_epCDn);

    TH1D *h_Edep_single_epFDn = new TH1D("Edep_single_epFDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_epFDn);
    TH2D *h_P_n_VS_Edep_single_epFDn = new TH2D("P_n_VS_Edep_single_epFDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_epFDn);
    TH2D *h_theta_n_VS_Edep_single_epFDn = new TH2D("theta_n_VS_Edep_single_epFDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_epFDn);
    TH2D *h_phi_n_VS_Edep_single_epFDn = new TH2D("phi_n_VS_Edep_single_epFDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_epFDn);
    TH2D *h_P_miss_VS_Edep_single_epFDn = new TH2D("P_miss_VS_Edep_single_epFDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_epFDn);
    TH2D *h_theta_miss_VS_Edep_single_epFDn = new TH2D("theta_miss_VS_Edep_single_epFDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_epFDn);
    TH2D *h_phi_miss_VS_Edep_single_epFDn = new TH2D("phi_miss_VS_Edep_single_epFDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_epFDn);
    TH2D *h_dpp_VS_Edep_single_epFDn = new TH2D("dpp_VS_Edep_single_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_epFDn);
    TH2D *h_beta_n_VS_Edep_single_epFDn = new TH2D("beta_n_VS_Edep_single_epFDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_epFDn);
    TH2D *h_E_miss_VS_Edep_single_epFDn = new TH2D("E_miss_VS_Edep_single_epFDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_epFDn);
    TH2D *h_M_miss_VS_Edep_single_epFDn = new TH2D("M_miss_VS_Edep_single_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_epFDn);
    TH2D *h_path_VS_Edep_single_epFDn = new TH2D("path_VS_Edep_single_epFDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_single_epFDn = new TH2D("theta_n_miss_VS_Edep_single_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_epFDn);
    TH2D *h_ToF_VS_Edep_single_epFDn = new TH2D("ToF_VS_Edep_single_epFDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_epFDn);
    TH2D *h_nSector_VS_Edep_single_epFDn = new TH2D("nSector_VS_Edep_single_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_epFDn);

    TH1D *h_Edep_CND1_epCDn = new TH1D("Edep_CND1_epCDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_epCDn);
    TH2D *h_P_n_VS_Edep_CND1_epCDn = new TH2D("P_n_VS_Edep_CND1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND1_epCDn = new TH2D("theta_n_VS_Edep_CND1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND1_epCDn = new TH2D("phi_n_VS_Edep_CND1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND1_epCDn = new TH2D("P_miss_VS_Edep_CND1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND1_epCDn = new TH2D("theta_miss_VS_Edep_CND1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND1_epCDn = new TH2D("phi_miss_VS_Edep_CND1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_epCDn);
    TH2D *h_dpp_VS_Edep_CND1_epCDn = new TH2D("dpp_VS_Edep_CND1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND1_epCDn = new TH2D("beta_n_VS_Edep_CND1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND1_epCDn = new TH2D("E_miss_VS_Edep_CND1_epCDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND1_epCDn = new TH2D("M_miss_VS_Edep_CND1_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_epCDn);
    TH2D *h_path_VS_Edep_CND1_epCDn = new TH2D("path_VS_Edep_CND1_epCDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_epCDn);
    TH2D *h_ToF_VS_Edep_CND1_epCDn = new TH2D("ToF_VS_Edep_CND1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_epCDn);
    TH2D *h_nSector_VS_Edep_CND1_epCDn = new TH2D("nSector_VS_Edep_CND1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_epCDn = new TH2D("Edep_CND2_VS_Edep_CND1_epCDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_epCDn = new TH2D("Edep_CND3_VS_Edep_CND1_epCDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_epCDn);

    TH1D *h_Edep_CND1_epFDn = new TH1D("Edep_CND1_epFDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_epFDn);
    TH2D *h_P_n_VS_Edep_CND1_epFDn = new TH2D("P_n_VS_Edep_CND1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND1_epFDn = new TH2D("theta_n_VS_Edep_CND1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND1_epFDn = new TH2D("phi_n_VS_Edep_CND1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND1_epFDn = new TH2D("P_miss_VS_Edep_CND1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND1_epFDn = new TH2D("theta_miss_VS_Edep_CND1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND1_epFDn = new TH2D("phi_miss_VS_Edep_CND1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_epFDn);
    TH2D *h_dpp_VS_Edep_CND1_epFDn = new TH2D("dpp_VS_Edep_CND1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND1_epFDn = new TH2D("beta_n_VS_Edep_CND1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND1_epFDn = new TH2D("E_miss_VS_Edep_CND1_epFDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND1_epFDn = new TH2D("M_miss_VS_Edep_CND1_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_epFDn);
    TH2D *h_path_VS_Edep_CND1_epFDn = new TH2D("path_VS_Edep_CND1_epFDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_epFDn);
    TH2D *h_ToF_VS_Edep_CND1_epFDn = new TH2D("ToF_VS_Edep_CND1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_epFDn);
    TH2D *h_nSector_VS_Edep_CND1_epFDn = new TH2D("nSector_VS_Edep_CND1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_epFDn = new TH2D("Edep_CND2_VS_Edep_CND1_epFDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_epFDn = new TH2D("Edep_CND3_VS_Edep_CND1_epFDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_epFDn);

    TH1D *h_Edep_CND2_epCDn = new TH1D("Edep_CND2_epCDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_epCDn);
    TH2D *h_P_n_VS_Edep_CND2_epCDn = new TH2D("P_n_VS_Edep_CND2_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_epCDn);
    TH2D *h_theta_n_VS_Edep_CND2_epCDn = new TH2D("theta_n_VS_Edep_CND2_epCDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_epCDn);
    TH2D *h_phi_n_VS_Edep_CND2_epCDn = new TH2D("phi_n_VS_Edep_CND2_epCDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_epCDn);
    TH2D *h_P_miss_VS_Edep_CND2_epCDn = new TH2D("P_miss_VS_Edep_CND2_epCDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND2_epCDn = new TH2D("theta_miss_VS_Edep_CND2_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND2_epCDn = new TH2D("phi_miss_VS_Edep_CND2_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_epCDn);
    TH2D *h_dpp_VS_Edep_CND2_epCDn = new TH2D("dpp_VS_Edep_CND2_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_epCDn);
    TH2D *h_beta_n_VS_Edep_CND2_epCDn = new TH2D("beta_n_VS_Edep_CND2_epCDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_epCDn);
    TH2D *h_E_miss_VS_Edep_CND2_epCDn = new TH2D("E_miss_VS_Edep_CND2_epCDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_epCDn);
    TH2D *h_M_miss_VS_Edep_CND2_epCDn = new TH2D("M_miss_VS_Edep_CND2_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_epCDn);
    TH2D *h_path_VS_Edep_CND2_epCDn = new TH2D("path_VS_Edep_CND2_epCDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_epCDn = new TH2D("theta_n_miss_VS_Edep_CND2_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_epCDn);
    TH2D *h_ToF_VS_Edep_CND2_epCDn = new TH2D("ToF_VS_Edep_CND2_epCDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_epCDn);
    TH2D *h_nSector_VS_Edep_CND2_epCDn = new TH2D("nSector_VS_Edep_CND2_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_epCDn = new TH2D("Edep_CND3_VS_Edep_CND2_epCDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_epCDn);

    TH1D *h_Edep_CND2_epFDn = new TH1D("Edep_CND2_epFDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_epFDn);
    TH2D *h_P_n_VS_Edep_CND2_epFDn = new TH2D("P_n_VS_Edep_CND2_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_epFDn);
    TH2D *h_theta_n_VS_Edep_CND2_epFDn = new TH2D("theta_n_VS_Edep_CND2_epFDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_epFDn);
    TH2D *h_phi_n_VS_Edep_CND2_epFDn = new TH2D("phi_n_VS_Edep_CND2_epFDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_epFDn);
    TH2D *h_P_miss_VS_Edep_CND2_epFDn = new TH2D("P_miss_VS_Edep_CND2_epFDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND2_epFDn = new TH2D("theta_miss_VS_Edep_CND2_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND2_epFDn = new TH2D("phi_miss_VS_Edep_CND2_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_epFDn);
    TH2D *h_dpp_VS_Edep_CND2_epFDn = new TH2D("dpp_VS_Edep_CND2_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_epFDn);
    TH2D *h_beta_n_VS_Edep_CND2_epFDn = new TH2D("beta_n_VS_Edep_CND2_epFDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_epFDn);
    TH2D *h_E_miss_VS_Edep_CND2_epFDn = new TH2D("E_miss_VS_Edep_CND2_epFDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_epFDn);
    TH2D *h_M_miss_VS_Edep_CND2_epFDn = new TH2D("M_miss_VS_Edep_CND2_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_epFDn);
    TH2D *h_path_VS_Edep_CND2_epFDn = new TH2D("path_VS_Edep_CND2_epFDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_epFDn = new TH2D("theta_n_miss_VS_Edep_CND2_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_epFDn);
    TH2D *h_ToF_VS_Edep_CND2_epFDn = new TH2D("ToF_VS_Edep_CND2_epFDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_epFDn);
    TH2D *h_nSector_VS_Edep_CND2_epFDn = new TH2D("nSector_VS_Edep_CND2_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_epFDn = new TH2D("Edep_CND3_VS_Edep_CND2_epFDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_epFDn);

    TH1D *h_Edep_CND3_epCDn = new TH1D("Edep_CND3_epCDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_epCDn);
    TH2D *h_P_n_VS_Edep_CND3_epCDn = new TH2D("P_n_VS_Edep_CND3_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_epCDn);
    TH2D *h_theta_n_VS_Edep_CND3_epCDn = new TH2D("theta_n_VS_Edep_CND3_epCDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_epCDn);
    TH2D *h_phi_n_VS_Edep_CND3_epCDn = new TH2D("phi_n_VS_Edep_CND3_epCDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_epCDn);
    TH2D *h_P_miss_VS_Edep_CND3_epCDn = new TH2D("P_miss_VS_Edep_CND3_epCDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND3_epCDn = new TH2D("theta_miss_VS_Edep_CND3_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND3_epCDn = new TH2D("phi_miss_VS_Edep_CND3_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_epCDn);
    TH2D *h_dpp_VS_Edep_CND3_epCDn = new TH2D("dpp_VS_Edep_CND3_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_epCDn);
    TH2D *h_beta_n_VS_Edep_CND3_epCDn = new TH2D("beta_n_VS_Edep_CND3_epCDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_epCDn);
    TH2D *h_E_miss_VS_Edep_CND3_epCDn = new TH2D("E_miss_VS_Edep_CND3_epCDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_epCDn);
    TH2D *h_M_miss_VS_Edep_CND3_epCDn = new TH2D("M_miss_VS_Edep_CND3_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_epCDn);
    TH2D *h_path_VS_Edep_CND3_epCDn = new TH2D("path_VS_Edep_CND3_epCDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_epCDn = new TH2D("theta_n_miss_VS_Edep_CND3_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_epCDn);
    TH2D *h_ToF_VS_Edep_CND3_epCDn = new TH2D("ToF_VS_Edep_CND3_epCDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_epCDn);
    TH2D *h_nSector_VS_Edep_CND3_epCDn = new TH2D("nSector_VS_Edep_CND3_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_epCDn);

    TH1D *h_Edep_CND3_epFDn = new TH1D("Edep_CND3_epFDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_epFDn);
    TH2D *h_P_n_VS_Edep_CND3_epFDn = new TH2D("P_n_VS_Edep_CND3_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_epFDn);
    TH2D *h_theta_n_VS_Edep_CND3_epFDn = new TH2D("theta_n_VS_Edep_CND3_epFDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_epFDn);
    TH2D *h_phi_n_VS_Edep_CND3_epFDn = new TH2D("phi_n_VS_Edep_CND3_epFDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_epFDn);
    TH2D *h_P_miss_VS_Edep_CND3_epFDn = new TH2D("P_miss_VS_Edep_CND3_epFDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND3_epFDn = new TH2D("theta_miss_VS_Edep_CND3_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND3_epFDn = new TH2D("phi_miss_VS_Edep_CND3_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_epFDn);
    TH2D *h_dpp_VS_Edep_CND3_epFDn = new TH2D("dpp_VS_Edep_CND3_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_epFDn);
    TH2D *h_beta_n_VS_Edep_CND3_epFDn = new TH2D("beta_n_VS_Edep_CND3_epFDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_epFDn);
    TH2D *h_E_miss_VS_Edep_CND3_epFDn = new TH2D("E_miss_VS_Edep_CND3_epFDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_epFDn);
    TH2D *h_M_miss_VS_Edep_CND3_epFDn = new TH2D("M_miss_VS_Edep_CND3_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_epFDn);
    TH2D *h_path_VS_Edep_CND3_epFDn = new TH2D("path_VS_Edep_CND3_epFDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_epFDn = new TH2D("theta_n_miss_VS_Edep_CND3_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_epFDn);
    TH2D *h_ToF_VS_Edep_CND3_epFDn = new TH2D("ToF_VS_Edep_CND3_epFDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_epFDn);
    TH2D *h_nSector_VS_Edep_CND3_epFDn = new TH2D("nSector_VS_Edep_CND3_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_epFDn);

    TH2D *h_Size_CND1_VS_Size_CND2_epCDn = new TH2D("Size_CND1_VS_Size_CND2_epCDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_epCDn);
    TH2D *h_Size_CND1_VS_Size_CND3_epCDn = new TH2D("Size_CND1_VS_Size_CND3_epCDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_epCDn);
    TH2D *h_Size_CND2_VS_Size_CND3_epCDn = new TH2D("Size_CND2_VS_Size_CND3_epCDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_epCDn);

    TH2D *h_Size_CND1_VS_Size_CND2_epFDn = new TH2D("Size_CND1_VS_Size_CND2_epFDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_epFDn);
    TH2D *h_Size_CND1_VS_Size_CND3_epFDn = new TH2D("Size_CND1_VS_Size_CND3_epFDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_epFDn);
    TH2D *h_Size_CND2_VS_Size_CND3_epFDn = new TH2D("Size_CND2_VS_Size_CND3_epFDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_epFDn);

    TH1D *h_ToF_epCDn = new TH1D("ToF_epCDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_epCDn);
    TH2D *h_P_n_VS_ToF_epCDn = new TH2D("P_n_VS_ToF_epCDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_epCDn);
    TH2D *h_theta_n_VS_ToF_epCDn = new TH2D("theta_n_VS_ToF_epCDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_epCDn);
    TH2D *h_phi_n_VS_ToF_epCDn = new TH2D("phi_n_VS_ToF_epCDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_epCDn);
    TH2D *h_P_miss_VS_ToF_epCDn = new TH2D("P_miss_VS_ToF_epCDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_epCDn);
    TH2D *h_theta_miss_VS_ToF_epCDn = new TH2D("theta_miss_VS_ToF_epCDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_epCDn);
    TH2D *h_phi_miss_VS_ToF_epCDn = new TH2D("phi_miss_VS_ToF_epCDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_epCDn);
    TH2D *h_dpp_VS_ToF_epCDn = new TH2D("dpp_VS_ToF_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_epCDn);
    TH2D *h_beta_n_VS_ToF_epCDn = new TH2D("beta_n_VS_ToF_epCDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_epCDn);
    TH2D *h_E_miss_VS_ToF_epCDn = new TH2D("E_miss_VS_ToF_epCDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_epCDn);
    TH2D *h_M_miss_VS_ToF_epCDn = new TH2D("M_miss_VS_ToF_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_epCDn);
    TH2D *h_path_VS_ToF_epCDn = new TH2D("path_VS_ToF_epCDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_epCDn);
    TH2D *h_theta_n_miss_VS_ToF_epCDn = new TH2D("theta_n_miss_VS_ToF_epCDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_epCDn);
    TH2D *h_nSector_VS_ToF_epCDn = new TH2D("nSector_VS_ToF_epCDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_epCDn);

    TH1D *h_ToF_epFDn = new TH1D("ToF_epFDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_epFDn);
    TH2D *h_P_n_VS_ToF_epFDn = new TH2D("P_n_VS_ToF_epFDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_epFDn);
    TH2D *h_theta_n_VS_ToF_epFDn = new TH2D("theta_n_VS_ToF_epFDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_epFDn);
    TH2D *h_phi_n_VS_ToF_epFDn = new TH2D("phi_n_VS_ToF_epFDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_epFDn);
    TH2D *h_P_miss_VS_ToF_epFDn = new TH2D("P_miss_VS_ToF_epFDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_epFDn);
    TH2D *h_theta_miss_VS_ToF_epFDn = new TH2D("theta_miss_VS_ToF_epFDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_epFDn);
    TH2D *h_phi_miss_VS_ToF_epFDn = new TH2D("phi_miss_VS_ToF_epFDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_epFDn);
    TH2D *h_dpp_VS_ToF_epFDn = new TH2D("dpp_VS_ToF_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_epFDn);
    TH2D *h_beta_n_VS_ToF_epFDn = new TH2D("beta_n_VS_ToF_epFDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_epFDn);
    TH2D *h_E_miss_VS_ToF_epFDn = new TH2D("E_miss_VS_ToF_epFDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_epFDn);
    TH2D *h_M_miss_VS_ToF_epFDn = new TH2D("M_miss_VS_ToF_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_epFDn);
    TH2D *h_path_VS_ToF_epFDn = new TH2D("path_VS_ToF_epFDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_epFDn);
    TH2D *h_theta_n_miss_VS_ToF_epFDn = new TH2D("theta_n_miss_VS_ToF_epFDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_epFDn);
    TH2D *h_nSector_VS_ToF_epFDn = new TH2D("nSector_VS_ToF_epFDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_epFDn);

    // Step Zero (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step0 = new TH1D("n_multiplicity_allN_epCDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epCDn_Step0);
    TH1D *h_n_multiplicity_goodN_epCDn_Step0 = new TH1D("n_multiplicity_goodN_epCDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epCDn_Step0);
    TH1D *h_n_multiplicity_badN_epCDn_Step0 = new TH1D("n_multiplicity_badN_epCDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epCDn_Step0);

    TH1D *h_n_multiplicity_allN_epFDn_Step0 = new TH1D("n_multiplicity_allN_epFDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epFDn_Step0);
    TH1D *h_n_multiplicity_goodN_epFDn_Step0 = new TH1D("n_multiplicity_goodN_epFDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epFDn_Step0);
    TH1D *h_n_multiplicity_badN_epFDn_Step0 = new TH1D("n_multiplicity_badN_epFDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epFDn_Step0);

    /* Step0 cuts */
    TH2D *h_dbeta_n_VS_P_n_BS0C_Step0_epCDn = new TH2D("dbeta_n_VS_P_n_BS0C_Step0_epCDn", "#Delta#beta_{n} vs Neutron Momentum (Before Step0 Cuts);P_{n} [GeV/c];#Delta#beta_{n}", 50, 0, 1.5, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_P_n_BS0C_Step0_epCDn);
    TH2D *h_dbeta_n_VS_ToF_BS0C_Step0_epCDn = new TH2D("dbeta_n_VS_ToF_BS0C_Step0_epCDn", "#Delta#beta_{n} vs Neutron ToF (Before Step0 Cuts);ToF [ns];#Delta#beta_{n}", 50, 0, 50, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_ToF_BS0C_Step0_epCDn);
    TH2D *h_dbeta_n_VS_P_n_AS0C_Step0_epCDn = new TH2D("dbeta_n_VS_P_n_AS0C_Step0_epCDn", "#Delta#beta_{n} vs Neutron Momentum (After Step0 Cuts);P_{n} [GeV/c];#Delta#beta_{n}", 50, 0, 1.5, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_P_n_AS0C_Step0_epCDn);
    TH2D *h_dbeta_n_VS_ToF_AS0C_Step0_epCDn = new TH2D("dbeta_n_VS_ToF_AS0C_Step0_epCDn", "#Delta#beta_{n} vs Neutron ToF (After Step0 Cuts);ToF [ns];#Delta#beta_{n}", 50, 0, 50, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_ToF_AS0C_Step0_epCDn);

    TH2D *h_dbeta_n_VS_P_n_BS0C_Step0_epFDn = new TH2D("dbeta_n_VS_P_n_BS0C_Step0_epFDn", "#Delta#beta_{n} vs Neutron Momentum (Before Step0 Cuts);P_{n} [GeV/c];#Delta#beta_{n}", 50, 0, 1.5, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_P_n_BS0C_Step0_epFDn);
    TH2D *h_dbeta_n_VS_ToF_BS0C_Step0_epFDn = new TH2D("dbeta_n_VS_ToF_BS0C_Step0_epFDn", "#Delta#beta_{n} vs Neutron ToF (Before Step0 Cuts);ToF [ns];#Delta#beta_{n}", 50, 0, 50, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_ToF_BS0C_Step0_epFDn);
    TH2D *h_dbeta_n_VS_P_n_AS0C_Step0_epFDn = new TH2D("dbeta_n_VS_P_n_AS0C_Step0_epFDn", "#Delta#beta_{n} vs Neutron Momentum (After Step0 Cuts);P_{n} [GeV/c];#Delta#beta_{n}", 50, 0, 1.5, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_P_n_AS0C_Step0_epFDn);
    TH2D *h_dbeta_n_VS_ToF_AS0C_Step0_epFDn = new TH2D("dbeta_n_VS_ToF_AS0C_Step0_epFDn", "#Delta#beta_{n} vs Neutron ToF (After Step0 Cuts);ToF [ns];#Delta#beta_{n}", 50, 0, 50, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_ToF_AS0C_Step0_epFDn);

    TH1D *h_Vhit_z_n_BS0C_Step0_epCDn = new TH1D("Vhit_z_n_BS0C_Step0_epCDn", "V_{hit,z}^{n} Distribution (Before Step0 Cuts);V_{hit,z}^{n} [cm]", 50, -50, 50);
    HistoList.push_back(h_Vhit_z_n_BS0C_Step0_epCDn);
    TH1D *h_Vhit_z_n_AS0C_Step0_epCDn = new TH1D("Vhit_z_n_AS0C_Step0_epCDn", "V_{hit,z}^{n} Distribution (After Step0 Cuts);V_{hit,z}^{n} [cm]", 50, -50, 50);
    HistoList.push_back(h_Vhit_z_n_AS0C_Step0_epCDn);

    TH1D *h_Vhit_z_n_BS0C_Step0_epFDn = new TH1D("Vhit_z_n_BS0C_Step0_epFDn", "V_{hit,z}^{n} Distribution (Before Step0 Cuts);V_{hit,z}^{n} [cm]", 50, -50, 50);
    HistoList.push_back(h_Vhit_z_n_BS0C_Step0_epFDn);
    TH1D *h_Vhit_z_n_AS0C_Step0_epFDn = new TH1D("Vhit_z_n_AS0C_Step0_epFDn", "V_{hit,z}^{n} Distribution (After Step0 Cuts);V_{hit,z}^{n} [cm]", 50, -50, 50);
    HistoList.push_back(h_Vhit_z_n_AS0C_Step0_epFDn);

    TH1D *h_ToF_n_BS0C_Step0_epCDn = new TH1D("ToF_n_BS0C_Step0_epCDn", "ToF Distribution (Before Step0 Cuts);ToF [ns]", 50, 0, 50);
    HistoList.push_back(h_ToF_n_BS0C_Step0_epCDn);
    TH1D *h_ToF_n_AS0C_Step0_epCDn = new TH1D("ToF_n_AS0C_Step0_epCDn", "ToF Distribution (After Step0 Cuts);ToF [ns]", 50, 0, 50);
    HistoList.push_back(h_ToF_n_AS0C_Step0_epCDn);

    TH1D *h_ToF_n_BS0C_Step0_epFDn = new TH1D("ToF_n_BS0C_Step0_epFDn", "ToF Distribution (Before Step0 Cuts);ToF [ns]", 50, 0, 50);
    HistoList.push_back(h_ToF_n_BS0C_Step0_epFDn);
    TH1D *h_ToF_n_AS0C_Step0_epFDn = new TH1D("ToF_n_AS0C_Step0_epFDn", "ToF Distribution (After Step0 Cuts);ToF [ns]", 50, 0, 50);
    HistoList.push_back(h_ToF_n_AS0C_Step0_epFDn);

    /* Kinematical variables */
    TH1D *h_theta_n_goodN_Step0_epCDn = new TH1D("theta_n_goodN_Step0_epCDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_goodN_Step0_epCDn);
    TH1D *h_theta_n_badN_Step0_epCDn = new TH1D("theta_n_badN_Step0_epCDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_badN_Step0_epCDn);
    TH1D *h_phi_n_goodN_Step0_epCDn = new TH1D("phi_n_goodN_Step0_epCDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_goodN_Step0_epCDn);
    TH1D *h_phi_n_badN_Step0_epCDn = new TH1D("phi_n_badN_Step0_epCDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_phi_n_goodN_Step0_epCDn = new TH2D("theta_n_VS_phi_n_goodN_Step0_epCDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_phi_n_badN_Step0_epCDn = new TH2D("theta_n_VS_phi_n_badN_Step0_epCDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_beta_n_goodN_Step0_epCDn = new TH2D("theta_VS_beta_goodN_Step0_epCDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_beta_n_badN_Step0_epCDn = new TH2D("theta_VS_beta_badN_Step0_epCDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_badN_Step0_epCDn);

    TH1D *h_theta_n_goodN_Step0_epFDn = new TH1D("theta_n_goodN_Step0_epFDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_goodN_Step0_epFDn);
    TH1D *h_theta_n_badN_Step0_epFDn = new TH1D("theta_n_badN_Step0_epFDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_badN_Step0_epFDn);
    TH1D *h_phi_n_goodN_Step0_epFDn = new TH1D("phi_n_goodN_Step0_epFDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_goodN_Step0_epFDn);
    TH1D *h_phi_n_badN_Step0_epFDn = new TH1D("phi_n_badN_Step0_epFDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_phi_n_goodN_Step0_epFDn = new TH2D("theta_n_VS_phi_n_goodN_Step0_epFDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_phi_n_badN_Step0_epFDn = new TH2D("theta_n_VS_phi_n_badN_Step0_epFDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_beta_n_goodN_Step0_epFDn = new TH2D("theta_VS_beta_goodN_Step0_epFDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_beta_n_badN_Step0_epFDn = new TH2D("theta_VS_beta_badN_Step0_epFDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_badN_Step0_epFDn);

    TH1D *h_P_n_goodN_Step0_epCDn = new TH1D("P_n_goodN_Step0_epCDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_goodN_Step0_epCDn);
    TH1D *h_P_n_badN_Step0_epCDn = new TH1D("P_n_badN_Step0_epCDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_badN_Step0_epCDn);
    TH2D *h_P_n_VS_theta_n_goodN_Step0_epCDn = new TH2D("P_n_VS_theta_n_goodN_Step0_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_theta_n_badN_Step0_epCDn = new TH2D("P_n_VS_theta_n_badN_Step0_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_badN_Step0_epCDn);

    TH1D *h_P_n_goodN_Step0_epFDn = new TH1D("P_n_goodN_Step0_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_goodN_Step0_epFDn);
    TH1D *h_P_n_badN_Step0_epFDn = new TH1D("P_n_badN_Step0_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_badN_Step0_epFDn);
    TH2D *h_P_n_VS_theta_n_goodN_Step0_epFDn = new TH2D("P_n_VS_theta_n_goodN_Step0_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_theta_n_badN_Step0_epFDn = new TH2D("P_n_VS_theta_n_badN_Step0_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_badN_Step0_epFDn);

    TH1D *h_P_miss_goodN_Step0_epCDn = new TH1D("P_miss_goodN_Step0_epCDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_goodN_Step0_epCDn);
    TH1D *h_P_miss_badN_Step0_epCDn = new TH1D("P_miss_badN_Step0_epCDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_theta_miss_goodN_Step0_epCDn = new TH2D("P_miss_VS_theta_miss_goodN_Step0_epCDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_theta_miss_badN_Step0_epCDn = new TH2D("P_miss_VS_theta_miss_badN_Step0_epCDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_phi_miss_goodN_Step0_epCDn = new TH2D("P_miss_VS_phi_miss_goodN_Step0_epCDn", "Missing Momentum vs #phi_{miss};#phi_{miss} [#circ];P_{miss} [GeV/c]", 48, -180, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_phi_miss_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_phi_miss_badN_Step0_epCDn = new TH2D("P_miss_VS_phi_miss_badN_Step0_epCDn", "Missing Momentum vs #phi_{miss};#phi_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, -180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_phi_miss_badN_Step0_epCDn);

    TH1D *h_P_miss_goodN_Step0_epFDn = new TH1D("P_miss_goodN_Step0_epFDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_goodN_Step0_epFDn);
    TH1D *h_P_miss_badN_Step0_epFDn = new TH1D("P_miss_badN_Step0_epFDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_theta_miss_goodN_Step0_epFDn = new TH2D("P_miss_VS_theta_miss_goodN_Step0_epFDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_theta_miss_badN_Step0_epFDn = new TH2D("P_miss_VS_theta_miss_badN_Step0_epFDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_phi_miss_goodN_Step0_epFDn = new TH2D("P_miss_VS_phi_miss_goodN_Step0_epFDn", "Missing Momentum vs #phi_{miss};#phi_{miss} [#circ];P_{miss} [GeV/c]", 48, -180, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_phi_miss_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_phi_miss_badN_Step0_epFDn = new TH2D("P_miss_VS_phi_miss_badN_Step0_epFDn", "Missing Momentum vs #phi_{miss};#phi_{miss} [#circ];P_{miss} [GeV/c]", 48, -180, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_phi_miss_badN_Step0_epFDn);

    TH1D *h_dpp_allN_Step0_epCDn = new TH1D("dpp_allN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_allN_Step0_epCDn);
    TH1D *h_dpp_goodN_Step0_epCDn = new TH1D("dpp_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_goodN_Step0_epCDn);
    TH1D *h_dpp_badN_Step0_epCDn = new TH1D("dpp_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_badN_Step0_epCDn);

    TH1D *h_dpp_allN_Step0_epFDn = new TH1D("dpp_allN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_allN_Step0_epFDn);
    TH1D *h_dpp_goodN_Step0_epFDn = new TH1D("dpp_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_goodN_Step0_epFDn);
    TH1D *h_dpp_badN_Step0_epFDn = new TH1D("dpp_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_badN_Step0_epFDn);

    TH1D *h_theta_n_miss_allN_Step0_epCDn = new TH1D("theta_n_miss_allN_Step0_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_allN_Step0_epCDn);
    TH1D *h_theta_n_miss_goodN_Step0_epCDn = new TH1D("theta_n_miss_goodN_Step0_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_goodN_Step0_epCDn);
    TH1D *h_theta_n_miss_badN_Step0_epCDn = new TH1D("theta_n_miss_badN_Step0_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_badN_Step0_epCDn);

    TH1D *h_theta_n_miss_allN_Step0_epFDn = new TH1D("theta_n_miss_allN_Step0_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_allN_Step0_epFDn);
    TH1D *h_theta_n_miss_goodN_Step0_epFDn = new TH1D("theta_n_miss_goodN_Step0_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_goodN_Step0_epFDn);
    TH1D *h_theta_n_miss_badN_Step0_epFDn = new TH1D("theta_n_miss_badN_Step0_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_badN_Step0_epFDn);

    TH2D *h_dpp_VS_theta_n_miss_allN_Step0_epCDn = new TH2D("dpp_VS_theta_n_miss_allN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -1.5, 1.5, 50, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_allN_Step0_epCDn);

    TH2D *h_dpp_VS_theta_n_miss_allN_Step0_epFDn = new TH2D("dpp_VS_theta_n_miss_allN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -1.5, 1.5, 50, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_allN_Step0_epFDn);

    TH1D *h_E_p_goodN_Step0_epCDn = new TH1D("E_p_goodN_Step0_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_goodN_Step0_epCDn);
    TH1D *h_E_p_badN_Step0_epCDn = new TH1D("E_p_badN_Step0_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_badN_Step0_epCDn);
    TH1D *h_E_miss_goodN_Step0_epCDn = new TH1D("E_miss_goodN_Step0_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_goodN_Step0_epCDn);
    TH1D *h_E_miss_badN_Step0_epCDn = new TH1D("E_miss_badN_Step0_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_badN_Step0_epCDn);
    TH1D *h_M_miss_goodN_Step0_epCDn = new TH1D("M_miss_goodN_Step0_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_goodN_Step0_epCDn);
    TH1D *h_M_miss_badN_Step0_epCDn = new TH1D("M_miss_badN_Step0_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_P_n_goodN_Step0_epCDn = new TH2D("M_miss_VS_P_n_goodN_Step0_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_P_n_badN_Step0_epCDn = new TH2D("M_miss_VS_P_n_badN_Step0_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_theta_n_goodN_Step0_epCDn = new TH2D("M_miss_VS_theta_n_goodN_Step0_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_theta_n_badN_Step0_epCDn = new TH2D("M_miss_VS_theta_n_badN_Step0_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_phi_n_goodN_Step0_epCDn = new TH2D("M_miss_VS_phi_n_goodN_Step0_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_phi_n_badN_Step0_epCDn = new TH2D("M_miss_VS_phi_n_badN_Step0_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_P_miss_goodN_Step0_epCDn = new TH2D("M_miss_VS_P_miss_goodN_Step0_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_P_miss_badN_Step0_epCDn = new TH2D("M_miss_VS_P_miss_badN_Step0_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_theta_miss_goodN_Step0_epCDn = new TH2D("M_miss_VS_theta_miss_goodN_Step0_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_theta_miss_badN_Step0_epCDn = new TH2D("M_miss_VS_theta_miss_badN_Step0_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_phi_miss_goodN_Step0_epCDn = new TH2D("M_miss_VS_phi_miss_goodN_Step0_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_phi_miss_badN_Step0_epCDn = new TH2D("M_miss_VS_phi_miss_badN_Step0_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_badN_Step0_epCDn);

    TH1D *h_E_p_goodN_Step0_epFDn = new TH1D("E_p_goodN_Step0_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_goodN_Step0_epFDn);
    TH1D *h_E_p_badN_Step0_epFDn = new TH1D("E_p_badN_Step0_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_badN_Step0_epFDn);
    TH1D *h_E_miss_goodN_Step0_epFDn = new TH1D("E_miss_goodN_Step0_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_goodN_Step0_epFDn);
    TH1D *h_E_miss_badN_Step0_epFDn = new TH1D("E_miss_badN_Step0_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_badN_Step0_epFDn);
    TH1D *h_M_miss_goodN_Step0_epFDn = new TH1D("M_miss_goodN_Step0_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_goodN_Step0_epFDn);
    TH1D *h_M_miss_badN_Step0_epFDn = new TH1D("M_miss_badN_Step0_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_P_n_goodN_Step0_epFDn = new TH2D("M_miss_VS_P_n_goodN_Step0_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_P_n_badN_Step0_epFDn = new TH2D("M_miss_VS_P_n_badN_Step0_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_theta_n_goodN_Step0_epFDn = new TH2D("M_miss_VS_theta_n_goodN_Step0_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_theta_n_badN_Step0_epFDn = new TH2D("M_miss_VS_theta_n_badN_Step0_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_phi_n_goodN_Step0_epFDn = new TH2D("M_miss_VS_phi_n_goodN_Step0_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_phi_n_badN_Step0_epFDn = new TH2D("M_miss_VS_phi_n_badN_Step0_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_P_miss_goodN_Step0_epFDn = new TH2D("M_miss_VS_P_miss_goodN_Step0_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_P_miss_badN_Step0_epFDn = new TH2D("M_miss_VS_P_miss_badN_Step0_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_theta_miss_goodN_Step0_epFDn = new TH2D("M_miss_VS_theta_miss_goodN_Step0_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_theta_miss_badN_Step0_epFDn = new TH2D("M_miss_VS_theta_miss_badN_Step0_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_phi_miss_goodN_Step0_epFDn = new TH2D("M_miss_VS_phi_miss_goodN_Step0_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_phi_miss_badN_Step0_epFDn = new TH2D("M_miss_VS_phi_miss_badN_Step0_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_badN_Step0_epFDn);

    TH1D *h_P_n_minus_P_miss_goodN_Step0_epCDn = new TH1D("P_n_minus_P_miss_goodN_Step0_epCDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_goodN_Step0_epCDn);
    TH1D *h_P_n_minus_P_miss_badN_Step0_epCDn = new TH1D("P_n_minus_P_miss_badN_Step0_epCDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_badN_Step0_epCDn);
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step0_epCDn = new TH1D("P_n_x_minus_P_miss_x_goodN_Step0_epCDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_goodN_Step0_epCDn);
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step0_epCDn = new TH1D("P_n_x_minus_P_miss_x_badN_Step0_epCDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_badN_Step0_epCDn);
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step0_epCDn = new TH1D("P_n_y_minus_P_miss_y_goodN_Step0_epCDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_goodN_Step0_epCDn);
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step0_epCDn = new TH1D("P_n_y_minus_P_miss_y_badN_Step0_epCDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_badN_Step0_epCDn);
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step0_epCDn = new TH1D("P_n_z_minus_P_miss_z_goodN_Step0_epCDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_goodN_Step0_epCDn);
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step0_epCDn = new TH1D("P_n_z_minus_P_miss_z_badN_Step0_epCDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_badN_Step0_epCDn);

    TH1D *h_P_n_minus_P_miss_goodN_Step0_epFDn = new TH1D("P_n_minus_P_miss_goodN_Step0_epFDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_goodN_Step0_epFDn);
    TH1D *h_P_n_minus_P_miss_badN_Step0_epFDn = new TH1D("P_n_minus_P_miss_badN_Step0_epFDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_badN_Step0_epFDn);
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step0_epFDn = new TH1D("P_n_x_minus_P_miss_x_goodN_Step0_epFDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_goodN_Step0_epFDn);
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step0_epFDn = new TH1D("P_n_x_minus_P_miss_x_badN_Step0_epFDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_badN_Step0_epFDn);
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step0_epFDn = new TH1D("P_n_y_minus_P_miss_y_goodN_Step0_epFDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_goodN_Step0_epFDn);
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step0_epFDn = new TH1D("P_n_y_minus_P_miss_y_badN_Step0_epFDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_badN_Step0_epFDn);
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step0_epFDn = new TH1D("P_n_z_minus_P_miss_z_goodN_Step0_epFDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_goodN_Step0_epFDn);
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step0_epFDn = new TH1D("P_n_z_minus_P_miss_z_badN_Step0_epFDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_badN_Step0_epFDn);

    TH2D *h_P_n_VS_P_miss_goodN_Step0_epCDn = new TH2D("P_n_VS_P_miss_goodN_Step0_epCDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_P_miss_badN_Step0_epCDn = new TH2D("P_n_VS_P_miss_badN_Step0_epCDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_badN_Step0_epCDn);
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step0_epCDn = new TH2D("P_n_x_VS_P_miss_x_goodN_Step0_epCDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_goodN_Step0_epCDn);
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step0_epCDn = new TH2D("P_n_x_VS_P_miss_x_badN_Step0_epCDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_badN_Step0_epCDn);
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step0_epCDn = new TH2D("P_n_y_VS_P_miss_y_goodN_Step0_epCDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_goodN_Step0_epCDn);
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step0_epCDn = new TH2D("P_n_y_VS_P_miss_y_badN_Step0_epCDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_badN_Step0_epCDn);
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step0_epCDn = new TH2D("P_n_z_VS_P_miss_z_goodN_Step0_epCDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_goodN_Step0_epCDn);
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step0_epCDn = new TH2D("P_n_z_VS_P_miss_z_badN_Step0_epCDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_badN_Step0_epCDn);

    TH2D *h_P_n_VS_P_miss_goodN_Step0_epFDn = new TH2D("P_n_VS_P_miss_goodN_Step0_epFDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_P_miss_badN_Step0_epFDn = new TH2D("P_n_VS_P_miss_badN_Step0_epFDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_badN_Step0_epFDn);
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step0_epFDn = new TH2D("P_n_x_VS_P_miss_x_goodN_Step0_epFDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_goodN_Step0_epFDn);
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step0_epFDn = new TH2D("P_n_x_VS_P_miss_x_badN_Step0_epFDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_badN_Step0_epFDn);
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step0_epFDn = new TH2D("P_n_y_VS_P_miss_y_goodN_Step0_epFDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_goodN_Step0_epFDn);
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step0_epFDn = new TH2D("P_n_y_VS_P_miss_y_badN_Step0_epFDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_badN_Step0_epFDn);
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step0_epFDn = new TH2D("P_n_z_VS_P_miss_z_goodN_Step0_epFDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_goodN_Step0_epFDn);
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step0_epFDn = new TH2D("P_n_z_VS_P_miss_z_badN_Step0_epFDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_badN_Step0_epFDn);

    TH1D *h_theta_n_p_goodN_Step0_epCDn = new TH1D("theta_n_p_goodN_Step0_epCDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_goodN_Step0_epCDn);
    TH1D *h_theta_n_p_badN_Step0_epCDn = new TH1D("theta_n_p_badN_Step0_epCDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_badN_Step0_epCDn);
    TH2D *h_theta_n_p_VS_P_p_goodN_Step0_epCDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_goodN_Step0_epCDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_goodN_Step0_epCDn);
    TH2D *h_theta_n_p_VS_P_p_badN_Step0_epCDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_badN_Step0_epCDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_badN_Step0_epCDn);

    TH1D *h_theta_n_p_goodN_Step0_epFDn = new TH1D("theta_n_p_goodN_Step0_epFDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_goodN_Step0_epFDn);
    TH1D *h_theta_n_p_badN_Step0_epFDn = new TH1D("theta_n_p_badN_Step0_epFDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_badN_Step0_epFDn);
    TH2D *h_theta_n_p_VS_P_p_goodN_Step0_epFDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_goodN_Step0_epFDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_goodN_Step0_epFDn);
    TH2D *h_theta_n_p_VS_P_p_badN_Step0_epFDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_badN_Step0_epFDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_badN_Step0_epFDn);

    TH1D *h_xB_goodN_Step0_epCDn = new TH1D("xB_goodN_Step0_epCDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_goodN_Step0_epCDn);
    TH1D *h_xB_badN_Step0_epCDn = new TH1D("xB_badN_Step0_epCDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_badN_Step0_epCDn);

    TH1D *h_xB_goodN_Step0_epFDn = new TH1D("xB_goodN_Step0_epFDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_goodN_Step0_epFDn);
    TH1D *h_xB_badN_Step0_epFDn = new TH1D("xB_badN_Step0_epFDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_badN_Step0_epFDn);

    /* Detector responses */
    TH1D *h_Edep_CND_goodN_Step0_epCDn = new TH1D("Edep_CND_goodN_Step0_epCDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_goodN_Step0_epCDn);
    TH1D *h_Edep_CND_badN_Step0_epCDn = new TH1D("Edep_CND_badN_Step0_epCDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_badN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("P_n_VS_Edep_CND_goodN_Step0_epCDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CND_badN_Step0_epCDn = new TH2D("P_n_VS_Edep_CND_badN_Step0_epCDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CND_goodN_Step0_epCDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CND_badN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CND_badN_Step0_epCDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CND_goodN_Step0_epCDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CND_badN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CND_badN_Step0_epCDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CND_goodN_Step0_epCDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CND_badN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CND_badN_Step0_epCDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CND_goodN_Step0_epCDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CND_badN_Step0_epCDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CND_goodN_Step0_epCDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CND_badN_Step0_epCDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("dpp_VS_Edep_CND_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CND_badN_Step0_epCDn = new TH2D("dpp_VS_Edep_CND_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CND_goodN_Step0_epCDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CND_badN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CND_badN_Step0_epCDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("E_p_VS_Edep_CND_goodN_Step0_epCDn", "E_{p} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CND_badN_Step0_epCDn = new TH2D("E_p_VS_Edep_CND_badN_Step0_epCDn", "E_{p} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CND_goodN_Step0_epCDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CND_badN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CND_badN_Step0_epCDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND_goodN_Step0_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND_badN_Step0_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("path_VS_Edep_CND_goodN_Step0_epCDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CND_badN_Step0_epCDn = new TH2D("path_VS_Edep_CND_badN_Step0_epCDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CND_goodN_Step0_epCDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CND_badN_Step0_epCDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("ToF_VS_Edep_CND_goodN_Step0_epCDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CND_badN_Step0_epCDn = new TH2D("ToF_VS_Edep_CND_badN_Step0_epCDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("nSector_VS_Edep_CND_goodN_Step0_epCDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CND_badN_Step0_epCDn = new TH2D("nSector_VS_Edep_CND_badN_Step0_epCDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("Edep_CND1_VS_Edep_CND_goodN_Step0_epCDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step0_epCDn = new TH2D("Edep_CND1_VS_Edep_CND_badN_Step0_epCDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("Edep_CND2_VS_Edep_CND_goodN_Step0_epCDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step0_epCDn = new TH2D("Edep_CND2_VS_Edep_CND_badN_Step0_epCDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_badN_Step0_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("Edep_CND3_VS_Edep_CND_goodN_Step0_epCDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step0_epCDn = new TH2D("Edep_CND3_VS_Edep_CND_badN_Step0_epCDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_badN_Step0_epCDn);

    TH1D *h_Edep_CND_goodN_Step0_epFDn = new TH1D("Edep_CND_goodN_Step0_epFDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_goodN_Step0_epFDn);
    TH1D *h_Edep_CND_badN_Step0_epFDn = new TH1D("Edep_CND_badN_Step0_epFDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_badN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("P_n_VS_Edep_CND_goodN_Step0_epFDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CND_badN_Step0_epFDn = new TH2D("P_n_VS_Edep_CND_badN_Step0_epFDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CND_goodN_Step0_epFDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CND_badN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CND_badN_Step0_epFDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CND_goodN_Step0_epFDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CND_badN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CND_badN_Step0_epFDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CND_goodN_Step0_epFDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CND_badN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CND_badN_Step0_epFDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CND_goodN_Step0_epFDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CND_badN_Step0_epFDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CND_goodN_Step0_epFDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CND_badN_Step0_epFDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("dpp_VS_Edep_CND_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CND_badN_Step0_epFDn = new TH2D("dpp_VS_Edep_CND_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CND_goodN_Step0_epFDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CND_badN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CND_badN_Step0_epFDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("E_p_VS_Edep_CND_goodN_Step0_epFDn", "E_{p} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CND_badN_Step0_epFDn = new TH2D("E_p_VS_Edep_CND_badN_Step0_epFDn", "E_{p} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CND_goodN_Step0_epFDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CND_badN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CND_badN_Step0_epFDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND_goodN_Step0_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND_badN_Step0_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("path_VS_Edep_CND_goodN_Step0_epFDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CND_badN_Step0_epFDn = new TH2D("path_VS_Edep_CND_badN_Step0_epFDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CND_goodN_Step0_epFDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CND_badN_Step0_epFDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("ToF_VS_Edep_CND_goodN_Step0_epFDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CND_badN_Step0_epFDn = new TH2D("ToF_VS_Edep_CND_badN_Step0_epFDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("nSector_VS_Edep_CND_goodN_Step0_epFDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CND_badN_Step0_epFDn = new TH2D("nSector_VS_Edep_CND_badN_Step0_epFDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("Edep_CND1_VS_Edep_CND_goodN_Step0_epFDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step0_epFDn = new TH2D("Edep_CND1_VS_Edep_CND_badN_Step0_epFDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("Edep_CND2_VS_Edep_CND_goodN_Step0_epFDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step0_epFDn = new TH2D("Edep_CND2_VS_Edep_CND_badN_Step0_epFDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_badN_Step0_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("Edep_CND3_VS_Edep_CND_goodN_Step0_epFDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step0_epFDn = new TH2D("Edep_CND3_VS_Edep_CND_badN_Step0_epFDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_badN_Step0_epFDn);

    TH1D *h_Edep_CTOF_goodN_Step0_epCDn = new TH1D("Edep_CTOF_goodN_Step0_epCDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_goodN_Step0_epCDn);
    TH1D *h_Edep_CTOF_badN_Step0_epCDn = new TH1D("Edep_CTOF_badN_Step0_epCDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("P_n_VS_Edep_CTOF_goodN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("P_n_VS_Edep_CTOF_badN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CTOF_goodN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CTOF_badN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CTOF_goodN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CTOF_badN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CTOF_goodN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CTOF_badN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CTOF_goodN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CTOF_badN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CTOF_goodN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CTOF_badN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("dpp_VS_Edep_CTOF_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("dpp_VS_Edep_CTOF_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CTOF_goodN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CTOF_badN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("E_p_VS_Edep_CTOF_goodN_Step0_epCDn", "E_{p} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("E_p_VS_Edep_CTOF_badN_Step0_epCDn", "E_{p} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CTOF_goodN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CTOF_badN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CTOF_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CTOF_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("path_VS_Edep_CTOF_goodN_Step0_epCDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("path_VS_Edep_CTOF_badN_Step0_epCDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CTOF_goodN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CTOF_badN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("ToF_VS_Edep_CTOF_goodN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("ToF_VS_Edep_CTOF_badN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("nSector_VS_Edep_CTOF_goodN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("nSector_VS_Edep_CTOF_badN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("Edep_CND1_VS_Edep_CTOF_goodN_Step0_epCDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("Edep_CND1_VS_Edep_CTOF_badN_Step0_epCDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("Edep_CND2_VS_Edep_CTOF_goodN_Step0_epCDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("Edep_CND2_VS_Edep_CTOF_badN_Step0_epCDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_badN_Step0_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("Edep_CND3_VS_Edep_CTOF_goodN_Step0_epCDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("Edep_CND3_VS_Edep_CTOF_badN_Step0_epCDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_badN_Step0_epCDn);

    TH1D *h_Edep_CTOF_goodN_Step0_epFDn = new TH1D("Edep_CTOF_goodN_Step0_epFDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_goodN_Step0_epFDn);
    TH1D *h_Edep_CTOF_badN_Step0_epFDn = new TH1D("Edep_CTOF_badN_Step0_epFDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("P_n_VS_Edep_CTOF_goodN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("P_n_VS_Edep_CTOF_badN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CTOF_goodN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CTOF_badN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CTOF_goodN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CTOF_badN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CTOF_goodN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CTOF_badN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CTOF_goodN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CTOF_badN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CTOF_goodN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CTOF_badN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("dpp_VS_Edep_CTOF_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("dpp_VS_Edep_CTOF_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CTOF_goodN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CTOF_badN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("E_p_VS_Edep_CTOF_goodN_Step0_epFDn", "E_{p} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("E_p_VS_Edep_CTOF_badN_Step0_epFDn", "E_{p} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CTOF_goodN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CTOF_badN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CTOF_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CTOF_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("path_VS_Edep_CTOF_goodN_Step0_epFDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("path_VS_Edep_CTOF_badN_Step0_epFDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CTOF_goodN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CTOF_badN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("ToF_VS_Edep_CTOF_goodN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("ToF_VS_Edep_CTOF_badN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("nSector_VS_Edep_CTOF_goodN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("nSector_VS_Edep_CTOF_badN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("Edep_CND1_VS_Edep_CTOF_goodN_Step0_epFDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("Edep_CND1_VS_Edep_CTOF_badN_Step0_epFDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("Edep_CND2_VS_Edep_CTOF_goodN_Step0_epFDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("Edep_CND2_VS_Edep_CTOF_badN_Step0_epFDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_badN_Step0_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("Edep_CND3_VS_Edep_CTOF_goodN_Step0_epFDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("Edep_CND3_VS_Edep_CTOF_badN_Step0_epFDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_badN_Step0_epFDn);

    TH1D *h_Edep_single_goodN_Step0_epCDn = new TH1D("Edep_single_goodN_Step0_epCDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_goodN_Step0_epCDn);
    TH1D *h_Edep_single_badN_Step0_epCDn = new TH1D("Edep_single_badN_Step0_epCDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_badN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_single_goodN_Step0_epCDn = new TH2D("P_n_VS_Edep_single_goodN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_single_badN_Step0_epCDn = new TH2D("P_n_VS_Edep_single_badN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_single_goodN_Step0_epCDn = new TH2D("theta_n_VS_Edep_single_goodN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_single_badN_Step0_epCDn = new TH2D("theta_n_VS_Edep_single_badN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_single_goodN_Step0_epCDn = new TH2D("phi_n_VS_Edep_single_goodN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_single_badN_Step0_epCDn = new TH2D("phi_n_VS_Edep_single_badN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_single_goodN_Step0_epCDn = new TH2D("P_miss_VS_Edep_single_goodN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_single_badN_Step0_epCDn = new TH2D("P_miss_VS_Edep_single_badN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_single_goodN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_single_badN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_single_badN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_single_goodN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_single_badN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_single_badN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_single_goodN_Step0_epCDn = new TH2D("dpp_VS_Edep_single_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_single_badN_Step0_epCDn = new TH2D("dpp_VS_Edep_single_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_single_goodN_Step0_epCDn = new TH2D("beta_n_VS_Edep_single_goodN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_single_badN_Step0_epCDn = new TH2D("beta_n_VS_Edep_single_badN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_single_goodN_Step0_epCDn = new TH2D("E_p_VS_Edep_single_goodN_Step0_epCDn", "E_{p} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_single_badN_Step0_epCDn = new TH2D("E_p_VS_Edep_single_badN_Step0_epCDn", "E_{p} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_single_goodN_Step0_epCDn = new TH2D("E_miss_VS_Edep_single_goodN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_single_badN_Step0_epCDn = new TH2D("E_miss_VS_Edep_single_badN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_single_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_single_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_single_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_single_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_path_VS_Edep_single_goodN_Step0_epCDn = new TH2D("path_VS_Edep_single_goodN_Step0_epCDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_path_VS_Edep_single_badN_Step0_epCDn = new TH2D("path_VS_Edep_single_badN_Step0_epCDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_single_goodN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_single_badN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_single_goodN_Step0_epCDn = new TH2D("ToF_VS_Edep_single_goodN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_single_badN_Step0_epCDn = new TH2D("ToF_VS_Edep_single_badN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_badN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_single_goodN_Step0_epCDn = new TH2D("nSector_VS_Edep_single_goodN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_single_badN_Step0_epCDn = new TH2D("nSector_VS_Edep_single_badN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_badN_Step0_epCDn);

    TH1D *h_Edep_single_goodN_Step0_epFDn = new TH1D("Edep_single_goodN_Step0_epFDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_goodN_Step0_epFDn);
    TH1D *h_Edep_single_badN_Step0_epFDn = new TH1D("Edep_single_badN_Step0_epFDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_badN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_single_goodN_Step0_epFDn = new TH2D("P_n_VS_Edep_single_goodN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_single_badN_Step0_epFDn = new TH2D("P_n_VS_Edep_single_badN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_single_goodN_Step0_epFDn = new TH2D("theta_n_VS_Edep_single_goodN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_single_badN_Step0_epFDn = new TH2D("theta_n_VS_Edep_single_badN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_single_goodN_Step0_epFDn = new TH2D("phi_n_VS_Edep_single_goodN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_single_badN_Step0_epFDn = new TH2D("phi_n_VS_Edep_single_badN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_single_goodN_Step0_epFDn = new TH2D("P_miss_VS_Edep_single_goodN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_single_badN_Step0_epFDn = new TH2D("P_miss_VS_Edep_single_badN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_single_goodN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_single_badN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_single_badN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_single_goodN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_single_badN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_single_badN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_single_goodN_Step0_epFDn = new TH2D("dpp_VS_Edep_single_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_single_badN_Step0_epFDn = new TH2D("dpp_VS_Edep_single_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_single_goodN_Step0_epFDn = new TH2D("beta_n_VS_Edep_single_goodN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_single_badN_Step0_epFDn = new TH2D("beta_n_VS_Edep_single_badN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_single_goodN_Step0_epFDn = new TH2D("E_p_VS_Edep_single_goodN_Step0_epFDn", "E_{P} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_single_badN_Step0_epFDn = new TH2D("E_p_VS_Edep_single_badN_Step0_epFDn", "E_{P} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_single_goodN_Step0_epFDn = new TH2D("E_miss_VS_Edep_single_goodN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_single_badN_Step0_epFDn = new TH2D("E_miss_VS_Edep_single_badN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_single_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_single_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_single_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_single_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_path_VS_Edep_single_goodN_Step0_epFDn = new TH2D("path_VS_Edep_single_goodN_Step0_epFDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_path_VS_Edep_single_badN_Step0_epFDn = new TH2D("path_VS_Edep_single_badN_Step0_epFDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_single_goodN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_single_badN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_single_goodN_Step0_epFDn = new TH2D("ToF_VS_Edep_single_goodN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_single_badN_Step0_epFDn = new TH2D("ToF_VS_Edep_single_badN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_badN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_single_goodN_Step0_epFDn = new TH2D("nSector_VS_Edep_single_goodN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_single_badN_Step0_epFDn = new TH2D("nSector_VS_Edep_single_badN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_badN_Step0_epFDn);

    TH1D *h_Edep_CND1_goodN_Step0_epCDn = new TH1D("Edep_CND1_goodN_Step0_epCDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_goodN_Step0_epCDn);
    TH1D *h_Edep_CND1_badN_Step0_epCDn = new TH1D("Edep_CND1_badN_Step0_epCDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("P_n_VS_Edep_CND1_goodN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("P_n_VS_Edep_CND1_badN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CND1_goodN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CND1_badN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CND1_goodN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CND1_badN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CND1_goodN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CND1_badN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CND1_goodN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CND1_badN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CND1_goodN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CND1_badN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("dpp_VS_Edep_CND1_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("dpp_VS_Edep_CND1_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CND1_goodN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CND1_badN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("E_p_VS_Edep_CND1_goodN_Step0_epCDn", "E_{P} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("E_p_VS_Edep_CND1_badN_Step0_epCDn", "E_{P} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CND1_goodN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CND1_badN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND1_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND1_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("path_VS_Edep_CND1_goodN_Step0_epCDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("path_VS_Edep_CND1_badN_Step0_epCDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CND1_goodN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CND1_badN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("ToF_VS_Edep_CND1_goodN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("ToF_VS_Edep_CND1_badN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("nSector_VS_Edep_CND1_goodN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("nSector_VS_Edep_CND1_badN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("Edep_CND2_VS_Edep_CND1_goodN_Step0_epCDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("Edep_CND2_VS_Edep_CND1_badN_Step0_epCDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_badN_Step0_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("Edep_CND3_VS_Edep_CND1_goodN_Step0_epCDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("Edep_CND3_VS_Edep_CND1_badN_Step0_epCDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_badN_Step0_epCDn);

    TH1D *h_Edep_CND1_goodN_Step0_epFDn = new TH1D("Edep_CND1_goodN_Step0_epFDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_goodN_Step0_epFDn);
    TH1D *h_Edep_CND1_badN_Step0_epFDn = new TH1D("Edep_CND1_badN_Step0_epFDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("P_n_VS_Edep_CND1_goodN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("P_n_VS_Edep_CND1_badN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CND1_goodN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CND1_badN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CND1_goodN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CND1_badN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CND1_goodN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CND1_badN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CND1_goodN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CND1_badN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CND1_goodN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CND1_badN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("dpp_VS_Edep_CND1_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("dpp_VS_Edep_CND1_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CND1_goodN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CND1_badN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("E_p_VS_Edep_CND1_goodN_Step0_epFDn", "E_{P} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("E_p_VS_Edep_CND1_badN_Step0_epFDn", "E_{P} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CND1_goodN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CND1_badN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND1_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND1_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("path_VS_Edep_CND1_goodN_Step0_epFDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("path_VS_Edep_CND1_badN_Step0_epFDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CND1_goodN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CND1_badN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("ToF_VS_Edep_CND1_goodN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("ToF_VS_Edep_CND1_badN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("nSector_VS_Edep_CND1_goodN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("nSector_VS_Edep_CND1_badN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("Edep_CND2_VS_Edep_CND1_goodN_Step0_epFDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("Edep_CND2_VS_Edep_CND1_badN_Step0_epFDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_badN_Step0_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("Edep_CND3_VS_Edep_CND1_goodN_Step0_epFDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("Edep_CND3_VS_Edep_CND1_badN_Step0_epFDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_badN_Step0_epFDn);

    TH1D *h_Edep_CND2_goodN_Step0_epCDn = new TH1D("Edep_CND2_goodN_Step0_epCDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_goodN_Step0_epCDn);
    TH1D *h_Edep_CND2_badN_Step0_epCDn = new TH1D("Edep_CND2_badN_Step0_epCDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("P_n_VS_Edep_CND2_goodN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("P_n_VS_Edep_CND2_badN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CND2_goodN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CND2_badN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CND2_goodN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CND2_badN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CND2_goodN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CND2_badN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CND2_goodN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CND2_badN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CND2_goodN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CND2_badN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("dpp_VS_Edep_CND2_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("dpp_VS_Edep_CND2_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CND2_goodN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CND2_badN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("E_p_VS_Edep_CND2_goodN_Step0_epCDn", "E_{P} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("E_p_VS_Edep_CND2_badN_Step0_epCDn", "E_{P} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CND2_goodN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CND2_badN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND2_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND2_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("path_VS_Edep_CND2_goodN_Step0_epCDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("path_VS_Edep_CND2_badN_Step0_epCDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CND2_goodN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CND2_badN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("ToF_VS_Edep_CND2_goodN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("ToF_VS_Edep_CND2_badN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("nSector_VS_Edep_CND2_goodN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("nSector_VS_Edep_CND2_badN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_badN_Step0_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("Edep_CND3_VS_Edep_CND2_goodN_Step0_epCDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("Edep_CND3_VS_Edep_CND2_badN_Step0_epCDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_badN_Step0_epCDn);

    TH1D *h_Edep_CND2_goodN_Step0_epFDn = new TH1D("Edep_CND2_goodN_Step0_epFDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_goodN_Step0_epFDn);
    TH1D *h_Edep_CND2_badN_Step0_epFDn = new TH1D("Edep_CND2_badN_Step0_epFDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("P_n_VS_Edep_CND2_goodN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("P_n_VS_Edep_CND2_badN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CND2_goodN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CND2_badN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CND2_goodN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CND2_badN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CND2_goodN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CND2_badN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CND2_goodN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CND2_badN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CND2_goodN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CND2_badN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("dpp_VS_Edep_CND2_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("dpp_VS_Edep_CND2_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CND2_goodN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CND2_badN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("E_p_VS_Edep_CND2_goodN_Step0_epFDn", "E_{P} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("E_p_VS_Edep_CND2_badN_Step0_epFDn", "E_{P} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CND2_goodN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CND2_badN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND2_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND2_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("path_VS_Edep_CND2_goodN_Step0_epFDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("path_VS_Edep_CND2_badN_Step0_epFDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CND2_goodN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CND2_badN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("ToF_VS_Edep_CND2_goodN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("ToF_VS_Edep_CND2_badN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("nSector_VS_Edep_CND2_goodN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("nSector_VS_Edep_CND2_badN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_badN_Step0_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("Edep_CND3_VS_Edep_CND2_goodN_Step0_epFDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("Edep_CND3_VS_Edep_CND2_badN_Step0_epFDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_badN_Step0_epFDn);

    TH1D *h_Edep_CND3_goodN_Step0_epCDn = new TH1D("Edep_CND3_goodN_Step0_epCDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_goodN_Step0_epCDn);
    TH1D *h_Edep_CND3_badN_Step0_epCDn = new TH1D("Edep_CND3_badN_Step0_epCDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("P_n_VS_Edep_CND3_goodN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("P_n_VS_Edep_CND3_badN_Step0_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CND3_goodN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("theta_n_VS_Edep_CND3_badN_Step0_epCDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CND3_goodN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("phi_n_VS_Edep_CND3_badN_Step0_epCDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CND3_goodN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("P_miss_VS_Edep_CND3_badN_Step0_epCDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CND3_goodN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("theta_miss_VS_Edep_CND3_badN_Step0_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CND3_goodN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("phi_miss_VS_Edep_CND3_badN_Step0_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("dpp_VS_Edep_CND3_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_dpp_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("dpp_VS_Edep_CND3_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CND3_goodN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("beta_n_VS_Edep_CND3_badN_Step0_epCDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("E_p_VS_Edep_CND3_goodN_Step0_epCDn", "E_{p} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_E_p_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("E_p_VS_Edep_CND3_badN_Step0_epCDn", "E_{p} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CND3_goodN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("E_miss_VS_Edep_CND3_badN_Step0_epCDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND3_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND3_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("path_VS_Edep_CND3_goodN_Step0_epCDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_path_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("path_VS_Edep_CND3_badN_Step0_epCDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CND3_goodN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("theta_n_miss_VS_Edep_CND3_badN_Step0_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("ToF_VS_Edep_CND3_goodN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_ToF_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("ToF_VS_Edep_CND3_badN_Step0_epCDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_badN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("nSector_VS_Edep_CND3_goodN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_nSector_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("nSector_VS_Edep_CND3_badN_Step0_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_badN_Step0_epCDn);

    TH1D *h_Edep_CND3_goodN_Step0_epFDn = new TH1D("Edep_CND3_goodN_Step0_epFDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_goodN_Step0_epFDn);
    TH1D *h_Edep_CND3_badN_Step0_epFDn = new TH1D("Edep_CND3_badN_Step0_epFDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("P_n_VS_Edep_CND3_goodN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("P_n_VS_Edep_CND3_badN_Step0_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CND3_goodN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("theta_n_VS_Edep_CND3_badN_Step0_epFDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CND3_goodN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("phi_n_VS_Edep_CND3_badN_Step0_epFDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CND3_goodN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("P_miss_VS_Edep_CND3_badN_Step0_epFDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CND3_goodN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("theta_miss_VS_Edep_CND3_badN_Step0_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CND3_goodN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("phi_miss_VS_Edep_CND3_badN_Step0_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("dpp_VS_Edep_CND3_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_dpp_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("dpp_VS_Edep_CND3_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CND3_goodN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("beta_n_VS_Edep_CND3_badN_Step0_epFDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("E_p_VS_Edep_CND3_goodN_Step0_epFDn", "E_{p} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_E_p_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("E_p_VS_Edep_CND3_badN_Step0_epFDn", "E_{p} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CND3_goodN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("E_miss_VS_Edep_CND3_badN_Step0_epFDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND3_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND3_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("path_VS_Edep_CND3_goodN_Step0_epFDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_path_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("path_VS_Edep_CND3_badN_Step0_epFDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CND3_goodN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("theta_n_miss_VS_Edep_CND3_badN_Step0_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("ToF_VS_Edep_CND3_goodN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_ToF_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("ToF_VS_Edep_CND3_badN_Step0_epFDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_badN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("nSector_VS_Edep_CND3_goodN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_nSector_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("nSector_VS_Edep_CND3_badN_Step0_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_badN_Step0_epFDn);

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step0_epCDn = new TH2D("Size_CND1_VS_Size_CND2_goodN_Step0_epCDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_goodN_Step0_epCDn);
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step0_epCDn = new TH2D("Size_CND1_VS_Size_CND2_badN_Step0_epCDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_badN_Step0_epCDn);
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step0_epCDn = new TH2D("Size_CND1_VS_Size_CND3_goodN_Step0_epCDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_goodN_Step0_epCDn);
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step0_epCDn = new TH2D("Size_CND1_VS_Size_CND3_badN_Step0_epCDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_badN_Step0_epCDn);
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step0_epCDn = new TH2D("Size_CND2_VS_Size_CND3_goodN_Step0_epCDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_goodN_Step0_epCDn);
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step0_epCDn = new TH2D("Size_CND2_VS_Size_CND3_badN_Step0_epCDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_badN_Step0_epCDn);

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step0_epFDn = new TH2D("Size_CND1_VS_Size_CND2_goodN_Step0_epFDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_goodN_Step0_epFDn);
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step0_epFDn = new TH2D("Size_CND1_VS_Size_CND2_badN_Step0_epFDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_badN_Step0_epFDn);
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step0_epFDn = new TH2D("Size_CND1_VS_Size_CND3_goodN_Step0_epFDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_goodN_Step0_epFDn);
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step0_epFDn = new TH2D("Size_CND1_VS_Size_CND3_badN_Step0_epFDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_badN_Step0_epFDn);
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step0_epFDn = new TH2D("Size_CND2_VS_Size_CND3_goodN_Step0_epFDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_goodN_Step0_epFDn);
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step0_epFDn = new TH2D("Size_CND2_VS_Size_CND3_badN_Step0_epFDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_badN_Step0_epFDn);

    TH1D *h_ToF_goodN_Step0_epCDn = new TH1D("ToF_goodN_Step0_epCDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_goodN_Step0_epCDn);
    TH1D *h_ToF_badN_Step0_epCDn = new TH1D("ToF_badN_Step0_epCDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_badN_Step0_epCDn);
    TH2D *h_P_n_VS_ToF_goodN_Step0_epCDn = new TH2D("P_n_VS_ToF_goodN_Step0_epCDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_ToF_badN_Step0_epCDn = new TH2D("P_n_VS_ToF_badN_Step0_epCDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_badN_Step0_epCDn);
    TH2D *h_theta_n_VS_ToF_goodN_Step0_epCDn = new TH2D("theta_n_VS_ToF_goodN_Step0_epCDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_theta_n_VS_ToF_badN_Step0_epCDn = new TH2D("theta_n_VS_ToF_badN_Step0_epCDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_badN_Step0_epCDn);
    TH2D *h_phi_n_VS_ToF_goodN_Step0_epCDn = new TH2D("phi_n_VS_ToF_goodN_Step0_epCDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_phi_n_VS_ToF_badN_Step0_epCDn = new TH2D("phi_n_VS_ToF_badN_Step0_epCDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_badN_Step0_epCDn);
    TH2D *h_P_miss_VS_ToF_goodN_Step0_epCDn = new TH2D("P_miss_VS_ToF_goodN_Step0_epCDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_P_miss_VS_ToF_badN_Step0_epCDn = new TH2D("P_miss_VS_ToF_badN_Step0_epCDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_badN_Step0_epCDn);
    TH2D *h_theta_miss_VS_ToF_goodN_Step0_epCDn = new TH2D("theta_miss_VS_ToF_goodN_Step0_epCDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_theta_miss_VS_ToF_badN_Step0_epCDn = new TH2D("theta_miss_VS_ToF_badN_Step0_epCDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_badN_Step0_epCDn);
    TH2D *h_phi_miss_VS_ToF_goodN_Step0_epCDn = new TH2D("phi_miss_VS_ToF_goodN_Step0_epCDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_phi_miss_VS_ToF_badN_Step0_epCDn = new TH2D("phi_miss_VS_ToF_badN_Step0_epCDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_badN_Step0_epCDn);
    TH2D *h_dpp_VS_ToF_goodN_Step0_epCDn = new TH2D("dpp_VS_ToF_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_dpp_VS_ToF_badN_Step0_epCDn = new TH2D("dpp_VS_ToF_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_badN_Step0_epCDn);
    TH2D *h_beta_n_VS_ToF_goodN_Step0_epCDn = new TH2D("beta_n_VS_ToF_goodN_Step0_epCDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_beta_n_VS_ToF_badN_Step0_epCDn = new TH2D("beta_n_VS_ToF_badN_Step0_epCDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_badN_Step0_epCDn);
    TH2D *h_E_p_VS_ToF_goodN_Step0_epCDn = new TH2D("E_p_VS_ToF_goodN_Step0_epCDn", "E_{p} vs ToF;ToF [ns];E_{p} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_E_p_VS_ToF_badN_Step0_epCDn = new TH2D("E_p_VS_ToF_badN_Step0_epCDn", "E_{p} vs ToF;ToF [ns];E_{p} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_ToF_badN_Step0_epCDn);
    TH2D *h_E_miss_VS_ToF_goodN_Step0_epCDn = new TH2D("E_miss_VS_ToF_goodN_Step0_epCDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_E_miss_VS_ToF_badN_Step0_epCDn = new TH2D("E_miss_VS_ToF_badN_Step0_epCDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_ToF_goodN_Step0_epCDn = new TH2D("M_miss_VS_ToF_goodN_Step0_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_ToF_badN_Step0_epCDn = new TH2D("M_miss_VS_ToF_badN_Step0_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_badN_Step0_epCDn);
    TH2D *h_path_VS_ToF_goodN_Step0_epCDn = new TH2D("path_VS_ToF_goodN_Step0_epCDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_path_VS_ToF_badN_Step0_epCDn = new TH2D("path_VS_ToF_badN_Step0_epCDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_badN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step0_epCDn = new TH2D("theta_n_miss_VS_ToF_goodN_Step0_epCDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_theta_n_miss_VS_ToF_badN_Step0_epCDn = new TH2D("theta_n_miss_VS_ToF_badN_Step0_epCDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_badN_Step0_epCDn);
    TH2D *h_nSector_VS_ToF_goodN_Step0_epCDn = new TH2D("nSector_VS_ToF_goodN_Step0_epCDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_nSector_VS_ToF_badN_Step0_epCDn = new TH2D("nSector_VS_ToF_badN_Step0_epCDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_badN_Step0_epCDn);

    TH1D *h_ToF_goodN_Step0_epFDn = new TH1D("ToF_goodN_Step0_epFDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_goodN_Step0_epFDn);
    TH1D *h_ToF_badN_Step0_epFDn = new TH1D("ToF_badN_Step0_epFDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_badN_Step0_epFDn);
    TH2D *h_P_n_VS_ToF_goodN_Step0_epFDn = new TH2D("P_n_VS_ToF_goodN_Step0_epFDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_ToF_badN_Step0_epFDn = new TH2D("P_n_VS_ToF_badN_Step0_epFDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_badN_Step0_epFDn);
    TH2D *h_theta_n_VS_ToF_goodN_Step0_epFDn = new TH2D("theta_n_VS_ToF_goodN_Step0_epFDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_theta_n_VS_ToF_badN_Step0_epFDn = new TH2D("theta_n_VS_ToF_badN_Step0_epFDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_badN_Step0_epFDn);
    TH2D *h_phi_n_VS_ToF_goodN_Step0_epFDn = new TH2D("phi_n_VS_ToF_goodN_Step0_epFDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_phi_n_VS_ToF_badN_Step0_epFDn = new TH2D("phi_n_VS_ToF_badN_Step0_epFDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_badN_Step0_epFDn);
    TH2D *h_P_miss_VS_ToF_goodN_Step0_epFDn = new TH2D("P_miss_VS_ToF_goodN_Step0_epFDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_P_miss_VS_ToF_badN_Step0_epFDn = new TH2D("P_miss_VS_ToF_badN_Step0_epFDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_badN_Step0_epFDn);
    TH2D *h_theta_miss_VS_ToF_goodN_Step0_epFDn = new TH2D("theta_miss_VS_ToF_goodN_Step0_epFDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_theta_miss_VS_ToF_badN_Step0_epFDn = new TH2D("theta_miss_VS_ToF_badN_Step0_epFDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_badN_Step0_epFDn);
    TH2D *h_phi_miss_VS_ToF_goodN_Step0_epFDn = new TH2D("phi_miss_VS_ToF_goodN_Step0_epFDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_phi_miss_VS_ToF_badN_Step0_epFDn = new TH2D("phi_miss_VS_ToF_badN_Step0_epFDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_badN_Step0_epFDn);
    TH2D *h_dpp_VS_ToF_goodN_Step0_epFDn = new TH2D("dpp_VS_ToF_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_dpp_VS_ToF_badN_Step0_epFDn = new TH2D("dpp_VS_ToF_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_badN_Step0_epFDn);
    TH2D *h_beta_n_VS_ToF_goodN_Step0_epFDn = new TH2D("beta_n_VS_ToF_goodN_Step0_epFDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_beta_n_VS_ToF_badN_Step0_epFDn = new TH2D("beta_n_VS_ToF_badN_Step0_epFDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_badN_Step0_epFDn);
    TH2D *h_E_p_VS_ToF_goodN_Step0_epFDn = new TH2D("E_p_VS_ToF_goodN_Step0_epFDn", "E_{p} vs ToF;ToF [ns];E_{p} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_E_p_VS_ToF_badN_Step0_epFDn = new TH2D("E_p_VS_ToF_badN_Step0_epFDn", "E_{p} vs ToF;ToF [ns];E_{p} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_ToF_badN_Step0_epFDn);
    TH2D *h_E_miss_VS_ToF_goodN_Step0_epFDn = new TH2D("E_miss_VS_ToF_goodN_Step0_epFDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_E_miss_VS_ToF_badN_Step0_epFDn = new TH2D("E_miss_VS_ToF_badN_Step0_epFDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_ToF_goodN_Step0_epFDn = new TH2D("M_miss_VS_ToF_goodN_Step0_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_ToF_badN_Step0_epFDn = new TH2D("M_miss_VS_ToF_badN_Step0_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_badN_Step0_epFDn);
    TH2D *h_path_VS_ToF_goodN_Step0_epFDn = new TH2D("path_VS_ToF_goodN_Step0_epFDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_path_VS_ToF_badN_Step0_epFDn = new TH2D("path_VS_ToF_badN_Step0_epFDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_badN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step0_epFDn = new TH2D("theta_n_miss_VS_ToF_goodN_Step0_epFDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_theta_n_miss_VS_ToF_badN_Step0_epFDn = new TH2D("theta_n_miss_VS_ToF_badN_Step0_epFDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_badN_Step0_epFDn);
    TH2D *h_nSector_VS_ToF_goodN_Step0_epFDn = new TH2D("nSector_VS_ToF_goodN_Step0_epFDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_nSector_VS_ToF_badN_Step0_epFDn = new TH2D("nSector_VS_ToF_badN_Step0_epFDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_badN_Step0_epFDn);

    TH1D *h_beta_n_goodN_Step0_epCDn = new TH1D("beta_n_goodN_Step0_epCDn", "#beta_{n} of CND Neutrons;#beta_{n};Counts", 50, 0, 1.1);
    HistoList.push_back(h_beta_n_goodN_Step0_epCDn);
    TH1D *h_beta_n_badN_Step0_epCDn = new TH1D("beta_n_badN_Step0_epCDn", "#beta_{n} of CND Neutrons;#beta_{n};Counts", 50, 0, 1.1);
    HistoList.push_back(h_beta_n_badN_Step0_epCDn);

    TH1D *h_beta_n_goodN_Step0_epFDn = new TH1D("beta_n_goodN_Step0_epFDn", "#beta_{n} of CND Neutrons;#beta_{n};Counts", 50, 0, 1.1);
    HistoList.push_back(h_beta_n_goodN_Step0_epFDn);
    TH1D *h_beta_n_badN_Step0_epFDn = new TH1D("beta_n_badN_Step0_epFDn", "#beta_{n} of CND Neutrons;#beta_{n};Counts", 50, 0, 1.1);
    HistoList.push_back(h_beta_n_badN_Step0_epFDn);

    // TODO: add these to code below
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step0_epCDn = new TH1D("neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step0_epCDn", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 50, 0, 5);
    HistoList.push_back(h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step0_epCDn);
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step0_epCDn = new TH1D("neut_Edep_CND_over_pos_Edep_CTOF_badN_Step0_epCDn", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 50, 0, 5);
    HistoList.push_back(h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step0_epCDn);

    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step0_epFDn = new TH1D("neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step0_epFDn", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 50, 0, 5);
    HistoList.push_back(h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step0_epFDn);
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step0_epFDn = new TH1D("neut_Edep_CND_over_pos_Edep_CTOF_badN_Step0_epFDn", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 50, 0, 5);
    HistoList.push_back(h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step0_epFDn);

    TH1D *h_Edep_CND_goodN_withNearbyPos_Step0_epCDn = new TH1D("Edep_CND_goodN_withNearbyPos_Step0_epCDn", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 50);
    HistoList.push_back(h_Edep_CND_goodN_withNearbyPos_Step0_epCDn);
    TH1D *h_Edep_CND_badN_withNearbyPos_Step0_epCDn = new TH1D("Edep_CND_badN_withNearbyPos_Step0_epCDn", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 50);
    HistoList.push_back(h_Edep_CND_badN_withNearbyPos_Step0_epCDn);

    TH1D *h_Edep_CND_goodN_withNearbyPos_Step0_epFDn = new TH1D("Edep_CND_goodN_withNearbyPos_Step0_epFDn", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 50);
    HistoList.push_back(h_Edep_CND_goodN_withNearbyPos_Step0_epFDn);
    TH1D *h_Edep_CND_badN_withNearbyPos_Step0_epFDn = new TH1D("Edep_CND_badN_withNearbyPos_Step0_epFDn", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 50);
    HistoList.push_back(h_Edep_CND_badN_withNearbyPos_Step0_epFDn);

    TH1D *h_sdiff_pos_goodN_Step0_layer_epCDn[7];
    TH1D *h_sdiff_pos_badN_Step0_layer_epCDn[7];

    TH1D *h_sdiff_pos_goodN_Step0_layer_epFDn[7];
    TH1D *h_sdiff_pos_badN_Step0_layer_epFDn[7];

    TH2D *h_sdiff_pos_mom_goodN_Step0_layer_epCDn[7];
    TH2D *h_sdiff_pos_mom_badN_Step0_layer_epCDn[7];

    TH2D *h_sdiff_pos_mom_goodN_Step0_layer_epFDn[7];
    TH2D *h_sdiff_pos_mom_badN_Step0_layer_epFDn[7];

    TH2D *h_sdiff_pos_z_goodN_Step0_layer_epCDn[7];
    TH2D *h_sdiff_pos_z_badN_Step0_layer_epCDn[7];

    TH2D *h_sdiff_pos_z_goodN_Step0_layer_epFDn[7];
    TH2D *h_sdiff_pos_z_badN_Step0_layer_epFDn[7];

    TH2D *h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epCDn[7];
    TH2D *h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epCDn[7];

    TH2D *h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epFDn[7];
    TH2D *h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epFDn[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name, "sdiff_pos_goodN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step0_layer_epCDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_badN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step0_layer_epCDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_goodN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step0_layer_epFDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_badN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step0_layer_epFDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step0_layer_epFDn[k]);

        sprintf(temp_name, "sdiff_pos_mom_goodN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step0_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_mom_badN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step0_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_mom_goodN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step0_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_mom_badN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step0_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step0_layer_epFDn[k]);

        sprintf(temp_name, "sdiff_pos_z_goodN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step0_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_z_badN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step0_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_z_goodN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step0_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_z_badN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step0_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_badN_Step0_layer_epFDn[k]);

        sprintf(temp_name, "sdiff_pos_diff_ToFc_z_goodN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_diff_ToFc_z_badN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_diff_ToFc_z_goodN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_diff_ToFc_z_badN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epFDn[k]);
    }

    TH2D *h_diff_ToFc_z_VS_Edep_noNear_goodN_Step0_epCDn = new TH2D("diff_ToFc_z_VS_Edep_noNear_goodN_Step0_epCDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_noNear_goodN_Step0_epCDn);
    TH2D *h_diff_ToFc_z_VS_Edep_noNear_badN_Step0_epCDn = new TH2D("diff_ToFc_z_VS_Edep_noNear_badN_Step0_epCDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_noNear_badN_Step0_epCDn);

    TH2D *h_diff_ToFc_z_VS_Edep_noNear_goodN_Step0_epFDn = new TH2D("diff_ToFc_z_VS_Edep_noNear_goodN_Step0_epFDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_noNear_goodN_Step0_epFDn);
    TH2D *h_diff_ToFc_z_VS_Edep_noNear_badN_Step0_epFDn = new TH2D("diff_ToFc_z_VS_Edep_noNear_badN_Step0_epFDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_noNear_badN_Step0_epFDn);

    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step0_epCDn = new TH2D("diff_ToFc_z_VS_Edep_yesNear_goodN_Step0_epCDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step0_epCDn);
    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_badN_Step0_epCDn = new TH2D("diff_ToFc_z_VS_Edep_yesNear_badN_Step0_epCDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_yesNear_badN_Step0_epCDn);

    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step0_epFDn = new TH2D("diff_ToFc_z_VS_Edep_yesNear_goodN_Step0_epFDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step0_epFDn);
    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_badN_Step0_epFDn = new TH2D("diff_ToFc_z_VS_Edep_yesNear_badN_Step0_epFDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_yesNear_badN_Step0_epFDn);

    TH2D *h_diff_ToFc_z_Edep_goodN_Step0_layer_epCDn[3];
    TH2D *h_diff_ToFc_z_Edep_badN_Step0_layer_epCDn[3];

    TH2D *h_diff_ToFc_z_Edep_goodN_Step0_layer_epFDn[3];
    TH2D *h_diff_ToFc_z_Edep_badN_Step0_layer_epFDn[3];

    for (int k = 0; k < 3; k++)
    {
        sprintf(temp_name, "diff_ToFc_z_goodN_Step0_layer_%d_epCDn", k + 1);
        sprintf(temp_title, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step0_layer_epCDn[k] = new TH2D(temp_name, temp_title, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name, "diff_ToFc_z_badN_Step0_layer_%d_epCDn", k + 1);
        sprintf(temp_title, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step0_layer_epCDn[k] = new TH2D(temp_name, temp_title, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name, "diff_ToFc_z_goodN_Step0_layer_%d_epFDn", k + 1);
        sprintf(temp_title, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step0_layer_epFDn[k] = new TH2D(temp_name, temp_title, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name, "diff_ToFc_z_badN_Step0_layer_%d_epFDn", k + 1);
        sprintf(temp_title, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step0_layer_epFDn[k] = new TH2D(temp_name, temp_title, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_badN_Step0_layer_epFDn[k]);
    }

    TH2D *h_sdiff_ldiff_allhit_goodN_Step0_epCDn = new TH2D("sdiff_ldiff_allhit_goodN_Step0_epCDn", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    HistoList.push_back(h_sdiff_ldiff_allhit_goodN_Step0_epCDn);
    TH2D *h_sdiff_ldiff_allhit_badN_Step0_epCDn = new TH2D("sdiff_ldiff_allhit_badN_Step0_epCDn", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    HistoList.push_back(h_sdiff_ldiff_allhit_badN_Step0_epCDn);

    TH2D *h_sdiff_ldiff_allhit_goodN_Step0_epFDn = new TH2D("sdiff_ldiff_allhit_goodN_Step0_epFDn", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    HistoList.push_back(h_sdiff_ldiff_allhit_goodN_Step0_epFDn);
    TH2D *h_sdiff_ldiff_allhit_badN_Step0_epFDn = new TH2D("sdiff_ldiff_allhit_badN_Step0_epFDn", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    HistoList.push_back(h_sdiff_ldiff_allhit_badN_Step0_epFDn);

    TH1D *h_numberNearby_goodN_Step0_epCDn = new TH1D("numberNearby_goodN_Step0_epCDn", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    HistoList.push_back(h_numberNearby_goodN_Step0_epCDn);
    TH1D *h_numberNearby_badN_Step0_epCDn = new TH1D("numberNearby_badN_Step0_epCDn", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    HistoList.push_back(h_numberNearby_badN_Step0_epCDn);

    TH1D *h_numberNearby_goodN_Step0_epFDn = new TH1D("numberNearby_goodN_Step0_epFDn", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    HistoList.push_back(h_numberNearby_goodN_Step0_epFDn);
    TH1D *h_numberNearby_badN_Step0_epFDn = new TH1D("numberNearby_badN_Step0_epFDn", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    HistoList.push_back(h_numberNearby_badN_Step0_epFDn);

    TH2D *h_numberNearby_momN_goodN_Step0_epCDn = new TH2D("numberNearby_momN_goodN_Step0_epCDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.5);
    HistoList.push_back(h_numberNearby_momN_goodN_Step0_epCDn);
    TH2D *h_numberNearby_momN_badN_Step0_epCDn = new TH2D("numberNearby_momN_badN_Step0_epCDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.5);
    HistoList.push_back(h_numberNearby_momN_badN_Step0_epCDn);

    TH2D *h_numberNearby_momN_goodN_Step0_epFDn = new TH2D("numberNearby_momN_goodN_Step0_epFDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.5);
    HistoList.push_back(h_numberNearby_momN_goodN_Step0_epFDn);
    TH2D *h_numberNearby_momN_badN_Step0_epFDn = new TH2D("numberNearby_momN_badN_Step0_epFDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.5);
    HistoList.push_back(h_numberNearby_momN_badN_Step0_epFDn);

    TH1D *h_NearbyEdep_goodN_Step0_epCDn = new TH1D("NearbyEdep_goodN_Step0_epCDn", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_NearbyEdep_goodN_Step0_epCDn);
    TH1D *h_NearbyEdep_badN_Step0_epCDn = new TH1D("NearbyEdep_badN_Step0_epCDn", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_NearbyEdep_badN_Step0_epCDn);

    TH1D *h_NearbyEdep_goodN_Step0_epFDn = new TH1D("NearbyEdep_goodN_Step0_epFDn", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_NearbyEdep_goodN_Step0_epFDn);
    TH1D *h_NearbyEdep_badN_Step0_epFDn = new TH1D("NearbyEdep_badN_Step0_epFDn", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_NearbyEdep_badN_Step0_epFDn);

    TH1D *h_nsector_goodN_Step0_epCDn = new TH1D("nsector_goodN_Step0_epCDn", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    HistoList.push_back(h_nsector_goodN_Step0_epCDn);
    TH1D *h_nsector_badN_Step0_epCDn = new TH1D("nsector_badN_Step0_epCDn", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    HistoList.push_back(h_nsector_badN_Step0_epCDn);

    TH1D *h_nsector_goodN_Step0_epFDn = new TH1D("nsector_goodN_Step0_epFDn", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    HistoList.push_back(h_nsector_goodN_Step0_epFDn);
    TH1D *h_nsector_badN_Step0_epFDn = new TH1D("nsector_badN_Step0_epFDn", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    HistoList.push_back(h_nsector_badN_Step0_epFDn);

    TH1D *h_phidiff_en_goodN_Step0_epCDn = new TH1D("phidiff_en_goodN_Step0_epCDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 50, 0, 180);
    HistoList.push_back(h_phidiff_en_goodN_Step0_epCDn);
    TH1D *h_phidiff_en_badN_Step0_epCDn = new TH1D("phidiff_en_badN_Step0_epCDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 50, 0, 180);
    HistoList.push_back(h_phidiff_en_badN_Step0_epCDn);

    TH1D *h_phidiff_en_goodN_Step0_epFDn = new TH1D("phidiff_en_goodN_Step0_epFDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 50, 0, 180);
    HistoList.push_back(h_phidiff_en_goodN_Step0_epFDn);
    TH1D *h_phidiff_en_badN_Step0_epFDn = new TH1D("phidiff_en_badN_Step0_epFDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 50, 0, 180);
    HistoList.push_back(h_phidiff_en_badN_Step0_epFDn);

    TH1D *h_TP_goodN_Step0_epCDn = new TH1D("TP_goodN_Step0_epCDn", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    HistoList.push_back(h_TP_goodN_Step0_epCDn);
    TH1D *h_TP_badN_Step0_epCDn = new TH1D("TP_badN_Step0_epCDn", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    HistoList.push_back(h_TP_badN_Step0_epCDn);

    TH1D *h_TP_goodN_Step0_epFDn = new TH1D("TP_goodN_Step0_epFDn", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    HistoList.push_back(h_TP_goodN_Step0_epFDn);
    TH1D *h_TP_badN_Step0_epFDn = new TH1D("TP_badN_Step0_epFDn", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    HistoList.push_back(h_TP_badN_Step0_epFDn);

    TH1D *h_Z_goodN_Step0_epCDn = new TH1D("Z_goodN_Step0_epCDn", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    HistoList.push_back(h_Z_goodN_Step0_epCDn);
    TH1D *h_Z_badN_Step0_epCDn = new TH1D("Z_badN_Step0_epCDn", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    HistoList.push_back(h_Z_badN_Step0_epCDn);

    TH1D *h_Z_goodN_Step0_epFDn = new TH1D("Z_goodN_Step0_epFDn", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    HistoList.push_back(h_Z_goodN_Step0_epFDn);
    TH1D *h_Z_badN_Step0_epFDn = new TH1D("Z_badN_Step0_epFDn", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    HistoList.push_back(h_Z_badN_Step0_epFDn);

    // Step One (After Beta Cut) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step1 = new TH1D("n_multiplicity_allN_epCDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epCDn_Step1);
    TH1D *h_n_multiplicity_goodN_epCDn_Step1 = new TH1D("n_multiplicity_goodN_epCDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epCDn_Step1);
    TH1D *h_n_multiplicity_badN_epCDn_Step1 = new TH1D("n_multiplicity_badN_epCDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epCDn_Step1);

    TH1D *h_n_multiplicity_allN_epFDn_Step1 = new TH1D("n_multiplicity_allN_epFDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epFDn_Step1);
    TH1D *h_n_multiplicity_goodN_epFDn_Step1 = new TH1D("n_multiplicity_goodN_epFDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epFDn_Step1);
    TH1D *h_n_multiplicity_badN_epFDn_Step1 = new TH1D("n_multiplicity_badN_epFDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epFDn_Step1);

    /*  */
    // TODO: add to code below!
    TH2D *h_dbeta_n_VS_P_n_goodN_Step1_epCDn = new TH2D("dbeta_n_VS_P_n_goodN_Step1_epCDn", "#Delta#beta_{n} vs Neutron Momentum;P_{n} [GeV/c];#Delta#beta_{n}", 50, 0, 1.5, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_P_n_goodN_Step1_epCDn);
    TH2D *h_dbeta_n_VS_ToF_goodN_Step1_epCDn = new TH2D("dbeta_n_VS_ToF_goodN_Step1_epCDn", "#Delta#beta_{n} vs Neutron ToF;ToF [ns];#Delta#beta_{n}", 50, 0, 50, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_dbeta_n_VS_P_n_badN_Step1_epCDn = new TH2D("dbeta_n_VS_P_n_badN_Step1_epCDn", "#Delta#beta_{n} vs Neutron Momentum;P_{n} [GeV/c];#Delta#beta_{n}", 50, 0, 1.5, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_P_n_badN_Step1_epCDn);
    TH2D *h_dbeta_n_VS_ToF_badN_Step1_epCDn = new TH2D("dbeta_n_VS_ToF_badN_Step1_epCDn", "#Delta#beta_{n} vs Neutron ToF;ToF [ns];#Delta#beta_{n}", 50, 0, 50, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_ToF_badN_Step1_epCDn);

    TH2D *h_dbeta_n_VS_P_n_goodN_Step1_epFDn = new TH2D("dbeta_n_VS_P_n_goodN_Step1_epFDn", "#Delta#beta_{n} vs Neutron Momentum;P_{n} [GeV/c];#Delta#beta_{n}", 50, 0, 1.5, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_P_n_goodN_Step1_epFDn);
    TH2D *h_dbeta_n_VS_ToF_goodN_Step1_epFDn = new TH2D("dbeta_n_VS_ToF_goodN_Step1_epFDn", "#Delta#beta_{n} vs Neutron ToF;ToF [ns];#Delta#beta_{n}", 50, 0, 50, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_dbeta_n_VS_P_n_badN_Step1_epFDn = new TH2D("dbeta_n_VS_P_n_badN_Step1_epFDn", "#Delta#beta_{n} vs Neutron Momentum;P_{n} [GeV/c];#Delta#beta_{n}", 50, 0, 1.5, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_P_n_badN_Step1_epFDn);
    TH2D *h_dbeta_n_VS_ToF_badN_Step1_epFDn = new TH2D("dbeta_n_VS_ToF_badN_Step1_epFDn", "#Delta#beta_{n} vs Neutron ToF;ToF [ns];#Delta#beta_{n}", 50, 0, 50, 50, -0.2, 0.2);
    HistoList.push_back(h_dbeta_n_VS_ToF_badN_Step1_epFDn);

    TH1D *h_Vhit_z_n_goodN_Step1_epCDn = new TH1D("Vhit_z_n_goodN_Step1_epCDn", "V_{hit,z}^{n} Distribution;V_{hit,z}^{n} [cm]", 50, -50, 50);
    HistoList.push_back(h_Vhit_z_n_goodN_Step1_epCDn);
    TH1D *h_Vhit_z_n_badN_Step1_epCDn = new TH1D("Vhit_z_n_badN_Step1_epCDn", "V_{hit,z}^{n} Distribution;V_{hit,z}^{n} [cm]", 50, -50, 50);
    HistoList.push_back(h_Vhit_z_n_badN_Step1_epCDn);

    TH1D *h_Vhit_z_n_goodN_Step1_epFDn = new TH1D("Vhit_z_n_goodN_Step1_epFDn", "V_{hit,z}^{n} Distribution;V_{hit,z}^{n} [cm]", 50, -50, 50);
    HistoList.push_back(h_Vhit_z_n_goodN_Step1_epFDn);
    TH1D *h_Vhit_z_n_badN_Step1_epFDn = new TH1D("Vhit_z_n_badN_Step1_epFDn", "V_{hit,z}^{n} Distribution;V_{hit,z}^{n} [cm]", 50, -50, 50);
    HistoList.push_back(h_Vhit_z_n_badN_Step1_epFDn);

    TH1D *h_ToF_n_goodN_Step1_epCDn = new TH1D("ToF_n_goodN_Step1_epCDn", "ToF Distribution;ToF [ns]", 50, 0, 50);
    HistoList.push_back(h_ToF_n_goodN_Step1_epCDn);
    TH1D *h_ToF_n_badN_Step1_epCDn = new TH1D("ToF_n_badN_Step1_epCDn", "ToF Distribution;ToF [ns]", 50, 0, 50);
    HistoList.push_back(h_ToF_n_badN_Step1_epCDn);

    TH1D *h_ToF_n_goodN_Step1_epFDn = new TH1D("ToF_n_goodN_Step1_epFDn", "ToF Distribution;ToF [ns]", 50, 0, 50);
    HistoList.push_back(h_ToF_n_goodN_Step1_epFDn);
    TH1D *h_ToF_n_badN_Step1_epFDn = new TH1D("ToF_n_badN_Step1_epFDn", "ToF Distribution;ToF [ns]", 50, 0, 50);
    HistoList.push_back(h_ToF_n_badN_Step1_epFDn);

    /* Kinematical variables */
    TH1D *h_theta_n_goodN_Step1_epCDn = new TH1D("theta_n_goodN_Step1_epCDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_goodN_Step1_epCDn);
    TH1D *h_theta_n_badN_Step1_epCDn = new TH1D("theta_n_badN_Step1_epCDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_badN_Step1_epCDn);
    TH1D *h_phi_n_goodN_Step1_epCDn = new TH1D("phi_n_goodN_Step1_epCDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_goodN_Step1_epCDn);
    TH1D *h_phi_n_badN_Step1_epCDn = new TH1D("phi_n_badN_Step1_epCDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_phi_n_goodN_Step1_epCDn = new TH2D("theta_n_VS_phi_n_goodN_Step1_epCDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_phi_n_badN_Step1_epCDn = new TH2D("theta_n_VS_phi_n_badN_Step1_epCDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_beta_n_goodN_Step1_epCDn = new TH2D("theta_VS_beta_goodN_Step1_epCDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_beta_n_badN_Step1_epCDn = new TH2D("theta_VS_beta_badN_Step1_epCDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_badN_Step1_epCDn);

    TH1D *h_theta_n_goodN_Step1_epFDn = new TH1D("theta_n_goodN_Step1_epFDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_goodN_Step1_epFDn);
    TH1D *h_theta_n_badN_Step1_epFDn = new TH1D("theta_n_badN_Step1_epFDn", "Neutron Polar Angle Distribution;#theta_{n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_badN_Step1_epFDn);
    TH1D *h_phi_n_goodN_Step1_epFDn = new TH1D("phi_n_goodN_Step1_epFDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_goodN_Step1_epFDn);
    TH1D *h_phi_n_badN_Step1_epFDn = new TH1D("phi_n_badN_Step1_epFDn", "Neutron Azimuthal Angle Distribution;#phi_{n} [#circ]", 48, -180, 180);
    HistoList.push_back(h_phi_n_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_phi_n_goodN_Step1_epFDn = new TH2D("theta_n_VS_phi_n_goodN_Step1_epFDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_phi_n_badN_Step1_epFDn = new TH2D("theta_n_VS_phi_n_badN_Step1_epFDn", "Neutron Angular Distribution;#phi_{n} [#circ];#theta_{n} [#circ]", 48, -180, 180, 50, 0, 180);
    HistoList.push_back(h_theta_n_VS_phi_n_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_beta_n_goodN_Step1_epFDn = new TH2D("theta_VS_beta_goodN_Step1_epFDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_beta_n_badN_Step1_epFDn = new TH2D("theta_VS_beta_badN_Step1_epFDn", "Neutron theta vs beta;#beta;#theta [#circ]", 50, -0.1, 1.1, 55, 0, 180);
    HistoList.push_back(h_theta_n_VS_beta_n_badN_Step1_epFDn);

    TH1D *h_P_n_goodN_Step1_epCDn = new TH1D("P_n_goodN_Step1_epCDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_goodN_Step1_epCDn);
    TH1D *h_P_n_badN_Step1_epCDn = new TH1D("P_n_badN_Step1_epCDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_badN_Step1_epCDn);
    TH2D *h_P_n_VS_theta_n_goodN_Step1_epCDn = new TH2D("P_n_VS_theta_n_goodN_Step1_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_theta_n_badN_Step1_epCDn = new TH2D("P_n_VS_theta_n_badN_Step1_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_badN_Step1_epCDn);

    TH1D *h_P_n_goodN_Step1_epFDn = new TH1D("P_n_goodN_Step1_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_goodN_Step1_epFDn);
    TH1D *h_P_n_badN_Step1_epFDn = new TH1D("P_n_badN_Step1_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_badN_Step1_epFDn);
    TH2D *h_P_n_VS_theta_n_goodN_Step1_epFDn = new TH2D("P_n_VS_theta_n_goodN_Step1_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_theta_n_badN_Step1_epFDn = new TH2D("P_n_VS_theta_n_badN_Step1_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_theta_n_badN_Step1_epFDn);

    TH1D *h_P_miss_goodN_Step1_epCDn = new TH1D("P_miss_goodN_Step1_epCDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_goodN_Step1_epCDn);
    TH1D *h_P_miss_badN_Step1_epCDn = new TH1D("P_miss_badN_Step1_epCDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_theta_miss_goodN_Step1_epCDn = new TH2D("P_miss_VS_theta_miss_goodN_Step1_epCDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_theta_miss_badN_Step1_epCDn = new TH2D("P_miss_VS_theta_miss_badN_Step1_epCDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_phi_miss_goodN_Step1_epCDn = new TH2D("P_miss_VS_phi_miss_goodN_Step1_epCDn", "Missing Momentum vs #phi_{miss};#phi_{miss} [#circ];P_{miss} [GeV/c]", 48, -180, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_phi_miss_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_phi_miss_badN_Step1_epCDn = new TH2D("P_miss_VS_phi_miss_badN_Step1_epCDn", "Missing Momentum vs #phi_{miss};#phi_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, -180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_phi_miss_badN_Step1_epCDn);

    TH1D *h_P_miss_goodN_Step1_epFDn = new TH1D("P_miss_goodN_Step1_epFDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_goodN_Step1_epFDn);
    TH1D *h_P_miss_badN_Step1_epFDn = new TH1D("P_miss_badN_Step1_epFDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_theta_miss_goodN_Step1_epFDn = new TH2D("P_miss_VS_theta_miss_goodN_Step1_epFDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_theta_miss_badN_Step1_epFDn = new TH2D("P_miss_VS_theta_miss_badN_Step1_epFDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_phi_miss_goodN_Step1_epFDn = new TH2D("P_miss_VS_phi_miss_goodN_Step1_epFDn", "Missing Momentum vs #phi_{miss};#phi_{miss} [#circ];P_{miss} [GeV/c]", 48, -180, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_phi_miss_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_phi_miss_badN_Step1_epFDn = new TH2D("P_miss_VS_phi_miss_badN_Step1_epFDn", "Missing Momentum vs #phi_{miss};#phi_{miss} [#circ];P_{miss} [GeV/c]", 48, -180, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_phi_miss_badN_Step1_epFDn);

    TH1D *h_dpp_allN_Step1_epCDn = new TH1D("dpp_allN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_allN_Step1_epCDn);
    TH1D *h_dpp_goodN_Step1_epCDn = new TH1D("dpp_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_goodN_Step1_epCDn);
    TH1D *h_dpp_badN_Step1_epCDn = new TH1D("dpp_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_badN_Step1_epCDn);

    TH1D *h_dpp_allN_Step1_epFDn = new TH1D("dpp_allN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_allN_Step1_epFDn);
    TH1D *h_dpp_goodN_Step1_epFDn = new TH1D("dpp_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_goodN_Step1_epFDn);
    TH1D *h_dpp_badN_Step1_epFDn = new TH1D("dpp_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_badN_Step1_epFDn);

    TH1D *h_theta_n_miss_allN_Step1_epCDn = new TH1D("theta_n_miss_allN_Step1_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_allN_Step1_epCDn);
    TH1D *h_theta_n_miss_goodN_Step1_epCDn = new TH1D("theta_n_miss_goodN_Step1_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_goodN_Step1_epCDn);
    TH1D *h_theta_n_miss_badN_Step1_epCDn = new TH1D("theta_n_miss_badN_Step1_epCDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_badN_Step1_epCDn);

    TH1D *h_theta_n_miss_allN_Step1_epFDn = new TH1D("theta_n_miss_allN_Step1_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_allN_Step1_epFDn);
    TH1D *h_theta_n_miss_goodN_Step1_epFDn = new TH1D("theta_n_miss_goodN_Step1_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_goodN_Step1_epFDn);
    TH1D *h_theta_n_miss_badN_Step1_epFDn = new TH1D("theta_n_miss_badN_Step1_epFDn", "#theta_{n,miss} Distribution;#theta_{n,miss} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_miss_badN_Step1_epFDn);

    TH2D *h_dpp_VS_theta_n_miss_allN_Step1_epCDn = new TH2D("dpp_VS_theta_n_miss_allN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -1.5, 1.5, 50, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_allN_Step1_epCDn);

    TH2D *h_dpp_VS_theta_n_miss_allN_Step1_epFDn = new TH2D("dpp_VS_theta_n_miss_allN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -1.5, 1.5, 50, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_allN_Step1_epFDn);

    TH1D *h_E_p_goodN_Step1_epCDn = new TH1D("E_p_goodN_Step1_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_goodN_Step1_epCDn);
    TH1D *h_E_p_badN_Step1_epCDn = new TH1D("E_p_badN_Step1_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_badN_Step1_epCDn);
    TH1D *h_E_miss_goodN_Step1_epCDn = new TH1D("E_miss_goodN_Step1_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_goodN_Step1_epCDn);
    TH1D *h_E_miss_badN_Step1_epCDn = new TH1D("E_miss_badN_Step1_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_badN_Step1_epCDn);
    TH1D *h_M_miss_goodN_Step1_epCDn = new TH1D("M_miss_goodN_Step1_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_goodN_Step1_epCDn);
    TH1D *h_M_miss_badN_Step1_epCDn = new TH1D("M_miss_badN_Step1_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_P_n_goodN_Step1_epCDn = new TH2D("M_miss_VS_P_n_goodN_Step1_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_P_n_badN_Step1_epCDn = new TH2D("M_miss_VS_P_n_badN_Step1_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_theta_n_goodN_Step1_epCDn = new TH2D("M_miss_VS_theta_n_goodN_Step1_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_theta_n_badN_Step1_epCDn = new TH2D("M_miss_VS_theta_n_badN_Step1_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_phi_n_goodN_Step1_epCDn = new TH2D("M_miss_VS_phi_n_goodN_Step1_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_phi_n_badN_Step1_epCDn = new TH2D("M_miss_VS_phi_n_badN_Step1_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_P_miss_goodN_Step1_epCDn = new TH2D("M_miss_VS_P_miss_goodN_Step1_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_P_miss_badN_Step1_epCDn = new TH2D("M_miss_VS_P_miss_badN_Step1_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_theta_miss_goodN_Step1_epCDn = new TH2D("M_miss_VS_theta_miss_goodN_Step1_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_theta_miss_badN_Step1_epCDn = new TH2D("M_miss_VS_theta_miss_badN_Step1_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_phi_miss_goodN_Step1_epCDn = new TH2D("M_miss_VS_phi_miss_goodN_Step1_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_phi_miss_badN_Step1_epCDn = new TH2D("M_miss_VS_phi_miss_badN_Step1_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_badN_Step1_epCDn);

    TH1D *h_E_p_goodN_Step1_epFDn = new TH1D("E_p_goodN_Step1_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_goodN_Step1_epFDn);
    TH1D *h_E_p_badN_Step1_epFDn = new TH1D("E_p_badN_Step1_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_badN_Step1_epFDn);
    TH1D *h_E_miss_goodN_Step1_epFDn = new TH1D("E_miss_goodN_Step1_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_goodN_Step1_epFDn);
    TH1D *h_E_miss_badN_Step1_epFDn = new TH1D("E_miss_badN_Step1_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_badN_Step1_epFDn);
    TH1D *h_M_miss_goodN_Step1_epFDn = new TH1D("M_miss_goodN_Step1_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_goodN_Step1_epFDn);
    TH1D *h_M_miss_badN_Step1_epFDn = new TH1D("M_miss_badN_Step1_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_P_n_goodN_Step1_epFDn = new TH2D("M_miss_VS_P_n_goodN_Step1_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_P_n_badN_Step1_epFDn = new TH2D("M_miss_VS_P_n_badN_Step1_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_n_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_theta_n_goodN_Step1_epFDn = new TH2D("M_miss_VS_theta_n_goodN_Step1_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_theta_n_badN_Step1_epFDn = new TH2D("M_miss_VS_theta_n_badN_Step1_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_n_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_phi_n_goodN_Step1_epFDn = new TH2D("M_miss_VS_phi_n_goodN_Step1_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_phi_n_badN_Step1_epFDn = new TH2D("M_miss_VS_phi_n_badN_Step1_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_n_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_P_miss_goodN_Step1_epFDn = new TH2D("M_miss_VS_P_miss_goodN_Step1_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_P_miss_badN_Step1_epFDn = new TH2D("M_miss_VS_P_miss_badN_Step1_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.5, 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_P_miss_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_theta_miss_goodN_Step1_epFDn = new TH2D("M_miss_VS_theta_miss_goodN_Step1_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_theta_miss_badN_Step1_epFDn = new TH2D("M_miss_VS_theta_miss_badN_Step1_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_theta_miss_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_phi_miss_goodN_Step1_epFDn = new TH2D("M_miss_VS_phi_miss_goodN_Step1_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_phi_miss_badN_Step1_epFDn = new TH2D("M_miss_VS_phi_miss_badN_Step1_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.7);
    HistoList.push_back(h_M_miss_VS_phi_miss_badN_Step1_epFDn);

    TH1D *h_P_n_minus_P_miss_goodN_Step1_epCDn = new TH1D("P_n_minus_P_miss_goodN_Step1_epCDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_goodN_Step1_epCDn);
    TH1D *h_P_n_minus_P_miss_badN_Step1_epCDn = new TH1D("P_n_minus_P_miss_badN_Step1_epCDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_badN_Step1_epCDn);
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step1_epCDn = new TH1D("P_n_x_minus_P_miss_x_goodN_Step1_epCDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_goodN_Step1_epCDn);
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step1_epCDn = new TH1D("P_n_x_minus_P_miss_x_badN_Step1_epCDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_badN_Step1_epCDn);
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step1_epCDn = new TH1D("P_n_y_minus_P_miss_y_goodN_Step1_epCDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_goodN_Step1_epCDn);
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step1_epCDn = new TH1D("P_n_y_minus_P_miss_y_badN_Step1_epCDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_badN_Step1_epCDn);
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step1_epCDn = new TH1D("P_n_z_minus_P_miss_z_goodN_Step1_epCDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_goodN_Step1_epCDn);
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step1_epCDn = new TH1D("P_n_z_minus_P_miss_z_badN_Step1_epCDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_badN_Step1_epCDn);

    TH1D *h_P_n_minus_P_miss_goodN_Step1_epFDn = new TH1D("P_n_minus_P_miss_goodN_Step1_epFDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_goodN_Step1_epFDn);
    TH1D *h_P_n_minus_P_miss_badN_Step1_epFDn = new TH1D("P_n_minus_P_miss_badN_Step1_epFDn", "P_{n}-P_{miss} Distribution;P_{n}-P_{miss} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_minus_P_miss_badN_Step1_epFDn);
    TH1D *h_P_n_x_minus_P_miss_x_goodN_Step1_epFDn = new TH1D("P_n_x_minus_P_miss_x_goodN_Step1_epFDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_goodN_Step1_epFDn);
    TH1D *h_P_n_x_minus_P_miss_x_badN_Step1_epFDn = new TH1D("P_n_x_minus_P_miss_x_badN_Step1_epFDn", "P_{n,x}-P_{miss,x} Distribution; P_{n,x}-P_{miss,x} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_x_minus_P_miss_x_badN_Step1_epFDn);
    TH1D *h_P_n_y_minus_P_miss_y_goodN_Step1_epFDn = new TH1D("P_n_y_minus_P_miss_y_goodN_Step1_epFDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_goodN_Step1_epFDn);
    TH1D *h_P_n_y_minus_P_miss_y_badN_Step1_epFDn = new TH1D("P_n_y_minus_P_miss_y_badN_Step1_epFDn", "P_{n,y}-P_{miss,y} Distribution; P_{n,y}-P_{miss,y} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_y_minus_P_miss_y_badN_Step1_epFDn);
    TH1D *h_P_n_z_minus_P_miss_z_goodN_Step1_epFDn = new TH1D("P_n_z_minus_P_miss_z_goodN_Step1_epFDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_goodN_Step1_epFDn);
    TH1D *h_P_n_z_minus_P_miss_z_badN_Step1_epFDn = new TH1D("P_n_z_minus_P_miss_z_badN_Step1_epFDn", "P_{n,z}-P_{miss,z} Distribution; P_{n,z}-P_{miss,z} [GeV/c];Counts", 50, -1.5, 1.5);
    HistoList.push_back(h_P_n_z_minus_P_miss_z_badN_Step1_epFDn);

    TH2D *h_P_n_VS_P_miss_goodN_Step1_epCDn = new TH2D("P_n_VS_P_miss_goodN_Step1_epCDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_P_miss_badN_Step1_epCDn = new TH2D("P_n_VS_P_miss_badN_Step1_epCDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_badN_Step1_epCDn);
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step1_epCDn = new TH2D("P_n_x_VS_P_miss_x_goodN_Step1_epCDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_goodN_Step1_epCDn);
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step1_epCDn = new TH2D("P_n_x_VS_P_miss_x_badN_Step1_epCDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_badN_Step1_epCDn);
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step1_epCDn = new TH2D("P_n_y_VS_P_miss_y_goodN_Step1_epCDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_goodN_Step1_epCDn);
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step1_epCDn = new TH2D("P_n_y_VS_P_miss_y_badN_Step1_epCDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_badN_Step1_epCDn);
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step1_epCDn = new TH2D("P_n_z_VS_P_miss_z_goodN_Step1_epCDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_goodN_Step1_epCDn);
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step1_epCDn = new TH2D("P_n_z_VS_P_miss_z_badN_Step1_epCDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_badN_Step1_epCDn);

    TH2D *h_P_n_VS_P_miss_goodN_Step1_epFDn = new TH2D("P_n_VS_P_miss_goodN_Step1_epFDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_P_miss_badN_Step1_epFDn = new TH2D("P_n_VS_P_miss_badN_Step1_epFDn", "P_{n} vs P_{miss} Distribution;P_{n} [GeV/c];P_{miss} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_P_miss_badN_Step1_epFDn);
    TH2D *h_P_n_x_VS_P_miss_x_goodN_Step1_epFDn = new TH2D("P_n_x_VS_P_miss_x_goodN_Step1_epFDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_goodN_Step1_epFDn);
    TH2D *h_P_n_x_VS_P_miss_x_badN_Step1_epFDn = new TH2D("P_n_x_VS_P_miss_x_badN_Step1_epFDn", "P_{n,x} vs P_{miss,x} Distribution;P_{n,x} [GeV/c];P_{miss,x} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_x_VS_P_miss_x_badN_Step1_epFDn);
    TH2D *h_P_n_y_VS_P_miss_y_goodN_Step1_epFDn = new TH2D("P_n_y_VS_P_miss_y_goodN_Step1_epFDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_goodN_Step1_epFDn);
    TH2D *h_P_n_y_VS_P_miss_y_badN_Step1_epFDn = new TH2D("P_n_y_VS_P_miss_y_badN_Step1_epFDn", "P_{n,y} vs P_{miss,y} Distribution;P_{n,y} [GeV/c];P_{miss,y} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_y_VS_P_miss_y_badN_Step1_epFDn);
    TH2D *h_P_n_z_VS_P_miss_z_goodN_Step1_epFDn = new TH2D("P_n_z_VS_P_miss_z_goodN_Step1_epFDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_goodN_Step1_epFDn);
    TH2D *h_P_n_z_VS_P_miss_z_badN_Step1_epFDn = new TH2D("P_n_z_VS_P_miss_z_badN_Step1_epFDn", "P_{n,z} vs P_{miss,z} Distribution;P_{n,z} [GeV/c];P_{miss,z} [GeV/c]", 50, 0, 1.5, 50, 0, 1.5);
    HistoList.push_back(h_P_n_z_VS_P_miss_z_badN_Step1_epFDn);

    TH1D *h_theta_n_p_goodN_Step1_epCDn = new TH1D("theta_n_p_goodN_Step1_epCDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_goodN_Step1_epCDn);
    TH1D *h_theta_n_p_badN_Step1_epCDn = new TH1D("theta_n_p_badN_Step1_epCDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_badN_Step1_epCDn);
    TH2D *h_theta_n_p_VS_P_p_goodN_Step1_epCDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_goodN_Step1_epCDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_goodN_Step1_epCDn);
    TH2D *h_theta_n_p_VS_P_p_badN_Step1_epCDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_badN_Step1_epCDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_badN_Step1_epCDn);

    TH1D *h_theta_n_p_goodN_Step1_epFDn = new TH1D("theta_n_p_goodN_Step1_epFDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_goodN_Step1_epFDn);
    TH1D *h_theta_n_p_badN_Step1_epFDn = new TH1D("theta_n_p_badN_Step1_epFDn", "#theta_{p,n} Distribution;#theta_{p,n} [#circ]", 50, 0, 180);
    HistoList.push_back(h_theta_n_p_badN_Step1_epFDn);
    TH2D *h_theta_n_p_VS_P_p_goodN_Step1_epFDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_goodN_Step1_epFDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_goodN_Step1_epFDn);
    TH2D *h_theta_n_p_VS_P_p_badN_Step1_epFDn = new TH2D("theta_n_p_VS_P_p_VS_P_p_badN_Step1_epFDn", "#theta_{p,n} vs P_{p};P_{p} [GeV/c];#theta_{p,n} [#circ]", 50, 0, 1.5, 50, 0, 180);
    HistoList.push_back(h_theta_n_p_VS_P_p_badN_Step1_epFDn);

    TH1D *h_xB_goodN_Step1_epCDn = new TH1D("xB_goodN_Step1_epCDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_goodN_Step1_epCDn);
    TH1D *h_xB_badN_Step1_epCDn = new TH1D("xB_badN_Step1_epCDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_badN_Step1_epCDn);

    TH1D *h_xB_goodN_Step1_epFDn = new TH1D("xB_goodN_Step1_epFDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_goodN_Step1_epFDn);
    TH1D *h_xB_badN_Step1_epFDn = new TH1D("xB_badN_Step1_epFDn", "x_{B} Distribution;x_{B};Counts", 50, 0, 2);
    HistoList.push_back(h_xB_badN_Step1_epFDn);

    /* Detector responses */
    TH1D *h_Edep_CND_goodN_Step1_epCDn = new TH1D("Edep_CND_goodN_Step1_epCDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_goodN_Step1_epCDn);
    TH1D *h_Edep_CND_badN_Step1_epCDn = new TH1D("Edep_CND_badN_Step1_epCDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_badN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("P_n_VS_Edep_CND_goodN_Step1_epCDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CND_badN_Step1_epCDn = new TH2D("P_n_VS_Edep_CND_badN_Step1_epCDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CND_goodN_Step1_epCDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND_badN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CND_badN_Step1_epCDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CND_goodN_Step1_epCDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND_badN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CND_badN_Step1_epCDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CND_goodN_Step1_epCDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND_badN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CND_badN_Step1_epCDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CND_goodN_Step1_epCDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CND_badN_Step1_epCDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CND_goodN_Step1_epCDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CND_badN_Step1_epCDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("dpp_VS_Edep_CND_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CND_badN_Step1_epCDn = new TH2D("dpp_VS_Edep_CND_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CND_goodN_Step1_epCDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND_badN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CND_badN_Step1_epCDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("E_p_VS_Edep_CND_goodN_Step1_epCDn", "E_{p} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CND_badN_Step1_epCDn = new TH2D("E_p_VS_Edep_CND_badN_Step1_epCDn", "E_{p} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CND_goodN_Step1_epCDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND_badN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CND_badN_Step1_epCDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND_goodN_Step1_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND_badN_Step1_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("path_VS_Edep_CND_goodN_Step1_epCDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CND_badN_Step1_epCDn = new TH2D("path_VS_Edep_CND_badN_Step1_epCDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND_goodN_Step1_epCDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND_badN_Step1_epCDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("ToF_VS_Edep_CND_goodN_Step1_epCDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CND_badN_Step1_epCDn = new TH2D("ToF_VS_Edep_CND_badN_Step1_epCDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("nSector_VS_Edep_CND_goodN_Step1_epCDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CND_badN_Step1_epCDn = new TH2D("nSector_VS_Edep_CND_badN_Step1_epCDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("Edep_CND1_VS_Edep_CND_goodN_Step1_epCDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step1_epCDn = new TH2D("Edep_CND1_VS_Edep_CND_badN_Step1_epCDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("Edep_CND2_VS_Edep_CND_goodN_Step1_epCDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step1_epCDn = new TH2D("Edep_CND2_VS_Edep_CND_badN_Step1_epCDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_badN_Step1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("Edep_CND3_VS_Edep_CND_goodN_Step1_epCDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step1_epCDn = new TH2D("Edep_CND3_VS_Edep_CND_badN_Step1_epCDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_badN_Step1_epCDn);

    TH1D *h_Edep_CND_goodN_Step1_epFDn = new TH1D("Edep_CND_goodN_Step1_epFDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_goodN_Step1_epFDn);
    TH1D *h_Edep_CND_badN_Step1_epFDn = new TH1D("Edep_CND_badN_Step1_epFDn", "Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND_badN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("P_n_VS_Edep_CND_goodN_Step1_epFDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CND_badN_Step1_epFDn = new TH2D("P_n_VS_Edep_CND_badN_Step1_epFDn", "Neutron Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CND_goodN_Step1_epFDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND_badN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CND_badN_Step1_epFDn", "#theta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CND_goodN_Step1_epFDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND_badN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CND_badN_Step1_epFDn", "#phi_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CND_goodN_Step1_epFDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND_badN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CND_badN_Step1_epFDn", "Missing Momentum vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CND_goodN_Step1_epFDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND_badN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CND_badN_Step1_epFDn", "#theta_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CND_goodN_Step1_epFDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND_badN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CND_badN_Step1_epFDn", "#phi_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("dpp_VS_Edep_CND_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CND_badN_Step1_epFDn = new TH2D("dpp_VS_Edep_CND_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CND_goodN_Step1_epFDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND_badN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CND_badN_Step1_epFDn", "#beta_{n} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("E_p_VS_Edep_CND_goodN_Step1_epFDn", "E_{p} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CND_badN_Step1_epFDn = new TH2D("E_p_VS_Edep_CND_badN_Step1_epFDn", "E_{p} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CND_goodN_Step1_epFDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND_badN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CND_badN_Step1_epFDn", "E_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND_goodN_Step1_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND_badN_Step1_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("path_VS_Edep_CND_goodN_Step1_epFDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CND_badN_Step1_epFDn = new TH2D("path_VS_Edep_CND_badN_Step1_epFDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND_goodN_Step1_epFDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_badN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND_badN_Step1_epFDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("ToF_VS_Edep_CND_goodN_Step1_epFDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CND_badN_Step1_epFDn = new TH2D("ToF_VS_Edep_CND_badN_Step1_epFDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("nSector_VS_Edep_CND_goodN_Step1_epFDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CND_badN_Step1_epFDn = new TH2D("nSector_VS_Edep_CND_badN_Step1_epFDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("Edep_CND1_VS_Edep_CND_goodN_Step1_epFDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CND_badN_Step1_epFDn = new TH2D("Edep_CND1_VS_Edep_CND_badN_Step1_epFDn", "E^{CND,1}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("Edep_CND2_VS_Edep_CND_goodN_Step1_epFDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND_badN_Step1_epFDn = new TH2D("Edep_CND2_VS_Edep_CND_badN_Step1_epFDn", "E^{CND,2}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND_badN_Step1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("Edep_CND3_VS_Edep_CND_goodN_Step1_epFDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND_badN_Step1_epFDn = new TH2D("Edep_CND3_VS_Edep_CND_badN_Step1_epFDn", "E^{CND,3}_{dep} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND_badN_Step1_epFDn);

    TH1D *h_Edep_CTOF_goodN_Step1_epCDn = new TH1D("Edep_CTOF_goodN_Step1_epCDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_goodN_Step1_epCDn);
    TH1D *h_Edep_CTOF_badN_Step1_epCDn = new TH1D("Edep_CTOF_badN_Step1_epCDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("P_n_VS_Edep_CTOF_goodN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("P_n_VS_Edep_CTOF_badN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CTOF_goodN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CTOF_badN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CTOF_goodN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CTOF_badN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CTOF_goodN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CTOF_badN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CTOF_goodN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CTOF_badN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CTOF_goodN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CTOF_badN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("dpp_VS_Edep_CTOF_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("dpp_VS_Edep_CTOF_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CTOF_goodN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CTOF_badN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("E_p_VS_Edep_CTOF_goodN_Step1_epCDn", "E_{p} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("E_p_VS_Edep_CTOF_badN_Step1_epCDn", "E_{p} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CTOF_goodN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CTOF_badN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CTOF_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CTOF_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("path_VS_Edep_CTOF_goodN_Step1_epCDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("path_VS_Edep_CTOF_badN_Step1_epCDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CTOF_goodN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CTOF_badN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("ToF_VS_Edep_CTOF_goodN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("ToF_VS_Edep_CTOF_badN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("nSector_VS_Edep_CTOF_goodN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("nSector_VS_Edep_CTOF_badN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("Edep_CND1_VS_Edep_CTOF_goodN_Step1_epCDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("Edep_CND1_VS_Edep_CTOF_badN_Step1_epCDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("Edep_CND2_VS_Edep_CTOF_goodN_Step1_epCDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("Edep_CND2_VS_Edep_CTOF_badN_Step1_epCDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_badN_Step1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("Edep_CND3_VS_Edep_CTOF_goodN_Step1_epCDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("Edep_CND3_VS_Edep_CTOF_badN_Step1_epCDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_badN_Step1_epCDn);

    TH1D *h_Edep_CTOF_goodN_Step1_epFDn = new TH1D("Edep_CTOF_goodN_Step1_epFDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_goodN_Step1_epFDn);
    TH1D *h_Edep_CTOF_badN_Step1_epFDn = new TH1D("Edep_CTOF_badN_Step1_epFDn", "Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("P_n_VS_Edep_CTOF_goodN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("P_n_VS_Edep_CTOF_badN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CTOF_goodN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CTOF_badN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CTOF_goodN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CTOF_badN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CTOF_goodN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CTOF_badN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CTOF_goodN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CTOF_badN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CTOF_goodN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CTOF_badN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("dpp_VS_Edep_CTOF_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("dpp_VS_Edep_CTOF_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CTOF_goodN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CTOF_badN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("E_p_VS_Edep_CTOF_goodN_Step1_epFDn", "E_{p} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("E_p_VS_Edep_CTOF_badN_Step1_epFDn", "E_{p} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CTOF_goodN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CTOF_badN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CTOF_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CTOF_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("path_VS_Edep_CTOF_goodN_Step1_epFDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("path_VS_Edep_CTOF_badN_Step1_epFDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CTOF_goodN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CTOF_badN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("ToF_VS_Edep_CTOF_goodN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("ToF_VS_Edep_CTOF_badN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("nSector_VS_Edep_CTOF_goodN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("nSector_VS_Edep_CTOF_badN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("Edep_CND1_VS_Edep_CTOF_goodN_Step1_epFDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_Edep_CND1_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("Edep_CND1_VS_Edep_CTOF_badN_Step1_epFDn", "E^{CND,1}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,1}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("Edep_CND2_VS_Edep_CTOF_goodN_Step1_epFDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("Edep_CND2_VS_Edep_CTOF_badN_Step1_epFDn", "E^{CND,2}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CTOF_badN_Step1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("Edep_CND3_VS_Edep_CTOF_goodN_Step1_epFDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("Edep_CND3_VS_Edep_CTOF_badN_Step1_epFDn", "E^{CND,3}_{dep} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CTOF_badN_Step1_epFDn);

    TH1D *h_Edep_single_goodN_Step1_epCDn = new TH1D("Edep_single_goodN_Step1_epCDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_goodN_Step1_epCDn);
    TH1D *h_Edep_single_badN_Step1_epCDn = new TH1D("Edep_single_badN_Step1_epCDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_badN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_single_goodN_Step1_epCDn = new TH2D("P_n_VS_Edep_single_goodN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_single_badN_Step1_epCDn = new TH2D("P_n_VS_Edep_single_badN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_single_goodN_Step1_epCDn = new TH2D("theta_n_VS_Edep_single_goodN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_single_badN_Step1_epCDn = new TH2D("theta_n_VS_Edep_single_badN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_single_goodN_Step1_epCDn = new TH2D("phi_n_VS_Edep_single_goodN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_single_badN_Step1_epCDn = new TH2D("phi_n_VS_Edep_single_badN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_single_goodN_Step1_epCDn = new TH2D("P_miss_VS_Edep_single_goodN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_single_badN_Step1_epCDn = new TH2D("P_miss_VS_Edep_single_badN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_single_goodN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_single_badN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_single_badN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_single_goodN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_single_badN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_single_badN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_single_goodN_Step1_epCDn = new TH2D("dpp_VS_Edep_single_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_single_badN_Step1_epCDn = new TH2D("dpp_VS_Edep_single_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_single_goodN_Step1_epCDn = new TH2D("beta_n_VS_Edep_single_goodN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_single_badN_Step1_epCDn = new TH2D("beta_n_VS_Edep_single_badN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_single_goodN_Step1_epCDn = new TH2D("E_p_VS_Edep_single_goodN_Step1_epCDn", "E_{p} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_single_badN_Step1_epCDn = new TH2D("E_p_VS_Edep_single_badN_Step1_epCDn", "E_{p} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_single_goodN_Step1_epCDn = new TH2D("E_miss_VS_Edep_single_goodN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_single_badN_Step1_epCDn = new TH2D("E_miss_VS_Edep_single_badN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_single_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_single_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_single_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_single_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_path_VS_Edep_single_goodN_Step1_epCDn = new TH2D("path_VS_Edep_single_goodN_Step1_epCDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_path_VS_Edep_single_badN_Step1_epCDn = new TH2D("path_VS_Edep_single_badN_Step1_epCDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_single_goodN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_single_badN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_single_goodN_Step1_epCDn = new TH2D("ToF_VS_Edep_single_goodN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_single_badN_Step1_epCDn = new TH2D("ToF_VS_Edep_single_badN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_badN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_single_goodN_Step1_epCDn = new TH2D("nSector_VS_Edep_single_goodN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_single_badN_Step1_epCDn = new TH2D("nSector_VS_Edep_single_badN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_badN_Step1_epCDn);

    TH1D *h_Edep_single_goodN_Step1_epFDn = new TH1D("Edep_single_goodN_Step1_epFDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_goodN_Step1_epFDn);
    TH1D *h_Edep_single_badN_Step1_epFDn = new TH1D("Edep_single_badN_Step1_epFDn", "Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_single_badN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_single_goodN_Step1_epFDn = new TH2D("P_n_VS_Edep_single_goodN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_single_badN_Step1_epFDn = new TH2D("P_n_VS_Edep_single_badN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_single_goodN_Step1_epFDn = new TH2D("theta_n_VS_Edep_single_goodN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_single_badN_Step1_epFDn = new TH2D("theta_n_VS_Edep_single_badN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_single_goodN_Step1_epFDn = new TH2D("phi_n_VS_Edep_single_goodN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_single_badN_Step1_epFDn = new TH2D("phi_n_VS_Edep_single_badN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_single_goodN_Step1_epFDn = new TH2D("P_miss_VS_Edep_single_goodN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_single_badN_Step1_epFDn = new TH2D("P_miss_VS_Edep_single_badN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_single_goodN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_single_goodN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_single_badN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_single_badN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_single_goodN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_single_goodN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_single_badN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_single_badN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_single_goodN_Step1_epFDn = new TH2D("dpp_VS_Edep_single_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_single_badN_Step1_epFDn = new TH2D("dpp_VS_Edep_single_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_single_goodN_Step1_epFDn = new TH2D("beta_n_VS_Edep_single_goodN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_single_badN_Step1_epFDn = new TH2D("beta_n_VS_Edep_single_badN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_single_goodN_Step1_epFDn = new TH2D("E_p_VS_Edep_single_goodN_Step1_epFDn", "E_{P} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_single_badN_Step1_epFDn = new TH2D("E_p_VS_Edep_single_badN_Step1_epFDn", "E_{P} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_single_goodN_Step1_epFDn = new TH2D("E_miss_VS_Edep_single_goodN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_single_badN_Step1_epFDn = new TH2D("E_miss_VS_Edep_single_badN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_single_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_single_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_single_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_single_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_path_VS_Edep_single_goodN_Step1_epFDn = new TH2D("path_VS_Edep_single_goodN_Step1_epFDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_path_VS_Edep_single_badN_Step1_epFDn = new TH2D("path_VS_Edep_single_badN_Step1_epFDn", "Path length vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_single_goodN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_single_goodN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_single_badN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_single_badN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_single_goodN_Step1_epFDn = new TH2D("ToF_VS_Edep_single_goodN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_single_badN_Step1_epFDn = new TH2D("ToF_VS_Edep_single_badN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_single_badN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_single_goodN_Step1_epFDn = new TH2D("nSector_VS_Edep_single_goodN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_single_badN_Step1_epFDn = new TH2D("nSector_VS_Edep_single_badN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_single_badN_Step1_epFDn);

    TH1D *h_Edep_CND1_goodN_Step1_epCDn = new TH1D("Edep_CND1_goodN_Step1_epCDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_goodN_Step1_epCDn);
    TH1D *h_Edep_CND1_badN_Step1_epCDn = new TH1D("Edep_CND1_badN_Step1_epCDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("P_n_VS_Edep_CND1_goodN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("P_n_VS_Edep_CND1_badN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CND1_goodN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CND1_badN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CND1_goodN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CND1_badN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CND1_goodN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CND1_badN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CND1_goodN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CND1_badN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CND1_goodN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CND1_badN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("dpp_VS_Edep_CND1_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("dpp_VS_Edep_CND1_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CND1_goodN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CND1_badN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("E_p_VS_Edep_CND1_goodN_Step1_epCDn", "E_{P} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("E_p_VS_Edep_CND1_badN_Step1_epCDn", "E_{P} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CND1_goodN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CND1_badN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND1_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND1_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("path_VS_Edep_CND1_goodN_Step1_epCDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("path_VS_Edep_CND1_badN_Step1_epCDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND1_goodN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND1_badN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("ToF_VS_Edep_CND1_goodN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("ToF_VS_Edep_CND1_badN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("nSector_VS_Edep_CND1_goodN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("nSector_VS_Edep_CND1_badN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("Edep_CND2_VS_Edep_CND1_goodN_Step1_epCDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("Edep_CND2_VS_Edep_CND1_badN_Step1_epCDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_badN_Step1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("Edep_CND3_VS_Edep_CND1_goodN_Step1_epCDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("Edep_CND3_VS_Edep_CND1_badN_Step1_epCDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_badN_Step1_epCDn);

    TH1D *h_Edep_CND1_goodN_Step1_epFDn = new TH1D("Edep_CND1_goodN_Step1_epFDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_goodN_Step1_epFDn);
    TH1D *h_Edep_CND1_badN_Step1_epFDn = new TH1D("Edep_CND1_badN_Step1_epFDn", "Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("P_n_VS_Edep_CND1_goodN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("P_n_VS_Edep_CND1_badN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CND1_goodN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CND1_badN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CND1_goodN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CND1_badN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CND1_goodN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CND1_badN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CND1_goodN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CND1_badN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CND1_goodN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CND1_badN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("dpp_VS_Edep_CND1_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("dpp_VS_Edep_CND1_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CND1_goodN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CND1_badN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("E_p_VS_Edep_CND1_goodN_Step1_epFDn", "E_{P} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("E_p_VS_Edep_CND1_badN_Step1_epFDn", "E_{P} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CND1_goodN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CND1_badN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND1_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND1_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("path_VS_Edep_CND1_goodN_Step1_epFDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("path_VS_Edep_CND1_badN_Step1_epFDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND1_goodN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND1_badN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("ToF_VS_Edep_CND1_goodN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("ToF_VS_Edep_CND1_badN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("nSector_VS_Edep_CND1_goodN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("nSector_VS_Edep_CND1_badN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("Edep_CND2_VS_Edep_CND1_goodN_Step1_epFDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_Edep_CND2_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("Edep_CND2_VS_Edep_CND1_badN_Step1_epFDn", "E^{CND,2}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,2}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_VS_Edep_CND1_badN_Step1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("Edep_CND3_VS_Edep_CND1_goodN_Step1_epFDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("Edep_CND3_VS_Edep_CND1_badN_Step1_epFDn", "E^{CND,3}_{dep} vs E^{CND,1}_{dep};E^{CND,1}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND1_badN_Step1_epFDn);

    TH1D *h_Edep_CND2_goodN_Step1_epCDn = new TH1D("Edep_CND2_goodN_Step1_epCDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_goodN_Step1_epCDn);
    TH1D *h_Edep_CND2_badN_Step1_epCDn = new TH1D("Edep_CND2_badN_Step1_epCDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("P_n_VS_Edep_CND2_goodN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("P_n_VS_Edep_CND2_badN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CND2_goodN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CND2_badN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CND2_goodN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CND2_badN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CND2_goodN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CND2_badN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CND2_goodN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CND2_badN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CND2_goodN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CND2_badN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("dpp_VS_Edep_CND2_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("dpp_VS_Edep_CND2_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CND2_goodN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CND2_badN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("E_p_VS_Edep_CND2_goodN_Step1_epCDn", "E_{P} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("E_p_VS_Edep_CND2_badN_Step1_epCDn", "E_{P} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CND2_goodN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CND2_badN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND2_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND2_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("path_VS_Edep_CND2_goodN_Step1_epCDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("path_VS_Edep_CND2_badN_Step1_epCDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND2_goodN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND2_badN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("ToF_VS_Edep_CND2_goodN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("ToF_VS_Edep_CND2_badN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("nSector_VS_Edep_CND2_goodN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("nSector_VS_Edep_CND2_badN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_badN_Step1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("Edep_CND3_VS_Edep_CND2_goodN_Step1_epCDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("Edep_CND3_VS_Edep_CND2_badN_Step1_epCDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_badN_Step1_epCDn);

    TH1D *h_Edep_CND2_goodN_Step1_epFDn = new TH1D("Edep_CND2_goodN_Step1_epFDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_goodN_Step1_epFDn);
    TH1D *h_Edep_CND2_badN_Step1_epFDn = new TH1D("Edep_CND2_badN_Step1_epFDn", "Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("P_n_VS_Edep_CND2_goodN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("P_n_VS_Edep_CND2_badN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CND2_goodN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CND2_badN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CND2_goodN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CND2_badN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CND2_goodN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CND2_badN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CND2_goodN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CND2_badN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CND2_goodN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CND2_badN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("dpp_VS_Edep_CND2_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("dpp_VS_Edep_CND2_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CND2_goodN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CND2_badN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("E_p_VS_Edep_CND2_goodN_Step1_epFDn", "E_{P} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("E_p_VS_Edep_CND2_badN_Step1_epFDn", "E_{P} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{P} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CND2_goodN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CND2_badN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND2_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND2_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("path_VS_Edep_CND2_goodN_Step1_epFDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("path_VS_Edep_CND2_badN_Step1_epFDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND2_goodN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND2_badN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("ToF_VS_Edep_CND2_goodN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("ToF_VS_Edep_CND2_badN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("nSector_VS_Edep_CND2_goodN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("nSector_VS_Edep_CND2_badN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_badN_Step1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("Edep_CND3_VS_Edep_CND2_goodN_Step1_epFDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_Edep_CND3_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("Edep_CND3_VS_Edep_CND2_badN_Step1_epFDn", "E^{CND,3}_{dep} vs E^{CND,2}_{dep};E^{CND,2}_{dep} [MeV];E^{CND,3}_{dep}", 50, 0, 100, 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_VS_Edep_CND2_badN_Step1_epFDn);

    TH1D *h_Edep_CND3_goodN_Step1_epCDn = new TH1D("Edep_CND3_goodN_Step1_epCDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_goodN_Step1_epCDn);
    TH1D *h_Edep_CND3_badN_Step1_epCDn = new TH1D("Edep_CND3_badN_Step1_epCDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("P_n_VS_Edep_CND3_goodN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("P_n_VS_Edep_CND3_badN_Step1_epCDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CND3_goodN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("theta_n_VS_Edep_CND3_badN_Step1_epCDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CND3_goodN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("phi_n_VS_Edep_CND3_badN_Step1_epCDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CND3_goodN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("P_miss_VS_Edep_CND3_badN_Step1_epCDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CND3_goodN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("theta_miss_VS_Edep_CND3_badN_Step1_epCDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CND3_goodN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("phi_miss_VS_Edep_CND3_badN_Step1_epCDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("dpp_VS_Edep_CND3_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_dpp_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("dpp_VS_Edep_CND3_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CND3_goodN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("beta_n_VS_Edep_CND3_badN_Step1_epCDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("E_p_VS_Edep_CND3_goodN_Step1_epCDn", "E_{p} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_E_p_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("E_p_VS_Edep_CND3_badN_Step1_epCDn", "E_{p} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CND3_goodN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("E_miss_VS_Edep_CND3_badN_Step1_epCDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND3_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND3_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("path_VS_Edep_CND3_goodN_Step1_epCDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_path_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("path_VS_Edep_CND3_badN_Step1_epCDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND3_goodN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND3_badN_Step1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("ToF_VS_Edep_CND3_goodN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_ToF_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("ToF_VS_Edep_CND3_badN_Step1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_badN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("nSector_VS_Edep_CND3_goodN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_nSector_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("nSector_VS_Edep_CND3_badN_Step1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_badN_Step1_epCDn);

    TH1D *h_Edep_CND3_goodN_Step1_epFDn = new TH1D("Edep_CND3_goodN_Step1_epFDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_goodN_Step1_epFDn);
    TH1D *h_Edep_CND3_badN_Step1_epFDn = new TH1D("Edep_CND3_badN_Step1_epFDn", "Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("P_n_VS_Edep_CND3_goodN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("P_n_VS_Edep_CND3_badN_Step1_epFDn", "Neutron Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{n} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CND3_goodN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("theta_n_VS_Edep_CND3_badN_Step1_epFDn", "#theta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CND3_goodN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_phi_n_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("phi_n_VS_Edep_CND3_badN_Step1_epFDn", "#phi_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{n} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CND3_goodN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("P_miss_VS_Edep_CND3_badN_Step1_epFDn", "Missing Momentum vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];P_{miss} [GeV/c]", 50, 0, 100, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CND3_goodN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_theta_miss_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("theta_miss_VS_Edep_CND3_badN_Step1_epFDn", "#theta_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CND3_goodN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_phi_miss_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("phi_miss_VS_Edep_CND3_badN_Step1_epFDn", "#phi_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#phi_{miss} [#circ]", 50, 0, 100, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("dpp_VS_Edep_CND3_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_dpp_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("dpp_VS_Edep_CND3_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 100, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CND3_goodN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_beta_n_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("beta_n_VS_Edep_CND3_badN_Step1_epFDn", "#beta_{n} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#beta_{n}", 50, 0, 100, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("E_p_VS_Edep_CND3_goodN_Step1_epFDn", "E_{p} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_E_p_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("E_p_VS_Edep_CND3_badN_Step1_epFDn", "E_{p} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{p} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CND3_goodN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_E_miss_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("E_miss_VS_Edep_CND3_badN_Step1_epFDn", "E_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];E_{miss} [GeV]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND3_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND3_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("path_VS_Edep_CND3_goodN_Step1_epFDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_path_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("path_VS_Edep_CND3_badN_Step1_epFDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND3_goodN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND3_badN_Step1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("ToF_VS_Edep_CND3_goodN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_ToF_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("ToF_VS_Edep_CND3_badN_Step1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_badN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("nSector_VS_Edep_CND3_goodN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_nSector_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("nSector_VS_Edep_CND3_badN_Step1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_badN_Step1_epFDn);

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step1_epCDn = new TH2D("Size_CND1_VS_Size_CND2_goodN_Step1_epCDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_goodN_Step1_epCDn);
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step1_epCDn = new TH2D("Size_CND1_VS_Size_CND2_badN_Step1_epCDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_badN_Step1_epCDn);
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step1_epCDn = new TH2D("Size_CND1_VS_Size_CND3_goodN_Step1_epCDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_goodN_Step1_epCDn);
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step1_epCDn = new TH2D("Size_CND1_VS_Size_CND3_badN_Step1_epCDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_badN_Step1_epCDn);
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step1_epCDn = new TH2D("Size_CND2_VS_Size_CND3_goodN_Step1_epCDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_goodN_Step1_epCDn);
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step1_epCDn = new TH2D("Size_CND2_VS_Size_CND3_badN_Step1_epCDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_badN_Step1_epCDn);

    TH2D *h_Size_CND1_VS_Size_CND2_goodN_Step1_epFDn = new TH2D("Size_CND1_VS_Size_CND2_goodN_Step1_epFDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_goodN_Step1_epFDn);
    TH2D *h_Size_CND1_VS_Size_CND2_badN_Step1_epFDn = new TH2D("Size_CND1_VS_Size_CND2_badN_Step1_epFDn", "Size(CND1) vs Size(CND2);Size(CND1);Size(CND2)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND2_badN_Step1_epFDn);
    TH2D *h_Size_CND1_VS_Size_CND3_goodN_Step1_epFDn = new TH2D("Size_CND1_VS_Size_CND3_goodN_Step1_epFDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_goodN_Step1_epFDn);
    TH2D *h_Size_CND1_VS_Size_CND3_badN_Step1_epFDn = new TH2D("Size_CND1_VS_Size_CND3_badN_Step1_epFDn", "Size(CND1) vs Size(CND3);Size(CND1);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND1_VS_Size_CND3_badN_Step1_epFDn);
    TH2D *h_Size_CND2_VS_Size_CND3_goodN_Step1_epFDn = new TH2D("Size_CND2_VS_Size_CND3_goodN_Step1_epFDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_goodN_Step1_epFDn);
    TH2D *h_Size_CND2_VS_Size_CND3_badN_Step1_epFDn = new TH2D("Size_CND2_VS_Size_CND3_badN_Step1_epFDn", "Size(CND2) vs Size(CND3);Size(CND2);Size(CND3)", 50, 0, 10, 50, 0, 10);
    HistoList.push_back(h_Size_CND2_VS_Size_CND3_badN_Step1_epFDn);

    TH1D *h_ToF_goodN_Step1_epCDn = new TH1D("ToF_goodN_Step1_epCDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_goodN_Step1_epCDn);
    TH1D *h_ToF_badN_Step1_epCDn = new TH1D("ToF_badN_Step1_epCDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_badN_Step1_epCDn);
    TH2D *h_P_n_VS_ToF_goodN_Step1_epCDn = new TH2D("P_n_VS_ToF_goodN_Step1_epCDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_ToF_badN_Step1_epCDn = new TH2D("P_n_VS_ToF_badN_Step1_epCDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_badN_Step1_epCDn);
    TH2D *h_theta_n_VS_ToF_goodN_Step1_epCDn = new TH2D("theta_n_VS_ToF_goodN_Step1_epCDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_theta_n_VS_ToF_badN_Step1_epCDn = new TH2D("theta_n_VS_ToF_badN_Step1_epCDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_badN_Step1_epCDn);
    TH2D *h_phi_n_VS_ToF_goodN_Step1_epCDn = new TH2D("phi_n_VS_ToF_goodN_Step1_epCDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_phi_n_VS_ToF_badN_Step1_epCDn = new TH2D("phi_n_VS_ToF_badN_Step1_epCDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_badN_Step1_epCDn);
    TH2D *h_P_miss_VS_ToF_goodN_Step1_epCDn = new TH2D("P_miss_VS_ToF_goodN_Step1_epCDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_P_miss_VS_ToF_badN_Step1_epCDn = new TH2D("P_miss_VS_ToF_badN_Step1_epCDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_badN_Step1_epCDn);
    TH2D *h_theta_miss_VS_ToF_goodN_Step1_epCDn = new TH2D("theta_miss_VS_ToF_goodN_Step1_epCDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_theta_miss_VS_ToF_badN_Step1_epCDn = new TH2D("theta_miss_VS_ToF_badN_Step1_epCDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_badN_Step1_epCDn);
    TH2D *h_phi_miss_VS_ToF_goodN_Step1_epCDn = new TH2D("phi_miss_VS_ToF_goodN_Step1_epCDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_phi_miss_VS_ToF_badN_Step1_epCDn = new TH2D("phi_miss_VS_ToF_badN_Step1_epCDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_badN_Step1_epCDn);
    TH2D *h_dpp_VS_ToF_goodN_Step1_epCDn = new TH2D("dpp_VS_ToF_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_dpp_VS_ToF_badN_Step1_epCDn = new TH2D("dpp_VS_ToF_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_badN_Step1_epCDn);
    TH2D *h_beta_n_VS_ToF_goodN_Step1_epCDn = new TH2D("beta_n_VS_ToF_goodN_Step1_epCDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_beta_n_VS_ToF_badN_Step1_epCDn = new TH2D("beta_n_VS_ToF_badN_Step1_epCDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_badN_Step1_epCDn);
    TH2D *h_E_p_VS_ToF_goodN_Step1_epCDn = new TH2D("E_p_VS_ToF_goodN_Step1_epCDn", "E_{p} vs ToF;ToF [ns];E_{p} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_E_p_VS_ToF_badN_Step1_epCDn = new TH2D("E_p_VS_ToF_badN_Step1_epCDn", "E_{p} vs ToF;ToF [ns];E_{p} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_ToF_badN_Step1_epCDn);
    TH2D *h_E_miss_VS_ToF_goodN_Step1_epCDn = new TH2D("E_miss_VS_ToF_goodN_Step1_epCDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_E_miss_VS_ToF_badN_Step1_epCDn = new TH2D("E_miss_VS_ToF_badN_Step1_epCDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_ToF_goodN_Step1_epCDn = new TH2D("M_miss_VS_ToF_goodN_Step1_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_ToF_badN_Step1_epCDn = new TH2D("M_miss_VS_ToF_badN_Step1_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_badN_Step1_epCDn);
    TH2D *h_path_VS_ToF_goodN_Step1_epCDn = new TH2D("path_VS_ToF_goodN_Step1_epCDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_path_VS_ToF_badN_Step1_epCDn = new TH2D("path_VS_ToF_badN_Step1_epCDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_badN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step1_epCDn = new TH2D("theta_n_miss_VS_ToF_goodN_Step1_epCDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_theta_n_miss_VS_ToF_badN_Step1_epCDn = new TH2D("theta_n_miss_VS_ToF_badN_Step1_epCDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_badN_Step1_epCDn);
    TH2D *h_nSector_VS_ToF_goodN_Step1_epCDn = new TH2D("nSector_VS_ToF_goodN_Step1_epCDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_nSector_VS_ToF_badN_Step1_epCDn = new TH2D("nSector_VS_ToF_badN_Step1_epCDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_badN_Step1_epCDn);

    TH1D *h_ToF_goodN_Step1_epFDn = new TH1D("ToF_goodN_Step1_epFDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_goodN_Step1_epFDn);
    TH1D *h_ToF_badN_Step1_epFDn = new TH1D("ToF_badN_Step1_epFDn", "Neutron ToF;ToF [ns];Counts", 50, 0, 20);
    HistoList.push_back(h_ToF_badN_Step1_epFDn);
    TH2D *h_P_n_VS_ToF_goodN_Step1_epFDn = new TH2D("P_n_VS_ToF_goodN_Step1_epFDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_ToF_badN_Step1_epFDn = new TH2D("P_n_VS_ToF_badN_Step1_epFDn", "Neutron Momentum vs ToF;ToF [ns];P_{n} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_n_VS_ToF_badN_Step1_epFDn);
    TH2D *h_theta_n_VS_ToF_goodN_Step1_epFDn = new TH2D("theta_n_VS_ToF_goodN_Step1_epFDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_theta_n_VS_ToF_badN_Step1_epFDn = new TH2D("theta_n_VS_ToF_badN_Step1_epFDn", "#theta_{n} vs ToF;ToF [ns];#theta_{n} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_VS_ToF_badN_Step1_epFDn);
    TH2D *h_phi_n_VS_ToF_goodN_Step1_epFDn = new TH2D("phi_n_VS_ToF_goodN_Step1_epFDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_phi_n_VS_ToF_badN_Step1_epFDn = new TH2D("phi_n_VS_ToF_badN_Step1_epFDn", "#phi_{n} vs ToF;ToF [ns];#phi_{n} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_n_VS_ToF_badN_Step1_epFDn);
    TH2D *h_P_miss_VS_ToF_goodN_Step1_epFDn = new TH2D("P_miss_VS_ToF_goodN_Step1_epFDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_P_miss_VS_ToF_badN_Step1_epFDn = new TH2D("P_miss_VS_ToF_badN_Step1_epFDn", "Missing Momentum vs ToF;ToF [ns];P_{miss} [GeV/c]", 50, 0, 20, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_ToF_badN_Step1_epFDn);
    TH2D *h_theta_miss_VS_ToF_goodN_Step1_epFDn = new TH2D("theta_miss_VS_ToF_goodN_Step1_epFDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_theta_miss_VS_ToF_badN_Step1_epFDn = new TH2D("theta_miss_VS_ToF_badN_Step1_epFDn", "#theta_{miss} vs ToF;ToF [ns];#theta_{miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_miss_VS_ToF_badN_Step1_epFDn);
    TH2D *h_phi_miss_VS_ToF_goodN_Step1_epFDn = new TH2D("phi_miss_VS_ToF_goodN_Step1_epFDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_phi_miss_VS_ToF_badN_Step1_epFDn = new TH2D("phi_miss_VS_ToF_badN_Step1_epFDn", "#phi_{miss} vs ToF;ToF [ns];#phi_{miss} [#circ]", 50, 0, 20, 50, -180., 180.);
    HistoList.push_back(h_phi_miss_VS_ToF_badN_Step1_epFDn);
    TH2D *h_dpp_VS_ToF_goodN_Step1_epFDn = new TH2D("dpp_VS_ToF_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_dpp_VS_ToF_badN_Step1_epFDn = new TH2D("dpp_VS_ToF_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs ToF;ToF [ns];|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 20, 50, -1.5, 1.5);
    HistoList.push_back(h_dpp_VS_ToF_badN_Step1_epFDn);
    TH2D *h_beta_n_VS_ToF_goodN_Step1_epFDn = new TH2D("beta_n_VS_ToF_goodN_Step1_epFDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_beta_n_VS_ToF_badN_Step1_epFDn = new TH2D("beta_n_VS_ToF_badN_Step1_epFDn", "#beta_{n} vs ToF;ToF [ns];#beta_{n}", 50, 0, 20, 50, -0.1, 1.1);
    HistoList.push_back(h_beta_n_VS_ToF_badN_Step1_epFDn);
    TH2D *h_E_p_VS_ToF_goodN_Step1_epFDn = new TH2D("E_p_VS_ToF_goodN_Step1_epFDn", "E_{p} vs ToF;ToF [ns];E_{p} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_E_p_VS_ToF_badN_Step1_epFDn = new TH2D("E_p_VS_ToF_badN_Step1_epFDn", "E_{p} vs ToF;ToF [ns];E_{p} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_p_VS_ToF_badN_Step1_epFDn);
    TH2D *h_E_miss_VS_ToF_goodN_Step1_epFDn = new TH2D("E_miss_VS_ToF_goodN_Step1_epFDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_E_miss_VS_ToF_badN_Step1_epFDn = new TH2D("E_miss_VS_ToF_badN_Step1_epFDn", "E_{miss} vs ToF;ToF [ns];E_{miss} [GeV]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_VS_ToF_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_ToF_goodN_Step1_epFDn = new TH2D("M_miss_VS_ToF_goodN_Step1_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_ToF_badN_Step1_epFDn = new TH2D("M_miss_VS_ToF_badN_Step1_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.7);
    HistoList.push_back(h_M_miss_VS_ToF_badN_Step1_epFDn);
    TH2D *h_path_VS_ToF_goodN_Step1_epFDn = new TH2D("path_VS_ToF_goodN_Step1_epFDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_path_VS_ToF_badN_Step1_epFDn = new TH2D("path_VS_ToF_badN_Step1_epFDn", "Path length vs ToF;ToF [ns];Path length [cm]", 50, 0, 20, 50, 0., 100.);
    HistoList.push_back(h_path_VS_ToF_badN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_ToF_goodN_Step1_epFDn = new TH2D("theta_n_miss_VS_ToF_goodN_Step1_epFDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_theta_n_miss_VS_ToF_badN_Step1_epFDn = new TH2D("theta_n_miss_VS_ToF_badN_Step1_epFDn", "#theta_{n,miss} vs ToF;ToF [ns];#theta_{n,miss} [#circ]", 50, 0, 20, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_ToF_badN_Step1_epFDn);
    TH2D *h_nSector_VS_ToF_goodN_Step1_epFDn = new TH2D("nSector_VS_ToF_goodN_Step1_epFDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_nSector_VS_ToF_badN_Step1_epFDn = new TH2D("nSector_VS_ToF_badN_Step1_epFDn", "Neutron Sector Number vs ToF;ToF [ns];Sector Number", 50, 0, 20, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_ToF_badN_Step1_epFDn);

    TH1D *h_beta_n_goodN_Step1_epCDn = new TH1D("beta_n_goodN_Step1_epCDn", "#beta_{n} of CND Neutrons;#beta_{n};Counts", 50, 0, 1.1);
    HistoList.push_back(h_beta_n_goodN_Step1_epCDn);
    TH1D *h_beta_n_badN_Step1_epCDn = new TH1D("beta_n_badN_Step1_epCDn", "#beta_{n} of CND Neutrons;#beta_{n};Counts", 50, 0, 1.1);
    HistoList.push_back(h_beta_n_badN_Step1_epCDn);

    TH1D *h_beta_n_goodN_Step1_epFDn = new TH1D("beta_n_goodN_Step1_epFDn", "#beta_{n} of CND Neutrons;#beta_{n};Counts", 50, 0, 1.1);
    HistoList.push_back(h_beta_n_goodN_Step1_epFDn);
    TH1D *h_beta_n_badN_Step1_epFDn = new TH1D("beta_n_badN_Step1_epFDn", "#beta_{n} of CND Neutrons;#beta_{n};Counts", 50, 0, 1.1);
    HistoList.push_back(h_beta_n_badN_Step1_epFDn);

    // TODO: add these to code below
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step1_epCDn = new TH1D("neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step1_epCDn", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 50, 0, 5);
    HistoList.push_back(h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step1_epCDn);
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step1_epCDn = new TH1D("neut_Edep_CND_over_pos_Edep_CTOF_badN_Step1_epCDn", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 50, 0, 5);
    HistoList.push_back(h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step1_epCDn);

    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step1_epFDn = new TH1D("neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step1_epFDn", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 50, 0, 5);
    HistoList.push_back(h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step1_epFDn);
    TH1D *h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step1_epFDn = new TH1D("neut_Edep_CND_over_pos_Edep_CTOF_badN_Step1_epFDn", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 50, 0, 5);
    HistoList.push_back(h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step1_epFDn);

    TH1D *h_Edep_CND_goodN_withNearbyPos_Step1_epCDn = new TH1D("Edep_CND_goodN_withNearbyPos_Step1_epCDn", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 50);
    HistoList.push_back(h_Edep_CND_goodN_withNearbyPos_Step1_epCDn);
    TH1D *h_Edep_CND_badN_withNearbyPos_Step1_epCDn = new TH1D("Edep_CND_badN_withNearbyPos_Step1_epCDn", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 50);
    HistoList.push_back(h_Edep_CND_badN_withNearbyPos_Step1_epCDn);

    TH1D *h_Edep_CND_goodN_withNearbyPos_Step1_epFDn = new TH1D("Edep_CND_goodN_withNearbyPos_Step1_epFDn", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 50);
    HistoList.push_back(h_Edep_CND_goodN_withNearbyPos_Step1_epFDn);
    TH1D *h_Edep_CND_badN_withNearbyPos_Step1_epFDn = new TH1D("Edep_CND_badN_withNearbyPos_Step1_epFDn", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 50);
    HistoList.push_back(h_Edep_CND_badN_withNearbyPos_Step1_epFDn);

    TH1D *h_sdiff_pos_goodN_Step1_layer_epCDn[7];
    TH1D *h_sdiff_pos_badN_Step1_layer_epCDn[7];

    TH1D *h_sdiff_pos_goodN_Step1_layer_epFDn[7];
    TH1D *h_sdiff_pos_badN_Step1_layer_epFDn[7];

    TH2D *h_sdiff_pos_mom_goodN_Step1_layer_epCDn[7];
    TH2D *h_sdiff_pos_mom_badN_Step1_layer_epCDn[7];

    TH2D *h_sdiff_pos_mom_goodN_Step1_layer_epFDn[7];
    TH2D *h_sdiff_pos_mom_badN_Step1_layer_epFDn[7];

    TH2D *h_sdiff_pos_z_goodN_Step1_layer_epCDn[7];
    TH2D *h_sdiff_pos_z_badN_Step1_layer_epCDn[7];

    TH2D *h_sdiff_pos_z_goodN_Step1_layer_epFDn[7];
    TH2D *h_sdiff_pos_z_badN_Step1_layer_epFDn[7];

    TH2D *h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epCDn[7];
    TH2D *h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epCDn[7];

    TH2D *h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epFDn[7];
    TH2D *h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epFDn[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name, "sdiff_pos_goodN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step1_layer_epCDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_badN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step1_layer_epCDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_goodN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step1_layer_epFDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_badN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step1_layer_epFDn[k] = new TH1D(temp_name, temp_title, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step1_layer_epFDn[k]);

        sprintf(temp_name, "sdiff_pos_mom_goodN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step1_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_mom_badN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step1_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_mom_goodN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step1_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_mom_badN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step1_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step1_layer_epFDn[k]);

        sprintf(temp_name, "sdiff_pos_z_goodN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step1_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_z_badN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step1_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_z_goodN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step1_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_z_badN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step1_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_badN_Step1_layer_epFDn[k]);

        sprintf(temp_name, "sdiff_pos_diff_ToFc_z_goodN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name, "sdiff_pos_diff_ToFc_z_badN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epCDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name, "sdiff_pos_diff_ToFc_z_goodN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name, "sdiff_pos_diff_ToFc_z_badN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epFDn[k] = new TH2D(temp_name, temp_title, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epFDn[k]);
    }

    TH2D *h_diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epCDn = new TH2D("diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epCDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epCDn);
    TH2D *h_diff_ToFc_z_VS_Edep_noNear_badN_Step1_epCDn = new TH2D("diff_ToFc_z_VS_Edep_noNear_badN_Step1_epCDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_noNear_badN_Step1_epCDn);

    TH2D *h_diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epFDn = new TH2D("diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epFDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epFDn);
    TH2D *h_diff_ToFc_z_VS_Edep_noNear_badN_Step1_epFDn = new TH2D("diff_ToFc_z_VS_Edep_noNear_badN_Step1_epFDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_noNear_badN_Step1_epFDn);

    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epCDn = new TH2D("diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epCDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epCDn);
    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epCDn = new TH2D("diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epCDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epCDn);

    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epFDn = new TH2D("diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epFDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epFDn);
    TH2D *h_diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epFDn = new TH2D("diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epFDn", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    HistoList.push_back(h_diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epFDn);

    TH2D *h_diff_ToFc_z_Edep_goodN_Step1_layer_epCDn[3];
    TH2D *h_diff_ToFc_z_Edep_badN_Step1_layer_epCDn[3];

    TH2D *h_diff_ToFc_z_Edep_goodN_Step1_layer_epFDn[3];
    TH2D *h_diff_ToFc_z_Edep_badN_Step1_layer_epFDn[3];

    for (int k = 0; k < 3; k++)
    {
        sprintf(temp_name, "diff_ToFc_z_goodN_Step1_layer_%d_epCDn", k + 1);
        sprintf(temp_title, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step1_layer_epCDn[k] = new TH2D(temp_name, temp_title, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name, "diff_ToFc_z_badN_Step1_layer_%d_epCDn", k + 1);
        sprintf(temp_title, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step1_layer_epCDn[k] = new TH2D(temp_name, temp_title, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name, "diff_ToFc_z_goodN_Step1_layer_%d_epFDn", k + 1);
        sprintf(temp_title, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step1_layer_epFDn[k] = new TH2D(temp_name, temp_title, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name, "diff_ToFc_z_badN_Step1_layer_%d_epFDn", k + 1);
        sprintf(temp_title, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step1_layer_epFDn[k] = new TH2D(temp_name, temp_title, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_badN_Step1_layer_epFDn[k]);
    }

    TH2D *h_sdiff_ldiff_allhit_goodN_Step1_epCDn = new TH2D("sdiff_ldiff_allhit_goodN_Step1_epCDn", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    HistoList.push_back(h_sdiff_ldiff_allhit_goodN_Step1_epCDn);
    TH2D *h_sdiff_ldiff_allhit_badN_Step1_epCDn = new TH2D("sdiff_ldiff_allhit_badN_Step1_epCDn", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    HistoList.push_back(h_sdiff_ldiff_allhit_badN_Step1_epCDn);

    TH2D *h_sdiff_ldiff_allhit_goodN_Step1_epFDn = new TH2D("sdiff_ldiff_allhit_goodN_Step1_epFDn", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    HistoList.push_back(h_sdiff_ldiff_allhit_goodN_Step1_epFDn);
    TH2D *h_sdiff_ldiff_allhit_badN_Step1_epFDn = new TH2D("sdiff_ldiff_allhit_badN_Step1_epFDn", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    HistoList.push_back(h_sdiff_ldiff_allhit_badN_Step1_epFDn);

    TH1D *h_numberNearby_goodN_Step1_epCDn = new TH1D("numberNearby_goodN_Step1_epCDn", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    HistoList.push_back(h_numberNearby_goodN_Step1_epCDn);
    TH1D *h_numberNearby_badN_Step1_epCDn = new TH1D("numberNearby_badN_Step1_epCDn", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    HistoList.push_back(h_numberNearby_badN_Step1_epCDn);

    TH1D *h_numberNearby_goodN_Step1_epFDn = new TH1D("numberNearby_goodN_Step1_epFDn", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    HistoList.push_back(h_numberNearby_goodN_Step1_epFDn);
    TH1D *h_numberNearby_badN_Step1_epFDn = new TH1D("numberNearby_badN_Step1_epFDn", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    HistoList.push_back(h_numberNearby_badN_Step1_epFDn);

    TH2D *h_numberNearby_momN_goodN_Step1_epCDn = new TH2D("numberNearby_momN_goodN_Step1_epCDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.5);
    HistoList.push_back(h_numberNearby_momN_goodN_Step1_epCDn);
    TH2D *h_numberNearby_momN_badN_Step1_epCDn = new TH2D("numberNearby_momN_badN_Step1_epCDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.5);
    HistoList.push_back(h_numberNearby_momN_badN_Step1_epCDn);

    TH2D *h_numberNearby_momN_goodN_Step1_epFDn = new TH2D("numberNearby_momN_goodN_Step1_epFDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.5);
    HistoList.push_back(h_numberNearby_momN_goodN_Step1_epFDn);
    TH2D *h_numberNearby_momN_badN_Step1_epFDn = new TH2D("numberNearby_momN_badN_Step1_epFDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.5);
    HistoList.push_back(h_numberNearby_momN_badN_Step1_epFDn);

    TH1D *h_NearbyEdep_goodN_Step1_epCDn = new TH1D("NearbyEdep_goodN_Step1_epCDn", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_NearbyEdep_goodN_Step1_epCDn);
    TH1D *h_NearbyEdep_badN_Step1_epCDn = new TH1D("NearbyEdep_badN_Step1_epCDn", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_NearbyEdep_badN_Step1_epCDn);

    TH1D *h_NearbyEdep_goodN_Step1_epFDn = new TH1D("NearbyEdep_goodN_Step1_epFDn", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_NearbyEdep_goodN_Step1_epFDn);
    TH1D *h_NearbyEdep_badN_Step1_epFDn = new TH1D("NearbyEdep_badN_Step1_epFDn", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    HistoList.push_back(h_NearbyEdep_badN_Step1_epFDn);

    TH1D *h_nsector_goodN_Step1_epCDn = new TH1D("nsector_goodN_Step1_epCDn", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    HistoList.push_back(h_nsector_goodN_Step1_epCDn);
    TH1D *h_nsector_badN_Step1_epCDn = new TH1D("nsector_badN_Step1_epCDn", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    HistoList.push_back(h_nsector_badN_Step1_epCDn);

    TH1D *h_nsector_goodN_Step1_epFDn = new TH1D("nsector_goodN_Step1_epFDn", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    HistoList.push_back(h_nsector_goodN_Step1_epFDn);
    TH1D *h_nsector_badN_Step1_epFDn = new TH1D("nsector_badN_Step1_epFDn", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    HistoList.push_back(h_nsector_badN_Step1_epFDn);

    TH1D *h_phidiff_en_goodN_Step1_epCDn = new TH1D("phidiff_en_goodN_Step1_epCDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 50, 0, 180);
    HistoList.push_back(h_phidiff_en_goodN_Step1_epCDn);
    TH1D *h_phidiff_en_badN_Step1_epCDn = new TH1D("phidiff_en_badN_Step1_epCDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 50, 0, 180);
    HistoList.push_back(h_phidiff_en_badN_Step1_epCDn);

    TH1D *h_phidiff_en_goodN_Step1_epFDn = new TH1D("phidiff_en_goodN_Step1_epFDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 50, 0, 180);
    HistoList.push_back(h_phidiff_en_goodN_Step1_epFDn);
    TH1D *h_phidiff_en_badN_Step1_epFDn = new TH1D("phidiff_en_badN_Step1_epFDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 50, 0, 180);
    HistoList.push_back(h_phidiff_en_badN_Step1_epFDn);

    TH1D *h_TP_goodN_Step1_epCDn = new TH1D("TP_goodN_Step1_epCDn", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    HistoList.push_back(h_TP_goodN_Step1_epCDn);
    TH1D *h_TP_badN_Step1_epCDn = new TH1D("TP_badN_Step1_epCDn", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    HistoList.push_back(h_TP_badN_Step1_epCDn);

    TH1D *h_TP_goodN_Step1_epFDn = new TH1D("TP_goodN_Step1_epFDn", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    HistoList.push_back(h_TP_goodN_Step1_epFDn);
    TH1D *h_TP_badN_Step1_epFDn = new TH1D("TP_badN_Step1_epFDn", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    HistoList.push_back(h_TP_badN_Step1_epFDn);

    TH1D *h_Z_goodN_Step1_epCDn = new TH1D("Z_goodN_Step1_epCDn", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    HistoList.push_back(h_Z_goodN_Step1_epCDn);
    TH1D *h_Z_badN_Step1_epCDn = new TH1D("Z_badN_Step1_epCDn", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    HistoList.push_back(h_Z_badN_Step1_epCDn);

    TH1D *h_Z_goodN_Step1_epFDn = new TH1D("Z_goodN_Step1_epFDn", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    HistoList.push_back(h_Z_goodN_Step1_epFDn);
    TH1D *h_Z_badN_Step1_epFDn = new TH1D("Z_badN_Step1_epFDn", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    HistoList.push_back(h_Z_badN_Step1_epFDn);

    // Step Two (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step2 = new TH1D("n_multiplicity_allN_epCDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epCDn_Step2);
    TH1D *h_n_multiplicity_goodN_epCDn_Step2 = new TH1D("n_multiplicity_goodN_epCDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epCDn_Step2);
    TH1D *h_n_multiplicity_badN_epCDn_Step2 = new TH1D("n_multiplicity_badN_epCDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epCDn_Step2);

    TH1D *h_n_multiplicity_allN_epFDn_Step2 = new TH1D("n_multiplicity_allN_epFDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epFDn_Step2);
    TH1D *h_n_multiplicity_goodN_epFDn_Step2 = new TH1D("n_multiplicity_goodN_epFDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epFDn_Step2);
    TH1D *h_n_multiplicity_badN_epFDn_Step2 = new TH1D("n_multiplicity_badN_epFDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epFDn_Step2);

    // Step Three (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step3 = new TH1D("n_multiplicity_allN_epCDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epCDn_Step3);
    TH1D *h_n_multiplicity_goodN_epCDn_Step3 = new TH1D("n_multiplicity_goodN_epCDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epCDn_Step3);
    TH1D *h_n_multiplicity_badN_epCDn_Step3 = new TH1D("n_multiplicity_badN_epCDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epCDn_Step3);

    TH1D *h_n_multiplicity_allN_epFDn_Step3 = new TH1D("n_multiplicity_allN_epFDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epFDn_Step3);
    TH1D *h_n_multiplicity_goodN_epFDn_Step3 = new TH1D("n_multiplicity_goodN_epFDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epFDn_Step3);
    TH1D *h_n_multiplicity_badN_epFDn_Step3 = new TH1D("n_multiplicity_badN_epFDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epFDn_Step3);

    // Step Four (After applying Phi Diff CND hit cut) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step4 = new TH1D("n_multiplicity_allN_epCDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epCDn_Step4);
    TH1D *h_n_multiplicity_goodN_epCDn_Step4 = new TH1D("n_multiplicity_goodN_epCDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epCDn_Step4);
    TH1D *h_n_multiplicity_badN_epCDn_Step4 = new TH1D("n_multiplicity_badN_epCDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epCDn_Step4);

    TH1D *h_n_multiplicity_allN_epFDn_Step4 = new TH1D("n_multiplicity_allN_epFDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epFDn_Step4);
    TH1D *h_n_multiplicity_goodN_epFDn_Step4 = new TH1D("n_multiplicity_goodN_epFDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epFDn_Step4);
    TH1D *h_n_multiplicity_badN_epFDn_Step4 = new TH1D("n_multiplicity_badN_epFDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epFDn_Step4);

    // Step Five (After event selection cuts) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step5 = new TH1D("n_multiplicity_allN_epCDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epCDn_Step5);
    TH1D *h_n_multiplicity_goodN_epCDn_Step5 = new TH1D("n_multiplicity_goodN_epCDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epCDn_Step5);
    TH1D *h_n_multiplicity_badN_epCDn_Step5 = new TH1D("n_multiplicity_badN_epCDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epCDn_Step5);

    TH1D *h_n_multiplicity_allN_epFDn_Step5 = new TH1D("n_multiplicity_allN_epFDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_allN_epFDn_Step5);
    TH1D *h_n_multiplicity_goodN_epFDn_Step5 = new TH1D("n_multiplicity_goodN_epFDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_goodN_epFDn_Step5);
    TH1D *h_n_multiplicity_badN_epFDn_Step5 = new TH1D("n_multiplicity_badN_epFDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    HistoList.push_back(h_n_multiplicity_badN_epFDn_Step5);

    for (int i = 0; i < HistoList.size(); i++)
    {
        if (HistoList[i]->InheritsFrom("TH1D"))
        {
            HistoList[i]->Sumw2();
        }

        HistoList[i]->GetXaxis()->CenterTitle();
        HistoList[i]->GetYaxis()->CenterTitle();
    }

#pragma endregion /* Veto histograms - end */

    // ======================================================================================================================================================================
    // Chain loop
    // ======================================================================================================================================================================

#pragma region /* Chain loop - start */

    int counter_A = 0; /* From Andrew */

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

    while (chain.Next())
    {
        // Display completed (from Andrew)
        counter_A++;

        if ((counter_A % 1000000) == 0)
        {
            cerr << "\n\n";
            cerr << "\033[33m" << counter_A / 1000000 << " million completed\033[0m\n\n";
        }

#pragma region /* PID & variable definitions - start */

        // PID
        // ===================================================================================================================================================================

        // PID (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        clasAna->Run(c12);

        auto Electrons = clasAna->getByPid(11);
        auto Protons = clasAna->getByPid(2212);
        auto Neutrons = clasAna->getByPid(2112);

        auto AllParticles = c12->getDetParticles();

        // Event selection
        // ===================================================================================================================================================================

        // Event selection (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if (Electrons.size() != 1) // One electron in event
        {
            continue;
        }

        if (Protons.size() != 1) // One proton in event
        {
            continue;
        }

        if (Neutrons.size() < 1) // At least one neutron in event
        {
            continue;
        }

        // Reject particles with the wrong PID
        bool trash = 0;

        for (int i = 0; i < AllParticles.size(); i++)
        {
            int pid = AllParticles[i]->par()->getPid();

            if (pid != 2112 && pid != 11 && pid != 2212 && pid != 0 && pid != 22)
            {
                trash = 1;
            }
        }

        if (trash == 1)
        {
            continue;
        }

        ++counter_epXn;
        numevent = numevent + 1;

        // Variable definitions
        // ===================================================================================================================================================================

        // Variable definitions (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        event = c12->runconfig()->getEvent() << '\n';

        double starttime = c12->event()->getStartTime();

        TVector3 P_b_3v(0, 0, Ebeam);

        // Variable definitions (from Andrew)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        double weight = 1;

        if (isMC)
        {
            weight = c12->mcevent()->getWeight();
        }

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
        double QSq = P_q_3v.Mag2() - (nu * nu); // 4-momentum transfer squared
        double xB = QSq / (2 * mN * nu);        // x Bjorken

        // Electrons (from Andrew)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        double EoP_e = (Electrons[0]->cal(PCAL)->getEnergy() + Electrons[0]->cal(ECIN)->getEnergy() + Electrons[0]->cal(ECOUT)->getEnergy()) / P_e_3v.Mag();
        int nphe = Electrons[0]->che(HTCC)->getNphe();

        int e_sector = Electrons[0]->getSector();

        double theta_q = P_q_3v.Theta() * 180 / M_PI;
        double WSq = (mN * mN) - QSq + (2 * nu * mN); // Hadronic mass
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
        for (int i = 0; i < Protons.size(); i++)
        {
            // define quantities
            P_p_3v.SetMagThetaPhi(Protons[i]->getP(), Protons[i]->getTheta(), Protons[i]->getPhi());
            double dbeta = Protons[i]->par()->getBeta() - P_p_3v.Mag() / sqrt(P_p_3v.Mag2() + mP * mP);
            double p_theta = P_p_3v.Theta() * 180. / M_PI;
            double Vz_p = Protons[i]->par()->getVz();
            double chipid = Protons[i]->par()->getChi2Pid();

            if (Protons[i]->getRegion() == FD)
            {
                ++counter_pFD_multiplicity_BPID;

                h_theta_p_VS_phi_p_BPID_epFD->Fill(P_p_3v.Phi() * 180. / M_PI, p_theta, weight);
                h_P_p_BPID_epFD->Fill(P_p_3v.Mag(), weight);
                h_dbeta_p_VS_P_p_BPID_epFD->Fill(P_p_3v.Mag(), dbeta, weight);
                h_dVz_p_BPID_epFD->Fill(Vz_p - Vz_e, weight);
                h_Chi2pid_p_BPID_epFD->Fill(chipid, weight);

                if (fabs(Vz_p - Vz_e) > 5)
                {
                    continue;
                }

                if (P_p_3v.Mag() < 0.5 || P_p_3v.Mag() > 3.0)
                {
                    continue;
                }

                if (fabs(dbeta) > 0.03)
                {
                    continue;
                }

                ++counter_pFD_multiplicity_APID;

                h_theta_p_VS_phi_p_APID_epFD->Fill(P_p_3v.Phi() * 180. / M_PI, p_theta);
                h_P_p_APID_epFD->Fill(P_p_3v.Mag(), weight);
                h_dbeta_p_VS_P_p_APID_epFD->Fill(P_p_3v.Mag(), dbeta, weight);
                h_dVz_p_APID_epFD->Fill(Vz_p - Vz_e, weight);
                h_Chi2pid_p_APID_epFD->Fill(chipid, weight);
            }
            else if (Protons[i]->getRegion() == CD)
            {
                ++counter_pCD_multiplicity_BPID;

                h_theta_p_VS_phi_p_BPID_epCD->Fill(P_p_3v.Phi() * 180. / M_PI, p_theta);
                h_P_p_BPID_epCD->Fill(P_p_3v.Mag(), weight);
                h_dbeta_p_VS_P_p_BPID_epCD->Fill(P_p_3v.Mag(), dbeta, weight);
                h_dVz_p_BPID_epCD->Fill(Vz_p - Vz_e, weight);
                h_Chi2pid_p_BPID_epCD->Fill(chipid, weight);

                if (fabs(Vz_p - Vz_e) > 4)
                {
                    continue;
                }

                if (P_p_3v.Mag() < 0.3 || P_p_3v.Mag() > 1.5)
                {
                    continue;
                }

                if (fabs(dbeta) > 0.05)
                {
                    continue;
                }

                ++counter_pCD_multiplicity_APID;

                h_theta_p_VS_phi_p_APID_epCD->Fill(P_p_3v.Phi() * 180. / M_PI, p_theta);
                h_P_p_APID_epCD->Fill(P_p_3v.Mag(), weight);
                h_dbeta_p_VS_P_p_APID_epCD->Fill(P_p_3v.Mag(), dbeta, weight);
                h_dVz_p_APID_epCD->Fill(Vz_p - Vz_e, weight);
                h_Chi2pid_p_APID_epCD->Fill(chipid, weight);
            }

            p_index = i;
        }

        h_p_multiplicity_BPID_epCD->Fill(counter_pCD_multiplicity_BPID, weight);
        h_p_multiplicity_BPID_epFD->Fill(counter_pFD_multiplicity_BPID, weight);

        if (p_index < 0)
        {
            continue;
        }

        h_p_multiplicity_APID_epCD->Fill(counter_pCD_multiplicity_APID, weight);
        h_p_multiplicity_APID_epFD->Fill(counter_pFD_multiplicity_APID, weight);

        P_p_3v.SetMagThetaPhi(Protons[p_index]->getP(), Protons[p_index]->getTheta(), Protons[p_index]->getPhi());

        // Determin where is the proton. Moved from angle cuts to getRegion() by the advice of Andrew.
        bool pInFD = (Protons[p_index]->getRegion() == FD); // My addition
        bool pInCD = (Protons[p_index]->getRegion() == CD); // My addition

#pragma endregion /* Protons - end */

#pragma region /* Missing momentum - start */

        // Missing momentum (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Missing momentum, energy, mass
        TVector3 P_miss_3v = P_q_3v - P_p_3v; // TODO: checkout difference from Andrew - he uses leading SRC proton here!

        momentum = P_miss_3v.Mag();

        double E_p = sqrt(mN * mN + P_p_3v.Mag2());
        double E_miss = Ebeam + mD - P_e_3v.Mag() - E_p;
        double M_miss = sqrt((E_miss * E_miss) - P_miss_3v.Mag2());

#pragma endregion /* Missing momentum - end */

        // ==================================================================================================================================================================
        // Andrew's manual work
        // ==================================================================================================================================================================

#pragma region /* Andrew's manual work */

#pragma region /* Missing momentum cuts (Andrew) - start */

        if (pInCD)
        {
            h_P_miss_BmissC_epCD->Fill(P_miss_3v.Mag(), weight);
            h_theta_miss_BmissC_epCD->Fill(P_miss_3v.Theta() * 180 / M_PI, weight);
            h_P_miss_VS_theta_miss_BmissC_epCD->Fill(P_miss_3v.Theta() * 180 / M_PI, P_miss_3v.Mag(), weight);
            h_E_p_BmissC_epCD->Fill(E_p, weight);
            h_E_miss_BmissC_epCD->Fill(E_miss, weight);
            h_M_miss_BmissC_epCD->Fill(M_miss, weight);
            h_xB_BmissC_epCD->Fill(xB, weight);
            h_xB_VS_M_miss_BmissC_epCD->Fill(xB, M_miss, weight);
        }
        else if (pInFD)
        {
            h_P_miss_BmissC_epFD->Fill(P_miss_3v.Mag(), weight);
            h_theta_miss_BmissC_epFD->Fill(P_miss_3v.Theta() * 180 / M_PI, weight);
            h_P_miss_VS_theta_miss_BmissC_epFD->Fill(P_miss_3v.Theta() * 180 / M_PI, P_miss_3v.Mag(), weight);
            h_E_p_BmissC_epFD->Fill(E_p, weight);
            h_E_miss_BmissC_epFD->Fill(E_miss, weight);
            h_M_miss_BmissC_epFD->Fill(M_miss, weight);
            h_xB_BmissC_epFD->Fill(xB, weight);
            h_xB_VS_M_miss_BmissC_epFD->Fill(xB, M_miss, weight);
        }

        if (P_miss_3v.Theta() * 180 / M_PI < 40 || P_miss_3v.Theta() * 180 / M_PI > 135)
        {
            continue;
        }

        if (P_miss_3v.Mag() < 0.2 || P_miss_3v.Mag() > 1.5)
        {
            continue;
        }

        if (M_miss < 0.7 || M_miss > 1.2)
        {
            continue;
        }

        if (pInCD)
        {
            h_P_miss_AmissC_epCD->Fill(P_miss_3v.Mag(), weight);
            h_theta_miss_AmissC_epCD->Fill(P_miss_3v.Theta() * 180 / M_PI, weight);
            h_P_miss_VS_theta_miss_AmissC_epCD->Fill(P_miss_3v.Theta() * 180 / M_PI, P_miss_3v.Mag(), weight);
            h_E_p_AmissC_epCD->Fill(E_p, weight);
            h_E_miss_AmissC_epCD->Fill(E_miss, weight);
            h_M_miss_AmissC_epCD->Fill(M_miss, weight);
            h_xB_AmissC_epCD->Fill(xB, weight);
            h_xB_VS_M_miss_AmissC_epCD->Fill(xB, M_miss, weight);
        }
        else if (pInFD)
        {
            h_P_miss_AmissC_epFD->Fill(P_miss_3v.Mag(), weight);
            h_theta_miss_AmissC_epFD->Fill(P_miss_3v.Theta() * 180 / M_PI, weight);
            h_P_miss_VS_theta_miss_AmissC_epFD->Fill(P_miss_3v.Theta() * 180 / M_PI, P_miss_3v.Mag(), weight);
            h_E_p_AmissC_epFD->Fill(E_p, weight);
            h_E_miss_AmissC_epFD->Fill(E_miss, weight);
            h_M_miss_AmissC_epFD->Fill(M_miss, weight);
            h_xB_AmissC_epFD->Fill(xB, weight);
            h_xB_VS_M_miss_AmissC_epFD->Fill(xB, M_miss, weight);
        }

#pragma endregion /* Missing momentum cuts (Andrew) - end */

#pragma region /* Neutrons (Andrew) */

        bool pass_step0_cuts = false, pass_step1_cuts = false, pass_step2_cuts = false, pass_step3_cuts = false, pass_step4_cuts = false, pass_step5_cuts = false;

        /////////////////////////////////////
        // Lead Neutron Checks
        /////////////////////////////////////
        for (int itr1 = 0; itr1 < AllParticles.size(); itr1++)
        {
            if (AllParticles[itr1]->par()->getCharge() != 0) // Cut out charged particles
            {
                continue;
            }

            // TODO: Confirm that these actually working! Try to move the Erin's variabels?
            bool CT = (AllParticles[itr1]->sci(clas12::CTOF)->getDetector() == 4);
            bool C1 = (AllParticles[itr1]->sci(clas12::CND1)->getDetector() == 3);
            bool C2 = (AllParticles[itr1]->sci(clas12::CND2)->getDetector() == 3);
            bool C3 = (AllParticles[itr1]->sci(clas12::CND3)->getDetector() == 3);

            if (!(C1 || C2 || C3)) // Cut out neutrons without a CND hit in one of it's layers
            {
                continue;
            }

            // Why this cut? reco code bug. Neutrons in this angle range are in the BAND and appear in the CND.
            // This bug is probobly fixed, yet the cut is still applied to mak sure.
            if (AllParticles[itr1]->getTheta() * 180 / M_PI > 160)
            {
                continue;
            }

            // Explicit calculation of the neutron's momentum (to bypass cases where P_n is E_dep)
            double theta = AllParticles[itr1]->getTheta() * 180 / M_PI;
            double beta = AllParticles[itr1]->par()->getBeta();
            double gamma = 1 / sqrt(1 - (beta * beta));
            double mom = gamma * beta * mN;
            double ToF = AllParticles[itr1]->getTime() - starttime;

            int detINTlayer = C1 ? 1 : C2 ? 2
                                          : 3;
            auto detlayer = C1 ? CND1 : C2 ? CND2
                                           : CND3; // CND layer with hit
            double Edep_CND1 = AllParticles[itr1]->sci(CND1)->getEnergy();
            double Edep_CND2 = AllParticles[itr1]->sci(CND2)->getEnergy();
            double Edep_CND3 = AllParticles[itr1]->sci(CND3)->getEnergy();
            double Edep_CND = Edep_CND1 + Edep_CND2 + Edep_CND3;
            double Edep_CTOF = AllParticles[itr1]->sci(CTOF)->getEnergy();
            double Edep_single = AllParticles[itr1]->sci(detlayer)->getEnergy();

            double Size_CND1 = AllParticles[itr1]->sci(CND1)->getSize();
            double Size_CND2 = AllParticles[itr1]->sci(CND2)->getSize();
            double Size_CND3 = AllParticles[itr1]->sci(CND3)->getSize();

            double nvtx_x = AllParticles[itr1]->par()->getVx();
            double nvtx_y = AllParticles[itr1]->par()->getVy();
            double nvtx_z = AllParticles[itr1]->par()->getVz();
            TVector3 v_nvtx_3v(nvtx_x, nvtx_y, nvtx_z); // Neutron's vertex location

            TVector3 v_hit_3v; // Neutron's hit location in CND
            v_hit_3v.SetXYZ(AllParticles[itr1]->sci(detlayer)->getX(), AllParticles[itr1]->sci(detlayer)->getY(), AllParticles[itr1]->sci(detlayer)->getZ());

            TVector3 v_path_3v = v_hit_3v - v_nvtx_3v; // Direct calculation of neutron's path (in vector form)
            TVector3 P_n_3v;
            P_n_3v.SetMagThetaPhi(mom, v_path_3v.Theta(), v_path_3v.Phi()); // Direct calculation of neutron momentum?
                                                                            // TODO: check with Andrew why he calculated this explicitly

            // Why "v_path_3v.Mag() / 100"? unit conversion.
            // TODO: check if this unit conversion is needed!
            double path = v_path_3v.Mag() / 100;
            // double path = v_path_3v.Mag();
            double theta_n_miss = P_n_3v.Angle(P_miss_3v) * 180 / M_PI; // Opening angle between calculated neutron's momentum and predicted neutron momentum (= missing momentum)
            double dpp = (P_miss_3v.Mag() - P_n_3v.Mag()) / P_miss_3v.Mag();
            int nSector = AllParticles[itr1]->sci(detlayer)->getSector(); // Number of CND sector with a neutron hit in the layer detlayer

            // Check to see if there is a good neutron
            bool isGN = false;

            // Why this cut? reco code bug. Neutrons in this angle range are in the BAND and appear in the CND.
            // This bug is probobly fixed, yet the cut is still applied to mak sure.
            if (P_n_3v.Theta() * 180. / M_PI > 160)
            {
                continue;
            }

            if ((theta_n_miss < 25.) && (dpp > -0.3) && (dpp < 0.3)) // Good neutron definition
            {
                isGN = true;
            }

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn, counter_n_multiplicity_goodN_epCDn, counter_n_multiplicity_badN_epCDn,
                               counter_n_multiplicity_allN_epFDn, counter_n_multiplicity_goodN_epFDn, counter_n_multiplicity_badN_epFDn);
            // SetNeutronCounters(isGN, counter_n_multiplicity_allN, counter_n_multiplicity_goodN, counter_n_multiplicity_badN);

            // FILL HISTOS FOR NEUTRON CANDIDATES
            if (pInCD)
            {
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

                if (isGN)
                {
                    h_dpp_goodN_epCDn->Fill(dpp, weight);
                    h_theta_n_miss_goodN_epCDn->Fill(theta_n_miss, weight);
                }
                else
                {
                    h_dpp_badN_epCDn->Fill(dpp, weight);
                    h_theta_n_miss_badN_epCDn->Fill(theta_n_miss, weight);
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
                h_theta_n_p_VS_P_p_epCDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

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

                h_Edep_single_epCDn->Fill(Edep_single, weight);
                h_P_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n_3v.Mag(), weight);
                h_theta_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n_3v.Theta() * 180. / M_PI, weight);
                h_phi_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n_3v.Phi() * 180. / M_PI, weight);
                h_P_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss_3v.Mag(), weight);
                h_theta_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss_3v.Theta() * 180. / M_PI, weight);
                h_phi_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss_3v.Phi() * 180. / M_PI, weight);
                h_dpp_VS_Edep_single_epCDn->Fill(Edep_single, dpp, weight);
                h_beta_n_VS_Edep_single_epCDn->Fill(Edep_single, beta, weight);
                h_E_miss_VS_Edep_single_epCDn->Fill(Edep_single, E_miss, weight);
                h_M_miss_VS_Edep_single_epCDn->Fill(Edep_single, M_miss, weight);
                h_path_VS_Edep_single_epCDn->Fill(Edep_single, path * 100, weight);
                h_theta_n_miss_VS_Edep_single_epCDn->Fill(Edep_single, theta_n_miss, weight);
                h_ToF_VS_Edep_single_epCDn->Fill(Edep_single, ToF, weight);
                h_nSector_VS_Edep_single_epCDn->Fill(Edep_single, nSector, weight);

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

                h_Size_CND1_VS_Size_CND2_epCDn->Fill(Size_CND1, Size_CND2, weight);
                h_Size_CND1_VS_Size_CND3_epCDn->Fill(Size_CND1, Size_CND3, weight);
                h_Size_CND2_VS_Size_CND3_epCDn->Fill(Size_CND2, Size_CND3, weight);

                h_ToF_epCDn->Fill(ToF, weight);
                h_P_n_VS_ToF_epCDn->Fill(ToF, P_n_3v.Mag(), weight);
                h_theta_n_VS_ToF_epCDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
                h_phi_n_VS_ToF_epCDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
                h_P_miss_VS_ToF_epCDn->Fill(ToF, P_miss_3v.Mag(), weight);
                h_theta_miss_VS_ToF_epCDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
                h_phi_miss_VS_ToF_epCDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
                h_dpp_VS_ToF_epCDn->Fill(ToF, dpp, weight);
                h_beta_n_VS_ToF_epCDn->Fill(ToF, beta, weight);
                h_E_miss_VS_ToF_epCDn->Fill(ToF, E_miss, weight);
                h_M_miss_VS_ToF_epCDn->Fill(ToF, M_miss, weight);
                h_path_VS_ToF_epCDn->Fill(ToF, path * 100, weight);
                h_theta_n_miss_VS_ToF_epCDn->Fill(ToF, theta_n_miss, weight);
                h_nSector_VS_ToF_epCDn->Fill(ToF, nSector, weight);
            }
            else if (pInFD)
            {
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

                if (isGN)
                {
                    h_dpp_goodN_epFDn->Fill(dpp, weight);
                    h_theta_n_miss_goodN_epFDn->Fill(theta_n_miss, weight);
                }
                else
                {
                    h_dpp_badN_epFDn->Fill(dpp, weight);
                    h_theta_n_miss_badN_epFDn->Fill(theta_n_miss, weight);
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
                h_theta_n_p_VS_P_p_epFDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

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

                h_Edep_single_epFDn->Fill(Edep_single, weight);
                h_P_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n_3v.Mag(), weight);
                h_theta_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n_3v.Theta() * 180. / M_PI, weight);
                h_phi_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n_3v.Phi() * 180. / M_PI, weight);
                h_P_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss_3v.Mag(), weight);
                h_theta_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss_3v.Theta() * 180. / M_PI, weight);
                h_phi_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss_3v.Phi() * 180. / M_PI, weight);
                h_dpp_VS_Edep_single_epFDn->Fill(Edep_single, dpp, weight);
                h_beta_n_VS_Edep_single_epFDn->Fill(Edep_single, beta, weight);
                h_E_miss_VS_Edep_single_epFDn->Fill(Edep_single, E_miss, weight);
                h_M_miss_VS_Edep_single_epFDn->Fill(Edep_single, M_miss, weight);
                h_path_VS_Edep_single_epFDn->Fill(Edep_single, path * 100, weight);
                h_theta_n_miss_VS_Edep_single_epFDn->Fill(Edep_single, theta_n_miss, weight);
                h_ToF_VS_Edep_single_epFDn->Fill(Edep_single, ToF, weight);
                h_nSector_VS_Edep_single_epFDn->Fill(Edep_single, nSector, weight);

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

                h_Size_CND1_VS_Size_CND2_epFDn->Fill(Size_CND1, Size_CND2, weight);
                h_Size_CND1_VS_Size_CND3_epFDn->Fill(Size_CND1, Size_CND3, weight);
                h_Size_CND2_VS_Size_CND3_epFDn->Fill(Size_CND2, Size_CND3, weight);

                h_ToF_epFDn->Fill(ToF, weight);
                h_P_n_VS_ToF_epFDn->Fill(ToF, P_n_3v.Mag(), weight);
                h_theta_n_VS_ToF_epFDn->Fill(ToF, P_n_3v.Theta() * 180. / M_PI, weight);
                h_phi_n_VS_ToF_epFDn->Fill(ToF, P_n_3v.Phi() * 180. / M_PI, weight);
                h_P_miss_VS_ToF_epFDn->Fill(ToF, P_miss_3v.Mag(), weight);
                h_theta_miss_VS_ToF_epFDn->Fill(ToF, P_miss_3v.Theta() * 180. / M_PI, weight);
                h_phi_miss_VS_ToF_epFDn->Fill(ToF, P_miss_3v.Phi() * 180. / M_PI, weight);
                h_dpp_VS_ToF_epFDn->Fill(ToF, dpp, weight);
                h_beta_n_VS_ToF_epFDn->Fill(ToF, beta, weight);
                h_E_miss_VS_ToF_epFDn->Fill(ToF, E_miss, weight);
                h_M_miss_VS_ToF_epFDn->Fill(ToF, M_miss, weight);
                h_path_VS_ToF_epFDn->Fill(ToF, path * 100, weight);
                h_theta_n_miss_VS_ToF_epFDn->Fill(ToF, theta_n_miss, weight);
                h_nSector_VS_ToF_epFDn->Fill(ToF, nSector, weight);
            }

            //////////////////////////////////////////////
            // Step Zero
            //////////////////////////////////////////////

#pragma region /* Step Zero - start */

            if (pInCD)
            {
                h_dbeta_n_VS_P_n_BS0C_Step0_epCDn->Fill(P_n_3v.Mag(), beta - (path * 100) / (ToF * c), weight);
                h_dbeta_n_VS_ToF_BS0C_Step0_epCDn->Fill(ToF, beta - (path * 100) / (ToF * c), weight);

                h_Vhit_z_n_BS0C_Step0_epCDn->Fill(v_hit_3v.Z(), weight);

                h_ToF_n_BS0C_Step0_epCDn->Fill(ToF, weight);
            }
            else if (pInFD)
            {
                h_dbeta_n_VS_P_n_BS0C_Step0_epFDn->Fill(P_n_3v.Mag(), beta - (path * 100) / (ToF * c), weight);
                h_dbeta_n_VS_ToF_BS0C_Step0_epFDn->Fill(ToF, beta - (path * 100) / (ToF * c), weight);

                h_Vhit_z_n_BS0C_Step0_epFDn->Fill(v_hit_3v.Z(), weight);

                h_ToF_n_BS0C_Step0_epFDn->Fill(ToF, weight);
            }

            // Why "path * 100"? unit conversion. Path is in cm; tof is in ns.
            // TODO: check if this unit conversion is needed!
            if (fabs(beta - (path * 100) / (ToF * c)) > 0.01) // A cut on delta beta
            // if (fabs(beta - path / (ToF * c)) > 0.01) // A cut on delta beta
            {
                continue;
            }

            // A cut on the z-component of the CND hit
            // This is a fiducial cut on the range that the CND can reach on the z-axis
            if (v_hit_3v.Z() > 45 || v_hit_3v.Z() < -40)
            {
                continue;
            }

            if (ToF < 0 || ToF > 20)
            {
                continue;
            }

            pass_step0_cuts = true;

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step0, counter_n_multiplicity_goodN_epCDn_Step0, counter_n_multiplicity_badN_epCDn_Step0,
                               counter_n_multiplicity_allN_epFDn_Step0, counter_n_multiplicity_goodN_epFDn_Step0, counter_n_multiplicity_badN_epFDn_Step0);
            // SetNeutronCounters(isGN, counter_n_multiplicity_allN_Step0, counter_n_multiplicity_goodN_Step0, counter_n_multiplicity_badN_Step0);

            if (pInCD)
            {
                h_dbeta_n_VS_P_n_AS0C_Step0_epCDn->Fill(P_n_3v.Mag(), beta - (path * 100) / (ToF * c), weight);
                h_dbeta_n_VS_ToF_AS0C_Step0_epCDn->Fill(ToF, beta - (path * 100) / (ToF * c), weight);

                h_Vhit_z_n_AS0C_Step0_epCDn->Fill(v_hit_3v.Z(), weight);

                h_ToF_n_AS0C_Step0_epCDn->Fill(ToF, weight);
            }
            else if (pInFD)
            {
                h_dbeta_n_VS_P_n_AS0C_Step0_epFDn->Fill(P_n_3v.Mag(), beta - (path * 100) / (ToF * c), weight);
                h_dbeta_n_VS_ToF_AS0C_Step0_epFDn->Fill(ToF, beta - (path * 100) / (ToF * c), weight);

                h_Vhit_z_n_AS0C_Step0_epFDn->Fill(v_hit_3v.Z(), weight);

                h_ToF_n_AS0C_Step0_epFDn->Fill(ToF, weight);
            }

            if (pInCD)
            {
                h_dpp_allN_Step0_epCDn->Fill(dpp, weight);
                h_theta_n_miss_allN_Step0_epCDn->Fill(theta_n_miss, weight);
                h_dpp_VS_theta_n_miss_allN_Step0_epCDn->Fill(dpp, theta_n_miss, weight);

                if (isGN)
                {
                    h_theta_n_goodN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                    h_phi_n_goodN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_goodN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_goodN_Step0_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                    h_P_n_goodN_Step0_epCDn->Fill(P_n_3v.Mag(), weight);
                    h_P_n_VS_theta_n_goodN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                    h_P_miss_goodN_Step0_epCDn->Fill(P_miss_3v.Mag(), weight);
                    h_P_miss_VS_theta_miss_goodN_Step0_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
                    h_P_miss_VS_phi_miss_goodN_Step0_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

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
                    h_theta_n_p_VS_P_p_goodN_Step0_epCDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

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
                    h_theta_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                }
                else
                {
                    h_theta_n_badN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                    h_phi_n_badN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_badN_Step0_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_badN_Step0_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                    h_P_n_badN_Step0_epCDn->Fill(P_n_3v.Mag(), weight);
                    h_P_n_VS_theta_n_badN_Step0_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                    h_P_miss_badN_Step0_epCDn->Fill(P_miss_3v.Mag(), weight);
                    h_P_miss_VS_theta_miss_badN_Step0_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
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
                    h_theta_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
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
            }
            else if (pInFD)
            {
                h_dpp_allN_Step0_epFDn->Fill(dpp, weight);
                h_theta_n_miss_allN_Step0_epFDn->Fill(theta_n_miss, weight);
                h_dpp_VS_theta_n_miss_allN_Step0_epFDn->Fill(dpp, theta_n_miss, weight);

                if (isGN)
                {
                    h_theta_n_goodN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                    h_phi_n_goodN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_goodN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_goodN_Step0_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                    h_P_n_goodN_Step0_epFDn->Fill(P_n_3v.Mag(), weight);
                    h_P_n_VS_theta_n_goodN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                    h_P_miss_goodN_Step0_epFDn->Fill(P_miss_3v.Mag(), weight);
                    h_P_miss_VS_theta_miss_goodN_Step0_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
                    h_P_miss_VS_phi_miss_goodN_Step0_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

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
                    h_theta_n_p_VS_P_p_goodN_Step0_epFDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

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
                    h_theta_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                }
                else
                {
                    h_theta_n_badN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                    h_phi_n_badN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_badN_Step0_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_badN_Step0_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                    h_P_n_badN_Step0_epFDn->Fill(P_n_3v.Mag(), weight);
                    h_P_n_VS_theta_n_badN_Step0_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                    h_P_miss_badN_Step0_epFDn->Fill(P_miss_3v.Mag(), weight);
                    h_P_miss_VS_theta_miss_badN_Step0_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
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
                    h_theta_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
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

#pragma endregion /* Step Zero - end */

            //////////////////////////////////////////////
            // Step One
            //////////////////////////////////////////////

#pragma region /* Step One - start */

            // Step One = Beta cut & Dep. energy cut

            // Beta cut:
            // Upper: beta > 0.8 -> cut out photons
            // Lower: beta < 0.15 ->
            if (beta < 0.15 || beta > 0.8)
            {
                continue;
            }

            // Total deposited energy in CND cut:
            // Upper: Edep_CND > (gamma - 1) * mN * 1000 -> the neutron's deposited energy should not exceed its relativistic kinematic energy. Factor 1000 -> convert GeV to MeV!
            // Lower: Edep_CND < 5 ->
            // TODO: add lower Edep_CND cut?
            if (Edep_CND < 5 || Edep_CND > (gamma - 1) * mN * 1000)
            // if (Edep_CND < 5)
            {
                continue;
            }

            pass_step1_cuts = true;

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step1, counter_n_multiplicity_goodN_epCDn_Step1, counter_n_multiplicity_badN_epCDn_Step1,
                               counter_n_multiplicity_allN_epFDn_Step1, counter_n_multiplicity_goodN_epFDn_Step1, counter_n_multiplicity_badN_epFDn_Step1);
            // SetNeutronCounters(isGN, counter_n_multiplicity_allN_Step1, counter_n_multiplicity_goodN_Step1, counter_n_multiplicity_badN_Step1);

            if (pInCD)
            {
                h_dpp_allN_Step1_epCDn->Fill(dpp, weight);
                h_theta_n_miss_allN_Step1_epCDn->Fill(theta_n_miss, weight);
                h_dpp_VS_theta_n_miss_allN_Step1_epCDn->Fill(dpp, theta_n_miss, weight);

                if (isGN)
                {
                    h_theta_n_goodN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                    h_phi_n_goodN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_goodN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_goodN_Step1_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                    h_P_n_goodN_Step1_epCDn->Fill(P_n_3v.Mag(), weight);
                    h_P_n_VS_theta_n_goodN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                    h_P_miss_goodN_Step1_epCDn->Fill(P_miss_3v.Mag(), weight);
                    h_P_miss_VS_theta_miss_goodN_Step1_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
                    h_P_miss_VS_phi_miss_goodN_Step1_epCDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

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
                    h_theta_n_p_VS_P_p_goodN_Step1_epCDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

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
                    h_theta_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                }
                else
                {
                    h_theta_n_badN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                    h_phi_n_badN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_badN_Step1_epCDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_badN_Step1_epCDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                    h_P_n_badN_Step1_epCDn->Fill(P_n_3v.Mag(), weight);
                    h_P_n_VS_theta_n_badN_Step1_epCDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                    h_P_miss_badN_Step1_epCDn->Fill(P_miss_3v.Mag(), weight);
                    h_P_miss_VS_theta_miss_badN_Step1_epCDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
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
                    h_theta_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
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
            }
            else if (pInFD)
            {
                h_dpp_allN_Step1_epFDn->Fill(dpp, weight);
                h_theta_n_miss_allN_Step1_epFDn->Fill(theta_n_miss, weight);
                h_dpp_VS_theta_n_miss_allN_Step1_epFDn->Fill(dpp, theta_n_miss, weight);

                if (isGN)
                {
                    h_theta_n_goodN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                    h_phi_n_goodN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_goodN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_goodN_Step1_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                    h_P_n_goodN_Step1_epFDn->Fill(P_n_3v.Mag(), weight);
                    h_P_n_VS_theta_n_goodN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                    h_P_miss_goodN_Step1_epFDn->Fill(P_miss_3v.Mag(), weight);
                    h_P_miss_VS_theta_miss_goodN_Step1_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
                    h_P_miss_VS_phi_miss_goodN_Step1_epFDn->Fill(P_miss_3v.Phi() * 180. / M_PI, P_miss_3v.Mag(), weight);

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
                    h_theta_n_p_VS_P_p_goodN_Step1_epFDn->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

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
                    h_theta_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                }
                else
                {
                    h_theta_n_badN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                    h_phi_n_badN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_badN_Step1_epFDn->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_badN_Step1_epFDn->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                    h_P_n_badN_Step1_epFDn->Fill(P_n_3v.Mag(), weight);
                    h_P_n_VS_theta_n_badN_Step1_epFDn->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                    h_P_miss_badN_Step1_epFDn->Fill(P_miss_3v.Mag(), weight);
                    h_P_miss_VS_theta_miss_badN_Step1_epFDn->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);
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
                    h_theta_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_miss_3v.Theta() * 180. / M_PI, weight);
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
                    h_theta_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_miss_3v.Theta() * 180. / M_PI, weight);
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

            bool CNDVeto = false;

            if (ToF * c - v_hit_3v.Z() < 70) // TODO: find a way to check what is this cut
            {
                if (pInCD)
                {
                    if (isGN)
                    {
                        h_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, weight);
                    }
                    else
                    {
                        h_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, weight);
                    }
                }
                else if (pInFD)
                {
                    if (isGN)
                    {
                        h_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, weight);
                    }
                    else
                    {
                        h_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, weight);
                    }
                }

                for (int itr2 = 0; itr2 < AllParticles.size(); itr2++)
                {
                    if (itr2 == 0) // Why skip itr2 == 0? it is the electron
                    {
                        continue;
                    }

                    if (itr2 == itr1)
                    {
                        continue;
                    }

                    // Cut negatively charged particles
                    // TODO: Maybe it is good to keep the nagativly charged particles in the future.
                    if (AllParticles[itr2]->par()->getCharge() <= 0)
                    {
                        continue;
                    }

                    // Why this cut? because the background (protons) have high probability of hitting the CTOF? all charged particles supposed to have a CTOF hit at the rime of writing the code
                    if (AllParticles[itr2]->sci(CTOF)->getDetector() == 0) // Cut out particles WITHOUT a CTOF hit
                    {
                        continue;
                    }

                    // TODO: what is this? check for sectors with proton hits in any of the layers of the CND and CTOF?
                    int vetoSectorbyLayer[4] = {(AllParticles[itr2]->sci(CTOF)->getComponent() + 1) / 2, // Normalizes CTOF components to CND sectors (since vetoSectorbyLayer is an array if integers)
                                                AllParticles[itr2]->sci(CND1)->getSector(),
                                                AllParticles[itr2]->sci(CND2)->getSector(),
                                                AllParticles[itr2]->sci(CND3)->getSector()};

                    TVector3 p_C_3v; // Momentum of the charged particle in the itr2-th entry of AllParticles
                    p_C_3v.SetMagThetaPhi(AllParticles[itr2]->getP(), AllParticles[itr2]->getTheta(), AllParticles[itr2]->getPhi());

                    double edep_pos = AllParticles[itr2]->sci(clas12::CTOF)->getEnergy(); // E_dep of positivly charged particle

                    for (int itr3 = 0; itr3 < 4; itr3++) //
                    {
                        if (vetoSectorbyLayer[itr3] == 0) // TODO: why this cut? no hit in the itr3-th layer?
                        {
                            continue;
                        }

                        int sdiff = nSector - vetoSectorbyLayer[itr3];

                        // sdiff normalization
                        if (sdiff <= -12)
                        {
                            sdiff += 24;
                        }
                        else if (sdiff > 12)
                        {
                            sdiff -= 24;
                        }

                        int ldiff = detINTlayer - itr3;

                        if (pInCD)
                        {
                            if (isGN) // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step1_layer
                            {
                                h_sdiff_pos_goodN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, weight);
                                h_sdiff_pos_mom_goodN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
                                h_sdiff_pos_z_goodN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, v_hit_3v.Z(), weight);
                                h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, ToF * c - v_hit_3v.Z(), weight);
                            }
                            else
                            {
                                h_sdiff_pos_badN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, weight);
                                h_sdiff_pos_mom_badN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
                                h_sdiff_pos_z_badN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, v_hit_3v.Z(), weight);
                                h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, ToF * c - v_hit_3v.Z(), weight);
                            }
                        }
                        else if (pInFD)
                        {
                            if (isGN) // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step1_layer
                            {
                                h_sdiff_pos_goodN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, weight);
                                h_sdiff_pos_mom_goodN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
                                h_sdiff_pos_z_goodN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, v_hit_3v.Z(), weight);
                                h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, ToF * c - v_hit_3v.Z(), weight);
                            }
                            else
                            {
                                h_sdiff_pos_badN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, weight);
                                h_sdiff_pos_mom_badN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, p_C_3v.Perp(), weight);
                                h_sdiff_pos_z_badN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, v_hit_3v.Z(), weight);
                                h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, ToF * c - v_hit_3v.Z(), weight);
                            }
                        }

                        if (isPosNear(sdiff, ldiff))
                        {
                            CNDVeto = true;
                        }
                    } // End of loop over vetoSectorbyLayer

                    if (CNDVeto)
                    {
                        if (pInCD)
                        {
                            if (isGN)
                            {
                                h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CND / edep_pos, weight);
                            }
                            else
                            {
                                h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CND / edep_pos, weight);
                            }
                        }
                        else if (pInFD)
                        {
                            if (isGN)
                            {
                                h_neut_Edep_CND_over_pos_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CND / edep_pos, weight);
                            }
                            else
                            {
                                h_neut_Edep_CND_over_pos_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CND / edep_pos, weight);
                            }
                        }
                    }
                } // End of second loop over AllParticles (step 1)

                if (CNDVeto)
                {
                    if (pInCD)
                    {
                        if (isGN)
                        {
                            h_Edep_CND_goodN_withNearbyPos_Step1_epCDn->Fill(Edep_CND, weight);
                        }
                        else
                        {
                            h_Edep_CND_badN_withNearbyPos_Step1_epCDn->Fill(Edep_CND, weight);
                        }
                    }
                    else if (pInFD)
                    {
                        if (isGN)
                        {
                            h_Edep_CND_goodN_withNearbyPos_Step1_epFDn->Fill(Edep_CND, weight);
                        }
                        else
                        {
                            h_Edep_CND_badN_withNearbyPos_Step1_epFDn->Fill(Edep_CND, weight);
                        }
                    }
                }

                if (pInCD)
                {
                    if (isGN)
                    {
                        if (!CNDVeto)
                            h_diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epCDn->Fill(ToF * c - v_hit_3v.Z(), Edep_CND, weight);
                        else
                        {
                            h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epCDn->Fill(ToF * c - v_hit_3v.Z(), Edep_CND, weight);
                        }
                    }
                    else
                    {
                        if (!CNDVeto)
                            h_diff_ToFc_z_VS_Edep_noNear_badN_Step1_epCDn->Fill(ToF * c - v_hit_3v.Z(), Edep_CND, weight);
                        else
                        {
                            h_diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epCDn->Fill(ToF * c - v_hit_3v.Z(), Edep_CND, weight);
                        }
                    }
                }
                else if (pInFD)
                {
                    if (isGN)
                    {
                        if (!CNDVeto)
                            h_diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epFDn->Fill(ToF * c - v_hit_3v.Z(), Edep_CND, weight);
                        else
                        {
                            h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epFDn->Fill(ToF * c - v_hit_3v.Z(), Edep_CND, weight);
                        }
                    }
                    else
                    {
                        if (!CNDVeto)
                            h_diff_ToFc_z_VS_Edep_noNear_badN_Step1_epFDn->Fill(ToF * c - v_hit_3v.Z(), Edep_CND, weight);
                        else
                        {
                            h_diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epFDn->Fill(ToF * c - v_hit_3v.Z(), Edep_CND, weight);
                        }
                    }
                }
            }

#pragma endregion /* Step One - end */

            //////////////////////////////////////////////
            // Step Two
            //////////////////////////////////////////////

#pragma region /* Step Two - start */

            /*
            // Step two = cut/veto out neutrons with charged particles close by

            if (CNDVeto)
            {
                continue;
            }

            pass_step2_cuts = true;

            SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step2, counter_n_multiplicity_goodN_epCDn_Step2, counter_n_multiplicity_badN_epCDn_Step2,
                               counter_n_multiplicity_allN_epFDn_Step2, counter_n_multiplicity_goodN_epFDn_Step2, counter_n_multiplicity_badN_epFDn_Step2);
            SetNeutronCounters(isGN, counter_n_multiplicity_allN_Step2, counter_n_multiplicity_goodN_Step2, counter_n_multiplicity_badN_Step2);

            h_pnRes_theta_n_miss_Step2->Fill(dpp, theta_n_miss, weight);

            if (isGN)
            {
                h_theta_n_goodN_Step2->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                h_phi_n_goodN_Step2->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                h_theta_n_VS_phi_n_goodN_Step2->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                h_theta_n_VS_beta_n_goodN_Step2->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                h_P_n_goodN_Step2->Fill(P_n_3v.Mag(), weight);
                h_P_n_VS_theta_n_goodN_Step2->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                h_P_miss_goodN_Step2->Fill(P_miss_3v.Mag(), weight);
                h_P_miss_VS_theta_miss_goodN_Step2->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);

                h_P_n_minus_P_miss_goodN_Step2->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
                h_P_n_x_minus_P_miss_x_goodN_Step2->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
                h_P_n_y_minus_P_miss_y_goodN_Step2->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
                h_P_n_z_minus_P_miss_z_goodN_Step2->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

                h_P_n_VS_P_miss_goodN_Step2->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
                h_P_n_x_VS_P_miss_x_goodN_Step2->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
                h_P_n_y_VS_P_miss_y_goodN_Step2->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
                h_P_n_z_VS_P_miss_z_goodN_Step2->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

                if (pInCD)
                {
                    h_E_p_CD_goodN_Step2->Fill(E_p, weight);
                }
                else if (pInFD)
                {
                    h_E_p_FD_goodN_Step2->Fill(E_p, weight);
                }

                h_E_miss_goodN_Step2->Fill(E_miss, weight);
                h_M_miss_goodN_Step2->Fill(M_miss, weight);
                h_M_miss_VS_P_n_goodN_Step2->Fill(P_n_3v.Mag(), M_miss, weight);
                h_M_miss_VS_P_miss_goodN_Step2->Fill(P_miss_3v.Mag(), M_miss, weight);

                h_theta_n_p_VS_P_p_goodN_Step2->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

                h_xB_goodN_Step2->Fill(xB, weight);

                h_Edep_goodN_Step2->Fill(edep, weight);
                h_P_n_VS_Edep_goodN_Step2->Fill(edep, P_n_3v.Mag(), weight);
                h_P_miss_VS_Edep_goodN_Step2->Fill(edep, P_miss_3v.Mag(), weight);

                h_dpp_VS_Edep_goodN_Step2->Fill(edep, dpp, weight);

                h_ToF_goodN_Step2->Fill(ToF, weight);
                h_pmiss_goodN_Step2->Fill(P_miss_3v.Mag(), weight);
                h_Edep_ToF_goodN_Step2->Fill(ToF, edep, weight);

                h_ToF_goodN_Step2->Fill(ToF, weight);
                h_Edep_ToF_goodN_Step2->Fill(ToF, edep, weight);
            }
            else
            {
                h_theta_n_badN_Step2->Fill(P_n_3v.Theta() * 180. / M_PI, weight);
                h_phi_n_badN_Step2->Fill(P_n_3v.Phi() * 180. / M_PI, weight);
                h_theta_n_VS_phi_n_badN_Step2->Fill(P_n_3v.Phi() * 180. / M_PI, P_n_3v.Theta() * 180. / M_PI, weight);
                h_theta_n_VS_beta_n_badN_Step2->Fill(beta, P_n_3v.Theta() * 180. / M_PI, weight);

                h_P_n_badN_Step2->Fill(P_n_3v.Mag(), weight);
                h_P_n_VS_theta_n_badN_Step2->Fill(P_n_3v.Theta() * 180. / M_PI, P_n_3v.Mag(), weight);

                h_P_miss_badN_Step2->Fill(P_miss_3v.Mag(), weight);
                h_P_miss_VS_theta_miss_badN_Step2->Fill(P_miss_3v.Theta() * 180. / M_PI, P_miss_3v.Mag(), weight);

                h_P_n_minus_P_miss_badN_Step2->Fill(P_n_3v.Mag() - P_miss_3v.Mag(), weight);
                h_P_n_x_minus_P_miss_x_badN_Step2->Fill(P_n_3v.X() - P_miss_3v.X(), weight);
                h_P_n_y_minus_P_miss_y_badN_Step2->Fill(P_n_3v.Y() - P_miss_3v.Y(), weight);
                h_P_n_z_minus_P_miss_z_badN_Step2->Fill(P_n_3v.Z() - P_miss_3v.Z(), weight);

                h_P_n_VS_P_miss_badN_Step2->Fill(P_n_3v.Mag(), P_miss_3v.Mag(), weight);
                h_P_n_x_VS_P_miss_x_badN_Step2->Fill(P_n_3v.X(), P_miss_3v.X(), weight);
                h_P_n_y_VS_P_miss_y_badN_Step2->Fill(P_n_3v.Y(), P_miss_3v.Y(), weight);
                h_P_n_z_VS_P_miss_z_badN_Step2->Fill(P_n_3v.Z(), P_miss_3v.Z(), weight);

                if (pInCD)
                {
                    h_E_p_CD_badN_Step2->Fill(E_p, weight);
                }
                else if (pInFD)
                {
                    h_E_p_FD_badN_Step2->Fill(E_p, weight);
                }

                h_E_miss_badN_Step2->Fill(E_miss, weight);
                h_M_miss_badN_Step2->Fill(M_miss, weight);
                h_M_miss_VS_P_n_badN_Step2->Fill(P_n_3v.Mag(), M_miss, weight);
                h_M_miss_VS_P_miss_badN_Step2->Fill(P_miss_3v.Mag(), M_miss, weight);

                h_theta_n_p_VS_P_p_badN_Step2->Fill(P_p_3v.Mag(), P_p_3v.Angle(P_n_3v) * 180. / M_PI, weight);

                h_xB_badN_Step2->Fill(xB, weight);

                h_Edep_badN_Step2->Fill(edep, weight);
                h_P_n_VS_Edep_badN_Step2->Fill(edep, P_n_3v.Mag(), weight);
                h_P_miss_VS_Edep_badN_Step2->Fill(edep, P_miss_3v.Mag(), weight);

                h_dpp_VS_Edep_badN_Step2->Fill(edep, dpp, weight);

                h_ToF_badN_Step2->Fill(ToF, weight);
                h_Edep_ToF_badN_Step2->Fill(ToF, edep, weight);
            }

            for (int itr4 = 0; itr4 < AllParticles.size(); itr4++)
            {
                if (itr4 == 0) // Why skip itr4 == 0? it is the electron
                {
                    continue;
                }

                if (itr4 == itr1)
                {
                    continue;
                }

                // Cut negatively charged particles
                // TODO: Maybe it is good to keep the nagativly charged particles in the future.
                if (AllParticles[itr4]->par()->getCharge() <= 0)
                {
                    continue;
                }

                // Why this cut? because the background (protons) have high probability of hitting the CTOF? all charged particles supposed to have a CTOF hit at the rime of writing the code
                if (AllParticles[itr4]->sci(CTOF)->getDetector() == 0) // Cut out particles WITHOUT a CTOF hit
                {
                    continue;
                }

                int vetoSectorbyLayer[4] = {(AllParticles[itr4]->sci(CTOF)->getComponent() + 1) / 2,
                                            AllParticles[itr4]->sci(CND1)->getSector(),
                                            AllParticles[itr4]->sci(CND2)->getSector(),
                                            AllParticles[itr4]->sci(CND3)->getSector()};

                for (int itr5 = 0; itr5 < 4; itr5++)
                {
                    if (vetoSectorbyLayer[itr5] == 0) // TODO: why this cut? no hit in the itr5-th layer?
                    {
                        continue;
                    }

                    int sdiff = nSector - vetoSectorbyLayer[itr5];

                    // sdiff normalization
                    if (sdiff <= -12)
                    {
                        sdiff += 24;
                    }
                    else if (sdiff > 12)
                    {
                        sdiff -= 24;
                    }

                    int ldiff = detINTlayer - itr5;

                    if (isGN) // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step2_layer
                    {
                        h_sdiff_pos_goodN_Step2_layer[ldiff + 3]->Fill(sdiff, weight);
                    }
                    else
                    {
                        h_sdiff_pos_badN_Step2_layer[ldiff + 3]->Fill(sdiff, weight);
                    }
                }
            } // End of third loop over AllParticles (step 2)

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
                        else
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
                    else
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
                if (isGN)
                {
                    h_numberNearby_goodN_Step2->Fill(hitsNear, weight);
                    h_numberNearby_momN_goodN_Step2->Fill(hitsNear, mom, weight);
                    h_nsector_goodN_Step2->Fill(nSector, weight);
                }
                else
                {
                    h_numberNearby_badN_Step2->Fill(hitsNear, weight);
                    h_numberNearby_momN_badN_Step2->Fill(hitsNear, mom, weight);
                    h_nsector_badN_Step2->Fill(nSector, weight);
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
            else
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
                else
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
                else
                {
                    h_sdiff_ldiff_CTOFhit_badN_Step3->Fill(sdiff, ldiff, weight);
                }
            }

            if (isGN)
            {
                h_numberCTOF_goodN_Step3->Fill(hitsCTOF, weight);
                h_numberCTOF_momN_goodN_Step3->Fill(hitsCTOF, mom, weight);
            }
            else
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
            else
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
            else
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
                    else
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
                    else
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
                    else
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
                    else
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

        } // End of Andrew's loop over all particles

        // Fill neutron multiplicity plots
        if (pInCD)
        {
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
        }
        else if (pInFD)
        {
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

        // Count events passing steps
        if (pass_step0_cuts)
        {
            ++counter_pass_step0_cuts;
        }

        if (pass_step1_cuts)
        {
            ++counter_pass_step1_cuts;
        }

        if (pass_step2_cuts)
        {
            ++counter_pass_step2_cuts;
        }

        if (pass_step3_cuts)
        {
            ++counter_pass_step3_cuts;
        }

        if (pass_step4_cuts)
        {
            ++counter_pass_step4_cuts;
        }

        if (pass_step5_cuts)
        {
            ++counter_pass_step5_cuts;
        }

#pragma endregion /* Neutrons */

#pragma endregion /* Andrew's manual work */

    } // closes event loop

#pragma endregion /* Chain loop - end */

    // ======================================================================================================================================================================
    // Wrap up
    // ======================================================================================================================================================================

#pragma region /* Wrap up - start */

    HistPrinter(HistoList, PDFFile);

#pragma endregion /* Wrap up - end */

    // ======================================================================================================================================================================
    // Save log file
    // ======================================================================================================================================================================

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

    myLogFile << "Total #(events) in sample:\t" << counter_A << "\n\n";

    myLogFile << "Total #((e,e'pXn) events):\t" << counter_epXn << "\n\n";

    myLogFile << "Total #((e,e'pXn) events) pass_step0_cuts:\t" << counter_pass_step0_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step1_cuts:\t" << counter_pass_step1_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step2_cuts:\t" << counter_pass_step2_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step3_cuts:\t" << counter_pass_step3_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step4_cuts:\t" << counter_pass_step4_cuts << "\n";
    myLogFile << "Total #((e,e'pXn) events) pass_step5_cuts:\t" << counter_pass_step5_cuts << "\n\n";

    myLogFile.close();

#pragma endregion /* Save log file - end */

    // ======================================================================================================================================================================
    // Printouts
    // ======================================================================================================================================================================

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

    if (Elapsed_time_seconds.count() < 60)
    {
        std::cout << "\033[33mRunning time:\033[0m\t\t" << Elapsed_time_seconds.count() << " seconds\n\n";
    }
    else
    {
        std::cout << "\033[33mRunning time:\033[0m\t\t" << to_string_with_precision(Elapsed_time_minutes, 3) << " minutes\n\n";
    }

#pragma endregion /* Printouts 2 - end */

    return 0;

} // closes main function