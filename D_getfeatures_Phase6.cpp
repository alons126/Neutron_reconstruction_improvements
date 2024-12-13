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

#pragma region /* Erin main function - start */

int D_getfeatures_Phase6(                                                                             //
    const string OutDir, string output_pdf_Erin,                                                      // My arguments
    double Ebeam, bool keep_good, string output_root_Erin, string output_txt_Erin, string input_hipo, // Erin's arguments
    string PDFFile, int isMC = 0                                                                      // Andrew's arguments
)
// int main(int argc, char **argv)
{
    auto Code_start_time = std::chrono::system_clock::now(); // Start counting running time

    const bool Run_Erins_features = false;
    const bool Run_Andrews_work = true;

    // ======================================================================================================================================================================
    // Printouts
    // ======================================================================================================================================================================

#pragma region /* Printouts - start */

    cout << "\033[33m\n\033[0m";
    cout << "\033[33mHIPO_FILES:\033[0m\t\t" << gSystem->Getenv("HIPO_FILES") << "\n";
    cout << "\033[33minput_hipo:\033[0m\t\t" << input_hipo << "\n";
    cout << "\033[33m\n\033[0m";
    cout << "\033[33mOUTDIR:\033[0m\t\t\t" << gSystem->Getenv("OUTDIR") << "\n";
    cout << "\033[33mOutDir:\033[0m\t\t\t" << OutDir << "\n";
    cout << "\033[33moutput_pdf_Erin:\033[0m\t" << output_pdf_Erin << "\n";
    cout << "\033[33moutput_root_Erin:\033[0m\t" << output_root_Erin << "\n";
    cout << "\033[33moutput_txt_Erin:\033[0m\t" << output_txt_Erin << "\n";
    cout << "\033[33mPDFFile:\033[0m\t\t" << PDFFile << "\n\n";

#pragma endregion /* Printouts - end */

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

    // prepare histograms
    vector<TH1 *> hist_list_1;
    vector<TH2 *> hist_list_2;

    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);

    gStyle->SetTitleXOffset(0.8);
    gStyle->SetTitleYOffset(0.8);

    char temp_name[100];
    char temp_title[100];

    // Set up root tree for TMVA
    Int_t nhits;
    double px, py, pz, momentum;
    Int_t sec[100] = {-1};
    Int_t lay[100] = {-1};
    int event;
    double energy, cnd_energy, ctof_energy, angle_diff;
    int layermult, size, cnd_hits, ctof_hits;
    bool is_CTOF, is_CND1, is_CND2, is_CND3;

    // TODO: add Andrew's variables to this TTree
    ntree->Branch("momentum", &momentum, "momentum/D");
    ntree->Branch("energy", &energy, "energy/D");
    ntree->Branch("layermult", &layermult, "layermult/I");
    ntree->Branch("size", &size, "size/I");
    ntree->Branch("cnd_hits", &cnd_hits, "cnd_hits/I");
    ntree->Branch("cnd_energy", &cnd_energy, "cnd_energy/D");
    ntree->Branch("ctof_energy", &ctof_energy, "ctof_energy/D");
    ntree->Branch("ctof_hits", &ctof_hits, "ctof_hits/I");
    ntree->Branch("angle_diff", &angle_diff, "angle_diff/D");

    int counter = 0;
    cout << endl;

    // set up instance of clas12ana
    clas12ana *clasAna = new clas12ana();

    clasAna->readEcalSFPar("src/cuts/paramsSF_LD2_x2.dat"); // TODO: check if applied
    clasAna->readEcalPPar("src/cuts/paramsPI_LD2_x2.dat");  // TODO: check if applied

    // clasAna->printParams();

    clasAna->setProtonPidCuts(true);

#pragma endregion /* Initial setup - end */

    // ======================================================================================================================================================================
    // Erin's histograms
    // ======================================================================================================================================================================

#pragma region /* Erin's histograms - start */

    // Proton histograms (Erin)
    // ======================================================================================================================================================================
    TH1D *h_psize = new TH1D("psize", "Number of Protons in Event", 10, 0, 10);
    hist_list_1.push_back(h_psize);
    TH2D *h_pangles = new TH2D("pangles", "Proton Angular Distribution;#phi_{p} (deg);#theta_{p} (deg)", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_pangles);

    TH2D *h_dbeta_p_cd = new TH2D("dbeta_p_cd", "#Delta #beta vs proton momentum (CD);p_{p} (GeV/c);#Delta#beta", 50, 0, 3, 50, -0.2, 0.2);
    hist_list_2.push_back(h_dbeta_p_cd);
    TH1D *h_vzp_cd = new TH1D("vzp_cd", "Vertex difference between proton and electron (CD);Proton Vertex z - Electron Vertex z (cm);Counts", 100, -8, 8);
    hist_list_1.push_back(h_vzp_cd);
    TH1D *h_chipid_cd = new TH1D("chipid_cd", "Proton #chi^{2}_{PID} (CD);#chi^{2}_{PID};Counts", 100, -6, 6);
    hist_list_1.push_back(h_chipid_cd);

    TH2D *h_dbeta_p_fd = new TH2D("dbeta_p_fd", "#Delta #beta vs proton momentum (FD);p_{p} (GeV/c);#Delta#beta", 50, 0, 3, 50, -0.2, 0.2);
    hist_list_2.push_back(h_dbeta_p_fd);
    TH1D *h_vzp_fd = new TH1D("vzp_fd", "Vertex difference between proton and electron (FD);Proton Vertex z - Electron Vertex z (cm);Counts", 100, -8, 8);
    hist_list_1.push_back(h_vzp_fd);
    TH1D *h_chipid_fd = new TH1D("chipid_fd", "Proton #chi^{2}_{PID} (FD);#chi^{2}_{PID};Counts", 100, -6, 6);
    hist_list_1.push_back(h_chipid_fd);

    // Neutron histograms (Erin)
    // ======================================================================================================================================================================
    TH1D *h_n_multiplicity = new TH1D("n_multiplicity", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity);

    // Reconstructed momentum (Erin)
    // ======================================================================================================================================================================
    TH2D *h_nangles = new TH2D("nangles", "Neutron Angular Distribution;#phi_{n};#theta_{n}", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_nangles);
    TH1D *h_pxminuspx = new TH1D("pxminuspx", "p_{n,x}-p_{x,pred};p_{n,x}-p_{x,pred} (GeV/c);Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pxminuspx);
    TH1D *h_pyminuspy = new TH1D("pyminuspy", "p_{n,y}-p_{y,pred};p_{n,y}-p_{y,pred} (GeV/c);Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pyminuspy);
    TH1D *h_pzminuspz = new TH1D("pzminuspz", "p_{n,z}-p_{z,pred};p_{n,z}-p_{z,pred} (GeV/c);Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pzminuspz);
    TH1D *h_pminusp = new TH1D("pminusp", "P_{n}-p_{pred};P_{n}-p_{pred} (GeV/c);Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pminusp);
    TH2D *h_pvsp = new TH2D("pvsp", "Momentum Resolution;p_{pred} (GeV/c);p_{measured} (GeV/c)", 100, 0, 1, 100, 0, 1);
    hist_list_2.push_back(h_pvsp);
    TH1D *h_cos0 = new TH1D("cos0", "Cosine of angle between generated and reconstructed p", 50, -1.1, 1.1);
    hist_list_1.push_back(h_cos0);
    TH2D *h_dpp = new TH2D("dpp", "Momentum Resolution;p_{generated} (GeV/c);#Delta p/p", 100, 0, 1, 100, -0.4, 0.4);
    hist_list_2.push_back(h_dpp);
    TH1D *h_energy = new TH1D("energy", "Neutron Energy Deposition;Energy (MeV);Counts", 100, 0, 25);
    hist_list_1.push_back(h_energy);
    TH2D *h_sec_phi = new TH2D("sec_phi", "Sector vs Phi of CND hits;#phi (deg);Sector", 90, 0, 360, 25, 0, 25);
    hist_list_2.push_back(h_sec_phi);
    TH1D *h_mmiss = new TH1D("mmiss", "Missing Mass", 50, 0.5, 1.5);
    hist_list_1.push_back(h_mmiss);
    TH2D *h_mmiss_pn = new TH2D("mmiss_pn", "Missing Mass vs Measured Neutron Momentum;P_{n} (GeV/c);Missing Mass (GeV/c^{2})", 50, 0.25, 1., 50, 0.5, 1.5);
    hist_list_2.push_back(h_mmiss_pn);
    TH2D *h_mmiss_pmiss = new TH2D("mmiss_pmiss", "Missing Mass vs Expected Neutron Momentum;p_{pred} (GeV/c);Missing Mass (GeV/c^{2})", 50, 0, 1.25, 50, 0.5, 1.5);
    hist_list_2.push_back(h_mmiss_pmiss);
    TH2D *h_mmiss_xb = new TH2D("mmiss_xb", "Missing Mass vs x_{B};x_{B};Missing Mass (GeV/c^{2})", 50, 0, 1.25, 50, 0.5, 1.5);
    hist_list_2.push_back(h_mmiss_xb);

    TH2D *h_theta_beta = new TH2D("theta_beta", "Neutron theta vs beta;#beta;#theta", 50, -0.1, 1.1, 55, 35, 145);
    hist_list_2.push_back(h_theta_beta);
    TH2D *h_p_theta = new TH2D("p_theta", "Neutron Momentum vs Theta;#theta;p (GeV/c)", 55, 35, 145, 50, 0, 1.2);
    hist_list_2.push_back(h_p_theta);
    TH2D *h_pmiss_thetamiss = new TH2D("pmiss_thetamiss", "Missing Momentum vs #theta_{pred};#theta_{pred} (deg);p_{pred} (GeV/c)", 90, 0, 180, 50, 0, 1.3);
    hist_list_2.push_back(h_pmiss_thetamiss);
    TH2D *h_thetapn_pp = new TH2D("thetapn_pp", "#theta_{pn} vs p_{p};p_{p} (GeV/c);#theta_{pn}", 40, 0, 1, 40, 0, 180);
    hist_list_2.push_back(h_thetapn_pp);
    TH1D *h_tof = new TH1D("tof", "Time of Flight;TOF (ns);Counts", 200, -10, 50);
    hist_list_1.push_back(h_tof);
    TH2D *h_compare = new TH2D("compare", "(p_{pred}-P_{n})/p_{pred} vs #theta_{n,pred} (deg);(p_{pred}-P_{n})/p_{pred};#theta_{n,pred}", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_compare);
    TH2D *h_Edep_beta = new TH2D("Edep_beta", "Energy deposition vs #beta;#beta;E_{dep}", 50, 0, 1, 50, 0, 100);
    hist_list_2.push_back(h_Edep_beta);
    TH1D *h_p_all = new TH1D("p_all", "Momentum", 100, 0, 1.2);
    hist_list_1.push_back(h_p_all);

    TH2D *h_dpp_edep = new TH2D("dpp_edep", "#Delta p/p vs Energy Deposition;E_{dep} (MeVee);#Delta p/p", 50, 0, 40, 50, -0.4, 0.4);
    hist_list_2.push_back(h_dpp_edep);

    // good n / bad n set (Erin)
    // ======================================================================================================================================================================
    TH2D *h_nangles2 = new TH2D("nangles2", "Neutron Angles;#phi;theta", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_nangles2);
    TH1D *h_pxminuspx2 = new TH1D("pxminuspx2", "p_{n,x}-p_{x,pred};p_{n,x}-p_{x,pred};Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pxminuspx2);
    TH1D *h_pyminuspy2 = new TH1D("pyminuspy2", "p_{n,y}-p_{y,pred};p_{n,y}-p_{y,pred};Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pyminuspy2);
    TH1D *h_pzminuspz2 = new TH1D("pzminuspz2", "p_{n,z}-p_{z,pred};p_{n,z}-p_{z,pred};Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pzminuspz2);
    TH1D *h_pminusp2 = new TH1D("pminusp2", "P_{n}-p_{gen};P_{n}-p_{gen};Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pminusp2);
    TH2D *h_pvsp2 = new TH2D("pvsp2", "Momentum Resolution;p_{pred} (GeV/c);p_{measured} (GeV/c)", 100, 0, 1, 100, 0, 1);
    hist_list_2.push_back(h_pvsp2);
    TH1D *h_cos02 = new TH1D("cos02", "Cosine of angle between generated and reconstructed p", 50, -1.1, 1.1);
    hist_list_1.push_back(h_cos02);
    TH2D *h_dpp2 = new TH2D("dpp2", "Momentum Resolution;p_{generated} (GeV/c);#Delta p/p", 100, 0, 1, 100, -0.4, 0.4);
    hist_list_2.push_back(h_dpp2);
    TH1D *h_mmiss2 = new TH1D("mmiss2", "Missing Mass", 50, 0.5, 1.5);
    hist_list_1.push_back(h_mmiss2);
    TH2D *h_mmiss_pn2 = new TH2D("mmiss_pn2", "Missing Mass vs Measured Neutron Momentum", 50, 0, 1.25, 50, 0.5, 1.5);
    hist_list_2.push_back(h_mmiss_pn2);
    TH1D *h_energy2 = new TH1D("energy2", "Neutron Energy Deposition;Energy (MeV);Counts", 100, 0, 25);
    hist_list_1.push_back(h_energy2);
    TH2D *h_theta_beta2 = new TH2D("theta_beta2", "Neutron theta vs beta;#beta;#theta", 50, -0.1, 1.1, 55, 35, 145);
    hist_list_2.push_back(h_theta_beta2);
    TH2D *h_p_theta2 = new TH2D("p_theta2", "Neutron Momentum vs Theta;#theta;p (GeV/c)", 55, 35, 145, 50, 0, 1.2);
    hist_list_2.push_back(h_p_theta2);
    TH2D *h_pmiss_thetamiss2 = new TH2D("pmiss_thetamiss2", "p_{pred} vs #theta_{pred};#theta_{pred};p_{pred} (GeV/c)", 90, 0, 180, 50, 0, 1.3);
    hist_list_2.push_back(h_pmiss_thetamiss2);
    TH2D *h_thetapn_pp2 = new TH2D("thetapn_pp2", "#theta_{pn} vs p_{p};p_{p} (GeV/c);#theta_{pn}", 40, 0, 1, 40, 0, 180);
    hist_list_2.push_back(h_thetapn_pp2);
    TH1D *h_tof2 = new TH1D("tof2", "Time of Flight;TOF (ns);Counts", 200, -10, 50);
    hist_list_1.push_back(h_tof2);
    TH2D *h_compare2 = new TH2D("compare2", "(p_{pred}-P_{n})/p_{pred} vs #theta_{n,pred};(p_{pred}-P_{n})/p_{pred};#theta_{n,pred} (deg)", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_compare2);
    TH2D *h_Edep_beta2 = new TH2D("Edep_beta2", "Energy deposition vs #beta;#beta;E_{dep}", 50, 0, 1, 50, 0, 100);
    hist_list_2.push_back(h_Edep_beta2);
    TH1D *h_p_cut = new TH1D("p_cut", "Momentum", 100, 0, 1.2);
    hist_list_1.push_back(h_p_cut);

    TH2D *h_thetapn_dpp = new TH2D("thetapn_dpp", "#theta_{pn} vs #Delta p/p;(p_{pred}-P_{n})/p_{pred};#theta_{pn} (deg)", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_thetapn_dpp);
    TH2D *h_thetapn_dpp1 = new TH2D("thetapn_dpp1", "#theta_{pn} vs #Delta p/p;(p_{pred}-P_{n})/p_{pred};#theta_{pn} (deg)", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_thetapn_dpp1);
    TH2D *h_thetapn_dpp2 = new TH2D("thetapn_dpp2", "#theta_{pn} vs #Delta p/p;(p_{pred}-P_{n})/p_{pred};#theta_{pn} (deg)", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_thetapn_dpp2);

    TH1D *h_anglediff = new TH1D("angle_diff", "Angle Diff between CVT hit and CND hit", 180, 0, 180);
    hist_list_1.push_back(h_anglediff);
    TH1D *h_anglediff2 = new TH1D("angle_diff2", "Angle Diff between CVT hit and CND hit", 180, 0, 180);
    hist_list_1.push_back(h_anglediff2);

    TH2D *h_ptheta_pred = new TH2D("ptheta_pred", "Predicted Momentum vs Angle for Final Signal Sample;#theta_{pred} (deg);p_{pred} (GeV/c)", 110, 35, 145, 100, 0.2, 1.3);
    hist_list_2.push_back(h_ptheta_pred);
    TH2D *h_ptheta = new TH2D("ptheta", "Measured Momentum vs Angle for Final Signal Sample;#theta_{n} (deg);P_{n} (GeV/c)", 110, 35, 145, 100, 0.2, 1.3);
    hist_list_2.push_back(h_ptheta);

    // ML features - all neutron candidates (Erin)
    // ======================================================================================================================================================================
    TH1D *h_energy_1 = new TH1D("f_energy_1", "Neutron Energy", 100, 0, 100);
    hist_list_1.push_back(h_energy_1);
    TH1D *h_layermult_1 = new TH1D("f_layermult_1", "CND Layer Mult", 4, 0, 4);
    hist_list_1.push_back(h_layermult_1);
    TH1D *h_size_1 = new TH1D("f_size_1", "Cluster Size", 5, 0, 5);
    hist_list_1.push_back(h_size_1);
    TH1D *h_cnd_hits_1 = new TH1D("f_cnd_hits_1", "Nearby CND Hits", 10, 0, 10);
    hist_list_1.push_back(h_cnd_hits_1);
    TH1D *h_cnd_energy_1 = new TH1D("f_cnd_energy_1", "Nearby CND Energy", 100, 0, 100);
    hist_list_1.push_back(h_cnd_energy_1);
    TH1D *h_ctof_energy_1 = new TH1D("f_ctof_energy_1", "Nearby CTOF Energy", 100, 0, 100);
    hist_list_1.push_back(h_ctof_energy_1);
    TH1D *h_ctof_hits_1 = new TH1D("f_ctof_hits_1", "Nearby CTOF Hits", 10, 0, 10);
    hist_list_1.push_back(h_ctof_hits_1);
    TH1D *h_anglediff_1 = new TH1D("f_anglediff_1", "CVT Angle Diff", 200, 0, 200);
    hist_list_1.push_back(h_anglediff_1);

    // ML features - signal/background only (Erin)
    // ======================================================================================================================================================================
    TH1D *h_energy_2 = new TH1D("f_energy_2", "Neutron Energy", 100, 0, 100);
    hist_list_1.push_back(h_energy_2);
    TH1D *h_layermult_2 = new TH1D("f_layermult_2", "CND Layer Mult", 4, 0, 4);
    hist_list_1.push_back(h_layermult_2);
    TH1D *h_size_2 = new TH1D("f_size_2", "Cluster Size", 5, 0, 5);
    hist_list_1.push_back(h_size_2);
    TH1D *h_cnd_hits_2 = new TH1D("f_cnd_hits_2", "Nearby CND Hits", 10, 0, 10);
    hist_list_1.push_back(h_cnd_hits_2);
    TH1D *h_cnd_energy_2 = new TH1D("f_cnd_energy_2", "Nearby CND Energy", 100, 0, 100);
    hist_list_1.push_back(h_cnd_energy_2);
    TH1D *h_ctof_energy_2 = new TH1D("f_ctof_energy_2", "Nearby CTOF Energy", 100, 0, 100);
    hist_list_1.push_back(h_ctof_energy_2);
    TH1D *h_ctof_hits_2 = new TH1D("f_ctof_hits_2", "Nearby CTOF Hits", 10, 0, 10);
    hist_list_1.push_back(h_ctof_hits_2);
    TH1D *h_anglediff_2 = new TH1D("f_anglediff_2", "CVT Angle Diff", 200, 0, 200);
    hist_list_1.push_back(h_anglediff_2);

#pragma endregion /* Erin's histograms - end */

    // ======================================================================================================================================================================
    // Andrew's histograms
    // ======================================================================================================================================================================

#pragma region /* Andrew's histograms - start */

    /////////////////////////////////////
    // Prepare histograms
    /////////////////////////////////////

    vector<TH1 *> HistoList;

    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);

    gStyle->SetTitleXOffset(0.8);
    gStyle->SetTitleYOffset(0.8);

    char temp_name_A[100];
    char temp_title_A[100];

    // (e,e'p) plots
    // ======================================================================================================================================================================

    /* Proton histograms (from Erin) */
    TH1D *h_p_multiplicity_BPID_epCD = new TH1D("p_multiplicity_BPID_epCD", "Number of CD Protons in Event (Before PID);CD Proton multiplicity;Counts", 10, 0, 10);
    HistoList.push_back(h_p_multiplicity_BPID_epCD);
    TH1D *h_P_p_BPID_epCD = new TH1D("P_p_BPID_epCD", "CD Proton momentum (Before PID);P_{p} [GeV/C];Counts", 50, 0., 3.);
    HistoList.push_back(h_P_p_BPID_epCD);
    TH2D *h_theta_p_VS_phi_p_BPID_epCD = new TH2D("theta_p_VS_phi_p_BPID_epCD", "#theta_{p} vs #phi_{p} of CD proton (Before PID);#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_theta_p_VS_phi_p_BPID_epCD);
    TH1D *h_p_multiplicity_APID_epCD = new TH1D("p_multiplicity_APID_epCD", "Number of CD Protons in Event (After PID);CD Proton multiplicity;Counts", 10, 0, 10);
    HistoList.push_back(h_p_multiplicity_APID_epCD);
    TH1D *h_P_p_APID_epCD = new TH1D("P_p_APID_epCD", "CD Proton momentum (After PID);P_{p} [GeV/C];Counts", 50, 0., 3.);
    HistoList.push_back(h_P_p_APID_epCD);
    TH2D *h_theta_p_VS_phi_p_APID_epCD = new TH2D("theta_p_VS_phi_p_APID_epCD", "#theta_{p} vs #phi_{p} of CD proton (After PID);#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_theta_p_VS_phi_p_APID_epCD);

    TH1D *h_p_multiplicity_BPID_epFD = new TH1D("p_multiplicity_BPID_epFD", "Number of FD Protons in Event (Before PID);FD Proton multiplicity;Counts", 10, 0, 10);
    HistoList.push_back(h_p_multiplicity_BPID_epFD);
    TH1D *h_P_p_BPID_epFD = new TH1D("P_p_BPID_epFD", "FD Proton momentum (Before PID);P_{p} [GeV/C];Counts", 50, 0., 3.);
    HistoList.push_back(h_P_p_BPID_epFD);
    TH2D *h_theta_p_VS_phi_p_BPID_epFD = new TH2D("theta_p_VS_phi_p_BPID_epFD", "#theta_{p} vs #phi_{p} of FD proton (Before PID);#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_theta_p_VS_phi_p_BPID_epFD);
    TH1D *h_p_multiplicity_APID_epFD = new TH1D("p_multiplicity_APID_epFD", "Number of FD Protons in Even (After PID);FD Proton multiplicity;Counts", 10, 0, 10);
    HistoList.push_back(h_p_multiplicity_APID_epFD);
    TH1D *h_P_p_APID_epFD = new TH1D("P_p_APID_epFD", "FD Proton momentum (After PID);P_{p} [GeV/C];Counts", 50, 0., 3.);
    HistoList.push_back(h_P_p_APID_epFD);
    TH2D *h_theta_p_VS_phi_p_APID_epFD = new TH2D("theta_p_VS_phi_p_APID_epFD", "#theta_{p} vs #phi_{p} of FD proton (After PID);#phi_{p} [#circ];#theta_{p} [#circ]", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_theta_p_VS_phi_p_APID_epFD);

    TH2D *h_dbeta_p_BPID_epCD = new TH2D("dbeta_p_BPID_epCD", "#Delta #beta vs CD proton momentum (Before PID);P_{p} [GeV/c];#Delta#beta", 50, 0, 3, 50, -0.2, 0.2);
    hist_list_2.push_back(h_dbeta_p_BPID_epCD);
    TH1D *h_dVz_p_BPID_epCD = new TH1D("dVz_p_BPID_epCD", "Vertex correlation between CD proton and electron (Before PID);dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
    HistoList.push_back(h_dVz_p_BPID_epCD);
    TH1D *h_Chi2pid_p_BPID_epCD = new TH1D("Chi2pid_p_BPID_epCD", "CD Proton #chi^{2}_{p} (Before PID);#chi^{2}_{p};Counts", 50, -6, 6);
    HistoList.push_back(h_Chi2pid_p_BPID_epCD);
    TH2D *h_dbeta_p_APID_epCD = new TH2D("dbeta_p_APID_epCD", "#Delta #beta vs CD proton momentum (After PID);P_{p} [GeV/c];#Delta#beta", 50, 0, 3, 50, -0.2, 0.2);
    hist_list_2.push_back(h_dbeta_p_APID_epCD);
    TH1D *h_dVz_p_APID_epCD = new TH1D("dVz_p_APID_epCD", "Vertex correlation between CD proton and electron (After PID);dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
    HistoList.push_back(h_dVz_p_APID_epCD);
    TH1D *h_Chi2pid_p_APID_epCD = new TH1D("Chi2pid_p_APID_epCD", "CD Proton #chi^{2}_{p} (After PID);#chi^{2}_{p};Counts", 50, -6, 6);
    HistoList.push_back(h_Chi2pid_p_APID_epCD);

    TH2D *h_dbeta_p_BPID_epFD = new TH2D("dbeta_p_BPID_epFD", "#Delta #beta vs FD proton momentum (Before PID);P_{p} [GeV/c];#Delta#beta", 50, 0, 3, 50, -0.2, 0.2);
    hist_list_2.push_back(h_dbeta_p_BPID_epFD);
    TH1D *h_dVz_p_BPID_epFD = new TH1D("dVz_p_BPID_epFD", "Vertex correlation between FD proton and electron (Before PID);dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
    HistoList.push_back(h_dVz_p_BPID_epFD);
    TH1D *h_Chi2pid_p_BPID_epFD = new TH1D("Chi2pid_p_BPID_epFD", "FD Proton #chi^{2}_{p} (Before PID);#chi^{2}_{p};Counts", 50, -6, 6);
    HistoList.push_back(h_Chi2pid_p_BPID_epFD);
    TH2D *h_dbeta_p_APID_epFD = new TH2D("dbeta_p_APID_epFD", "#Delta #beta vs FD proton momentum (After PID);P_{p} [GeV/c];#Delta#beta", 50, 0, 3, 50, -0.2, 0.2);
    hist_list_2.push_back(h_dbeta_p_APID_epFD);
    TH1D *h_dVz_p_APID_epFD = new TH1D("dVz_p_APID_epFD", "Vertex correlation between FD proton and electron (After PID);dV^{p}_{z}=V^{p}_{z}-V^{e}_{z} [cm];Counts", 50, -8, 8);
    HistoList.push_back(h_dVz_p_APID_epFD);
    TH1D *h_Chi2pid_p_APID_epFD = new TH1D("Chi2pid_p_APID_epFD", "FD Proton #chi^{2}_{p} (After PID);#chi^{2}_{p};Counts", 50, -6, 6);
    HistoList.push_back(h_Chi2pid_p_APID_epFD);

    /* Missing variabels */
    TH1D *h_P_miss_BmissC_epCD = new TH1D("P_miss_BmissC_epCD", "Missing Momentum (Before M_{miss},P_{miss} Cuts);P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_BmissC_epCD);
    TH1D *h_theta_miss_BmissC_epCD = new TH1D("theta_miss_BmissC_epCD", "Missing Momentum (Before M_{miss},P_{miss} Cuts);theta_{miss} [#circ]", 48, -180, 180);
    HistoList.push_back(h_theta_miss_BmissC_epCD);
    TH2D *h_P_miss_VS_theta_miss_BmissC_epCD = new TH2D("P_miss_VS_theta_miss_BmissC_epCD", "Missing Momentum vs #theta_{miss} (Before M_{miss},P_{miss} Cuts);#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_BmissC_epCD);
    TH1D *h_P_miss_AmissC_epCD = new TH1D("P_miss_AmissC_epCD", "Missing Momentum (After M_{miss},P_{miss} Cuts);P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_AmissC_epCD);
    TH1D *h_theta_miss_AmissC_epCD = new TH1D("theta_miss_AmissC_epCD", "Missing Momentum (After M_{miss},P_{miss} Cuts);theta_{miss} [#circ]", 48, -180, 180);
    HistoList.push_back(h_theta_miss_AmissC_epCD);
    TH2D *h_P_miss_VS_theta_miss_AmissC_epCD = new TH2D("P_miss_VS_theta_miss_AmissC_epCD", "Missing Momentum vs #theta_{miss} (After M_{miss},P_{miss} Cuts);#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_AmissC_epCD);

    TH1D *h_P_miss_BmissC_epFD = new TH1D("P_miss_BmissC_epFD", "Missing Momentum (Before M_{miss},P_{miss} Cuts);P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_BmissC_epFD);
    TH1D *h_theta_miss_BmissC_epFD = new TH1D("theta_miss_BmissC_epFD", "Missing Momentum (Before M_{miss},P_{miss} Cuts);theta_{miss} [#circ]", 48, -180, 180);
    HistoList.push_back(h_theta_miss_BmissC_epFD);
    TH2D *h_P_miss_VS_theta_miss_BmissC_epFD = new TH2D("P_miss_VS_theta_miss_BmissC_epFD", "Missing Momentum vs #theta_{miss} (Before M_{miss},P_{miss} Cuts);#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_BmissC_epFD);
    TH1D *h_P_miss_AmissC_epFD = new TH1D("P_miss_AmissC_epFD", "Missing Momentum (After M_{miss},P_{miss} Cuts);P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_AmissC_epFD);
    TH1D *h_theta_miss_AmissC_epFD = new TH1D("theta_miss_AmissC_epFD", "Missing Momentum (After M_{miss},P_{miss} Cuts);theta_{miss} [#circ]", 48, -180, 180);
    HistoList.push_back(h_theta_miss_AmissC_epFD);
    TH2D *h_P_miss_VS_theta_miss_AmissC_epFD = new TH2D("P_miss_VS_theta_miss_AmissC_epFD", "Missing Momentum vs #theta_{miss} (After M_{miss},P_{miss} Cuts);#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_AmissC_epFD);

    TH1D *h_E_p_BmissC_epCD = new TH1D("E_p_BmissC_epCD", "CD Proton Energy (Before M_{miss},P_{miss} Cuts);E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_BmissC_epCD);
    TH1D *h_E_miss_BmissC_epCD = new TH1D("E_miss_BmissC_epCD", "Missing Energy (Before M_{miss},P_{miss} Cuts);E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_BmissC_epCD);
    TH1D *h_M_miss_BmissC_epCD = new TH1D("M_miss_BmissC_epCD", "Missing Mass (Before M_{miss},P_{miss} Cuts);M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_BmissC_epCD);
    TH1D *h_E_p_AmissC_epCD = new TH1D("E_p_AmissC_epCD", "CD Proton Energy (After M_{miss},P_{miss} Cuts);E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_AmissC_epCD);
    TH1D *h_E_miss_AmissC_epCD = new TH1D("E_miss_AmissC_epCD", "Missing Energy (After M_{miss},P_{miss} Cuts);E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_AmissC_epCD);
    TH1D *h_M_miss_AmissC_epCD = new TH1D("M_miss_AmissC_epCD", "Missing Mass (After M_{miss},P_{miss} Cuts);M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_AmissC_epCD);

    TH1D *h_E_p_BmissC_epFD = new TH1D("E_p_BmissC_epFD", "FD Proton Energy (Before M_{miss},P_{miss} Cuts);E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_BmissC_epFD);
    TH1D *h_E_miss_BmissC_epFD = new TH1D("E_miss_BmissC_epFD", "Missing Energy (Before M_{miss},P_{miss} Cuts);E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_BmissC_epFD);
    TH1D *h_M_miss_BmissC_epFD = new TH1D("M_miss_BmissC_epFD", "Missing Mass (Before M_{miss},P_{miss} Cuts);M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_BmissC_epFD);
    TH1D *h_E_p_AmissC_epFD = new TH1D("E_p_AmissC_epFD", "FD Proton Energy (After M_{miss},P_{miss} Cuts);E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_AmissC_epFD);
    TH1D *h_E_miss_AmissC_epFD = new TH1D("E_miss_AmissC_epFD", "Missing Energy (After M_{miss},P_{miss} Cuts);E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_AmissC_epFD);
    TH1D *h_M_miss_AmissC_epFD = new TH1D("M_miss_AmissC_epFD", "Missing Mass (After M_{miss},P_{miss} Cuts);M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_AmissC_epFD);

    /* Checks on which events have neutrons (Andrew) */
    TH1D *h_xB_BmissC_epCD = new TH1D("xB_BmissC_epCD", "x_{B} Distribution (Before M_{miss},P_{miss} Cuts);x_{B}", 100, 0.0, 2.0);
    HistoList.push_back(h_xB_BmissC_epCD);
    TH2D *h_xB_VS_M_miss_BmissC_epCD = new TH2D("xB_VS_M_miss_BmissC_epCD", "x_{B} vs. M_{miss} (Before M_{miss},P_{miss} Cuts);x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_BmissC_epCD);
    TH1D *h_xB_AmissC_epCD = new TH1D("xB_AmissC_epCD", "x_{B} Distribution (After M_{miss},P_{miss} Cuts);x_{B}", 100, 0.0, 2.0);
    HistoList.push_back(h_xB_AmissC_epCD);
    TH2D *h_xB_VS_M_miss_AmissC_epCD = new TH2D("xB_VS_M_miss_AmissC_epCD", "x_{B} vs. M_{miss} (After M_{miss},P_{miss} Cuts);x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_AmissC_epCD);

    TH1D *h_xB_BmissC_epFD = new TH1D("xB_BmissC_epFD", "x_{B} Distribution (Before M_{miss},P_{miss} Cuts);x_{B}", 100, 0.0, 2.0);
    HistoList.push_back(h_xB_BmissC_epFD);
    TH2D *h_xB_VS_M_miss_BmissC_epFD = new TH2D("xB_VS_M_miss_BmissC_epFD", "x_{B} vs. M_{miss} (Before M_{miss},P_{miss} Cuts);x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_BmissC_epFD);
    TH1D *h_xB_AmissC_epFD = new TH1D("xB_AmissC_epFD", "x_{B} Distribution (After M_{miss},P_{miss} Cuts);x_{B}", 100, 0.0, 2.0);
    HistoList.push_back(h_xB_AmissC_epFD);
    TH2D *h_xB_VS_M_miss_AmissC_epFD = new TH2D("xB_VS_M_miss_AmissC_epFD", "x_{B} vs. M_{miss} (After M_{miss},P_{miss} Cuts);x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_AmissC_epFD);

    TH2D *h_xB_VS_M_miss_epCDn = new TH2D("xB_VS_M_miss_epCDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_epCDn);
    TH2D *h_xB_VS_M_miss_epFDn = new TH2D("xB_VS_M_miss_epFDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_epFDn);

    TH2D *h_xB_VS_M_miss_goodN_epCDn = new TH2D("xB_VS_M_miss_goodN_epCDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_goodN_epCDn);
    TH2D *h_xB_VS_M_miss_badN_epCDn = new TH2D("xB_VS_M_miss_badN_epCDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_badN_epCDn);

    TH2D *h_xB_VS_M_miss_goodN_epFDn = new TH2D("xB_VS_M_miss_goodN_epFDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
    HistoList.push_back(h_xB_VS_M_miss_goodN_epFDn);
    TH2D *h_xB_VS_M_miss_badN_epFDn = new TH2D("xB_VS_M_miss_badN_epFDn", "x_{B} vs. M_{miss};x_{B};M_{miss} [GeV/c^{2}]", 100, 0.0, 2.0, 100, 0.5, 1.5);
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
    TH2D *h_P_n_VS_theta_n_epCDn = new TH2D("P_n_VS_theta_n_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
    HistoList.push_back(h_P_n_VS_theta_n_epCDn);

    TH1D *h_P_n_epFDn = new TH1D("P_n_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_epFDn);
    TH2D *h_P_n_VS_theta_n_epFDn = new TH2D("P_n_VS_theta_n_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
    HistoList.push_back(h_P_n_VS_theta_n_epFDn);

    TH1D *h_P_miss_epCDn = new TH1D("P_miss_epCDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_epCDn);
    TH2D *h_P_miss_VS_theta_miss_epCDn = new TH2D("P_miss_VS_theta_miss_epCDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_epCDn);

    TH1D *h_P_miss_epFDn = new TH1D("P_miss_epFDn", "Missing Momentum;P_{miss} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_miss_epFDn);
    TH2D *h_P_miss_VS_theta_miss_epFDn = new TH2D("P_miss_VS_theta_miss_epFDn", "Missing Momentum vs #theta_{miss};#theta_{miss} [#circ];P_{miss} [GeV/c]", 50, 0, 180, 50, 0, 1.5);
    HistoList.push_back(h_P_miss_VS_theta_miss_epFDn);

    TH1D *h_dpp_allN_epCDn = new TH1D("dpp_allN_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_allN_epCDn);
    TH1D *h_dpp_goodN_epCDn = new TH1D("dpp_goodN_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_goodN_epCDn);
    TH1D *h_dpp_badN_epCDn = new TH1D("dpp_badN_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_badN_epCDn);

    TH1D *h_dpp_allN_epFDn = new TH1D("dpp_allN_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_allN_epFDn);
    TH1D *h_dpp_goodN_epFDn = new TH1D("dpp_goodN_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_goodN_epFDn);
    TH1D *h_dpp_badN_epFDn = new TH1D("dpp_badN_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
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

    TH2D *h_dpp_VS_theta_n_miss_epCDn = new TH2D("dpp_VS_theta_n_miss_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_epCDn);

    TH2D *h_dpp_VS_theta_n_miss_epFDn = new TH2D("dpp_VS_theta_n_miss_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_epFDn);

    TH1D *h_E_p_epCDn = new TH1D("E_p_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_epCDn);
    TH1D *h_E_miss_epCDn = new TH1D("E_miss_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_epCDn);
    TH1D *h_M_miss_epCDn = new TH1D("M_miss_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_epCDn);
    TH2D *h_M_miss_VS_P_n_epCDn = new TH2D("M_miss_VS_P_n_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_epCDn);
    TH2D *h_M_miss_VS_theta_n_epCDn = new TH2D("M_miss_VS_theta_n_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_epCDn);
    TH2D *h_M_miss_VS_phi_n_epCDn = new TH2D("M_miss_VS_phi_n_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_epCDn);
    TH2D *h_M_miss_VS_P_miss_epCDn = new TH2D("M_miss_VS_P_miss_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_epCDn);
    TH2D *h_M_miss_VS_theta_miss_epCDn = new TH2D("M_miss_VS_theta_miss_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_epCDn);
    TH2D *h_M_miss_VS_phi_miss_epCDn = new TH2D("M_miss_VS_phi_miss_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_miss_epCDn);

    TH1D *h_E_p_epFDn = new TH1D("E_p_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_epFDn);
    TH1D *h_E_miss_epFDn = new TH1D("E_miss_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_epFDn);
    TH1D *h_M_miss_epFDn = new TH1D("M_miss_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_epFDn);
    TH2D *h_M_miss_VS_P_n_epFDn = new TH2D("M_miss_VS_P_n_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_epFDn);
    TH2D *h_M_miss_VS_theta_n_epFDn = new TH2D("M_miss_VS_theta_n_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_epFDn);
    TH2D *h_M_miss_VS_phi_n_epFDn = new TH2D("M_miss_VS_phi_n_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_epFDn);
    TH2D *h_M_miss_VS_P_miss_epFDn = new TH2D("M_miss_VS_P_miss_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_epFDn);
    TH2D *h_M_miss_VS_theta_miss_epFDn = new TH2D("M_miss_VS_theta_miss_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_epFDn);
    TH2D *h_M_miss_VS_phi_miss_epFDn = new TH2D("M_miss_VS_phi_miss_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND_epCDn = new TH2D("M_miss_VS_Edep_CND_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND_epCDn);
    TH2D *h_path_VS_Edep_CND_epCDn = new TH2D("path_VS_Edep_CND_epCDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_epCDn = new TH2D("theta_n_miss_VS_Edep_CND_epCDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_epCDn);
    TH2D *h_ToF_VS_Edep_CND_epCDn = new TH2D("ToF_VS_Edep_CND_epCDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_epCDn);
    TH2D *h_nSector_VS_Edep_CND_epCDn = new TH2D("nSector_VS_Edep_CND_epCDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_epCDn);

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
    TH2D *h_M_miss_VS_Edep_CND_epFDn = new TH2D("M_miss_VS_Edep_CND_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND_epFDn);
    TH2D *h_path_VS_Edep_CND_epFDn = new TH2D("path_VS_Edep_CND_epFDn", "Path length vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND_epFDn = new TH2D("theta_n_miss_VS_Edep_CND_epFDn", "#theta_{n,miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND_epFDn);
    TH2D *h_ToF_VS_Edep_CND_epFDn = new TH2D("ToF_VS_Edep_CND_epFDn", "Neutron ToF vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND_epFDn);
    TH2D *h_nSector_VS_Edep_CND_epFDn = new TH2D("nSector_VS_Edep_CND_epFDn", "Neutron Sector Number vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND_epFDn);

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
    TH2D *h_M_miss_VS_Edep_CTOF_epCDn = new TH2D("M_miss_VS_Edep_CTOF_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_epCDn);
    TH2D *h_path_VS_Edep_CTOF_epCDn = new TH2D("path_VS_Edep_CTOF_epCDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_epCDn = new TH2D("theta_n_miss_VS_Edep_CTOF_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_epCDn);
    TH2D *h_ToF_VS_Edep_CTOF_epCDn = new TH2D("ToF_VS_Edep_CTOF_epCDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_epCDn);
    TH2D *h_nSector_VS_Edep_CTOF_epCDn = new TH2D("nSector_VS_Edep_CTOF_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_epCDn);

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
    TH2D *h_M_miss_VS_Edep_CTOF_epFDn = new TH2D("M_miss_VS_Edep_CTOF_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_epFDn);
    TH2D *h_path_VS_Edep_CTOF_epFDn = new TH2D("path_VS_Edep_CTOF_epFDn", "Path length vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CTOF_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CTOF_epFDn = new TH2D("theta_n_miss_VS_Edep_CTOF_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CTOF_epFDn);
    TH2D *h_ToF_VS_Edep_CTOF_epFDn = new TH2D("ToF_VS_Edep_CTOF_epFDn", "Neutron ToF vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CTOF_epFDn);
    TH2D *h_nSector_VS_Edep_CTOF_epFDn = new TH2D("nSector_VS_Edep_CTOF_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CTOF_epFDn);

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
    TH2D *h_M_miss_VS_Edep_single_epCDn = new TH2D("M_miss_VS_Edep_single_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_single_epFDn = new TH2D("M_miss_VS_Edep_single_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND1_epCDn = new TH2D("M_miss_VS_Edep_CND1_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_epCDn);
    TH2D *h_path_VS_Edep_CND1_epCDn = new TH2D("path_VS_Edep_CND1_epCDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_epCDn = new TH2D("theta_n_miss_VS_Edep_CND1_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_epCDn);
    TH2D *h_ToF_VS_Edep_CND1_epCDn = new TH2D("ToF_VS_Edep_CND1_epCDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_epCDn);
    TH2D *h_nSector_VS_Edep_CND1_epCDn = new TH2D("nSector_VS_Edep_CND1_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_epCDn);

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
    TH2D *h_M_miss_VS_Edep_CND1_epFDn = new TH2D("M_miss_VS_Edep_CND1_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_epFDn);
    TH2D *h_path_VS_Edep_CND1_epFDn = new TH2D("path_VS_Edep_CND1_epFDn", "Path length vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND1_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND1_epFDn = new TH2D("theta_n_miss_VS_Edep_CND1_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND1_epFDn);
    TH2D *h_ToF_VS_Edep_CND1_epFDn = new TH2D("ToF_VS_Edep_CND1_epFDn", "Neutron ToF vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND1_epFDn);
    TH2D *h_nSector_VS_Edep_CND1_epFDn = new TH2D("nSector_VS_Edep_CND1_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND1_epFDn);

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
    TH2D *h_M_miss_VS_Edep_CND2_epCDn = new TH2D("M_miss_VS_Edep_CND2_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_epCDn);
    TH2D *h_path_VS_Edep_CND2_epCDn = new TH2D("path_VS_Edep_CND2_epCDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_epCDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_epCDn = new TH2D("theta_n_miss_VS_Edep_CND2_epCDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_epCDn);
    TH2D *h_ToF_VS_Edep_CND2_epCDn = new TH2D("ToF_VS_Edep_CND2_epCDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_epCDn);
    TH2D *h_nSector_VS_Edep_CND2_epCDn = new TH2D("nSector_VS_Edep_CND2_epCDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_epCDn);

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
    TH2D *h_M_miss_VS_Edep_CND2_epFDn = new TH2D("M_miss_VS_Edep_CND2_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_epFDn);
    TH2D *h_path_VS_Edep_CND2_epFDn = new TH2D("path_VS_Edep_CND2_epFDn", "Path length vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND2_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND2_epFDn = new TH2D("theta_n_miss_VS_Edep_CND2_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND2_epFDn);
    TH2D *h_ToF_VS_Edep_CND2_epFDn = new TH2D("ToF_VS_Edep_CND2_epFDn", "Neutron ToF vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND2_epFDn);
    TH2D *h_nSector_VS_Edep_CND2_epFDn = new TH2D("nSector_VS_Edep_CND2_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND2_epFDn);

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
    TH2D *h_M_miss_VS_Edep_CND3_epCDn = new TH2D("M_miss_VS_Edep_CND3_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND3_epFDn = new TH2D("M_miss_VS_Edep_CND3_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_epFDn);
    TH2D *h_path_VS_Edep_CND3_epFDn = new TH2D("path_VS_Edep_CND3_epFDn", "Path length vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Path length [cm]", 50, 0, 100, 50, 0., 100.);
    HistoList.push_back(h_path_VS_Edep_CND3_epFDn);
    TH2D *h_theta_n_miss_VS_Edep_CND3_epFDn = new TH2D("theta_n_miss_VS_Edep_CND3_epFDn", "#theta_{n,miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];#theta_{n,miss} [#circ]", 50, 0, 100, 50, 0., 180.);
    HistoList.push_back(h_theta_n_miss_VS_Edep_CND3_epFDn);
    TH2D *h_ToF_VS_Edep_CND3_epFDn = new TH2D("ToF_VS_Edep_CND3_epFDn", "Neutron ToF vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];ToF [ns]", 50, 0, 100, 50, 0., 20.);
    HistoList.push_back(h_ToF_VS_Edep_CND3_epFDn);
    TH2D *h_nSector_VS_Edep_CND3_epFDn = new TH2D("nSector_VS_Edep_CND3_epFDn", "Neutron Sector Number vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];Sector Number", 50, 0, 100, 24, 0., 24.);
    HistoList.push_back(h_nSector_VS_Edep_CND3_epFDn);

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
    TH2D *h_M_miss_VS_ToF_epCDn = new TH2D("M_miss_VS_ToF_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_ToF_epFDn = new TH2D("M_miss_VS_ToF_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
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
    hist_list_1.push_back(h_n_multiplicity_allN_epCDn_Step0);
    TH1D *h_n_multiplicity_goodN_epCDn_Step0 = new TH1D("n_multiplicity_goodN_epCDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epCDn_Step0);
    TH1D *h_n_multiplicity_badN_epCDn_Step0 = new TH1D("n_multiplicity_badN_epCDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epCDn_Step0);

    TH1D *h_n_multiplicity_allN_epFDn_Step0 = new TH1D("n_multiplicity_allN_epFDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epFDn_Step0);
    TH1D *h_n_multiplicity_goodN_epFDn_Step0 = new TH1D("n_multiplicity_goodN_epFDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epFDn_Step0);
    TH1D *h_n_multiplicity_badN_epFDn_Step0 = new TH1D("n_multiplicity_badN_epFDn_Step0", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epFDn_Step0);

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
    TH2D *h_P_n_VS_theta_n_goodN_Step0_epCDn = new TH2D("P_n_VS_theta_n_goodN_Step0_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
    HistoList.push_back(h_P_n_VS_theta_n_goodN_Step0_epCDn);
    TH2D *h_P_n_VS_theta_n_badN_Step0_epCDn = new TH2D("P_n_VS_theta_n_badN_Step0_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
    HistoList.push_back(h_P_n_VS_theta_n_badN_Step0_epCDn);

    TH1D *h_P_n_goodN_Step0_epFDn = new TH1D("P_n_goodN_Step0_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_goodN_Step0_epFDn);
    TH1D *h_P_n_badN_Step0_epFDn = new TH1D("P_n_badN_Step0_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_badN_Step0_epFDn);
    TH2D *h_P_n_VS_theta_n_goodN_Step0_epFDn = new TH2D("P_n_VS_theta_n_goodN_Step0_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
    HistoList.push_back(h_P_n_VS_theta_n_goodN_Step0_epFDn);
    TH2D *h_P_n_VS_theta_n_badN_Step0_epFDn = new TH2D("P_n_VS_theta_n_badN_Step0_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
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

    TH1D *h_dpp_allN_Step0_epCDn = new TH1D("dpp_allN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_allN_Step0_epCDn);
    TH1D *h_dpp_goodN_Step0_epCDn = new TH1D("dpp_goodN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_goodN_Step0_epCDn);
    TH1D *h_dpp_badN_Step0_epCDn = new TH1D("dpp_badN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_badN_Step0_epCDn);

    TH1D *h_dpp_allN_Step0_epFDn = new TH1D("dpp_allN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_allN_Step0_epFDn);
    TH1D *h_dpp_goodN_Step0_epFDn = new TH1D("dpp_goodN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_goodN_Step0_epFDn);
    TH1D *h_dpp_badN_Step0_epFDn = new TH1D("dpp_badN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
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

    TH2D *h_dpp_VS_theta_n_miss_allN_Step0_epCDn = new TH2D("dpp_VS_theta_n_miss_allN_Step0_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_allN_Step0_epCDn);

    TH2D *h_dpp_VS_theta_n_miss_allN_Step0_epFDn = new TH2D("dpp_VS_theta_n_miss_allN_Step0_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_allN_Step0_epFDn);

    TH1D *h_E_p_goodN_Step0_epCDn = new TH1D("E_p_goodN_Step0_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_goodN_Step0_epCDn);
    TH1D *h_E_p_badN_Step0_epCDn = new TH1D("E_p_badN_Step0_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_badN_Step0_epCDn);
    TH1D *h_E_miss_goodN_Step0_epCDn = new TH1D("E_miss_goodN_Step0_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_goodN_Step0_epCDn);
    TH1D *h_E_miss_badN_Step0_epCDn = new TH1D("E_miss_badN_Step0_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_badN_Step0_epCDn);
    TH1D *h_M_miss_goodN_Step0_epCDn = new TH1D("M_miss_goodN_Step0_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_goodN_Step0_epCDn);
    TH1D *h_M_miss_badN_Step0_epCDn = new TH1D("M_miss_badN_Step0_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_P_n_goodN_Step0_epCDn = new TH2D("M_miss_VS_P_n_goodN_Step0_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_P_n_badN_Step0_epCDn = new TH2D("M_miss_VS_P_n_badN_Step0_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_theta_n_goodN_Step0_epCDn = new TH2D("M_miss_VS_theta_n_goodN_Step0_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_theta_n_badN_Step0_epCDn = new TH2D("M_miss_VS_theta_n_badN_Step0_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_phi_n_goodN_Step0_epCDn = new TH2D("M_miss_VS_phi_n_goodN_Step0_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_phi_n_badN_Step0_epCDn = new TH2D("M_miss_VS_phi_n_badN_Step0_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_P_miss_goodN_Step0_epCDn = new TH2D("M_miss_VS_P_miss_goodN_Step0_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_P_miss_badN_Step0_epCDn = new TH2D("M_miss_VS_P_miss_badN_Step0_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_theta_miss_goodN_Step0_epCDn = new TH2D("M_miss_VS_theta_miss_goodN_Step0_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_theta_miss_badN_Step0_epCDn = new TH2D("M_miss_VS_theta_miss_badN_Step0_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_badN_Step0_epCDn);
    TH2D *h_M_miss_VS_phi_miss_goodN_Step0_epCDn = new TH2D("M_miss_VS_phi_miss_goodN_Step0_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_miss_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_phi_miss_badN_Step0_epCDn = new TH2D("M_miss_VS_phi_miss_badN_Step0_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_miss_badN_Step0_epCDn);

    TH1D *h_E_p_goodN_Step0_epFDn = new TH1D("E_p_goodN_Step0_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_goodN_Step0_epFDn);
    TH1D *h_E_p_badN_Step0_epFDn = new TH1D("E_p_badN_Step0_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_badN_Step0_epFDn);
    TH1D *h_E_miss_goodN_Step0_epFDn = new TH1D("E_miss_goodN_Step0_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_goodN_Step0_epFDn);
    TH1D *h_E_miss_badN_Step0_epFDn = new TH1D("E_miss_badN_Step0_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_badN_Step0_epFDn);
    TH1D *h_M_miss_goodN_Step0_epFDn = new TH1D("M_miss_goodN_Step0_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_goodN_Step0_epFDn);
    TH1D *h_M_miss_badN_Step0_epFDn = new TH1D("M_miss_badN_Step0_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_P_n_goodN_Step0_epFDn = new TH2D("M_miss_VS_P_n_goodN_Step0_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_P_n_badN_Step0_epFDn = new TH2D("M_miss_VS_P_n_badN_Step0_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_theta_n_goodN_Step0_epFDn = new TH2D("M_miss_VS_theta_n_goodN_Step0_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_theta_n_badN_Step0_epFDn = new TH2D("M_miss_VS_theta_n_badN_Step0_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_phi_n_goodN_Step0_epFDn = new TH2D("M_miss_VS_phi_n_goodN_Step0_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_phi_n_badN_Step0_epFDn = new TH2D("M_miss_VS_phi_n_badN_Step0_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_P_miss_goodN_Step0_epFDn = new TH2D("M_miss_VS_P_miss_goodN_Step0_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_P_miss_badN_Step0_epFDn = new TH2D("M_miss_VS_P_miss_badN_Step0_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_theta_miss_goodN_Step0_epFDn = new TH2D("M_miss_VS_theta_miss_goodN_Step0_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_theta_miss_badN_Step0_epFDn = new TH2D("M_miss_VS_theta_miss_badN_Step0_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_badN_Step0_epFDn);
    TH2D *h_M_miss_VS_phi_miss_goodN_Step0_epFDn = new TH2D("M_miss_VS_phi_miss_goodN_Step0_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_miss_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_phi_miss_badN_Step0_epFDn = new TH2D("M_miss_VS_phi_miss_badN_Step0_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND_goodN_Step0_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND_badN_Step0_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND_goodN_Step0_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND_badN_Step0_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CTOF_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CTOF_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CTOF_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CTOF_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_single_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_single_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_single_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_single_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_single_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_single_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_single_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_single_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_single_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_single_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND1_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND1_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND1_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND1_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND2_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND2_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND2_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND2_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND3_goodN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step0_epCDn = new TH2D("M_miss_VS_Edep_CND3_badN_Step0_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND3_goodN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step0_epFDn = new TH2D("M_miss_VS_Edep_CND3_badN_Step0_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_ToF_goodN_Step0_epCDn = new TH2D("M_miss_VS_ToF_goodN_Step0_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_ToF_goodN_Step0_epCDn);
    TH2D *h_M_miss_VS_ToF_badN_Step0_epCDn = new TH2D("M_miss_VS_ToF_badN_Step0_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_ToF_goodN_Step0_epFDn = new TH2D("M_miss_VS_ToF_goodN_Step0_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_ToF_goodN_Step0_epFDn);
    TH2D *h_M_miss_VS_ToF_badN_Step0_epFDn = new TH2D("M_miss_VS_ToF_badN_Step0_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
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
        sprintf(temp_name_A, "sdiff_pos_goodN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step0_layer_epCDn[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name_A, "sdiff_pos_badN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step0_layer_epCDn[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name_A, "sdiff_pos_goodN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step0_layer_epFDn[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name_A, "sdiff_pos_badN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step0_layer_epFDn[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step0_layer_epFDn[k]);

        sprintf(temp_name_A, "sdiff_pos_mom_goodN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step0_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name_A, "sdiff_pos_mom_badN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step0_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name_A, "sdiff_pos_mom_goodN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step0_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name_A, "sdiff_pos_mom_badN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step0_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step0_layer_epFDn[k]);

        sprintf(temp_name_A, "sdiff_pos_z_goodN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step0_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name_A, "sdiff_pos_z_badN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step0_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name_A, "sdiff_pos_z_goodN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step0_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name_A, "sdiff_pos_z_badN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step0_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_badN_Step0_layer_epFDn[k]);

        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_goodN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_badN_Step0_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_goodN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_badN_Step0_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step0_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
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
        sprintf(temp_name_A, "diff_ToFc_z_goodN_Step0_layer_%d_epCDn", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step0_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_goodN_Step0_layer_epCDn[k]);
        sprintf(temp_name_A, "diff_ToFc_z_badN_Step0_layer_%d_epCDn", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step0_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_badN_Step0_layer_epCDn[k]);

        sprintf(temp_name_A, "diff_ToFc_z_goodN_Step0_layer_%d_epFDn", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step0_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_goodN_Step0_layer_epFDn[k]);
        sprintf(temp_name_A, "diff_ToFc_z_badN_Step0_layer_%d_epFDn", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step0_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
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

    TH2D *h_numberNearby_momN_goodN_Step0_epCDn = new TH2D("numberNearby_momN_goodN_Step0_epCDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    HistoList.push_back(h_numberNearby_momN_goodN_Step0_epCDn);
    TH2D *h_numberNearby_momN_badN_Step0_epCDn = new TH2D("numberNearby_momN_badN_Step0_epCDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    HistoList.push_back(h_numberNearby_momN_badN_Step0_epCDn);

    TH2D *h_numberNearby_momN_goodN_Step0_epFDn = new TH2D("numberNearby_momN_goodN_Step0_epFDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    HistoList.push_back(h_numberNearby_momN_goodN_Step0_epFDn);
    TH2D *h_numberNearby_momN_badN_Step0_epFDn = new TH2D("numberNearby_momN_badN_Step0_epFDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
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

    TH1D *h_phidiff_en_goodN_Step0_epCDn = new TH1D("phidiff_en_goodN_Step0_epCDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
    HistoList.push_back(h_phidiff_en_goodN_Step0_epCDn);
    TH1D *h_phidiff_en_badN_Step0_epCDn = new TH1D("phidiff_en_badN_Step0_epCDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
    HistoList.push_back(h_phidiff_en_badN_Step0_epCDn);

    TH1D *h_phidiff_en_goodN_Step0_epFDn = new TH1D("phidiff_en_goodN_Step0_epFDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
    HistoList.push_back(h_phidiff_en_goodN_Step0_epFDn);
    TH1D *h_phidiff_en_badN_Step0_epFDn = new TH1D("phidiff_en_badN_Step0_epFDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
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
    hist_list_1.push_back(h_n_multiplicity_allN_epCDn_Step1);
    TH1D *h_n_multiplicity_goodN_epCDn_Step1 = new TH1D("n_multiplicity_goodN_epCDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epCDn_Step1);
    TH1D *h_n_multiplicity_badN_epCDn_Step1 = new TH1D("n_multiplicity_badN_epCDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epCDn_Step1);

    TH1D *h_n_multiplicity_allN_epFDn_Step1 = new TH1D("n_multiplicity_allN_epFDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epFDn_Step1);
    TH1D *h_n_multiplicity_goodN_epFDn_Step1 = new TH1D("n_multiplicity_goodN_epFDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epFDn_Step1);
    TH1D *h_n_multiplicity_badN_epFDn_Step1 = new TH1D("n_multiplicity_badN_epFDn_Step1", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epFDn_Step1);

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
    TH2D *h_P_n_VS_theta_n_goodN_Step1_epCDn = new TH2D("P_n_VS_theta_n_goodN_Step1_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
    HistoList.push_back(h_P_n_VS_theta_n_goodN_Step1_epCDn);
    TH2D *h_P_n_VS_theta_n_badN_Step1_epCDn = new TH2D("P_n_VS_theta_n_badN_Step1_epCDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
    HistoList.push_back(h_P_n_VS_theta_n_badN_Step1_epCDn);

    TH1D *h_P_n_goodN_Step1_epFDn = new TH1D("P_n_goodN_Step1_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_goodN_Step1_epFDn);
    TH1D *h_P_n_badN_Step1_epFDn = new TH1D("P_n_badN_Step1_epFDn", "Neutron Momentum;P_{n} [GeV/c]", 50, 0, 1.5);
    HistoList.push_back(h_P_n_badN_Step1_epFDn);
    TH2D *h_P_n_VS_theta_n_goodN_Step1_epFDn = new TH2D("P_n_VS_theta_n_goodN_Step1_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
    HistoList.push_back(h_P_n_VS_theta_n_goodN_Step1_epFDn);
    TH2D *h_P_n_VS_theta_n_badN_Step1_epFDn = new TH2D("P_n_VS_theta_n_badN_Step1_epFDn", "Neutron Momentum vs Theta;#theta_{n} [#circ];P_{n} [GeV/c]", 55, 35, 145, 50, 0, 1.2);
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

    TH1D *h_dpp_allN_Step1_epCDn = new TH1D("dpp_allN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_allN_Step1_epCDn);
    TH1D *h_dpp_goodN_Step1_epCDn = new TH1D("dpp_goodN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_goodN_Step1_epCDn);
    TH1D *h_dpp_badN_Step1_epCDn = new TH1D("dpp_badN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_badN_Step1_epCDn);

    TH1D *h_dpp_allN_Step1_epFDn = new TH1D("dpp_allN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_allN_Step1_epFDn);
    TH1D *h_dpp_goodN_Step1_epFDn = new TH1D("dpp_goodN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
    HistoList.push_back(h_dpp_goodN_Step1_epFDn);
    TH1D *h_dpp_badN_Step1_epFDn = new TH1D("dpp_badN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} Distribution;|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss}", 50, 0, 180);
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

    TH2D *h_dpp_VS_theta_n_miss_allN_Step1_epCDn = new TH2D("dpp_VS_theta_n_miss_allN_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_allN_Step1_epCDn);

    TH2D *h_dpp_VS_theta_n_miss_allN_Step1_epFDn = new TH2D("dpp_VS_theta_n_miss_allN_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_allN_Step1_epFDn);

    TH1D *h_E_p_goodN_Step1_epCDn = new TH1D("E_p_goodN_Step1_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_goodN_Step1_epCDn);
    TH1D *h_E_p_badN_Step1_epCDn = new TH1D("E_p_badN_Step1_epCDn", "CD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_badN_Step1_epCDn);
    TH1D *h_E_miss_goodN_Step1_epCDn = new TH1D("E_miss_goodN_Step1_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_goodN_Step1_epCDn);
    TH1D *h_E_miss_badN_Step1_epCDn = new TH1D("E_miss_badN_Step1_epCDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_badN_Step1_epCDn);
    TH1D *h_M_miss_goodN_Step1_epCDn = new TH1D("M_miss_goodN_Step1_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_goodN_Step1_epCDn);
    TH1D *h_M_miss_badN_Step1_epCDn = new TH1D("M_miss_badN_Step1_epCDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_P_n_goodN_Step1_epCDn = new TH2D("M_miss_VS_P_n_goodN_Step1_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_P_n_badN_Step1_epCDn = new TH2D("M_miss_VS_P_n_badN_Step1_epCDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_theta_n_goodN_Step1_epCDn = new TH2D("M_miss_VS_theta_n_goodN_Step1_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_theta_n_badN_Step1_epCDn = new TH2D("M_miss_VS_theta_n_badN_Step1_epCDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_phi_n_goodN_Step1_epCDn = new TH2D("M_miss_VS_phi_n_goodN_Step1_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_phi_n_badN_Step1_epCDn = new TH2D("M_miss_VS_phi_n_badN_Step1_epCDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_P_miss_goodN_Step1_epCDn = new TH2D("M_miss_VS_P_miss_goodN_Step1_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_P_miss_badN_Step1_epCDn = new TH2D("M_miss_VS_P_miss_badN_Step1_epCDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_theta_miss_goodN_Step1_epCDn = new TH2D("M_miss_VS_theta_miss_goodN_Step1_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_theta_miss_badN_Step1_epCDn = new TH2D("M_miss_VS_theta_miss_badN_Step1_epCDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_badN_Step1_epCDn);
    TH2D *h_M_miss_VS_phi_miss_goodN_Step1_epCDn = new TH2D("M_miss_VS_phi_miss_goodN_Step1_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_miss_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_phi_miss_badN_Step1_epCDn = new TH2D("M_miss_VS_phi_miss_badN_Step1_epCDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_miss_badN_Step1_epCDn);

    TH1D *h_E_p_goodN_Step1_epFDn = new TH1D("E_p_goodN_Step1_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_goodN_Step1_epFDn);
    TH1D *h_E_p_badN_Step1_epFDn = new TH1D("E_p_badN_Step1_epFDn", "FD Proton Energy;E_{p} [GeV]", 50, 0, 3.);
    HistoList.push_back(h_E_p_badN_Step1_epFDn);
    TH1D *h_E_miss_goodN_Step1_epFDn = new TH1D("E_miss_goodN_Step1_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_goodN_Step1_epFDn);
    TH1D *h_E_miss_badN_Step1_epFDn = new TH1D("E_miss_badN_Step1_epFDn", "Missing Energy;E_{miss} [GeV]", 50, 0.5, 1.5);
    HistoList.push_back(h_E_miss_badN_Step1_epFDn);
    TH1D *h_M_miss_goodN_Step1_epFDn = new TH1D("M_miss_goodN_Step1_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_goodN_Step1_epFDn);
    TH1D *h_M_miss_badN_Step1_epFDn = new TH1D("M_miss_badN_Step1_epFDn", "Missing Mass;M_{miss} [GeV/c^{2}]", 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_P_n_goodN_Step1_epFDn = new TH2D("M_miss_VS_P_n_goodN_Step1_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_P_n_badN_Step1_epFDn = new TH2D("M_miss_VS_P_n_badN_Step1_epFDn", "Missing Mass vs Measured Neutron Momentum;P_{n} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0.25, 1., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_n_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_theta_n_goodN_Step1_epFDn = new TH2D("M_miss_VS_theta_n_goodN_Step1_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_theta_n_badN_Step1_epFDn = new TH2D("M_miss_VS_theta_n_badN_Step1_epFDn", "Missing Mass vs Measured #theta_{n};#theta_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_n_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_phi_n_goodN_Step1_epFDn = new TH2D("M_miss_VS_phi_n_goodN_Step1_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_phi_n_badN_Step1_epFDn = new TH2D("M_miss_VS_phi_n_badN_Step1_epFDn", "Missing Mass vs Measured #phi_{n};#phi_{n} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_n_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_P_miss_goodN_Step1_epFDn = new TH2D("M_miss_VS_P_miss_goodN_Step1_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_P_miss_badN_Step1_epFDn = new TH2D("M_miss_VS_P_miss_badN_Step1_epFDn", "Missing Mass vs Missing Momentum;P_{miss} [GeV/c];M_{miss} [GeV/c^{2}]", 50, 0, 1.25, 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_P_miss_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_theta_miss_goodN_Step1_epFDn = new TH2D("M_miss_VS_theta_miss_goodN_Step1_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_theta_miss_badN_Step1_epFDn = new TH2D("M_miss_VS_theta_miss_badN_Step1_epFDn", "Missing Mass vs #theta_{miss};#theta_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, 0., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_theta_miss_badN_Step1_epFDn);
    TH2D *h_M_miss_VS_phi_miss_goodN_Step1_epFDn = new TH2D("M_miss_VS_phi_miss_goodN_Step1_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
    HistoList.push_back(h_M_miss_VS_phi_miss_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_phi_miss_badN_Step1_epFDn = new TH2D("M_miss_VS_phi_miss_badN_Step1_epFDn", "Missing Mass vs #phi_{miss};#phi_{miss} [#circ];M_{miss} [GeV/c^{2}]", 50, -180., 180., 50, 0, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND_goodN_Step1_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND_badN_Step1_epCDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND_goodN_Step1_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND_badN_Step1_epFDn", "M_{miss} vs Total Neutron Energy Deposition in the CND;E^{CND}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CTOF_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CTOF_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CTOF_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CTOF_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CTOF_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CTOF_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CTOF_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the CTOF;E^{CTOF}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_single_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_single_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_single_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_single_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_single_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_single_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_single_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_single_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_single_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_single_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in a single CND layer;E^{CND,i}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND1_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND1_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND1_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND1_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND1_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND1_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND1_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the first CND layer;E^{CND,1}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND2_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND2_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND2_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND2_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND2_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND2_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND2_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the second CND layer;E^{CND,2}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND3_goodN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step1_epCDn = new TH2D("M_miss_VS_Edep_CND3_badN_Step1_epCDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_Edep_CND3_goodN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND3_goodN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_Edep_CND3_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_Edep_CND3_badN_Step1_epFDn = new TH2D("M_miss_VS_Edep_CND3_badN_Step1_epFDn", "M_{miss} vs Neutron Energy Deposition in the third CND layer;E^{CND,3}_{dep} [MeV];M_{miss} [GeV/c^{2}]", 50, 0, 100, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_ToF_goodN_Step1_epCDn = new TH2D("M_miss_VS_ToF_goodN_Step1_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_ToF_goodN_Step1_epCDn);
    TH2D *h_M_miss_VS_ToF_badN_Step1_epCDn = new TH2D("M_miss_VS_ToF_badN_Step1_epCDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
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
    TH2D *h_M_miss_VS_ToF_goodN_Step1_epFDn = new TH2D("M_miss_VS_ToF_goodN_Step1_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
    HistoList.push_back(h_M_miss_VS_ToF_goodN_Step1_epFDn);
    TH2D *h_M_miss_VS_ToF_badN_Step1_epFDn = new TH2D("M_miss_VS_ToF_badN_Step1_epFDn", "M_{miss} vs ToF;ToF [ns];M_{miss} [GeV/c^{2}]", 50, 0, 20, 50, 0.5, 1.5);
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

    TH2D *h_dpp_VS_theta_n_miss_Step1_epCDn = new TH2D("dpp_VS_theta_n_miss_Step1_epCDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_Step1_epCDn);

    TH2D *h_dpp_VS_theta_n_miss_Step1_epFDn = new TH2D("dpp_VS_theta_n_miss_Step1_epFDn", "|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss} vs. #theta_{n,miss};|#vec{P}_{miss}|-|#vec{P}_{n}|/P_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    HistoList.push_back(h_dpp_VS_theta_n_miss_Step1_epFDn);

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
        sprintf(temp_name_A, "sdiff_pos_goodN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step1_layer_epCDn[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name_A, "sdiff_pos_badN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step1_layer_epCDn[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name_A, "sdiff_pos_goodN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step1_layer_epFDn[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name_A, "sdiff_pos_badN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step1_layer_epFDn[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        HistoList.push_back(h_sdiff_pos_badN_Step1_layer_epFDn[k]);

        sprintf(temp_name_A, "sdiff_pos_mom_goodN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step1_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name_A, "sdiff_pos_mom_badN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step1_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name_A, "sdiff_pos_mom_goodN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step1_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name_A, "sdiff_pos_mom_badN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step1_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        HistoList.push_back(h_sdiff_pos_mom_badN_Step1_layer_epFDn[k]);

        sprintf(temp_name_A, "sdiff_pos_z_goodN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step1_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name_A, "sdiff_pos_z_badN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step1_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name_A, "sdiff_pos_z_goodN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step1_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name_A, "sdiff_pos_z_badN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step1_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        HistoList.push_back(h_sdiff_pos_z_badN_Step1_layer_epFDn[k]);

        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_goodN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_badN_Step1_layer_%d_epCDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_goodN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
        HistoList.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_badN_Step1_layer_%d_epFDn", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
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
        sprintf(temp_name_A, "diff_ToFc_z_goodN_Step1_layer_%d_epCDn", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step1_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_goodN_Step1_layer_epCDn[k]);
        sprintf(temp_name_A, "diff_ToFc_z_badN_Step1_layer_%d_epCDn", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step1_layer_epCDn[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_badN_Step1_layer_epCDn[k]);

        sprintf(temp_name_A, "diff_ToFc_z_goodN_Step1_layer_%d_epFDn", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step1_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        HistoList.push_back(h_diff_ToFc_z_Edep_goodN_Step1_layer_epFDn[k]);
        sprintf(temp_name_A, "diff_ToFc_z_badN_Step1_layer_%d_epFDn", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step1_layer_epFDn[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
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

    TH2D *h_numberNearby_momN_goodN_Step1_epCDn = new TH2D("numberNearby_momN_goodN_Step1_epCDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    HistoList.push_back(h_numberNearby_momN_goodN_Step1_epCDn);
    TH2D *h_numberNearby_momN_badN_Step1_epCDn = new TH2D("numberNearby_momN_badN_Step1_epCDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    HistoList.push_back(h_numberNearby_momN_badN_Step1_epCDn);

    TH2D *h_numberNearby_momN_goodN_Step1_epFDn = new TH2D("numberNearby_momN_goodN_Step1_epFDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    HistoList.push_back(h_numberNearby_momN_goodN_Step1_epFDn);
    TH2D *h_numberNearby_momN_badN_Step1_epFDn = new TH2D("numberNearby_momN_badN_Step1_epFDn", "Number of Nearby Hits vs. P_{n} for CND Neutrons;# Hits;P_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
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

    TH1D *h_phidiff_en_goodN_Step1_epCDn = new TH1D("phidiff_en_goodN_Step1_epCDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
    HistoList.push_back(h_phidiff_en_goodN_Step1_epCDn);
    TH1D *h_phidiff_en_badN_Step1_epCDn = new TH1D("phidiff_en_badN_Step1_epCDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
    HistoList.push_back(h_phidiff_en_badN_Step1_epCDn);

    TH1D *h_phidiff_en_goodN_Step1_epFDn = new TH1D("phidiff_en_goodN_Step1_epFDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
    HistoList.push_back(h_phidiff_en_goodN_Step1_epFDn);
    TH1D *h_phidiff_en_badN_Step1_epFDn = new TH1D("phidiff_en_badN_Step1_epFDn", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
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
    hist_list_1.push_back(h_n_multiplicity_allN_epCDn_Step2);
    TH1D *h_n_multiplicity_goodN_epCDn_Step2 = new TH1D("n_multiplicity_goodN_epCDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epCDn_Step2);
    TH1D *h_n_multiplicity_badN_epCDn_Step2 = new TH1D("n_multiplicity_badN_epCDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epCDn_Step2);

    TH1D *h_n_multiplicity_allN_epFDn_Step2 = new TH1D("n_multiplicity_allN_epFDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epFDn_Step2);
    TH1D *h_n_multiplicity_goodN_epFDn_Step2 = new TH1D("n_multiplicity_goodN_epFDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epFDn_Step2);
    TH1D *h_n_multiplicity_badN_epFDn_Step2 = new TH1D("n_multiplicity_badN_epFDn_Step2", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epFDn_Step2);

    // Step Three (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step3 = new TH1D("n_multiplicity_allN_epCDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epCDn_Step3);
    TH1D *h_n_multiplicity_goodN_epCDn_Step3 = new TH1D("n_multiplicity_goodN_epCDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epCDn_Step3);
    TH1D *h_n_multiplicity_badN_epCDn_Step3 = new TH1D("n_multiplicity_badN_epCDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epCDn_Step3);

    TH1D *h_n_multiplicity_allN_epFDn_Step3 = new TH1D("n_multiplicity_allN_epFDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epFDn_Step3);
    TH1D *h_n_multiplicity_goodN_epFDn_Step3 = new TH1D("n_multiplicity_goodN_epFDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epFDn_Step3);
    TH1D *h_n_multiplicity_badN_epFDn_Step3 = new TH1D("n_multiplicity_badN_epFDn_Step3", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epFDn_Step3);

    // Step Four (After applying Phi Diff CND hit cut) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step4 = new TH1D("n_multiplicity_allN_epCDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epCDn_Step4);
    TH1D *h_n_multiplicity_goodN_epCDn_Step4 = new TH1D("n_multiplicity_goodN_epCDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epCDn_Step4);
    TH1D *h_n_multiplicity_badN_epCDn_Step4 = new TH1D("n_multiplicity_badN_epCDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epCDn_Step4);

    TH1D *h_n_multiplicity_allN_epFDn_Step4 = new TH1D("n_multiplicity_allN_epFDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epFDn_Step4);
    TH1D *h_n_multiplicity_goodN_epFDn_Step4 = new TH1D("n_multiplicity_goodN_epFDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epFDn_Step4);
    TH1D *h_n_multiplicity_badN_epFDn_Step4 = new TH1D("n_multiplicity_badN_epFDn_Step4", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epFDn_Step4);

    // Step Five (After event selection cuts) (Andrew)
    // ======================================================================================================================================================================

    /* Neutron histograms (from Erin) */
    TH1D *h_n_multiplicity_allN_epCDn_Step5 = new TH1D("n_multiplicity_allN_epCDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epCDn_Step5);
    TH1D *h_n_multiplicity_goodN_epCDn_Step5 = new TH1D("n_multiplicity_goodN_epCDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epCDn_Step5);
    TH1D *h_n_multiplicity_badN_epCDn_Step5 = new TH1D("n_multiplicity_badN_epCDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epCDn_Step5);

    TH1D *h_n_multiplicity_allN_epFDn_Step5 = new TH1D("n_multiplicity_allN_epFDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_allN_epFDn_Step5);
    TH1D *h_n_multiplicity_goodN_epFDn_Step5 = new TH1D("n_multiplicity_goodN_epFDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_goodN_epFDn_Step5);
    TH1D *h_n_multiplicity_badN_epFDn_Step5 = new TH1D("n_multiplicity_badN_epFDn_Step5", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_n_multiplicity_badN_epFDn_Step5);

    for (int i = 0; i < HistoList.size(); i++)
    {
        if (HistoList[i]->InheritsFrom("TH1D"))
        {
            HistoList[i]->Sumw2();
        }

        HistoList[i]->GetXaxis()->CenterTitle();
        HistoList[i]->GetYaxis()->CenterTitle();
    }

#pragma endregion /* Andrew's histograms - end */

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

        /*
                if ((counter_A % 100000) == 0)
                {
                    cerr << "\033[33m.\033[0m";
                }
        */

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

        TVector3 P_b(0, 0, Ebeam);

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

        TVector3 P_e(0., 0., 0.);

        double P_e_x = Electrons[0]->par()->getPx();
        double P_e_y = Electrons[0]->par()->getPy();
        double P_e_z = Electrons[0]->par()->getPz();

        P_e.SetXYZ(P_e_x, P_e_y, P_e_z);

        double Vz_e = Electrons[0]->par()->getVz();

        TVector3 P_q = P_b - P_e;            // 3-momentum transfer
        double nu = Ebeam - P_e.Mag();       // Energy transfer
        double QSq = P_q.Mag2() - (nu * nu); // 4-momentum transfer squared
        double xB = QSq / (2 * mN * nu);     // x Bjorken

        // Electrons (from Andrew)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        double EoP_e = (Electrons[0]->cal(PCAL)->getEnergy() + Electrons[0]->cal(ECIN)->getEnergy() + Electrons[0]->cal(ECOUT)->getEnergy()) / P_e.Mag();
        int nphe = Electrons[0]->che(HTCC)->getNphe();

        int e_sector = Electrons[0]->getSector();

        double theta_q = P_q.Theta() * 180 / M_PI;
        double WSq = (mN * mN) - QSq + (2 * nu * mN); // Hadronic mass
        double theta_e = P_e.Theta() * 180 / M_PI;

#pragma endregion /* Electrons - end */

#pragma region /* Protons - start */

        // Protons (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        h_psize->Fill(Protons.size());

        int counter_pCD_multiplicity_BPID = 0, counter_pFD_multiplicity_BPID = 0;
        int counter_pCD_multiplicity_APID = 0, counter_pFD_multiplicity_APID = 0;

        int p_index = -1;

        TVector3 P_p(0., 0., 0.);

        // Technically not optimized - this doesn't address what happens if there are two protons passing cuts
        // TODO: recheck this!
        for (int i = 0; i < Protons.size(); i++)
        {
            // define quantities
            P_p.SetMagThetaPhi(Protons[i]->getP(), Protons[i]->getTheta(), Protons[i]->getPhi());
            double dbeta = Protons[i]->par()->getBeta() - P_p.Mag() / sqrt(P_p.Mag2() + mP * mP);
            double p_theta = P_p.Theta() * 180. / M_PI;
            double Vz_p = Protons[i]->par()->getVz();
            double chipid = Protons[i]->par()->getChi2Pid();

            // fill histos
            h_pangles->Fill(P_p.Phi() * 180. / M_PI, p_theta);

            if (Protons[i]->getRegion() == FD)
            {
                ++counter_pFD_multiplicity_BPID;

                h_theta_p_VS_phi_p_BPID_epFD->Fill(P_p.Phi() * 180. / M_PI, p_theta, weight);
                h_P_p_BPID_epFD->Fill(P_p.Mag(), weight);
                h_dbeta_p_BPID_epFD->Fill(P_p.Mag(), dbeta, weight);
                h_dVz_p_BPID_epFD->Fill(Vz_p - Vz_e, weight);
                h_Chi2pid_p_BPID_epFD->Fill(chipid, weight);

                h_vzp_fd->Fill(Vz_p - Vz_e);

                if (fabs(Vz_p - Vz_e) > 5)
                {
                    continue;
                }

                h_chipid_fd->Fill(chipid);
                h_dbeta_p_fd->Fill(P_p.Mag(), dbeta);

                if (P_p.Mag() < 0.5 || P_p.Mag() > 3.0)
                {
                    continue;
                }

                if (fabs(dbeta) > 0.03)
                {
                    continue;
                }

                ++counter_pFD_multiplicity_APID;

                h_theta_p_VS_phi_p_APID_epFD->Fill(P_p.Phi() * 180. / M_PI, p_theta);
                h_P_p_APID_epFD->Fill(P_p.Mag(), weight);
                h_dbeta_p_APID_epFD->Fill(P_p.Mag(), dbeta, weight);
                h_dVz_p_APID_epFD->Fill(Vz_p - Vz_e, weight);
                h_Chi2pid_p_APID_epFD->Fill(chipid, weight);
            }
            else if (Protons[i]->getRegion() == CD)
            {
                ++counter_pCD_multiplicity_BPID;

                h_theta_p_VS_phi_p_BPID_epCD->Fill(P_p.Phi() * 180. / M_PI, p_theta);
                h_P_p_BPID_epCD->Fill(P_p.Mag(), weight);
                h_dbeta_p_BPID_epCD->Fill(P_p.Mag(), dbeta, weight);
                h_dVz_p_BPID_epCD->Fill(Vz_p - Vz_e, weight);
                h_Chi2pid_p_BPID_epCD->Fill(chipid, weight);

                h_vzp_cd->Fill(Vz_p - Vz_e);

                if (fabs(Vz_p - Vz_e) > 4)
                {
                    continue;
                }

                h_chipid_cd->Fill(chipid);
                h_dbeta_p_cd->Fill(P_p.Mag(), dbeta);

                if (P_p.Mag() < 0.3 || P_p.Mag() > 1.5)
                {
                    continue;
                }

                if (fabs(dbeta) > 0.05)
                {
                    continue;
                }

                ++counter_pCD_multiplicity_APID;

                h_theta_p_VS_phi_p_APID_epCD->Fill(P_p.Phi() * 180. / M_PI, p_theta);
                h_P_p_APID_epCD->Fill(P_p.Mag(), weight);
                h_dbeta_p_APID_epCD->Fill(P_p.Mag(), dbeta, weight);
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

        P_p.SetMagThetaPhi(Protons[p_index]->getP(), Protons[p_index]->getTheta(), Protons[p_index]->getPhi());

        // Determin where is the proton. Moved from angle cuts to getRegion() by the advice of Andrew.
        bool pInFD = (Protons[p_index]->getRegion() == FD); // My addition
        bool pInCD = (Protons[p_index]->getRegion() == CD); // My addition

#pragma endregion /* Protons - end */

#pragma region /* Missing momentum - start */

        // Missing momentum (from Erin)
        // -------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Missing momentum, energy, mass
        TVector3 P_miss = P_q - P_p; // TODO: checkout difference from Andrew - he uses leading SRC proton here!

        momentum = P_miss.Mag();

        double E_p = sqrt(mN * mN + P_p.Mag2());
        double E_miss = Ebeam + mD - P_e.Mag() - E_p;
        double M_miss = sqrt((E_miss * E_miss) - P_miss.Mag2());

#pragma endregion /* Missing momentum - end */

        // ==================================================================================================================================================================
        // Erin's features
        // ==================================================================================================================================================================

#pragma region /* Erin's features - start */

        if (Run_Erins_features)
        {
            // initialize features
            energy = 0;
            cnd_energy = 0;
            ctof_energy = 0;
            angle_diff = 180;
            layermult = 0;
            size = 0;
            cnd_hits = 0;
            ctof_hits = 0;
            is_CTOF = 0;
            is_CND1 = 0;
            is_CND2 = 0;
            is_CND3 = 0;

            /*
             // Particle PID
            // ===================================================================================================================================================================

            clasAna->Run(c12);

            auto elec = clasAna->getByPid(11);
            auto prot = clasAna->getByPid(2212);
            auto neut = clasAna->getByPid(2112);

            auto AllParticles = c12->getDetParticles();

            if (elec.size() != 1) // One electron in event
            {
                continue;
            }

            if (prot.size() != 1) // One proton in event
            {
                continue;
            }

            if (neut.size() < 1) // At least one neutron in event
            {
                continue;
            }

            event = c12->runconfig()->getEvent() << '\n';

            // reject particles with the wrong PID
            bool trash = 0;

            for (int i = 0; i < AllParticles.size(); i++) // TODO: ask Larry if I should use these from the event selection for Andrew's work
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

            numevent = numevent + 1;
            */

#pragma region /* Electrons (Erin) - start */

            /*
            //////////////////////////
            /////    ELECTRONS   /////
            //////////////////////////
            TVector3 pe(0., 0., 0.);

            double pe_x = elec[0]->par()->getPx();
            double pe_y = elec[0]->par()->getPy();
            double pe_z = elec[0]->par()->getPz();

            pe.SetXYZ(pe_x, pe_y, pe_z);

            double vze = elec[0]->par()->getVz();

            TVector3 P_b(0, 0, Ebeam);
            TVector3 pq = P_b - pe; // 3-momentum transfer

            double nu = Ebeam - pe.Mag();       // Energy transfer
            double QSq = pq.Mag2() - (nu * nu); // 4-momentum transfer squared
            double xB = QSq / (2 * mN * nu);    // x Bjorken
            */

#pragma endregion /* Electrons (Erin) - end */

#pragma region /* Protons - start */

            /*
            //////////////////////////
            /////     PROTONS    /////
            //////////////////////////
            h_psize->Fill(prot.size());
            int p_index = -1;
            TVector3 pp(0., 0., 0.);

            // technically not optimized - this doesn't address what happens if there are two protons passing cuts
            // TODO: recheck this!
            for (int i = 0; i < prot.size(); i++)
            {
                // define quantities
                pp.SetMagThetaPhi(prot[i]->getP(), prot[i]->getTheta(), prot[i]->getPhi());
                double dbeta = prot[i]->par()->getBeta() - pp.Mag() / sqrt(pp.Mag2() + mP * mP);
                double p_theta = pp.Theta() * 180. / M_PI;
                double vzp = prot[i]->par()->getVz();
                double chipid = prot[i]->par()->getChi2Pid();

                // fill histos
                h_pangles->Fill(pp.Phi() * 180. / M_PI, p_theta);

                if (prot[i]->getRegion() == FD)
                {
                    h_vzp_fd->Fill(vzp - vze);

                    if (fabs(vzp - vze) > 5)
                    // if (abs(vzp - vze) > 5) // Erin's original
                    {
                        continue;
                    }

                    h_chipid_fd->Fill(chipid);
                    h_dbeta_p_fd->Fill(pp.Mag(), dbeta);

                    if (pp.Mag() < 0.5)
                    {
                        continue;
                    }

                    if (pp.Mag() > 3.0)
                    {
                        continue;
                    }

                    if (fabs(dbeta) > 0.03)
                    // if (abs(dbeta) > 0.03) // Erin's original
                    {
                        continue;
                    }
                }
                else if (prot[i]->getRegion() == CD)
                {
                    h_vzp_cd->Fill(vzp - vze);

                    if (fabs(vzp - vze) > 4)
                    // if (abs(vzp - vze) > 4) // Erin's original
                    {
                        continue;
                    }

                    h_chipid_cd->Fill(chipid);
                    // if (abs(chipid)>4) {continue;}
                    h_dbeta_p_cd->Fill(pp.Mag(), dbeta);

                    if (pp.Mag() < 0.3)
                    {
                        continue;
                    }

                    if (pp.Mag() > 1.5)
                    {
                        continue;
                    }

                    if (fabs(dbeta) > 0.05)
                    // if (abs(dbeta) > 0.05) // Erin's original
                    {
                        continue;
                    }
                }

                p_index = i;
            }

            if (p_index < 0)
            {
                continue;
            }

            pp.SetMagThetaPhi(prot[p_index]->getP(), prot[p_index]->getTheta(), prot[p_index]->getPhi());

            if (pp.Theta() * 180. / M_PI < 40 || pp.Theta() * 180. / M_PI > 140) // p goes to CD
            {
                continue;
            }
            // if (pp.Theta()*180./M_PI>40) {continue;}  // p goes to FD
            */

#pragma endregion /* Protons - end */

#pragma region /* Missing momentum - start */

            /*
            //////////////////////////
            //  MISSING MOMENTUM    //
            //////////////////////////

            // missing momentum, energy, mass
            TVector3 pmiss = pq - pp; // TODO: checkout difference from Andrew - he uses leading SRC proton here!

            momentum = pmiss.Mag();

            double Ep = sqrt(mN * mN + pp.Mag2());
            double Emiss = Ebeam + mD - pe.Mag() - Ep;
            double mmiss = sqrt((Emiss * Emiss) - pmiss.Mag2());
            */

#pragma endregion /* Missing momentum - end */

#pragma region /* Neutrons - start */

            //////////////////////////
            ////     NEUTRONS    /////
            //////////////////////////

            // LOOP OVER NEUTRONS
            h_n_multiplicity->Fill(Neutrons.size());

            for (int i = 0; i < Neutrons.size(); i++)
            {
                // GET NEUTRON INFORMATION

                // get neutron momentums
                double P_n_x = Neutrons[i]->par()->getPx();
                double P_n_y = Neutrons[i]->par()->getPy();
                double P_n_z = Neutrons[i]->par()->getPz();

                TVector3 P_n;
                P_n.SetXYZ(P_n_x, P_n_y, P_n_z);

                double dpp = (P_miss.Mag() - P_n.Mag()) / P_miss.Mag();

                // figure out what layer the hit is in
                is_CND1 = (Neutrons[i]->sci(CND1)->getLayer() == 1);
                is_CND2 = (Neutrons[i]->sci(CND2)->getLayer() == 2);
                is_CND3 = (Neutrons[i]->sci(CND3)->getLayer() == 3);
                is_CTOF = Neutrons[i]->sci(CTOF)->getDetector() == 4;

                // put REC::Scintillator information
                double time; // Neutron TOF

                int status = 0;

                double beta = Neutrons[i]->par()->getBeta();

                if (is_CND1)
                {
                    time = Neutrons[i]->sci(CND1)->getTime() - starttime;
                    status = status + Neutrons[i]->sci(CND1)->getStatus();
                }

                if (is_CND3)
                {
                    time = Neutrons[i]->sci(CND3)->getTime() - starttime;
                    status = status + Neutrons[i]->sci(CND3)->getStatus();
                }

                if (is_CND2)
                {
                    time = Neutrons[i]->sci(CND2)->getTime() - starttime;
                    status = status + Neutrons[i]->sci(CND2)->getStatus();
                }

                // PROBLEM: this gives preference to 2nd-layer hits
                // TODO: recheck this!
                if (is_CTOF)
                {
                    time = Neutrons[i]->sci(CTOF)->getTime() - starttime;
                }

                double cos0 = P_miss.Dot(P_n) / (P_miss.Mag() * P_n.Mag());

                if (status != 0) // Cutting out neutrons suspected to have double-hit CND hits
                {
                    continue;
                }

                // GET ML FEATURES FOR THIS NEUTRON
                Struct ninfo = getFeatures(Neutrons, AllParticles, i);
                cnd_hits = ninfo.cnd_hits;
                ctof_hits = ninfo.ctof_hits;
                cnd_energy = ninfo.cnd_energy;
                ctof_energy = ninfo.ctof_energy;
                layermult = ninfo.layermult;
                energy = ninfo.energy;
                size = ninfo.size;
                angle_diff = ninfo.angle_diff;

                if (cnd_energy > 1000) // TODO: why this cut?
                {
                    continue;
                }

                // ESSENTIAL NEUTRONS CUTS
                h_tof->Fill(time);

                if (P_n_x == 0 || P_n_y == 0 || P_n_z == 0) // Negative TOF cut
                {
                    continue;
                }

                if (time > 10) // TODO: why this cut?
                {
                    continue;
                }

                h_pvsp->Fill(P_miss.Mag(), P_n.Mag());

                // select neutrons in momentum and angle accepted by CND
                double theta_n = P_n.Theta() * 180. / M_PI;

                if (P_n.Mag() < 0.25 || P_n.Mag() > 1)
                {
                    continue;
                }

                if (theta_n < 45 || theta_n > 140)
                {
                    continue;
                }

                h_mmiss_pn->Fill(P_n.Mag(), M_miss);
                h_mmiss_pmiss->Fill(P_miss.Mag(), M_miss);
                h_mmiss_xb->Fill(xB, M_miss);

                h_mmiss->Fill(M_miss);

                if (M_miss > 1.) // Missing mass cut
                {
                    continue;
                }

                h_pmiss_thetamiss->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag());
                h_thetapn_pp->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI);

                if (P_miss.Mag() < 0.25 || P_miss.Mag() > 1.) // Missing momentum cut
                {
                    continue;
                }

                if (P_miss.Theta() * 180. / M_PI < 45 || P_miss.Theta() * 180. / M_PI > 140) // Missing momentum theta cut
                {
                    continue;
                }

                h_dpp_edep->Fill(energy, dpp);

                if (energy < 5)
                {
                    continue;
                }

                // FILL HISTOS FOR NEUTRON CANDIDATES
                h_nangles->Fill(P_n.Phi() * 180. / M_PI, theta_n);
                h_energy->Fill(energy);
                h_Edep_beta->Fill(Neutrons[i]->getBeta(), energy);

                h_cos0->Fill(P_miss.Dot(P_n) / (P_miss.Mag() * P_n.Mag()));
                h_pxminuspx->Fill(P_n_x - P_miss.X());
                h_pyminuspy->Fill(P_n_y - P_miss.Y());
                h_pzminuspz->Fill(P_n_z - P_miss.Z());
                h_pxminuspx->Fill(P_n_x - P_miss.X(), weight);
                h_pyminuspy->Fill(P_n_y - P_miss.Y(), weight);
                h_pzminuspz->Fill(P_n_z - P_miss.Z(), weight);
                h_pminusp->Fill(P_n.Mag() - P_miss.Mag());

                h_dpp->Fill(P_miss.Mag(), (P_miss.Mag() - P_n.Mag()) / P_miss.Mag());
                h_theta_beta->Fill(beta, theta_n);
                h_p_theta->Fill(theta_n, P_n.Mag());
                h_p_all->Fill(P_miss.Mag());
                h_anglediff->Fill(angle_diff);

                h_compare->Fill((P_miss.Mag() - P_n.Mag()) / P_miss.Mag(), P_n.Angle(P_miss) * 180. / M_PI);

                if ((fabs(P_miss.Mag() - P_n.Mag()) / P_miss.Mag()) > 0.2) // Relative momentum difference cut
                // if ((abs(P_miss.Mag() - pn.Mag()) / P_miss.Mag()) > 0.2) // Erin's original
                {
                    continue;
                }

                h_thetapn_dpp->Fill((P_miss.Mag() - P_n.Mag()) / P_miss.Mag(), P_n.Angle(P_p) * 180. / M_PI);
                h_thetapn_dpp1->Fill((P_miss.Mag() - P_n.Mag()) / P_miss.Mag(), P_n.Angle(P_p) * 180. / M_PI);

                if (P_n.Angle(P_miss) * 180. / M_PI > 20) // pn close to P_miss cut
                {
                    continue;
                }

                // ML features
                h_energy_1->Fill(energy);
                h_layermult_1->Fill(layermult);
                h_size_1->Fill(size);
                h_cnd_hits_1->Fill(cnd_hits);
                h_cnd_energy_1->Fill(cnd_energy);
                h_ctof_energy_1->Fill(ctof_energy);
                h_ctof_hits_1->Fill(ctof_hits);
                h_anglediff_1->Fill(angle_diff);

                // physics cuts - not being used
                // if (angle_diff<30) {continue;}
                // if (cnd_hits>2) {continue;}
                // if (size>1) {continue;}
                // if (ctof_hits>2) {continue;}
                // if (ctof_hits>0) {continue;}

                //////////////////////////
                /////     SORT       /////
                ////   GOOD / BAD    /////
                ////    NEUTRONS     /////
                //////////////////////////

                // Determine whether to write to "good (signal) neutron" or "bad (background) neutron" file

                bool good_N = (P_n.Angle(P_miss) * 180. / M_PI < 20) &&
                              (fabs((P_miss.Mag() - P_n.Mag()) / P_miss.Mag()) < 0.2) &&
                              //   abs((P_miss.Mag() - pn.Mag()) / P_miss.Mag()) < 0.2 && // Erin's original
                              (cnd_energy < 1000) &&
                              (P_p.Angle(P_n) * 180. / M_PI > 60) &&
                              (P_miss.Mag() > 0.25 && P_miss.Mag() < 1.) &&
                              (P_miss.Theta() * 180. / M_PI > 45 && P_miss.Theta() * 180. / M_PI < 140);
                // bool good_N = pn.Angle(P_miss) * 180. / M_PI < 20 &&
                //               fabs((P_miss.Mag() - pn.Mag()) / P_miss.Mag()) < 0.2 &&
                //               //   abs((P_miss.Mag() - pn.Mag()) / P_miss.Mag()) < 0.2 && // Erin's original
                //               cnd_energy < 1000 &&
                //               pp.Angle(pn) * 180. / M_PI > 60 &&
                //               (P_miss.Mag() > 0.25 &&
                //                P_miss.Mag() < 1.) &&
                //               (P_miss.Theta() * 180. / M_PI > 45 &&
                //                P_miss.Theta() * 180. / M_PI < 140);

                bool bad_N = ((P_n.Angle(P_miss) * 180. / M_PI > 50) ||
                              (fabs((P_miss.Mag() - P_n.Mag()) / P_miss.Mag()) > 0.6)) &&
                             //   (abs((P_miss.Mag() - pn.Mag()) / P_miss.Mag()) > 0.6)) && // Erin's original
                             cnd_energy < 1000; // && (pp.Angle(pn)*180./M_PI<60);

                bool keep_this_one = keep_good ? good_N : bad_N;

                if (keep_this_one)
                {
                    // all neutrons - print features
                    outtxt << P_miss.Mag() << ' ';
                    cout << P_miss.Mag() << ' ';
                    outtxt << energy << ' ';
                    cout << energy << ' ';
                    outtxt << layermult << ' ';
                    cout << layermult << ' ';
                    outtxt << size << ' ';
                    cout << size << ' ';
                    outtxt << cnd_hits << ' ';
                    cout << cnd_hits << ' ';
                    outtxt << cnd_energy << ' ';
                    cout << cnd_energy << ' ';
                    outtxt << ctof_energy << ' ';
                    cout << ctof_energy << ' ';
                    outtxt << ctof_hits << ' ';
                    cout << ctof_hits << ' ';
                    outtxt << angle_diff << ' ';
                    cout << angle_diff << ' ';
                    outtxt << '\n';

                    // FILL HISTOS FOR SIGNAL/BACKGROUND EVENTS
                    h_nangles2->Fill(P_n.Phi() * 180. / M_PI, theta_n);
                    h_cos02->Fill(P_miss.Dot(P_n) / (P_miss.Mag() * P_n.Mag()));
                    h_pxminuspx2->Fill(P_n_x - P_miss.X());
                    h_pyminuspy2->Fill(P_n_y - P_miss.Y());
                    h_pzminuspz2->Fill(P_n_z - P_miss.Z());
                    h_pminusp2->Fill(P_n.Mag() - P_miss.Mag());
                    h_pvsp2->Fill(P_miss.Mag(), P_n.Mag());
                    h_dpp2->Fill(P_miss.Mag(), (P_miss.Mag() - P_n.Mag()) / P_miss.Mag());
                    h_mmiss2->Fill(M_miss);
                    h_mmiss_pn2->Fill(P_n.Mag(), M_miss);
                    h_energy2->Fill(energy);
                    h_theta_beta2->Fill(beta, theta_n);
                    h_p_theta2->Fill(theta_n, P_n.Mag());
                    h_pmiss_thetamiss2->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag());
                    h_thetapn_pp2->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI);
                    h_tof2->Fill(time);
                    h_compare2->Fill((P_miss.Mag() - P_n.Mag()) / P_miss.Mag(), P_n.Angle(P_miss) * 180. / M_PI);
                    h_Edep_beta2->Fill(Neutrons[i]->getBeta(), energy);
                    h_p_cut->Fill(P_miss.Mag());
                    h_anglediff2->Fill(angle_diff);
                    h_thetapn_dpp2->Fill((P_miss.Mag() - P_n.Mag()) / P_miss.Mag(), P_n.Angle(P_p) * 180. / M_PI);

                    h_ptheta_pred->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag());
                    h_ptheta->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag());

                    // ML features
                    h_energy_2->Fill(energy);
                    h_layermult_2->Fill(layermult);
                    h_size_2->Fill(size);
                    h_cnd_hits_2->Fill(cnd_hits);
                    h_cnd_energy_2->Fill(cnd_energy);
                    h_ctof_energy_2->Fill(ctof_energy);
                    h_ctof_hits_2->Fill(ctof_hits);
                    h_anglediff_2->Fill(angle_diff);

                    // write events to tree
                    ntree->Fill();

                } // closes condition for good/bad neutron

            } // closes neutron loop

            // chain.WriteEvent();
            counter++;

#pragma endregion /* Neutrons - end */
        }

#pragma endregion /* Erin's features - end */

        // ==================================================================================================================================================================
        // Andrew's manual work
        // ==================================================================================================================================================================

#pragma region /* Andrew's manual work */

        if (Run_Andrews_work)
        {
#pragma region /* Missing momentum cuts (Andrew) - start */

            if (pInCD)
            {
                h_P_miss_BmissC_epCD->Fill(P_miss.Mag(), weight);
                h_theta_miss_BmissC_epCD->Fill(P_miss.Theta() * 180 / M_PI, weight);
                h_P_miss_VS_theta_miss_BmissC_epCD->Fill(P_miss.Theta() * 180, P_miss.Mag(), weight);
                h_E_p_BmissC_epCD->Fill(E_p, weight);
                h_E_miss_BmissC_epCD->Fill(E_miss, weight);
                h_M_miss_BmissC_epCD->Fill(M_miss, weight);
                h_xB_BmissC_epCD->Fill(xB, weight);
                h_xB_VS_M_miss_BmissC_epCD->Fill(xB, M_miss, weight);
            }
            else if (pInFD)
            {
                h_P_miss_BmissC_epFD->Fill(P_miss.Mag(), weight);
                h_theta_miss_BmissC_epFD->Fill(P_miss.Theta() * 180 / M_PI, weight);
                h_P_miss_VS_theta_miss_BmissC_epFD->Fill(P_miss.Theta() * 180, P_miss.Mag(), weight);
                h_E_p_BmissC_epFD->Fill(E_p, weight);
                h_E_miss_BmissC_epFD->Fill(E_miss, weight);
                h_M_miss_BmissC_epFD->Fill(M_miss, weight);
                h_xB_BmissC_epFD->Fill(xB, weight);
                h_xB_VS_M_miss_BmissC_epFD->Fill(xB, M_miss, weight);
            }

            if (P_miss.Theta() * 180 / M_PI < 40 || P_miss.Theta() * 180 / M_PI > 135)
            {
                continue;
            }

            if (P_miss.Mag() < 0.2 || P_miss.Mag() > 1.5)
            {
                continue;
            }

            if (M_miss < 0.7 || M_miss > 1.2)
            {
                continue;
            }

            if (pInCD)
            {
                h_P_miss_AmissC_epCD->Fill(P_miss.Mag(), weight);
                h_theta_miss_AmissC_epCD->Fill(P_miss.Theta() * 180 / M_PI, weight);
                h_P_miss_VS_theta_miss_AmissC_epCD->Fill(P_miss.Theta() * 180, P_miss.Mag(), weight);
                h_E_p_AmissC_epCD->Fill(E_p, weight);
                h_E_miss_AmissC_epCD->Fill(E_miss, weight);
                h_M_miss_AmissC_epCD->Fill(M_miss, weight);
                h_xB_AmissC_epCD->Fill(xB, weight);
                h_xB_VS_M_miss_AmissC_epCD->Fill(xB, M_miss, weight);
            }
            else if (pInFD)
            {
                h_P_miss_AmissC_epFD->Fill(P_miss.Mag(), weight);
                h_theta_miss_AmissC_epFD->Fill(P_miss.Theta() * 180 / M_PI, weight);
                h_P_miss_VS_theta_miss_AmissC_epFD->Fill(P_miss.Theta() * 180, P_miss.Mag(), weight);
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

                double nvtx_x = AllParticles[itr1]->par()->getVx();
                double nvtx_y = AllParticles[itr1]->par()->getVy();
                double nvtx_z = AllParticles[itr1]->par()->getVz();
                TVector3 v_nvtx(nvtx_x, nvtx_y, nvtx_z); // Neutron's vertex location

                TVector3 v_hit; // Neutron's hit location in CND
                v_hit.SetXYZ(AllParticles[itr1]->sci(detlayer)->getX(), AllParticles[itr1]->sci(detlayer)->getY(), AllParticles[itr1]->sci(detlayer)->getZ());

                TVector3 v_path = v_hit - v_nvtx; // Direct calculation of neutron's path (in vector form)
                TVector3 P_n;
                P_n.SetMagThetaPhi(mom, v_path.Theta(), v_path.Phi()); // Direct calculation of neutron momentum?
                                                                       // TODO: check with Andrew why he calculated this explicitly

                // Why "v_path.Mag() / 100"? unit conversion.
                // TODO: check if this unit conversion is needed!
                double path = v_path.Mag() / 100;
                // double path = v_path.Mag();
                double theta_n_miss = P_n.Angle(P_miss) * 180 / M_PI; // Opening angle between calculated neutron's momentum and predicted neutron momentum (= missing momentum)
                TVector3 dp = P_miss - P_n;
                // double dpp = dp.Mag() / P_miss.Mag();
                double dpp = (P_miss.Mag() - P_n.Mag()) / P_miss.Mag();
                int nSector = AllParticles[itr1]->sci(detlayer)->getSector(); // Number of CND sector with a neutron hit in the layer detlayer

                // Check to see if there is a good neutron
                bool isGN = false;

                // Why this cut? reco code bug. Neutrons in this angle range are in the BAND and appear in the CND.
                // This bug is probobly fixed, yet the cut is still applied to mak sure.
                if (P_n.Theta() * 180. / M_PI > 160)
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
                    h_theta_n_epCDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_epCDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_epCDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_epCDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                    h_P_n_epCDn->Fill(P_n.Mag(), weight);
                    h_P_n_VS_theta_n_epCDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                    h_P_miss_epCDn->Fill(P_miss.Mag(), weight);
                    h_P_miss_VS_theta_miss_epCDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);

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
                    h_M_miss_VS_P_n_epCDn->Fill(P_n.Mag(), M_miss, weight);
                    h_M_miss_VS_theta_n_epCDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                    h_M_miss_VS_phi_n_epCDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                    h_M_miss_VS_P_miss_epCDn->Fill(P_miss.Mag(), M_miss, weight);
                    h_M_miss_VS_theta_miss_epCDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                    h_M_miss_VS_phi_miss_epCDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                    h_P_n_minus_P_miss_epCDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                    h_P_n_x_minus_P_miss_x_epCDn->Fill(P_n.X() - P_miss.X(), weight);
                    h_P_n_y_minus_P_miss_y_epCDn->Fill(P_n.Y() - P_miss.Y(), weight);
                    h_P_n_z_minus_P_miss_z_epCDn->Fill(P_n.Z() - P_miss.Z(), weight);

                    h_P_n_VS_P_miss_epCDn->Fill(P_miss.Mag(), P_n.Mag(), weight);
                    h_P_n_x_VS_P_miss_x_epCDn->Fill(P_n.X() - P_miss.X(), weight);
                    h_P_n_y_VS_P_miss_y_epCDn->Fill(P_n.Y() - P_miss.Y(), weight);
                    h_P_n_z_VS_P_miss_z_epCDn->Fill(P_n.Z() - P_miss.Z(), weight);

                    h_theta_n_p_epCDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                    h_theta_n_p_VS_P_p_epCDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                    h_xB_epCDn->Fill(xB, weight);

                    h_Edep_CND_epCDn->Fill(Edep_CND, weight);
                    h_P_n_VS_Edep_CND_epCDn->Fill(Edep_CND, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CND_epCDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CND_epCDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CND_epCDn->Fill(Edep_CND, dpp, weight);
                    h_beta_n_VS_Edep_CND_epCDn->Fill(Edep_CND, beta, weight);
                    h_E_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, E_miss, weight);
                    h_M_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, M_miss, weight);
                    h_path_VS_Edep_CND_epCDn->Fill(Edep_CND, path, weight);
                    h_theta_n_miss_VS_Edep_CND_epCDn->Fill(Edep_CND, theta_n_miss, weight);
                    h_ToF_VS_Edep_CND_epCDn->Fill(Edep_CND, ToF, weight);
                    h_nSector_VS_Edep_CND_epCDn->Fill(Edep_CND, nSector, weight);

                    h_Edep_CTOF_epCDn->Fill(Edep_CTOF, weight);
                    h_P_n_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, dpp, weight);
                    h_beta_n_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, beta, weight);
                    h_E_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, E_miss, weight);
                    h_M_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, M_miss, weight);
                    h_path_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, path, weight);
                    h_theta_n_miss_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
                    h_ToF_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, ToF, weight);
                    h_nSector_VS_Edep_CTOF_epCDn->Fill(Edep_CTOF, nSector, weight);

                    h_Edep_single_epCDn->Fill(Edep_single, weight);
                    h_P_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_single_epCDn->Fill(Edep_single, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_single_epCDn->Fill(Edep_single, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_single_epCDn->Fill(Edep_single, dpp, weight);
                    h_beta_n_VS_Edep_single_epCDn->Fill(Edep_single, beta, weight);
                    h_E_miss_VS_Edep_single_epCDn->Fill(Edep_single, E_miss, weight);
                    h_M_miss_VS_Edep_single_epCDn->Fill(Edep_single, M_miss, weight);
                    h_path_VS_Edep_single_epCDn->Fill(Edep_single, path, weight);
                    h_theta_n_miss_VS_Edep_single_epCDn->Fill(Edep_single, theta_n_miss, weight);
                    h_ToF_VS_Edep_single_epCDn->Fill(Edep_single, ToF, weight);
                    h_nSector_VS_Edep_single_epCDn->Fill(Edep_single, nSector, weight);

                    h_Edep_CND1_epCDn->Fill(Edep_CND1, weight);
                    h_P_n_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CND1_epCDn->Fill(Edep_CND1, dpp, weight);
                    h_beta_n_VS_Edep_CND1_epCDn->Fill(Edep_CND1, beta, weight);
                    h_E_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, E_miss, weight);
                    h_M_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, M_miss, weight);
                    h_path_VS_Edep_CND1_epCDn->Fill(Edep_CND1, path, weight);
                    h_theta_n_miss_VS_Edep_CND1_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
                    h_ToF_VS_Edep_CND1_epCDn->Fill(Edep_CND1, ToF, weight);
                    h_nSector_VS_Edep_CND1_epCDn->Fill(Edep_CND1, nSector, weight);

                    h_Edep_CND2_epCDn->Fill(Edep_CND2, weight);
                    h_P_n_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CND2_epCDn->Fill(Edep_CND2, dpp, weight);
                    h_beta_n_VS_Edep_CND2_epCDn->Fill(Edep_CND2, beta, weight);
                    h_E_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, E_miss, weight);
                    h_M_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, M_miss, weight);
                    h_path_VS_Edep_CND2_epCDn->Fill(Edep_CND2, path, weight);
                    h_theta_n_miss_VS_Edep_CND2_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
                    h_ToF_VS_Edep_CND2_epCDn->Fill(Edep_CND2, ToF, weight);
                    h_nSector_VS_Edep_CND2_epCDn->Fill(Edep_CND2, nSector, weight);

                    h_Edep_CND3_epCDn->Fill(Edep_CND3, weight);
                    h_P_n_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CND3_epCDn->Fill(Edep_CND3, dpp, weight);
                    h_beta_n_VS_Edep_CND3_epCDn->Fill(Edep_CND3, beta, weight);
                    h_E_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, E_miss, weight);
                    h_M_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, M_miss, weight);
                    h_path_VS_Edep_CND3_epCDn->Fill(Edep_CND3, path, weight);
                    h_theta_n_miss_VS_Edep_CND3_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
                    h_ToF_VS_Edep_CND3_epCDn->Fill(Edep_CND3, ToF, weight);
                    h_nSector_VS_Edep_CND3_epCDn->Fill(Edep_CND3, nSector, weight);

                    h_ToF_epCDn->Fill(ToF, weight);
                    h_P_n_VS_ToF_epCDn->Fill(ToF, P_n.Mag(), weight);
                    h_theta_n_VS_ToF_epCDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_ToF_epCDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_ToF_epCDn->Fill(ToF, P_miss.Mag(), weight);
                    h_theta_miss_VS_ToF_epCDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_ToF_epCDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_ToF_epCDn->Fill(ToF, dpp, weight);
                    h_beta_n_VS_ToF_epCDn->Fill(ToF, beta, weight);
                    h_E_miss_VS_ToF_epCDn->Fill(ToF, E_miss, weight);
                    h_M_miss_VS_ToF_epCDn->Fill(ToF, M_miss, weight);
                    h_path_VS_ToF_epCDn->Fill(ToF, path, weight);
                    h_theta_n_miss_VS_ToF_epCDn->Fill(ToF, theta_n_miss, weight);
                    h_nSector_VS_ToF_epCDn->Fill(ToF, nSector, weight);
                }
                else if (pInFD)
                {
                    h_theta_n_epFDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_epFDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_epFDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_epFDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                    h_P_n_epFDn->Fill(P_n.Mag(), weight);
                    h_P_n_VS_theta_n_epFDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                    h_P_miss_epFDn->Fill(P_miss.Mag(), weight);
                    h_P_miss_VS_theta_miss_epFDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);

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
                    h_M_miss_VS_P_n_epFDn->Fill(P_n.Mag(), M_miss, weight);
                    h_M_miss_VS_theta_n_epFDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                    h_M_miss_VS_phi_n_epFDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                    h_M_miss_VS_P_miss_epFDn->Fill(P_miss.Mag(), M_miss, weight);
                    h_M_miss_VS_theta_miss_epFDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                    h_M_miss_VS_phi_miss_epFDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                    h_P_n_minus_P_miss_epFDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                    h_P_n_x_minus_P_miss_x_epFDn->Fill(P_n.X() - P_miss.X(), weight);
                    h_P_n_y_minus_P_miss_y_epFDn->Fill(P_n.Y() - P_miss.Y(), weight);
                    h_P_n_z_minus_P_miss_z_epFDn->Fill(P_n.Z() - P_miss.Z(), weight);

                    h_P_n_VS_P_miss_epFDn->Fill(P_miss.Mag(), P_n.Mag(), weight);
                    h_P_n_x_VS_P_miss_x_epFDn->Fill(P_n.X() - P_miss.X(), weight);
                    h_P_n_y_VS_P_miss_y_epFDn->Fill(P_n.Y() - P_miss.Y(), weight);
                    h_P_n_z_VS_P_miss_z_epFDn->Fill(P_n.Z() - P_miss.Z(), weight);

                    h_theta_n_p_epFDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                    h_theta_n_p_VS_P_p_epFDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                    h_xB_epFDn->Fill(xB, weight);

                    h_Edep_CND_epFDn->Fill(Edep_CND, weight);
                    h_P_n_VS_Edep_CND_epFDn->Fill(Edep_CND, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CND_epFDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CND_epFDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CND_epFDn->Fill(Edep_CND, dpp, weight);
                    h_beta_n_VS_Edep_CND_epFDn->Fill(Edep_CND, beta, weight);
                    h_E_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, E_miss, weight);
                    h_M_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, M_miss, weight);
                    h_path_VS_Edep_CND_epFDn->Fill(Edep_CND, path, weight);
                    h_theta_n_miss_VS_Edep_CND_epFDn->Fill(Edep_CND, theta_n_miss, weight);
                    h_ToF_VS_Edep_CND_epFDn->Fill(Edep_CND, ToF, weight);
                    h_nSector_VS_Edep_CND_epFDn->Fill(Edep_CND, nSector, weight);

                    h_Edep_CTOF_epFDn->Fill(Edep_CTOF, weight);
                    h_P_n_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, dpp, weight);
                    h_beta_n_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, beta, weight);
                    h_E_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, E_miss, weight);
                    h_M_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, M_miss, weight);
                    h_path_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, path, weight);
                    h_theta_n_miss_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
                    h_ToF_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, ToF, weight);
                    h_nSector_VS_Edep_CTOF_epFDn->Fill(Edep_CTOF, nSector, weight);

                    h_Edep_single_epFDn->Fill(Edep_single, weight);
                    h_P_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_single_epFDn->Fill(Edep_single, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_single_epFDn->Fill(Edep_single, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_single_epFDn->Fill(Edep_single, dpp, weight);
                    h_beta_n_VS_Edep_single_epFDn->Fill(Edep_single, beta, weight);
                    h_E_miss_VS_Edep_single_epFDn->Fill(Edep_single, E_miss, weight);
                    h_M_miss_VS_Edep_single_epFDn->Fill(Edep_single, M_miss, weight);
                    h_path_VS_Edep_single_epFDn->Fill(Edep_single, path, weight);
                    h_theta_n_miss_VS_Edep_single_epFDn->Fill(Edep_single, theta_n_miss, weight);
                    h_ToF_VS_Edep_single_epFDn->Fill(Edep_single, ToF, weight);
                    h_nSector_VS_Edep_single_epFDn->Fill(Edep_single, nSector, weight);

                    h_Edep_CND1_epFDn->Fill(Edep_CND1, weight);
                    h_P_n_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CND1_epFDn->Fill(Edep_CND1, dpp, weight);
                    h_beta_n_VS_Edep_CND1_epFDn->Fill(Edep_CND1, beta, weight);
                    h_E_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, E_miss, weight);
                    h_M_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, M_miss, weight);
                    h_path_VS_Edep_CND1_epFDn->Fill(Edep_CND1, path, weight);
                    h_theta_n_miss_VS_Edep_CND1_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
                    h_ToF_VS_Edep_CND1_epFDn->Fill(Edep_CND1, ToF, weight);
                    h_nSector_VS_Edep_CND1_epFDn->Fill(Edep_CND1, nSector, weight);

                    h_Edep_CND2_epFDn->Fill(Edep_CND2, weight);
                    h_P_n_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CND2_epFDn->Fill(Edep_CND2, dpp, weight);
                    h_beta_n_VS_Edep_CND2_epFDn->Fill(Edep_CND2, beta, weight);
                    h_E_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, E_miss, weight);
                    h_M_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, M_miss, weight);
                    h_path_VS_Edep_CND2_epFDn->Fill(Edep_CND2, path, weight);
                    h_theta_n_miss_VS_Edep_CND2_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
                    h_ToF_VS_Edep_CND2_epFDn->Fill(Edep_CND2, ToF, weight);
                    h_nSector_VS_Edep_CND2_epFDn->Fill(Edep_CND2, nSector, weight);

                    h_Edep_CND3_epFDn->Fill(Edep_CND3, weight);
                    h_P_n_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_n.Mag(), weight);
                    h_theta_n_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                    h_theta_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_Edep_CND3_epFDn->Fill(Edep_CND3, dpp, weight);
                    h_beta_n_VS_Edep_CND3_epFDn->Fill(Edep_CND3, beta, weight);
                    h_E_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, E_miss, weight);
                    h_M_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, M_miss, weight);
                    h_path_VS_Edep_CND3_epFDn->Fill(Edep_CND3, path, weight);
                    h_theta_n_miss_VS_Edep_CND3_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
                    h_ToF_VS_Edep_CND3_epFDn->Fill(Edep_CND3, ToF, weight);
                    h_nSector_VS_Edep_CND3_epFDn->Fill(Edep_CND3, nSector, weight);

                    h_ToF_epFDn->Fill(ToF, weight);
                    h_P_n_VS_ToF_epFDn->Fill(ToF, P_n.Mag(), weight);
                    h_theta_n_VS_ToF_epFDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_VS_ToF_epFDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                    h_P_miss_VS_ToF_epFDn->Fill(ToF, P_miss.Mag(), weight);
                    h_theta_miss_VS_ToF_epFDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                    h_phi_miss_VS_ToF_epFDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
                    h_dpp_VS_ToF_epFDn->Fill(ToF, dpp, weight);
                    h_beta_n_VS_ToF_epFDn->Fill(ToF, beta, weight);
                    h_E_miss_VS_ToF_epFDn->Fill(ToF, E_miss, weight);
                    h_M_miss_VS_ToF_epFDn->Fill(ToF, M_miss, weight);
                    h_path_VS_ToF_epFDn->Fill(ToF, path, weight);
                    h_theta_n_miss_VS_ToF_epFDn->Fill(ToF, theta_n_miss, weight);
                    h_nSector_VS_ToF_epFDn->Fill(ToF, nSector, weight);
                }

                //////////////////////////////////////////////
                // Step Zero
                //////////////////////////////////////////////

#pragma region /* Step Zero - start */

                // Why "path * 100"? unit conversion. Path is in cm; tof is in ns.
                // TODO: check if this unit conversion is needed!
                if (fabs(beta - (path * 100) / (ToF * c)) > 0.01) // A cut on delta beta
                // if (fabs(beta - path / (ToF * c)) > 0.01) // A cut on delta beta
                {
                    continue;
                }

                // A cut on the z-component of the CND hit
                // This is a fiducial cut on the range that the CND can reach on the z-axis
                if (v_hit.Z() > 45 || v_hit.Z() < -40)
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
                    h_xB_VS_M_miss_epCDn->Fill(xB, M_miss, weight);

                    h_dpp_allN_Step0_epCDn->Fill(dpp, weight);
                    h_theta_n_miss_allN_Step0_epCDn->Fill(theta_n_miss, weight);
                    h_dpp_VS_theta_n_miss_allN_Step0_epCDn->Fill(dpp, theta_n_miss, weight);

                    if (isGN)
                    {
                        h_theta_n_goodN_Step0_epCDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_goodN_Step0_epCDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                        h_theta_n_VS_phi_n_goodN_Step0_epCDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                        h_theta_n_VS_beta_n_goodN_Step0_epCDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                        h_P_n_goodN_Step0_epCDn->Fill(P_n.Mag(), weight);
                        h_P_n_VS_theta_n_goodN_Step0_epCDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                        h_P_miss_goodN_Step0_epCDn->Fill(P_miss.Mag(), weight);
                        h_P_miss_VS_theta_miss_goodN_Step0_epCDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);
                        h_P_miss_VS_phi_miss_goodN_Step0_epCDn->Fill(P_miss.Phi() * 180. / M_PI, P_miss.Mag(), weight);

                        h_dpp_goodN_Step0_epCDn->Fill(dpp, weight);
                        h_theta_n_miss_goodN_Step0_epCDn->Fill(theta_n_miss, weight);

                        h_E_p_goodN_Step0_epCDn->Fill(E_p, weight);
                        h_E_miss_goodN_Step0_epCDn->Fill(E_miss, weight);
                        h_M_miss_goodN_Step0_epCDn->Fill(M_miss, weight);
                        h_M_miss_VS_P_n_goodN_Step0_epCDn->Fill(P_n.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_n_goodN_Step0_epCDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_n_goodN_Step0_epCDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_P_miss_goodN_Step0_epCDn->Fill(P_miss.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_miss_goodN_Step0_epCDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_miss_goodN_Step0_epCDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                        h_P_n_minus_P_miss_goodN_Step0_epCDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                        h_P_n_x_minus_P_miss_x_goodN_Step0_epCDn->Fill(P_n.X() - P_miss.X(), weight);
                        h_P_n_y_minus_P_miss_y_goodN_Step0_epCDn->Fill(P_n.Y() - P_miss.Y(), weight);
                        h_P_n_z_minus_P_miss_z_goodN_Step0_epCDn->Fill(P_n.Z() - P_miss.Z(), weight);

                        h_P_n_VS_P_miss_goodN_Step0_epCDn->Fill(P_n.Mag(), P_miss.Mag(), weight);
                        h_P_n_x_VS_P_miss_x_goodN_Step0_epCDn->Fill(P_n.X(), P_miss.X(), weight);
                        h_P_n_y_VS_P_miss_y_goodN_Step0_epCDn->Fill(P_n.Y(), P_miss.Y(), weight);
                        h_P_n_z_VS_P_miss_z_goodN_Step0_epCDn->Fill(P_n.Z(), P_miss.Z(), weight);

                        h_theta_n_p_goodN_Step0_epCDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                        h_theta_n_p_VS_P_p_goodN_Step0_epCDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                        h_xB_goodN_Step0_epCDn->Fill(xB, weight);

                        h_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, weight);
                        h_P_n_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, dpp, weight);
                        h_beta_n_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, beta, weight);
                        h_E_p_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, E_p, weight);
                        h_E_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, E_miss, weight);
                        h_M_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, M_miss, weight);
                        h_path_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, ToF, weight);
                        h_nSector_VS_Edep_CND_goodN_Step0_epCDn->Fill(Edep_CND, nSector, weight);

                        h_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, weight);
                        h_P_n_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, dpp, weight);
                        h_beta_n_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, beta, weight);
                        h_E_p_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, E_p, weight);
                        h_E_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, E_miss, weight);
                        h_M_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, M_miss, weight);
                        h_path_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
                        h_ToF_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, ToF, weight);
                        h_nSector_VS_Edep_CTOF_goodN_Step0_epCDn->Fill(Edep_CTOF, nSector, weight);

                        h_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, weight);
                        h_P_n_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, dpp, weight);
                        h_beta_n_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, beta, weight);
                        h_E_p_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, E_p, weight);
                        h_E_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, E_miss, weight);
                        h_M_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, M_miss, weight);
                        h_path_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, ToF, weight);
                        h_nSector_VS_Edep_CND1_goodN_Step0_epCDn->Fill(Edep_CND1, nSector, weight);

                        h_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, weight);
                        h_P_n_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, dpp, weight);
                        h_beta_n_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, beta, weight);
                        h_E_p_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, E_p, weight);
                        h_E_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, E_miss, weight);
                        h_M_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, M_miss, weight);
                        h_path_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, ToF, weight);
                        h_nSector_VS_Edep_CND2_goodN_Step0_epCDn->Fill(Edep_CND2, nSector, weight);

                        h_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, weight);
                        h_P_n_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, dpp, weight);
                        h_beta_n_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, beta, weight);
                        h_E_p_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, E_p, weight);
                        h_E_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, E_miss, weight);
                        h_M_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, M_miss, weight);
                        h_path_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, ToF, weight);
                        h_nSector_VS_Edep_CND3_goodN_Step0_epCDn->Fill(Edep_CND3, nSector, weight);

                        h_ToF_goodN_Step0_epCDn->Fill(ToF, weight);
                        h_P_n_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_n.Mag(), weight);
                        h_theta_n_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_miss.Mag(), weight);
                        h_theta_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_ToF_goodN_Step0_epCDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
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
                        h_theta_n_badN_Step0_epCDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_badN_Step0_epCDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                        h_theta_n_VS_phi_n_badN_Step0_epCDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                        h_theta_n_VS_beta_n_badN_Step0_epCDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                        h_P_n_badN_Step0_epCDn->Fill(P_n.Mag(), weight);
                        h_P_n_VS_theta_n_badN_Step0_epCDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                        h_P_miss_badN_Step0_epCDn->Fill(P_miss.Mag(), weight);
                        h_P_miss_VS_theta_miss_badN_Step0_epCDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);
                        h_P_miss_VS_phi_miss_badN_Step0_epCDn->Fill(P_miss.Phi() * 180. / M_PI, P_miss.Mag(), weight);

                        h_dpp_badN_Step0_epCDn->Fill(dpp, weight);
                        h_theta_n_miss_badN_Step0_epCDn->Fill(theta_n_miss, weight);

                        h_E_p_badN_Step0_epCDn->Fill(E_p, weight);
                        h_E_miss_badN_Step0_epCDn->Fill(E_miss, weight);
                        h_M_miss_badN_Step0_epCDn->Fill(M_miss, weight);
                        h_M_miss_VS_P_n_badN_Step0_epCDn->Fill(P_n.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_n_badN_Step0_epCDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_n_badN_Step0_epCDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_P_miss_badN_Step0_epCDn->Fill(P_miss.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_miss_badN_Step0_epCDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_miss_badN_Step0_epCDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                        h_P_n_minus_P_miss_badN_Step0_epCDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                        h_P_n_x_minus_P_miss_x_badN_Step0_epCDn->Fill(P_n.X() - P_miss.X(), weight);
                        h_P_n_y_minus_P_miss_y_badN_Step0_epCDn->Fill(P_n.Y() - P_miss.Y(), weight);
                        h_P_n_z_minus_P_miss_z_badN_Step0_epCDn->Fill(P_n.Z() - P_miss.Z(), weight);

                        h_P_n_VS_P_miss_badN_Step0_epCDn->Fill(P_n.Mag(), P_miss.Mag(), weight);
                        h_P_n_x_VS_P_miss_x_badN_Step0_epCDn->Fill(P_n.X(), P_miss.X(), weight);
                        h_P_n_y_VS_P_miss_y_badN_Step0_epCDn->Fill(P_n.Y(), P_miss.Y(), weight);
                        h_P_n_z_VS_P_miss_z_badN_Step0_epCDn->Fill(P_n.Z(), P_miss.Z(), weight);

                        h_theta_n_p_badN_Step0_epCDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                        h_theta_n_p_VS_P_p_badN_Step0_epCDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                        h_xB_badN_Step0_epCDn->Fill(xB, weight);

                        h_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, weight);
                        h_P_n_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, dpp, weight);
                        h_beta_n_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, beta, weight);
                        h_E_p_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, E_p, weight);
                        h_E_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, E_miss, weight);
                        h_M_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, M_miss, weight);
                        h_path_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, ToF, weight);
                        h_nSector_VS_Edep_CND_badN_Step0_epCDn->Fill(Edep_CND, nSector, weight);

                        h_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, weight);
                        h_P_n_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, dpp, weight);
                        h_beta_n_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, beta, weight);
                        h_E_p_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, E_p, weight);
                        h_E_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, E_miss, weight);
                        h_M_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, M_miss, weight);
                        h_path_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
                        h_ToF_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, ToF, weight);
                        h_nSector_VS_Edep_CTOF_badN_Step0_epCDn->Fill(Edep_CTOF, nSector, weight);

                        h_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, weight);
                        h_P_n_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, dpp, weight);
                        h_beta_n_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, beta, weight);
                        h_E_p_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, E_p, weight);
                        h_E_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, E_miss, weight);
                        h_M_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, M_miss, weight);
                        h_path_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, ToF, weight);
                        h_nSector_VS_Edep_CND1_badN_Step0_epCDn->Fill(Edep_CND1, nSector, weight);

                        h_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, weight);
                        h_P_n_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, dpp, weight);
                        h_beta_n_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, beta, weight);
                        h_E_p_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, E_p, weight);
                        h_E_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, E_miss, weight);
                        h_M_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, M_miss, weight);
                        h_path_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, ToF, weight);
                        h_nSector_VS_Edep_CND2_badN_Step0_epCDn->Fill(Edep_CND2, nSector, weight);

                        h_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, weight);
                        h_P_n_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, dpp, weight);
                        h_beta_n_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, beta, weight);
                        h_E_p_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, E_p, weight);
                        h_E_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, E_miss, weight);
                        h_M_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, M_miss, weight);
                        h_path_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, ToF, weight);
                        h_nSector_VS_Edep_CND3_badN_Step0_epCDn->Fill(Edep_CND3, nSector, weight);

                        h_ToF_badN_Step0_epCDn->Fill(ToF, weight);
                        h_P_n_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_n.Mag(), weight);
                        h_theta_n_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_miss.Mag(), weight);
                        h_theta_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_ToF_badN_Step0_epCDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
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
                    h_xB_VS_M_miss_epFDn->Fill(xB, M_miss, weight);

                    h_dpp_allN_Step0_epFDn->Fill(dpp, weight);
                    h_theta_n_miss_allN_Step0_epFDn->Fill(theta_n_miss, weight);
                    h_dpp_VS_theta_n_miss_allN_Step0_epFDn->Fill(dpp, theta_n_miss, weight);

                    if (isGN)
                    {
                        h_theta_n_goodN_Step0_epFDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_goodN_Step0_epFDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                        h_theta_n_VS_phi_n_goodN_Step0_epFDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                        h_theta_n_VS_beta_n_goodN_Step0_epFDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                        h_P_n_goodN_Step0_epFDn->Fill(P_n.Mag(), weight);
                        h_P_n_VS_theta_n_goodN_Step0_epFDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                        h_P_miss_goodN_Step0_epFDn->Fill(P_miss.Mag(), weight);
                        h_P_miss_VS_theta_miss_goodN_Step0_epFDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);
                        h_P_miss_VS_phi_miss_goodN_Step0_epFDn->Fill(P_miss.Phi() * 180. / M_PI, P_miss.Mag(), weight);

                        h_dpp_goodN_Step0_epFDn->Fill(dpp, weight);
                        h_theta_n_miss_goodN_Step0_epFDn->Fill(theta_n_miss, weight);

                        h_E_p_goodN_Step0_epFDn->Fill(E_p, weight);
                        h_E_miss_goodN_Step0_epFDn->Fill(E_miss, weight);
                        h_M_miss_goodN_Step0_epFDn->Fill(M_miss, weight);
                        h_M_miss_VS_P_n_goodN_Step0_epFDn->Fill(P_n.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_n_goodN_Step0_epFDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_n_goodN_Step0_epFDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_P_miss_goodN_Step0_epFDn->Fill(P_miss.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_miss_goodN_Step0_epFDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_miss_goodN_Step0_epFDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                        h_P_n_minus_P_miss_goodN_Step0_epFDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                        h_P_n_x_minus_P_miss_x_goodN_Step0_epFDn->Fill(P_n.X() - P_miss.X(), weight);
                        h_P_n_y_minus_P_miss_y_goodN_Step0_epFDn->Fill(P_n.Y() - P_miss.Y(), weight);
                        h_P_n_z_minus_P_miss_z_goodN_Step0_epFDn->Fill(P_n.Z() - P_miss.Z(), weight);

                        h_P_n_VS_P_miss_goodN_Step0_epFDn->Fill(P_n.Mag(), P_miss.Mag(), weight);
                        h_P_n_x_VS_P_miss_x_goodN_Step0_epFDn->Fill(P_n.X(), P_miss.X(), weight);
                        h_P_n_y_VS_P_miss_y_goodN_Step0_epFDn->Fill(P_n.Y(), P_miss.Y(), weight);
                        h_P_n_z_VS_P_miss_z_goodN_Step0_epFDn->Fill(P_n.Z(), P_miss.Z(), weight);

                        h_theta_n_p_goodN_Step0_epFDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                        h_theta_n_p_VS_P_p_goodN_Step0_epFDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                        h_xB_goodN_Step0_epFDn->Fill(xB, weight);

                        h_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, weight);
                        h_P_n_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, dpp, weight);
                        h_beta_n_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, beta, weight);
                        h_E_p_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, E_p, weight);
                        h_E_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, E_miss, weight);
                        h_M_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, M_miss, weight);
                        h_path_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, ToF, weight);
                        h_nSector_VS_Edep_CND_goodN_Step0_epFDn->Fill(Edep_CND, nSector, weight);

                        h_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, weight);
                        h_P_n_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, dpp, weight);
                        h_beta_n_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, beta, weight);
                        h_E_p_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, E_p, weight);
                        h_E_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, E_miss, weight);
                        h_M_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, M_miss, weight);
                        h_path_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
                        h_ToF_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, ToF, weight);
                        h_nSector_VS_Edep_CTOF_goodN_Step0_epFDn->Fill(Edep_CTOF, nSector, weight);

                        h_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, weight);
                        h_P_n_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, dpp, weight);
                        h_beta_n_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, beta, weight);
                        h_E_p_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, E_p, weight);
                        h_E_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, E_miss, weight);
                        h_M_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, M_miss, weight);
                        h_path_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, ToF, weight);
                        h_nSector_VS_Edep_CND1_goodN_Step0_epFDn->Fill(Edep_CND1, nSector, weight);

                        h_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, weight);
                        h_P_n_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, dpp, weight);
                        h_beta_n_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, beta, weight);
                        h_E_p_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, E_p, weight);
                        h_E_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, E_miss, weight);
                        h_M_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, M_miss, weight);
                        h_path_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, ToF, weight);
                        h_nSector_VS_Edep_CND2_goodN_Step0_epFDn->Fill(Edep_CND2, nSector, weight);

                        h_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, weight);
                        h_P_n_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, dpp, weight);
                        h_beta_n_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, beta, weight);
                        h_E_p_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, E_p, weight);
                        h_E_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, E_miss, weight);
                        h_M_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, M_miss, weight);
                        h_path_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, ToF, weight);
                        h_nSector_VS_Edep_CND3_goodN_Step0_epFDn->Fill(Edep_CND3, nSector, weight);

                        h_ToF_goodN_Step0_epFDn->Fill(ToF, weight);
                        h_P_n_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_n.Mag(), weight);
                        h_theta_n_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_miss.Mag(), weight);
                        h_theta_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_ToF_goodN_Step0_epFDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
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
                        h_theta_n_badN_Step0_epFDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_badN_Step0_epFDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                        h_theta_n_VS_phi_n_badN_Step0_epFDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                        h_theta_n_VS_beta_n_badN_Step0_epFDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                        h_P_n_badN_Step0_epFDn->Fill(P_n.Mag(), weight);
                        h_P_n_VS_theta_n_badN_Step0_epFDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                        h_P_miss_badN_Step0_epFDn->Fill(P_miss.Mag(), weight);
                        h_P_miss_VS_theta_miss_badN_Step0_epFDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);
                        h_P_miss_VS_phi_miss_badN_Step0_epFDn->Fill(P_miss.Phi() * 180. / M_PI, P_miss.Mag(), weight);

                        h_dpp_badN_Step0_epFDn->Fill(dpp, weight);
                        h_theta_n_miss_badN_Step0_epFDn->Fill(theta_n_miss, weight);

                        h_E_p_badN_Step0_epFDn->Fill(E_p, weight);
                        h_E_miss_badN_Step0_epFDn->Fill(E_miss, weight);
                        h_M_miss_badN_Step0_epFDn->Fill(M_miss, weight);
                        h_M_miss_VS_P_n_badN_Step0_epFDn->Fill(P_n.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_n_badN_Step0_epFDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_n_badN_Step0_epFDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_P_miss_badN_Step0_epFDn->Fill(P_miss.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_miss_badN_Step0_epFDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_miss_badN_Step0_epFDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                        h_P_n_minus_P_miss_badN_Step0_epFDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                        h_P_n_x_minus_P_miss_x_badN_Step0_epFDn->Fill(P_n.X() - P_miss.X(), weight);
                        h_P_n_y_minus_P_miss_y_badN_Step0_epFDn->Fill(P_n.Y() - P_miss.Y(), weight);
                        h_P_n_z_minus_P_miss_z_badN_Step0_epFDn->Fill(P_n.Z() - P_miss.Z(), weight);

                        h_P_n_VS_P_miss_badN_Step0_epFDn->Fill(P_n.Mag(), P_miss.Mag(), weight);
                        h_P_n_x_VS_P_miss_x_badN_Step0_epFDn->Fill(P_n.X(), P_miss.X(), weight);
                        h_P_n_y_VS_P_miss_y_badN_Step0_epFDn->Fill(P_n.Y(), P_miss.Y(), weight);
                        h_P_n_z_VS_P_miss_z_badN_Step0_epFDn->Fill(P_n.Z(), P_miss.Z(), weight);

                        h_theta_n_p_badN_Step0_epFDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                        h_theta_n_p_VS_P_p_badN_Step0_epFDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                        h_xB_badN_Step0_epFDn->Fill(xB, weight);

                        h_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, weight);
                        h_P_n_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, dpp, weight);
                        h_beta_n_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, beta, weight);
                        h_E_p_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, E_p, weight);
                        h_E_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, E_miss, weight);
                        h_M_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, M_miss, weight);
                        h_path_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, ToF, weight);
                        h_nSector_VS_Edep_CND_badN_Step0_epFDn->Fill(Edep_CND, nSector, weight);

                        h_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, weight);
                        h_P_n_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, dpp, weight);
                        h_beta_n_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, beta, weight);
                        h_E_p_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, E_p, weight);
                        h_E_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, E_miss, weight);
                        h_M_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, M_miss, weight);
                        h_path_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
                        h_ToF_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, ToF, weight);
                        h_nSector_VS_Edep_CTOF_badN_Step0_epFDn->Fill(Edep_CTOF, nSector, weight);

                        h_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, weight);
                        h_P_n_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, dpp, weight);
                        h_beta_n_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, beta, weight);
                        h_E_p_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, E_p, weight);
                        h_E_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, E_miss, weight);
                        h_M_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, M_miss, weight);
                        h_path_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, ToF, weight);
                        h_nSector_VS_Edep_CND1_badN_Step0_epFDn->Fill(Edep_CND1, nSector, weight);

                        h_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, weight);
                        h_P_n_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, dpp, weight);
                        h_beta_n_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, beta, weight);
                        h_E_p_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, E_p, weight);
                        h_E_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, E_miss, weight);
                        h_M_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, M_miss, weight);
                        h_path_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, ToF, weight);
                        h_nSector_VS_Edep_CND2_badN_Step0_epFDn->Fill(Edep_CND2, nSector, weight);

                        h_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, weight);
                        h_P_n_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, dpp, weight);
                        h_beta_n_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, beta, weight);
                        h_E_p_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, E_p, weight);
                        h_E_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, E_miss, weight);
                        h_M_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, M_miss, weight);
                        h_path_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, ToF, weight);
                        h_nSector_VS_Edep_CND3_badN_Step0_epFDn->Fill(Edep_CND3, nSector, weight);

                        h_ToF_badN_Step0_epFDn->Fill(ToF, weight);
                        h_P_n_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_n.Mag(), weight);
                        h_theta_n_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_miss.Mag(), weight);
                        h_theta_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_ToF_badN_Step0_epFDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
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

                if (beta > 0.8) // Beta cut
                {
                    continue;
                }

                // Dep. energy cut
                // TODO: check if should be 12 MeV
                if (Edep_CND < 5)
                {
                    continue;
                }

                pass_step1_cuts = true;

                SetNeutronCounters(pInCD, pInFD, isGN, counter_n_multiplicity_allN_epCDn_Step1, counter_n_multiplicity_goodN_epCDn_Step1, counter_n_multiplicity_badN_epCDn_Step1,
                                   counter_n_multiplicity_allN_epFDn_Step1, counter_n_multiplicity_goodN_epFDn_Step1, counter_n_multiplicity_badN_epFDn_Step1);
                // SetNeutronCounters(isGN, counter_n_multiplicity_allN_Step1, counter_n_multiplicity_goodN_Step1, counter_n_multiplicity_badN_Step1);

                if (pInCD)
                {
                    h_xB_VS_M_miss_epCDn->Fill(xB, M_miss, weight);

                    h_dpp_allN_Step1_epCDn->Fill(dpp, weight);
                    h_theta_n_miss_allN_Step1_epCDn->Fill(theta_n_miss, weight);
                    h_dpp_VS_theta_n_miss_allN_Step1_epCDn->Fill(dpp, theta_n_miss, weight);

                    if (isGN)
                    {
                        h_theta_n_goodN_Step1_epCDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_goodN_Step1_epCDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                        h_theta_n_VS_phi_n_goodN_Step1_epCDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                        h_theta_n_VS_beta_n_goodN_Step1_epCDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                        h_P_n_goodN_Step1_epCDn->Fill(P_n.Mag(), weight);
                        h_P_n_VS_theta_n_goodN_Step1_epCDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                        h_P_miss_goodN_Step1_epCDn->Fill(P_miss.Mag(), weight);
                        h_P_miss_VS_theta_miss_goodN_Step1_epCDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);
                        h_P_miss_VS_phi_miss_goodN_Step1_epCDn->Fill(P_miss.Phi() * 180. / M_PI, P_miss.Mag(), weight);

                        h_dpp_goodN_Step1_epCDn->Fill(dpp, weight);
                        h_theta_n_miss_goodN_Step1_epCDn->Fill(theta_n_miss, weight);

                        h_E_p_goodN_Step1_epCDn->Fill(E_p, weight);
                        h_E_miss_goodN_Step1_epCDn->Fill(E_miss, weight);
                        h_M_miss_goodN_Step1_epCDn->Fill(M_miss, weight);
                        h_M_miss_VS_P_n_goodN_Step1_epCDn->Fill(P_n.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_n_goodN_Step1_epCDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_n_goodN_Step1_epCDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_P_miss_goodN_Step1_epCDn->Fill(P_miss.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_miss_goodN_Step1_epCDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_miss_goodN_Step1_epCDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                        h_P_n_minus_P_miss_goodN_Step1_epCDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                        h_P_n_x_minus_P_miss_x_goodN_Step1_epCDn->Fill(P_n.X() - P_miss.X(), weight);
                        h_P_n_y_minus_P_miss_y_goodN_Step1_epCDn->Fill(P_n.Y() - P_miss.Y(), weight);
                        h_P_n_z_minus_P_miss_z_goodN_Step1_epCDn->Fill(P_n.Z() - P_miss.Z(), weight);

                        h_P_n_VS_P_miss_goodN_Step1_epCDn->Fill(P_n.Mag(), P_miss.Mag(), weight);
                        h_P_n_x_VS_P_miss_x_goodN_Step1_epCDn->Fill(P_n.X(), P_miss.X(), weight);
                        h_P_n_y_VS_P_miss_y_goodN_Step1_epCDn->Fill(P_n.Y(), P_miss.Y(), weight);
                        h_P_n_z_VS_P_miss_z_goodN_Step1_epCDn->Fill(P_n.Z(), P_miss.Z(), weight);

                        h_theta_n_p_goodN_Step1_epCDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                        h_theta_n_p_VS_P_p_goodN_Step1_epCDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                        h_xB_goodN_Step1_epCDn->Fill(xB, weight);

                        h_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, weight);
                        h_P_n_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, dpp, weight);
                        h_beta_n_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, beta, weight);
                        h_E_p_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, E_p, weight);
                        h_E_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, E_miss, weight);
                        h_M_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, M_miss, weight);
                        h_path_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, ToF, weight);
                        h_nSector_VS_Edep_CND_goodN_Step1_epCDn->Fill(Edep_CND, nSector, weight);

                        h_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, weight);
                        h_P_n_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, dpp, weight);
                        h_beta_n_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, beta, weight);
                        h_E_p_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, E_p, weight);
                        h_E_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, E_miss, weight);
                        h_M_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, M_miss, weight);
                        h_path_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
                        h_ToF_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, ToF, weight);
                        h_nSector_VS_Edep_CTOF_goodN_Step1_epCDn->Fill(Edep_CTOF, nSector, weight);

                        h_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, weight);
                        h_P_n_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, dpp, weight);
                        h_beta_n_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, beta, weight);
                        h_E_p_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, E_p, weight);
                        h_E_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, E_miss, weight);
                        h_M_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, M_miss, weight);
                        h_path_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, ToF, weight);
                        h_nSector_VS_Edep_CND1_goodN_Step1_epCDn->Fill(Edep_CND1, nSector, weight);

                        h_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, weight);
                        h_P_n_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, dpp, weight);
                        h_beta_n_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, beta, weight);
                        h_E_p_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, E_p, weight);
                        h_E_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, E_miss, weight);
                        h_M_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, M_miss, weight);
                        h_path_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, ToF, weight);
                        h_nSector_VS_Edep_CND2_goodN_Step1_epCDn->Fill(Edep_CND2, nSector, weight);

                        h_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, weight);
                        h_P_n_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, dpp, weight);
                        h_beta_n_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, beta, weight);
                        h_E_p_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, E_p, weight);
                        h_E_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, E_miss, weight);
                        h_M_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, M_miss, weight);
                        h_path_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, ToF, weight);
                        h_nSector_VS_Edep_CND3_goodN_Step1_epCDn->Fill(Edep_CND3, nSector, weight);

                        h_ToF_goodN_Step1_epCDn->Fill(ToF, weight);
                        h_P_n_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_n.Mag(), weight);
                        h_theta_n_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_miss.Mag(), weight);
                        h_theta_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_ToF_goodN_Step1_epCDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
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
                        h_theta_n_badN_Step1_epCDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_badN_Step1_epCDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                        h_theta_n_VS_phi_n_badN_Step1_epCDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                        h_theta_n_VS_beta_n_badN_Step1_epCDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                        h_P_n_badN_Step1_epCDn->Fill(P_n.Mag(), weight);
                        h_P_n_VS_theta_n_badN_Step1_epCDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                        h_P_miss_badN_Step1_epCDn->Fill(P_miss.Mag(), weight);
                        h_P_miss_VS_theta_miss_badN_Step1_epCDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);
                        h_P_miss_VS_phi_miss_badN_Step1_epCDn->Fill(P_miss.Phi() * 180. / M_PI, P_miss.Mag(), weight);

                        h_dpp_badN_Step1_epCDn->Fill(dpp, weight);
                        h_theta_n_miss_badN_Step1_epCDn->Fill(theta_n_miss, weight);

                        h_E_p_badN_Step1_epCDn->Fill(E_p, weight);
                        h_E_miss_badN_Step1_epCDn->Fill(E_miss, weight);
                        h_M_miss_badN_Step1_epCDn->Fill(M_miss, weight);
                        h_M_miss_VS_P_n_badN_Step1_epCDn->Fill(P_n.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_n_badN_Step1_epCDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_n_badN_Step1_epCDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_P_miss_badN_Step1_epCDn->Fill(P_miss.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_miss_badN_Step1_epCDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_miss_badN_Step1_epCDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                        h_P_n_minus_P_miss_badN_Step1_epCDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                        h_P_n_x_minus_P_miss_x_badN_Step1_epCDn->Fill(P_n.X() - P_miss.X(), weight);
                        h_P_n_y_minus_P_miss_y_badN_Step1_epCDn->Fill(P_n.Y() - P_miss.Y(), weight);
                        h_P_n_z_minus_P_miss_z_badN_Step1_epCDn->Fill(P_n.Z() - P_miss.Z(), weight);

                        h_P_n_VS_P_miss_badN_Step1_epCDn->Fill(P_n.Mag(), P_miss.Mag(), weight);
                        h_P_n_x_VS_P_miss_x_badN_Step1_epCDn->Fill(P_n.X(), P_miss.X(), weight);
                        h_P_n_y_VS_P_miss_y_badN_Step1_epCDn->Fill(P_n.Y(), P_miss.Y(), weight);
                        h_P_n_z_VS_P_miss_z_badN_Step1_epCDn->Fill(P_n.Z(), P_miss.Z(), weight);

                        h_theta_n_p_badN_Step1_epCDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                        h_theta_n_p_VS_P_p_badN_Step1_epCDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                        h_xB_badN_Step1_epCDn->Fill(xB, weight);

                        h_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, weight);
                        h_P_n_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, dpp, weight);
                        h_beta_n_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, beta, weight);
                        h_E_p_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, E_p, weight);
                        h_E_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, E_miss, weight);
                        h_M_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, M_miss, weight);
                        h_path_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, ToF, weight);
                        h_nSector_VS_Edep_CND_badN_Step1_epCDn->Fill(Edep_CND, nSector, weight);

                        h_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, weight);
                        h_P_n_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, dpp, weight);
                        h_beta_n_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, beta, weight);
                        h_E_p_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, E_p, weight);
                        h_E_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, E_miss, weight);
                        h_M_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, M_miss, weight);
                        h_path_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, theta_n_miss, weight);
                        h_ToF_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, ToF, weight);
                        h_nSector_VS_Edep_CTOF_badN_Step1_epCDn->Fill(Edep_CTOF, nSector, weight);

                        h_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, weight);
                        h_P_n_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, dpp, weight);
                        h_beta_n_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, beta, weight);
                        h_E_p_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, E_p, weight);
                        h_E_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, E_miss, weight);
                        h_M_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, M_miss, weight);
                        h_path_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, ToF, weight);
                        h_nSector_VS_Edep_CND1_badN_Step1_epCDn->Fill(Edep_CND1, nSector, weight);

                        h_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, weight);
                        h_P_n_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, dpp, weight);
                        h_beta_n_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, beta, weight);
                        h_E_p_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, E_p, weight);
                        h_E_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, E_miss, weight);
                        h_M_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, M_miss, weight);
                        h_path_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, ToF, weight);
                        h_nSector_VS_Edep_CND2_badN_Step1_epCDn->Fill(Edep_CND2, nSector, weight);

                        h_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, weight);
                        h_P_n_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, dpp, weight);
                        h_beta_n_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, beta, weight);
                        h_E_p_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, E_p, weight);
                        h_E_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, E_miss, weight);
                        h_M_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, M_miss, weight);
                        h_path_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, ToF, weight);
                        h_nSector_VS_Edep_CND3_badN_Step1_epCDn->Fill(Edep_CND3, nSector, weight);

                        h_ToF_badN_Step1_epCDn->Fill(ToF, weight);
                        h_P_n_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_n.Mag(), weight);
                        h_theta_n_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_miss.Mag(), weight);
                        h_theta_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_ToF_badN_Step1_epCDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
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
                    h_xB_VS_M_miss_epFDn->Fill(xB, M_miss, weight);

                    h_dpp_allN_Step1_epFDn->Fill(dpp, weight);
                    h_theta_n_miss_allN_Step1_epFDn->Fill(theta_n_miss, weight);
                    h_dpp_VS_theta_n_miss_allN_Step1_epFDn->Fill(dpp, theta_n_miss, weight);

                    if (isGN)
                    {
                        h_theta_n_goodN_Step1_epFDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_goodN_Step1_epFDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                        h_theta_n_VS_phi_n_goodN_Step1_epFDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                        h_theta_n_VS_beta_n_goodN_Step1_epFDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                        h_P_n_goodN_Step1_epFDn->Fill(P_n.Mag(), weight);
                        h_P_n_VS_theta_n_goodN_Step1_epFDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                        h_P_miss_goodN_Step1_epFDn->Fill(P_miss.Mag(), weight);
                        h_P_miss_VS_theta_miss_goodN_Step1_epFDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);
                        h_P_miss_VS_phi_miss_goodN_Step1_epFDn->Fill(P_miss.Phi() * 180. / M_PI, P_miss.Mag(), weight);

                        h_dpp_goodN_Step1_epFDn->Fill(dpp, weight);
                        h_theta_n_miss_goodN_Step1_epFDn->Fill(theta_n_miss, weight);

                        h_E_p_goodN_Step1_epFDn->Fill(E_p, weight);
                        h_E_miss_goodN_Step1_epFDn->Fill(E_miss, weight);
                        h_M_miss_goodN_Step1_epFDn->Fill(M_miss, weight);
                        h_M_miss_VS_P_n_goodN_Step1_epFDn->Fill(P_n.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_n_goodN_Step1_epFDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_n_goodN_Step1_epFDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_P_miss_goodN_Step1_epFDn->Fill(P_miss.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_miss_goodN_Step1_epFDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_miss_goodN_Step1_epFDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                        h_P_n_minus_P_miss_goodN_Step1_epFDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                        h_P_n_x_minus_P_miss_x_goodN_Step1_epFDn->Fill(P_n.X() - P_miss.X(), weight);
                        h_P_n_y_minus_P_miss_y_goodN_Step1_epFDn->Fill(P_n.Y() - P_miss.Y(), weight);
                        h_P_n_z_minus_P_miss_z_goodN_Step1_epFDn->Fill(P_n.Z() - P_miss.Z(), weight);

                        h_P_n_VS_P_miss_goodN_Step1_epFDn->Fill(P_n.Mag(), P_miss.Mag(), weight);
                        h_P_n_x_VS_P_miss_x_goodN_Step1_epFDn->Fill(P_n.X(), P_miss.X(), weight);
                        h_P_n_y_VS_P_miss_y_goodN_Step1_epFDn->Fill(P_n.Y(), P_miss.Y(), weight);
                        h_P_n_z_VS_P_miss_z_goodN_Step1_epFDn->Fill(P_n.Z(), P_miss.Z(), weight);

                        h_theta_n_p_goodN_Step1_epFDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                        h_theta_n_p_VS_P_p_goodN_Step1_epFDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                        h_xB_goodN_Step1_epFDn->Fill(xB, weight);

                        h_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, weight);
                        h_P_n_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, dpp, weight);
                        h_beta_n_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, beta, weight);
                        h_E_p_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, E_p, weight);
                        h_E_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, E_miss, weight);
                        h_M_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, M_miss, weight);
                        h_path_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, ToF, weight);
                        h_nSector_VS_Edep_CND_goodN_Step1_epFDn->Fill(Edep_CND, nSector, weight);

                        h_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, weight);
                        h_P_n_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, dpp, weight);
                        h_beta_n_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, beta, weight);
                        h_E_p_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, E_p, weight);
                        h_E_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, E_miss, weight);
                        h_M_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, M_miss, weight);
                        h_path_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
                        h_ToF_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, ToF, weight);
                        h_nSector_VS_Edep_CTOF_goodN_Step1_epFDn->Fill(Edep_CTOF, nSector, weight);

                        h_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, weight);
                        h_P_n_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, dpp, weight);
                        h_beta_n_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, beta, weight);
                        h_E_p_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, E_p, weight);
                        h_E_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, E_miss, weight);
                        h_M_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, M_miss, weight);
                        h_path_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, ToF, weight);
                        h_nSector_VS_Edep_CND1_goodN_Step1_epFDn->Fill(Edep_CND1, nSector, weight);

                        h_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, weight);
                        h_P_n_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, dpp, weight);
                        h_beta_n_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, beta, weight);
                        h_E_p_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, E_p, weight);
                        h_E_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, E_miss, weight);
                        h_M_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, M_miss, weight);
                        h_path_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, ToF, weight);
                        h_nSector_VS_Edep_CND2_goodN_Step1_epFDn->Fill(Edep_CND2, nSector, weight);

                        h_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, weight);
                        h_P_n_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, dpp, weight);
                        h_beta_n_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, beta, weight);
                        h_E_p_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, E_p, weight);
                        h_E_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, E_miss, weight);
                        h_M_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, M_miss, weight);
                        h_path_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, ToF, weight);
                        h_nSector_VS_Edep_CND3_goodN_Step1_epFDn->Fill(Edep_CND3, nSector, weight);

                        h_ToF_goodN_Step1_epFDn->Fill(ToF, weight);
                        h_P_n_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_n.Mag(), weight);
                        h_theta_n_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_miss.Mag(), weight);
                        h_theta_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_ToF_goodN_Step1_epFDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
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
                        h_theta_n_badN_Step1_epFDn->Fill(P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_badN_Step1_epFDn->Fill(P_n.Phi() * 180. / M_PI, weight);
                        h_theta_n_VS_phi_n_badN_Step1_epFDn->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                        h_theta_n_VS_beta_n_badN_Step1_epFDn->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                        h_P_n_badN_Step1_epFDn->Fill(P_n.Mag(), weight);
                        h_P_n_VS_theta_n_badN_Step1_epFDn->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                        h_P_miss_badN_Step1_epFDn->Fill(P_miss.Mag(), weight);
                        h_P_miss_VS_theta_miss_badN_Step1_epFDn->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);
                        h_P_miss_VS_phi_miss_badN_Step1_epFDn->Fill(P_miss.Phi() * 180. / M_PI, P_miss.Mag(), weight);

                        h_dpp_badN_Step1_epFDn->Fill(dpp, weight);
                        h_theta_n_miss_badN_Step1_epFDn->Fill(theta_n_miss, weight);

                        h_E_p_badN_Step1_epFDn->Fill(E_p, weight);
                        h_E_miss_badN_Step1_epFDn->Fill(E_miss, weight);
                        h_M_miss_badN_Step1_epFDn->Fill(M_miss, weight);
                        h_M_miss_VS_P_n_badN_Step1_epFDn->Fill(P_n.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_n_badN_Step1_epFDn->Fill(P_n.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_n_badN_Step1_epFDn->Fill(P_n.Phi() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_P_miss_badN_Step1_epFDn->Fill(P_miss.Mag(), M_miss, weight);
                        h_M_miss_VS_theta_miss_badN_Step1_epFDn->Fill(P_miss.Theta() * 180. / M_PI, M_miss, weight);
                        h_M_miss_VS_phi_miss_badN_Step1_epFDn->Fill(P_miss.Phi() * 180. / M_PI, M_miss, weight);

                        h_P_n_minus_P_miss_badN_Step1_epFDn->Fill(P_n.Mag() - P_miss.Mag(), weight);
                        h_P_n_x_minus_P_miss_x_badN_Step1_epFDn->Fill(P_n.X() - P_miss.X(), weight);
                        h_P_n_y_minus_P_miss_y_badN_Step1_epFDn->Fill(P_n.Y() - P_miss.Y(), weight);
                        h_P_n_z_minus_P_miss_z_badN_Step1_epFDn->Fill(P_n.Z() - P_miss.Z(), weight);

                        h_P_n_VS_P_miss_badN_Step1_epFDn->Fill(P_n.Mag(), P_miss.Mag(), weight);
                        h_P_n_x_VS_P_miss_x_badN_Step1_epFDn->Fill(P_n.X(), P_miss.X(), weight);
                        h_P_n_y_VS_P_miss_y_badN_Step1_epFDn->Fill(P_n.Y(), P_miss.Y(), weight);
                        h_P_n_z_VS_P_miss_z_badN_Step1_epFDn->Fill(P_n.Z(), P_miss.Z(), weight);

                        h_theta_n_p_badN_Step1_epFDn->Fill(P_p.Angle(P_n) * 180. / M_PI, weight);
                        h_theta_n_p_VS_P_p_badN_Step1_epFDn->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                        h_xB_badN_Step1_epFDn->Fill(xB, weight);

                        h_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, weight);
                        h_P_n_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, dpp, weight);
                        h_beta_n_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, beta, weight);
                        h_E_p_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, E_p, weight);
                        h_E_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, E_miss, weight);
                        h_M_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, M_miss, weight);
                        h_path_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, ToF, weight);
                        h_nSector_VS_Edep_CND_badN_Step1_epFDn->Fill(Edep_CND, nSector, weight);

                        h_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, weight);
                        h_P_n_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, dpp, weight);
                        h_beta_n_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, beta, weight);
                        h_E_p_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, E_p, weight);
                        h_E_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, E_miss, weight);
                        h_M_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, M_miss, weight);
                        h_path_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, theta_n_miss, weight);
                        h_ToF_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, ToF, weight);
                        h_nSector_VS_Edep_CTOF_badN_Step1_epFDn->Fill(Edep_CTOF, nSector, weight);

                        h_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, weight);
                        h_P_n_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, dpp, weight);
                        h_beta_n_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, beta, weight);
                        h_E_p_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, E_p, weight);
                        h_E_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, E_miss, weight);
                        h_M_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, M_miss, weight);
                        h_path_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, ToF, weight);
                        h_nSector_VS_Edep_CND1_badN_Step1_epFDn->Fill(Edep_CND1, nSector, weight);

                        h_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, weight);
                        h_P_n_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, dpp, weight);
                        h_beta_n_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, beta, weight);
                        h_E_p_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, E_p, weight);
                        h_E_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, E_miss, weight);
                        h_M_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, M_miss, weight);
                        h_path_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, ToF, weight);
                        h_nSector_VS_Edep_CND2_badN_Step1_epFDn->Fill(Edep_CND2, nSector, weight);

                        h_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, weight);
                        h_P_n_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_n.Mag(), weight);
                        h_theta_n_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_miss.Mag(), weight);
                        h_theta_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, P_miss.Phi() * 180. / M_PI, weight);
                        h_dpp_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, dpp, weight);
                        h_beta_n_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, beta, weight);
                        h_E_p_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, E_p, weight);
                        h_E_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, E_miss, weight);
                        h_M_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, M_miss, weight);
                        h_path_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, path * 100, weight);
                        h_theta_n_miss_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, theta_n_miss, weight);
                        h_ToF_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, ToF, weight);
                        h_nSector_VS_Edep_CND3_badN_Step1_epFDn->Fill(Edep_CND3, nSector, weight);

                        h_ToF_badN_Step1_epFDn->Fill(ToF, weight);
                        h_P_n_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_n.Mag(), weight);
                        h_theta_n_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_n.Theta() * 180. / M_PI, weight);
                        h_phi_n_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_n.Phi() * 180. / M_PI, weight);
                        h_P_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_miss.Mag(), weight);
                        h_theta_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_miss.Theta() * 180. / M_PI, weight);
                        h_phi_miss_VS_ToF_badN_Step1_epFDn->Fill(ToF, P_miss.Phi() * 180. / M_PI, weight);
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

                if (ToF * c - v_hit.Z() < 70) // TODO: find a way to check what is this cut
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

                        TVector3 p_C; // Momentum of the charged particle in the itr2-th entry of AllParticles
                        p_C.SetMagThetaPhi(AllParticles[itr2]->getP(), AllParticles[itr2]->getTheta(), AllParticles[itr2]->getPhi());

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
                                    h_sdiff_pos_mom_goodN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, p_C.Perp(), weight);
                                    h_sdiff_pos_z_goodN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, v_hit.Z(), weight);
                                    h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, ToF * c - v_hit.Z(), weight);
                                }
                                else
                                {
                                    h_sdiff_pos_badN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, weight);
                                    h_sdiff_pos_mom_badN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, p_C.Perp(), weight);
                                    h_sdiff_pos_z_badN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, v_hit.Z(), weight);
                                    h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epCDn[ldiff + 3]->Fill(sdiff, ToF * c - v_hit.Z(), weight);
                                }
                            }
                            else if (pInFD)
                            {
                                if (isGN) // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step1_layer
                                {
                                    h_sdiff_pos_goodN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, weight);
                                    h_sdiff_pos_mom_goodN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, p_C.Perp(), weight);
                                    h_sdiff_pos_z_goodN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, v_hit.Z(), weight);
                                    h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, ToF * c - v_hit.Z(), weight);
                                }
                                else
                                {
                                    h_sdiff_pos_badN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, weight);
                                    h_sdiff_pos_mom_badN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, p_C.Perp(), weight);
                                    h_sdiff_pos_z_badN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, v_hit.Z(), weight);
                                    h_sdiff_pos_diff_ToFc_z_badN_Step1_layer_epFDn[ldiff + 3]->Fill(sdiff, ToF * c - v_hit.Z(), weight);
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
                                h_diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epCDn->Fill(ToF * c - v_hit.Z(), Edep_CND, weight);
                            else
                            {
                                h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epCDn->Fill(ToF * c - v_hit.Z(), Edep_CND, weight);
                            }
                        }
                        else
                        {
                            if (!CNDVeto)
                                h_diff_ToFc_z_VS_Edep_noNear_badN_Step1_epCDn->Fill(ToF * c - v_hit.Z(), Edep_CND, weight);
                            else
                            {
                                h_diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epCDn->Fill(ToF * c - v_hit.Z(), Edep_CND, weight);
                            }
                        }
                    }
                    else if (pInFD)
                    {
                        if (isGN)
                        {
                            if (!CNDVeto)
                                h_diff_ToFc_z_VS_Edep_noNear_goodN_Step1_epFDn->Fill(ToF * c - v_hit.Z(), Edep_CND, weight);
                            else
                            {
                                h_diff_ToFc_z_VS_Edep_yesNear_goodN_Step1_epFDn->Fill(ToF * c - v_hit.Z(), Edep_CND, weight);
                            }
                        }
                        else
                        {
                            if (!CNDVeto)
                                h_diff_ToFc_z_VS_Edep_noNear_badN_Step1_epFDn->Fill(ToF * c - v_hit.Z(), Edep_CND, weight);
                            else
                            {
                                h_diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epFDn->Fill(ToF * c - v_hit.Z(), Edep_CND, weight);
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
                    h_theta_n_goodN_Step2->Fill(P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_goodN_Step2->Fill(P_n.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_goodN_Step2->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_goodN_Step2->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                    h_P_n_goodN_Step2->Fill(P_n.Mag(), weight);
                    h_P_n_VS_theta_n_goodN_Step2->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                    h_P_miss_goodN_Step2->Fill(P_miss.Mag(), weight);
                    h_P_miss_VS_theta_miss_goodN_Step2->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);

                    h_P_n_minus_P_miss_goodN_Step2->Fill(P_n.Mag() - P_miss.Mag(), weight);
                    h_P_n_x_minus_P_miss_x_goodN_Step2->Fill(P_n.X() - P_miss.X(), weight);
                    h_P_n_y_minus_P_miss_y_goodN_Step2->Fill(P_n.Y() - P_miss.Y(), weight);
                    h_P_n_z_minus_P_miss_z_goodN_Step2->Fill(P_n.Z() - P_miss.Z(), weight);

                    h_P_n_VS_P_miss_goodN_Step2->Fill(P_miss.Mag(), P_n.Mag(), weight);
                    h_P_n_x_VS_P_miss_x_goodN_Step2->Fill(P_n.X() - P_miss.X(), weight);
                    h_P_n_y_VS_P_miss_y_goodN_Step2->Fill(P_n.Y() - P_miss.Y(), weight);
                    h_P_n_z_VS_P_miss_z_goodN_Step2->Fill(P_n.Z() - P_miss.Z(), weight);

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
                    h_M_miss_VS_P_n_goodN_Step2->Fill(P_n.Mag(), M_miss, weight);
                    h_M_miss_VS_P_miss_goodN_Step2->Fill(P_miss.Mag(), M_miss, weight);

                    h_theta_n_p_VS_P_p_goodN_Step2->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                    h_xB_goodN_Step2->Fill(xB, weight);

                    h_Edep_goodN_Step2->Fill(edep, weight);
                    h_P_n_VS_Edep_goodN_Step2->Fill(edep, P_n.Mag(), weight);
                    h_P_miss_VS_Edep_goodN_Step2->Fill(edep, P_miss.Mag(), weight);

                    h_dpp_VS_Edep_goodN_Step2->Fill(edep, dpp, weight);

                    h_ToF_goodN_Step2->Fill(ToF, weight);
                    h_pmiss_goodN_Step2->Fill(P_miss.Mag(), weight);
                    h_Edep_ToF_goodN_Step2->Fill(ToF, edep, weight);

                    h_ToF_goodN_Step2->Fill(ToF, weight);
                    h_Edep_ToF_goodN_Step2->Fill(ToF, edep, weight);
                }
                else
                {
                    h_theta_n_badN_Step2->Fill(P_n.Theta() * 180. / M_PI, weight);
                    h_phi_n_badN_Step2->Fill(P_n.Phi() * 180. / M_PI, weight);
                    h_theta_n_VS_phi_n_badN_Step2->Fill(P_n.Phi() * 180. / M_PI, P_n.Theta() * 180. / M_PI, weight);
                    h_theta_n_VS_beta_n_badN_Step2->Fill(beta, P_n.Theta() * 180. / M_PI, weight);

                    h_P_n_badN_Step2->Fill(P_n.Mag(), weight);
                    h_P_n_VS_theta_n_badN_Step2->Fill(P_n.Theta() * 180. / M_PI, P_n.Mag(), weight);

                    h_P_miss_badN_Step2->Fill(P_miss.Mag(), weight);
                    h_P_miss_VS_theta_miss_badN_Step2->Fill(P_miss.Theta() * 180. / M_PI, P_miss.Mag(), weight);

                    h_P_n_minus_P_miss_badN_Step2->Fill(P_n.Mag() - P_miss.Mag(), weight);
                    h_P_n_x_minus_P_miss_x_badN_Step2->Fill(P_n.X() - P_miss.X(), weight);
                    h_P_n_y_minus_P_miss_y_badN_Step2->Fill(P_n.Y() - P_miss.Y(), weight);
                    h_P_n_z_minus_P_miss_z_badN_Step2->Fill(P_n.Z() - P_miss.Z(), weight);

                    h_P_n_VS_P_miss_badN_Step2->Fill(P_miss.Mag(), P_n.Mag(), weight);
                    h_P_n_x_VS_P_miss_x_badN_Step2->Fill(P_n.X() - P_miss.X(), weight);
                    h_P_n_y_VS_P_miss_y_badN_Step2->Fill(P_n.Y() - P_miss.Y(), weight);
                    h_P_n_z_VS_P_miss_z_badN_Step2->Fill(P_n.Z() - P_miss.Z(), weight);

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
                    h_M_miss_VS_P_n_badN_Step2->Fill(P_n.Mag(), M_miss, weight);
                    h_M_miss_VS_P_miss_badN_Step2->Fill(P_miss.Mag(), M_miss, weight);

                    h_theta_n_p_VS_P_p_badN_Step2->Fill(P_p.Mag(), P_p.Angle(P_n) * 180. / M_PI, weight);

                    h_xB_badN_Step2->Fill(xB, weight);

                    h_Edep_badN_Step2->Fill(edep, weight);
                    h_P_n_VS_Edep_badN_Step2->Fill(edep, P_n.Mag(), weight);
                    h_P_miss_VS_Edep_badN_Step2->Fill(edep, P_miss.Mag(), weight);

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

                    // if(ToF*c-v_hit.Z() < 70){
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
                // if(ToF*c-v_hit.Z() < 70){
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
                h_pmiss_allN_Step5->Fill(p_miss.Mag(), weight);
                if (isGN)
                {
                    h_ToF_goodN_Step5->Fill(ToF, weight);
                    h_edep_ToF_goodN_Step5->Fill(ToF, edep, weight);
                    h_pmiss_goodN_Step5->Fill(p_miss.Mag(), weight);
                    h_diff_ToFc_z_Edep_goodN_Step5->Fill(ToF * c - v_hit.Z(), edep, weight);
                    h_diff_ToFc_z_Edep_goodN_Step5_layer_epCDn[detINTlayer - 1]->Fill(ToF * c - v_hit.Z(), edep, weight);
                    h_phidiff_en_goodN_Step5->Fill(get_phi_diff(p_e, P_n), weight);
                    h_TP_goodN_Step5->Fill(ToF / path, weight);
                    h_Z_goodN_Step5->Fill(v_hit.Z(), weight);
                    h_beta_Edep_goodN_Step5->Fill(beta, edep, weight);

                    h_ToF_Edep_goodN_Step5->Fill(ToF, edep, weight);
                    h_TP_Edep_goodN_Step5->Fill(ToF / path, edep, weight);
                }
                else
                {
                    h_ToF_badN_Step5->Fill(ToF, weight);
                    h_edep_ToF_badN_Step5->Fill(ToF, edep, weight);
                    // if(ToF<8){
                    h_diff_ToFc_z_Edep_badN_Step5->Fill(ToF * c - v_hit.Z(), edep, weight);
                    h_diff_ToFc_z_Edep_badN_Step5_layer_epCDn[detINTlayer - 1]->Fill(ToF * c - v_hit.Z(), edep, weight);
                    h_phidiff_en_badN_Step5->Fill(get_phi_diff(p_e, P_n), weight);
                    h_TP_badN_Step5->Fill(ToF / path, weight);
                    h_Z_badN_Step5->Fill(v_hit.Z(), weight);
                    h_beta_Edep_badN_Step5->Fill(beta, edep, weight);

                    h_ToF_Edep_badN_Step5->Fill(ToF, edep, weight);
                    h_TP_Edep_badN_Step5->Fill(ToF / path, edep, weight);
                    //}
                    if (ToF < 5)
                    {
                        if (edep > 20)
                        {
                            // cerr<<"Event="<<c12->runconfig()->getEvent()<<endl;
                            // cerr<<"Neutron Sector = "<<nSector<<endl;
                            // cerr<<"Neutron Z Hit = "<<v_hit.Z()<<endl;
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
                    h_ToF_z_badN->Fill(ToF, v_hit.Z(), weight);
                    if (ToF < 10)
                    {
                        h_TM_badN->Fill(ToF / path, weight);
                        h_beta_badN->Fill(beta, weight);
                        h_mom_badN->Fill(P_n.Mag(), weight);
                        h_Edep_z_badN->Fill(Edep_single, v_hit.Z(), weight);
                        h_Edep_ToF_badN->Fill(Edep_single, ToF, weight);
                        h_beta_z_badN->Fill(beta, v_hit.Z(), weight);
                        if (C1)
                        {
                            h_ToFc_z_1_badN->Fill(ToF * c, v_hit.Z(), weight);
                            h_diff_ToFc_z_1_badN->Fill(ToF * c - v_hit.Z(), weight);
                            h_diff_ToFc_z_Edep_1_badN->Fill(ToF * c - v_hit.Z(), edep, weight);
                            h_Edep_z_1_badN->Fill(Edep_single, v_hit.Z(), weight);
                        }
                        if (C2)
                        {
                            h_ToFc_z_2_badN->Fill(ToF * c, v_hit.Z(), weight);
                            h_diff_ToFc_z_2_badN->Fill(ToF * c - v_hit.Z(), weight);
                            h_diff_ToFc_z_Edep_2_badN->Fill(ToF * c - v_hit.Z(), edep, weight);
                            h_Edep_z_2_badN->Fill(Edep_single, v_hit.Z(), weight);
                        }
                        if (C3)
                        {
                            h_ToFc_z_3_badN->Fill(ToF * c, v_hit.Z(), weight);
                            h_diff_ToFc_z_3_badN->Fill(ToF * c - v_hit.Z(), weight);
                            h_diff_ToFc_z_Edep_3_badN->Fill(ToF * c - v_hit.Z(), edep, weight);
                            h_Edep_z_3_badN->Fill(Edep_single, v_hit.Z(), weight);
                        }

                        h_Edep_mom_badN->Fill(edep, mom, weight);

                        if (v_hit.Z() > 10)
                        {
                            cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
                            cerr << "Run=" << c12->runconfig()->getRun() << endl;
                            cerr << "Event=" << c12->runconfig()->getEvent() << endl;
                            cerr << "Neutron Sector = " << nSector << endl;
                            cerr << "Neutron Z Hit = " << v_hit.Z() << endl;
                            cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
                                 << endl
                                 << endl;
                        }
                    }
                }
                else
                {
                    h_ToF_goodN->Fill(ToF, weight);
                    h_ToF_z_goodN->Fill(ToF, v_hit.Z(), weight);
                    if (ToF < 10)
                    {
                        h_TM_goodN->Fill(ToF / path, weight);
                        h_beta_goodN->Fill(beta, weight);
                        h_mom_goodN->Fill(P_n.Mag(), weight);
                        h_Edep_z_goodN->Fill(Edep_single, v_hit.Z(), weight);
                        h_Edep_ToF_goodN->Fill(Edep_single, ToF, weight);
                        h_beta_z_goodN->Fill(beta, v_hit.Z(), weight);
                        if (C1)
                        {
                            h_ToFc_z_1_goodN->Fill(ToF * c, v_hit.Z(), weight);
                            h_diff_ToFc_z_1_goodN->Fill(ToF * c - v_hit.Z(), weight);
                            h_diff_ToFc_z_Edep_1_goodN->Fill(ToF * c - v_hit.Z(), edep, weight);
                            h_Edep_z_1_goodN->Fill(Edep_single, v_hit.Z(), weight);
                        }
                        if (C2)
                        {
                            h_ToFc_z_2_goodN->Fill(ToF * c, v_hit.Z(), weight);
                            h_diff_ToFc_z_2_goodN->Fill(ToF * c - v_hit.Z(), weight);
                            h_diff_ToFc_z_Edep_2_goodN->Fill(ToF * c - v_hit.Z(), edep, weight);
                            h_Edep_z_2_goodN->Fill(Edep_single, v_hit.Z(), weight);
                        }
                        if (C3)
                        {
                            h_ToFc_z_3_goodN->Fill(ToF * c, v_hit.Z(), weight);
                            h_diff_ToFc_z_3_goodN->Fill(ToF * c - v_hit.Z(), weight);
                            h_diff_ToFc_z_Edep_3_goodN->Fill(ToF * c - v_hit.Z(), edep, weight);
                            h_Edep_z_3_goodN->Fill(Edep_single, v_hit.Z(), weight);
                        }
                    }
                    h_Edep_mom_goodN->Fill(edep, mom, weight);
                }

                if (edep < 12.5)
                {
                    continue;
                }
                if ((edep < -(40.0 / 110.0) * ((ToF * c - v_hit.Z()) - 110)) && C1)
                {
                    continue;
                }
                if ((edep < -(32.0 / 110.0) * ((ToF * c - v_hit.Z()) - 110)) && C2)
                {
                    continue;
                }
                if ((edep < -(26.0 / 110.0) * ((ToF * c - v_hit.Z()) - 110)) && C3)
                {
                    continue;
                }
                */

                // if (C3 && (v_hit.Z() > 25))
                // {
                //     continue;
                // }
                // else if (C2 && (v_hit.Z() > 20))
                // {
                //     continue;
                // }
                // else if (C1 && (v_hit.Z() > 10))
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
        }

#pragma endregion /* Andrew's manual work */

    } // closes event loop

#pragma endregion /* Chain loop - end */

    // ======================================================================================================================================================================
    // Andrew's wrap up
    // ======================================================================================================================================================================

    HistPrinter(HistoList, PDFFile);

#pragma region /* Andrew's wrap up - start */

    //     /////////////////////////////////////////////////////
    //     // Now create the output PDFs
    //     /////////////////////////////////////////////////////

    //     int pixelx = 1980, pixely = 1530;
    //     // int pixelx = 1000 * 4 * 1.5 * 2, pixely = 750 * 3 * 1.5 * 2;

    //     TCanvas *myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);
    //     TCanvas *myText = new TCanvas("myText", "myText", pixelx, pixely);

    // #pragma region /* Saving all plots - start */

    //     TLatex text;
    //     text.SetTextSize(0.05);

    //     const char *pdfFile = PDFFile.c_str();

    //     char fileName[100];
    //     sprintf(fileName, "%s[", pdfFile);
    //     myText->SaveAs(fileName);
    //     sprintf(fileName, "%s", pdfFile);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text.DrawLatex(0.2, 0.7, "Neutrons in CND");

    //     myText->Print(fileName, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     // myCanvas->SetGrid();
    //     myCanvas->Divide(4, 3);
    //     // myCanvas->SetGrid(), myCanvas->cd()->SetBottomMargin(0.14), myCanvas->cd()->SetLeftMargin(0.16), myCanvas->cd()->SetRightMargin(0.16), myCanvas->cd()->SetTopMargin(0.12);

    //     double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    //     double diplayTextSize = 0.1;

    //     // int canvas_ind = 1;

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         int canvas_ind = (i % 12) + 1; // Determine the pad number (1 to 12)

    //         myCanvas->cd(canvas_ind);
    //         gPad->SetGrid();

    //         HistoList[i]->SetLineWidth(2);
    //         HistoList[i]->SetLineColor(kBlue);

    //         if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //         {
    //             TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //             displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //             HistoList[i]->Draw(), displayText->Draw("same");
    //         }
    //         else
    //         {
    //             HistoList[i]->Draw();
    //         }

    //         // Save the canvas to a PDF page after filling 12 pads or processing the last histogram
    //         if (canvas_ind == 12 || i == HistoList.size() - 1)
    //         {
    //             myCanvas->Print(fileName); // Save the current page
    //             if (i != HistoList.size() - 1)
    //             {
    //                 myCanvas->Clear();      // Clear the canvas for the next page
    //                 myCanvas->Divide(4, 3); // Reset the grid layout
    //             }
    //         }
    //     }

    //     // for (int i = 0; i < HistoList.size(); i++)
    //     // {
    //     //     // myCanvas->cd(canvas_ind);
    //     //     // gPad->SetGrid();

    //     //     HistoList[i]->SetLineWidth(2);
    //     //     HistoList[i]->SetLineColor(kBlue);

    //     //     if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //     //     {
    //     //         TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //     //         displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //     //         HistoList[i]->Draw(), displayText->Draw("same");
    //     //     }
    //     //     else
    //     //     {
    //     //         HistoList[i]->Draw();
    //     //     }

    //     //     myCanvas->Print(fileName, "pdf");
    //     //     myCanvas->Clear();

    //     //     // ++canvas_ind;

    //     //     // if (i > 12 && 12 % i == 0)
    //     //     // {
    //     //     //     myCanvas->Print(fileName, "pdf");
    //     //     //     myCanvas->Clear();
    //     //     //     myCanvas->Divide(4, 3);

    //     //     //     canvas_ind = 1;
    //     //     // }
    //     // }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         myCanvas->cd(1);

    //         if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //         {
    //             TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //             displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //             HistoList[i]->Draw("colz"), displayText->Draw("same");
    //         }
    //         else
    //         {
    //             HistoList[i]->Draw("colz");
    //         }

    //         myCanvas->Print(fileName, "pdf");
    //         myCanvas->Clear();
    //     }

    //     sprintf(fileName, "%s]", pdfFile);
    //     myCanvas->Print(fileName, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving all plots - end */

    // #pragma region /* Saving only CD proton plots - start */

    //     TLatex text_CD;
    //     text_CD.SetTextSize(0.05);

    //     string pdfFile_CD_0 = ConfigOutPutName(PDFFile, "pCD_only").c_str();
    //     const char *pdfFile_CD = pdfFile_CD_0.c_str();

    //     char fileName_CD[100];
    //     sprintf(fileName_CD, "%s[", pdfFile_CD);
    //     myText->SaveAs(fileName_CD);
    //     sprintf(fileName_CD, "%s", pdfFile_CD);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_CD.DrawLatex(0.2, 0.9, "(e,e'pCD) Cuts:");
    //     text_CD.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_CD.DrawLatex(0.2, 0.7, "Neutrons in CND");

    //     myText->Print(fileName_CD, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();
    //     // myCanvas->Divide(4, 3);
    //     // myCanvas->SetGrid(), myCanvas->cd()->SetBottomMargin(0.14), myCanvas->cd()->SetLeftMargin(0.16), myCanvas->cd()->SetRightMargin(0.16), myCanvas->cd()->SetTopMargin(0.12);

    //     // double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    //     // double diplayTextSize = 0.1;

    //     // int canvas_ind = 1;

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "CD"))
    //         {
    //             // myCanvas->cd(canvas_ind);
    //             // gPad->SetGrid();

    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_CD, "pdf");
    //             myCanvas->Clear();

    //             // ++canvas_ind;

    //             // if (i > 12 && 12 % i == 0)
    //             // {
    //             //     myCanvas->Print(fileName_CD, "pdf");
    //             //     myCanvas->Clear();
    //             //     myCanvas->Divide(4, 3);

    //             //     canvas_ind = 1;
    //             // }
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "CD"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_CD, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_CD, "%s]", pdfFile_CD);
    //     myCanvas->Print(fileName_CD, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving only CD proton plots - end */

    // #pragma region /* Saving Step0 plots - start */

    //     TLatex text_Step0;
    //     text_Step0.SetTextSize(0.05);

    //     string pdfFile_Step0_0 = ConfigOutPutName(PDFFile, "Step0").c_str();
    //     const char *pdfFile_Step0 = pdfFile_Step0_0.c_str();

    //     cout << "\nPDFFile = " << PDFFile << "\n";
    //     cout << "\npdfFile_Step0_0 = " << pdfFile_Step0_0 << "\n";
    //     cout << "\npdfFile_Step0 = " << pdfFile_Step0 << "\n";

    //     char fileName_Step0[100];
    //     sprintf(fileName_Step0, "%s[", pdfFile_Step0);
    //     myText->SaveAs(fileName_Step0);
    //     sprintf(fileName_Step0, "%s", pdfFile_Step0);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step0.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step0.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step0.DrawLatex(0.2, 0.7, "Neutrons in CND - step 0");

    //     myText->Print(fileName_Step0, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step0"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step0, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step0"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step0, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step0, "%s]", pdfFile_Step0);
    //     myCanvas->Print(fileName_Step0, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step0 plots - end */

    // #pragma region /* Saving Step1 plots - start */

    //     TLatex text_Step1;
    //     text_Step1.SetTextSize(0.05);

    //     string pdfFile_Step1_0 = ConfigOutPutName(PDFFile, "Step1").c_str();
    //     const char *pdfFile_Step1 = pdfFile_Step1_0.c_str();

    //     char fileName_Step1[100];
    //     sprintf(fileName_Step1, "%s[", pdfFile_Step1);
    //     myText->SaveAs(fileName_Step1);
    //     sprintf(fileName_Step1, "%s", pdfFile_Step1);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step1.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step1.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step1.DrawLatex(0.2, 0.7, "Neutrons in CND - step 1");

    //     myText->Print(fileName_Step1, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step1"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step1, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step1"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step1, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step1, "%s]", pdfFile_Step1);
    //     myCanvas->Print(fileName_Step1, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step1 plots - end */

    // #pragma region /* Saving Step2 plots - start */

    //     TLatex text_Step2;
    //     text_Step2.SetTextSize(0.05);

    //     string pdfFile_Step2_0 = ConfigOutPutName(PDFFile, "Step2").c_str();
    //     const char *pdfFile_Step2 = pdfFile_Step2_0.c_str();

    //     char fileName_Step2[100];
    //     sprintf(fileName_Step2, "%s[", pdfFile_Step2);
    //     myText->SaveAs(fileName_Step2);
    //     sprintf(fileName_Step2, "%s", pdfFile_Step2);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step2.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step2.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step2.DrawLatex(0.2, 0.7, "Neutrons in CND - step 2");

    //     myText->Print(fileName_Step2, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step2"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step2, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step2"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step2, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step2, "%s]", pdfFile_Step2);
    //     myCanvas->Print(fileName_Step2, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step2 plots - end */

    // #pragma region /* Saving Step3 plots - start */

    //     TLatex text_Step3;
    //     text_Step3.SetTextSize(0.05);

    //     string pdfFile_Step3_0 = ConfigOutPutName(PDFFile, "Step3").c_str();
    //     const char *pdfFile_Step3 = pdfFile_Step3_0.c_str();

    //     char fileName_Step3[100];
    //     sprintf(fileName_Step3, "%s[", pdfFile_Step3);
    //     myText->SaveAs(fileName_Step3);
    //     sprintf(fileName_Step3, "%s", pdfFile_Step3);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step3.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step3.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step3.DrawLatex(0.2, 0.7, "Neutrons in CND - step 3");

    //     myText->Print(fileName_Step3, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step3"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step3, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step3"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step3, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step3, "%s]", pdfFile_Step3);
    //     myCanvas->Print(fileName_Step3, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step3 plots - end */

    // #pragma region /* Saving Step4 plots - start */

    //     TLatex text_Step4;
    //     text_Step4.SetTextSize(0.05);

    //     string pdfFile_Step4_0 = ConfigOutPutName(PDFFile, "Step4").c_str();
    //     const char *pdfFile_Step4 = pdfFile_Step4_0.c_str();

    //     char fileName_Step4[100];
    //     sprintf(fileName_Step4, "%s[", pdfFile_Step4);
    //     myText->SaveAs(fileName_Step4);
    //     sprintf(fileName_Step4, "%s", pdfFile_Step4);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step4.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step4.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step4.DrawLatex(0.2, 0.7, "Neutrons in CND - step 4");

    //     myText->Print(fileName_Step4, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step4"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step4, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step4"))
    //         {
    //             myCanvas->cd(1);
    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step4, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step4, "%s]", pdfFile_Step4);
    //     myCanvas->Print(fileName_Step4, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step4 plots - end */

    // #pragma region /* Saving Step5 plots - start */

    //     TLatex text_Step5;
    //     text_Step5.SetTextSize(0.05);

    //     string pdfFile_Step5_0 = ConfigOutPutName(PDFFile, "Step5").c_str();
    //     const char *pdfFile_Step5 = pdfFile_Step5_0.c_str();

    //     char fileName_Step5[100];
    //     sprintf(fileName_Step5, "%s[", pdfFile_Step5);
    //     myText->SaveAs(fileName_Step5);
    //     sprintf(fileName_Step5, "%s", pdfFile_Step5);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step5.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step5.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step5.DrawLatex(0.2, 0.7, "Neutrons in CND - step 5");

    //     myText->Print(fileName_Step5, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step5"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step5, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step5"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step5, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step5, "%s]", pdfFile_Step5);
    //     myCanvas->Print(fileName_Step5, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step5 plots - end */

#pragma endregion /* Andrew's wrap up - end */

    // ======================================================================================================================================================================
    // Erin's wrap up
    // ======================================================================================================================================================================

#pragma region /* Erin's wrap up - start */

    // ======================================================================================================================================================================
    // Erin's wrap up
    // ======================================================================================================================================================================

    delete clasAna;

    f->cd();

    for (int i = 0; i < hist_list_1.size(); i++)
    {
        hist_list_1[i]->Write();
    }

    for (int i = 0; i < hist_list_2.size(); i++)
    {
        hist_list_2[i]->SetOption("colz");
        hist_list_2[i]->Write();
    }

    cout << "\n";
    cout << counter << " events counted!\n\n";

    // wrap it up
    outtxt.close();
    ntree->Write();
    f->Close();

#pragma endregion /* Erin's wrap up - end */

    // ======================================================================================================================================================================
    // Save log file
    // ======================================================================================================================================================================

#pragma region /* Save log file - start */

    // TODO: does't work - fix this!

    // Saving setup to log file ------------------------------------------------------------------------------------------------------------------------------------------

    // Saving setup to log file
    ofstream myLogFile;
    myLogFile.open(("./" + OutDir + "/Log_file.txt").c_str());

    myLogFile << "///////////////////////////////////////////////////////////////////////////\n";
    myLogFile << "// Input file was " << input_hipo << "\n";
    myLogFile << "// Beam energy was" << Ebeam << "\n";
    myLogFile << "///////////////////////////////////////////////////////////////////////////\n\n";

    myLogFile << "Run_Erins_features:\t" << Run_Erins_features << "\n";
    myLogFile << "Run_Andrews_work:\t" << Run_Andrews_work << "\n\n";

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

#pragma region /* Printouts - start */

    cout << "\033[33m\n\033[0m";
    cout << "\033[33mHIPO_FILES:\033[0m\t\t" << gSystem->Getenv("HIPO_FILES") << "\n";
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

#pragma endregion /* Printouts - end */

    return 0;

} // closes main function

#pragma endregion /* Erin main function - end */