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
#include "src/functions/GeneralFunctions.h"
#include "src/functions/neutron-veto/veto_functions.cpp"
#include "src/functions/Andrews_functions/Andrews_functions.cpp"
#include "src/functions/HipoChain_config.cpp"
#include "src/classes/clas12ana/clas12ana.cpp"

using namespace std;
using namespace clas12;

#pragma region /* Erin main function - start */

int D_getfeatures_Phase5(                                                                             //
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
    cout << "\033[33minput_hipo:\033[0m\t\t" << input_hipo << "\n";
    cout << "\033[33m\n\033[0m";
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
    cout << "\nClearing '" << OutDir << "'\n";
    system(("rm -r " + OutDir).c_str());
    cout << "\n";

    // Remake old output folder
    cout << "\nRemaking '" << OutDir << "'\n";
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

    const double mP = 0.93828;
    const double mN = 0.939;
    const double mD = 1.8756;
    double c = 29.9792; // cm/ns

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
    TH1D *h_nsize = new TH1D("nsize", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_nsize);

    // Reconstructed momentum (Erin)
    // ======================================================================================================================================================================
    TH2D *h_nangles = new TH2D("nangles", "Neutron Angular Distribution;#phi_{n};#theta_{n}", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_nangles);
    TH1D *h_pxminuspx = new TH1D("pxminuspx", "p_{x,n}-p_{x,pred};p_{x,n}-p_{x,pred} (GeV/c);Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pxminuspx);
    TH1D *h_pyminuspy = new TH1D("pyminuspy", "p_{y,n}-p_{y,pred};p_{y,n}-p_{y,pred} (GeV/c);Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pyminuspy);
    TH1D *h_pzminuspz = new TH1D("pzminuspz", "p_{z,n}-p_{z,pred};p_{z,n}-p_{z,pred} (GeV/c);Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pzminuspz);
    TH1D *h_pminusp = new TH1D("pminusp", "p_{n}-p_{pred};p_{n}-p_{pred} (GeV/c);Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pminusp);
    TH2D *h_pvsp = new TH2D("pvsp", "Momentum Resolution;p_{pred} (GeV/c);p_{measured} (GeV/c)", 100, 0, 1, 100, 0, 1);
    hist_list_2.push_back(h_pvsp);
    TH1D *h_cos0 = new TH1D("cos0", "Cosine of angle between generated and reconstructed p", 50, -1.1, 1.1);
    hist_list_1.push_back(h_cos0);
    TH2D *h_dpp = new TH2D("dpp", "Momentum Resolution;p_{generated} (GeV/c);#Delta p/p", 100, 0, 1, 100, -0.4, 0.4);
    hist_list_2.push_back(h_dpp);
    TH1D *h_energy = new TH1D("energy", "Neutron Energy Deposition;Energy (MeV);Counts", 100, 0, 25);
    hist_list_1.push_back(h_energy);
    TH2D *h_sec_phi = new TH2D("sec_phi", "Sector vs Phi of CND hits;phi (deg);Sector", 90, 0, 360, 25, 0, 25);
    hist_list_2.push_back(h_sec_phi);
    TH1D *h_mmiss = new TH1D("mmiss", "Missing Mass", 50, 0.5, 1.5);
    hist_list_1.push_back(h_mmiss);
    TH2D *h_mmiss_pn = new TH2D("mmiss_pn", "Missing Mass vs Measured Neutron Momentum;p_{n} (GeV/c);Missing Mass (GeV/c^{2})", 50, 0.25, 1., 50, 0.5, 1.5);
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
    TH2D *h_compare = new TH2D("compare", "(p_{pred}-p_{n})/p_{pred} vs #theta_{n,pred} (deg);(p_{pred}-p_{n})/p_{pred};#theta_{n,pred}", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_compare);
    TH2D *h_Edep_beta = new TH2D("Edep_beta", "Energy deposition vs #beta;#beta;E_{dep}", 50, 0, 1, 50, 0, 100);
    hist_list_2.push_back(h_Edep_beta);
    TH1D *h_p_all = new TH1D("p_all", "Momentum", 100, 0, 1.2);
    hist_list_1.push_back(h_p_all);

    TH2D *h_dpp_edep = new TH2D("dpp_edep", "#Delta p/p vs Energy Deposition;E_{dep} (MeVee);#Delta p/p", 50, 0, 40, 50, -0.4, 0.4);
    hist_list_2.push_back(h_dpp_edep);

    // good n / bad n set (Erin)
    // ======================================================================================================================================================================
    TH2D *h_nangles2 = new TH2D("nangles2", "Neutron Angles;phi;theta", 48, -180, 180, 45, 0, 180);
    hist_list_2.push_back(h_nangles2);
    TH1D *h_pxminuspx2 = new TH1D("pxminuspx2", "p_{x,n}-p_{x,pred};p_{x,n}-p_{x,pred};Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pxminuspx2);
    TH1D *h_pyminuspy2 = new TH1D("pyminuspy2", "p_{y,n}-p_{y,pred};p_{y,n}-p_{y,pred};Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pyminuspy2);
    TH1D *h_pzminuspz2 = new TH1D("pzminuspz2", "p_{z,n}-p_{z,pred};p_{z,n}-p_{z,pred};Counts", 100, -0.5, 0.5);
    hist_list_1.push_back(h_pzminuspz2);
    TH1D *h_pminusp2 = new TH1D("pminusp2", "p_{n}-p_{gen};p_{n}-p_{gen};Counts", 100, -0.5, 0.5);
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
    TH2D *h_compare2 = new TH2D("compare2", "(p_{pred}-p_{n})/p_{pred} vs #theta_{n,pred};(p_{pred}-p_{n})/p_{pred};#theta_{n,pred} (deg)", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_compare2);
    TH2D *h_Edep_beta2 = new TH2D("Edep_beta2", "Energy deposition vs #beta;#beta;E_{dep}", 50, 0, 1, 50, 0, 100);
    hist_list_2.push_back(h_Edep_beta2);
    TH1D *h_p_cut = new TH1D("p_cut", "Momentum", 100, 0, 1.2);
    hist_list_1.push_back(h_p_cut);

    TH2D *h_thetapn_dpp = new TH2D("thetapn_dpp", "#theta_{pn} vs #Delta p/p;(p_{pred}-p_{n})/p_{pred};#theta_{pn} (deg)", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_thetapn_dpp);
    TH2D *h_thetapn_dpp1 = new TH2D("thetapn_dpp1", "#theta_{pn} vs #Delta p/p;(p_{pred}-p_{n})/p_{pred};#theta_{pn} (deg)", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_thetapn_dpp1);
    TH2D *h_thetapn_dpp2 = new TH2D("thetapn_dpp2", "#theta_{pn} vs #Delta p/p;(p_{pred}-p_{n})/p_{pred};#theta_{pn} (deg)", 100, -2, 2, 90, 0, 180);
    hist_list_2.push_back(h_thetapn_dpp2);

    TH1D *h_anglediff = new TH1D("angle_diff", "Angle Diff between CVT hit and CND hit", 180, 0, 180);
    hist_list_1.push_back(h_anglediff);
    TH1D *h_anglediff2 = new TH1D("angle_diff2", "Angle Diff between CVT hit and CND hit", 180, 0, 180);
    hist_list_1.push_back(h_anglediff2);

    TH2D *h_ptheta_pred = new TH2D("ptheta_pred", "Predicted Momentum vs Angle for Final Signal Sample;#theta_{pred} (deg);p_{pred} (GeV/c)", 110, 35, 145, 100, 0.2, 1.3);
    hist_list_2.push_back(h_ptheta_pred);
    TH2D *h_ptheta = new TH2D("ptheta", "Measured Momentum vs Angle for Final Signal Sample;#theta_{n} (deg);p_{n} (GeV/c)", 110, 35, 145, 100, 0.2, 1.3);
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

    vector<TH1 *> hist_list_1_A;
    vector<TH1 *> hist_list_1_Step0_A;
    vector<TH1 *> hist_list_1_Step1_A;
    vector<TH1 *> hist_list_1_Step2_A;
    vector<TH1 *> hist_list_1_Step3_A;
    vector<TH1 *> hist_list_1_Step4_A;
    vector<TH1 *> hist_list_1_Step5_A;
    vector<TH2 *> hist_list_2_A;
    vector<TH2 *> hist_list_2_Step0_A;
    vector<TH2 *> hist_list_2_Step1_A;
    vector<TH2 *> hist_list_2_Step2_A;
    vector<TH2 *> hist_list_2_Step3_A;
    vector<TH2 *> hist_list_2_Step4_A;
    vector<TH2 *> hist_list_2_Step5_A;

    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);

    gStyle->SetTitleXOffset(0.8);
    gStyle->SetTitleYOffset(0.8);

    char temp_name_A[100];
    char temp_title_A[100];

    // Checks on which events have neutrons (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_xB_mmiss_epFD = new TH2D("xB_mmiss_epFD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epFD);
    TH2D *h_xB_mmiss_epnFD = new TH2D("xB_mmiss_epnFD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epnFD);
    TH2D *h_xB_mmiss_epn_goodN_pFD = new TH2D("xB_mmiss_epngoodFD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epn_goodN_pFD);
    TH2D *h_xB_mmiss_epCD = new TH2D("xB_mmiss_epCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epCD);
    TH2D *h_xB_mmiss_epnCD = new TH2D("xB_mmiss_epnCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epnCD);
    TH2D *h_xB_mmiss_epn_goodN_pCD = new TH2D("xB_mmiss_epngoodCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epn_goodN_pCD);

    TH1D *h_pmiss_ep = new TH1D("pmiss_ep", "p_{miss} ep;p_{miss} [GeV/c];Counts", 25, 0.25, 1.0);
    hist_list_1_A.push_back(h_pmiss_ep);

    // Step Zero (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step0 = new TH2D("pnRes_theta_nmiss_Step0", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step0);
    hist_list_2_Step0_A.push_back(h_pnRes_theta_nmiss_Step0);

    TH1D *h_ToF_goodN_Step0 = new TH1D("ToF_goodN_Step0", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step0);
    hist_list_1_Step0_A.push_back(h_ToF_goodN_Step0);
    TH1D *h_beta_goodN_Step0 = new TH1D("beta_goodN_Step0", "#beta of CND Neutrons;#beta;Counts", 50, 0, 1.1);
    hist_list_1_A.push_back(h_beta_goodN_Step0);
    hist_list_1_Step0_A.push_back(h_beta_goodN_Step0);
    TH1D *h_Edep_goodN_Step0 = new TH1D("Edep_goodN_Step0", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_goodN_Step0);
    hist_list_1_Step0_A.push_back(h_Edep_goodN_Step0);
    TH2D *h_beta_Edep_goodN_Step0 = new TH2D("Edep_beta_goodN", "#beta vs. E_{dep} of CND Neutrons;#beta;E_{dep} [MeV]", 50, 0, 1.1, 50, 0, 100);
    hist_list_2_A.push_back(h_beta_Edep_goodN_Step0);
    hist_list_2_Step0_A.push_back(h_beta_Edep_goodN_Step0);
    TH2D *h_Edep_ToF_goodN_Step0 = new TH2D("Edep_ToF_goodN_Step0", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_goodN_Step0);
    hist_list_2_Step0_A.push_back(h_Edep_ToF_goodN_Step0);

    TH1D *h_ToF_badN_Step0 = new TH1D("ToF_badN_Step0", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step0);
    hist_list_1_Step0_A.push_back(h_ToF_badN_Step0);
    TH1D *h_beta_badN_Step0 = new TH1D("beta_badN_Step0", "#beta of CND Neutrons;#beta;Counts", 50, 0, 1.1);
    hist_list_1_A.push_back(h_beta_badN_Step0);
    hist_list_1_Step0_A.push_back(h_beta_badN_Step0);
    TH1D *h_Edep_badN_Step0 = new TH1D("Edep_badN_Step0", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_badN_Step0);
    hist_list_1_Step0_A.push_back(h_Edep_badN_Step0);
    TH2D *h_beta_Edep_badN_Step0 = new TH2D("Edep_beta_badN", "#beta vs. E_{dep} of CND Neutrons;#beta;E_{dep} [MeV]", 50, 0, 1.1, 50, 0, 100);
    hist_list_2_A.push_back(h_beta_Edep_badN_Step0);
    hist_list_2_Step0_A.push_back(h_beta_Edep_badN_Step0);
    TH2D *h_Edep_ToF_badN_Step0 = new TH2D("Edep_ToF_badN_Step0", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_badN_Step0);
    hist_list_2_Step0_A.push_back(h_Edep_ToF_badN_Step0);

    // Step One (After Beta Cut) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step1 = new TH2D("pnRes_theta_nmiss_Step1", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step1);
    hist_list_2_Step1_A.push_back(h_pnRes_theta_nmiss_Step1);

    TH1D *h_pmiss_goodN_Step1 = new TH1D("pmiss_goodN_Step1", "p_{miss} Step1;p_{miss} [GeV/c];Counts", 25, 0.25, 1.0);
    hist_list_1_A.push_back(h_pmiss_goodN_Step1);
    hist_list_1_Step1_A.push_back(h_pmiss_goodN_Step1);

    TH1D *h_ToF_goodN_Step1 = new TH1D("ToF_goodN_Step1", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step1);
    hist_list_1_Step1_A.push_back(h_ToF_goodN_Step1);
    TH1D *h_edep_goodN_Step1 = new TH1D("edep_goodN_Step1", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 100, 0, 50);
    hist_list_1_A.push_back(h_edep_goodN_Step1);
    hist_list_1_Step1_A.push_back(h_edep_goodN_Step1);
    TH2D *h_Edep_ToF_goodN_Step1 = new TH2D("Edep_ToF_goodN_Step1", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_goodN_Step1);
    hist_list_2_Step1_A.push_back(h_Edep_ToF_goodN_Step1);

    TH1D *h_ToF_badN_Step1 = new TH1D("ToF_badN_Step1", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step1);
    hist_list_1_Step1_A.push_back(h_ToF_badN_Step1);
    TH1D *h_edep_badN_Step1 = new TH1D("edep_badN_Step1", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 100, 0, 50);
    hist_list_1_A.push_back(h_edep_badN_Step1);
    hist_list_1_Step1_A.push_back(h_edep_badN_Step1);
    TH2D *h_Edep_ToF_badN_Step1 = new TH2D("Edep_ToF_badN_Step1", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_badN_Step1);
    hist_list_2_Step1_A.push_back(h_Edep_ToF_badN_Step1);

    TH1D *h_edep_over_edepCTOT_goodN_Step1 = new TH1D("edep_over_edepCTOT_goodN_Step1", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 100, 0, 5);
    hist_list_1_A.push_back(h_edep_over_edepCTOT_goodN_Step1);
    hist_list_1_Step1_A.push_back(h_edep_over_edepCTOT_goodN_Step1);
    TH1D *h_edep_over_edepCTOT_badN_Step1 = new TH1D("edep_over_edepCTOT_badN_Step1", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 100, 0, 5);
    hist_list_1_A.push_back(h_edep_over_edepCTOT_badN_Step1);
    hist_list_1_Step1_A.push_back(h_edep_over_edepCTOT_badN_Step1);

    TH1D *h_edep_goodN_withNearbyPos_Step1 = new TH1D("edep_goodN_withNearbyPos_Step1", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 100, 0, 50);
    hist_list_1_A.push_back(h_edep_goodN_withNearbyPos_Step1);
    hist_list_1_Step1_A.push_back(h_edep_goodN_withNearbyPos_Step1);
    TH1D *h_edep_badN_withNearbyPos_Step1 = new TH1D("edep_badN_withNearbyPos_Step1", "E_{dep} of CND Neutrons;E_{dep} [MeF];Counts", 100, 0, 50);
    hist_list_1_A.push_back(h_edep_badN_withNearbyPos_Step1);
    hist_list_1_Step1_A.push_back(h_edep_badN_withNearbyPos_Step1);

    TH1D *h_sdiff_pos_goodN_Step1_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_goodN_Step1_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step1_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        hist_list_1_A.push_back(h_sdiff_pos_goodN_Step1_layer[k]);
        hist_list_1_Step1_A.push_back(h_sdiff_pos_goodN_Step1_layer[k]);
    }

    TH1D *h_sdiff_pos_badN_Step1_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_badN_Step1_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step1_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        hist_list_1_A.push_back(h_sdiff_pos_badN_Step1_layer[k]);
        hist_list_1_Step1_A.push_back(h_sdiff_pos_badN_Step1_layer[k]);
    }

    TH2D *h_sdiff_pos_mom_goodN_Step1_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_mom_goodN_Step1_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_goodN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        hist_list_2_A.push_back(h_sdiff_pos_mom_goodN_Step1_layer[k]);
        hist_list_2_Step1_A.push_back(h_sdiff_pos_mom_goodN_Step1_layer[k]);
    }

    TH2D *h_sdiff_pos_mom_badN_Step1_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_mom_badN_Step1_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)", k - 3);
        h_sdiff_pos_mom_badN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 4);
        hist_list_2_A.push_back(h_sdiff_pos_mom_badN_Step1_layer[k]);
        hist_list_2_Step1_A.push_back(h_sdiff_pos_mom_badN_Step1_layer[k]);
    }

    TH2D *h_sdiff_pos_z_goodN_Step1_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_z_goodN_Step1_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_goodN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        hist_list_2_A.push_back(h_sdiff_pos_z_goodN_Step1_layer[k]);
        hist_list_2_Step1_A.push_back(h_sdiff_pos_z_goodN_Step1_layer[k]);
    }

    TH2D *h_sdiff_pos_z_badN_Step1_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_z_badN_Step1_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_z_badN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, -40.0, 40.0);
        hist_list_2_A.push_back(h_sdiff_pos_z_badN_Step1_layer[k]);
        hist_list_2_Step1_A.push_back(h_sdiff_pos_z_badN_Step1_layer[k]);
    }

    TH2D *h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_goodN_Step1_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
        hist_list_2_A.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[k]);
        hist_list_2_Step1_A.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[k]);
    }

    TH2D *h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_diff_ToFc_z_badN_Step1_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)", k - 3);
        h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[k] = new TH2D(temp_name_A, temp_title_A, 24, -11.5, 12.5, 50, 0, 300);
        hist_list_2_A.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[k]);
        hist_list_2_Step1_A.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[k]);
    }

    TH2D *h_diff_ToFc_z_Edep_noNear_goodN_Step1 = new TH2D("diff_ToFc_z_Edep_noNear_goodN_Step1", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_noNear_goodN_Step1);
    hist_list_2_Step1_A.push_back(h_diff_ToFc_z_Edep_noNear_goodN_Step1);
    TH2D *h_diff_ToFc_z_Edep_yesNear_goodN_Step1 = new TH2D("diff_ToFc_z_Edep_yesNear_goodN_Step1", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_yesNear_goodN_Step1);
    hist_list_2_Step1_A.push_back(h_diff_ToFc_z_Edep_yesNear_goodN_Step1);
    TH2D *h_diff_ToFc_z_Edep_noNear_badN_Step1 = new TH2D("diff_ToFc_z_Edep_noNear_badN_Step1", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_noNear_badN_Step1);
    hist_list_2_Step1_A.push_back(h_diff_ToFc_z_Edep_noNear_badN_Step1);
    TH2D *h_diff_ToFc_z_Edep_yesNear_badN_Step1 = new TH2D("diff_ToFc_z_Edep_yesNear_badN_Step1", "ToF*c - z vs. E_{dep} of CND Neutrons with no Nearby Tracks;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_yesNear_badN_Step1);
    hist_list_2_Step1_A.push_back(h_diff_ToFc_z_Edep_yesNear_badN_Step1);

    // Step Two (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step2 = new TH2D("pnRes_theta_nmiss_Step2", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step2);
    hist_list_2_Step2_A.push_back(h_pnRes_theta_nmiss_Step2);

    TH1D *h_ToF_goodN_Step2 = new TH1D("ToF_goodN_Step2", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step2);
    hist_list_1_Step2_A.push_back(h_ToF_goodN_Step2);
    TH2D *h_Edep_ToF_goodN_Step2 = new TH2D("Edep_ToF_goodN_Step2", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_goodN_Step2);
    hist_list_2_Step2_A.push_back(h_Edep_ToF_goodN_Step2);

    TH1D *h_ToF_badN_Step2 = new TH1D("ToF_badN_Step2", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step2);
    hist_list_1_Step2_A.push_back(h_ToF_badN_Step2);
    TH2D *h_Edep_ToF_badN_Step2 = new TH2D("Edep_ToF_badN_Step2", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_badN_Step2);
    hist_list_2_Step2_A.push_back(h_Edep_ToF_badN_Step2);

    TH1D *h_sdiff_pos_goodN_Step2_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_goodN_Step2_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_goodN_Step2_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        hist_list_1_A.push_back(h_sdiff_pos_goodN_Step2_layer[k]);
        hist_list_1_Step2_A.push_back(h_sdiff_pos_goodN_Step2_layer[k]);
    }

    TH1D *h_sdiff_pos_badN_Step2_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_pos_badN_Step2_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)", k - 3);
        h_sdiff_pos_badN_Step2_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        hist_list_1_A.push_back(h_sdiff_pos_badN_Step2_layer[k]);
        hist_list_1_Step2_A.push_back(h_sdiff_pos_badN_Step2_layer[k]);
    }

    TH1D *h_sdiff_allhit_goodN_Step2_layer[7];

    for (int k = 0; k < 7; k++)
    {
        sprintf(temp_name_A, "sdiff_allhit_goodN_Step2_layer_%d", k - 3);
        sprintf(temp_title_A, "Nuetral Sector minus Random Hit Sector (Layer Difference = %d)", k - 3);
        h_sdiff_allhit_goodN_Step2_layer[k] = new TH1D(temp_name_A, temp_title_A, 24, -11.5, 12.5);
        hist_list_1_A.push_back(h_sdiff_allhit_goodN_Step2_layer[k]);
        hist_list_1_Step2_A.push_back(h_sdiff_allhit_goodN_Step2_layer[k]);
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
        hist_list_1_Step2_A.push_back(h_sdiff_allhit_badN_Step2_layer[k]);
    }

    TH2D *h_sdiff_ldiff_allhit_badN_Step2 = new TH2D("sdiff_ldiff_allhit_badN_Step2", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    hist_list_2_A.push_back(h_sdiff_ldiff_allhit_badN_Step2);
    hist_list_2_Step2_A.push_back(h_sdiff_ldiff_allhit_badN_Step2);

    TH1D *h_numberNearby_goodN_Step2 = new TH1D("numberNearby_goodN_Step2", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    hist_list_1_A.push_back(h_numberNearby_goodN_Step2);
    hist_list_1_Step2_A.push_back(h_numberNearby_goodN_Step2);
    TH2D *h_numberNearby_momN_goodN_Step2 = new TH2D("numberNearby_momN_goodN_Step2", "Number of Nearby Hits vs. p_{n} for CND Neutrons;# Hits;p_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    hist_list_2_A.push_back(h_numberNearby_momN_goodN_Step2);
    hist_list_2_Step2_A.push_back(h_numberNearby_momN_goodN_Step2);
    TH1D *h_numberNearby_badN_Step2 = new TH1D("numberNearby_badN_Step2", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    hist_list_1_A.push_back(h_numberNearby_badN_Step2);
    hist_list_1_Step2_A.push_back(h_numberNearby_badN_Step2);
    TH2D *h_numberNearby_momN_badN_Step2 = new TH2D("numberNearby_momN_badN_Step2", "Number of Nearby Hits vs. p_{n} for CND Neutrons;# Hits;p_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    hist_list_2_A.push_back(h_numberNearby_momN_badN_Step2);
    hist_list_2_Step2_A.push_back(h_numberNearby_momN_badN_Step2);

    TH1D *h_NearbyEdep_goodN_Step2 = new TH1D("NearbyEdep_goodN_Step2", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_NearbyEdep_goodN_Step2);
    hist_list_1_Step2_A.push_back(h_NearbyEdep_goodN_Step2);
    TH1D *h_NearbyEdep_badN_Step2 = new TH1D("NearbyEdep_badN_Step2", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_NearbyEdep_badN_Step2);
    hist_list_1_Step2_A.push_back(h_NearbyEdep_badN_Step2);

    TH1D *h_nsector_goodN_Step2 = new TH1D("nsector_goodN_Step2", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    hist_list_1_A.push_back(h_nsector_goodN_Step2);
    TH1D *h_nsector_badN_Step2 = new TH1D("nsector_badN_Step2", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    hist_list_1_A.push_back(h_nsector_badN_Step2);
    hist_list_1_Step2_A.push_back(h_nsector_badN_Step2);

    // Step Three (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step3 = new TH2D("pnRes_theta_nmiss_Step3", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step3);
    hist_list_2_Step3_A.push_back(h_pnRes_theta_nmiss_Step3);

    TH1D *h_ToF_goodN_Step3 = new TH1D("ToF_goodN_Step3", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step3);
    hist_list_1_Step3_A.push_back(h_ToF_goodN_Step3);
    TH2D *h_Edep_ToF_goodN_Step3 = new TH2D("Edep_ToF_goodN_Step3", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_goodN_Step3);
    hist_list_2_Step3_A.push_back(h_Edep_ToF_goodN_Step3);

    TH1D *h_ToF_badN_Step3 = new TH1D("ToF_badN_Step3", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step3);
    hist_list_1_Step3_A.push_back(h_ToF_badN_Step3);
    TH2D *h_Edep_ToF_badN_Step3 = new TH2D("Edep_ToF_badN_Step3", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_badN_Step3);
    hist_list_2_Step3_A.push_back(h_Edep_ToF_badN_Step3);

    TH2D *h_sdiff_ldiff_allhit_goodN_Step3 = new TH2D("sdiff_ldiff_allhit_goodN_Step3", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    hist_list_2_A.push_back(h_sdiff_ldiff_allhit_goodN_Step3);
    hist_list_2_Step3_A.push_back(h_sdiff_ldiff_allhit_goodN_Step3);
    TH2D *h_sdiff_ldiff_allhit_badN_Step3 = new TH2D("sdiff_ldiff_allhit_badN_Step3", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 7, -3.5, 3.5);
    hist_list_2_A.push_back(h_sdiff_ldiff_allhit_badN_Step3);
    hist_list_2_Step3_A.push_back(h_sdiff_ldiff_allhit_badN_Step3);

    TH2D *h_sdiff_ldiff_CTOFhit_goodN_Step3 = new TH2D("sdiff_ldiff_CTOFhit_goodN_Step3", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 3, 0.5, 3.5);
    hist_list_2_A.push_back(h_sdiff_ldiff_CTOFhit_goodN_Step3);
    hist_list_2_Step3_A.push_back(h_sdiff_ldiff_CTOFhit_goodN_Step3);
    TH2D *h_sdiff_ldiff_CTOFhit_badN_Step3 = new TH2D("sdiff_ldiff_CTOFhit_badN_Step3", "Sector Difference vs. Layer Difference;Sector Difference;Layer Difference", 24, -11.5, 12.5, 3, 0.5, 3.5);
    hist_list_2_A.push_back(h_sdiff_ldiff_CTOFhit_badN_Step3);
    hist_list_2_Step3_A.push_back(h_sdiff_ldiff_CTOFhit_badN_Step3);

    TH1D *h_numberCTOF_goodN_Step3 = new TH1D("numberCTOF_goodN_Step3", "Number of CTOF Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    hist_list_1_A.push_back(h_numberCTOF_goodN_Step3);
    hist_list_1_Step3_A.push_back(h_numberCTOF_goodN_Step3);
    TH2D *h_numberCTOF_momN_goodN_Step3 = new TH2D("numberCTOF_momN_goodN_Step3", "Number of CTOF Hits vs. p_{n} for CND Neutrons;# Hits;p_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    hist_list_2_A.push_back(h_numberCTOF_momN_goodN_Step3);
    hist_list_2_Step3_A.push_back(h_numberCTOF_momN_goodN_Step3);
    TH1D *h_numberCTOF_badN_Step3 = new TH1D("numberCTOF_badN_Step3", "Number of CTOF Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    hist_list_1_A.push_back(h_numberCTOF_badN_Step3);
    hist_list_1_Step3_A.push_back(h_numberCTOF_badN_Step3);
    TH2D *h_numberCTOF_momN_badN_Step3 = new TH2D("numberCTOF_momN_badN_Step3", "Number of CTOF Hits vs. p_{n} for CND Neutrons;# Hits;p_{n} [GeV/c]", 9, -0.5, 8.5, 50, 0, 1.3);
    hist_list_2_A.push_back(h_numberCTOF_momN_badN_Step3);
    hist_list_2_Step3_A.push_back(h_numberCTOF_momN_badN_Step3);

    // Step Four (After applying Phi Diff CND hit cut) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step4 = new TH2D("pnRes_theta_nmiss_Step4", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step4);
    hist_list_2_Step4_A.push_back(h_pnRes_theta_nmiss_Step4);

    TH1D *h_ToF_goodN_Step4 = new TH1D("ToF_goodN_Step4", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step4);
    hist_list_1_Step4_A.push_back(h_ToF_goodN_Step4);
    TH2D *h_Edep_ToF_goodN_Step4 = new TH2D("Edep_ToF_goodN_Step4", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_goodN_Step4);
    hist_list_2_Step4_A.push_back(h_Edep_ToF_goodN_Step4);

    TH1D *h_ToF_badN_Step4 = new TH1D("ToF_badN_Step4", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step4);
    hist_list_1_Step4_A.push_back(h_ToF_badN_Step4);
    TH2D *h_Edep_ToF_badN_Step4 = new TH2D("Edep_ToF_badN_Step4", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_badN_Step4);
    hist_list_2_Step4_A.push_back(h_Edep_ToF_badN_Step4);

    // Step Five (After event selection cuts) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step5 = new TH2D("pnRes_theta_nmiss_Step5", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss} [#circ]", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step5);
    hist_list_2_Step5_A.push_back(h_pnRes_theta_nmiss_Step5);

    TH1D *h_ToF_goodN_Step5 = new TH1D("ToF_goodN_Step5", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step5);
    hist_list_1_Step5_A.push_back(h_ToF_goodN_Step5);
    TH2D *h_Edep_ToF_goodN_Step5 = new TH2D("Edep_ToF_goodN_Step5", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_goodN_Step5);
    hist_list_2_Step5_A.push_back(h_Edep_ToF_goodN_Step5);

    TH1D *h_ToF_badN_Step5 = new TH1D("ToF_badN_Step5", "ToF of CND Neutrons;ToF [ns];Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step5);
    hist_list_1_Step5_A.push_back(h_ToF_badN_Step5);
    TH2D *h_Edep_ToF_badN_Step5 = new TH2D("Edep_ToF_badN_Step5", "E_{dep} vs. ToF of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 100, 0, 50);
    hist_list_2_A.push_back(h_Edep_ToF_badN_Step5);
    hist_list_2_Step5_A.push_back(h_Edep_ToF_badN_Step5);

    TH1D *h_pmiss_goodN_Step5 = new TH1D("pmiss_goodN_Step5", "p_{miss} good N Step5;p_{miss} [GeV/c];Counts", 25, 0.25, 1.0);
    hist_list_1_A.push_back(h_pmiss_goodN_Step5);
    hist_list_1_Step5_A.push_back(h_pmiss_goodN_Step5);
    TH1D *h_pmiss_allN_Step5 = new TH1D("pmiss_allN_Step5", "p_{miss} all N Step5;p_{miss} [GeV/c];Counts", 25, 0.25, 1.0);
    hist_list_1_A.push_back(h_pmiss_allN_Step5);
    hist_list_1_Step5_A.push_back(h_pmiss_allN_Step5);

    TH1D *h_Edep_infront_goodN_Step5 = new TH1D("Edep_infront_goodN_Step5", "E_{dep} of Hit infront CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_infront_goodN_Step5);
    hist_list_1_Step5_A.push_back(h_Edep_infront_goodN_Step5);
    TH1D *h_Edep_behind_goodN_Step5 = new TH1D("Edep_behind_goodN_Step5", "E_{dep} of Hit behind CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_behind_goodN_Step5);
    hist_list_1_Step5_A.push_back(h_Edep_behind_goodN_Step5);

    TH1D *h_Edep_infront_badN_Step5 = new TH1D("Edep_infront_badN_Step5", "E_{dep} of Hit infront CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_infront_badN_Step5);
    hist_list_1_Step5_A.push_back(h_Edep_infront_badN_Step5);
    TH1D *h_Edep_behind_badN_Step5 = new TH1D("Edep_behind_badN_Step5", "E_{dep} of Hit behind CND Neutrons;E_{dep} [MeV];Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_behind_badN_Step5);
    hist_list_1_Step5_A.push_back(h_Edep_behind_badN_Step5);

    TH2D *h_diff_ToFc_z_Edep_goodN_Step5 = new TH2D("diff_ToFc_z_Edep_goodN_Step5", "ToF*c - z vs. E_{dep} of CND Neutrons;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_goodN_Step5);
    hist_list_2_Step5_A.push_back(h_diff_ToFc_z_Edep_goodN_Step5);
    TH2D *h_diff_ToFc_z_Edep_badN_Step5 = new TH2D("diff_ToFc_z_Edep_badN_Step5", "ToF*c - z vs. E_{dep} of CND Neutrons;ToF*c-z [cm];E_{dep} [MeV]", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_badN_Step5);
    hist_list_2_Step5_A.push_back(h_diff_ToFc_z_Edep_badN_Step5);

    TH2D *h_diff_ToFc_z_Edep_goodN_Step5_layer[3];

    for (int k = 0; k < 3; k++)
    {
        sprintf(temp_name_A, "diff_ToFc_z_goodN_Step5_layer_%d", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step5_layer[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        hist_list_2_A.push_back(h_diff_ToFc_z_Edep_goodN_Step5_layer[k]);
        hist_list_2_Step5_A.push_back(h_diff_ToFc_z_Edep_goodN_Step5_layer[k]);
    }

    TH2D *h_diff_ToFc_z_Edep_badN_Step5_layer[3];

    for (int k = 0; k < 3; k++)
    {
        sprintf(temp_name_A, "diff_ToFc_z_badN_Step5_layer_%d", k + 1);
        sprintf(temp_title_A, "ToF*c - z vs. E_{dep} of CND Neutrons (Layer Difference = %d);ToF*c-z [cm];E_{dep} [MeV]", k + 1);
        h_diff_ToFc_z_Edep_badN_Step5_layer[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        hist_list_2_A.push_back(h_diff_ToFc_z_Edep_badN_Step5_layer[k]);
        hist_list_2_Step5_A.push_back(h_diff_ToFc_z_Edep_badN_Step5_layer[k]);
    }

    TH1D *h_phidiff_en_goodN_Step5 = new TH1D("phidiff_en_goodN_Step5", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
    hist_list_1_A.push_back(h_phidiff_en_goodN_Step5);
    hist_list_1_Step5_A.push_back(h_phidiff_en_goodN_Step5);

    TH1D *h_phidiff_en_badN_Step5 = new TH1D("phidiff_en_badN_Step5", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}| [#circ];Counts", 90, 0, 180);
    hist_list_1_A.push_back(h_phidiff_en_badN_Step5);
    hist_list_1_Step5_A.push_back(h_phidiff_en_badN_Step5);

    TH1D *h_TP_goodN_Step5 = new TH1D("TP_goodN_Step5", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    hist_list_1_A.push_back(h_TP_goodN_Step5);
    hist_list_1_Step5_A.push_back(h_TP_goodN_Step5);
    TH1D *h_TP_badN_Step5 = new TH1D("TP_badN_Step5", "ToF/path of CND Neutrons;ToF/path [ns/m];Counts", 150, 0, 50);
    hist_list_1_A.push_back(h_TP_badN_Step5);
    hist_list_1_Step5_A.push_back(h_TP_badN_Step5);

    TH1D *h_Z_goodN_Step5 = new TH1D("Z_goodN_Step5", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    hist_list_1_A.push_back(h_Z_goodN_Step5);
    hist_list_1_Step5_A.push_back(h_Z_goodN_Step5);
    TH1D *h_Z_badN_Step5 = new TH1D("Z_badN_Step5", "Z of CND Neutrons;Z [cm];Counts", 100, -60, 60);
    hist_list_1_A.push_back(h_Z_badN_Step5);
    hist_list_1_Step5_A.push_back(h_Z_badN_Step5);

    TH2D *h_beta_Edep_goodN_Step5 = new TH2D("Edep_beta_goodN", "#beta vs. E_{dep} of CND Neutrons;#beta;E_{dep} [MeV]", 50, 0, 1.1, 50, 0, 100);
    hist_list_2_A.push_back(h_beta_Edep_goodN_Step5);
    hist_list_2_Step5_A.push_back(h_beta_Edep_goodN_Step5);
    TH2D *h_beta_Edep_badN_Step5 = new TH2D("Edep_beta_badN", "#beta vs. E_{dep} of CND Neutrons;#beta;E_{dep} [MeV]", 50, 0, 1.1, 50, 0, 100);
    hist_list_2_A.push_back(h_beta_Edep_badN_Step5);
    hist_list_2_Step5_A.push_back(h_beta_Edep_badN_Step5);

    TH2D *h_ToF_Edep_goodN_Step5 = new TH2D("ToF_Edep_goodN_Step5", "ToF vs. E_{dep} of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
    hist_list_2_A.push_back(h_ToF_Edep_goodN_Step5);
    hist_list_2_Step5_A.push_back(h_ToF_Edep_goodN_Step5);
    TH2D *h_ToF_Edep_badN_Step5 = new TH2D("ToF_Edep_badN_Step5", "ToF vs. E_{dep} of CND Neutrons;ToF [ns];E_{dep} [MeV]", 100, 0, 20, 50, 0, 100);
    hist_list_2_A.push_back(h_ToF_Edep_badN_Step5);
    hist_list_2_Step5_A.push_back(h_ToF_Edep_badN_Step5);

    TH2D *h_TP_Edep_goodN_Step5 = new TH2D("TP_Edep_goodN_Step5", "TP vs. E_{dep} of CND Neutrons;TP [ns/m];E_{dep} [MeV]", 150, 0, 50, 50, 0, 100);
    hist_list_2_A.push_back(h_TP_Edep_goodN_Step5);
    hist_list_2_Step5_A.push_back(h_TP_Edep_goodN_Step5);
    TH2D *h_TP_Edep_badN_Step5 = new TH2D("TP_Edep_badN_Step5", "TP vs. E_{dep} of CND Neutrons;TP [ns/m];E_{dep} [MeV]", 150, 0, 50, 50, 0, 100);
    hist_list_2_A.push_back(h_TP_Edep_badN_Step5);
    hist_list_2_Step5_A.push_back(h_TP_Edep_badN_Step5);

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        hist_list_1_A[i]->Sumw2();
        hist_list_1_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_A[i]->GetYaxis()->CenterTitle();
        hist_list_1_Step0_A[i]->Sumw2();
        hist_list_1_Step0_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_Step0_A[i]->GetYaxis()->CenterTitle();
        hist_list_1_Step1_A[i]->Sumw2();
        hist_list_1_Step1_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_Step1_A[i]->GetYaxis()->CenterTitle();
        hist_list_1_Step3_A[i]->Sumw2();
        hist_list_1_Step3_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_Step3_A[i]->GetYaxis()->CenterTitle();
        hist_list_1_Step3_A[i]->Sumw2();
        hist_list_1_Step3_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_Step3_A[i]->GetYaxis()->CenterTitle();
        hist_list_1_Step4_A[i]->Sumw2();
        hist_list_1_Step4_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_Step4_A[i]->GetYaxis()->CenterTitle();
        hist_list_1_Step5_A[i]->Sumw2();
        hist_list_1_Step5_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_Step5_A[i]->GetYaxis()->CenterTitle();
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        hist_list_2_A[i]->Sumw2();
        hist_list_2_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_A[i]->GetYaxis()->CenterTitle();
        hist_list_2_Step0_A[i]->Sumw2();
        hist_list_2_Step0_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_Step0_A[i]->GetYaxis()->CenterTitle();
        hist_list_2_Step1_A[i]->Sumw2();
        hist_list_2_Step1_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_Step1_A[i]->GetYaxis()->CenterTitle();
        hist_list_2_Step3_A[i]->Sumw2();
        hist_list_2_Step3_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_Step3_A[i]->GetYaxis()->CenterTitle();
        hist_list_2_Step3_A[i]->Sumw2();
        hist_list_2_Step3_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_Step3_A[i]->GetYaxis()->CenterTitle();
        hist_list_2_Step4_A[i]->Sumw2();
        hist_list_2_Step4_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_Step4_A[i]->GetYaxis()->CenterTitle();
        hist_list_2_Step5_A[i]->Sumw2();
        hist_list_2_Step5_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_Step5_A[i]->GetYaxis()->CenterTitle();
    }

#pragma endregion /* Andrew's histograms - end */

    // ======================================================================================================================================================================
    // Chain loop
    // ======================================================================================================================================================================

#pragma region /* Chain loop - start */

    int counter_A = 0; /* From Andrew */

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
                h_vzp_fd->Fill(Vz_p - Vz_e);

                if (fabs(Vz_p - Vz_e) > 5)
                // if (abs(vzp - vze) > 5) // Erin's original
                {
                    continue;
                }

                h_chipid_fd->Fill(chipid);
                h_dbeta_p_fd->Fill(P_p.Mag(), dbeta);

                if (P_p.Mag() < 0.5)
                {
                    continue;
                }

                if (P_p.Mag() > 3.0)
                {
                    continue;
                }

                if (fabs(dbeta) > 0.03)
                // if (abs(dbeta) > 0.03) // Erin's original
                {
                    continue;
                }
            }
            else if (Protons[i]->getRegion() == CD)
            {
                h_vzp_cd->Fill(Vz_p - Vz_e);

                if (fabs(Vz_p - Vz_e) > 4)
                // if (abs(vzp - vze) > 4) // Erin's original
                {
                    continue;
                }

                h_chipid_cd->Fill(chipid);
                // if (abs(chipid)>4) {continue;}
                h_dbeta_p_cd->Fill(P_p.Mag(), dbeta);

                if (P_p.Mag() < 0.3)
                {
                    continue;
                }

                if (P_p.Mag() > 1.5)
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

        P_p.SetMagThetaPhi(Protons[p_index]->getP(), Protons[p_index]->getTheta(), Protons[p_index]->getPhi());

        // Determin where is the proton. Moved from angle cuts to getRegion() by the advice of Andrew.
        bool pInFD = (Protons[p_index]->getRegion() == FD); // My addition
        bool pInCD = (Protons[p_index]->getRegion() == CD); // My addition

        // Proton goes to the CD only.
        // TODO: include FD protons in the future?
        if (!pInCD)
        {
            continue;
        }

        /*
        // todo: use getdregion - much better!
        bool pInFD = (P_p.Theta() * 180. / M_PI < 40);                                      // My addition
        bool pInCD = (P_p.Theta() * 180. / M_PI >= 40 && P_p.Theta() * 180. / M_PI <= 140); // My addition

        if (P_p.Theta() * 180. / M_PI < 40 || P_p.Theta() * 180. / M_PI > 140) // p goes to CD
        {
            continue;
        }
        */

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

#pragma region /* Erin's features */

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
            h_nsize->Fill(Neutrons.size());

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

#pragma endregion /* Erin's features */

        // ==================================================================================================================================================================
        // Andrew's manual work
        // ==================================================================================================================================================================

#pragma region /* Andrew's manual work */

        if (Run_Andrews_work)
        {

            /*
            auto electrons = clasAna->getByPid(11); // From Erin's code
            auto protons = clasAna->getByPid(2212); // From Erin's code

            auto AllParticles = c12->getDetParticles();

            // Andrew's original - commented out!
            // auto AllParticles = c12->getDetParticles();
            // auto electrons = c12->getByID(11);

            double weight = 1;

            if (isMC)
            {
                weight = c12->mcevent()->getWeight();
            }

            TVector3 p_b(0, 0, Ebeam);
             */

#pragma region /* Electrons - start */

            /*
            if (electrons.size() != 1)
            {
                continue;
            }

            TVector3 p_e;
            p_e.SetMagThetaPhi(electrons[0]->getP(), electrons[0]->getTheta(), electrons[0]->getPhi());

            double EoP_e = (electrons[0]->cal(PCAL)->getEnergy() + electrons[0]->cal(ECIN)->getEnergy() + electrons[0]->cal(ECOUT)->getEnergy()) / p_e.Mag();
            int nphe = electrons[0]->che(HTCC)->getNphe();
            double vtz_e = electrons[0]->par()->getVz();

            // Andrew's original - commented out!
            // if (!myCut.electroncut(c12))
            // {
            //     continue;
            // }

            int esector = electrons[0]->getSector();

            /////////////////////////////////////
            // Electron Kinematics
            /////////////////////////////////////
            TVector3 p_q = P_b - p_e; // 3-momentum transfer (same as Erin's code)
            double theta_q = p_q.Theta() * 180 / M_PI;
            double nu = Ebeam - p_e.Mag();       // Energy transfer (same as Erin's code)
            double QSq = p_q.Mag2() - (nu * nu); // 4-momentum transfer squared (same as Erin's code)
            double xB = QSq / (2 * mN * nu);     // x Bjorken (same as Erin's code)
            double WSq = (mN * mN) - QSq + (2 * nu * mN);
            double theta_e = p_e.Theta() * 180 / M_PI;
            */

#pragma endregion /* Electrons - end */

#pragma region /* Protons - start */

            /*
            // From Erin's code
            if (protons.size() != 1) // One proton in event
            {
                continue;
            }

            int p_index = -1;
            TVector3 p_p(0., 0., 0.);

            // technically not optimized - this doesn't address what happens if there are two protons passing cuts
            // TODO: recheck this!
            for (int i = 0; i < protons.size(); i++)
            {
                // define quantities
                p_p.SetMagThetaPhi(protons[i]->getP(), protons[i]->getTheta(), protons[i]->getPhi());
                double dbeta = protons[i]->par()->getBeta() - p_p.Mag() / sqrt(p_p.Mag2() + mP * mP);
                double p_theta = p_p.Theta() * 180. / M_PI;
                double vzp = protons[i]->par()->getVz();
                double chipid = protons[i]->par()->getChi2Pid();

                // fill histos
                h_pangles->Fill(p_p.Phi() * 180. / M_PI, p_theta);

                if (protons[i]->getRegion() == FD)
                {
                    h_vzp_fd->Fill(vzp - vze);

                    if (fabs(vzp - vze) > 5)
                    // if (abs(vzp - vze) > 5) // Erin's original
                    {
                        continue;
                    }

                    h_chipid_fd->Fill(chipid);
                    h_dbeta_p_fd->Fill(p_p.Mag(), dbeta);

                    if (p_p.Mag() < 0.5)
                    {
                        continue;
                    }

                    if (p_p.Mag() > 3.0)
                    {
                        continue;
                    }

                    if (fabs(dbeta) > 0.03)
                    // if (abs(dbeta) > 0.03) // Erin's original
                    {
                        continue;
                    }
                }
                else if (protons[i]->getRegion() == CD)
                {
                    h_vzp_cd->Fill(vzp - vze);

                    if (fabs(vzp - vze) > 4)
                    // if (abs(vzp - vze) > 4) // Erin's original
                    {
                        continue;
                    }

                    h_chipid_cd->Fill(chipid);
                    // if (abs(chipid)>4) {continue;}
                    h_dbeta_p_cd->Fill(p_p.Mag(), dbeta);

                    if (p_p.Mag() < 0.3)
                    {
                        continue;
                    }

                    if (p_p.Mag() > 1.5)
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

            p_p.SetMagThetaPhi(protons[p_index]->getP(), protons[p_index]->getTheta(), protons[p_index]->getPhi());

            if (p_p.Theta() * 180. / M_PI < 40 || p_p.Theta() * 180. / M_PI > 140) // p goes to CD (CD proton cut)
            {
                continue;
            }
            // if (p_p.Theta()*180./M_PI>40) {continue;}  // p goes to FD

            // Andrew's original - commented out!
            // // Lead Proton
            // int num_L = 0;
            // int index_L = -1;

            // for (int j = 0; j < AllParticles.size(); j++)
            // {
            //     if ((LeadFDProton_Cut(c12, Ebeam, j)) || (LeadCDProton_Cut(c12, Ebeam, j)))
            //     {
            //         num_L++;
            //         index_L = j;
            //     }
            // }

            // if (num_L != 1)
            // {
            //     continue;
            // }

            // bool LeadCD = LeadCDProton_Cut(c12, Ebeam, index_L);
            // bool LeadFD = LeadFDProton_Cut(c12, Ebeam, index_L);

            // if (LeadCD && LeadFD)
            // {
            //     cout << "Problem!\n";
            // }

            // Andrew's original (p_L) was replaced by p_p from Erin's code
            // TVector3 p_L; // Andrew's original
            // p_L.SetMagThetaPhi(AllParticles[index_L]->getP(), AllParticles[index_L]->getTheta(), AllParticles[index_L]->getPhi()); // Andrew's original
            */

#pragma endregion /* Protons - end */

#pragma region /* Missing momentum - start */

            /*
            TVector3 p_miss = p_q - p_L;
            double Ep = sqrt(mN * mN + pp.Mag2());
            double Emiss = Ebeam + mD - pe.Mag() - Ep;
            double mmiss = sqrt((Emiss * Emiss) - p_miss.Mag2());
            // double mmiss = get_mmiss(p_b, p_e, p_L); // Andrew's original - commented out!
            */

#pragma region /* Missing momentum cuts (Andrew) - start */

            if (P_miss.Theta() * 180 / M_PI < 40)
            {
                continue;
            }

            if (P_miss.Theta() * 180 / M_PI > 135)
            {
                continue;
            }

            if (P_miss.Mag() < 0.2)
            {
                continue;
            }

            if (P_miss.Mag() > 1.25)
            {
                continue;
            }

            if (M_miss < 0.7)
            {
                continue;
            }

            if (M_miss > 1.2)
            {
                continue;
            }

#pragma endregion /* Missing momentum cuts (Andrew) - end */

            //////////////////////////////////////////////////
            // For after checking the hipo banks
            //////////////////////////////////////////////////

            // TODO: Do we need this?
            // int num_Charge = 0;

            // for (int j = 0; j < AllParticles.size(); j++)
            // {
            //     if (j == 0)
            //     {
            //         continue;
            //     }

            //     if (j == index_L)
            //     {
            //         continue;
            //     }

            //     // if(j==index_Rp1){continue;}
            //     if (AllParticles[j]->par()->getCharge() == 0)
            //     {
            //         continue;
            //     }

            //     num_Charge++;
            // }

            // TODO: need this?
            // if (num_Charge > 0)
            // {
            //     continue;
            // }

            if (pInFD)
            {
                h_xB_mmiss_epFD->Fill(xB, M_miss, weight);
            }
            else if (pInCD)
            {
                h_xB_mmiss_epCD->Fill(xB, M_miss, weight);
            }

            h_pmiss_ep->Fill(P_miss.Mag(), weight);

            if (M_miss > 1.05) // TODO: Missing mass cut; why is this cut being applied twice?
            {
                continue;
            }

            // TODO: check if works!
            if (pInCD && (xB < 1.1)) // Cutting out CD leading protons with xB < 1.1
            {
                continue;
            }

            // TODO: check if works!
            if (P_p.Theta() * 180. / M_PI < 40) // Cutting out FD protons
            {
                continue;
            }
            // if(LeadFD && (xB<0.8)){continue;}

#pragma endregion /* Missing momentum - end */

#pragma region /* Neutrons */

            /////////////////////////////////////
            // Lead Neutron Checks
            /////////////////////////////////////
            for (int itr1 = 0; itr1 < AllParticles.size(); itr1++)
            {
                if (AllParticles[itr1]->par()->getCharge() != 0) // Cut out charged particles
                {
                    continue;
                }

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

                double theta = AllParticles[itr1]->getTheta() * 180 / M_PI;
                double beta = AllParticles[itr1]->par()->getBeta();
                double gamma = 1 / sqrt(1 - (beta * beta));
                double mom = gamma * beta * mN;
                double ToF = AllParticles[itr1]->getTime() - starttime;

                int detINTlayer = C1 ? 1 : C2 ? 2
                                              : 3;
                auto detlayer = C1 ? CND1 : C2 ? CND2
                                               : CND3; // CND layer with hit
                double edep = AllParticles[itr1]->sci(CND1)->getEnergy() + AllParticles[itr1]->sci(CND2)->getEnergy() + AllParticles[itr1]->sci(CND3)->getEnergy();
                double edep_CTOF = AllParticles[itr1]->sci(CTOF)->getEnergy();
                double edep_single = AllParticles[itr1]->sci(detlayer)->getEnergy();

                double nvtx_x = AllParticles[itr1]->par()->getVx();
                double nvtx_y = AllParticles[itr1]->par()->getVy();
                double nvtx_z = AllParticles[itr1]->par()->getVz();
                TVector3 v_nvtx(nvtx_x, nvtx_y, nvtx_z); // Neutron's vertex location

                TVector3 v_hit;
                v_hit.SetXYZ(AllParticles[itr1]->sci(detlayer)->getX(), AllParticles[itr1]->sci(detlayer)->getY(), AllParticles[itr1]->sci(detlayer)->getZ()); // Neutron's hit location in CND

                TVector3 v_path = v_hit - v_nvtx; // Direct calculation of neutron's path (in vector form)
                TVector3 v_n;
                v_n.SetMagThetaPhi(mom, v_path.Theta(), v_path.Phi()); // Direct calculation of neutron momentum?
                // TODO: check with Andrew why he calculated this explicitly

                double path = v_path.Mag() / 100;
                double theta_nmiss = v_n.Angle(P_miss) * 180 / M_PI; // Opening angle between calculated neutron's momentum and predicted neutron momentum (= missing momentum)
                double dm_nmiss = (P_miss.Mag() - v_n.Mag()) / P_miss.Mag();
                int nSector = AllParticles[itr1]->sci(detlayer)->getSector(); // Number of CND sector with a neutron hit in the layer detlayer

                // Check to see if there is a good neutron
                bool isGN = false;

                if ((theta_nmiss < 40) && (dm_nmiss > -0.5) && (dm_nmiss < 0.5)) // Good neutron definition
                {
                    isGN = true;
                }

                //////////////////////////////////////////////
                // Step Zero
                //////////////////////////////////////////////
                // Why "path * 100"? unit conversion. Path is in cm; tof is in ns.
                if (fabs(beta - (path * 100) / (ToF * c)) > 0.01) // A cut on delta beta
                {
                    continue;
                }

                // Andrew's original
                // if (beta - (path * 100) / (ToF * c) < -0.01)
                // {
                //     continue;
                // }

                // Andrew's original
                // if (beta - (path * 100) / (ToF * c) > 0.01)
                // {
                //     continue;
                // }

                // A cut on the z-component of the CND hit
                // This is a fiducial cut on the range that the CND can reach on the z-axis
                if (v_hit.Z() > 45 || v_hit.Z() < -40)
                {
                    continue;
                }

                // Andrew's original
                // if (v_hit.Z() > 45)
                // {
                //     continue;
                // }

                // Andrew's original
                // if (v_hit.Z() < -40)
                // {
                //     continue;
                // }

                if (ToF < 0 || ToF > 20)
                {
                    continue;
                }

                // Andrew's original
                // if (ToF < 0)
                // {
                //     continue;
                // }

                // Andrew's original
                // if (ToF > 20)
                // {
                //     continue;
                // }

                if (pInFD)
                {
                    h_xB_mmiss_epnFD->Fill(xB, M_miss, weight);
                }
                else if (pInCD)
                {
                    h_xB_mmiss_epnCD->Fill(xB, M_miss, weight);
                }

                h_pnRes_theta_nmiss_Step0->Fill(dm_nmiss, theta_nmiss, weight);

                if (isGN)
                {
                    if (pInFD)
                    {
                        h_xB_mmiss_epn_goodN_pFD->Fill(xB, M_miss, weight);
                    }
                    else if (pInFD)
                    {
                        h_xB_mmiss_epn_goodN_pCD->Fill(xB, M_miss, weight);
                    }

                    h_ToF_goodN_Step0->Fill(ToF, weight);
                    h_beta_goodN_Step0->Fill(beta, weight);
                    h_Edep_goodN_Step0->Fill(edep, weight);
                    h_beta_Edep_goodN_Step0->Fill(beta, edep, weight);
                    h_Edep_ToF_goodN_Step0->Fill(ToF, edep, weight);
                }
                else
                {
                    h_ToF_badN_Step0->Fill(ToF, weight);
                    h_beta_badN_Step0->Fill(beta, weight);
                    h_Edep_badN_Step0->Fill(edep, weight);
                    h_beta_Edep_badN_Step0->Fill(beta, edep, weight);
                    h_Edep_ToF_badN_Step0->Fill(ToF, edep, weight);
                }

                //////////////////////////////////////////////
                // Step One
                //////////////////////////////////////////////

                // Step One = Beta cut & Dep. energy cut

                if (beta > 0.8) // Beta cut
                {
                    continue;
                }

                // Dep. energy cut
                // TODO: check if should be 12 MeV
                if (edep < 5)
                {
                    continue;
                }

                h_pnRes_theta_nmiss_Step1->Fill(dm_nmiss, theta_nmiss, weight);

                if (isGN) // Surviving good neutrons after step one
                {
                    h_ToF_goodN_Step1->Fill(ToF, weight);
                    h_pmiss_goodN_Step1->Fill(P_miss.Mag(), weight);
                    h_Edep_ToF_goodN_Step1->Fill(ToF, edep, weight);
                }
                else
                {
                    h_ToF_badN_Step1->Fill(ToF, weight);
                    h_Edep_ToF_badN_Step1->Fill(ToF, edep, weight);
                }

                bool CNDVeto = false;

                if (ToF * c - v_hit.Z() < 70) // TODO: ask Andrew why this cut
                {

                    if (isGN)
                    {
                        h_edep_goodN_Step1->Fill(edep, weight);
                    }
                    else
                    {
                        h_edep_badN_Step1->Fill(edep, weight);
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
                        int vetoSectorbyLayer[4] = {(AllParticles[itr2]->sci(CTOF)->getComponent() + 1) / 2,
                                                    AllParticles[itr2]->sci(CND1)->getSector(),
                                                    AllParticles[itr2]->sci(CND2)->getSector(),
                                                    AllParticles[itr2]->sci(CND3)->getSector()};

                        TVector3 p_C;
                        p_C.SetMagThetaPhi(AllParticles[itr2]->getP(), AllParticles[itr2]->getTheta(), AllParticles[itr2]->getPhi());

                        double edep_pos = AllParticles[itr2]->sci(clas12::CTOF)->getEnergy(); // E_dep of positivly charged particle

                        for (int itr3 = 0; itr3 < 4; itr3++) //
                        {
                            if (vetoSectorbyLayer[itr3] == 0) // TODO: why this cut? no hit in the itr3-th layer?
                            {
                                continue;
                            }

                            int sdiff = nSector - vetoSectorbyLayer[itr3];

                            if (sdiff <= -12)
                            {
                                sdiff += 24;
                            }
                            else if (sdiff > 12)
                            {
                                sdiff -= 24;
                            }

                            int ldiff = detINTlayer - itr3;

                            if (isGN) // ldiff + 3 == 0 -> first element in h_sdiff_pos_goodN_Step1_layer
                            {
                                h_sdiff_pos_goodN_Step1_layer[ldiff + 3]->Fill(sdiff, weight);
                                h_sdiff_pos_mom_goodN_Step1_layer[ldiff + 3]->Fill(sdiff, p_C.Perp(), weight);
                                h_sdiff_pos_z_goodN_Step1_layer[ldiff + 3]->Fill(sdiff, v_hit.Z(), weight);
                                h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[ldiff + 3]->Fill(sdiff, ToF * c - v_hit.Z(), weight);
                            }
                            else
                            {
                                h_sdiff_pos_badN_Step1_layer[ldiff + 3]->Fill(sdiff, weight);
                                h_sdiff_pos_mom_badN_Step1_layer[ldiff + 3]->Fill(sdiff, p_C.Perp(), weight);
                                h_sdiff_pos_z_badN_Step1_layer[ldiff + 3]->Fill(sdiff, v_hit.Z(), weight);
                                h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[ldiff + 3]->Fill(sdiff, ToF * c - v_hit.Z(), weight);
                            }

                            if (isPosNear(sdiff, ldiff))
                            {
                                CNDVeto = true;
                            }
                        } // End of loop over vetoSectorbyLayer

                        if (CNDVeto)
                        {
                            if (isGN)
                            {
                                h_edep_over_edepCTOT_goodN_Step1->Fill(edep / edep_pos, weight);
                            }
                            else
                            {
                                h_edep_over_edepCTOT_badN_Step1->Fill(edep / edep_pos, weight);
                            }
                        }
                    } // End of second loop over AllParticles (step 1)

                    if (CNDVeto)
                    {
                        if (isGN)
                        {
                            h_edep_goodN_withNearbyPos_Step1->Fill(edep, weight);
                        }
                        else
                        {
                            h_edep_badN_withNearbyPos_Step1->Fill(edep, weight);
                        }
                    }

                    if (isGN)
                    {
                        if (!CNDVeto)
                            h_diff_ToFc_z_Edep_noNear_goodN_Step1->Fill(ToF * c - v_hit.Z(), edep, weight);
                        else
                        {
                            h_diff_ToFc_z_Edep_yesNear_goodN_Step1->Fill(ToF * c - v_hit.Z(), edep, weight);
                        }
                    }
                    else
                    {
                        if (!CNDVeto)
                            h_diff_ToFc_z_Edep_noNear_badN_Step1->Fill(ToF * c - v_hit.Z(), edep, weight);
                        else
                        {
                            h_diff_ToFc_z_Edep_yesNear_badN_Step1->Fill(ToF * c - v_hit.Z(), edep, weight);
                        }
                    }
                }

                //////////////////////////////////////////////
                // Step Two
                //////////////////////////////////////////////

                // Step two = cut/veto out neutrons with charged particles close by

                if (CNDVeto)
                {
                    continue;
                }

                h_pnRes_theta_nmiss_Step2->Fill(dm_nmiss, theta_nmiss, weight);

                if (isGN)
                {
                    h_ToF_goodN_Step2->Fill(ToF, weight);
                    h_Edep_ToF_goodN_Step2->Fill(ToF, edep, weight);
                }
                else
                {
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

                /*
                bool AllHitVeto = false;
                int hitsNear=0;

                for( int row = 0; row < c12->getBank(cnd_hits)->getRows();row++){
                  int hit_sector = c12->getBank(cnd_hits)->getInt(cnd_hit_sector,row);
                  int hit_layer = c12->getBank(cnd_hits)->getInt(cnd_hit_layer,row);
                  double hit_energy = c12->getBank(cnd_hits)->getFloat(cnd_hit_energy,row);

                  int sdiff = nSector - hit_sector;
                  if(sdiff<=-12){sdiff+=24;}
                  else if(sdiff>12){sdiff-=24;}
                  int ldiff = detINTlayer - hit_layer;

                  if((ldiff==0) && (sdiff==0)){continue;}
                  if(isNear(sdiff,ldiff)){
                    if(hit_energy>20){
                    if(isGN){
                      h_NearbyEdep_goodN_Step2->Fill(hit_energy,weight);
                    }
                    else{
                      h_NearbyEdep_badN_Step2->Fill(hit_energy,weight);
                    }
                      hitsNear++;
                    }
                  }

                  //if(ToF*c-v_hit.Z() < 70){
                if(ToF<6){
                  if(isGN){
                    h_sdiff_allhit_goodN_Step2_layer[ldiff+3]->Fill(sdiff,weight);
                    h_sdiff_ldiff_allhit_goodN_Step2->Fill(sdiff,ldiff,weight);
                  }
                  else{
                    h_sdiff_allhit_badN_Step2_layer[ldiff+3]->Fill(sdiff,weight);
                    h_sdiff_ldiff_allhit_badN_Step2->Fill(sdiff,ldiff,weight);
                  }
                }
                //}
                }

                //Now that we have nearby hits, look at how many we have that are not off time
                //if(ToF*c-v_hit.Z() < 70){
                if(ToF<6){
                if(isGN){
                  h_numberNearby_goodN_Step2->Fill(hitsNear,weight);
                  h_numberNearby_momN_goodN_Step2->Fill(hitsNear,mom,weight);
                  h_nsector_goodN_Step2->Fill(nSector,weight);
                }
                else{
                  h_numberNearby_badN_Step2->Fill(hitsNear,weight);
                  h_numberNearby_momN_badN_Step2->Fill(hitsNear,mom,weight);
                  h_nsector_badN_Step2->Fill(nSector,weight);
                }
                }
                //}
                if(hitsNear>=1){AllHitVeto=true;}

                //////////////////////////////////////////////
                //Step 3
                //////////////////////////////////////////////
                if(AllHitVeto){continue;}
                bool CTOFHitVeto = false;
                int hitsCTOF=0;
                h_pnRes_theta_nmiss_Step3->Fill(dm_nmiss,theta_nmiss,weight);
                if(isGN){
                  h_ToF_goodN_Step3->Fill(ToF,weight);
                  h_Edep_ToF_goodN_Step3->Fill(ToF,edep,weight);
                }
                else{
                  h_ToF_badN_Step3->Fill(ToF,weight);
                  h_Edep_ToF_badN_Step3->Fill(ToF,edep,weight);
                }

                for( int row = 0; row < c12->getBank(cnd_hits)->getRows();row++){
                  int hit_sector = c12->getBank(cnd_hits)->getInt(cnd_hit_sector,row);
                  int hit_layer = c12->getBank(cnd_hits)->getInt(cnd_hit_layer,row);
                  int sdiff = nSector - hit_sector;
                  if(sdiff<=-12){sdiff+=24;}
                  else if(sdiff>12){sdiff-=24;}
                  int ldiff = detINTlayer - hit_layer;
                  if((ldiff==0) && (sdiff==0)){continue;}
                  if(isGN){h_sdiff_ldiff_allhit_goodN_Step3->Fill(sdiff,ldiff,weight);}
                  else{h_sdiff_ldiff_allhit_badN_Step3->Fill(sdiff,ldiff,weight);}
                }

                for( int row = 0; row < c12->getBank(ctof_hits)->getRows();row++){
                  int hit_sector = (c12->getBank(ctof_hits)->getInt(ctof_hit_component,row)+1)/2;
                  double hit_energy = c12->getBank(ctof_hits)->getFloat(ctof_hit_energy,row);

                  int sdiff = nSector - hit_sector;
                  if(sdiff<=-12){sdiff+=24;}
                  else if(sdiff>12){sdiff-=24;}
                  int ldiff = detINTlayer;

                  if((ldiff==0) && (sdiff==0)){continue;}
                  if(isNearCTOF(sdiff,ldiff)){
                    //if(isGN){h_NearbyEdep_goodN_Step3->Fill(hit_energy,weight);}
                    //else{h_NearbyEdep_badN_Step3->Fill(hit_energy,weight);}
                    if(hit_energy>5){hitsCTOF++;}
                  }
                  if(isGN){
                    h_sdiff_ldiff_CTOFhit_goodN_Step3->Fill(sdiff,ldiff,weight);
                  }
                  else{
                    h_sdiff_ldiff_CTOFhit_badN_Step3->Fill(sdiff,ldiff,weight);
                  }
                }

                if(isGN){
                  h_numberCTOF_goodN_Step3->Fill(hitsCTOF,weight);
                  h_numberCTOF_momN_goodN_Step3->Fill(hitsCTOF,mom,weight);
                }
                else{
                  h_numberCTOF_badN_Step3->Fill(hitsCTOF,weight);
                  h_numberCTOF_momN_badN_Step3->Fill(hitsCTOF,mom,weight);
                }
                if(hitsCTOF>=1){CTOFHitVeto=true;}

                //////////////////////////////////////////////
                //Step 4
                //////////////////////////////////////////////
                if(CTOFHitVeto){continue;}

                h_pnRes_theta_nmiss_Step4->Fill(dm_nmiss,theta_nmiss,weight);
                if(isGN){
                  h_ToF_goodN_Step4->Fill(ToF,weight);
                  h_edep_ToF_goodN_Step4->Fill(ToF,edep,weight);
                }
                else{
                  h_ToF_badN_Step4->Fill(ToF,weight);
                  h_Edep_ToF_badN_Step3->Fill(ToF,edep,weight);
                }

                ///////////////////
                h_pnRes_theta_nmiss_Step5->Fill(dm_nmiss,theta_nmiss,weight);
                h_pmiss_allN_Step5->Fill(p_miss.Mag(),weight);
                if(isGN){
                  h_ToF_goodN_Step5->Fill(ToF,weight);
                  h_edep_ToF_goodN_Step5->Fill(ToF,edep,weight);
                  h_pmiss_goodN_Step5->Fill(p_miss.Mag(),weight);
                  h_diff_ToFc_z_Edep_goodN_Step5->Fill(ToF*c-v_hit.Z(),edep,weight);
                  h_diff_ToFc_z_Edep_goodN_Step5_layer[detINTlayer-1]->Fill(ToF*c-v_hit.Z(),edep,weight);
                  h_phidiff_en_goodN_Step5->Fill(get_phi_diff(p_e,v_n),weight);
                  h_TP_goodN_Step5->Fill(ToF/path,weight);
                  h_Z_goodN_Step5->Fill(v_hit.Z(),weight);
                  h_beta_Edep_goodN_Step5->Fill(beta,edep,weight);

                  h_ToF_Edep_goodN_Step5->Fill(ToF,edep,weight);
                  h_TP_Edep_goodN_Step5->Fill(ToF/path,edep,weight);
                }
                else{
                  h_ToF_badN_Step5->Fill(ToF,weight);
                  h_edep_ToF_badN_Step5->Fill(ToF,edep,weight);
                  //if(ToF<8){
                  h_diff_ToFc_z_Edep_badN_Step5->Fill(ToF*c-v_hit.Z(),edep,weight);
                  h_diff_ToFc_z_Edep_badN_Step5_layer[detINTlayer-1]->Fill(ToF*c-v_hit.Z(),edep,weight);
                  h_phidiff_en_badN_Step5->Fill(get_phi_diff(p_e,v_n),weight);
                  h_TP_badN_Step5->Fill(ToF/path,weight);
                  h_Z_badN_Step5->Fill(v_hit.Z(),weight);
                  h_beta_Edep_badN_Step5->Fill(beta,edep,weight);

                  h_ToF_Edep_badN_Step5->Fill(ToF,edep,weight);
                  h_TP_Edep_badN_Step5->Fill(ToF/path,edep,weight);
                  //}
                  if(ToF<5){
                    if(edep>20){
                      //cerr<<"Event="<<c12->runconfig()->getEvent()<<endl;
                      //cerr<<"Neutron Sector = "<<nSector<<endl;
                      //cerr<<"Neutron Z Hit = "<<v_hit.Z()<<endl;
                    }
                  }
                }

                for( int row = 0; row < c12->getBank(cnd_hits)->getRows();row++){
                  int hit_sector = c12->getBank(cnd_hits)->getInt(cnd_hit_sector,row);
                  int hit_layer = c12->getBank(cnd_hits)->getInt(cnd_hit_layer,row);
                  double hit_energy = c12->getBank(cnd_hits)->getFloat(cnd_hit_energy,row);

                  int sdiff = nSector - hit_sector;
                  if(sdiff<=-12){sdiff+=24;}
                  else if(sdiff>12){sdiff-=24;}
                  int ldiff = detINTlayer - hit_layer;

                  if((ldiff==1) && (sdiff==0)){
                    if(isGN){
                      h_Edep_infront_goodN_Step5->Fill(hit_energy,weight);
                    }
                    else{
                      h_Edep_infront_badN_Step5->Fill(hit_energy,weight);
                    }
                  }

                  else if((ldiff==-1) && (sdiff==0)){
                    if(isGN){
                      h_Edep_behind_goodN_Step5->Fill(hit_energy,weight);
                    }
                    else{
                      h_Edep_behind_badN_Step5->Fill(hit_energy,weight);
                    }
                  }
                }
                for( int row = 0; row < c12->getBank(ctof_hits)->getRows();row++){
                  int hit_sector = (c12->getBank(ctof_hits)->getInt(ctof_hit_component,row)+1)/2;
                  double hit_energy = c12->getBank(ctof_hits)->getFloat(ctof_hit_energy,row);

                  int sdiff = nSector - hit_sector;
                  if(sdiff<=-12){sdiff+=24;}
                  else if(sdiff>12){sdiff-=24;}
                  int ldiff = detINTlayer;

                  if((ldiff==1) && (sdiff==0)){
                    if(isGN){
                      h_Edep_infront_goodN_Step5->Fill(hit_energy,weight);
                    }
                    else{
                      h_Edep_infront_badN_Step5->Fill(hit_energy,weight);
                    }
                  }

                  else if((ldiff==-1) && (sdiff==0)){
                    if(isGN){
                      h_Edep_behind_goodN_Step5->Fill(hit_energy,weight);
                    }
                    else{
                      h_Edep_behind_badN_Step5->Fill(hit_energy,weight);
                    }
                  }

                }
                */
                /*
                if(!isGN){
                  h_ToF_badN->Fill(ToF,weight);
                  h_ToF_z_badN->Fill(ToF,v_hit.Z(),weight);
                  if(ToF<10){
                    h_TM_badN->Fill(ToF/path,weight);
                    h_beta_badN->Fill(beta,weight);
                    h_mom_badN->Fill(v_n.Mag(),weight);
                    h_Edep_z_badN->Fill(edep_single,v_hit.Z(),weight);
                    h_Edep_ToF_badN->Fill(edep_single,ToF,weight);
                    h_beta_z_badN->Fill(beta,v_hit.Z(),weight);
                    if(C1){
                      h_ToFc_z_1_badN->Fill(ToF*c,v_hit.Z(),weight);
                      h_diff_ToFc_z_1_badN->Fill(ToF*c-v_hit.Z(),weight);
                      h_diff_ToFc_z_Edep_1_badN->Fill(ToF*c-v_hit.Z(),edep,weight);
                      h_Edep_z_1_badN->Fill(edep_single,v_hit.Z(),weight);}
                    if(C2){
                      h_ToFc_z_2_badN->Fill(ToF*c,v_hit.Z(),weight);
                      h_diff_ToFc_z_2_badN->Fill(ToF*c-v_hit.Z(),weight);
                      h_diff_ToFc_z_Edep_2_badN->Fill(ToF*c-v_hit.Z(),edep,weight);
                      h_Edep_z_2_badN->Fill(edep_single,v_hit.Z(),weight);}
                    if(C3){
                      h_ToFc_z_3_badN->Fill(ToF*c,v_hit.Z(),weight);
                      h_diff_ToFc_z_3_badN->Fill(ToF*c-v_hit.Z(),weight);
                      h_diff_ToFc_z_Edep_3_badN->Fill(ToF*c-v_hit.Z(),edep,weight);
                      h_Edep_z_3_badN->Fill(edep_single,v_hit.Z(),weight);}

                    h_Edep_mom_badN->Fill(edep,mom,weight);

                      if(v_hit.Z()>10){
                    cerr<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
                    cerr<<"Run="<<c12->runconfig()->getRun()<<endl;
                    cerr<<"Event="<<c12->runconfig()->getEvent()<<endl;
                    cerr<<"Neutron Sector = "<<nSector<<endl;
                    cerr<<"Neutron Z Hit = "<<v_hit.Z()<<endl;
                    cerr<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl<<endl<<endl;
                   }

                  }
                }
                else{
                  h_ToF_goodN->Fill(ToF,weight);
                  h_ToF_z_goodN->Fill(ToF,v_hit.Z(),weight);
                  if(ToF<10){
                    h_TM_goodN->Fill(ToF/path,weight);
                    h_beta_goodN->Fill(beta,weight);
                    h_mom_goodN->Fill(v_n.Mag(),weight);
                    h_Edep_z_goodN->Fill(edep_single,v_hit.Z(),weight);
                    h_Edep_ToF_goodN->Fill(edep_single,ToF,weight);
                    h_beta_z_goodN->Fill(beta,v_hit.Z(),weight);
                    if(C1){
                      h_ToFc_z_1_goodN->Fill(ToF*c,v_hit.Z(),weight);
                      h_diff_ToFc_z_1_goodN->Fill(ToF*c-v_hit.Z(),weight);
                      h_diff_ToFc_z_Edep_1_goodN->Fill(ToF*c-v_hit.Z(),edep,weight);
                      h_Edep_z_1_goodN->Fill(edep_single,v_hit.Z(),weight);}
                    if(C2){
                      h_ToFc_z_2_goodN->Fill(ToF*c,v_hit.Z(),weight);
                      h_diff_ToFc_z_2_goodN->Fill(ToF*c-v_hit.Z(),weight);
                      h_diff_ToFc_z_Edep_2_goodN->Fill(ToF*c-v_hit.Z(),edep,weight);
                      h_Edep_z_2_goodN->Fill(edep_single,v_hit.Z(),weight);}
                    if(C3){
                      h_ToFc_z_3_goodN->Fill(ToF*c,v_hit.Z(),weight);
                      h_diff_ToFc_z_3_goodN->Fill(ToF*c-v_hit.Z(),weight);
                      h_diff_ToFc_z_Edep_3_goodN->Fill(ToF*c-v_hit.Z(),edep,weight);
                      h_Edep_z_3_goodN->Fill(edep_single,v_hit.Z(),weight);}
                  }
                  h_Edep_mom_goodN->Fill(edep,mom,weight);
                }

                if(edep<12.5){continue;}
                  if((edep<-(40.0/110.0)*((ToF*c-v_hit.Z())-110)) && C1){continue;}
                  if((edep<-(32.0/110.0)*((ToF*c-v_hit.Z())-110)) && C2){continue;}
                  if((edep<-(26.0/110.0)*((ToF*c-v_hit.Z())-110)) && C3){continue;}

                */

                // if(C3 && (v_hit.Z()>25)){continue;}
                // else if(C2 && (v_hit.Z()>20)){continue;}
                // else if(C1 && (v_hit.Z()>10)){continue;}
            } // End of Andrew's loop over all particles

#pragma endregion /* Neutrons */
        }

#pragma endregion /* Andrew's manual work */

    } // closes event loop

#pragma endregion /* Chain loop - end */

    // ======================================================================================================================================================================
    // Andrew's wrap up
    // ======================================================================================================================================================================

#pragma region /* Andrew's wrap up - start */

    /////////////////////////////////////////////////////
    // Now create the output PDFs
    /////////////////////////////////////////////////////

    int pixelx = 1980;
    int pixely = 1530;

    TCanvas *myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);
    TCanvas *myText = new TCanvas("myText", "myText", pixelx, pixely);

#pragma region /* Saving all plots - start */

    TLatex text;
    text.SetTextSize(0.05);

    const char *pdfFile = PDFFile.c_str();

    char fileName[100];
    sprintf(fileName, "%s[", pdfFile);
    myText->SaveAs(fileName);
    sprintf(fileName, "%s", pdfFile);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text.DrawLatex(0.2, 0.7, "Neutrons in CND");

    myText->Print(fileName, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();

    double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    double diplayTextSize = 0.1;

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        myCanvas->cd(1);
        hist_list_1_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_A[i]->GetXaxis()->SetTitleSize(0.06);
        hist_list_1_A[i]->GetXaxis()->SetLabelSize(0.0425);
        hist_list_1_A[i]->GetYaxis()->CenterTitle();
        hist_list_1_A[i]->GetYaxis()->SetTitleSize(0.06);
        hist_list_1_A[i]->GetYaxis()->SetLabelSize(0.0425);
        hist_list_1_A[i]->SetLineWidth(2);
        hist_list_1_A[i]->SetLineColor(kBlue);

        if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
        {
            TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
            displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
            hist_list_1_A[i]->Draw(), displayText->Draw("same");
        }
        else
        {
            hist_list_1_A[i]->Draw();
        }

        myCanvas->Print(fileName, "pdf");
        myCanvas->Clear();
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        myCanvas->cd(1);
        hist_list_2_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_A[i]->GetXaxis()->SetTitleSize(0.06);
        hist_list_2_A[i]->GetXaxis()->SetLabelSize(0.0425);
        hist_list_2_A[i]->GetYaxis()->CenterTitle();
        hist_list_2_A[i]->GetYaxis()->SetTitleSize(0.06);
        hist_list_2_A[i]->GetYaxis()->SetLabelSize(0.0425);

        if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
        {
            TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
            displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
            hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
        }
        else
        {
            hist_list_2_A[i]->Draw("colz");
        }

        myCanvas->Print(fileName, "pdf");
        myCanvas->Clear();
    }

    sprintf(fileName, "%s]", pdfFile);
    myCanvas->Print(fileName, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving all plots - end */

#pragma region /* Saving Step0 plots - start */

    TLatex text_Step0;
    text_Step0.SetTextSize(0.05);

    const char *pdfFile_Step0 = (ConfigOutPutName(PDFFile, "Step0")).c_str();

    char fileName_Step0[100];
    sprintf(fileName_Step0, "%s[", pdfFile_Step0);
    myText->SaveAs(fileName_Step0);
    sprintf(fileName_Step0, "%s", pdfFile_Step0);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text_Step0.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text_Step0.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text_Step0.DrawLatex(0.2, 0.7, "Neutrons in CND - step 0");

    myText->Print(fileName_Step0, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();

    for (int i = 0; i < hist_list_1_Step0_A.size(); i++)
    {
        myCanvas->cd(1);
        hist_list_1_Step0_A[i]->GetXaxis()->CenterTitle();
        hist_list_1_Step0_A[i]->GetXaxis()->SetTitleSize(0.06);
        hist_list_1_Step0_A[i]->GetXaxis()->SetLabelSize(0.0425);
        hist_list_1_Step0_A[i]->GetYaxis()->CenterTitle();
        hist_list_1_Step0_A[i]->GetYaxis()->SetTitleSize(0.06);
        hist_list_1_Step0_A[i]->GetYaxis()->SetLabelSize(0.0425);
        hist_list_1_Step0_A[i]->SetLineWidth(2);
        hist_list_1_Step0_A[i]->SetLineColor(kBlue);

        if (hist_list_1_Step0_A[i]->GetEntries() == 0 || hist_list_1_Step0_A[i]->Integral() == 0)
        {
            TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
            displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
            hist_list_1_Step0_A[i]->Draw(), displayText->Draw("same");
        }
        else
        {
            hist_list_1_Step0_A[i]->Draw();
        }

        myCanvas->Print(fileName_Step0, "pdf");
        myCanvas->Clear();
    }

    for (int i = 0; i < hist_list_2_Step0_A.size(); i++)
    {
        myCanvas->cd(1);
        hist_list_2_Step0_A[i]->GetXaxis()->CenterTitle();
        hist_list_2_Step0_A[i]->GetXaxis()->SetTitleSize(0.06);
        hist_list_2_Step0_A[i]->GetXaxis()->SetLabelSize(0.0425);
        hist_list_2_Step0_A[i]->GetYaxis()->CenterTitle();
        hist_list_2_Step0_A[i]->GetYaxis()->SetTitleSize(0.06);
        hist_list_2_Step0_A[i]->GetYaxis()->SetLabelSize(0.0425);

        if (hist_list_2_Step0_A[i]->GetEntries() == 0 || hist_list_2_Step0_A[i]->Integral() == 0)
        {
            TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
            displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
            hist_list_2_Step0_A[i]->Draw("colz"), displayText->Draw("same");
        }
        else
        {
            hist_list_2_Step0_A[i]->Draw("colz");
        }

        myCanvas->Print(fileName_Step0, "pdf");
        myCanvas->Clear();
    }

    sprintf(fileName_Step0, "%s]", pdfFile_Step0);
    myCanvas->Print(fileName_Step0, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving Step0 plots - end */

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
    // Printouts
    // ======================================================================================================================================================================

#pragma region /* Printouts - start */

    cout << "\033[33m\n\033[0m";
    cout << "\033[33minput_hipo:\033[0m\t\t" << input_hipo << "\n";
    cout << "\033[33m\n\033[0m";
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