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
#include "src/functions/neutron-veto/veto_functions.cpp"
#include "src/functions/Andrews_functions/Andrews_functions.cpp"
#include "src/functions/HipoChain_config.cpp"
#include "src/classes/clas12ana/clas12ana.cpp"

using namespace std;
using namespace clas12;

#pragma region /* Erin main function */

int D_getfeatures_Phase3(double Ebeam, bool keep_good, string output_root, string output_txt, string input_hipo, /* Erin's arguments*/
                         string PDFFile, int isMC = 0 /* Andrew's arguments*/)
// int main(int argc, char **argv)
{
    const bool Run_Erins_features = true;
    const bool Run_Andrews_work = true;

    cout << endl;

    // ======================================================================================================================================================================
    // Initial setup
    // ======================================================================================================================================================================

#pragma region /* Initial setup */

    // arg 1: beam energy

    // arg 2: keep good

    // args 3-4: output file names
    TFile *f = new TFile(output_root.c_str(), "RECREATE");
    TTree *ntree = new TTree("T", "NeutronTree");
    std::ofstream outtxt(output_txt);

    // arg 5+: input hipo file
    clas12root::HipoChain chain;
    // TODO: add all run files in folder to the chain
    // for (int k = 5; k < argc; k++)
    // {
    //     std::cout << "Input file " << argv[k] << std::endl;
    //     chain.Add(argv[k]);
    // }

    chain.Add(input_hipo);

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

    clasAna->readEcalSFPar("src/cuts/paramsSF_LD2_x2.dat");
    clasAna->readEcalPPar("src/cuts/paramsPI_LD2_x2.dat");

    clasAna->printParams();

    clasAna->setProtonPidCuts(true);

#pragma endregion /* Initial setup */

    // ======================================================================================================================================================================
    // Erin's histograms
    // ======================================================================================================================================================================

#pragma region /* Erin's histograms */

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

#pragma endregion /* Erin's histograms */

    // ======================================================================================================================================================================
    // Andrew's histograms
    // ======================================================================================================================================================================

#pragma region /* Andrew's histograms */

    /////////////////////////////////////
    // Prepare histograms
    /////////////////////////////////////

    system(("rm -r " + PDFFile).c_str()); // Delete old output folder

    vector<TH1 *> hist_list_1_A;
    vector<TH2 *> hist_list_2_A;

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
    TH2D *h_xB_mmiss_epngoodFD = new TH2D("xB_mmiss_epngoodFD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epngoodFD);
    TH2D *h_xB_mmiss_epCD = new TH2D("xB_mmiss_epCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epCD);
    TH2D *h_xB_mmiss_epnCD = new TH2D("xB_mmiss_epnCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epnCD);
    TH2D *h_xB_mmiss_epngoodCD = new TH2D("xB_mmiss_epngoodCD", "x_{B} vs. m_{miss};x_{B};m_{miss}", 100, 0.0, 2.0, 100, 0.5, 1.5);
    hist_list_2_A.push_back(h_xB_mmiss_epngoodCD);

    TH1D *h_pmiss_ep = new TH1D("pmiss_ep", "p_{miss} ep;p_{miss};Counts", 25, 0.25, 1.0);
    hist_list_1_A.push_back(h_pmiss_ep);

    // Step Zero (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step0 = new TH2D("pnRes_theta_nmiss_Step0", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step0);

    TH1D *h_ToF_goodN_Step0 = new TH1D("ToF_goodN_Step0", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step0);
    TH1D *h_ToF_badN_Step0 = new TH1D("ToF_badN_Step0", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step0);

    TH1D *h_beta_goodN_Step0 = new TH1D("beta_goodN_Step0", "#beta of CND Neutrons;#beta;Counts", 50, 0, 1.1);
    hist_list_1_A.push_back(h_beta_goodN_Step0);
    TH1D *h_Edep_goodN_Step0 = new TH1D("Edep_goodN_Step0", "E_{dep} [MeF] of CND Neutrons;E_{dep};Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_goodN_Step0);
    TH2D *h_beta_Edep_goodN_Step0 = new TH2D("Edep_beta_goodN", "#beta vs. E_{dep} [MeV] of CND Neutrons;#beta;E_{dep}", 50, 0, 1.1, 50, 0, 100);
    hist_list_2_A.push_back(h_beta_Edep_goodN_Step0);

    TH1D *h_beta_badN_Step0 = new TH1D("beta_badN_Step0", "#beta of CND Neutrons;#beta;Counts", 50, 0, 1.1);
    hist_list_1_A.push_back(h_beta_badN_Step0);
    TH1D *h_Edep_badN_Step0 = new TH1D("Edep_badN_Step0", "E_{dep} [MeF] of CND Neutrons;E_{dep};Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_badN_Step0);
    TH2D *h_beta_Edep_badN_Step0 = new TH2D("Edep_beta_badN", "#beta vs. E_{dep} [MeV] of CND Neutrons;#beta;E_{dep}", 50, 0, 1.1, 50, 0, 100);
    hist_list_2_A.push_back(h_beta_Edep_badN_Step0);

    // Step One (After Beta Cut) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step1 = new TH2D("pnRes_theta_nmiss_Step1", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step1);

    TH1D *h_pmiss_goodN_Step1 = new TH1D("pmiss_goodN_Step1", "p_{miss} Step1;p_{miss};Counts", 25, 0.25, 1.0);
    hist_list_1_A.push_back(h_pmiss_goodN_Step1);
    TH1D *h_ToF_goodN_Step1 = new TH1D("ToF_goodN_Step1", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step1);
    TH1D *h_ToF_badN_Step1 = new TH1D("ToF_badN_Step1", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step1);

    TH1D *h_edep_goodN_Step1 = new TH1D("edep_goodN_Step1", "edep [MeV] of CND Neutrons;edep;Counts", 100, 0, 50);
    hist_list_1_A.push_back(h_edep_goodN_Step1);
    TH1D *h_edep_badN_Step1 = new TH1D("edep_badN_Step1", "edep [MeV] of CND Neutrons;edep;Counts", 100, 0, 50);
    hist_list_1_A.push_back(h_edep_badN_Step1);

    TH1D *h_edep_over_edepCTOT_goodN_Step1 = new TH1D("edep_over_edepCTOT_goodN_Step1", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 100, 0, 5);
    hist_list_1_A.push_back(h_edep_over_edepCTOT_goodN_Step1);
    TH1D *h_edep_over_edepCTOT_badN_Step1 = new TH1D("edep_over_edepCTOT_badN_Step1", "E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts", 100, 0, 5);
    hist_list_1_A.push_back(h_edep_over_edepCTOT_badN_Step1);

    TH1D *h_edep_goodN_withNearbyPos_Step1 = new TH1D("edep_goodN_withNearbyPos_Step1", "edep [MeV] of CND Neutrons;edep;Counts", 100, 0, 50);
    hist_list_1_A.push_back(h_edep_goodN_withNearbyPos_Step1);
    TH1D *h_edep_badN_withNearbyPos_Step1 = new TH1D("edep_badN_withNearbyPos_Step1", "edep [MeV] of CND Neutrons;edep;Counts", 100, 0, 50);
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

    TH2D *h_diff_ToFc_z_Edep_noNear_goodN_Step1 = new TH2D("diff_ToFc_z_Edep_noNear_goodN_Step1", "ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons with no Nearby Tracks;ToF*c-z;E_{dep}", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_noNear_goodN_Step1);
    TH2D *h_diff_ToFc_z_Edep_yesNear_goodN_Step1 = new TH2D("diff_ToFc_z_Edep_yesNear_goodN_Step1", "ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons with no Nearby Tracks;ToF*c-z;E_{dep}", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_yesNear_goodN_Step1);
    TH2D *h_diff_ToFc_z_Edep_noNear_badN_Step1 = new TH2D("diff_ToFc_z_Edep_noNear_badN_Step1", "ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons with no Nearby Tracks;ToF*c-z;E_{dep}", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_noNear_badN_Step1);
    TH2D *h_diff_ToFc_z_Edep_yesNear_badN_Step1 = new TH2D("diff_ToFc_z_Edep_yesNear_badN_Step1", "ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons with no Nearby Tracks;ToF*c-z;E_{dep}", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_yesNear_badN_Step1);

    // Step Two (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step2 = new TH2D("pnRes_theta_nmiss_Step2", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step2);

    TH1D *h_ToF_goodN_Step2 = new TH1D("ToF_goodN_Step2", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step2);
    TH1D *h_ToF_badN_Step2 = new TH1D("ToF_badN_Step2", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step2);

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
    TH2D *h_numberNearby_momN_goodN_Step2 = new TH2D("numberNearby_momN_goodN_Step2", "Number of Nearby Hits vs. p_{n} for CND Neutrons;# Hits;p_{n}", 9, -0.5, 8.5, 50, 0, 1.3);
    hist_list_2_A.push_back(h_numberNearby_momN_goodN_Step2);
    TH1D *h_numberNearby_badN_Step2 = new TH1D("numberNearby_badN_Step2", "Number of Nearby Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    hist_list_1_A.push_back(h_numberNearby_badN_Step2);
    TH2D *h_numberNearby_momN_badN_Step2 = new TH2D("numberNearby_momN_badN_Step2", "Number of Nearby Hits vs. p_{n} for CND Neutrons;# Hits;p_{n}", 9, -0.5, 8.5, 50, 0, 1.3);
    hist_list_2_A.push_back(h_numberNearby_momN_badN_Step2);

    TH1D *h_NearbyEdep_goodN_Step2 = new TH1D("NearbyEdep_goodN_Step2", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep};Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_NearbyEdep_goodN_Step2);
    TH1D *h_NearbyEdep_badN_Step2 = new TH1D("NearbyEdep_badN_Step2", "E_{dep} of Nearby Hits for CND Neutrons;E_{dep};Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_NearbyEdep_badN_Step2);

    TH1D *h_nsector_goodN_Step2 = new TH1D("nsector_goodN_Step2", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    hist_list_1_A.push_back(h_nsector_goodN_Step2);
    TH1D *h_nsector_badN_Step2 = new TH1D("nsector_badN_Step2", "Neutron Sector for CND Neutrons;Neutron Sector;Counts", 24, 0.5, 24.5);
    hist_list_1_A.push_back(h_nsector_badN_Step2);

    // Step Three (After applying Phi Diff Charge Track cut) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step3 = new TH2D("pnRes_theta_nmiss_Step3", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step3);

    TH1D *h_ToF_goodN_Step3 = new TH1D("ToF_goodN_Step3", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step3);
    TH1D *h_ToF_badN_Step3 = new TH1D("ToF_badN_Step3", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step3);

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
    TH2D *h_numberCTOF_momN_goodN_Step3 = new TH2D("numberCTOF_momN_goodN_Step3", "Number of CTOF Hits vs. p_{n} for CND Neutrons;# Hits;p_{n}", 9, -0.5, 8.5, 50, 0, 1.3);
    hist_list_2_A.push_back(h_numberCTOF_momN_goodN_Step3);
    TH1D *h_numberCTOF_badN_Step3 = new TH1D("numberCTOF_badN_Step3", "Number of CTOF Hits for CND Neutrons;# Hits;Counts", 9, -0.5, 8.5);
    hist_list_1_A.push_back(h_numberCTOF_badN_Step3);
    TH2D *h_numberCTOF_momN_badN_Step3 = new TH2D("numberCTOF_momN_badN_Step3", "Number of CTOF Hits vs. p_{n} for CND Neutrons;# Hits;p_{n}", 9, -0.5, 8.5, 50, 0, 1.3);
    hist_list_2_A.push_back(h_numberCTOF_momN_badN_Step3);

    // Step Four (After applying Phi Diff CND hit cut) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step4 = new TH2D("pnRes_theta_nmiss_Step4", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step4);

    TH1D *h_ToF_goodN_Step4 = new TH1D("ToF_goodN_Step4", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step4);
    TH1D *h_ToF_badN_Step4 = new TH1D("ToF_badN_Step4", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step4);

    // Step Five (After event selection cuts) (Andrew)
    // ======================================================================================================================================================================
    TH2D *h_pnRes_theta_nmiss_Step5 = new TH2D("pnRes_theta_nmiss_Step5", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step5);

    TH1D *h_ToF_goodN_Step5 = new TH1D("ToF_goodN_Step5", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step5);
    TH1D *h_ToF_badN_Step5 = new TH1D("ToF_badN_Step5", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step5);
    TH1D *h_pmiss_goodN_Step5 = new TH1D("pmiss_goodN_Step5", "p_{miss} good N Step5;p_{miss};Counts", 25, 0.25, 1.0);
    hist_list_1_A.push_back(h_pmiss_goodN_Step5);
    TH1D *h_pmiss_allN_Step5 = new TH1D("pmiss_allN_Step5", "p_{miss} all N Step5;p_{miss};Counts", 25, 0.25, 1.0);
    hist_list_1_A.push_back(h_pmiss_allN_Step5);

    TH1D *h_Edep_infront_goodN_Step5 = new TH1D("Edep_infront_goodN_Step5", "E_{dep} [MeV] of Hit infront CND Neutrons;E_{dep};Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_infront_goodN_Step5);
    TH1D *h_Edep_behind_goodN_Step5 = new TH1D("Edep_behind_goodN_Step5", "E_{dep} [MeV] of Hit behind CND Neutrons;E_{dep};Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_behind_goodN_Step5);

    TH1D *h_Edep_infront_badN_Step5 = new TH1D("Edep_infront_badN_Step5", "E_{dep} [MeV] of Hit infront CND Neutrons;E_{dep};Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_infront_badN_Step5);
    TH1D *h_Edep_behind_badN_Step5 = new TH1D("Edep_behind_badN_Step5", "E_{dep} [MeV] of Hit behind CND Neutrons;E_{dep};Counts", 50, 0, 100);
    hist_list_1_A.push_back(h_Edep_behind_badN_Step5);

    TH2D *h_diff_ToFc_z_Edep_goodN_Step5 = new TH2D("diff_ToFc_z_Edep_goodN_Step5", "ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_goodN_Step5);
    TH2D *h_diff_ToFc_z_Edep_badN_Step5 = new TH2D("diff_ToFc_z_Edep_badN_Step5", "ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}", 50, 0, 300, 50, 0, 100);
    hist_list_2_A.push_back(h_diff_ToFc_z_Edep_badN_Step5);

    TH2D *h_diff_ToFc_z_Edep_goodN_Step5_layer[3];

    for (int k = 0; k < 3; k++)
    {
        sprintf(temp_name_A, "diff_ToFc_z_goodN_Step5_layer_%d", k + 1);
        sprintf(temp_title_A, "ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons (Layer Difference = %d);ToF*c-z;E_{dep}", k + 1);
        h_diff_ToFc_z_Edep_goodN_Step5_layer[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        hist_list_2_A.push_back(h_diff_ToFc_z_Edep_goodN_Step5_layer[k]);
    }

    TH2D *h_diff_ToFc_z_Edep_badN_Step5_layer[3];

    for (int k = 0; k < 3; k++)
    {
        sprintf(temp_name_A, "diff_ToFc_z_badN_Step5_layer_%d", k + 1);
        sprintf(temp_title_A, "ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons (Layer Difference = %d);ToF*c-z;E_{dep}", k + 1);
        h_diff_ToFc_z_Edep_badN_Step5_layer[k] = new TH2D(temp_name_A, temp_title_A, 50, 0, 300, 50, 0, 100);
        hist_list_2_A.push_back(h_diff_ToFc_z_Edep_badN_Step5_layer[k]);
    }

    TH1D *h_phidiff_en_goodN_Step5 = new TH1D("phidiff_en_goodN_Step5", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}|;Counts", 90, 0, 180);
    hist_list_1_A.push_back(h_phidiff_en_goodN_Step5);

    TH1D *h_phidiff_en_badN_Step5 = new TH1D("phidiff_en_badN_Step5", "|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}|;Counts", 90, 0, 180);
    hist_list_1_A.push_back(h_phidiff_en_badN_Step5);

    TH1D *h_TP_goodN_Step5 = new TH1D("TP_goodN_Step5", "ToF/path [ns/m] of CND Neutrons;ToF/path;Counts", 150, 0, 50);
    hist_list_1_A.push_back(h_TP_goodN_Step5);
    TH1D *h_TP_badN_Step5 = new TH1D("TP_badN_Step5", "ToF/path [ns/m] of CND Neutrons;ToF/path;Counts", 150, 0, 50);
    hist_list_1_A.push_back(h_TP_badN_Step5);

    TH1D *h_Z_goodN_Step5 = new TH1D("Z_goodN_Step5", "Z [cm] of CND Neutrons;Z;Counts", 100, -60, 60);
    hist_list_1_A.push_back(h_Z_goodN_Step5);
    TH1D *h_Z_badN_Step5 = new TH1D("Z_badN_Step5", "Z [cm] of CND Neutrons;Z;Counts", 100, -60, 60);
    hist_list_1_A.push_back(h_Z_badN_Step5);

    TH2D *h_beta_Edep_goodN_Step5 = new TH2D("Edep_beta_goodN", "#beta vs. E_{dep} [MeV] of CND Neutrons;#beta;E_{dep}", 50, 0, 1.1, 50, 0, 100);
    hist_list_2_A.push_back(h_beta_Edep_goodN_Step5);
    TH2D *h_beta_Edep_badN_Step5 = new TH2D("Edep_beta_badN", "#beta vs. E_{dep} [MeV] of CND Neutrons;#beta;E_{dep}", 50, 0, 1.1, 50, 0, 100);
    hist_list_2_A.push_back(h_beta_Edep_badN_Step5);

    TH2D *h_ToF_Edep_goodN_Step5 = new TH2D("ToF_Edep_goodN_Step5", "ToF [ns] vs. E_{dep} [MeV] of CND Neutrons;ToF;E_{dep} MeV", 100, 0, 20, 50, 0, 100);
    hist_list_2_A.push_back(h_ToF_Edep_goodN_Step5);
    TH2D *h_ToF_Edep_badN_Step5 = new TH2D("ToF_Edep_badN_Step5", "ToF [ns] vs. E_{dep} [MeV] of CND Neutrons;ToF;E_{dep} MeV", 100, 0, 20, 50, 0, 100);
    hist_list_2_A.push_back(h_ToF_Edep_badN_Step5);

    TH2D *h_TP_Edep_goodN_Step5 = new TH2D("TP_Edep_goodN_Step5", "TP [ns/m] vs. E_{dep} [MeV] of CND Neutrons;TP;E_{dep} MeV", 150, 0, 50, 50, 0, 100);
    hist_list_2_A.push_back(h_TP_Edep_goodN_Step5);
    TH2D *h_TP_Edep_badN_Step5 = new TH2D("TP_Edep_badN_Step5", "TP [ns/m] vs. E_{dep} [MeV] of CND Neutrons;TP;E_{dep} MeV", 150, 0, 50, 50, 0, 100);
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

#pragma endregion /* Andrew's histograms */

    // ======================================================================================================================================================================
    // Chain loop
    // ======================================================================================================================================================================

#pragma region /* Chain loop */

    int counter_A = 0; /* From Andrew */

    while (chain.Next())
    {
        // Display completed (from Andrew)
        counter_A++;
        if ((counter_A % 1000000) == 0)
        {
            cerr << "\n"
                 << counter_A / 1000000 << " million completed";
        }
        if ((counter_A % 100000) == 0)
        {
            cerr << ".";
        }

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

            // Particle PID
            // ===================================================================================================================================================================

            clasAna->Run(c12);

            auto elec = clasAna->getByPid(11);
            auto prot = clasAna->getByPid(2212);
            auto neut = clasAna->getByPid(2112);

            auto allParticles = c12->getDetParticles();

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

            for (int i = 0; i < allParticles.size(); i++)
            {
                int pid = allParticles[i]->par()->getPid();

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
            double starttime = c12->event()->getStartTime();

#pragma region /* Electrons */

            //////////////////////////
            /////    ELECTRONS   /////
            //////////////////////////
            TVector3 pe(0., 0., 0.);

            double pe_x = elec[0]->par()->getPx();
            double pe_y = elec[0]->par()->getPy();
            double pe_z = elec[0]->par()->getPz();

            pe.SetXYZ(pe_x, pe_y, pe_z);

            double vze = elec[0]->par()->getVz();

            TVector3 pb(0, 0, Ebeam);
            TVector3 pq = pb - pe; // 3-momentum transfer

            double nu = Ebeam - pe.Mag();       // Energy transfer
            double QSq = pq.Mag2() - (nu * nu); // 4-momentum transfer squared
            double xB = QSq / (2 * mN * nu);    // x Bjorken

#pragma endregion /* Electrons */

#pragma region /* Protons */

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

#pragma endregion /* Protons */

#pragma region /* Missing momentum */

            //////////////////////////
            //  MISSING MOMENTUM    //
            //////////////////////////

            // missing momentum, energy, mass
            TVector3 pmiss = pq - pp; // TODO: checkout difference from Andrew - he uses leading SRC proton here!

            momentum = pmiss.Mag();

            double Ep = sqrt(mN * mN + pp.Mag2());
            double Emiss = Ebeam + mD - pe.Mag() - Ep;
            double mmiss = sqrt((Emiss * Emiss) - pmiss.Mag2());

#pragma endregion /* Missing momentum */

#pragma region /* Neutrons */
            //////////////////////////
            ////     NEUTRONS    /////
            //////////////////////////

            // LOOP OVER NEUTRONS
            h_nsize->Fill(neut.size());

            for (int i = 0; i < neut.size(); i++)
            {
                // GET NEUTRON INFORMATION

                // get neutron momentums
                double pn_x = neut[i]->par()->getPx();
                double pn_y = neut[i]->par()->getPy();
                double pn_z = neut[i]->par()->getPz();

                TVector3 pn;
                pn.SetXYZ(pn_x, pn_y, pn_z);

                double dpp = (pmiss.Mag() - pn.Mag()) / pmiss.Mag();

                // figure out what layer the hit is in
                is_CND1 = (neut[i]->sci(CND1)->getLayer() == 1);
                is_CND2 = (neut[i]->sci(CND2)->getLayer() == 2);
                is_CND3 = (neut[i]->sci(CND3)->getLayer() == 3);
                is_CTOF = neut[i]->sci(CTOF)->getDetector() == 4;

                // put REC::Scintillator information
                double time;

                int status = 0;

                double beta = neut[i]->par()->getBeta();

                if (is_CND1)
                {
                    time = neut[i]->sci(CND1)->getTime() - starttime;
                    status = status + neut[i]->sci(CND1)->getStatus();
                }

                if (is_CND3)
                {
                    time = neut[i]->sci(CND3)->getTime() - starttime;
                    status = status + neut[i]->sci(CND3)->getStatus();
                }

                if (is_CND2)
                {
                    time = neut[i]->sci(CND2)->getTime() - starttime;
                    status = status + neut[i]->sci(CND2)->getStatus();
                }

                // PROBLEM: this gives preference to 2nd-layer hits
                // TODO: recheck this!
                if (is_CTOF)
                {
                    time = neut[i]->sci(CTOF)->getTime() - starttime;
                }

                double cos0 = pmiss.Dot(pn) / (pmiss.Mag() * pn.Mag());

                if (status != 0) // Cutting out neutrons suspected to have double-hit CND hits
                {
                    continue;
                }

                // GET ML FEATURES FOR THIS NEUTRON
                Struct ninfo = getFeatures(neut, allParticles, i);
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

                if (pn_x == 0 || pn_y == 0 || pn_z == 0)
                {
                    continue;
                }

                if (time > 10) // TODO: why this cut?
                {
                    continue;
                }

                h_pvsp->Fill(pmiss.Mag(), pn.Mag());

                // select neutrons in momentum and angle accepted by CND
                double n_theta = pn.Theta() * 180. / M_PI;

                if (pn.Mag() < 0.25 || pn.Mag() > 1)
                {
                    continue;
                }

                if (n_theta < 45 || n_theta > 140)
                {
                    continue;
                }

                h_mmiss_pn->Fill(pn.Mag(), mmiss);
                h_mmiss_pmiss->Fill(pmiss.Mag(), mmiss);
                h_mmiss_xb->Fill(xB, mmiss);

                h_mmiss->Fill(mmiss);

                if (mmiss > 1.) // Missing mass cut
                {
                    continue;
                }

                h_pmiss_thetamiss->Fill(pmiss.Theta() * 180. / M_PI, pmiss.Mag());
                h_thetapn_pp->Fill(pp.Mag(), pp.Angle(pn) * 180. / M_PI);

                if (pmiss.Mag() < 0.25 || pmiss.Mag() > 1.) // Missing momentum cut
                {
                    continue;
                }

                if (pmiss.Theta() * 180. / M_PI < 45 || pmiss.Theta() * 180. / M_PI > 140) // Missing momentum theta cut
                {
                    continue;
                }

                h_dpp_edep->Fill(energy, dpp);

                if (energy < 5)
                {
                    continue;
                }

                // FILL HISTOS FOR NEUTRON CANDIDATES
                h_nangles->Fill(pn.Phi() * 180. / M_PI, n_theta);
                h_energy->Fill(energy);
                h_Edep_beta->Fill(neut[i]->getBeta(), energy);

                h_cos0->Fill(pmiss.Dot(pn) / (pmiss.Mag() * pn.Mag()));
                h_pxminuspx->Fill(pn_x - pmiss.X());
                h_pyminuspy->Fill(pn_y - pmiss.Y());
                h_pzminuspz->Fill(pn_z - pmiss.Z());
                h_pminusp->Fill(pn.Mag() - pmiss.Mag());

                h_dpp->Fill(pmiss.Mag(), (pmiss.Mag() - pn.Mag()) / pmiss.Mag());
                h_theta_beta->Fill(beta, n_theta);
                h_p_theta->Fill(n_theta, pn.Mag());
                h_p_all->Fill(pmiss.Mag());
                h_anglediff->Fill(angle_diff);

                h_compare->Fill((pmiss.Mag() - pn.Mag()) / pmiss.Mag(), pn.Angle(pmiss) * 180. / M_PI);

                if ((fabs(pmiss.Mag() - pn.Mag()) / pmiss.Mag()) > 0.2) // Relative momentum difference cut
                // if ((abs(pmiss.Mag() - pn.Mag()) / pmiss.Mag()) > 0.2) // Erin's original
                {
                    continue;
                }

                h_thetapn_dpp->Fill((pmiss.Mag() - pn.Mag()) / pmiss.Mag(), pn.Angle(pp) * 180. / M_PI);
                h_thetapn_dpp1->Fill((pmiss.Mag() - pn.Mag()) / pmiss.Mag(), pn.Angle(pp) * 180. / M_PI);

                if (pn.Angle(pmiss) * 180. / M_PI > 20) // pn close to pmiss cut
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

                bool good_N = pn.Angle(pmiss) * 180. / M_PI < 20 &&
                              fabs((pmiss.Mag() - pn.Mag()) / pmiss.Mag()) < 0.2 &&
                              //   abs((pmiss.Mag() - pn.Mag()) / pmiss.Mag()) < 0.2 && // Erin's original
                              cnd_energy < 1000 &&
                              pp.Angle(pn) * 180. / M_PI > 60 &&
                              (pmiss.Mag() > 0.25 &&
                               pmiss.Mag() < 1.) &&
                              (pmiss.Theta() * 180. / M_PI > 45 &&
                               pmiss.Theta() * 180. / M_PI < 140);

                bool bad_N = (pn.Angle(pmiss) * 180. / M_PI > 50 ||
                              fabs((pmiss.Mag() - pn.Mag()) / pmiss.Mag()) > 0.6) && // Erin's original
                             cnd_energy < 1000;                                      // && (pp.Angle(pn)*180./M_PI<60);

                bool keep_this_one = keep_good ? good_N : bad_N;

                if (keep_this_one)
                {
                    // all neutrons - print features
                    outtxt << pmiss.Mag() << ' ';
                    cout << pmiss.Mag() << ' ';
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
                    h_nangles2->Fill(pn.Phi() * 180. / M_PI, n_theta);
                    h_cos02->Fill(pmiss.Dot(pn) / (pmiss.Mag() * pn.Mag()));
                    h_pxminuspx2->Fill(pn_x - pmiss.X());
                    h_pyminuspy2->Fill(pn_y - pmiss.Y());
                    h_pzminuspz2->Fill(pn_z - pmiss.Z());
                    h_pminusp2->Fill(pn.Mag() - pmiss.Mag());
                    h_pvsp2->Fill(pmiss.Mag(), pn.Mag());
                    h_dpp2->Fill(pmiss.Mag(), (pmiss.Mag() - pn.Mag()) / pmiss.Mag());
                    h_mmiss2->Fill(mmiss);
                    h_mmiss_pn2->Fill(pn.Mag(), mmiss);
                    h_energy2->Fill(energy);
                    h_theta_beta2->Fill(beta, n_theta);
                    h_p_theta2->Fill(n_theta, pn.Mag());
                    h_pmiss_thetamiss2->Fill(pmiss.Theta() * 180. / M_PI, pmiss.Mag());
                    h_thetapn_pp2->Fill(pp.Mag(), pp.Angle(pn) * 180. / M_PI);
                    h_tof2->Fill(time);
                    h_compare2->Fill((pmiss.Mag() - pn.Mag()) / pmiss.Mag(), pn.Angle(pmiss) * 180. / M_PI);
                    h_Edep_beta2->Fill(neut[i]->getBeta(), energy);
                    h_p_cut->Fill(pmiss.Mag());
                    h_anglediff2->Fill(angle_diff);
                    h_thetapn_dpp2->Fill((pmiss.Mag() - pn.Mag()) / pmiss.Mag(), pn.Angle(pp) * 180. / M_PI);

                    h_ptheta_pred->Fill(pmiss.Theta() * 180. / M_PI, pmiss.Mag());
                    h_ptheta->Fill(pn.Theta() * 180. / M_PI, pn.Mag());

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

#pragma endregion /* Neutrons */

            // chain.WriteEvent();
            counter++;
                }

#pragma endregion /* Erin's features */

        // ==================================================================================================================================================================
        // Andrew's manual work
        // ==================================================================================================================================================================

#pragma region /* Andrew's manual work */

        if (Run_Andrews_work)
        {

            auto electrons = clasAna->getByPid(11); /* From Erin's code */
            auto protons = clasAna->getByPid(2212); /* From Erin's code */

            auto allParticles = c12->getDetParticles();

            /*
            // Andrew's original - commented out!
            auto allParticles = c12->getDetParticles();
            auto electrons = c12->getByID(11);
            */

            double weight = 1;

            if (isMC)
            {
                weight = c12->mcevent()->getWeight();
            }

            TVector3 p_b(0, 0, Ebeam);

#pragma region /* Electrons */

            if (electrons.size() != 1)
            {
                continue;
            }

            TVector3 p_e;
            p_e.SetMagThetaPhi(electrons[0]->getP(), electrons[0]->getTheta(), electrons[0]->getPhi());

            double EoP_e = (electrons[0]->cal(PCAL)->getEnergy() + electrons[0]->cal(ECIN)->getEnergy() + electrons[0]->cal(ECOUT)->getEnergy()) / p_e.Mag();
            int nphe = electrons[0]->che(HTCC)->getNphe();
            double vtz_e = electrons[0]->par()->getVz();

            /*
            // Andrew's original - commented out!
            if (!myCut.electroncut(c12))
            {
                continue;
            }
            */

            int esector = electrons[0]->getSector();

            /////////////////////////////////////
            // Electron Kinematics
            /////////////////////////////////////
            TVector3 p_q = p_b - p_e; // 3-momentum transfer (same as Erin's code)
            double theta_q = p_q.Theta() * 180 / M_PI;
            double nu = Ebeam - p_e.Mag();       // Energy transfer (same as Erin's code)
            double QSq = p_q.Mag2() - (nu * nu); // 4-momentum transfer squared (same as Erin's code)
            double xB = QSq / (2 * mN * nu);     // x Bjorken (same as Erin's code)
            double WSq = (mN * mN) - QSq + (2 * nu * mN);
            double theta_e = p_e.Theta() * 180 / M_PI;

#pragma endregion /* Electrons */

#pragma region /* Protons */

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

            /*
            // Andrew's original - commented out!
            // Lead Proton
            int num_L = 0;
            int index_L = -1;

            for (int j = 0; j < allParticles.size(); j++)
            {
                if ((LeadFDProton_Cut(c12, Ebeam, j)) || (LeadCDProton_Cut(c12, Ebeam, j)))
                {
                    num_L++;
                    index_L = j;
                }
            }

            if (num_L != 1)
            {
                continue;
            }

            bool LeadCD = LeadCDProton_Cut(c12, Ebeam, index_L);
            bool LeadFD = LeadFDProton_Cut(c12, Ebeam, index_L);

            if (LeadCD && LeadFD)
            {
                cout << "Problem!\n";
            }
            */

            // Andrew's original (p_L) was replaced by p_p from Erin's code
            // TVector3 p_L; // Andrew's original
            // p_L.SetMagThetaPhi(allParticles[index_L]->getP(), allParticles[index_L]->getTheta(), allParticles[index_L]->getPhi()); // Andrew's original

#pragma endregion

#pragma region /* Missing momentum */

            TVector3 p_miss = p_q - p_p;
            double Ep = sqrt(mN * mN + pp.Mag2());
            double Emiss = Ebeam + mD - pe.Mag() - Ep;
            double mmiss = sqrt((Emiss * Emiss) - p_miss.Mag2());
            // double mmiss = get_mmiss(p_b, p_e, p_L); // Andrew's original - commented out!

            if (p_miss.Theta() * 180 / M_PI < 40)
            {
                continue;
            }

            if (p_miss.Theta() * 180 / M_PI > 135)
            {
                continue;
            }

            if (p_miss.Mag() < 0.2)
            {
                continue;
            }

            if (p_miss.Mag() > 1.25)
            {
                continue;
            }

            if (mmiss < 0.7)
            {
                continue;
            }

            if (mmiss > 1.2)
            {
                continue;
            }

            bool match = false;

            //////////////////////////////////////////////////
            // For after checking the hipo banks
            //////////////////////////////////////////////////
            int num_Charge = 0;

            for (int j = 0; j < allParticles.size(); j++)
            {
                if (j == 0)
                {
                    continue;
                }

                if (j == index_L)
                {
                    continue;
                }

                // if(j==index_Rp1){continue;}
                if (allParticles[j]->par()->getCharge() == 0)
                {
                    continue;
                }

                num_Charge++;
            }

            if (num_Charge > 0)
            {
                continue;
            }

            if (LeadFD)
            {
                h_xB_mmiss_epFD->Fill(xB, mmiss, weight);
            }
            else if (LeadCD)
            {
                h_xB_mmiss_epCD->Fill(xB, mmiss, weight);
            }

            h_pmiss_ep->Fill(p_miss.Mag(), weight);

            if (mmiss > 1.05) // Missing mass cut
            {
                continue;
            }

            if (LeadCD && (xB < 1.1)) // Cutting out CD leading protons with xB < 1.1
            {
                continue;
            }

            if (LeadFD) // Cutting out FD leading protons
            {
                continue;
            }

            // if(LeadFD && (xB<0.8)){continue;}

#pragma endregion /* Missing momentum */

            /////////////////////////////////////
            // Lead Neutron Checks
            /////////////////////////////////////
            for (int j = 0; j < allParticles.size(); j++)
            {
                if (allParticles[j]->par()->getCharge() != 0)
                {
                    continue;
                }

                bool CT = (allParticles[j]->sci(clas12::CTOF)->getDetector() == 4);
                bool C1 = (allParticles[j]->sci(clas12::CND1)->getDetector() == 3);
                bool C2 = (allParticles[j]->sci(clas12::CND2)->getDetector() == 3);
                bool C3 = (allParticles[j]->sci(clas12::CND3)->getDetector() == 3);

                if (!(C1 || C2 || C3))
                {
                    continue;
                }
                if (allParticles[j]->getTheta() * 180 / M_PI > 160)
                {
                    continue;
                }
                double theta = allParticles[j]->getTheta() * 180 / M_PI;
                double beta = allParticles[j]->par()->getBeta();
                double gamma = 1 / sqrt(1 - (beta * beta));
                double mom = gamma * beta * mN;
                double ToF = allParticles[j]->getTime() - c12->event()->getStartTime();

                int detINTlayer = C1 ? 1 : C2 ? 2
                                              : 3;
                auto detlayer = C1 ? CND1 : C2 ? CND2
                                               : CND3;
                double edep = allParticles[j]->sci(CND1)->getEnergy() + allParticles[j]->sci(CND2)->getEnergy() + allParticles[j]->sci(CND3)->getEnergy();
                double edep_CTOF = allParticles[j]->sci(CTOF)->getEnergy();
                double edep_single = allParticles[j]->sci(detlayer)->getEnergy();

                double nvtx_x = allParticles[j]->par()->getVx();
                double nvtx_y = allParticles[j]->par()->getVy();
                double nvtx_z = allParticles[j]->par()->getVz();
                TVector3 v_nvtx(nvtx_x, nvtx_y, nvtx_z);

                TVector3 v_hit;
                v_hit.SetXYZ(allParticles[j]->sci(detlayer)->getX(), allParticles[j]->sci(detlayer)->getY(), allParticles[j]->sci(detlayer)->getZ());

                TVector3 v_path = v_hit - v_nvtx;
                TVector3 v_n;
                v_n.SetMagThetaPhi(mom, v_path.Theta(), v_path.Phi());

                double path = v_path.Mag() / 100;
                double theta_nmiss = v_n.Angle(p_miss) * 180 / M_PI;
                double dm_nmiss = (p_miss.Mag() - v_n.Mag()) / p_miss.Mag();
                int nSector = allParticles[j]->sci(detlayer)->getSector();

                // Check to see if there is a good neutron
                bool isGN = false;

                if ((theta_nmiss < 40) && (dm_nmiss > -0.5) && (dm_nmiss < 0.5))
                {
                    isGN = true;
                }

                //////////////////////////////////////////////
                // Step Zero
                //////////////////////////////////////////////
                if (beta - (path * 100) / (ToF * c) < -0.01)
                {
                    continue;
                }

                if (beta - (path * 100) / (ToF * c) > 0.01)
                {
                    continue;
                }

                if (v_hit.Z() > 45)
                {
                    continue;
                }

                if (v_hit.Z() < -40)
                {
                    continue;
                }

                if (ToF < 0)
                {
                    continue;
                }

                if (ToF > 20)
                {
                    continue;
                }

                if (LeadFD)
                {
                    h_xB_mmiss_epnFD->Fill(xB, mmiss, weight);
                }

                else if (LeadCD)
                {
                    h_xB_mmiss_epnCD->Fill(xB, mmiss, weight);
                }

                h_pnRes_theta_nmiss_Step0->Fill(dm_nmiss, theta_nmiss, weight);

                if (isGN)
                {
                    if (LeadFD)
                    {
                        h_xB_mmiss_epngoodFD->Fill(xB, mmiss, weight);
                    }
                    else if (LeadCD)
                    {
                        h_xB_mmiss_epngoodCD->Fill(xB, mmiss, weight);
                    }

                    h_ToF_goodN_Step0->Fill(ToF, weight);
                    h_beta_goodN_Step0->Fill(beta, weight);
                    h_Edep_goodN_Step0->Fill(edep, weight);
                    h_beta_Edep_goodN_Step0->Fill(beta, edep, weight);
                }
                else
                {
                    h_ToF_badN_Step0->Fill(ToF, weight);
                    h_beta_badN_Step0->Fill(beta, weight);
                    h_Edep_badN_Step0->Fill(edep, weight);
                    h_beta_Edep_badN_Step0->Fill(beta, edep, weight);
                }

                //////////////////////////////////////////////
                // Step One
                //////////////////////////////////////////////
                if (beta > 0.8) // Beta cut
                {
                    continue;
                }

                if (edep < 5) // Dep. energy cut
                {
                    continue;
                }

                h_pnRes_theta_nmiss_Step1->Fill(dm_nmiss, theta_nmiss, weight);

                if (isGN)
                {
                    h_ToF_goodN_Step1->Fill(ToF, weight);
                    h_pmiss_goodN_Step1->Fill(p_miss.Mag(), weight);
                }
                else
                {
                    h_ToF_badN_Step1->Fill(ToF, weight);
                }

                bool CNDVeto = false;

                if (ToF * c - v_hit.Z() < 70) // TODO: why this cut?
                {

                    if (isGN)
                    {
                        h_edep_goodN_Step1->Fill(edep, weight);
                    }
                    else
                    {
                        h_edep_badN_Step1->Fill(edep, weight);
                    }

                    for (int k = 0; k < allParticles.size(); k++)
                    {
                        if (k == 0)
                        {
                            continue;
                        }

                        if (k == j)
                        {
                            continue;
                        }

                        if (allParticles[k]->par()->getCharge() <= 0)
                        {
                            continue;
                        }

                        if (allParticles[k]->sci(CTOF)->getDetector() == 0)
                        {
                            continue;
                        }

                        // TODO: what is this?
                        int vetoSectorbyLayer[4] = {(allParticles[k]->sci(CTOF)->getComponent() + 1) / 2, allParticles[k]->sci(CND1)->getSector(), allParticles[k]->sci(CND2)->getSector(), allParticles[k]->sci(CND3)->getSector()};

                        TVector3 p_C;
                        p_C.SetMagThetaPhi(allParticles[k]->getP(), allParticles[k]->getTheta(), allParticles[k]->getPhi());

                        double edep_pos = allParticles[k]->sci(clas12::CTOF)->getEnergy();

                        for (int k = 0; k < 4; k++)
                        {
                            if (vetoSectorbyLayer[k] == 0)
                            {
                                continue;
                            }

                            int sdiff = nSector - vetoSectorbyLayer[k];

                            if (sdiff <= -12)
                            {
                                sdiff += 24;
                            }
                            else if (sdiff > 12)
                            {
                                sdiff -= 24;
                            }

                            int ldiff = detINTlayer - k;

                            if (isGN)
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
                        }

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
                    }

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
                if (CNDVeto)
                {
                    continue;
                }

                h_pnRes_theta_nmiss_Step2->Fill(dm_nmiss, theta_nmiss, weight);

                if (isGN)
                {
                    h_ToF_goodN_Step2->Fill(ToF, weight);
                }
                else
                {
                    h_ToF_badN_Step2->Fill(ToF, weight);
                }

                for (int k = 0; k < allParticles.size(); k++)
                {
                    if (k == 0)
                    {
                        continue;
                    }

                    if (k == j)
                    {
                        continue;
                    }

                    if (allParticles[k]->par()->getCharge() <= 0)
                    {
                        continue;
                    }

                    if (allParticles[k]->sci(CTOF)->getDetector() == 0)
                    {
                        continue;
                    }

                    int vetoSectorbyLayer[4] = {(allParticles[k]->sci(CTOF)->getComponent() + 1) / 2, allParticles[k]->sci(CND1)->getSector(), allParticles[k]->sci(CND2)->getSector(), allParticles[k]->sci(CND3)->getSector()};

                    for (int k = 0; k < 4; k++)
                    {
                        if (vetoSectorbyLayer[k] == 0)
                        {
                            continue;
                        }

                        int sdiff = nSector - vetoSectorbyLayer[k];

                        if (sdiff <= -12)
                        {
                            sdiff += 24;
                        }
                        else if (sdiff > 12)
                        {
                            sdiff -= 24;
                        }

                        int ldiff = detINTlayer - k;

                        if (isGN)
                        {
                            h_sdiff_pos_goodN_Step2_layer[ldiff + 3]->Fill(sdiff, weight);
                        }
                        else
                        {
                            h_sdiff_pos_badN_Step2_layer[ldiff + 3]->Fill(sdiff, weight);
                        }
                    }
                }
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
                }
                else{
                  h_ToF_badN_Step3->Fill(ToF,weight);
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
                    //if(isGN){h_NearbyEdep_goodN_Step2->Fill(hit_energy,weight);}
                    //else{h_NearbyEdep_badN_Step2->Fill(hit_energy,weight);}
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
                }
                else{
                  h_ToF_badN_Step4->Fill(ToF,weight);
                }

                ///////////////////
                h_pnRes_theta_nmiss_Step5->Fill(dm_nmiss,theta_nmiss,weight);
                h_pmiss_allN_Step5->Fill(p_miss.Mag(),weight);
                if(isGN){
                  h_ToF_goodN_Step5->Fill(ToF,weight);
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
            }
        }

#pragma endregion /* Andrew's manual work */

    } // closes event loop

#pragma endregion /* Chain loop */

    // ======================================================================================================================================================================
    // Andrew's wrap up
    // ======================================================================================================================================================================

#pragma region /* Andrew's wrap up */

    /////////////////////////////////////////////////////
    // Now create the output PDFs
    /////////////////////////////////////////////////////

    int pixelx = 1980;
    int pixely = 1530;

    TCanvas *myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);
    TCanvas *myText = new TCanvas("myText", "myText", pixelx, pixely);

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

    /*

         std::cout << "\nAndrew's wrap up 1\n\n";
        // std::cout << "\hist_list_1_A.size() = " << hist_list_1_A.size() << endl;

        // TH1D *hhhhh123 = new TH1D("pmiss_ep", "p_{miss} ep;p_{miss}", 25, 0.25, 1.0);
        // // TH1D *hhhhh123 = new TH1D("h_pmiss_ep_11", "p_{miss} ep;p_{miss}", 25, 0.25, 1.0);
        // // TH1D *hhhhh123 = new TH1D("hhhhh123", "p_{miss} ep;p_{miss}", 25, 0.25, 1.0);
        // // TH1D *hhhhh123 = new TH1D("hhhhh123", "h title;var", 25, 0.25, 1.0);

        // h_pmiss_ep->Draw();
        hhhhh123->Draw();
        // h_pmiss_ep_1->Draw();
        // myCanvas->Print(fileName, "pdf");
        // myCanvas->Clear();

        std::cout << "\nAndrew's wrap up 1a\n\n";

     */

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        myCanvas->cd(1);
        // hist_list_1_A[i]->GetXaxis()->CenterTitle();
        // hist_list_1_A[i]->GetXaxis()->SetTitleSize(0.06);
        // hist_list_1_A[i]->GetXaxis()->SetLabelSize(0.0425);
        // hist_list_1_A[i]->GetYaxis()->CenterTitle();
        // hist_list_1_A[i]->GetYaxis()->SetTitleSize(0.06);
        // hist_list_1_A[i]->GetYaxis()->SetLabelSize(0.0425);
        hist_list_1_A[i]->Draw();
        myCanvas->Print(fileName, "pdf");
        myCanvas->Clear();
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        myCanvas->cd(1);
        // hist_list_2_A[i]->GetXaxis()->CenterTitle();
        // hist_list_2_A[i]->GetXaxis()->SetTitleSize(0.06);
        // hist_list_2_A[i]->GetXaxis()->SetLabelSize(0.0425);
        // hist_list_2_A[i]->GetYaxis()->CenterTitle();
        // hist_list_2_A[i]->GetYaxis()->SetTitleSize(0.06);
        // hist_list_2_A[i]->GetYaxis()->SetLabelSize(0.0425);
        hist_list_2_A[i]->Draw("colz");
        myCanvas->Print(fileName, "pdf");
        myCanvas->Clear();
    }

    sprintf(fileName, "%s]", pdfFile);
    myCanvas->Print(fileName, "pdf");

#pragma endregion

    // ======================================================================================================================================================================
    // Erin's wrap up
    // ======================================================================================================================================================================

#pragma region /* Erin's wrap up */

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

    std::cout << '\n'
              << counter << " events counted!\n\n";

    // wrap it up
    outtxt.close();
    ntree->Write();
    f->Close();

#pragma endregion

    return 0;

} // closes main function

#pragma endregion /* Erin main function */