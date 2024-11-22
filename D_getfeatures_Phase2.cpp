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
#include "src/classes/clas12ana/clas12ana.cpp"

using namespace std;
using namespace clas12;

// ===========================================================================================================================================================================
// Erins's main function
// ===========================================================================================================================================================================

int D_getfeatures_Phase2(double Ebeam, bool keep_good, string output_root, string output_txt, string input_hipo,  /* Erin's arguments*/
                         string pdfFile /* Andrew's arguments*/)
// int main(int argc, char **argv)
{
    cout << endl;

    // ======================================================================================================================================================================
    // initial setup
    // ======================================================================================================================================================================

    // arg 1: beam energy

    // arg 2: keep good

    // args 3-4: output file names
    TFile *f = new TFile(output_root.c_str(), "RECREATE");
    TTree *ntree = new TTree("T", "NeutronTree");
    std::ofstream outtxt(output_txt);

    // arg 5+: input hipo file
    clas12root::HipoChain chain;
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

    // ======================================================================================================================================================================
    // Erin's histograms
    // ======================================================================================================================================================================

    // Proton histograms
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

    // Neutron histograms
    // ======================================================================================================================================================================
    TH1D *h_nsize = new TH1D("nsize", "Number of Neutrons in Event", 10, 0, 10);
    hist_list_1.push_back(h_nsize);

    // reconstructed momentum

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

    // good n / bad n set
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

    // ML features - all neutron candidates
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

    // ML features - signal/background only
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

    const double mP = 0.93828;
    const double mN = 0.939;
    const double mD = 1.8756;

    int numevent = 0;

    // ======================================================================================================================================================================
    // Andrew's histograms
    // ======================================================================================================================================================================

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

    // Checks on which events have neutrons
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

    // Step Zero
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

    // Step One (After Beta Cut)
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

    // Step Two (After applying Phi Diff Charge Track cut)
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

    // Step Three (After applying Phi Diff Charge Track cut)
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

    // Step Four (After applying Phi Diff CND hit cut)
    TH2D *h_pnRes_theta_nmiss_Step4 = new TH2D("pnRes_theta_nmiss_Step4", "(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}", 50, -3.0, 1.0, 90, 0, 180);
    hist_list_2_A.push_back(h_pnRes_theta_nmiss_Step4);

    TH1D *h_ToF_goodN_Step4 = new TH1D("ToF_goodN_Step4", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_goodN_Step4);
    TH1D *h_ToF_badN_Step4 = new TH1D("ToF_badN_Step4", "ToF [ns] of CND Neutrons;ToF;Counts", 100, 0, 20);
    hist_list_1_A.push_back(h_ToF_badN_Step4);

    // Step Five (After event selection cuts)
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

    // ======================================================================================================================================================================
    // Chain loop
    // ======================================================================================================================================================================

    while (chain.Next())
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
        TVector3 pq = pb - pe;

        double nu = Ebeam - pe.Mag();
        double QSq = pq.Mag2() - (nu * nu);
        double xB = QSq / (2 * mN * nu);

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

                if (abs(vzp - vze) > 5)
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

                if (abs(dbeta) > 0.03)
                {
                    continue;
                }
            }
            else if (prot[i]->getRegion() == CD)
            {
                h_vzp_cd->Fill(vzp - vze);

                if (abs(vzp - vze) > 4)
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

                if (abs(dbeta) > 0.05)
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

        //////////////////////////
        //  MISSING MOMENTUM    //
        //////////////////////////

        // missing momentum, energy, mass
        TVector3 pmiss = pq - pp;

        momentum = pmiss.Mag();

        double Ep = sqrt(mN * mN + pp.Mag2());
        double Emiss = Ebeam + mD - pe.Mag() - Ep;
        double mmiss = sqrt((Emiss * Emiss) - pmiss.Mag2());

        //////////////////////////
        ////     NEUTRONS    /////
        //////////////////////////

        // LOOP OVER NEUTRONS
        h_nsize->Fill(neut.size());

        for (int i = 0; i < neut.size(); i++)
        {
            // GET NEUTRON INFORMATION

            // get neutron momentum
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

            if (cnd_energy > 1000)
            {
                continue;
            }

            // ESSENTIAL NEUTRONS CUTS
            h_tof->Fill(time);

            if (pn_x == 0 || pn_y == 0 || pn_z == 0)
            {
                continue;
            }

            if (time > 10)
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

            if (mmiss > 1.)
            {
                continue;
            }

            h_pmiss_thetamiss->Fill(pmiss.Theta() * 180. / M_PI, pmiss.Mag());
            h_thetapn_pp->Fill(pp.Mag(), pp.Angle(pn) * 180. / M_PI);

            if (pmiss.Mag() < 0.25 || pmiss.Mag() > 1.)
            {
                continue;
            }

            if (pmiss.Theta() * 180. / M_PI < 45 || pmiss.Theta() * 180. / M_PI > 140)
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

            if ((fabs(pmiss.Mag() - pn.Mag()) / pmiss.Mag()) > 0.2)
            // if ((abs(pmiss.Mag() - pn.Mag()) / pmiss.Mag()) > 0.2) // Erin's original
            {
                continue;
            }

            h_thetapn_dpp->Fill((pmiss.Mag() - pn.Mag()) / pmiss.Mag(), pn.Angle(pp) * 180. / M_PI);
            h_thetapn_dpp1->Fill((pmiss.Mag() - pn.Mag()) / pmiss.Mag(), pn.Angle(pp) * 180. / M_PI);

            if (pn.Angle(pmiss) * 180. / M_PI > 20)
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

        // chain.WriteEvent();
        counter++;

    } // closes event loop

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

    // ======================================================================================================================================================================
    // Andrew's wrap up
    // ======================================================================================================================================================================

    std::cout << "\nAndrew's wrap up 1\n\n";

    /////////////////////////////////////////////////////
    // Now create the output PDFs
    /////////////////////////////////////////////////////

    int pixelx = 1980;
    int pixely = 1530;

    TCanvas *myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);
    TCanvas *myText = new TCanvas("myText", "myText", pixelx, pixely);

    TLatex text;
    text.SetTextSize(0.05);


    std::cout << "\nAndrew's wrap up 2\n\n";



    char fileName[100];
    sprintf(fileName, "%s[", pdfFile.c_str());
    myText->SaveAs(fileName);
    sprintf(fileName, "%s", pdfFile.c_str());

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text.DrawLatex(0.2, 0.7, "Neutrons in CND");

    myText->Print(fileName, "pdf");
    myText->Clear();

    myCanvas->Divide(1, 1);
    
    std::cout << "\nAndrew's wrap up 3\n\n";

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        myCanvas->cd(1);
        hist_list_1_A[i]->Draw();
        myCanvas->Print(fileName, "pdf");
        myCanvas->Clear();
    }

    std::cout << "\nAndrew's wrap up 4\n\n";

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        myCanvas->cd(1);
        hist_list_2_A[i]->Draw("colz");
        myCanvas->Print(fileName, "pdf");
        myCanvas->Clear();
    }

    std::cout << "\nAndrew's wrap up 5\n\n";

    sprintf(fileName, "%s]", pdfFile.c_str());
    myCanvas->Print(fileName, "pdf");

    return 0;

} // closes main function