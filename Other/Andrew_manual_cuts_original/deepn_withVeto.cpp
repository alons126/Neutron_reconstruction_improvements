#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "clas12reader.h"
#include "HipoChain.h"
#include "eventcut.h"
#include "functions.h"

using namespace std;
using namespace clas12;

void printProgress(double percentage);
void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> <path/to/cutfile.txt> <path/to/input.hipo> \n";
}

bool isPosNear(int sdiff, int ldiff){
  if((ldiff==-2) && (sdiff>=-1) && (sdiff<=0)){return true;}
  if((ldiff==-1) && (sdiff>=-1) && (sdiff<=2)){return true;}
  if((ldiff== 0) && (sdiff>=-1) && (sdiff<=2)){return true;}
  if((ldiff== 1) && (sdiff>=-1) && (sdiff<=2)){return true;}
  if((ldiff== 2) && (sdiff>=-1) && (sdiff<=2)){return true;}
  if((ldiff== 3) && (sdiff>=-1) && (sdiff<=2)){return true;}
  return false;
}


bool isNear(int sdiff, int ldiff){
  /*
  //if((ldiff== 2) && (sdiff==-2)){return true;}
  //if((ldiff== 2) && (sdiff==-1)){return true;}
  if((ldiff== 2) && (sdiff== 0)){return true;}
  if((ldiff== 2) && (sdiff== 1)){return true;}
  if((ldiff== 2) && (sdiff== 2)){return true;}

  //if((ldiff== 1) && (sdiff== 1)){return true;}
  if((ldiff== 1) && (sdiff== 2)){return true;}

  //if((ldiff== 0) && (sdiff== 1)){return true;}
  if((ldiff== 0) && (sdiff== 2)){return true;}

  if((ldiff==-1) && (sdiff== -1)){return true;}
  */
  
  if((ldiff==-2) && (sdiff==-2)){return true;}
  if((ldiff==-2) && (sdiff==-1)){return true;}
  if((ldiff==-2) && (sdiff== 0)){return true;}
  if((ldiff==-2) && (sdiff== 1)){return true;}
  if((ldiff==-2) && (sdiff== 2)){return true;}

  if((ldiff==-1) && (sdiff==-2)){return true;}
  if((ldiff==-1) && (sdiff==-1)){return true;}
  //if((ldiff==-1) && (sdiff== 0)){return true;}
  if((ldiff==-1) && (sdiff== 1)){return true;}
  if((ldiff==-1) && (sdiff== 2)){return true;}

  if((ldiff== 0) && (sdiff==-2)){return true;}
  //if((ldiff== 0) && (sdiff==-1)){return true;}
  //if((ldiff== 0) && (sdiff== 0)){return true;}
  //if((ldiff== 0) && (sdiff== 1)){return true;}
  if((ldiff== 0) && (sdiff== 2)){return true;}

  if((ldiff== 1) && (sdiff==-2)){return true;}
  if((ldiff== 1) && (sdiff==-1)){return true;}
  //if((ldiff== 1) && (sdiff== 0)){return true;}
  if((ldiff== 1) && (sdiff== 1)){return true;}
  if((ldiff== 1) && (sdiff== 2)){return true;}

  if((ldiff== 2) && (sdiff==-2)){return true;}
  if((ldiff== 2) && (sdiff==-1)){return true;}
  if((ldiff== 2) && (sdiff== 0)){return true;}
  if((ldiff== 2) && (sdiff== 1)){return true;}
  if((ldiff== 2) && (sdiff== 2)){return true;}
  
  /*
  //if((ldiff==-1) && (sdiff==-2)){return true;}
  if((ldiff==-1) && (sdiff==-1)){return true;}
  if((ldiff==-1) && (sdiff== 0)){return true;}
  if((ldiff==-1) && (sdiff== 1)){return true;}
  //if((ldiff==-1) && (sdiff== 2)){return true;}

  //if((ldiff== 0) && (sdiff==-4)){return true;}
  //if((ldiff== 0) && (sdiff==-3)){return true;}
  //if((ldiff== 0) && (sdiff==-2)){return true;}
  if((ldiff== 0) && (sdiff==-1)){return true;}
  if((ldiff== 0) && (sdiff== 1)){return true;}
  //if((ldiff== 0) && (sdiff== 2)){return true;}
  //if((ldiff== 0) && (sdiff== 3)){return true;}

  //if((ldiff== 1) && (sdiff==-2)){return true;}
  if((ldiff== 1) && (sdiff==-1)){return true;}
  if((ldiff== 1) && (sdiff== 0)){return true;}
  if((ldiff== 1) && (sdiff== 1)){return true;}
  //if((ldiff== 1) && (sdiff== 2)){return true;}
  */
  return false;
}

bool isNearCTOF(int sdiff, int ldiff){
  if((ldiff== 1) && (sdiff==-3)){return true;}
  if((ldiff== 1) && (sdiff==-2)){return true;}
  if((ldiff== 1) && (sdiff==-1)){return true;}
  if((ldiff== 1) && (sdiff== 1)){return true;}
  if((ldiff== 1) && (sdiff== 2)){return true;}
  if((ldiff== 1) && (sdiff== 3)){return true;}

  if((ldiff== 2) && (sdiff==-3)){return true;}
  if((ldiff== 2) && (sdiff==-2)){return true;}
  if((ldiff== 2) && (sdiff==-1)){return true;}
  if((ldiff== 2) && (sdiff== 0)){return true;}
  if((ldiff== 2) && (sdiff== 1)){return true;}
  if((ldiff== 2) && (sdiff== 2)){return true;}
  if((ldiff== 2) && (sdiff== 3)){return true;}

  if((ldiff== 3) && (sdiff==-3)){return true;}
  if((ldiff== 3) && (sdiff==-2)){return true;}
  if((ldiff== 3) && (sdiff==-1)){return true;}
  if((ldiff== 3) && (sdiff== 0)){return true;}
  if((ldiff== 3) && (sdiff== 1)){return true;}
  if((ldiff== 3) && (sdiff== 2)){return true;}
  if((ldiff== 3) && (sdiff== 3)){return true;}

  return false;
}

int main(int argc, char ** argv)
{

  if(argc < 7)
    {
      std::cerr<<"Wrong number of arguments.\n";
      Usage();
      return -1;
    }

  /////////////////////////////////////
  
  bool isMC = false;
  if(atoi(argv[1]) == 1){isMC=true;}
  TRandom3 myRand(0);

  double Ebeam = atof(argv[2]);
  
  TFile * outFile = new TFile(argv[3],"RECREATE");
  char * pdfFile = argv[4];
  eventcut myCut(Ebeam,argv[5]);
  myCut.print_cuts();
  clas12root::HipoChain chain;
  for(int k = 6; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  //chain.NextFile();//Need to open a file to get entry orders
  auto config_c12=chain.GetC12Reader(); 
  //Get the important CND information
  auto cnd_hits = config_c12->addBank("CND::hits");
  auto cnd_hit_layer = config_c12->getBankOrder(cnd_hits,"layer");  
  auto cnd_hit_sector = config_c12->getBankOrder(cnd_hits,"sector");  
  auto cnd_hit_energy = config_c12->getBankOrder(cnd_hits,"energy");  
  auto cnd_hit_time = config_c12->getBankOrder(cnd_hits,"time");  
  auto cnd_hit_z = config_c12->getBankOrder(cnd_hits,"z");  
  //Get the important CTOF information  
  auto ctof_hits = config_c12->addBank("CTOF::hits");
  auto ctof_hit_component = config_c12->getBankOrder(ctof_hits,"component");  
  auto ctof_hit_energy = config_c12->getBankOrder(ctof_hits,"energy");  
  auto ctof_hit_time = config_c12->getBankOrder(ctof_hits,"time"); 
  auto ctof_hit_z = config_c12->getBankOrder(ctof_hits,"z");  
  //Now define the object you will call
  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();
  /////////////////////////////////////
  //Prepare histograms
  /////////////////////////////////////

  vector<TH1*> hist_list_1;
  vector<TH2*> hist_list_2;

  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);

  gStyle->SetTitleXOffset(0.8);
  gStyle->SetTitleYOffset(0.8);

  char temp_name[100];
  char temp_title[100];

  //Checks on which events have neutrons
  TH2D * h_xB_mmiss_epFD = new TH2D("xB_mmiss_epFD","x_{B} vs. m_{miss};x_{B};m_{miss}",100,0.0,2.0,100,0.5,1.5);
  hist_list_2.push_back(h_xB_mmiss_epFD);
  TH2D * h_xB_mmiss_epnFD = new TH2D("xB_mmiss_epnFD","x_{B} vs. m_{miss};x_{B};m_{miss}",100,0.0,2.0,100,0.5,1.5);
  hist_list_2.push_back(h_xB_mmiss_epnFD);
  TH2D * h_xB_mmiss_epngoodFD = new TH2D("xB_mmiss_epngoodFD","x_{B} vs. m_{miss};x_{B};m_{miss}",100,0.0,2.0,100,0.5,1.5);
  hist_list_2.push_back(h_xB_mmiss_epngoodFD);
  TH2D * h_xB_mmiss_epCD = new TH2D("xB_mmiss_epCD","x_{B} vs. m_{miss};x_{B};m_{miss}",100,0.0,2.0,100,0.5,1.5);
  hist_list_2.push_back(h_xB_mmiss_epCD);
  TH2D * h_xB_mmiss_epnCD = new TH2D("xB_mmiss_epnCD","x_{B} vs. m_{miss};x_{B};m_{miss}",100,0.0,2.0,100,0.5,1.5);
  hist_list_2.push_back(h_xB_mmiss_epnCD);
  TH2D * h_xB_mmiss_epngoodCD = new TH2D("xB_mmiss_epngoodCD","x_{B} vs. m_{miss};x_{B};m_{miss}",100,0.0,2.0,100,0.5,1.5);
  hist_list_2.push_back(h_xB_mmiss_epngoodCD);

  TH1D * h_pmiss_ep = new TH1D("pmiss_ep","p_{miss} ep;p_{miss};Counts",25,0.25,1.0);
  hist_list_1.push_back(h_pmiss_ep);

  //Step Zero
  TH2D * h_pnRes_theta_nmiss_Step0 = new TH2D("pnRes_theta_nmiss_Step0","(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}",50,-3.0,1.0,90,0,180);
  hist_list_2.push_back(h_pnRes_theta_nmiss_Step0);

  TH1D * h_ToF_goodN_Step0 = new TH1D("ToF_goodN_Step0","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_goodN_Step0);
  TH1D * h_ToF_badN_Step0 = new TH1D("ToF_badN_Step0","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_badN_Step0);

  TH1D * h_beta_goodN_Step0 = new TH1D("beta_goodN_Step0","#beta of CND Neutrons;#beta;Counts",50,0,1.1);
  hist_list_1.push_back(h_beta_goodN_Step0);
  TH1D * h_Edep_goodN_Step0 = new TH1D("Edep_goodN_Step0","E_{dep} [MeF] of CND Neutrons;E_{dep};Counts",50,0,100);
  hist_list_1.push_back(h_Edep_goodN_Step0);
  TH2D * h_beta_Edep_goodN_Step0 = new TH2D("Edep_beta_goodN","#beta vs. E_{dep} [MeV] of CND Neutrons;#beta;E_{dep}",50,0,1.1,50,0,100);
  hist_list_2.push_back(h_beta_Edep_goodN_Step0);

  TH1D * h_beta_badN_Step0 = new TH1D("beta_badN_Step0","#beta of CND Neutrons;#beta;Counts",50,0,1.1);
  hist_list_1.push_back(h_beta_badN_Step0);
  TH1D * h_Edep_badN_Step0 = new TH1D("Edep_badN_Step0","E_{dep} [MeF] of CND Neutrons;E_{dep};Counts",50,0,100);
  hist_list_1.push_back(h_Edep_badN_Step0);
  TH2D * h_beta_Edep_badN_Step0 = new TH2D("Edep_beta_badN","#beta vs. E_{dep} [MeV] of CND Neutrons;#beta;E_{dep}",50,0,1.1,50,0,100);
  hist_list_2.push_back(h_beta_Edep_badN_Step0);


  //Step One (After Beta Cut)
  TH2D * h_pnRes_theta_nmiss_Step1 = new TH2D("pnRes_theta_nmiss_Step1","(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}",50,-3.0,1.0,90,0,180);
  hist_list_2.push_back(h_pnRes_theta_nmiss_Step1);

  TH1D * h_pmiss_goodN_Step1 = new TH1D("pmiss_goodN_Step1","p_{miss} Step1;p_{miss};Counts",25,0.25,1.0);
  hist_list_1.push_back(h_pmiss_goodN_Step1);
  TH1D * h_ToF_goodN_Step1 = new TH1D("ToF_goodN_Step1","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_goodN_Step1);
  TH1D * h_ToF_badN_Step1 = new TH1D("ToF_badN_Step1","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_badN_Step1);

  TH1D * h_edep_goodN_Step1 = new TH1D("edep_goodN_Step1","edep [MeV] of CND Neutrons;edep;Counts",100,0,50);
  hist_list_1.push_back(h_edep_goodN_Step1);
  TH1D * h_edep_badN_Step1 = new TH1D("edep_badN_Step1","edep [MeV] of CND Neutrons;edep;Counts",100,0,50);
  hist_list_1.push_back(h_edep_badN_Step1);

  TH1D * h_edep_over_edepCTOT_goodN_Step1 = new TH1D("edep_over_edepCTOT_goodN_Step1","E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts",100,0,5);
  hist_list_1.push_back(h_edep_over_edepCTOT_goodN_Step1);
  TH1D * h_edep_over_edepCTOT_badN_Step1 = new TH1D("edep_over_edepCTOT_badN_Step1","E_{dep,N}/E_{dep,pos};E_{dep,N}/E_{dep,pos};Counts",100,0,5);
  hist_list_1.push_back(h_edep_over_edepCTOT_badN_Step1);

  TH1D * h_edep_goodN_withNearbyPos_Step1 = new TH1D("edep_goodN_withNearbyPos_Step1","edep [MeV] of CND Neutrons;edep;Counts",100,0,50);
  hist_list_1.push_back(h_edep_goodN_withNearbyPos_Step1);
  TH1D * h_edep_badN_withNearbyPos_Step1 = new TH1D("edep_badN_withNearbyPos_Step1","edep [MeV] of CND Neutrons;edep;Counts",100,0,50);
  hist_list_1.push_back(h_edep_badN_withNearbyPos_Step1);

  TH1D * h_sdiff_pos_goodN_Step1_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_goodN_Step1_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)",k-3);
    h_sdiff_pos_goodN_Step1_layer[k] = new TH1D(temp_name,temp_title,24,-11.5,12.5);
    hist_list_1.push_back(h_sdiff_pos_goodN_Step1_layer[k]);
  }

  TH1D * h_sdiff_pos_badN_Step1_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_badN_Step1_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)",k-3);
    h_sdiff_pos_badN_Step1_layer[k] = new TH1D(temp_name,temp_title,24,-11.5,12.5);
    hist_list_1.push_back(h_sdiff_pos_badN_Step1_layer[k]);
  }

  TH2D * h_sdiff_pos_mom_goodN_Step1_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_mom_goodN_Step1_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)",k-3);
    h_sdiff_pos_mom_goodN_Step1_layer[k] = new TH2D(temp_name,temp_title,24,-11.5,12.5,50,0,4);
    hist_list_2.push_back(h_sdiff_pos_mom_goodN_Step1_layer[k]);
  }
  TH2D * h_sdiff_pos_mom_badN_Step1_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_mom_badN_Step1_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector vs. Momentum Proton (Layer Difference = %d)",k-3);
    h_sdiff_pos_mom_badN_Step1_layer[k] = new TH2D(temp_name,temp_title,24,-11.5,12.5,50,0,4);
    hist_list_2.push_back(h_sdiff_pos_mom_badN_Step1_layer[k]);
  }

  TH2D * h_sdiff_pos_z_goodN_Step1_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_z_goodN_Step1_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)",k-3);
    h_sdiff_pos_z_goodN_Step1_layer[k] = new TH2D(temp_name,temp_title,24,-11.5,12.5,50,-40.0,40.0);
    hist_list_2.push_back(h_sdiff_pos_z_goodN_Step1_layer[k]);
  }
  TH2D * h_sdiff_pos_z_badN_Step1_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_z_badN_Step1_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector vs. Z (Layer Difference = %d)",k-3);
    h_sdiff_pos_z_badN_Step1_layer[k] = new TH2D(temp_name,temp_title,24,-11.5,12.5,50,-40.0,40.0);
    hist_list_2.push_back(h_sdiff_pos_z_badN_Step1_layer[k]);
  }


  TH2D * h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_diff_ToFc_z_goodN_Step1_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)",k-3);
    h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[k] = new TH2D(temp_name,temp_title,24,-11.5,12.5,50,0,300);
    hist_list_2.push_back(h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[k]);
  }
  TH2D * h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_diff_ToFc_z_badN_Step1_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector vs. ToF*c-z (Layer Difference = %d)",k-3);
    h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[k] = new TH2D(temp_name,temp_title,24,-11.5,12.5,50,0,300);
    hist_list_2.push_back(h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[k]);
  }

  TH2D * h_diff_ToFc_z_Edep_noNear_goodN_Step1 = new TH2D("diff_ToFc_z_Edep_noNear_goodN_Step1","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons with no Nearby Tracks;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_noNear_goodN_Step1);
  TH2D * h_diff_ToFc_z_Edep_yesNear_goodN_Step1 = new TH2D("diff_ToFc_z_Edep_yesNear_goodN_Step1","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons with no Nearby Tracks;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_yesNear_goodN_Step1);
  TH2D * h_diff_ToFc_z_Edep_noNear_badN_Step1 = new TH2D("diff_ToFc_z_Edep_noNear_badN_Step1","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons with no Nearby Tracks;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_noNear_badN_Step1);
  TH2D * h_diff_ToFc_z_Edep_yesNear_badN_Step1 = new TH2D("diff_ToFc_z_Edep_yesNear_badN_Step1","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons with no Nearby Tracks;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_yesNear_badN_Step1);

  //Step Two (After applying Phi Diff Charge Track cut)
  TH2D * h_pnRes_theta_nmiss_Step2 = new TH2D("pnRes_theta_nmiss_Step2","(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}",50,-3.0,1.0,90,0,180);
  hist_list_2.push_back(h_pnRes_theta_nmiss_Step2);

  TH1D * h_ToF_goodN_Step2 = new TH1D("ToF_goodN_Step2","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_goodN_Step2);
  TH1D * h_ToF_badN_Step2 = new TH1D("ToF_badN_Step2","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_badN_Step2);

  TH1D * h_sdiff_pos_goodN_Step2_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_goodN_Step2_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)",k-3);
    h_sdiff_pos_goodN_Step2_layer[k] = new TH1D(temp_name,temp_title,24,-11.5,12.5);
    hist_list_1.push_back(h_sdiff_pos_goodN_Step2_layer[k]);
  }

  TH1D * h_sdiff_pos_badN_Step2_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_pos_badN_Step2_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus +Charge Particle Sector (Layer Difference = %d)",k-3);
    h_sdiff_pos_badN_Step2_layer[k] = new TH1D(temp_name,temp_title,24,-11.5,12.5);
    hist_list_1.push_back(h_sdiff_pos_badN_Step2_layer[k]);
  }

  TH1D * h_sdiff_allhit_goodN_Step2_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_allhit_goodN_Step2_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus Random Hit Sector (Layer Difference = %d)",k-3);
    h_sdiff_allhit_goodN_Step2_layer[k] = new TH1D(temp_name,temp_title,24,-11.5,12.5);
    hist_list_1.push_back(h_sdiff_allhit_goodN_Step2_layer[k]);
  }
  TH2D * h_sdiff_ldiff_allhit_goodN_Step2 = new TH2D("sdiff_ldiff_allhit_goodN_Step2","Sector Difference vs. Layer Difference;Sector Difference;Layer Difference",24,-11.5,12.5,7,-3.5,3.5);
  hist_list_2.push_back(h_sdiff_ldiff_allhit_goodN_Step2);

  TH1D * h_sdiff_allhit_badN_Step2_layer[7];
  for(int k = 0; k < 7; k++){
    sprintf(temp_name,"sdiff_allhit_badN_Step2_layer_%d",k-3);
    sprintf(temp_title,"Nuetral Sector minus Random Hit Sector (Layer Difference = %d)",k-3);
    h_sdiff_allhit_badN_Step2_layer[k] = new TH1D(temp_name,temp_title,24,-11.5,12.5);
    hist_list_1.push_back(h_sdiff_allhit_badN_Step2_layer[k]);
  }
  TH2D * h_sdiff_ldiff_allhit_badN_Step2 = new TH2D("sdiff_ldiff_allhit_badN_Step2","Sector Difference vs. Layer Difference;Sector Difference;Layer Difference",24,-11.5,12.5,7,-3.5,3.5);
  hist_list_2.push_back(h_sdiff_ldiff_allhit_badN_Step2);

  TH1D * h_numberNearby_goodN_Step2 = new TH1D("numberNearby_goodN_Step2","Number of Nearby Hits for CND Neutrons;# Hits;Counts",9,-0.5,8.5);
  hist_list_1.push_back(h_numberNearby_goodN_Step2);
  TH2D * h_numberNearby_momN_goodN_Step2 = new TH2D("numberNearby_momN_goodN_Step2","Number of Nearby Hits vs. p_{n} for CND Neutrons;# Hits;p_{n}",9,-0.5,8.5,50,0,1.3);
  hist_list_2.push_back(h_numberNearby_momN_goodN_Step2);
  TH1D * h_numberNearby_badN_Step2 = new TH1D("numberNearby_badN_Step2","Number of Nearby Hits for CND Neutrons;# Hits;Counts",9,-0.5,8.5);
  hist_list_1.push_back(h_numberNearby_badN_Step2);
  TH2D * h_numberNearby_momN_badN_Step2 = new TH2D("numberNearby_momN_badN_Step2","Number of Nearby Hits vs. p_{n} for CND Neutrons;# Hits;p_{n}",9,-0.5,8.5,50,0,1.3);
  hist_list_2.push_back(h_numberNearby_momN_badN_Step2);

  TH1D * h_NearbyEdep_goodN_Step2 = new TH1D("NearbyEdep_goodN_Step2","E_{dep} of Nearby Hits for CND Neutrons;E_{dep};Counts",50,0,100);
  hist_list_1.push_back(h_NearbyEdep_goodN_Step2);
  TH1D * h_NearbyEdep_badN_Step2 = new TH1D("NearbyEdep_badN_Step2","E_{dep} of Nearby Hits for CND Neutrons;E_{dep};Counts",50,0,100);
  hist_list_1.push_back(h_NearbyEdep_badN_Step2);


  TH1D * h_nsector_goodN_Step2 = new TH1D("nsector_goodN_Step2","Neutron Sector for CND Neutrons;Neutron Sector;Counts",24,0.5,24.5);
  hist_list_1.push_back(h_nsector_goodN_Step2);
  TH1D * h_nsector_badN_Step2 = new TH1D("nsector_badN_Step2","Neutron Sector for CND Neutrons;Neutron Sector;Counts",24,0.5,24.5);
  hist_list_1.push_back(h_nsector_badN_Step2);

  //Step Three (After applying Phi Diff Charge Track cut)
  TH2D * h_pnRes_theta_nmiss_Step3 = new TH2D("pnRes_theta_nmiss_Step3","(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}",50,-3.0,1.0,90,0,180);
  hist_list_2.push_back(h_pnRes_theta_nmiss_Step3);

  TH1D * h_ToF_goodN_Step3 = new TH1D("ToF_goodN_Step3","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_goodN_Step3);
  TH1D * h_ToF_badN_Step3 = new TH1D("ToF_badN_Step3","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_badN_Step3);

  TH2D * h_sdiff_ldiff_allhit_goodN_Step3 = new TH2D("sdiff_ldiff_allhit_goodN_Step3","Sector Difference vs. Layer Difference;Sector Difference;Layer Difference",24,-11.5,12.5,7,-3.5,3.5);
  hist_list_2.push_back(h_sdiff_ldiff_allhit_goodN_Step3);
  TH2D * h_sdiff_ldiff_allhit_badN_Step3 = new TH2D("sdiff_ldiff_allhit_badN_Step3","Sector Difference vs. Layer Difference;Sector Difference;Layer Difference",24,-11.5,12.5,7,-3.5,3.5);
  hist_list_2.push_back(h_sdiff_ldiff_allhit_badN_Step3);

  TH2D * h_sdiff_ldiff_CTOFhit_goodN_Step3 = new TH2D("sdiff_ldiff_CTOFhit_goodN_Step3","Sector Difference vs. Layer Difference;Sector Difference;Layer Difference",24,-11.5,12.5,3,0.5,3.5);
  hist_list_2.push_back(h_sdiff_ldiff_CTOFhit_goodN_Step3);
  TH2D * h_sdiff_ldiff_CTOFhit_badN_Step3 = new TH2D("sdiff_ldiff_CTOFhit_badN_Step3","Sector Difference vs. Layer Difference;Sector Difference;Layer Difference",24,-11.5,12.5,3,0.5,3.5);
  hist_list_2.push_back(h_sdiff_ldiff_CTOFhit_badN_Step3);

  TH1D * h_numberCTOF_goodN_Step3 = new TH1D("numberCTOF_goodN_Step3","Number of CTOF Hits for CND Neutrons;# Hits;Counts",9,-0.5,8.5);
  hist_list_1.push_back(h_numberCTOF_goodN_Step3);
  TH2D * h_numberCTOF_momN_goodN_Step3 = new TH2D("numberCTOF_momN_goodN_Step3","Number of CTOF Hits vs. p_{n} for CND Neutrons;# Hits;p_{n}",9,-0.5,8.5,50,0,1.3);
  hist_list_2.push_back(h_numberCTOF_momN_goodN_Step3);
  TH1D * h_numberCTOF_badN_Step3 = new TH1D("numberCTOF_badN_Step3","Number of CTOF Hits for CND Neutrons;# Hits;Counts",9,-0.5,8.5);
  hist_list_1.push_back(h_numberCTOF_badN_Step3);
  TH2D * h_numberCTOF_momN_badN_Step3 = new TH2D("numberCTOF_momN_badN_Step3","Number of CTOF Hits vs. p_{n} for CND Neutrons;# Hits;p_{n}",9,-0.5,8.5,50,0,1.3);
  hist_list_2.push_back(h_numberCTOF_momN_badN_Step3);


  //Step Four (After applying Phi Diff CND hit cut)
  TH2D * h_pnRes_theta_nmiss_Step4 = new TH2D("pnRes_theta_nmiss_Step4","(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}",50,-3.0,1.0,90,0,180);
  hist_list_2.push_back(h_pnRes_theta_nmiss_Step4);

  TH1D * h_ToF_goodN_Step4 = new TH1D("ToF_goodN_Step4","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_goodN_Step4);
  TH1D * h_ToF_badN_Step4 = new TH1D("ToF_badN_Step4","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_badN_Step4);


  //Step Five (After event selection cuts)
  TH2D * h_pnRes_theta_nmiss_Step5 = new TH2D("pnRes_theta_nmiss_Step5","(p_{miss}-p_{n})/p_{miss} vs. #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}",50,-3.0,1.0,90,0,180);
  hist_list_2.push_back(h_pnRes_theta_nmiss_Step5);

  TH1D * h_ToF_goodN_Step5 = new TH1D("ToF_goodN_Step5","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_goodN_Step5);
  TH1D * h_ToF_badN_Step5 = new TH1D("ToF_badN_Step5","ToF [ns] of CND Neutrons;ToF;Counts",100,0,20);
  hist_list_1.push_back(h_ToF_badN_Step5);
  TH1D * h_pmiss_goodN_Step5 = new TH1D("pmiss_goodN_Step5","p_{miss} good N Step5;p_{miss};Counts",25,0.25,1.0);
  hist_list_1.push_back(h_pmiss_goodN_Step5);
  TH1D * h_pmiss_allN_Step5 = new TH1D("pmiss_allN_Step5","p_{miss} all N Step5;p_{miss};Counts",25,0.25,1.0);
  hist_list_1.push_back(h_pmiss_allN_Step5);

  TH1D * h_Edep_infront_goodN_Step5 = new TH1D("Edep_infront_goodN_Step5","E_{dep} [MeV] of Hit infront CND Neutrons;E_{dep};Counts",50,0,100);
  hist_list_1.push_back(h_Edep_infront_goodN_Step5);
  TH1D * h_Edep_behind_goodN_Step5 = new TH1D("Edep_behind_goodN_Step5","E_{dep} [MeV] of Hit behind CND Neutrons;E_{dep};Counts",50,0,100);
  hist_list_1.push_back(h_Edep_behind_goodN_Step5);

  TH1D * h_Edep_infront_badN_Step5 = new TH1D("Edep_infront_badN_Step5","E_{dep} [MeV] of Hit infront CND Neutrons;E_{dep};Counts",50,0,100);
  hist_list_1.push_back(h_Edep_infront_badN_Step5);
  TH1D * h_Edep_behind_badN_Step5 = new TH1D("Edep_behind_badN_Step5","E_{dep} [MeV] of Hit behind CND Neutrons;E_{dep};Counts",50,0,100);
  hist_list_1.push_back(h_Edep_behind_badN_Step5);

  TH2D * h_diff_ToFc_z_Edep_goodN_Step5 = new TH2D("diff_ToFc_z_Edep_goodN_Step5","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_goodN_Step5);
  TH2D * h_diff_ToFc_z_Edep_badN_Step5 = new TH2D("diff_ToFc_z_Edep_badN_Step5","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_badN_Step5);

  TH2D * h_diff_ToFc_z_Edep_goodN_Step5_layer[3];
  for(int k = 0; k < 3; k++){
    sprintf(temp_name,"diff_ToFc_z_goodN_Step5_layer_%d",k+1);
    sprintf(temp_title,"ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons (Layer Difference = %d);ToF*c-z;E_{dep}",k+1);
    h_diff_ToFc_z_Edep_goodN_Step5_layer[k] = new TH2D(temp_name,temp_title,50,0,300,50,0,100);
    hist_list_2.push_back(h_diff_ToFc_z_Edep_goodN_Step5_layer[k]);
  }

  TH2D * h_diff_ToFc_z_Edep_badN_Step5_layer[3];
  for(int k = 0; k < 3; k++){
    sprintf(temp_name,"diff_ToFc_z_badN_Step5_layer_%d",k+1);
    sprintf(temp_title,"ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons (Layer Difference = %d);ToF*c-z;E_{dep}",k+1);
    h_diff_ToFc_z_Edep_badN_Step5_layer[k] = new TH2D(temp_name,temp_title,50,0,300,50,0,100);
    hist_list_2.push_back(h_diff_ToFc_z_Edep_badN_Step5_layer[k]);
  }
  
  TH1D * h_phidiff_en_goodN_Step5 = new TH1D("phidiff_en_goodN_Step5","|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}|;Counts",90,0,180);
  hist_list_1.push_back(h_phidiff_en_goodN_Step5);

  TH1D * h_phidiff_en_badN_Step5 = new TH1D("phidiff_en_badN_Step5","|#phi_{e}-#phi_{n}| of CND Neutrons;|#phi_{e}-#phi_{n}|;Counts",90,0,180);
  hist_list_1.push_back(h_phidiff_en_badN_Step5);

  TH1D * h_TP_goodN_Step5 = new TH1D("TP_goodN_Step5","ToF/path [ns/m] of CND Neutrons;ToF/path;Counts",150,0,50);
  hist_list_1.push_back(h_TP_goodN_Step5);
  TH1D * h_TP_badN_Step5 = new TH1D("TP_badN_Step5","ToF/path [ns/m] of CND Neutrons;ToF/path;Counts",150,0,50);
  hist_list_1.push_back(h_TP_badN_Step5);

  TH1D * h_Z_goodN_Step5 = new TH1D("Z_goodN_Step5","Z [cm] of CND Neutrons;Z;Counts",100,-60,60);
  hist_list_1.push_back(h_Z_goodN_Step5);
  TH1D * h_Z_badN_Step5 = new TH1D("Z_badN_Step5","Z [cm] of CND Neutrons;Z;Counts",100,-60,60);
  hist_list_1.push_back(h_Z_badN_Step5);

  TH2D * h_beta_Edep_goodN_Step5 = new TH2D("Edep_beta_goodN","#beta vs. E_{dep} [MeV] of CND Neutrons;#beta;E_{dep}",50,0,1.1,50,0,100);
  hist_list_2.push_back(h_beta_Edep_goodN_Step5);
  TH2D * h_beta_Edep_badN_Step5 = new TH2D("Edep_beta_badN","#beta vs. E_{dep} [MeV] of CND Neutrons;#beta;E_{dep}",50,0,1.1,50,0,100);
  hist_list_2.push_back(h_beta_Edep_badN_Step5);

  TH2D * h_ToF_Edep_goodN_Step5 = new TH2D("ToF_Edep_goodN_Step5","ToF [ns] vs. E_{dep} [MeV] of CND Neutrons;ToF;E_{dep} MeV",100,0,20,50,0,100);
  hist_list_2.push_back(h_ToF_Edep_goodN_Step5);
  TH2D * h_ToF_Edep_badN_Step5 = new TH2D("ToF_Edep_badN_Step5","ToF [ns] vs. E_{dep} [MeV] of CND Neutrons;ToF;E_{dep} MeV",100,0,20,50,0,100);
  hist_list_2.push_back(h_ToF_Edep_badN_Step5);

  TH2D * h_TP_Edep_goodN_Step5 = new TH2D("TP_Edep_goodN_Step5","TP [ns/m] vs. E_{dep} [MeV] of CND Neutrons;TP;E_{dep} MeV",150,0,50,50,0,100);
  hist_list_2.push_back(h_TP_Edep_goodN_Step5);
  TH2D * h_TP_Edep_badN_Step5 = new TH2D("TP_Edep_badN_Step5","TP [ns/m] vs. E_{dep} [MeV] of CND Neutrons;TP;E_{dep} MeV",150,0,50,50,0,100);
  hist_list_2.push_back(h_TP_Edep_badN_Step5);

  ////////  

  /*
  TH1D * h_ToF_bad = new TH1D("ToF_bad","ToF [ns] of CND Neutrons;ToF;Counts",500,0,15);
  hist_list_1.push_back(h_ToF_bad);
  TH2D * h_Edep_z_bad = new TH2D("Edep_z_bad","E_{dep} [MeV] vs. z [cm] of CND Neutrons;E_{dep};z",200,0,100,200,-60,60);
  hist_list_2.push_back(h_Edep_z_bad);
  TH2D * h_Edep_ToF_bad = new TH2D("Edep_ToF_bad","E_{dep} [MeV] vs. ToF [ns] of CND Neutrons;E_{dep};ToF",200,0,100,500,-4,10);
  hist_list_2.push_back(h_Edep_ToF_bad);
  TH2D * h_ToF_z_bad = new TH2D("ToF_z_bad","ToF [ns] vs. z [cm] of CND Neutrons;ToF;z",500,-4,10,200,-60,60);
  hist_list_2.push_back(h_ToF_z_bad);

  TH1D * h_ToF_goodN = new TH1D("ToF_goodN","ToF [ns] of CND Neutrons;ToF;Counts",100,0,50);
  hist_list_1.push_back(h_ToF_goodN);
  TH1D * h_TM_goodN = new TH1D("TM_goodN","ToF/path [ns/m] of CND Neutrons;ToF/path;Counts",20,0,20);
  hist_list_1.push_back(h_TM_goodN);
  TH1D * h_beta_goodN = new TH1D("beta_goodN","#beta of CND Neutrons;#beta;Counts",20,0,1);
  hist_list_1.push_back(h_beta_goodN);
  TH1D * h_mom_goodN = new TH1D("p_goodN","Momentum of CND Neutrons;p;Counts",20,0,2);
  hist_list_1.push_back(h_mom_goodN);
  TH2D * h_Edep_z_goodN = new TH2D("Edep_z_goodN","E_{dep} [MeV] vs. z [cm] of CND Neutrons;E_{dep};z",100,0,100,120,-60,60);
  hist_list_2.push_back(h_Edep_z_goodN);
  TH2D * h_Edep_ToF_goodN = new TH2D("Edep_ToF_goodN","E_{dep} [MeV] vs. ToF [ns] of CND Neutrons;E_{dep};ToF",100,0,100,50,0,30);
  hist_list_2.push_back(h_Edep_ToF_goodN);
  TH2D * h_beta_z_goodN = new TH2D("beta_z_goodN","#beta [ns] vs. z [cm] of CND Neutrons;ToF;z",50,0,1,120,-60,60);
  hist_list_2.push_back(h_beta_z_goodN);
  TH2D * h_Edep_z_1_goodN = new TH2D("Edep_z_1_goodN","E_{dep} [MeV] vs. z [cm] of CND1 Neutrons;E_{dep};z",100,0,100,120,-40,40);
  hist_list_2.push_back(h_Edep_z_1_goodN);
  TH2D * h_Edep_z_2_goodN = new TH2D("Edep_z_2_goodN","E_{dep} [MeV] vs. z [cm] of CND2 Neutrons;E_{dep};z",100,0,100,120,-40,40);
  hist_list_2.push_back(h_Edep_z_2_goodN);
  TH2D * h_Edep_z_3_goodN = new TH2D("Edep_z_3_goodN","E_{dep} [MeV] vs. z [cm] of CND3 Neutrons;E_{dep};z",100,0,100,120,-40,40);
  hist_list_2.push_back(h_Edep_z_3_goodN);
  TH2D * h_ToFc_z_1_goodN = new TH2D("ToFc_z_1_goodN","ToF*c [cm] vs. z [cm] of CND Neutrons;ToF*c;z",50,12,120,50,-40,40);
  hist_list_2.push_back(h_ToFc_z_1_goodN);
  TH2D * h_ToFc_z_2_goodN = new TH2D("ToFc_z_2_goodN","ToF*c [cm] vs. z [cm] of CND Neutrons;ToF*c;z",50,12,120,50,-40,40);
  hist_list_2.push_back(h_ToFc_z_2_goodN);
  TH2D * h_ToFc_z_3_goodN = new TH2D("ToFc_z_3_goodN","ToF*c [cm] vs. z [cm] of CND Neutrons;ToF*c;z",50,12,120,50,-40,40);
  hist_list_2.push_back(h_ToFc_z_3_goodN);
  TH1D * h_diff_ToFc_z_1_goodN = new TH1D("diff_ToFc_z_1_goodN","ToF*c - z [cm] of CND Neutrons;ToF*c-z",100,0,300);
  hist_list_1.push_back(h_diff_ToFc_z_1_goodN);
  TH1D * h_diff_ToFc_z_2_goodN = new TH1D("diff_ToFc_z_2_goodN","ToF*c - z [cm] of CND Neutrons;ToF*c-z",100,0,300);
  hist_list_1.push_back(h_diff_ToFc_z_2_goodN);
  TH1D * h_diff_ToFc_z_3_goodN = new TH1D("diff_ToFc_z_3_goodN","ToF*c - z [cm] of CND Neutrons;ToF*c-z",100,0,300);
  hist_list_1.push_back(h_diff_ToFc_z_3_goodN);
  TH2D * h_diff_ToFc_z_Edep_1_goodN = new TH2D("diff_ToFc_z_Edep_1_goodN","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_1_goodN);
  TH2D * h_diff_ToFc_z_Edep_2_goodN = new TH2D("diff_ToFc_z_Edep_2_goodN","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_2_goodN);
  TH2D * h_diff_ToFc_z_Edep_3_goodN = new TH2D("diff_ToFc_z_Edep_3_goodN","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_3_goodN);
  TH2D * h_Edep_mom_goodN = new TH2D("Edep_mom_goodN","E_{dep} [MeV] vs. Momentum [GeV] of CND1 Neutrons;E_{dep};Momentum",100,0,100,100,0,1.5);
  hist_list_2.push_back(h_Edep_mom_goodN);





  TH1D * h_ToF_badN = new TH1D("ToF_badN","ToF [ns] of CND Neutrons;ToF;Counts",100,0,50);
  hist_list_1.push_back(h_ToF_badN);
  TH1D * h_TM_badN = new TH1D("TM_badN","ToF/path [ns/m] of CND Neutrons;ToF/path;Counts",20,0,20);
  hist_list_1.push_back(h_TM_badN);
  TH1D * h_beta_badN = new TH1D("beta_badN","#beta of CND Neutrons;#beta;Counts",20,0,1);
  hist_list_1.push_back(h_beta_badN);
  TH1D * h_mom_badN = new TH1D("p_badN","Momentum of CND Neutrons;p;Counts",20,0,2);
  hist_list_1.push_back(h_mom_badN);
  TH2D * h_Edep_z_badN = new TH2D("Edep_z_badN","E_{dep} [MeV] vs. z [cm] of CND Neutrons;E_{dep};z",100,0,100,120,-60,60);
  hist_list_2.push_back(h_Edep_z_badN);
  TH2D * h_Edep_ToF_badN = new TH2D("Edep_ToF_badN","E_{dep} [MeV] vs. ToF [ns] of CND Neutrons;E_{dep};ToF",100,0,100,50,0,30);
  hist_list_2.push_back(h_Edep_ToF_badN);
  TH2D * h_beta_z_badN = new TH2D("beta_z_badN","#beta [ns] vs. z [cm] of CND Neutrons;ToF;z",50,0,1,120,-60,60);
  hist_list_2.push_back(h_beta_z_badN);
  TH2D * h_Edep_z_1_badN = new TH2D("Edep_z_1_badN","E_{dep} [MeV] vs. z [cm] of CND1 Neutrons;E_{dep};z",100,0,100,120,-40,40);
  hist_list_2.push_back(h_Edep_z_1_badN);
  TH2D * h_Edep_z_2_badN = new TH2D("Edep_z_2_badN","E_{dep} [MeV] vs. z [cm] of CND2 Neutrons;E_{dep};z",100,0,100,120,-40,40);
  hist_list_2.push_back(h_Edep_z_2_badN);
  TH2D * h_Edep_z_3_badN = new TH2D("Edep_z_3_badN","E_{dep} [MeV] vs. z [cm] of CND3 Neutrons;E_{dep};z",100,0,100,120,-40,40);
  hist_list_2.push_back(h_Edep_z_3_badN);
  TH2D * h_ToFc_z_1_badN = new TH2D("ToFc_z_1_badN","ToF*c [cm] vs. z [cm] of CND Neutrons;ToF*c;z",50,12,120,50,-40,40);
  hist_list_2.push_back(h_ToFc_z_1_badN);
  TH2D * h_ToFc_z_2_badN = new TH2D("ToFc_z_2_badN","ToF*c [cm] vs. z [cm] of CND Neutrons;ToF*c;z",50,12,120,50,-40,40);
  hist_list_2.push_back(h_ToFc_z_2_badN);
  TH2D * h_ToFc_z_3_badN = new TH2D("ToFc_z_3_badN","ToF*c [cm] vs. z [cm] of CND Neutrons;ToF*c;z",50,12,120,50,-40,40);
  hist_list_2.push_back(h_ToFc_z_3_badN);
  TH1D * h_diff_ToFc_z_1_badN = new TH1D("diff_ToFc_z_1_badN","ToF*c - z [cm] of CND Neutrons;ToF*c-z",100,0,300);
  hist_list_1.push_back(h_diff_ToFc_z_1_badN);
  TH1D * h_diff_ToFc_z_2_badN = new TH1D("diff_ToFc_z_2_badN","ToF*c - z [cm] of CND Neutrons;ToF*c-z",100,0,300);
  hist_list_1.push_back(h_diff_ToFc_z_2_badN);
  TH1D * h_diff_ToFc_z_3_badN = new TH1D("diff_ToFc_z_3_badN","ToF*c - z [cm] of CND Neutrons;ToF*c-z",100,0,300);
  hist_list_1.push_back(h_diff_ToFc_z_3_badN);
  TH2D * h_diff_ToFc_z_Edep_1_badN = new TH2D("diff_ToFc_z_Edep_1_badN","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_1_badN);
  TH2D * h_diff_ToFc_z_Edep_2_badN = new TH2D("diff_ToFc_z_Edep_2_badN","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_2_badN);
  TH2D * h_diff_ToFc_z_Edep_3_badN = new TH2D("diff_ToFc_z_Edep_3_badN","ToF*c - z [cm] vs. E_{dep} [MeV] of CND Neutrons;ToF*c-z;E_{dep}",50,0,300,50,0,100);
  hist_list_2.push_back(h_diff_ToFc_z_Edep_3_badN);
  TH2D * h_Edep_mom_badN = new TH2D("Edep_mom_badN","E_{dep} [MeV] vs. Momentum [GeV] of CND1 Neutrons;E_{dep};Momentum",100,0,100,100,0,1.5);
  hist_list_2.push_back(h_Edep_mom_badN);

  */
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Sumw2();
    hist_list_1[i]->GetXaxis()->CenterTitle();
    hist_list_1[i]->GetYaxis()->CenterTitle();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Sumw2();
    hist_list_2[i]->GetXaxis()->CenterTitle();
    hist_list_2[i]->GetYaxis()->CenterTitle();
  }


  int counter = 0;

  //Define cut class
  //while((chain.Next()==true) &&(counter<10000)){
  while(chain.Next()==true){
      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
      }    

      auto allParticles = c12->getDetParticles();
      auto electrons=c12->getByID(11);
      double weight = 1;
      if(isMC){weight=c12->mcevent()->getWeight();}
      TVector3 	p_b(0,0,Ebeam);
      if(electrons.size()!=1){ continue;}
      TVector3 p_e;
      p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
      double EoP_e =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy()) / p_e.Mag();
      int nphe = electrons[0]->che(HTCC)->getNphe();
      double vtz_e = electrons[0]->par()->getVz();
      if(!myCut.electroncut(c12)){continue;}      
      int esector = electrons[0]->getSector();
  /////////////////////////////////////
  //Electron Kinematics  
  /////////////////////////////////////
      TVector3	p_q = p_b - p_e;
      double theta_q =  p_q.Theta() * 180 / M_PI;
      double nu = Ebeam - p_e.Mag();
      double QSq = p_q.Mag2() - (nu*nu);
      double xB = QSq / (2 * mN * nu);
      double WSq = (mN*mN) - QSq + (2*nu*mN);
      double theta_e = p_e.Theta() * 180 / M_PI;
      //Lead Proton
      int num_L = 0;
      int index_L = -1;
      for(int j = 0; j < allParticles.size(); j ++){
	if((LeadFDProton_Cut(c12,Ebeam,j)) || (LeadCDProton_Cut(c12,Ebeam,j))){
	  num_L++;
	  index_L=j;
	}
      }
      if(num_L!=1){continue;}
      bool LeadCD = LeadCDProton_Cut(c12,Ebeam,index_L);
      bool LeadFD = LeadFDProton_Cut(c12,Ebeam,index_L);
      if(LeadCD&&LeadFD){
	cout<<"Problem!\n";
      }
      TVector3 p_L;
      p_L.SetMagThetaPhi(allParticles[index_L]->getP(),allParticles[index_L]->getTheta(),allParticles[index_L]->getPhi());
      TVector3 p_miss = p_q - p_L;
      double mmiss = get_mmiss(p_b,p_e,p_L);      
      if(p_miss.Theta()*180/M_PI< 40){continue;}
      if(p_miss.Theta()*180/M_PI>135){continue;}
      if(p_miss.Mag()<0.2){continue;}
      if(p_miss.Mag()>1.25){continue;}
      if(mmiss<0.7){continue;}
      if(mmiss>1.2){continue;}
      
      bool match = false;
      //////////////////////////////////////////////////
      //For after checking the hipo banks
      //////////////////////////////////////////////////      
      int num_Charge = 0;
      for(int j = 0; j < allParticles.size(); j++){
	if(j==0){continue;}
	if(j==index_L){continue;}
	//if(j==index_Rp1){continue;}
	if(allParticles[j]->par()->getCharge()==0){continue;}
	num_Charge++;
      }
      if(num_Charge>0){continue;}

      if(LeadFD){
	h_xB_mmiss_epFD->Fill(xB,mmiss,weight);
      }
      else if(LeadCD){
	h_xB_mmiss_epCD->Fill(xB,mmiss,weight);
      }
      h_pmiss_ep->Fill(p_miss.Mag(),weight);
      if(mmiss>1.05){continue;}
      if(LeadCD && (xB<1.2)){continue;}
      //if(LeadFD && (xB<0.8)){continue;}
  /////////////////////////////////////
  //Lead Neutron Checks
  /////////////////////////////////////
      for(int j = 0; j < allParticles.size(); j++){
	if(allParticles[j]->par()->getCharge()!=0){continue;}
	bool CT = (allParticles[j]->sci(clas12::CTOF)->getDetector() == 4);
	bool C1 = (allParticles[j]->sci(clas12::CND1)->getDetector() == 3);
	bool C2 = (allParticles[j]->sci(clas12::CND2)->getDetector() == 3);
	bool C3 = (allParticles[j]->sci(clas12::CND3)->getDetector() == 3);

	if(!(C1||C2||C3)){continue;}
	if(allParticles[j]->getTheta()*180/M_PI > 160){continue;}
	double theta = allParticles[j]->getTheta()*180/M_PI;
	double beta = allParticles[j]->par()->getBeta();
	double gamma = 1/sqrt(1-(beta*beta));
	double mom = gamma * beta * mN;
	double ToF = allParticles[j]->getTime() - c12->event()->getStartTime();

	int detINTlayer = C1?1:C2?2:3;
	auto detlayer = C1?CND1:C2?CND2:CND3;
	double edep = allParticles[j]->sci(CND1)->getEnergy() + allParticles[j]->sci(CND2)->getEnergy() + allParticles[j]->sci(CND3)->getEnergy();
	double edep_CTOF = allParticles[j]->sci(CTOF)->getEnergy();
	double edep_single = allParticles[j]->sci(detlayer)->getEnergy();

	double nvtx_x = allParticles[j]->par()->getVx();
	double nvtx_y = allParticles[j]->par()->getVy();
	double nvtx_z = allParticles[j]->par()->getVz();
	TVector3 v_nvtx(nvtx_x,nvtx_y,nvtx_z);
	TVector3 v_hit;	
	v_hit.SetXYZ(allParticles[j]->sci(detlayer)->getX(),allParticles[j]->sci(detlayer)->getY(),allParticles[j]->sci(detlayer)->getZ());

	TVector3 v_path = v_hit - v_nvtx;
	TVector3 v_n;
	v_n.SetMagThetaPhi(mom,v_path.Theta(),v_path.Phi());
	double path = v_path.Mag()/100;
	double theta_nmiss = v_n.Angle(p_miss)*180/M_PI;
	double dm_nmiss = (p_miss.Mag()-v_n.Mag())/p_miss.Mag();
	int nSector=allParticles[j]->sci(detlayer)->getSector();

	//Check to see if there is a good neutron
	bool isGN = false;
	if((theta_nmiss<40) && (dm_nmiss>-0.5) && (dm_nmiss< 0.5)){
	  isGN=true;
	}
	//////////////////////////////////////////////
	//Step Zero
	//////////////////////////////////////////////
	if(beta-(path*100)/(ToF*c)<-0.01){continue;}
	if(beta-(path*100)/(ToF*c)> 0.01){continue;}
	if(v_hit.Z()>45){continue;}	
	if(v_hit.Z()<-40){continue;}
	if(ToF<0){continue;}
	if(ToF>20){continue;}

	if(LeadFD){
	  h_xB_mmiss_epnFD->Fill(xB,mmiss,weight);
	}
	else if(LeadCD){
	  h_xB_mmiss_epnCD->Fill(xB,mmiss,weight);
	}

	h_pnRes_theta_nmiss_Step0->Fill(dm_nmiss,theta_nmiss,weight);
	
	if(isGN){	
	  if(LeadFD){
	    h_xB_mmiss_epngoodFD->Fill(xB,mmiss,weight);
	  }
	  else if(LeadCD){
	    h_xB_mmiss_epngoodCD->Fill(xB,mmiss,weight);
	  }
	  
	  h_ToF_goodN_Step0->Fill(ToF,weight);
	  h_beta_goodN_Step0->Fill(beta,weight);
	  h_Edep_goodN_Step0->Fill(edep,weight);
	  h_beta_Edep_goodN_Step0->Fill(beta,edep,weight);
	}
	else{
	  h_ToF_badN_Step0->Fill(ToF,weight);
	  h_beta_badN_Step0->Fill(beta,weight);	
	  h_Edep_badN_Step0->Fill(edep,weight);
	  h_beta_Edep_badN_Step0->Fill(beta,edep,weight);
	}

	//////////////////////////////////////////////
	//Step One
	//////////////////////////////////////////////
	if(beta>0.8){continue;}
	if(edep<5){continue;}

	h_pnRes_theta_nmiss_Step1->Fill(dm_nmiss,theta_nmiss,weight);	
	if(isGN){	
	  h_ToF_goodN_Step1->Fill(ToF,weight);
	  h_pmiss_goodN_Step1->Fill(p_miss.Mag(),weight);
	}
	else{
	  h_ToF_badN_Step1->Fill(ToF,weight);
	}

	bool CNDVeto = false;
	if(ToF*c-v_hit.Z() < 70){
	  
	if(isGN){	
	  h_edep_goodN_Step1->Fill(edep,weight);
	}
	else{
	  h_edep_badN_Step1->Fill(edep,weight);
	}

	for(int k = 0; k < allParticles.size(); k++){
	  if(k==0){continue;}
	  if(k==j){continue;}
	  if(allParticles[k]->par()->getCharge()<=0){continue;}
	  if(allParticles[k]->sci(CTOF)->getDetector() == 0){continue;}
	  int vetoSectorbyLayer[4] = {(allParticles[k]->sci(CTOF)->getComponent()+1)/2,allParticles[k]->sci(CND1)->getSector(),allParticles[k]->sci(CND2)->getSector(),allParticles[k]->sci(CND3)->getSector()};
	  double momP = allParticles[k]->getP();
	  double edep_pos = allParticles[k]->sci(clas12::CTOF)->getEnergy();

	  for(int k = 0; k < 4; k++){	    
	    if(vetoSectorbyLayer[k]==0){continue;}
	    int sdiff = nSector - vetoSectorbyLayer[k];
	    if(sdiff<=-12){sdiff+=24;}
	    else if(sdiff>12){sdiff-=24;}
	    int ldiff = detINTlayer - k;

	    if(isGN){	      
	      h_sdiff_pos_goodN_Step1_layer[ldiff+3]->Fill(sdiff,weight);
	      h_sdiff_pos_mom_goodN_Step1_layer[ldiff+3]->Fill(sdiff,momP,weight);
	      h_sdiff_pos_z_goodN_Step1_layer[ldiff+3]->Fill(sdiff,v_hit.Z(),weight);
	      h_sdiff_pos_diff_ToFc_z_goodN_Step1_layer[ldiff+3]->Fill(sdiff,ToF*c-v_hit.Z(),weight);
	    }
	    else{
	      h_sdiff_pos_badN_Step1_layer[ldiff+3]->Fill(sdiff,weight);
	      h_sdiff_pos_mom_badN_Step1_layer[ldiff+3]->Fill(sdiff,momP,weight);
	      h_sdiff_pos_z_badN_Step1_layer[ldiff+3]->Fill(sdiff,v_hit.Z(),weight);
	      h_sdiff_pos_diff_ToFc_z_badN_Step1_layer[ldiff+3]->Fill(sdiff,ToF*c-v_hit.Z(),weight);
	    }
	    if(isPosNear(sdiff,ldiff)){
	      CNDVeto=true;
	    }
	  }

	if(CNDVeto){
	if(isGN){	
	  h_edep_over_edepCTOT_goodN_Step1->Fill(edep/edep_pos,weight);
	}
	else{
	  h_edep_over_edepCTOT_badN_Step1->Fill(edep/edep_pos,weight);
	}
	}
	}

	if(CNDVeto){
	if(isGN){	
	  h_edep_goodN_withNearbyPos_Step1->Fill(edep,weight);
	}
	else{
	  h_edep_badN_withNearbyPos_Step1->Fill(edep,weight);
	}
	}

	if(isGN){	
	  if(!CNDVeto)
	    h_diff_ToFc_z_Edep_noNear_goodN_Step1->Fill(ToF*c-v_hit.Z(),edep,weight);
	  else{
	    h_diff_ToFc_z_Edep_yesNear_goodN_Step1->Fill(ToF*c-v_hit.Z(),edep,weight);	  
	  }
	}
	else{
	  if(!CNDVeto)
	    h_diff_ToFc_z_Edep_noNear_badN_Step1->Fill(ToF*c-v_hit.Z(),edep,weight);
	  else{
	    h_diff_ToFc_z_Edep_yesNear_badN_Step1->Fill(ToF*c-v_hit.Z(),edep,weight);	  
	  }
	}
	}
	//////////////////////////////////////////////
	//Step Two
	//////////////////////////////////////////////
	if(CNDVeto){continue;}
	h_pnRes_theta_nmiss_Step2->Fill(dm_nmiss,theta_nmiss,weight);	
	if(isGN){	
	  h_ToF_goodN_Step2->Fill(ToF,weight);
	}
	else{
	  h_ToF_badN_Step2->Fill(ToF,weight);
	}
	for(int k = 0; k < allParticles.size(); k++){
	  if(k==0){continue;}
	  if(k==j){continue;}
	  if(allParticles[k]->par()->getCharge()<=0){continue;}
	  if(allParticles[k]->sci(CTOF)->getDetector() == 0){continue;}
	  int vetoSectorbyLayer[4] = {(allParticles[k]->sci(CTOF)->getComponent()+1)/2,allParticles[k]->sci(CND1)->getSector(),allParticles[k]->sci(CND2)->getSector(),allParticles[k]->sci(CND3)->getSector()};

	  for(int k = 0; k < 4; k++){	    
	    if(vetoSectorbyLayer[k]==0){continue;}
	    int sdiff = nSector - vetoSectorbyLayer[k];
	    if(sdiff<=-12){sdiff+=24;}
	    else if(sdiff>12){sdiff-=24;}
	    int ldiff = detINTlayer - k;

	    if(isGN){
	      h_sdiff_pos_goodN_Step2_layer[ldiff+3]->Fill(sdiff,weight);
	    }
	    else{
	      h_sdiff_pos_badN_Step2_layer[ldiff+3]->Fill(sdiff,weight);
	    }
	  }	    	  
	}
	
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

	//if(C3 && (v_hit.Z()>25)){continue;}
	//else if(C2 && (v_hit.Z()>20)){continue;}
	//else if(C1 && (v_hit.Z()>10)){continue;}
      }
  }
  cout<<counter<<endl;

  TH1D * h_pmissrat_goodN_Step1 = (TH1D*)h_pmiss_goodN_Step1->Clone("pmissrat_goodN_Step1");
  h_pmissrat_goodN_Step1->Divide(h_pmiss_ep);
  hist_list_1.push_back(h_pmissrat_goodN_Step1);  

  TH1D * h_pmissrat_goodN_Step5 = (TH1D*)h_pmiss_goodN_Step5->Clone("pmissrat_goodN_Step5");
  h_pmissrat_goodN_Step5->Divide(h_pmiss_ep);
  hist_list_1.push_back(h_pmissrat_goodN_Step5);  
  
  TH1D * h_pmissrat_allN_Step5 = (TH1D*)h_pmiss_allN_Step5->Clone("pmissrat_allN_Step5");
  h_pmissrat_allN_Step5->Divide(h_pmiss_ep);
  hist_list_1.push_back(h_pmissrat_allN_Step5);  
  
  outFile->cd();
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Write();
  }


  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);

  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);


  /////////////////////////////////////
  //CND Neutron Information
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Neutrons in CND");
  myText->Print(fileName,"pdf");
  myText->Clear();
  
  myCanvas->Divide(1,1);
  for(int i = 0; i < hist_list_1.size(); i++){
    myCanvas->cd(1);
    hist_list_1[i]->Draw();
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  for(int i = 0; i < hist_list_2.size(); i++){
    myCanvas->cd(1);
    hist_list_2[i]->Draw("colz");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  outFile->Close();
}


void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}
