#include <TMVA/DataLoader.h>
#include <TMVA/Factory.h>

/* Below, anything in multiline comments is my addition! */

void trainBDT() {
	/* A file for some plots, such as ROC curve and classifier output distribution. */
	TFile* outputFile = new TFile("TMVABDT.root", "recreate");

	/* A factory class that handles training and testing of the classifiers. */
	TMVA::Factory *factory = new TMVA::Factory( "trainBDT", outputFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

	/* A class is responsible for handling the loading and preprocessing of data used in machine learning algorithms. */
	TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

	TFile* sigSrc = new TFile("mcDVCS_rgbNeutronForTMVA.root");
	TFile* bkgSrc = new TFile("mcDVCS_rgbProtonForTMVA.root");
	TTree* sigTree = (TTree*)sigSrc->Get("nDVCSTree");
	TTree* bkgTree = (TTree*)bkgSrc->Get("nDVCSTree");
	std::cout << sigTree->GetEntries() << std::endl;
	std::cout << bkgTree->GetEntries() << std::endl;

	Double_t sigWeight = 1.0;
	Double_t bkgWeight = 1.0;
	dataloader->AddVariable("nu_SciCD_energy[0]", 'F');
	dataloader->AddVariable("nu_SciCD_energy[1]", 'F');
	dataloader->AddVariable("nu_SciCD_energy[2]", 'F');
	dataloader->AddVariable("nu_SciCD_energy[3]", 'F');
	dataloader->AddVariable("nu_SciCD_size[0]", 'I');
	dataloader->AddVariable("nu_SciCD_size[1]", 'I');
	dataloader->AddVariable("nu_SciCD_size[2]", 'I');
	dataloader->AddVariable("nu_SciCD_size[3]", 'I');
	dataloader->AddVariable("CVT_nTrack", 'I');
	dataloader->AddVariable("CVT_deltaAngle", 'F');
	dataloader->AddVariable("CTOF_nCluster", 'I');
	dataloader->AddVariable("CTOF_deltaAngle", 'F');
	dataloader->AddSignalTree(sigTree, sigWeight);
	dataloader->AddBackgroundTree(bkgTree, bkgWeight);

	/* Cuts on variables from the TTrees. */
	TCut mycuts = "nu_status>3000";
	TCut mycutb = "nu_status>3000";

	dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "SplitMode=random:!V" );

	factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
                        "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();

	outputFile->Close();
}
