void trainBDT() {
	TFile* outputFile = new TFile("TMVABDT.root", "recreate");
	TMVA::Factory *factory = new TMVA::Factory( "trainBDT", outputFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
	TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

	TFile* sigSrc = new TFile("mcDVCS_rgbNeutronForTMVA.root");
	TFile* bkgSrc = new TFile("mcDVCS_rgbProtonForTMVA.root");
	TTree* sigTree = (TTree*)sigSrc->Get("nDVCSTree");
	TTree* bkgTree = (TTree*)bkgSrc->Get("nDVCSTree");
	cout << sigTree->GetEntries() << endl;
	cout << bkgTree->GetEntries() << endl;

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
