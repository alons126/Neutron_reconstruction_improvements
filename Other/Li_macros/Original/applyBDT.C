void applyBDT() {
	TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
	Float_t nu_SciCD_energy0, nu_SciCD_energy1, nu_SciCD_energy2, nu_SciCD_energy3;
	Float_t nu_SciCD_size0, nu_SciCD_size1, nu_SciCD_size2, nu_SciCD_size3;
	Float_t CVT_nTrack, CTOF_nCluster, CND_nTrack;
	Float_t CVT_deltaAngle, CTOF_deltaAngle;
	Float_t nu_SciCD_energy[4];
	Int_t nu_SciCD_size[4];
	Int_t CVTnTrack, CTOFnCluster, CNDnTrack;
	Float_t CVTDeltaAngle, CTOFDeltaAngle;
	Float_t responseBDT;

	reader->AddVariable("nu_SciCD_energy[0]", &nu_SciCD_energy0);
	reader->AddVariable("nu_SciCD_energy[1]", &nu_SciCD_energy1);
	reader->AddVariable("nu_SciCD_energy[2]", &nu_SciCD_energy2);
	reader->AddVariable("nu_SciCD_energy[3]", &nu_SciCD_energy3);
	reader->AddVariable("nu_SciCD_size[0]", &nu_SciCD_size0);
	reader->AddVariable("nu_SciCD_size[1]", &nu_SciCD_size1);
	reader->AddVariable("nu_SciCD_size[2]", &nu_SciCD_size2);
	reader->AddVariable("nu_SciCD_size[3]", &nu_SciCD_size3);
	reader->AddVariable("CVT_nTrack", &CVT_nTrack);
	reader->AddVariable("CVT_deltaAngle", &CVT_deltaAngle);
	reader->AddVariable("CTOF_nCluster", &CTOF_nCluster);
	reader->AddVariable("CTOF_deltaAngle", &CTOF_deltaAngle);
	reader->BookMVA("BDT", "dataset/weights/trainBDT_BDT.weights.xml");

	TFile* myFile = new TFile("dataDVCS_forTMVA.root");
	TTree* myTree = (TTree*) myFile->Get("nDVCSTree");
	myTree->SetBranchAddress("nu_SciCD_energy", &nu_SciCD_energy);
	myTree->SetBranchAddress("nu_SciCD_size", &nu_SciCD_size);
	myTree->SetBranchAddress("CVT_nTrack", &CVTnTrack);
	myTree->SetBranchAddress("CVT_deltaAngle", &CVTDeltaAngle);
	myTree->SetBranchAddress("CTOF_nCluster", &CTOFnCluster);
	myTree->SetBranchAddress("CTOF_deltaAngle", &CTOFDeltaAngle);

	TFile* newFile = new TFile("dataDVCS_BDT.root", "recreate");
	TTree* newTree = myTree->CopyTree("");
	TBranch* bResponseBDT = newTree->Branch("BDT_response", &responseBDT, "BDT_response/F");
	for (Long64_t ievt=0; ievt<myTree->GetEntries(); ievt++) {
		myTree->GetEntry(ievt);
		nu_SciCD_energy0 = nu_SciCD_energy[0];
		nu_SciCD_energy1 = nu_SciCD_energy[1];
		nu_SciCD_energy2 = nu_SciCD_energy[2];
		nu_SciCD_energy3 = nu_SciCD_energy[3];
		nu_SciCD_size0 = nu_SciCD_size[0];
		nu_SciCD_size1 = nu_SciCD_size[1];
		nu_SciCD_size2 = nu_SciCD_size[2];
		nu_SciCD_size3 = nu_SciCD_size[3];
		CVT_nTrack = CVTnTrack;
		CTOF_nCluster = CTOFnCluster;
		CVT_deltaAngle = CVTDeltaAngle;
		CTOF_deltaAngle = CTOFDeltaAngle;
		responseBDT = reader->EvaluateMVA("BDT");
		bResponseBDT->Fill();
	}
	newTree->Write();
	newFile->Close();
}
