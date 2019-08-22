void runStarlight() {
  
  TClonesArray *fParticles = new TClonesArray("TParticle", 200);
  TTree *fEventTree = new TTree("fEventTree", "fEventTree");
  fEventTree ->Branch("fParticles", &fParticles);

  //TStarlight
  gSystem->Load("libStarLight.so");
  gSystem->Load("libAliStarLight.so");
  
  TClonesArray *fSLparticles = new TClonesArray("TParticle", 200);
  TStarLight *fStarlight = new TStarLight("fStarlight","fStarlight","./slight.in");
  fStarlight->InitStarLight();
  TMap *fSLparams = new TMap();
  
  //Neutron generator
  gROOT->LoadMacro("NeutronGenerator.cxx+g");
  
  TClonesArray *fNGparticles = NULL;
  NeutronGenerator *fNeutronGen = new NeutronGenerator();
  fNeutronGen->SetHadronicInteractionModel(NeutronGenerator::kHardSphere);
  fNeutronGen->Initialize();
  fNeutronGen->LoadENDF(); 
  fNeutronGen->SetRunMode(NeutronGenerator::kInterface);
  fNeutronGen->SetStoreQA();
  fNeutronGen->SetStoreGeneratorFunctions();
  fNeutronGen->Setup();
  
  Double_t VMrapidity = 0;
  Double_t VMmass = 3.09;
  Double_t photonK = 0;
  Double_t fRapMin = 2.5;
  Double_t fRapMax = 4.0;
  
  //Run
  UInt_t nEvents = 10000;
  cout<<"Running production"<<endl; 
  for(Int_t iEvent = 0; iEvent<=nEvents; iEvent++){
    if(iEvent%(nEvents/10) == 0){
      if(iEvent != 0){ printf("\033[1A"); printf("\033[K");}   
      cout<<100*iEvent/nEvents<<"% "<<endl;
      }
    Int_t nTotalPart = 0;  
    
    Bool_t genOK = kFALSE;
    while(!genOK){  
      fStarlight->GenerateEvent();
      fStarlight->BoostEvent();
      fStarlight->ImportParticles(fSLparticles, "ALL");
      fStarlight->ImportEventInfo(fSLparams);
    
      const TParameter<double>* p = dynamic_cast<const TParameter<double>* >(fSLparams->GetValue("Egam"));
      photonK = p->GetVal();
    
      VMrapidity = TMath::Abs(TMath::Log(2*photonK/3.09));
      if(VMrapidity<fRapMax && VMrapidity>fRapMin) genOK = kTRUE;
      }
      
    for(Int_t i = 0; i<fSLparticles->GetEntriesFast(); i++){
      TParticle *part(dynamic_cast<TParticle*>(fSLparticles->At(i)));
      new((*fParticles)[nTotalPart++]) TParticle(*part);
      }
      
    fNeutronGen->GenerateEvent(photonK);
    fNGparticles = fNeutronGen->ImportParticles();
    for(Int_t i = 0; i<fNGparticles->GetEntriesFast(); i++){
      TParticle *part(dynamic_cast<TParticle*>(fNGparticles->At(i)));
      new((*fParticles)[nTotalPart++]) TParticle(*part);
      }
    fNeutronGen->FinishEvent();
    
    fEventTree->Fill();
    fParticles->Clear("C");
    }
  fNeutronGen->FinishProduction();
  
  TFile *fOutputFile = new TFile("SLoutput.root","RECREATE");
  fEventTree->Write();
  fOutputFile->Close();  
}
