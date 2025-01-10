// --- Standard library ---
#include <iostream>
#include <fstream>
#include <string>

// --- ROOT system ---
#include <TSystem.h>
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"

#include "NeutronGenerator.h"

ClassImp(NeutronGenerator)

using std::ios;
using std::cout;
using std::endl;

//_____________________________________________________________________________
NeutronGenerator::NeutronGenerator()
  : TObject()
  , fRunMode(kInterface)
  , fHadronicIntModel(kGlauber)
  , fProductionMode(kPhotonPomeron)
  , fNucleus(kPb208)
  , iEvent(0)
  , hInputRapidity(NULL)
  , hInputMass(NULL)
  , fRapMin(-666)
  , fRapMax(666)
  , fMassMin(1)
  , fMassMax(666)
  , lineString(0)
  , fDataPath("./Data/")
  , nFluxes(2+(maxNeutrons)*(maxNeutrons+1)/2)
  , nucleus_Z(82)
  , nucleus_A(208)
  , nucleus_R(6.624)
  , beamGamma(0)
  , gammaTarget(2.0*beamGamma*beamGamma-1.0)
  , neutronSepThr(0.0)
  , saturationEnergy(1e6)
  , gXsection(NULL)
  , hXsection(NULL)
  , gPhotonFluxTable(NULL)
  , gNucleusBreakupTable(NULL)
  , gTwoPhotonFluxTable(NULL)
  , gNuclearThickness(NULL)
  , gHadronicProbabilityTable(NULL)
  , gMeanNeutronN(NULL)
  , gWidthNeutronN(NULL)
  , fitMeanNeutronN(NULL)
  , fitWidthNeutronN(NULL)
  , hENDF_2D(NULL)
  , hENDF_1D(NULL)
  , fENDFFile(NULL)
  , hBranchingRatioMap(NULL)
  , hEventBreakupMap(new TH1D("hEventBreakupMap","",nFluxes,0,nFluxes))
  , gUnitaryLeak(new TGraph)
  , fEventTree(NULL)
  , fParticles(NULL)
  , fOutputFile(NULL)
  , fImpProfile(NULL) 
  , fQAhistList(NULL)
  , hNeutronMultiplicity(NULL)
  , hEnergyGen(NULL)
  , hKinEnergyGen(NULL)
  , hMomGen(NULL)
  , hPhiGen(NULL)
  , hThetaGen(NULL)
  , hEnergyBoosted(NULL)
  , hMomBoosted(NULL)
  , hPhiBoosted(NULL)
  , hThetaBoosted(NULL)
  , hNeutronRapidity(NULL)
  , hNeutronEta(NULL)
  , hEnergyBin(NULL)
  , hEnergyForNeutronMulti(NULL) 
  , hProbabilityXn(NULL) 
  , kStoreQA(kFALSE)
  , kStoreGen(kFALSE)
{
 // Default constructor
}
//_____________________________________________________________________________
NeutronGenerator::~NeutronGenerator()
{
  //destructor
  if(hInputRapidity) {delete hInputRapidity; hInputRapidity = NULL;}
  if(hInputMass) {delete hInputMass; hInputMass = NULL;}
  if(gPhotonFluxTable) {delete gPhotonFluxTable; gPhotonFluxTable = NULL;}
  if(gTwoPhotonFluxTable) {delete gTwoPhotonFluxTable; gTwoPhotonFluxTable = NULL;}
  if(gNuclearThickness) {delete gNuclearThickness; gNuclearThickness = NULL;}
  if(gHadronicProbabilityTable) {delete gHadronicProbabilityTable; gHadronicProbabilityTable = NULL;}
  if(gMeanNeutronN) {delete gMeanNeutronN; gMeanNeutronN = NULL;}
  if(gWidthNeutronN) {delete gWidthNeutronN; gWidthNeutronN = NULL;}
  if(fitMeanNeutronN) {delete fitMeanNeutronN; fitMeanNeutronN = NULL;}
  if(fitWidthNeutronN) {delete fitWidthNeutronN; fitWidthNeutronN = NULL;}
  if(hENDF_2D) {delete hENDF_2D; hENDF_2D = NULL;}
  if(hENDF_1D) {delete hENDF_1D; hENDF_1D = NULL;}
  if(fENDFFile) {delete fENDFFile; fENDFFile = NULL;} 
  if(fEventTree) {delete fEventTree; fEventTree = NULL;}
  if(fParticles) {delete fParticles; fParticles = NULL;}
  if(fOutputFile) {delete fOutputFile; fOutputFile = NULL;}
  if(fQAhistList){delete fQAhistList; fQAhistList = NULL;}
  
}
//______________________________________________________________________________
void NeutronGenerator::SetRunMode(RunMode_t mode,const char *filename, const char *name1, const char *name2)
{
  fRunMode = mode;
  
  if(fRunMode == kMassRapidity){
    TFile *inputFile = new TFile(filename,"READ");
    hInputMass = (TH1*)inputFile->Get(name1);
    hInputMass->SetDirectory(0);
    hInputRapidity = (TH1*)inputFile->Get(name2);
    hInputRapidity->SetDirectory(0);
    inputFile->Close();
    delete inputFile;
    }
    
  if(fRunMode == kStarlightAscii){
    //fInputStarlightAscii.open(filename);
    //std::string newLine;
    //getline(fInputStarlightAscii,newLine);
    //fOutputStarlightAscii.open(name1);
    //fOutputStarlightAscii<<newLine<<endl;
    }
   if(fRunMode == kFlatMultiplicity){}
   if(fRunMode == kInterface){}
   if(fRunMode == k1n1n){}
}
//______________________________________________________________________________
void NeutronGenerator::Setup(){

 fParticles = new TClonesArray("TParticle", 200);
  if(fRunMode == kMassRapidity || fRunMode == kFlatMultiplicity || fRunMode == k1n1n){
    fEventTree = new TTree("fEventTree", "fEventTree");
    fEventTree ->Branch("fParticles", &fParticles);
    fParticles ->BypassStreamer();
    }
  if(fRunMode == kStarlightAscii){}
  if(fRunMode == kInterface){}
  if(fRunMode == k1n1n){}
}
//______________________________________________________________________________
void NeutronGenerator::Run(const UInt_t nEvents)
{
  TDatabasePDG pdgData;
  if(fRunMode == kInterface){cout<<"This method is obsolette in interface mode"; return;}
  Double_t VMrapidity = 0;
  Double_t VMmass = 0;
  Double_t photonK = 0;
  //
  cout<<"Running production"<<endl; 
  for(iEvent = 0; iEvent<=nEvents; iEvent++){
    if(iEvent%(nEvents/10) == 0){
      if(iEvent != 0){ printf("\033[1A"); printf("\033[K");}   
      cout<<100*iEvent/nEvents<<"% "<<endl;
      }
      
    if(fRunMode == kMassRapidity){
      do VMmass = hInputMass->GetRandom(); while(VMmass>fMassMax || VMmass<fMassMin);
      do VMrapidity = hInputRapidity->GetRandom(); while(!((VMrapidity<fRapMax && VMrapidity>fRapMin) || (VMrapidity> -1*fRapMax && VMrapidity< -1*fRapMin)));
      hMassVM->Fill(VMmass);
      hRapidityVM->Fill(VMrapidity);
      if(VMrapidity<fRapMax && VMrapidity>fRapMin){ 
        photonK = 0.5*VMmass*TMath::Exp(-1*VMrapidity);
	}
      else if(VMrapidity> -1*fRapMax && VMrapidity< -1*fRapMin){
        photonK = 0.5*VMmass*TMath::Exp(VMrapidity);
	}
      }
    
    if(fRunMode == kFlatMultiplicity || fRunMode == k1n1n) photonK = -1;

    GenerateEvent(photonK);
    FinishEvent();
    }
    //
  FinishProduction();
}
//______________________________________________________________________________
void NeutronGenerator::GenerateEvent(const Double_t photonK1, const Double_t photonK2)
{
  if(fRunMode != kFlatMultiplicity && fRunMode != k1n1n){
    hPhotonK->Fill(photonK1);
    if(photonK2>0)hPhotonK->Fill(photonK2);
    }
  Int_t nNeutronsBeam1 = 0, nNeutronsBeam2 = 0;
  
  if(fRunMode == kFlatMultiplicity){
    nNeutronsBeam1 = maxNeutrons*gRandom->Rndm(); 
    nNeutronsBeam2 = maxNeutrons*gRandom->Rndm();
    }
  else if(fRunMode == k1n1n){
    nNeutronsBeam1 = 1; 
    nNeutronsBeam2 = 1;
    }
  else{
    if(fProductionMode == kPhotonPomeron){ 
      for(Int_t i = 0; i < maxNeutrons; i++){ 
        for(Int_t j = i; j < maxNeutrons; j++){
          if(i != j)hEventBreakupMap->SetBinContent(FromMatrixToVector(i,j)+1,2.0*GetBreakupProbability(photonK1,i,j));
	  if(i == j)hEventBreakupMap->SetBinContent(FromMatrixToVector(i,j)+1,GetBreakupProbability(photonK1,i,j));
	  }
        }
      }
    if(fProductionMode == kTwoPhoton){
      for(Int_t i = 0; i < maxNeutrons; i++){ 
        for(Int_t j = i; j < maxNeutrons; j++){
          if(i != j)hEventBreakupMap->SetBinContent(FromMatrixToVector(i,j)+1,2.0*GetBreakupProbability(photonK1,photonK2,i,j));
          if(i == j)hEventBreakupMap->SetBinContent(FromMatrixToVector(i,j)+1,GetBreakupProbability(photonK1,photonK2,i,j));
          }
        }
      }
      
    Int_t randBin = hEventBreakupMap->FindBin(hEventBreakupMap->GetRandom())-1;
    FromVectorToMatrix(randBin,nNeutronsBeam1,nNeutronsBeam2);
    }
    
  if(nNeutronsBeam1 != nNeutronsBeam2 && gRandom->Rndm()<0.5)CreateNeutrons(nNeutronsBeam2,nNeutronsBeam1);
  else CreateNeutrons(nNeutronsBeam1,nNeutronsBeam2);
}
//______________________________________________________________________________
Double_t NeutronGenerator::GetBreakupProbability(const Double_t photonK, const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2)
{
  Double_t probability = 0;  
  if(nNeutronsBeam1+nNeutronsBeam2 == -2) probability = gPhotonFluxTable[nFluxes-1].Eval(photonK);
  if(nNeutronsBeam1+nNeutronsBeam2 == -1) probability = 1-gPhotonFluxTable[nFluxes-1].Eval(photonK)-gPhotonFluxTable[FromMatrixToVector(0,0)+1].Eval(photonK);
  if(nNeutronsBeam1 >= 0 && nNeutronsBeam2 >= 0) probability = gPhotonFluxTable[FromMatrixToVector(nNeutronsBeam1,nNeutronsBeam2)+1].Eval(photonK);
  			   
  return probability;
}
//______________________________________________________________________________
Double_t NeutronGenerator::GetBreakupProbability(const Double_t photonK1, const Double_t photonK2, const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2, const Bool_t interpolate)
{

  Double_t probability = 0;
  if(interpolate){
    if(nNeutronsBeam1+nNeutronsBeam2 == -2)probability = gTwoPhotonFluxTable[nFluxes-1].Interpolate(photonK1,photonK2);
    if(nNeutronsBeam1+nNeutronsBeam2 == -1)probability = 1-gTwoPhotonFluxTable[nFluxes-1].Interpolate(photonK1,photonK2)-gTwoPhotonFluxTable[FromMatrixToVector(0,0)+1].Interpolate(photonK1,photonK2);
    if(nNeutronsBeam1 >= 0 && nNeutronsBeam2 >= 0) probability = gTwoPhotonFluxTable[FromMatrixToVector(nNeutronsBeam1,nNeutronsBeam2)+1].Interpolate(photonK1,photonK2);
    }
  else{
    Double_t* flux = TwoPhotonFlux(photonK1,photonK2);
    
    if(nNeutronsBeam1+nNeutronsBeam2 == -2) probability = flux[nFluxes-1];
    if(nNeutronsBeam1+nNeutronsBeam2 == -1) probability = 1-flux[nFluxes-1]-flux[FromMatrixToVector(0,0)+1];
    if(nNeutronsBeam1 >= 0 && nNeutronsBeam2 >= 0) probability = flux[FromMatrixToVector(nNeutronsBeam1,nNeutronsBeam2)+1];
    delete [] flux; 
    }
		   
  return probability;
}

//______________________________________________________________________________
Double_t NeutronGenerator::GetTotalFlux(const Double_t photonK1, const Double_t photonK2)
{
  if(fProductionMode == kPhotonPomeron) return gPhotonFluxTable[0].Eval(photonK1);
  if(fProductionMode == kTwoPhoton) return gTwoPhotonFluxTable[0].Interpolate(photonK1,photonK2);
  return -666;
}
//______________________________________________________________________________
void NeutronGenerator::FinishProduction(){

  if(fRunMode == kMassRapidity || fRunMode == kFlatMultiplicity || fRunMode == k1n1n || kStoreQA || kStoreGen)fOutputFile = new TFile("Output.root","RECREATE");
  if(kStoreQA)fQAhistList->Write();
  if(kStoreGen){
    if(gNuclearThickness)gNuclearThickness->Write();
    gHadronicProbabilityTable->Write();
    for(Int_t i=0; i<=10; i++){
      gXsection[i].SetName(Form("gXsection%d",i));
      gXsection[i].Write(); 
      }
    fitMeanNeutronN->Write();
    fitWidthNeutronN->Write();
    gNucleusBreakupTable[maxNeutrons].Write();
    for(Int_t i = 0; i<10; i++) gNucleusBreakupTable[i].Write();
    gUnitaryLeak->SetName("gUnitaryLeak");
    gUnitaryLeak->Write();
    hBranchingRatioMap->Write();
    
    for(Int_t i=0; i<=9; i++){
      fImpProfile[i].SetName(Form("fImpProfile%d",i));
      fImpProfile[i].Write(); 
      }
    
    }
  if(fRunMode == kMassRapidity || fRunMode == kFlatMultiplicity || fRunMode == k1n1n) fEventTree->Write();
  if(fOutputFile)fOutputFile->Close();  
}
//______________________________________________________________________________
void NeutronGenerator::FinishEvent(){
  
  if(fRunMode == kMassRapidity || fRunMode == kFlatMultiplicity || fRunMode == k1n1n)fEventTree ->Fill();
  if(fRunMode == kStarlightAscii){
    //for(Int_t i=0; i< fParticles->GetEntriesFast(); i++)fOutputStarlightAscii<<"TRACK: "<<i+11<<" "<<
    							      // ((TParticle*)fParticles->At(i))->Px()<<" "<<
    							      // ((TParticle*)fParticles->At(i))->Py()<<" "<<
							      // ((TParticle*)fParticles->At(i))->Pz()<<" "<<iEvent+1<<" 0 0 "<<"2112"<<endl;
    //fOutputStarlightAscii<<lineString.Data()<<endl;
    }
  fParticles->Clear("C");
}
//______________________________________________________________________________
void NeutronGenerator::CreateNeutrons(const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2)
{
  hNeutronMultiplicity->Fill(nNeutronsBeam1,nNeutronsBeam2);
  Double_t energyKin = 0;
  Double_t mom = 0,phi = 0, theta = 0;
  TLorentzVector vec;
  Int_t nNeutrons[2] = {nNeutronsBeam1,nNeutronsBeam2};
  Double_t energyPhoton = -1;
  Int_t energyBin = -1;
  
  Int_t nGenerated=0;
  for(Int_t side = 0; side<=1; side++){
  
    if(nNeutrons[side] == 0)continue;
    if(nNeutrons[side] <= 10){
      energyPhoton = hXsection[nNeutrons[side]].GetRandom();
      energyBin = hENDF_2D->GetXaxis()->FindBin(energyPhoton);
      if(energyPhoton>140) energyBin--;
      hENDF_1D = hENDF_2D->ProjectionY("hENDF_1D",energyBin,energyBin);
      hEnergyBin->Fill(energyBin);
      hEnergyForNeutronMulti->Fill(energyPhoton);
      }
    if(nNeutrons[side] > 10)hENDF_1D = hENDF_2D->ProjectionY("hENDF_1D",hENDF_2D->GetNbinsX(),hENDF_2D->GetNbinsX()); 
       
    for(Int_t i = 0; i<nNeutrons[side]; i++){
      energyKin = hENDF_1D->GetRandom();
      mom = TMath::Sqrt((energyKin + neutron_M)*(energyKin + neutron_M) - neutron_M*neutron_M);
      phi = 2*pi*gRandom->Rndm();
      theta = pi*gRandom->Rndm();

      hKinEnergyGen->Fill(energyKin);
      vec.SetXYZM(mom*TMath::Sin(theta)*TMath::Cos(phi),mom*TMath::Sin(theta)*TMath::Sin(phi), mom*TMath::Cos(theta), neutron_M);
      hEnergyGen->Fill(vec.Energy());
      hMomGen->Fill(vec.P());
      hPhiGen->Fill(vec.Phi());
      hThetaGen->Fill(vec.Theta());
      
      vec.Boost(0,0,TMath::Power(-1,side)*TMath::Sqrt(1.0-1.0/beamGamma/beamGamma));
      hEnergyBoosted->Fill(vec.Energy());
      hMomBoosted->Fill(vec.P());
      hPhiBoosted->Fill(vec.Phi());
      hThetaBoosted->Fill(vec.Theta());

      //TParticle units are GeV, up to now we were in MeV
      TParticle *part = (TParticle*) fParticles->ConstructedAt(nGenerated++);
      part->SetMomentum(vec.Px()*0.001, vec.Py()*0.001, vec.Pz()*0.001, vec.Energy()*0.001);
      part->SetPdgCode(2112);
      part->SetStatusCode(1);
      part->SetProductionVertex(0,0,0,0);
      part->SetMother(0,-1);
      part->SetMother(1,-1);
      part->SetDaughter(0,-1);
      part->SetDaughter(1,-1);
      
      hNeutronRapidity->Fill(part->Y());
      hNeutronEta->Fill(part->Eta());
      }
    }
}
//______________________________________________________________________________
void NeutronGenerator::BuildTwoPhotonFluxModulationTables()
{

  Double_t energy_min = 0.5*fMassMin*TMath::Exp(-TMath::Max(TMath::Abs(fRapMin), TMath::Abs(fRapMax))); 
  Double_t energy_max = 0.5*fMassMax*TMath::Exp(TMath::Max(TMath::Abs(fRapMin), TMath::Abs(fRapMax)));  
  Double_t energy_delta = TMath::Exp(TMath::Log(energy_max/energy_min)/Double_t(nSteps_GG));
  
  Int_t iStep = 0;
  
  gTwoPhotonFluxTable = new TGraph2D[nFluxes];
  gTwoPhotonFluxTable[0].SetName("gTwoPhotonFluxTable");
  for(Int_t i = 0; i < maxNeutrons; i++){
    for(Int_t j = i; j < maxNeutrons; j++){
      gTwoPhotonFluxTable[FromMatrixToVector(i,j)+1].SetName(TString::Format("gTwoPhotonFluxTable%d%d",i,j));
      }
    }
    
  cout<<"Building two photon flux modulation table in range Mass = ("<<TString::Format("%.1f",fMassMin)<<","<<TString::Format("%.1f",fMassMax)<<") Rapidity = ("<<TString::Format("%.1f",fRapMin)<<","<<TString::Format("%.1f",fRapMax)<<")"<<endl;  
  for(Double_t energy1_value = energy_min; energy1_value<=energy_max; energy1_value *= energy_delta){
    for(Double_t energy2_value = energy_min; energy2_value<=energy_max; energy2_value *= energy_delta){
      Double_t* flux = TwoPhotonFlux(energy1_value,energy2_value);
      gTwoPhotonFluxTable[0].SetPoint(iStep,energy1_value,energy2_value,flux[0]);
      //cout<<energy1_value<<" "<<energy2_value<<" "<<flux[nFluxes-1]<<endl;
      for(Int_t i = 0; i < maxNeutrons; i++){
        for(Int_t j = i; j < maxNeutrons; j++){
          gTwoPhotonFluxTable[FromMatrixToVector(i,j)+1].SetPoint(iStep,energy1_value,energy2_value,flux[FromMatrixToVector(i,j)+1]);
	  }
        }
      gTwoPhotonFluxTable[nFluxes-1].SetPoint(iStep,energy1_value,energy2_value,flux[nFluxes-1]);
      delete [] flux;
      iStep++;
      if(iStep%12 == 0){
        if(iStep > 12){ printf("\033[1A"); printf("\033[K");}   
        cout<<Int_t(iStep/1.2)<<"%"<<endl;
        }
      }
    }
  for(Int_t iFlux=0; iFlux<nFluxes; iFlux++)gTwoPhotonFluxTable[iFlux].GetHistogram("empty");
}
//______________________________________________________________________________
Double_t *NeutronGenerator::TwoPhotonFlux(const Double_t energyPhoton1, const Double_t energyPhoton2) 
{
  Double_t integral_Phi[nFluxes];
  Double_t integral_Beam2[nFluxes];
  Double_t diferential_Phi = 0;
  Double_t diferential_Density = 0;
  
  Double_t *flux_integral = new Double_t[nFluxes];
  for(Int_t iFlux=0; iFlux<nFluxes; iFlux++)flux_integral[iFlux] = 0.0;
  Double_t breakupProb[maxNeutrons+1];
  
  Double_t impactPar_min=0.8*nucleus_R;
  Double_t impactPar1_max=impactPar_min + 6.0*hbarc*beamGamma/energyPhoton1;  
  Double_t impactPar2_max=impactPar_min + 6.0*hbarc*beamGamma/energyPhoton2;
  Double_t impactPar1_delta=TMath::Exp(TMath::Log(impactPar1_max/impactPar_min)/Double_t(nSteps_GG));
  Double_t impactPar2_delta=TMath::Exp(TMath::Log(impactPar2_max/impactPar_min)/Double_t(nSteps_GG));
  
  const int ngi = 5;
  double weights[ngi] = {0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881};
  double abscissas[ngi] = {-0.1488743389816312, -0.4333953941292472, -0.6794095682990244, -0.8650633666889845, -0.9739065285171717};
  
  for(Double_t impactPar1_value = impactPar_min; impactPar1_value<=impactPar1_max; impactPar1_value *= impactPar1_delta){
    for(Int_t iFlux=0; iFlux<nFluxes; iFlux++)integral_Beam2[iFlux] = 0;
    
    for(Double_t impactPar2_value = impactPar_min; impactPar2_value<=impactPar2_max; impactPar2_value *= impactPar2_delta){
      for(Int_t iFlux=0; iFlux<nFluxes; iFlux++)integral_Phi[iFlux] = 0;
      
      for(Int_t k = 0; k < ngi; k++){
        Double_t impactPar_diff = TMath::Sqrt(impactPar1_value*impactPar1_value+impactPar2_value*impactPar2_value + 2.*impactPar1_value*impactPar2_value*TMath::Cos(pi*(abscissas[k]+1)));
	
	for(Int_t i = 0; i < maxNeutrons+1; i++) breakupProb[i] = gNucleusBreakupTable[i].Eval(impactPar_diff);
	
	diferential_Phi = weights[k]*gHadronicProbabilityTable->Eval(impactPar_diff);
	integral_Phi[0] += diferential_Phi;
	integral_Phi[nFluxes-1] += diferential_Phi*breakupProb[maxNeutrons]*breakupProb[maxNeutrons];
          for(Int_t i = 0; i < maxNeutrons; i++){
            for(Int_t j = i; j < maxNeutrons; j++){
              integral_Phi[FromMatrixToVector(i,j)+1] += diferential_Phi*breakupProb[i]*breakupProb[j];
	      }
	    }
        }
      diferential_Density = PhotonDensity(impactPar2_value,energyPhoton2)*2.0*pi*impactPar2_value*impactPar2_value*(1.0-1.0/impactPar2_delta);		      
      for(Int_t iFlux=0; iFlux<nFluxes; iFlux++) integral_Beam2[iFlux] += integral_Phi[iFlux]*diferential_Density;
      }
    diferential_Density = PhotonDensity(impactPar1_value,energyPhoton1)*2.0*pi*impactPar1_value*impactPar1_value*(1.0-1.0/impactPar1_delta);  
    for(Int_t iFlux=0; iFlux<nFluxes; iFlux++) flux_integral[iFlux] += integral_Beam2[iFlux]*diferential_Density;  
    }
    
  for(Int_t iFlux=1; iFlux<nFluxes; iFlux++) flux_integral[iFlux] /= flux_integral[0];
  return flux_integral; 
}
//______________________________________________________________________________
void NeutronGenerator::BuildPhotonFluxModulationTables()
{
  Double_t energy_min = 1e-5; 
  Double_t energy_max = 12.0 * beamGamma * hbarc/(2.0*nucleus_R);  
  Double_t energy_delta = TMath::Exp(TMath::Log(energy_max/energy_min)/Double_t(nSteps_Energy));
  energy_min *= energy_delta;
  Int_t iStep = 0;
  
  gPhotonFluxTable = new TGraph[nFluxes];
  gPhotonFluxTable[0].SetName("gPhotonFluxTable");
  for(Int_t i = 0; i < maxNeutrons; i++){
    for(Int_t j = i; j < maxNeutrons; j++){
      gPhotonFluxTable[FromMatrixToVector(i,j)+1].SetName(TString::Format("gPhotonFluxTable%d%d",i,j));
      }
    }
  cout<<"Building photon flux modulation tables"<<endl;  
  for(Double_t energy_value = energy_min; energy_value<=energy_max; energy_value *= energy_delta){
    Double_t* flux = PhotonFlux(energy_value,kFALSE);
    gPhotonFluxTable[0].SetPoint(iStep,energy_value,energy_value*flux[0]);
    for(Int_t i = 0; i < maxNeutrons; i++){
      for(Int_t j = i; j < maxNeutrons; j++){
        gPhotonFluxTable[FromMatrixToVector(i,j)+1].SetPoint(iStep,energy_value,flux[FromMatrixToVector(i,j)+1]);
	}
      }
    gPhotonFluxTable[nFluxes-1].SetPoint(iStep,energy_value,flux[nFluxes-1]);
    delete [] flux;
    iStep++;
    if(iStep%10 == 0){
      if(iStep > 10){ printf("\033[1A"); printf("\033[K");}   
       cout<<iStep<<"%"<<endl;
       }
   }
}
//______________________________________________________________________________
Double_t *NeutronGenerator::PhotonFlux(const Double_t energyPhoton, const Bool_t printout)
{
  Double_t flux_differential = 0.0;
  Double_t flux_differentialMod[nFluxes], flux_integralMod[nFluxes];
  Double_t *flux_integral = new Double_t[nFluxes];
  for(Int_t i=0; i<nFluxes; i++)flux_integral[i] = 0.0;
  Double_t dist;
  Double_t breakupProb[maxNeutrons+1];
 
  Double_t impactPar_min=1.6*nucleus_R;
  Double_t impactPar_max=impactPar_min + 6.0*hbarc*beamGamma/energyPhoton;  //6.0*adiabatic cutoff energy
  
  Double_t impactPar_delta=TMath::Exp(TMath::Log(impactPar_max/impactPar_min)/Double_t(nSteps_impactPar));
  Double_t R_delta = nucleus_R/Double_t(nSteps_R);
  Double_t Phi_delta=pi/Double_t(nSteps_Phi);
  
  impactPar_min*=impactPar_delta;
  Int_t nStep = 0;
  for(Double_t impactPar_value = impactPar_min; impactPar_value<=impactPar_max; impactPar_value *= impactPar_delta){
    //if(printout)cout<<impactPar_value<<", ";
    // When we get to b>20R_A change methods - just take the photon flux at the center of the nucleus
    if(impactPar_value < (10.0*nucleus_R)){
      //Integrate over nuclear surface. n.b. this assumes total shadowing - treat photons hitting the nucleus the same no matter where they strike
      flux_differential=0.0;
      for(Double_t R_value = 0.5*R_delta; R_value<=nucleus_R; R_value += R_delta){
	//use symmetry;  only integrate from 0 to pi (half circle)
	for(Double_t Phi_value = 0.5*Phi_delta; Phi_value<= pi; Phi_value += Phi_delta){
	  // dist is the distance from the center of the emitting nucleus to the point in question
	  dist=TMath::Sqrt((impactPar_value+R_value*TMath::Cos(Phi_value))*(impactPar_value+R_value*TMath::Cos(Phi_value))+(R_value*TMath::Sin(Phi_value))*(R_value*TMath::Sin(Phi_value)));
	  //The surface  element is 2.* delta phi* r * delta r. The '2' is because the phi integral only goes from 0 to pi
	  flux_differential += PhotonDensity(dist,energyPhoton)*2.0*Phi_delta*R_value*R_delta;
	  }
  	}
  	//average flux over the nuclear surface
  	flux_differential /= (pi*nucleus_R*nucleus_R);
      } 
    else{
      // We can do a single flux calculation, at the center of the nucleus Eq. 41 of Vidovic, Greiner and Soff, Phys.Rev.C47,2308(1993)
      flux_differential = PhotonDensity(impactPar_value,energyPhoton);
      }
    //multiply by volume element to get total flux in the volume element
    flux_differential *= 2.0*pi*impactPar_value*impactPar_value*(1.0-1.0/impactPar_delta);
    flux_differential *= gHadronicProbabilityTable->Eval(impactPar_value);
    
    //modulate by the probability of nuclear breakup as f(impactPar_value)
    for(Int_t i = 0; i < maxNeutrons+1; i++) breakupProb[i] = gNucleusBreakupTable[i].Eval(impactPar_value);
    
    breakupProb[2] = 0;
    for(Int_t i = 1; i < 3; i++) breakupProb[2] += gNucleusBreakupTable[i].Eval(impactPar_value); 
    
    breakupProb[3] = 0;
    for(Int_t i = 3; i < maxNeutrons; i++) breakupProb[3] += gNucleusBreakupTable[i].Eval(impactPar_value); 
    breakupProb[maxNeutrons - 1] = breakupProb[maxNeutrons];
      
    flux_differentialMod[nFluxes-1] = flux_differential*breakupProb[maxNeutrons]*breakupProb[maxNeutrons];
    flux_integral[nFluxes-1] += flux_differentialMod[nFluxes-1];
    for(Int_t i = 0; i < maxNeutrons; i++){
      for(Int_t j = i; j < maxNeutrons; j++){
        flux_differentialMod[FromMatrixToVector(i,j)] = flux_differential*breakupProb[i]*breakupProb[j];
	flux_integral[FromMatrixToVector(i,j)+1] += flux_differentialMod[FromMatrixToVector(i,j)];
	}
      }
    if(printout){
      fImpProfile[2].SetPoint(nStep,impactPar_value,(flux_differential-flux_differentialMod[FromMatrixToVector(0,0)]-flux_differentialMod[nFluxes-1])/(impactPar_value*(1.0-1.0/impactPar_delta)));
      fImpProfile[0].SetPoint(nStep,impactPar_value,flux_differentialMod[FromMatrixToVector(0,0)]/(impactPar_value*(1.0-1.0/impactPar_delta)));
      fImpProfile[1].SetPoint(nStep,impactPar_value,flux_differentialMod[nFluxes-1]/(impactPar_value*(1.0-1.0/impactPar_delta)));
      
      fImpProfile[3].SetPoint(nStep,impactPar_value,flux_differentialMod[FromMatrixToVector(0,1)]/(impactPar_value*(1.0-1.0/impactPar_delta)));
      fImpProfile[4].SetPoint(nStep,impactPar_value,flux_differentialMod[FromMatrixToVector(0,2)]/(impactPar_value*(1.0-1.0/impactPar_delta)));
      
      fImpProfile[5].SetPoint(nStep,impactPar_value,flux_differentialMod[FromMatrixToVector(maxNeutrons - 1,1)]/(impactPar_value*(1.0-1.0/impactPar_delta)));
      fImpProfile[6].SetPoint(nStep,impactPar_value,flux_differentialMod[FromMatrixToVector(maxNeutrons - 1,2)]/(impactPar_value*(1.0-1.0/impactPar_delta)));
      
      fImpProfile[7].SetPoint(nStep,impactPar_value,flux_differentialMod[FromMatrixToVector(0,3)]/(impactPar_value*(1.0-1.0/impactPar_delta)));
 
      }

 
    flux_integral[0] += flux_differential;
    nStep++;
    }
    
  if(printout){  
    //cout<<"0n12n = "<<energyPhoton*flux_integral[FromMatrixToVector(0,2)+1]<<endl;  
    //cout<<"0n3pn = "<<energyPhoton*flux_integral[FromMatrixToVector(0,3)+1]<<endl; 
    //cout<<"Xn1n = "<<energyPhoton*flux_integral[FromMatrixToVector(maxNeutrons - 1,1)+1]<<endl; 
    //cout<<"Xn2pn = "<<energyPhoton*flux_integral[FromMatrixToVector(maxNeutrons - 1,2)+1]<<endl;
    } 
  
  for(Int_t iFlux=1; iFlux<nFluxes; iFlux++) flux_integral[iFlux] /= flux_integral[0];
  
  if(printout){
    fOutputFile = new TFile("Output.root","RECREATE");
    for(Int_t i=0; i<=9; i++){
      fImpProfile[i].SetName(Form("fImpProfile%d",i));
      fImpProfile[i].Write(); 
      }
  }
  
  return flux_integral;
}
//______________________________________________________________________________
Double_t NeutronGenerator::PhotonDensity(const Double_t distance, const Double_t energyPhoton)
{
  Double_t Xvar=energyPhoton*distance/(hbarc*beamGamma); 
  return (nucleus_Z*nucleus_Z*alpha*energyPhoton)*(TMath::BesselK1(Xvar)*TMath::BesselK1(Xvar))/((pi*beamGamma*hbarc)*(pi*beamGamma*hbarc));
}
//______________________________________________________________________________
Double_t NeutronGenerator::PhotonDensity2(const Double_t distance, const Double_t energyPhoton)
{
  Double_t Xvar=energyPhoton*distance/(hbarcmev*gammaTarget); 
  return (nucleus_Z*nucleus_Z*alpha*energyPhoton)*(TMath::BesselK1(Xvar)*TMath::BesselK1(Xvar))/((pi*gammaTarget*hbarcmev)*(pi*gammaTarget*hbarcmev));
  
}
//______________________________________________________________________________
Double_t NeutronGenerator::GetEMDCrossSection()
{
  Double_t EMDxsection[maxNeutrons+1];
  for(Int_t i = 0; i < maxNeutrons+1; i++)EMDxsection[i] = 0;
  
  Double_t impactPar_min=1.8*nucleus_R;
  Double_t impactPar_max=1e9;
  Double_t nStepsEMD = 100000;
  
  Double_t impactPar_delta=TMath::Exp(TMath::Log(impactPar_max/impactPar_min)/Double_t(nStepsEMD));
  
  cout<<"Computing EMD cross section"<<endl; 
  Int_t iStep = 0;
  for(Double_t impactPar_value = impactPar_min; impactPar_value<=impactPar_max; impactPar_value *= impactPar_delta){
    Double_t* prob_Breakup = NucleusBreakupProbability(impactPar_value);
  
    EMDxsection[0] += gHadronicProbabilityTable->Eval(impactPar_value)*ComputeXnProbability(impactPar_value)*2*pi*impactPar_value*impactPar_value*(1.0-1.0/impactPar_delta);
    for(Int_t i = 1; i < 11; i++)EMDxsection[i] += gHadronicProbabilityTable->Eval(impactPar_value)*prob_Breakup[i]*2*pi*impactPar_value*impactPar_value*(1.0-1.0/impactPar_delta);
    delete [] prob_Breakup;
    //cout<<ComputeXnProbability(impactPar_value)<<endl;
    iStep++;
    if(iStep%10000 == 0){
      if(iStep > 10000){ printf("\033[1A"); printf("\033[K");}   
       cout<<iStep/1000<<"%"<<endl;
       }
    }
  
  for(Int_t i = 0; i < 11; i++)EMDxsection[i] /= 100; //barns 

  cout<<"EMD cross section total = "<<EMDxsection[0]<<" barn"<<endl;
  for(Int_t i = 1; i < 11; i++)cout<<"EMD cross section "<<i<<"n = "<<EMDxsection[i]<<" barn"<<endl;
  for(Int_t i = 1; i < 11; i++)cout<<"EMD cross section fraction "<<i<<"n = "<<EMDxsection[i]/EMDxsection[0]<<endl;
  return EMDxsection[maxNeutrons];
}
//______________________________________________________________________________
Double_t NeutronGenerator::ComputeXnProbability(const Double_t impactPar)
{
  //Double_t prob_Breakup = 0.0;
  Double_t prob_Xn=0.0;
  
  Double_t prob_Breakup = 0.0;

  //Maximum energy for GDR dissocation (in target frame, in MeV)
  Double_t photonEnergyLimit_Xn = 0.0;
  if (beamGamma > 500.){
      photonEnergyLimit_Xn=1.E10;
  }
  else{
      photonEnergyLimit_Xn=1.E7;
  }
  
   Double_t gk1=0,gk1m=0;
  //Compute the probabilities 
  //Xn dissociation
  Double_t maxPhotonEnergy = TMath::Min(photonEnergyLimit_Xn,4.0*gammaTarget*hbarcmev/impactPar);
  
  gk1m = PhotonDensity2(impactPar,gXsection[0].GetX()[0]);
  for(Int_t k = 1; gXsection[0].GetX()[k] < maxPhotonEnergy; k++){
    gk1 = PhotonDensity2(impactPar,gXsection[0].GetX()[k]);
    prob_Xn += (gXsection[0].GetX()[k]-gXsection[0].GetX()[k-1])*0.5*(gXsection[0].GetY()[k-1]*gk1m+gXsection[0].GetY()[k]*gk1);    
    gk1m = gk1;
    }
    
  prob_Breakup = 1-TMath::Exp(-1*prob_Xn);
  
  return prob_Breakup;
}
//______________________________________________________________________________
void NeutronGenerator::BuildNucleusBreakupProbabilityTable()
{
  gNucleusBreakupTable = new TGraph[maxNeutrons+1];
  for(Int_t i = 0; i < maxNeutrons; i++) gNucleusBreakupTable[i].SetName(TString::Format("gNucleusBreakupTable%d",i));
  gNucleusBreakupTable[maxNeutrons].SetName("gNucleusBreakupTableXn");

  cout<<"Building nucleus breakup probability tables"<<endl;  
  
  Double_t impactPar_min=1.6*nucleus_R;
  Double_t impactPar_max=impactPar_min + 6.0*hbarc*beamGamma; 
  Double_t impactPar_delta=TMath::Exp(TMath::Log(impactPar_max/impactPar_min)/Double_t(nSteps_impactPar*10));
  Int_t iStep = 0;
  
  for(Double_t impactPar_value = impactPar_min; impactPar_value<=impactPar_max; impactPar_value *= impactPar_delta){
    Double_t* prob_Breakup = NucleusBreakupProbability(impactPar_value);
    for(Int_t i = 0; i < maxNeutrons+1; i++)gNucleusBreakupTable[i].SetPoint(iStep,impactPar_value,prob_Breakup[i]);
    delete [] prob_Breakup;
    iStep++;
    if(iStep%120 == 0){
      if(iStep > 120){ printf("\033[1A"); printf("\033[K");}   
       cout<<iStep/12<<"%"<<endl;
       }
    }
}
//______________________________________________________________________________
Double_t *NeutronGenerator::NucleusBreakupProbability(const Double_t impactPar)
{
  //Double_t prob_Breakup = 0.0;
  Double_t prob_Xn=0.0;
  Double_t prob_Nn[maxNeutrons-1];
  for(Int_t i=0; i<maxNeutrons-1; i++)prob_Nn[i] = 0.0;
  
  Double_t *prob_Breakup = new Double_t[maxNeutrons+1];
  for(Int_t i=0; i<maxNeutrons+1; i++)prob_Breakup[i] = 0.0;

  //Maximum energy for GDR dissocation (in target frame, in MeV)
  Double_t photonEnergyLimit_Xn = 0.0;
  if (beamGamma > 500.){
      photonEnergyLimit_Xn=1.E10;
  }
  else{
      photonEnergyLimit_Xn=1.E7;
  }
  
   Double_t gk1=0,gk1m=0;
  //Compute the probabilities 
  //Xn dissociation
  Double_t maxPhotonEnergy = TMath::Min(photonEnergyLimit_Xn,4.0*gammaTarget*hbarcmev/impactPar);
  
  gk1m = PhotonDensity2(impactPar,gXsection[0].GetX()[0]);
  for(Int_t k = 1; gXsection[0].GetX()[k] < maxPhotonEnergy; k++){
    gk1 = PhotonDensity2(impactPar,gXsection[0].GetX()[k]);
    
    prob_Xn += (gXsection[0].GetX()[k]-gXsection[0].GetX()[k-1])*0.5*(gXsection[0].GetY()[k-1]*gk1m+gXsection[0].GetY()[k]*gk1);
    for(Int_t i=11; i<maxNeutrons; i++) 
      prob_Nn[i-1] += (gXsection[0].GetX()[k]-gXsection[0].GetX()[k-1])*0.5*(gXsection[0].GetY()[k-1]*GetBR(gXsection[0].GetX()[k-1],i)*gk1m+gXsection[0].GetY()[k]*GetBR(gXsection[0].GetX()[k],i)*gk1);
    
    gk1m = gk1;
    }
    
  prob_Breakup[maxNeutrons] = 1-TMath::Exp(-1*prob_Xn);
  prob_Breakup[0] = TMath::Exp(-1*prob_Xn);
  
  //Nn dissociation
  for(Int_t i=1; i<=10; i++){
    gk1m = PhotonDensity2(impactPar,gXsection[i].GetX()[0]);
    for(Int_t k = 1; gXsection[i].GetX()[k] < maxPhotonEnergy; k++){
      gk1 = PhotonDensity2(impactPar,gXsection[i].GetX()[k]);
      prob_Nn[i-1] += (gXsection[i].GetX()[k]-gXsection[i].GetX()[k-1])*0.5*(gXsection[i].GetY()[k-1]*gk1m+gXsection[i].GetY()[k]*gk1);
      gk1m = gk1;
      }
    }
  Double_t integralPoisson = 0;
  Int_t maxExcitation = 5;
  //if(impactPar<22) maxExcitation = 6;
  
  /* for(Int_t i=0; integralPoisson<0.999; i++){
    integralPoisson += TMath::PoissonI(i,prob_Xn);
    maxExcitation = i;
    } */
  
  Double_t norm_Xn = 0.0;  
  for(Int_t i=1; i<maxNeutrons; i++){ 
    prob_Breakup[i] += prob_Nn[i-1];
    if(maxExcitation<2) break;
    for(Int_t j=1; i+j<maxNeutrons; j++){
      prob_Breakup[i+j] += prob_Nn[i-1]*prob_Nn[j-1]/2;
      if(maxExcitation<3) break;
      for(Int_t k=1; i+j+k<maxNeutrons; k++){ 
        prob_Breakup[i+j+k] += prob_Nn[i-1]*prob_Nn[j-1]*prob_Nn[k-1]/6;
	if(maxExcitation<4) break;
        for(Int_t l=1; i+j+k+l<maxNeutrons; l++){
          prob_Breakup[i+j+k+l] += prob_Nn[i-1]*prob_Nn[j-1]*prob_Nn[k-1]*prob_Nn[l-1]/24;
	  if(maxExcitation<5) break;
	  for(Int_t m=1; i+j+k+l+m<maxNeutrons; m++){
            prob_Breakup[i+j+k+l+m] += prob_Nn[i-1]*prob_Nn[j-1]*prob_Nn[k-1]*prob_Nn[l-1]*prob_Nn[m-1]/120;
	    if(maxExcitation<6) break;
	    for(Int_t n=1; i+j+k+l+m+n<maxNeutrons; n++){
              prob_Breakup[i+j+k+l+m+n] += prob_Nn[i-1]*prob_Nn[j-1]*prob_Nn[k-1]*prob_Nn[l-1]*prob_Nn[m-1]*prob_Nn[n-1]/720;
	      if(maxExcitation<7) break;
	        for(Int_t p=1; i+j+k+l+m+n+p<maxNeutrons; p++){
                  prob_Breakup[i+j+k+l+m+n+p] += prob_Nn[i-1]*prob_Nn[j-1]*prob_Nn[k-1]*prob_Nn[l-1]*prob_Nn[m-1]*prob_Nn[n-1]*prob_Nn[p-1]/5040;
	        }
	      }
	    }
	  }
        }
      }
      prob_Breakup[i] *= TMath::Exp(-1*prob_Xn);
      norm_Xn += prob_Breakup[i];
    }
    
  for(Int_t i=1; i<maxNeutrons; i++) if(norm_Xn != 0)prob_Breakup[i] *= prob_Breakup[maxNeutrons]/norm_Xn;
  if(prob_Breakup[maxNeutrons] > 0)gUnitaryLeak->SetPoint(gUnitaryLeak->GetN(),impactPar,norm_Xn/prob_Breakup[maxNeutrons]);
  
  return prob_Breakup;  
}

//______________________________________________________________________________
Double_t NeutronGenerator::GetBR(const Double_t energyPhoton, const Int_t nNeutrons)
{
  return hBranchingRatioMap->GetBinContent(hBranchingRatioMap->FindBin(energyPhoton,nNeutrons));
} 
//______________________________________________________________________________
void NeutronGenerator::Initialize()
{
  gXsection = new TGraph[11];
  TGraph gXsectionF;
  
  if(fNucleus == kPb208){
    InsertDataset(gXsection[0],ReadXSFile("XS_Pb_208_Xn_Veyssiere.txt"));
    
    for(Int_t nNeutrons = 1; nNeutrons<=2; nNeutrons++)
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Pb_208_%dn_Veyssiere.txt",nNeutrons))); 
    } 
  if(fNucleus == kAu197){
    InsertDataset(gXsection[0],ReadXSFile("XS_Au_197_Xn_Veyssiere.txt"));
    
    for(Int_t nNeutrons = 1; nNeutrons<=2; nNeutrons++)
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Au_197_%dn_Veyssiere.txt",nNeutrons)));
    }
  if(fNucleus == kPb208 || fNucleus == kAu197){
    InsertDataset(gXsection[0],ReadXSFile("XS_Au_197_Xn_Michalowski.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_Pb_208_Xn_Caldwell.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_Pb_208_Xn_Lepretre.txt"));
    
    /*/ 
  Double_t *xsData = ReadXSFile(".XS_n_2_Armstrong.txt");
  Double_t *pData = ReadXSFile("XS_p_1_Armstrong.txt");
  for(Int_t i = 1; i<xsData[0]; i++)xsData[2*i] += pData[2*i];
  InsertDataset(gXsection[0],xsData);
  /*/
    
    InsertDataset(gXsection[0],ReadXSFile("XS_Pb_208_Xn_Bianchi.txt"));
    
    for(Int_t nNeutrons = 2; nNeutrons<=10; nNeutrons++)
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Pb_208_%dn_Lepretre.txt",nNeutrons)));
    }
    
  if(fNucleus == kU238){
    //InsertDataset(gXsection[0],ReadXSFile("XS_U_238_Xn_Caldwell.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_U_238_Xn_Michalowski.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_U_238_Xn_Veyssiere.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_U_238_Xn_Lepretre.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_U_238_Xn_Bergere.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_U_238_Xn_Bianchi.txt"));
    
    InsertDataset(gXsectionF,ReadXSFile("XS_U_238_Xn_Michalowski.txt"));
    InsertDataset(gXsectionF,ReadXSFile("XS_U_238_F_Veyssiere.txt"));
    InsertDataset(gXsectionF,ReadXSFile("XS_U_238_Xn_Lepretre.txt"));
    InsertDataset(gXsectionF,ReadXSFile("XS_U_238_Xn_Bergere.txt"));
    InsertDataset(gXsectionF,ReadXSFile("XS_U_238_Xn_Bianchi.txt"));
    InsertDataset(gXsectionF,ReadXSFile("XS_Pb_208_Xn_Carlos.txt"));
    InsertDataset(gXsectionF,MakeReggeParametrization(gXsectionF.GetX()[gXsectionF.GetN()-1], 150));
    
    for(Int_t nNeutrons = 1; nNeutrons<=2; nNeutrons++){
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_U_238_%dn_Veyssiere.txt",nNeutrons)));
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_U_238_%dn_Caldwell.txt",nNeutrons)));
      }
    for(Int_t nNeutrons = 2; nNeutrons<=10; nNeutrons++)
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Pb_208_%dn_Lepretre.txt",nNeutrons)));
    }
  
  if(fNucleus == kXe129){
    InsertDataset(gXsection[0],ReadXSFile("XS_I_127_Xn_Bergere.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_I_127_Xn_Bramblett.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_Ce_140_Xn_Lepretre.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_Sn_118_Xn_Lepretre.txt"));
    InsertDataset(gXsection[0],ReadXSFile("XS_Sn_118_Xn_Bianchi.txt"));
    
    for(Int_t nNeutrons = 1; nNeutrons<=2; nNeutrons++){
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_I_127_%dn_Bergere.txt",nNeutrons)));
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_I_127_%dn_Bramblett.txt",nNeutrons)));
      //InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Pb_208_%dn_Veyssiere.txt",nNeutrons)));
      }
    for(Int_t nNeutrons = 2; nNeutrons<=7; nNeutrons++){
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Ce_140_%dn_Lepretre.txt",nNeutrons)));
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Sn_118_%dn_Lepretre.txt",nNeutrons)));
      //InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Pb_208_%dn_Lepretre.txt",nNeutrons)));
      }
    InsertDataset(gXsection[8],ReadXSFile("XS_Ce_140_8n_Lepretre.txt"));
    for(Int_t nNeutrons = 9; nNeutrons<=10; nNeutrons++)
      InsertDataset(gXsection[nNeutrons],ReadXSFile(Form("XS_Pb_208_%dn_Lepretre.txt",nNeutrons)));   
    }

  InsertDataset(gXsection[0],ReadXSFile("XS_Pb_208_Xn_Carlos.txt"));

  InsertDataset(gXsection[0],MakeReggeParametrization(gXsection[0].GetX()[gXsection[0].GetN()-1], 150));
   
  //"Intepreted" ENDF files from https://www-nds.iaea.org/exfor/endf.htm
  if (fNucleus == kAu197) ReadENDF("ENDF_Au197.txt");
  if (fNucleus == kPb208) ReadENDF("ENDF_Pb208.txt");
  if (fNucleus == kU238) ReadENDF("ENDF_U238.txt");
  if (fNucleus == kXe129) ReadENDF("ENDF_Cs133.txt");
 
  fitMeanNeutronN = new TF1("fitMeanNeutronN", "[0]*TMath::Log(x)+[1]", neutronSepThr, 1e9);
  fitMeanNeutronN->SetParameters(2.16117,-4.47756);
 
  fitWidthNeutronN = new TF1("fitWidthNeutronN", "[0]*TMath::Log(x)+[1]", neutronSepThr, 1e9);
  fitWidthNeutronN->SetParameters(1.53,-4.68);
 
  TF1 *fitMeanFissionN = new TF1("fitMeanFissionN", "[0]+([1]-[0])*(1-TMath::Exp([2]*(x-[4])))",neutronSepThr,1e9);
  fitMeanFissionN->SetParameters(2.41,12.65,-0.015,4.8);

  TGraph gBranchingRatio;
  if(fNucleus == kU238)gBranchingRatio = gXsectionF;
  else gBranchingRatio = gXsection[0];

  Double_t mean; 
  Double_t width = 0;
  Double_t integral;
  Double_t energyValue;
  Double_t probValue;  
  hBranchingRatioMap = new TH2D("hBranchingRatioMap","hBranchingRatioMap",gBranchingRatio.GetN()-1,gBranchingRatio.GetX(),50,0,50);
  for(Int_t i = 0; i<hBranchingRatioMap->GetNbinsX(); i++){
    energyValue = hBranchingRatioMap->GetXaxis()->GetBinCenter(i+1);
    if(energyValue > saturationEnergy) energyValue = saturationEnergy;
    if(fNucleus == kU238)
      mean = fitMeanFissionN->Eval(energyValue);
    else{
      mean = fitMeanNeutronN->Eval(energyValue);
      width = fitWidthNeutronN->Eval(energyValue);
      }
    integral = 0;
    for(Int_t j = 1; j<50; j++){
      if(fNucleus == kU238)
        probValue = TMath::Poisson(j,mean);
      else{
        if(j == 1 && energyValue < 2*neutronSepThr) probValue = 1.0;
        else if(energyValue < j*neutronSepThr) probValue = 0.0;
        else probValue = TMath::Gaus(j,mean,width,kTRUE);
        }
      
      hBranchingRatioMap->SetBinContent(i+1,j+1,probValue);
      integral+=probValue;
      }
    for(Int_t j = 1; j<50; j++)if(integral != 0)hBranchingRatioMap->SetBinContent(i+1,j+1,hBranchingRatioMap->GetBinContent(i+1,j+1)/integral);
    }
  
  for(Int_t nNeutrons=1; nNeutrons<=10; nNeutrons++){  
    Int_t startingPoint = 0;
    if(gXsection[nNeutrons].GetN() == 0)startingPoint = 0;
    else{ 
      for(Int_t i = 0; gXsection[0].GetX()[i] < gXsection[nNeutrons].GetX()[gXsection[nNeutrons].GetN()-1]; i++) startingPoint = i+1;
      if(fNucleus == kU238)
        for(Int_t iPoint = 0; iPoint<gXsection[nNeutrons].GetN(); iPoint++)
	  gXsection[nNeutrons].SetPoint(iPoint,gXsection[nNeutrons].GetX()[iPoint],gXsection[nNeutrons].GetY()[iPoint]+gXsection[0].Eval(gXsection[nNeutrons].GetX()[iPoint]));
      }
    for(Int_t iPoint = startingPoint; iPoint<gXsection[0].GetN(); iPoint++)
      gXsection[nNeutrons].SetPoint(gXsection[nNeutrons].GetN(), gXsection[0].GetX()[iPoint], gXsection[0].GetY()[iPoint]*GetBR(gXsection[0].GetX()[iPoint],nNeutrons));
      
    }
  //hBranchingRatioMap->Draw("COLZ");
  
  delete fitMeanFissionN;
  
  hXsection = new TH1D[11]; 
  for(Int_t i = 0; i<=10; i++){
    hXsection[i] = TH1D(TString::Format("hXsection%d",i)," ",gXsection[i].GetN()-1,gXsection[i].GetX());
    for(Int_t iBin = 1; iBin <= hXsection[i].GetNbinsX(); iBin++)hXsection[i].SetBinContent(iBin,gXsection[i].Eval(hXsection[i].GetBinCenter(iBin)));
    hXsection[i].SetDirectory(0);
    }

  if(fRunMode != kFlatMultiplicity && fRunMode != k1n1n){
    BuildHadronicInteractionProbabilityTable();
    BuildNucleusBreakupProbabilityTable();
    //if(fProductionMode == kPhotonPomeron)
    BuildPhotonFluxModulationTables();
    if(fProductionMode == kTwoPhoton)
    BuildTwoPhotonFluxModulationTables();
    }
  InitQAhistograms();
  
}
//______________________________________________________________________________
void NeutronGenerator::BuildHadronicInteractionProbabilityTable()
{
  Int_t iStep = 0;
  gHadronicProbabilityTable = new TGraph();
  gHadronicProbabilityTable->SetName("gHadronicProbabilityTable");
  
  Double_t impactPar_delta=1.007;
  Double_t probValue = 1.0;
  TString hadrModel = fHadronicIntModel == kHardSphere ? "Hard sphere":"Glauber";
  cout<<"Building hadronic interaction probability table for "<<hadrModel.Data()<<" model"<<endl;

  for(Double_t impactPar_value = nucleus_R; impactPar_value<=25.0; impactPar_value *= impactPar_delta){
    probValue = TMath::Exp(-1.0*HadronicInteractionProbability(impactPar_value));
    gHadronicProbabilityTable->SetPoint(iStep++,impactPar_value,probValue);
    
    if(Int_t(iStep/1.91)%10 == 0){
      if(iStep != 1){ printf("\033[1A"); printf("\033[K");}   
       cout<<Int_t(iStep/1.91)<<"%"<<endl;
       }
    }
    //gHadronicProbabilityTable->Draw();
}
//______________________________________________________________________________
Double_t NeutronGenerator::HadronicInteractionProbability(const Double_t impactPar)
{
  Double_t prob_Handronic = 0.0;
  if(fHadronicIntModel == kHardSphere)prob_Handronic = impactPar<2*nucleus_R ? 1e20 : 0;

  if(fHadronicIntModel == kGlauber){
    Double_t nucleus_WSSD = 0.535; //Woods-saxon skin depth
    Double_t energyCMS = 2*beamGamma*0.938; 
    // This equation is from section 50 of the particle data book, the subsection on "Total Hadronic Cross-Sections, using the parameterization for sqrt{s} > 7 GeV.
    // only the first and second terms contribute significantly, but leave them all here for good measure
    Double_t xSection_IntNN = 0.1*(0.2838*TMath::Power(TMath::Log(energyCMS),2)+33.73+13.67*TMath::Power(energyCMS,-0.412)-7.77*TMath::Power(energyCMS,-0.5626));

    Double_t maxDistance = nucleus_R+5.0;	
    Double_t RZ_delta = 0.01;
    Double_t nuclearThickness_value = 0, nuclearThickness_norm = 0; 

    if(!gNuclearThickness){
      Int_t nPoint = 0;
      gNuclearThickness = new TGraph;
      gNuclearThickness->SetName("gNuclearThickness");
      // This calculates T_A(b) and stores it in TGraph 
      for (Double_t R_value = 0; R_value <= maxDistance; R_value+=RZ_delta){
        nuclearThickness_value = 0.0;
        for (Double_t Z_value = RZ_delta/2; Z_value <= maxDistance; Z_value+=RZ_delta) nuclearThickness_value += 1.0/(1.0+TMath::Exp((TMath::Sqrt(R_value*R_value+Z_value*Z_value)-nucleus_R)/nucleus_WSSD));
        nuclearThickness_value *= 2.0*RZ_delta;
	//cout<<"R = "<<R_value<<" thickness = "<<nuclearThickness_value<<endl;
	
        nuclearThickness_norm += (R_value+RZ_delta)*nuclearThickness_value*RZ_delta*2.0*pi;
        nuclearThickness_value *= nucleus_A;
        gNuclearThickness->SetPoint(nPoint++,R_value,nuclearThickness_value);
	
        }
       for (Int_t i=0; i<gNuclearThickness->GetN(); i++) gNuclearThickness->GetY()[i] *= 1.0/nuclearThickness_norm; 
      //gNuclearThickness->Draw();
      }
    prob_Handronic = 0.0;
    Double_t XY_delta = 0.05;

    if(impactPar > 25.0) return prob_Handronic;
    for(Double_t Y_value = XY_delta/2; Y_value <= maxDistance; Y_value += XY_delta){
      for(Double_t X_value = -maxDistance+XY_delta/2; X_value <= maxDistance; X_value += XY_delta){
        Double_t dist_bR = TMath::Sqrt((impactPar-X_value)*(impactPar-X_value)+Y_value*Y_value);
        Double_t dist_R = TMath::Sqrt(X_value*X_value+Y_value*Y_value);
        if(dist_bR >= maxDistance || dist_R >= maxDistance) continue;

        prob_Handronic += 2.0*gNuclearThickness->Eval(dist_bR)*(1.0-TMath::Exp(-xSection_IntNN*gNuclearThickness->Eval(dist_R)))*XY_delta*XY_delta;
        }
      }
    }
  return prob_Handronic;
}
//______________________________________________________________________________
Double_t *NeutronGenerator::MakeReggeParametrization(Double_t startEnergy, Int_t nPoints){

  Double_t *xsPoints = new Double_t[2*nPoints+1];
  xsPoints[0] = nPoints;

  //Regge parameters
  Double_t x=0,y=0,eps=0,eta=0,em=0,exx=0,s=0,ictr=0,pom=0,vec=0;
  x = 0.0677;
  y = 0.129;
  eps = 0.0808;
  eta = 0.4525;
  em = 0.94;
  exx = TMath::Power(10,0.05);

  //Regge model for high energy
  s = 0.002*em*startEnergy;
  for ( Int_t j = 1; j <= nPoints; j++ ) {
    s *= exx;
    pom = x*TMath::Power(s,eps);
    vec = y*TMath::Power(s,(-eta));
  
    xsPoints[2*j-1] = 1000.0*0.5*(s-em*em)/em;
    xsPoints[2*j] = 0.1*0.65*nucleus_A*(pom+vec);
    }
  return xsPoints;
}
//______________________________________________________________________________
void NeutronGenerator::InsertDataset(TGraph& graph, Double_t *dataset){

  for(Int_t i = 1; i<=dataset[0]; i++){
    InsertPoint(graph,dataset[2*i-1],dataset[2*i]);
  }
  delete [] dataset;
}
//______________________________________________________________________________
void NeutronGenerator::InsertPoint(TGraph& graph, Double_t x, Double_t y){

  Bool_t inserted = kFALSE;
  if(graph.GetN() != 0){
    if(graph.GetX()[0] > x){
      for(Int_t j = graph.GetN(); j>0; j--)graph.SetPoint(j,graph.GetX()[j-1],graph.GetY()[j-1]);
        graph.SetPoint(0,x,y);
        inserted = kTRUE;
      }
    else{
      for(Int_t i = 1; i<graph.GetN(); i++){
        if(graph.GetX()[i-1] < x && graph.GetX()[i] >= x){
          for(Int_t j = graph.GetN(); j>i; j--)graph.SetPoint(j,graph.GetX()[j-1],graph.GetY()[j-1]);
	  if(graph.GetX()[i] == x) x -= 1e-6;
          graph.SetPoint(i,x,y);
          inserted = kTRUE;
          break;
          }
        }
      }
    }
  if(!inserted)graph.SetPoint(graph.GetN(),x,y); 
}
//______________________________________________________________________________
Double_t *NeutronGenerator::ReadXSFile(const char *filename, Double_t scale){

  TString fullFileName = fDataPath;
  fullFileName += filename;

  std::string newLine;
  TString token;
  Ssiz_t from = 0;
  Ssiz_t pos = 0;
  UInt_t nPoints = 0; 
  
  ifstream infile(gSystem->ExpandPathName(fullFileName.Data()));
  
  while(getline(infile,newLine)){
    TString lineString = newLine;
    if(lineString.Contains("#") || lineString.IsWhitespace())continue;
    nPoints++;
    }
  infile.clear();
  infile.seekg(0, ios::beg);
  
  Double_t *xsPoints = new Double_t[2*nPoints+1];
  xsPoints[0] = nPoints;

  nPoints = 1;
  while(getline(infile,newLine)){
    TString lineString = newLine;
    if(lineString.Contains("#") || lineString.IsWhitespace())continue;
    from = 0;
    while (lineString.Tokenize(token, from, "[|]")){
      do{ 
        pos = token.Index("\t");
	if(pos >= 1)token.Replace(pos, 1, " ", 1);
	}
	while(pos >= 1);
        xsPoints[nPoints++] = token.Atof();
      }
  }
  Int_t fileNucleusA = 1;
  TString fileString = filename;
  //cout<<fileString.Data()<<endl;
  //for(Int_t i = 1; i<=xsPoints[0]; i++)cout<<xsPoints[2*i-1]<<" "<<xsPoints[2*i]<<endl;
  from = 0;
  while (fileString.Tokenize(token, from, "[_]")){
    if(token.IsFloat()){
      fileNucleusA = token.Atoi();
      break;
      }
    }
  //Scaling by A, conversion to fm
  for(Int_t i = 1; i<=xsPoints[0]; i++){
    xsPoints[2*i] *= scale;
    if(xsPoints[2*i] < 0)xsPoints[2*i-1] = 0;
    }
  if(fileNucleusA == 1) for(Int_t i = 1; i<=xsPoints[0]; i++)xsPoints[2*i] *= 0.1*nucleus_Z;
  else if(fileNucleusA == 2) for(Int_t i = 1; i<=xsPoints[0]; i++)xsPoints[2*i] *= 0.1*(nucleus_A-nucleus_Z);
  
  else for(Int_t i = 1; i<=xsPoints[0]; i++)xsPoints[2*i] *= 0.1*((Float_t)nucleus_A/(Float_t)fileNucleusA);
  
  return xsPoints;
}
//______________________________________________________________________________
void NeutronGenerator::InitQAhistograms(){
  
  fQAhistList = new TList();
  fQAhistList ->SetOwner();
    
  hNeutronMultiplicity = CreateHist2D("hNeutronMultiplicity","Neutron multiplicity",100,0,100,100,0,100,"Number of neutrons (beam 1)","Number of neutrons (beam 2)","Counts");
  fQAhistList->Add(hNeutronMultiplicity);
  
  hEnergyGen = CreateHist1D("hEnergyGen","Generated energy before boost",3000,900,1200,"Energy [MeV]","Counts");
  fQAhistList->Add(hEnergyGen);
  hKinEnergyGen = CreateHist1D("hKinEnergyGen","Generated kinetic energy before boost",3000,0,300,"Energy [MeV]","Counts");
  fQAhistList->Add(hKinEnergyGen);
  hMomGen = CreateHist1D("hMomGen","Generated total momentum before boost",1200,0,300,"Momentum [MeV/c]","Counts");
  fQAhistList->Add(hMomGen);
  hPhiGen = CreateHist1D("hPhiGen","Generated #phi before boost",100,0,pi,"#phi","Counts");
  fQAhistList->Add(hPhiGen);
  hThetaGen = CreateHist1D("hThetaGen","Generated #theta before boost",100,0,pi,"#theta","Counts");
  fQAhistList->Add(hThetaGen);
  
  hEnergyBoosted = CreateHist1D("hEnergyBoosted","Generated energy after boost",3000,1e6,5e6,"Energy [MeV]","Counts");
  fQAhistList->Add(hEnergyBoosted);
  hMomBoosted = CreateHist1D("hMomBoosted","Generated total momentum after boost",1200,1e6,5e6,"Momentum [MeV/c]","Counts");
  fQAhistList->Add(hMomBoosted);
  hPhiBoosted = CreateHist1D("hPhiBoosted","Generated #phi after boost",100,0,pi,"#phi","Counts");
  fQAhistList->Add(hPhiBoosted);
  hThetaBoosted = CreateHist1D("hThetaBoosted","Generated #theta after boost",1000,0,pi,"#theta","Counts");
  fQAhistList->Add(hThetaBoosted);
  
  hNeutronRapidity = CreateHist1D("hNeutronRapidity","Neutron rapidity",2000,-20,20,"y","Counts");
  fQAhistList->Add(hNeutronRapidity);
  hNeutronEta = CreateHist1D("hNeutronEta","Neutron #eta",2000,-20,20,"#eta","Counts");
  fQAhistList->Add(hNeutronEta);
  
  hEnergyBin = CreateHist1D("hEnergyBin","Bin of the ENDF 2D hist used for energy generation",140,0,140,"Bin","Counts");
  fQAhistList->Add(hEnergyBin);
  hEnergyForNeutronMulti = CreateHist1D("hEnergyForNeutronMulti","hEnergyForNeutronMulti",gXsection[0].GetN()-1,gXsection[0].GetX(),"Energy [MeV]","Counts");
  fQAhistList->Add(hEnergyForNeutronMulti);
  
  hRapidityVM = CreateHist1D("hRapidityVM","Rapidity of VM used in event generation",1000,-10,10,"Rapidity","Counts");
  fQAhistList->Add(hRapidityVM);
  hMassVM = CreateHist1D("hMassVM","Mass of VM used in event generation",1000,0,10,"Rapidity","Counts");
  fQAhistList->Add(hMassVM);
  
  if(fRunMode != kFlatMultiplicity && fRunMode != k1n1n){
    if(fProductionMode == kPhotonPomeron) 
      hPhotonK = CreateHist1D("hPhotonK","Virtual photon k=0.5*M_{VM}*exp(y_{VM}) used in event generation",nSteps_Energy-1,gPhotonFluxTable[0].GetX(),"k [GeV/c]","Counts");
    if(fProductionMode == kTwoPhoton){
      Double_t energy_min = 0.5*fMassMin*TMath::Exp(-TMath::Max(TMath::Abs(fRapMin), TMath::Abs(fRapMax))); 
      Double_t energy_max = 0.5*fMassMax*TMath::Exp(TMath::Max(TMath::Abs(fRapMin), TMath::Abs(fRapMax))); 
      hPhotonK = CreateHist1D("hPhotonK","Virtual photon k=0.5*M_{VM}*exp(y_{VM}) used in event generation",10*nSteps_GG,energy_min,energy_max,"k [GeV/c]","Counts");
      }
    fQAhistList->Add(hPhotonK);
    }
    
  hProbabilityXn = CreateHist1D("hProbabilityXn","Single side probability of Xn",1000,0,1,"Probability","Counts");
  fQAhistList->Add(hProbabilityXn);
  
  fImpProfile = new TGraph[10];

}
//______________________________________________________________________________
void NeutronGenerator::SetBeamParameters(Int_t nuclZ, Int_t nuclA, Double_t gamma)
{
  if(nuclZ == 82 && nuclA == 208)SetBeamParameters(kPb208, gamma);
  if(nuclZ == 79 && nuclA == 197)SetBeamParameters(kAu197, gamma);
  if(nuclZ == 92 && nuclA == 238)SetBeamParameters(kU238, gamma);
  if(nuclZ == 54 && nuclA == 129)SetBeamParameters(kXe129, gamma);
}
//______________________________________________________________________________
void NeutronGenerator::SetBeamParameters(Nucleus_t nucleus, Double_t gamma)
{
  beamGamma = gamma;
  gammaTarget = 2.0*beamGamma*beamGamma-1.0;
  fNucleus = nucleus;
  
  cout<<"nOOn setup for ";
  
  if(nucleus == kPb208){
    nucleus_Z = 82; 
    nucleus_A = 208;
    neutronSepThr = 7.4;
    nucleus_R = 6.624;
    cout<<"Pb208 beam with gamma = ";
    }
  if(nucleus == kAu197){
    nucleus_Z = 79; 
    nucleus_A = 197;
    neutronSepThr = 8.1;
    nucleus_R = 6.38;
    cout<<"Au197 beam with gamma = ";
    } 
  if(nucleus == kU238){
    nucleus_Z = 92; 
    nucleus_A = 238;
    neutronSepThr = 6.152;
    nucleus_R = 6.624;    //FIXME: This is from Pb208 
    cout<<"U238 beam with gamma = ";
    }
  if(nucleus == kXe129){
    nucleus_Z = 54; 
    nucleus_A = 129;
    neutronSepThr = 9.05;
    nucleus_R = 5.36;
    cout<<"Xe129 beam with gamma = ";
    }
  cout<<TString::Format("%.1f",beamGamma)<<endl;  
}
//______________________________________________________________________________
void NeutronGenerator::LoadENDF(const char *filename,const char *histname)
{
  fENDFFile = new TFile(filename,"READ");
  hENDF_2D = (TH2D*) fENDFFile->Get(histname);
  hENDF_2D->SetDirectory(0);
  fENDFFile->Close();
}
//______________________________________________________________________________
void NeutronGenerator::ReadENDF(const char *filename, Bool_t saveToFile)
{
  std::string newLine;
  TString token;
  Ssiz_t from = 0;
  Ssiz_t pos = 0;
  
  TString energy;
  Double_t energyPhoton = 0;
  Double_t energyNeutron = 0, probability = 0; 
  
  TString fullFileName = fDataPath;
  fullFileName += filename;
  
  ifstream infile(gSystem->ExpandPathName(fullFileName.Data()));
  
  Int_t counterPhE = 0;
  while(getline(infile,newLine)){
    TString lineString = newLine;
    if(lineString.Contains("Energy"))counterPhE++;
    }
  infile.clear();
  infile.seekg(0, ios::beg);
  
  const Int_t nBinsPhE = counterPhE;
  Double_t binsPhE[nBinsPhE+1];
  counterPhE = 1;
  binsPhE[0] = 0.0;
  while(getline(infile,newLine)){
    TString lineString = newLine;
    if(lineString.Contains("Energy")){
      from = 0;
      while (lineString.Tokenize(token, from, "[ ]"))if(!token.IsWhitespace()){
        pos = token.Index("-");
        if (pos >= 1) token.Replace(pos, 1, "E-", 2);
	pos = token.Index("+");
        if (pos >= 1) token.Replace(pos, 1, "E+", 2);
	if(token.IsFloat())energyPhoton = token.Atof();
        }
      binsPhE[counterPhE++] = energyPhoton/1e6;
      }
    }
  infile.clear();
  infile.seekg(0, ios::beg);
  
  hENDF_2D = new TH2D("hENDF_2D","hENDF_2D",nBinsPhE,binsPhE,1400,0,140);
  hENDF_2D->GetXaxis()->SetTitle("Photon energy [MeV]");
  hENDF_2D->GetYaxis()->SetTitle("Neutron energy [MeV]"); 
  hENDF_2D->SetStats(kFALSE);
  
  Int_t counterNuE = 0;
  TGraph *graphNuE = NULL;
  counterPhE = 0;

  while(getline(infile,newLine)){
    TString lineString = newLine;
    if(lineString.Contains("eV") || lineString.Contains("----------") || lineString.IsWhitespace())continue;
    //cout<<lineString.Data()<<endl;
    if(lineString.Contains("Energy")){
      if(graphNuE){
        Double_t maxX,maxY;
	graphNuE->GetPoint(graphNuE->GetN()-1,maxX,maxY);
        for(Int_t j = 1; j<=maxX*10; j++){
          Double_t value = graphNuE->Eval(hENDF_2D->GetYaxis()->GetBinCenter(j));
	  if(value<1e-11)value = 0;
          hENDF_2D->SetBinContent(counterPhE,j,value);
	  }
	}
      counterPhE++;
      delete graphNuE;
      counterNuE = 0;
      graphNuE = new TGraph();
      }
    else{
      UInt_t iNumber = 0;
      from = 0;
      while (lineString.Tokenize(token, from, "[ ]"))if(!token.IsWhitespace()){
        //energy prob r energy prob r
	iNumber++;
        pos = token.Index("-");
        if (pos >= 1) token.Replace(pos, 1, "E-", 2);
	pos = token.Index("+");
        if (pos >= 1) token.Replace(pos, 1, "E+", 2);
	if(token.IsFloat() && (iNumber == 1 || iNumber == 4))energyNeutron = token.Atof();
	if(token.IsFloat() && (iNumber == 2 || iNumber == 5))probability = token.Atof();
	if((iNumber == 2 || iNumber == 5) && probability > 1e-11){ 
	  if(probability>1e-11){
	    graphNuE->SetPoint(counterNuE,energyNeutron/1e6,probability);
	    counterNuE++;
	    }
	  }
        }
      }  
  }
  if(graphNuE){
    Double_t maxX,maxY;
    graphNuE->GetPoint(graphNuE->GetN()-1,maxX,maxY);
    for(Int_t j = 1; j<=maxX*10; j++){
      Double_t value = graphNuE->Eval(hENDF_2D->GetYaxis()->GetBinCenter(j));
      if(value<1e-11)value = 0;
      hENDF_2D->SetBinContent(counterPhE,j,value);
      }
    }
  delete graphNuE;
  //hENDF_2D->Draw();
  if(saveToFile){
    fENDFFile = new TFile("hENDF.root","RECREATE");
    hENDF_2D->Write();
    fENDFFile->Close();
    }
}

/*/
if(fRunMode == kStarlightAscii){
      TLorentzVector totgen, genpart;
      totgen.SetXYZM(0.0,0.0,0.0,0.0);
      std::string newLine;
      do{ 
        getline(fInputStarlightAscii,newLine);
	lineString = newLine;
	if(!lineString.Contains("EVENT"))fOutputStarlightAscii<<newLine<<endl;
	if(lineString.Contains("TRACK")){ 
          TObjArray *splitLine = lineString.Tokenize(" ");
	  genpart.SetXYZM(((TObjString*)splitLine->At(2))->String().Atof(), ((TObjString*)splitLine->At(3))->String().Atof(),
	  ((TObjString*)splitLine->At(4))->String().Atof(), pdgData.GetParticle(((TObjString*)splitLine->At(8))->String().Atoi())->Mass());
	  totgen += genpart;
	  } 
	}while(!lineString.Contains("EVENT"));
      VMmass = totgen.M();
      VMrapidity = totgen.Rapidity();
      }

/*/
//________________________________________________________________________
TH1D* NeutronGenerator::CreateHist1D(const char* name, const char* title,Int_t nBins, Double_t xMin, Double_t xMax, const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH1D* result = new TH1D(name, title, nBins, xMin, xMax);
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  return result;
}
//________________________________________________________________________
TH1D* NeutronGenerator::CreateHist1D(const char* name, const char* title,Int_t nBins, Double_t* xBins, const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH1D* result = new TH1D(name, title, nBins, xBins);
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  return result;
}
//________________________________________________________________________
TH2D* NeutronGenerator::CreateHist2D(const char* name, const char* title,Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t yMin, Double_t yMax, const char* xLabel, const char* yLabel, const char* zLabel)
{
  // create a histogram
  TH2D* result = new TH2D(name, title, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  result->SetOption("COLZ");
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  if (zLabel) result->GetZaxis()->SetTitle(zLabel);
  return result;
}
//________________________________________________________________________
Int_t NeutronGenerator::FromMatrixToVector(Int_t i, Int_t j)
{
  if (i <= j) return i * maxNeutrons - (i - 1) * i / 2 + j - i;
  else return j * maxNeutrons - (j - 1) * j / 2 + i - j;
}
//________________________________________________________________________
void NeutronGenerator::FromVectorToMatrix(Int_t index, Int_t &row, Int_t &col) 
{ 
  row = 0; 
  Int_t keyafter; 
  do 
  { 
      row++; 
      keyafter = row * maxNeutrons - (row - 1) * row / 2; 
  } while (index >= keyafter); 
  row--; 
  col = maxNeutrons - keyafter + index; 
}
