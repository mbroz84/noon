// --- Standard library ---
#include <iostream>
#include <fstream>
#include <string>

// --- ROOT system ---
#include "TGraph.h"
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

//_____________________________________________________________________________
NeutronGenerator::NeutronGenerator()
  : TObject()
  , fRunMode(kInterface)
  , fHadronicIntModel(kGlauber)
  , iEvent(0)
  , hInputRapidity(NULL)
  , hInputMass(NULL)
  , fRapMin(-666)
  , fRapMax(666)
  , fMassMin(0)
  , fMassMax(666)
  , lineString(0)
  , nFluxes(2+(maxNeutrons)*(maxNeutrons+1)/2)
  , nucleus_Z(82)
  , nucleus_A(208)
  , beamGamma(2942)
  , gammaTarget(2.0*beamGamma*beamGamma-1.0)
  , neutronSepThr(0.0)
  , saturationEnergy(1e6)
  , gSection_Nn(NULL)
  , hSection_Nn(NULL)
  , gPhotonFluxTable(NULL)
  , gNucleusBreakupTable(NULL)
  , hTwoPhotonFluxModulationTable(NULL)
  , gNuclearThickness(NULL)
  , gHadronicProbabilityTable(NULL)
  , hXsectionXn(NULL)
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
 for(Int_t i = 0; i<625; i++){
   energyGamma_Xn[i] = 0.0;
   xSection_Xn[i] = 0.0;
   }
}
//_____________________________________________________________________________
NeutronGenerator::~NeutronGenerator()
{
  //destructor
  if(hInputRapidity) {delete hInputRapidity; hInputRapidity = NULL;}
  if(hInputMass) {delete hInputMass; hInputMass = NULL;}
  if(gPhotonFluxTable) {delete gPhotonFluxTable; gPhotonFluxTable = NULL;}
  if(hTwoPhotonFluxModulationTable) {delete hTwoPhotonFluxModulationTable; hTwoPhotonFluxModulationTable = NULL;}
  if(gNuclearThickness) {delete gNuclearThickness; gNuclearThickness = NULL;}
  if(gHadronicProbabilityTable) {delete gHadronicProbabilityTable; gHadronicProbabilityTable = NULL;}
  if(hXsectionXn) {delete hXsectionXn; hXsectionXn = NULL;} 
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
    fInputStarlightAscii.open(filename);
    std::string newLine;
    getline(fInputStarlightAscii,newLine);
    fOutputStarlightAscii.open(name1);
    fOutputStarlightAscii<<newLine<<endl;
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
  FinishProduction();
}
//______________________________________________________________________________
void NeutronGenerator::GenerateEvent(const Double_t photonK)
{
  if(fRunMode != kFlatMultiplicity && fRunMode != k1n1n)hPhotonK->Fill(photonK);
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
    for(Int_t i = 0; i < maxNeutrons; i++){ 
      for(Int_t j = i; j < maxNeutrons; j++){
        if(i != j)hEventBreakupMap->SetBinContent(FromMatrixToVector(i,j)+1,2.0*GetBreakupProbability(photonK,i,j));
	if(i == j)hEventBreakupMap->SetBinContent(FromMatrixToVector(i,j)+1,GetBreakupProbability(photonK,i,j));
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
  //cout<<gPhotonFluxTable[0].Eval(photonK)<<endl;
  if(nNeutronsBeam1+nNeutronsBeam2 == -2)probability = gPhotonFluxTable[nFluxes-1].Eval(photonK)/gPhotonFluxTable[0].Eval(photonK);
  if(nNeutronsBeam1+nNeutronsBeam2 == -1)probability = 1-gPhotonFluxTable[nFluxes-1].Eval(photonK)/gPhotonFluxTable[0].Eval(photonK)-
  							 gPhotonFluxTable[FromMatrixToVector(0,0)+1].Eval(photonK)/gPhotonFluxTable[0].Eval(photonK);
  if(nNeutronsBeam1 >= 0 && nNeutronsBeam2 >= 0) probability = gPhotonFluxTable[FromMatrixToVector(nNeutronsBeam1,nNeutronsBeam2)+1].Eval(photonK)/gPhotonFluxTable[0].Eval(photonK);
  			   
  return probability;
}
//______________________________________________________________________________
Double_t NeutronGenerator::GetTotalFlux(const Double_t photonK)
{

  return gPhotonFluxTable[0].Eval(photonK);
}
//______________________________________________________________________________
void NeutronGenerator::FinishProduction(){

  if(fRunMode == kMassRapidity || fRunMode == kFlatMultiplicity || fRunMode == k1n1n || kStoreQA || kStoreGen)fOutputFile = new TFile("Output.root","RECREATE");
  if(kStoreQA)fQAhistList->Write();
  if(kStoreGen){
    if(gNuclearThickness)gNuclearThickness->Write();
    gHadronicProbabilityTable->Write();
    hXsectionXn->Write(); 
    gMeanNeutronN->Write();
    gWidthNeutronN->Write();
    fitMeanNeutronN->Write();
    fitWidthNeutronN->Write();
    gNucleusBreakupTable[maxNeutrons].Write();
    for(Int_t i = 0; i<10; i++) gNucleusBreakupTable[i].Write();
    gUnitaryLeak->SetName("gUnitaryLeak");
    gUnitaryLeak->Write();
    
    }
  if(fRunMode == kMassRapidity || fRunMode == kFlatMultiplicity || fRunMode == k1n1n) fEventTree->Write();
  if(fOutputFile)fOutputFile->Close();  
}
//______________________________________________________________________________
void NeutronGenerator::FinishEvent(){
  
  if(fRunMode == kMassRapidity || fRunMode == kFlatMultiplicity || fRunMode == k1n1n)fEventTree ->Fill();
  if(fRunMode == kStarlightAscii){
    for(Int_t i=0; i< fParticles->GetEntriesFast(); i++)fOutputStarlightAscii<<"TRACK: "<<i+11<<" "<<
    							       ((TParticle*)fParticles->At(i))->Px()<<" "<<
    							       ((TParticle*)fParticles->At(i))->Py()<<" "<<
							       ((TParticle*)fParticles->At(i))->Pz()<<" "<<iEvent+1<<" 0 0 "<<"2112"<<endl;
    fOutputStarlightAscii<<lineString.Data()<<endl;
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
      energyPhoton = hSection_Nn[nNeutrons[side]-1].GetRandom();
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
void NeutronGenerator::BuildPhotonFluxModulationTables()
{
  Double_t energy_min = 1e-5; 
  Double_t energy_max = 12.0 * beamGamma * hbarc/(2.0*nucleus_R);  
  Double_t energy_delta = TMath::Exp(TMath::Log(energy_max/energy_min)/Double_t(nSteps_Energy));
  
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
    Double_t* flux = PhotonFlux(energy_value);
    gPhotonFluxTable[0].SetPoint(iStep,energy_value,energy_value*flux[0]);
    for(Int_t i = 0; i < maxNeutrons; i++){
      for(Int_t j = i; j < maxNeutrons; j++){
        gPhotonFluxTable[FromMatrixToVector(i,j)+1].SetPoint(iStep,energy_value,energy_value*flux[FromMatrixToVector(i,j)+1]);
	}
      }
    gPhotonFluxTable[nFluxes-1].SetPoint(iStep,energy_value,energy_value*flux[nFluxes-1]);
    iStep++;
    if(iStep%10 == 0){
      if(iStep > 10){ printf("\033[1A"); printf("\033[K");}   
       cout<<iStep<<"%"<<endl;
       }
   }
}
//______________________________________________________________________________
Double_t *NeutronGenerator::PhotonFlux(const Double_t energyPhoton)
{
  Double_t flux_differential = 0.0;
  Double_t flux_differentialMod[nFluxes], flux_integralMod[nFluxes];
  Double_t *flux_integral = new Double_t[nFluxes];
  for(Int_t i=0; i<nFluxes; i++)flux_integral[i] = 0.0;
  Double_t dist;
  Double_t breakupProb[maxNeutrons+1];
 
  Double_t impactPar_min=1.8*nucleus_R;
  Double_t impactPar_max=impactPar_min + 6.0*hbarc*beamGamma/energyPhoton;  //6.0*adiabatic cutoff energy
  
  Double_t impactPar_delta=TMath::Exp(TMath::Log(impactPar_max/impactPar_min)/Double_t(nSteps_impactPar));
  Double_t R_delta = nucleus_R/Double_t(nSteps_R);
  Double_t Phi_delta=pi/Double_t(nSteps_Phi);

  for(Double_t impactPar_value = impactPar_min; impactPar_value<=impactPar_max; impactPar_value *= impactPar_delta){
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
    flux_differential *= 2.0*pi*impactPar_value*impactPar_value*(1.0-1.0/impactPar_delta)*gHadronicProbabilityTable->Eval(impactPar_value);
    //modulate by the probability of nuclear breakup as f(impactPar_value)
    for(Int_t i = 0; i < maxNeutrons+1; i++) breakupProb[i] = gNucleusBreakupTable[i].Eval(impactPar_value);    
    flux_differentialMod[nFluxes-1] = flux_differential*breakupProb[maxNeutrons]*breakupProb[maxNeutrons];
    flux_integral[nFluxes-1] += flux_differentialMod[nFluxes-1];
    for(Int_t i = 0; i < maxNeutrons; i++){
      for(Int_t j = i; j < maxNeutrons; j++){
        flux_differentialMod[FromMatrixToVector(i,j)] = flux_differential*breakupProb[i]*breakupProb[j];
	flux_integral[FromMatrixToVector(i,j)+1] += flux_differentialMod[FromMatrixToVector(i,j)];
	}
      }
    flux_integral[0] += flux_differential;
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
void NeutronGenerator::BuildNucleusBreakupProbabilityTable()
{
  gNucleusBreakupTable = new TGraph[maxNeutrons+1];
  for(Int_t i = 0; i < maxNeutrons; i++) gNucleusBreakupTable[i].SetName(TString::Format("gNucleusBreakupTable%d",i));
  gNucleusBreakupTable[maxNeutrons].SetName("gNucleusBreakupTableXn");

  cout<<"Building nucleus breakup probability tables"<<endl;  
  
  Double_t impactPar_min=1.8*nucleus_R;
  Double_t impactPar_max=impactPar_min + 6.0*hbarc*beamGamma/(1e-5); 
  Double_t impactPar_delta=TMath::Exp(TMath::Log(impactPar_max/impactPar_min)/Double_t(nSteps_impactPar*10));
  Int_t iStep = 0;
  
  for(Double_t impactPar_value = impactPar_min; impactPar_value<=impactPar_max; impactPar_value *= impactPar_delta){
    Double_t* prob_Breakup = NucleusBreakupProbability(impactPar_value);
    for(Int_t i = 0; i < maxNeutrons+1; i++)gNucleusBreakupTable[i].SetPoint(iStep,impactPar_value,prob_Breakup[i]);
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
  Double_t photonEnergyLimit_Nn = 140.0;
  if (beamGamma > 500.){
      photonEnergyLimit_Xn=1.E10;
  }
  else{
      photonEnergyLimit_Xn=1.E7;
  }
   Double_t gk1=0,gk1m=0,maxPhotonEnergy = 0;
   Double_t integralXn = 0, integralNn = 0;

  //Compute the probabilities 
  //Xn dissociation
  maxPhotonEnergy = TMath::Min(photonEnergyLimit_Xn,4.0*gammaTarget*hbarcmev/impactPar);
  gk1m = PhotonDensity2(impactPar,energyGamma_Xn[0]);
  for(Int_t k = 1; energyGamma_Xn[k] < maxPhotonEnergy; k++){
    gk1 = PhotonDensity2(impactPar,energyGamma_Xn[k]);
    prob_Xn += (energyGamma_Xn[k]-energyGamma_Xn[k-1])*0.5*(xSection_Xn[k-1]*gk1m+xSection_Xn[k]*gk1);
    for(Int_t i=1; i<maxNeutrons; i++) prob_Nn[i-1] += (energyGamma_Xn[k]-energyGamma_Xn[k-1])*0.5*(xSection_Xn[k-1]*GetBR(energyGamma_Xn[k-1],i)*gk1m+xSection_Xn[k]*GetBR(energyGamma_Xn[k],i)*gk1);
    gk1m = gk1;
    }
    
  prob_Breakup[maxNeutrons] = 1-TMath::Exp(-1*prob_Xn);
  prob_Breakup[0] = TMath::Exp(-1*prob_Xn);
  
  //Nn dissociation
  for(Int_t i=0; i<10; i++){
    photonEnergyLimit_Nn = gSection_Nn[i].GetX()[gSection_Nn[i].GetN()-1];
    maxPhotonEnergy = TMath::Min(photonEnergyLimit_Nn,4.0*gammaTarget*hbarcmev/impactPar); 
    gk1m = PhotonDensity2(impactPar,gSection_Nn[i].GetX()[0]);
    for(Int_t k = 1; gSection_Nn[i].GetX()[k] < maxPhotonEnergy; k++){
      gk1 = PhotonDensity2(impactPar,gSection_Nn[i].GetX()[k]);
      prob_Nn[i] += (gSection_Nn[i].GetX()[k]-gSection_Nn[i].GetX()[k-1])*0.5*(gSection_Nn[i].GetY()[k-1]*gk1m+gSection_Nn[i].GetY()[k]*gk1);
      gk1m = gk1;
      }
    }
  
  Double_t integralPoisson = 0;
  Int_t maxExcitation = 5;
  if(impactPar<22) maxExcitation = 6;
  /*/
  for(Int_t i=0; integralPoisson<0.999; i++){
    integralPoisson += TMath::PoissonI(i,prob_Xn);
    maxExcitation = i;
    }
/*/
  Double_t norm_Xn = 0.0;  
  for(Int_t i=1; i<maxNeutrons; i++){ 
    prob_Breakup[i] += prob_Nn[i-1];
    for(Int_t j=1; i+j<maxNeutrons; j++){
    if(maxExcitation<2) break;
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

  for(Int_t i=1; i<maxNeutrons; i++) prob_Breakup[i] *= prob_Breakup[maxNeutrons]/norm_Xn;
  gUnitaryLeak->SetPoint(gUnitaryLeak->GetN(),impactPar,norm_Xn/prob_Breakup[maxNeutrons]);
  
  return prob_Breakup;  
}

//______________________________________________________________________________
Double_t NeutronGenerator::GetBR(const Double_t energyPhoton, const Int_t nNeutrons)
{
  if(nNeutrons < 2) return 0;
  if(energyPhoton < nNeutrons*neutronSepThr) return 0;
  if(nNeutrons <= 10)if(energyPhoton < gSection_Nn[nNeutrons-1].GetX()[gSection_Nn[nNeutrons-1].GetN()-1]) return 0;

  return hBranchingRatioMap->GetBinContent(hBranchingRatioMap->FindBin(energyPhoton,nNeutrons));
} 
//______________________________________________________________________________
void NeutronGenerator::Initialize()
{
  //Data for Xn in 25GeV-103GeV from Pb Nucl. Phys. A367, 237 (1981) 		
  Double_t energyGamma_Lepretre[28]={26.,28.,30.,32.,34.,36.,38.,40.,44.,46.,48.,50.,52.,55.,57.,62.,64.,66.,69.,72.,74.,76.,79.,82.,86.,92.,98.,103.};
  Double_t xSection_Lepretre[28]={30.,21.5,22.5,18.5,17.5,15.,14.5,19.,17.5,16.,14.,20.,16.5,17.5,17.,15.5,18.,15.5,15.5,15.,13.5,18.,14.5,15.5,12.5,13.,13.,12.};

  //Data for Xn in 103-440 MeV Nucl. Phys. A431, 573 (1984)
  Double_t energyGamma_Carlos[22]= {103.,106.,112.,119.,127.,132.,145.,171.,199.,230.,235.,
  		  254.,280.,300.,320.,330.,333.,373.,390.,420.,426.,440.};
  Double_t xSection_Carlos[22]= {12.0,11.5,12.0,12.0,12.0,15.0,17.0,28.0,33.0,
  		  52.0,60.0,70.0,76.0,85.0,86.0,89.0,89.0,75.0,76.0,69.0,59.0,61.0};
		  
  //Data for Xn in 2-16.4 GeV from Phys. Rev. Lett. 39, 737 (1977) and Phys. Rev. D 7, 1362 (1973)		  
  Double_t energyGamma_MichaDwell[11]={2000.0,3270.0,4100.0,4810.0,6210.0,6600.0,7790.0,8400.0,9510.0,13600.0,16400.0};
  Double_t xSection_MichaDwell[11]={0.1266,0.1080,0.0805,0.1017,0.0942,0.0844,0.0841,0.0755,0.0827,0.0626,0.0740};
		 
  
  // gammay,p gamma,n of Armstrong begin at 265 incr 25
  Double_t xSection_p_Armstrong[160]={0.,.4245,.4870,.5269,.4778,.4066,.3341,.2444,.2245,.2005,
  		    .1783,.1769,.1869,.1940,.2117,.2226,.2327,.2395,.2646,.2790,.2756,
  		    .2607,.2447,.2211,.2063,.2137,.2088,.2017,.2050,.2015,.2121,.2175,
  		    .2152,.1917,.1911,.1747,.1650,.1587,.1622,.1496,.1486,.1438,.1556,
  		    .1468,.1536,.1544,.1536,.1468,.1535,.1442,.1515,.1559,.1541,.1461,
  		    .1388,.1565,.1502,.1503,.1454,.1389,.1445,.1425,.1415,.1424,.1432,
  		    .1486,.1539,.1354,.1480,.1443,.1435,.1491,.1435,.1380,.1317,.1445,
  		    .1375,.1449,.1359,.1383,.1390,.1361,.1286,.1359,.1395,.1327,.1387,
  		    .1431,.1403,.1404,.1389,.1410,.1304,.1363,.1241,.1284,.1299,.1325,
  		    .1343,.1387,.1328,.1444,.1334,.1362,.1302,.1338,.1339,.1304,.1314,
  		    .1287,.1404,.1383,.1292,.1436,.1280,.1326,.1321,.1268,.1278,.1243,
  		    .1239,.1271,.1213,.1338,.1287,.1343,.1231,.1317,.1214,.1370,.1232,
  		    .1301,.1348,.1294,.1278,.1227,.1218,.1198,.1193,.1342,.1323,.1248,
  		    .1220,.1139,.1271,.1224,.1347,.1249,.1163,.1362,.1236,.1462,.1356,
  		    .1198,.1419,.1324,.1288,.1336,.1335,.1266};


  Double_t xSection_n_Armstrong[160]={0.,.3125,.3930,.4401,.4582,.3774,.3329,.2996,.2715,.2165,
  		     .2297,.1861,.1551,.2020,.2073,.2064,.2193,.2275,.2384,.2150,.2494,
  		     .2133,.2023,.1969,.1797,.1693,.1642,.1463,.1280,.1555,.1489,.1435,
  		     .1398,.1573,.1479,.1493,.1417,.1403,.1258,.1354,.1394,.1420,.1364,
  		     .1325,.1455,.1326,.1397,.1286,.1260,.1314,.1378,.1353,.1264,.1471,
  		     .1650,.1311,.1261,.1348,.1277,.1518,.1297,.1452,.1453,.1598,.1323,
  		     .1234,.1212,.1333,.1434,.1380,.1330,.12,.12,.12,.12,.12,.12,.12,.12,
  		     .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
  		     .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
  		     .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
  		     .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
  		     .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12};

  Double_t BrtWgnr_xsection = 0, BrtWgnr_width = 0, BrtWgnr_mean = 0;
  Int_t	nStepsLowE = 0, iStepE =0;
  Double_t eStepLowE = 0.05;
  
  if (nucleus_Z == 79){
    //GDR Breit-Wigner fit parameters from Nucl. Phys. A159, 561 (1970)  
    BrtWgnr_xsection=540.0;
    BrtWgnr_width=4.75;
    BrtWgnr_mean=13.70;
    neutronSepThr=8.1;
      }
  else{
    //GDR Breit-Wigner fit parameters from Nucl. Phys. A159, 561 (1970)
    BrtWgnr_xsection=640.0;
    BrtWgnr_width=4.05;
    BrtWgnr_mean=13.42;
    neutronSepThr=7.4;
    }

  //Veyssiere et al. Nucl. Phys. A159, 561 (1970) - Lorentzian parametrization of GDR peak
  Double_t BrtWgnr_Con = 0.1*BrtWgnr_width*BrtWgnr_width*BrtWgnr_xsection;
  nStepsLowE=Int_t((25.0-neutronSepThr)/eStepLowE)+1;
  for ( Int_t i = 0; i < nStepsLowE; i++ ) {
    energyGamma_Xn[i] = i*eStepLowE + neutronSepThr;
    xSection_Xn[i] = 0.93*BrtWgnr_Con*energyGamma_Xn[i]*energyGamma_Xn[i]/(((BrtWgnr_mean*BrtWgnr_mean-energyGamma_Xn[i]*energyGamma_Xn[i])*(BrtWgnr_mean*BrtWgnr_mean-energyGamma_Xn[i]*energyGamma_Xn[i]))
  				+energyGamma_Xn[i]*energyGamma_Xn[i]*BrtWgnr_width*BrtWgnr_width);
  }
  
  iStepE = nStepsLowE;   
  //25-103 MeV, Lepretre, et al., Nucl. Phys. A367, 237 (1981)
  for ( Int_t j = 0; j < 27; j++ ) {
    energyGamma_Xn[iStepE] = energyGamma_Lepretre[j];
    xSection_Xn[iStepE] = 0.1*nucleus_A*xSection_Lepretre[j]/208.0;
    iStepE++; 
    }
  //103-440 MeV, Carlos, et al., Nucl. Phys. A431, 573 (1984)
  for ( Int_t j = 0; j < 22; j++ ) {
    energyGamma_Xn[iStepE] = energyGamma_Carlos[j];
    xSection_Xn[iStepE] = 0.1*nucleus_A*xSection_Carlos[j]/208.0;
    iStepE++;
    }
    
  //440 MeV-2 GeV Armstrong et al.
  for ( Int_t j = 9; j <= 70; j++) {
      energyGamma_Xn[iStepE] = energyGamma_Xn[iStepE-1]+25.0;
      xSection_Xn[iStepE] = 0.1*(nucleus_Z*xSection_p_Armstrong[j]+(nucleus_A-nucleus_Z)*xSection_n_Armstrong[j]);
      iStepE++;
      }
      
  //2-16.4 GeV Michalowski; Caldwell
  for ( Int_t j = 0; j < 11; j++) {
      energyGamma_Xn[iStepE] = energyGamma_MichaDwell[j];
      xSection_Xn[iStepE] = 0.1*nucleus_A*xSection_MichaDwell[j];
      iStepE++;
      }
  //Regge parameters
  Double_t x=0,y=0,eps=0,eta=0,em=0,exx=0,s=0,ictr=0,pom=0,vec=0;
  x = 0.0677;
  y = 0.129;
  eps = 0.0808;
  eta = 0.4525;
  em = 0.94;
  exx = pow(10,0.05);

  //Regge model for high energy
  s = 0.002*em*energyGamma_Xn[iStepE-1];
  //make sure we reach LHC energies
  ictr = 100;
  if ( gammaTarget > (2.*150.*150.)) ictr = 150;
  for ( Int_t j = 1; j <= ictr; j++ ) {
      s = s*exx;
      energyGamma_Xn[iStepE] = 1000.0*0.5*(s-em*em)/em;
      pom = x*pow(s,eps);
      vec = y*pow(s,(-eta));
      xSection_Xn[iStepE] = 0.1*0.65*nucleus_A*(pom+vec);
      iStepE++;
  }
  hXsectionXn = new TH1D("hXsectionXn","hXsectionXn",624,energyGamma_Xn);
  for(Int_t i = 0; i<624; i++)hXsectionXn->SetBinContent(i+1,xSection_Xn[i]);
  //cout<<"Hist = "<<hXsectionXn->Integral(0,hXsectionXn->FindBin(140))<<endl;
  
  Double_t eMeanN_Lepretre[34] = {24.38, 27.07, 28.96, 31.89, 33.51, 34.3, 35.39, 38.06, 40.18, 45.28, 46.09, 47.42, 51.17, 54.9, 58.1, 64.26, 62.67, 65.88, 70.16, 72.57, 75.51, 77.09, 79.24, 81.91, 85.96, 92.63, 98.25, 102.81, 106.29, 112.16, 119.12, 127.41, 132.45, 140.76};

  Double_t vMeanN_Lepretre[34] = {2.057, 2.515, 3.055, 2.948, 3.38, 3.165, 3.785, 3.543, 3.248, 3.708, 4.085, 3.978, 4.276, 3.873, 3.793, 4.065, 4.415, 4.497, 4.795, 4.903, 4.985, 4.42, 4.69, 4.503, 5.932, 5.477, 5.856, 6.316, 6.721, 6.427, 6.78, 7.133, 5.923, 6.734};

  Double_t eWidthN_Lepretre[31] = {24.05, 26.45, 29.12, 32.06, 34.2, 37.94, 39.54, 43.82, 45.15, 46.76, 49.96, 51.03, 54.24, 57.71, 61.45, 64.12, 65.19, 66.26, 71.87, 74.54, 76.68, 79.08, 81.76, 86.3, 92.18, 97.79, 102.33, 111.68, 118.89, 126.37, 131.72};

  Double_t vWidthN_Lepretre[31] = {0.308, 0.521, 0.215, 0.415, 0.402, 0.708, 1, 1.001, 1.081, 1.347, 1.413, 1.307, 1.719, 1.295, 1.694, 1.906, 1.707, 1.694, 1.92, 2.319, 2.2, 1.881, 2.121, 2.108, 2.122, 2.508, 2.229, 2.522, 2.111, 2.205, 3.308};

  Double_t eMeanN_Veyssiere[39] = {11.63, 11.97, 12.25, 12.53, 12.80, 13.10, 13.35, 13.65, 13.92, 14.18, 14.46, 14.74, 15.01, 15.28, 15.54, 15.82, 16.09, 16.37, 16.63, 16.91, 17.18, 17.17, 17.46, 17.71, 17.98, 18.26, 18.52, 18.78, 19.07, 19.35, 19.67, 19.92, 20.22, 20.48, 20.76, 21.02, 21.32, 21.58, 21.87}; 

  Double_t vMeanN_Veyssiere[39] = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.04, 1.09, 1.15, 1.20, 1.28, 1.31, 1.39, 1.45, 1.49, 1.55, 1.56, 1.61, 1.63, 1.64, 1.71, 1.68, 1.72, 1.72, 1.74, 1.77, 1.69, 1.70, 1.64, 1.76, 1.74, 1.73, 1.70, 1.66, 1.73};

  Double_t thresholdNeutronN[4] = {7.4, 14.1, 22.5, 30.0};

  Double_t eMeanNeutronN[74],vMeanNeutronN[74]; 
  Double_t eWidthNeutronN[34],vWidthNeutronN[34]; 
  //Double_t eWidthNeutronN[71],vWidthNeutronN[71];
  
  eMeanNeutronN[0] = neutronSepThr;
  vMeanNeutronN[0] = 1.0;
  iStepE = 1;
  for ( Int_t j = 0; j < 39; j++) {
    eMeanNeutronN[iStepE] = eMeanN_Veyssiere[j];
    vMeanNeutronN[iStepE] = vMeanN_Veyssiere[j];
    iStepE++;
    }
  for ( Int_t j = 0; j < 34; j++) {
    eMeanNeutronN[iStepE] = eMeanN_Lepretre[j];
    vMeanNeutronN[iStepE] = vMeanN_Lepretre[j];
    iStepE++;
    }
  
  //There are no data for widths from Veyssiere, we use uncertainty of the mean as the width in the 2neutron region. 
  //   
  eWidthNeutronN[0] = thresholdNeutronN[0];
  vWidthNeutronN[0] = 0.1;
  eWidthNeutronN[1] = thresholdNeutronN[1];
  vWidthNeutronN[1] = 0.1;
  eWidthNeutronN[2] = thresholdNeutronN[2];
  vWidthNeutronN[2] = 0.1;
  iStepE = 3;
  /*/
  for ( Int_t j = 0; j < 39; j++) {
    eWidthNeutronN[iStepE] = eMeanN_Veyssiere[j];
    vWidthNeutronN[iStepE] = TMath::Sqrt(vMeanN_Veyssiere[j]);
    iStepE++;
    }
  /*/  
  for ( Int_t j = 0; j < 31; j++) {
    eWidthNeutronN[iStepE] = eWidthN_Lepretre[j];
    vWidthNeutronN[iStepE] = vWidthN_Lepretre[j];
    iStepE++;
    }
  
 gMeanNeutronN = new TGraph(74,eMeanNeutronN,vMeanNeutronN);
 gMeanNeutronN->SetName("gMeanNeutronN");
 gWidthNeutronN = new TGraph(34,eWidthNeutronN,vWidthNeutronN);
 gWidthNeutronN->SetName("gWidthNeutronN");
 
 fitMeanNeutronN = new TF1("fitMeanNeutronN", "[0]*TMath::Log(x)+[1]", 140, 1e9);
 fitMeanNeutronN->SetParameters(2.16117,-4.47756);
 
 fitWidthNeutronN = new TF1("fitWidthNeutronN", "[0]*TMath::Log(x)+[1]", 140, 1e9);
 fitWidthNeutronN->SetParameters(1.53,-4.68);
 
 Double_t mean; 
 Double_t width;
 Double_t integral;
 Double_t energyValue;
 hBranchingRatioMap = new TH2D("hBranchingRatioMap","hBranchingRatioMap",624,energyGamma_Xn,50,0,50);
 for(Int_t i = 0; i<624; i++){
   if(hXsectionXn->GetBinCenter(i+1)<100)continue;
   energyValue = hXsectionXn->GetBinCenter(i+1);
   if(energyValue > saturationEnergy) energyValue = saturationEnergy;
   mean = fitMeanNeutronN->Eval(energyValue);
   width = fitWidthNeutronN->Eval(energyValue);
   integral = 0;
   for(Int_t j = 2; j<50; j++){
     hBranchingRatioMap->SetBinContent(i+1,j+1,TMath::Gaus(j,mean,width,kTRUE));
     integral+=TMath::Gaus(j,mean,width,kTRUE);
     }
   for(Int_t j = 2; j<50; j++) hBranchingRatioMap->SetBinContent(i+1,j+1,hBranchingRatioMap->GetBinContent(i+1,j+1)/integral);
   }
 
 gSection_Nn = new TGraph[10];
 
 Double_t energyGamma_1n[76] = {7.5, 7.57, 7.64, 7.71, 7.79, 7.88, 7.95, 8.05, 8.12, 8.19, 8.26, 8.33, 8.47, 8.56, 8.65, 8.74, 8.87, 8.98, 9.12, 9.23, 9.35, 9.49, 9.6, 9.74, 9.84, 9.96, 10.08, 10.22, 10.37, 10.49, 10.61, 10.73, 10.85, 11.0, 11.13, 11.24, 11.38, 11.51, 11.64, 11.76, 11.89, 12.0, 12.13, 12.26, 12.53, 12.68, 12.95, 13.12, 13.5, 13.79, 14.05, 14.28, 14.58, 14.8, 15.08, 15.37, 15.66, 15.95, 16.21, 16.53, 16.79, 17.05, 17.33, 17.6, 17.86, 18.13, 18.4, 18.66, 18.9, 19.2, 19.46, 19.74, 20.03, 20.28, 20.55, 20.85}; 

  Double_t xSection_1n[150] = {11.0, 21.0, 34.0, 31.0, 21.0, 19.0, 32.0, 36.0, 29.0, 25.0, 26.0, 35.0, 31.0, 34.0, 32.0, 38.0, 38.0, 46.0, 54.0, 60.0, 71.0, 73.0, 67.0, 80.0, 115.0, 134.0, 113.0, 101.0, 132.0, 164.0, 176.0, 183.0, 208.0, 241.0, 306.0, 311.0, 302.0, 321.0, 336.0, 386.0, 407.0, 418.0, 432.0, 459.0, 512.0, 534.0, 589.0, 641.0, 645.0, 625.0, 598.0, 560.0, 470.0, 418.0, 336.0, 287.0, 230.0, 184.0, 137.0, 116.0, 95.0, 82.0, 71.0, 59.0, 53.0, 46.0, 42.0, 35.0, 31.0, 27.0, 26.0, 25.0, 25.0, 23.0, 23.0, 19.0}; 

  gSection_Nn[0] = TGraph(76,energyGamma_1n,xSection_1n);
  for (Int_t i=0;i<gSection_Nn[0].GetN();i++) gSection_Nn[0].GetY()[i] *= 0.1;
  for (Int_t i=0;i<gSection_Nn[0].GetN();i++) gSection_Nn[0].GetY()[i] *= 0.93; //Correction from Phys.Rev. C 36, 1286 (1987)
  gSection_Nn[0].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[0].SetMarkerStyle(kFullCircle);
  gSection_Nn[0].SetMarkerColor(kBlue);
 
  Double_t energyGamma_2nV[50] = {14.587, 14.758, 14.915, 15.007, 15.092, 15.281, 15.386, 15.602, 15.836, 16.077, 16.316, 16.590, 16.801, 17.005, 17.226, 17.447, 17.682, 17.933, 18.115, 18.303, 18.504, 18.822, 19.067, 19.308, 19.545, 19.776, 20.116, 20.344, 20.575, 20.736, 20.951, 21.187, 21.471, 21.685, 21.930, 22.270, 22.500, 22.790, 23.030, 23.255, 23.486, 23.711, 23.919, 24.151, 24.375, 24.600, 24.808, 25.040, 25.475, 25.708};
  Double_t xSection_2nV[50] = {0.904, 1.809, 2.932, 3.551, 4.267,  5.114, 6.017, 6.850, 7.613, 7.988, 8.537, 8.664, 8.811, 8.637, 8.493, 8.210, 8.139, 7.693, 7.574, 7.471, 7.128, 6.798, 6.290, 6.124, 5.603, 5.313, 5.029, 5.037, 4.724, 4.521, 4.602, 4.532, 4.507, 4.435, 4.270, 4.247, 4.023, 3.897, 3.887, 3.793, 3.632, 3.402, 3.262, 3.122, 3.041, 2.842, 2.535, 2.244, 1.778, 1.589};

  Double_t energyGamma_2nL[40] = {26.001, 26.141, 26.149, 26.435, 27.091, 27.098, 27.885, 29.022, 29.927, 31.763, 34.364, 37.141, 42.241, 45.281, 48.289, 51.868, 55.422, 58.917, 62.668, 64.688, 68.408, 73.211, 76.919, 84.075, 87.430, 91.308, 92.320, 95.786, 99.933, 102.943, 106.520, 107.693, 111.102, 115.628, 120.750, 123.759, 127.338, 128.881, 131.890, 135.306};
  Double_t xSection_2nL[40] = { 30.110, 29.326, 27.974, 26.828, 25.854, 24.851, 23.790, 22.571, 21.159, 19.817, 18.783, 17.806, 16.979, 16.985, 16.244, 16.008, 15.672, 15.403, 15.350, 14.835, 14.566, 14.384, 14.240, 13.961, 13.602, 13.489, 13.044, 13.429, 13.338, 13.054, 13.032, 12.714, 12.784, 12.773, 12.817, 12.589, 12.285, 12.596, 12.424, 12.303};
  
  gSection_Nn[1] = TGraph(50,energyGamma_2nV,xSection_2nV);
  for (Int_t i=0;i<gSection_Nn[1].GetN();i++) gSection_Nn[1].GetY()[i] *= 0.93; //Correction from Phys.Rev. C 36, 1286 (1987)
  gSection_Nn[1].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[1].SetMarkerStyle(kFullCircle);
  gSection_Nn[1].SetMarkerColor(kGreen);
  
  for(Int_t i = 0; i< 40; i++)gSection_Nn[1].SetPoint(50+i,energyGamma_2nL[i],xSection_2nL[i]*0.1/2);

  Double_t energyGamma_3n[44] = {22.515, 22.994, 23.148, 23.305, 23.623, 23.834,  23.935, 24.255, 24.733, 25.308, 26.009, 26.813, 28.595, 30.444, 32.659, 34.611, 37.150, 37.875, 40.211, 43.416, 46.641, 49.378, 53.185, 56.672, 60.200, 63.925, 67.504, 70.432, 74.660, 77.751, 82.468, 85.720, 90.274, 93.852, 98.730, 101.750, 106.210, 109.535, 114.341, 117.593, 122.471, 125.561, 129.951, 133.366}; 
  Double_t xSection_3n[44] = {3.124, 4.637, 5.964, 6.951, 8.114, 9.123, 10.271, 11.182, 12.646, 14.889, 16.752, 18.192, 19.338, 20.086, 19.558, 18.536, 17.781, 17.352, 16.573, 15.312, 14.790, 13.715, 13.306, 12.812, 12.410, 11.949, 11.669, 11.510, 11.398, 11.123, 10.971, 10.946, 10.807, 10.716, 10.696, 10.655, 10.583, 10.581, 10.459, 10.508, 10.456, 10.483, 10.460, 10.457};  
 
  Double_t energyGamma_4n[35] = {31.473, 32.280, 33.087, 34.054, 35.344, 35.826, 36.469, 37.276, 38.830, 40.193, 45.066, 48.583, 52.221, 55.706, 59.377, 62.954, 66.393, 69.840, 73.689, 80.844, 84.421, 88.091, 91.576, 94.502, 98.556, 102.261, 105.979, 110.765, 114.927, 118.281, 121.124, 123.813, 127.902, 131.101, 134.206};
  Double_t xSection_4n[35] = {0.838, 1.718, 2.775, 4.111, 5.857, 6.894, 8.109, 9.103, 10.259, 10.645, 11.446, 11.615, 11.435, 11.333, 11.200, 11.156, 11.056, 10.838, 10.695, 10.645, 10.732, 10.687, 10.580, 10.654, 10.450, 10.129, 9.925, 10.244, 10.346, 10.204, 9.874, 10.085, 9.906, 9.961, 9.851};

  Double_t energyGamma_5n[30] = {39.137, 40.408, 41.704, 43.161, 45.105, 46.920, 48.867, 52.326, 55.655, 60.369, 63.782, 68.496, 71.423, 75.975, 79.551, 82.153, 86.705, 90.282, 95.160, 98.737, 102.856, 105.891, 109.288, 113.046, 115.973, 120.525, 123.858, 127.678, 131.256, 134.345}; 
  Double_t xSection_5n[30] = {0.699, 2.029, 2.950, 3.920, 5.131, 6.049, 6.819, 7.618, 8.208, 8.513, 8.683, 8.941, 9.056, 9.201, 9.292, 9.331, 9.456, 9.495, 9.542, 9.559, 9.546, 9.547, 9.563, 9.551, 9.553, 9.710, 9.725, 9.845, 9.779, 9.814};
 
  Double_t energyGamma_6n[26] = {46.139, 47.999, 50.162, 54.058, 57.885, 61.042, 65.158, 68.191, 71.674, 75.342, 78.268, 83.144, 86.558, 91.922, 95.498, 99.660, 102.651, 105.414, 110.129, 113.380, 118.258, 121.834, 125.690, 128.988, 132.686, 135.166}; 
  Double_t xSection_6n[26] = {0.363, 1.617, 2.388, 3.561, 4.416, 4.891, 5.348, 5.642, 5.942, 6.245, 6.460, 6.792, 6.952, 7.209, 7.404, 7.573, 7.668, 7.809, 7.977, 8.074, 8.197, 8.329, 8.371, 8.440, 8.559, 8.605};
    
  Double_t energyGamma_7n[20] = {58.466, 62.292, 65.342, 69.187, 73.225, 76.173, 80.885, 84.094, 88.036, 91.449, 96.487, 100.063, 105.264, 108.352, 113.066, 116.154, 120.998, 124.281, 129.321, 132.247}; 
  Double_t xSection_7n[20] = {0.599, 1.326, 1.871, 2.467, 3.056, 3.389, 3.957, 4.275, 4.530, 4.836, 5.259, 5.484, 5.818, 6.036, 6.301, 6.485, 6.762, 6.947, 7.157, 7.282};

  Double_t energyGamma_8n[19] = {66.435, 71.147, 75.264, 77.972, 82.522, 85.512, 89.673, 93.443, 96.824, 102.350, 105.113, 109.248, 112.264, 116.652, 119.578, 124.551, 127.380, 131.736, 133.718}; 
  Double_t xSection_8n[19] = {0.369, 0.959, 1.399, 1.639, 2.095, 2.379, 2.716, 3.038, 3.315, 3.730, 3.971, 4.309, 4.445, 4.779, 5.003, 5.356, 5.515, 5.802, 5.962};

  Double_t energyGamma_9n[18] = {74.077, 78.116, 80.740, 85.127, 88.432, 92.279, 95.854, 99.431, 103.192, 106.420, 111.135, 114.061, 118.288, 121.214, 125.441, 128.652, 132.594, 134.871};
  Double_t xSection_9n[18] = {0.388, 0.909, 1.121, 1.579, 1.788, 2.070, 2.303, 2.490, 2.698, 2.873, 3.054, 3.130, 3.281, 3.356, 3.489, 3.536, 3.653, 3.659};
  
  Double_t energyGamma_10n[16] = {81.558, 85.134, 88.709, 92.285, 95.861, 99.437, 102.986, 106.589, 110.165, 113.742, 117.318, 120.895, 124.471, 128.047, 131.624, 134.551}; 
  Double_t xSection_10n[18] = {0.186, 0.418, 0.760, 1.048, 1.295, 1.447, 1.642, 1.851, 2.036, 2.162, 2.269, 2.431, 2.607, 2.691, 2.775, 2.831};
  
  gSection_Nn[2] = TGraph(44,energyGamma_3n,xSection_3n);
  for (Int_t i=0;i<gSection_Nn[2].GetN();i++) gSection_Nn[2].GetY()[i] *= 0.1/3;
  gSection_Nn[2].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[2].SetMarkerStyle(kFullCircle);
  gSection_Nn[2].SetMarkerColor(kRed);
  
  gSection_Nn[3] = TGraph(35,energyGamma_4n,xSection_4n);
  for (Int_t i=0;i<gSection_Nn[3].GetN();i++) gSection_Nn[3].GetY()[i] *= 0.1/4;
  gSection_Nn[3].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[3].SetMarkerStyle(kOpenCircle);
  gSection_Nn[3].SetMarkerColor(kRed);

  gSection_Nn[4] = TGraph(30,energyGamma_5n,xSection_5n);
  for (Int_t i=0;i<gSection_Nn[4].GetN();i++) gSection_Nn[4].GetY()[i] *= 0.1/5;
  gSection_Nn[4].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[4].SetMarkerStyle(kOpenCircle);
  gSection_Nn[4].SetMarkerColor(kGreen);
  
  gSection_Nn[5] = TGraph(26,energyGamma_6n,xSection_6n);
  for (Int_t i=0;i<gSection_Nn[5].GetN();i++) gSection_Nn[5].GetY()[i] *= 0.1/6;
  gSection_Nn[5].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[5].SetMarkerStyle(kOpenCircle);
  gSection_Nn[5].SetMarkerColor(kBlue);
  
  gSection_Nn[6] = TGraph(20,energyGamma_7n,xSection_7n);
  for (Int_t i=0;i<gSection_Nn[6].GetN();i++) gSection_Nn[6].GetY()[i] *= 0.1/7;
  gSection_Nn[6].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[6].SetMarkerStyle(kOpenCircle);
  gSection_Nn[6].SetMarkerColor(kBlack);
  
  gSection_Nn[7] = TGraph(19,energyGamma_8n,xSection_8n);
  for (Int_t i=0;i<gSection_Nn[7].GetN();i++) gSection_Nn[7].GetY()[i] *= 0.1/8;
  gSection_Nn[7].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[7].SetMarkerStyle(kOpenSquare);
  gSection_Nn[7].SetMarkerColor(kBlack);

  gSection_Nn[8] = TGraph(18,energyGamma_9n,xSection_9n);
  for (Int_t i=0;i<gSection_Nn[8].GetN();i++) gSection_Nn[8].GetY()[i] *= 0.1/9;
  gSection_Nn[8].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[8].SetMarkerStyle(kOpenStar);
  gSection_Nn[8].SetMarkerColor(kBlack);
  
  gSection_Nn[9] = TGraph(16,energyGamma_10n,xSection_10n);
  for (Int_t i=0;i<gSection_Nn[9].GetN();i++) gSection_Nn[9].GetY()[i] *= 0.01;
  gSection_Nn[9].GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  gSection_Nn[9].SetMarkerStyle(kFullStar);
  gSection_Nn[9].SetMarkerColor(10);

  hSection_Nn = new TH1D[10];
  for(Int_t i = 0; i<10; i++){
    hSection_Nn[i] = TH1D(TString::Format("hSection_Nn%d",i)," ",gSection_Nn[i].GetN()-1,gSection_Nn[i].GetX());
    for(Int_t iBin = 1; iBin <= hSection_Nn[i].GetNbinsX(); iBin++)hSection_Nn[i].SetBinContent(iBin,gSection_Nn[i].Eval(hSection_Nn[i].GetBinCenter(iBin)));
    }

  if(fRunMode != kFlatMultiplicity && fRunMode != k1n1n){
    BuildHadronicInteractionProbabilityTable();
    BuildNucleusBreakupProbabilityTable();
    BuildPhotonFluxModulationTables();
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
    //cout<<impactPar_value<<"  "<<HadronicInteractionProbability(impactPar_value)<<endl;
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
  if(fHadronicIntModel == kHardSphere)prob_Handronic = impactPar<2*1.2*TMath::Power(nucleus_A,1.0/3) ? 1e20 : 0;

  if(fHadronicIntModel == kGlauber){
    Double_t nucleus_WSSD = 0.535; //Woods-saxon skin depth
    Double_t energyCMS = 2*beamGamma*0.938; 
    // This equation is from section 50 of the particle data book, the subsection on "Total Hadronic Cross-Sections, using the parameterization for sqrt{s} > 7 GeV.
    // only the first and second terms contribute significantly, but leave them all here for good measure
    Double_t xSection_IntNN = 0.1*0.2838*TMath::Power(TMath::Log(energyCMS),2)+33.73+13.67*TMath::Power(energyCMS,-0.412)-7.77*TMath::Power(energyCMS,-0.5626);

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
      for(Double_t X_value = -maxDistance; X_value <= maxDistance; X_value += XY_delta){
        Double_t dist_bR = TMath::Sqrt((impactPar-X_value)*(impactPar-X_value)+Y_value*Y_value);
        Double_t dist_R = TMath::Sqrt(X_value*X_value+Y_value*Y_value);
        if(dist_bR >= maxDistance || dist_R >= maxDistance) continue;
        //cout<<dist_bR<<" "<<gNuclearThickness->Eval(dist_bR)<<" "<<dist_R<<" "<<gNuclearThickness->Eval(dist_R)<<endl;

        prob_Handronic += 2.0*gNuclearThickness->Eval(dist_bR)*(1.0-TMath::Exp(-xSection_IntNN*gNuclearThickness->Eval(dist_R)))*XY_delta*XY_delta;
        }
      }
    }
  return prob_Handronic;
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
  hEnergyForNeutronMulti = CreateHist1D("hEnergyForNeutronMulti","hEnergyForNeutronMulti",624,energyGamma_Xn,"Energy [MeV]","Counts");
  fQAhistList->Add(hEnergyForNeutronMulti);
  
  hRapidityVM = CreateHist1D("hRapidityVM","Rapidity of VM used in event generation",1000,-10,10,"Rapidity","Counts");
  fQAhistList->Add(hRapidityVM);
  hMassVM = CreateHist1D("hMassVM","Mass of VM used in event generation",1000,0,10,"Rapidity","Counts");
  fQAhistList->Add(hMassVM);
  
  if(fRunMode != kFlatMultiplicity && fRunMode != k1n1n){
    hPhotonK = CreateHist1D("hPhotonK","Virtual photon k=0.5*M_{VM}*exp(y_{VM}) used in event generation",nSteps_Energy-1,gPhotonFluxTable[0].GetX(),"k [GeV/c]","Counts");
    fQAhistList->Add(hPhotonK);
    }
    
  hProbabilityXn = CreateHist1D("hProbabilityXn","Single side probability of Xn",1000,0,1,"Probability","Counts");
  fQAhistList->Add(hProbabilityXn);

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
void NeutronGenerator::ReadENDF(Bool_t saveToFile)
{
  std::string newLine;
  TString token;
  Ssiz_t from = 0;
  Ssiz_t pos = 0;
  
  TString energy;
  Double_t energyPhoton = 0;
  Double_t energyNeutron = 0, probability = 0; 
    
  //"Intepreted" ENDF file from https://www-nds.iaea.org/exfor/endf.htm
  ifstream infile("ENDF_Pb208.txt");
  
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
