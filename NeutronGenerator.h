#ifndef NEUTRONGENERATOR_H
#define NEUTRONGENERATOR_H

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

class NeutronGenerator : public TObject
{
public:
  enum RunMode_t {
    kMassRapidity,
    kStarlightAscii,
    kFlatMultiplicity,
    kInterface,
  };
  enum HadronicIntModel_t {
    kGlauber,
    kHardSphere
  };

  NeutronGenerator();
  virtual ~NeutronGenerator();
  
  void       	GenerateEvent(const Double_t photonK);
  TClonesArray* ImportParticles(){return fParticles;}
  
  void Initialize();
  void SetBeamParameters(Int_t nuclZ, Int_t nuclA, Double_t gamma){nucleus_Z = nuclZ; nucleus_A = nuclA; beamGamma = gamma;}
  void ReadENDF(Bool_t saveToFile);
  void LoadENDF(const char *filename = "hENDF.root",const char *histname = "hENDF_2D");
  void Setup();
  void Run(const UInt_t nEvents);
  
  void SetRunMode(RunMode_t mode,const char *filename = " ", const char *name1 = " ", const char *name2 = " ");
  void SetRapidityCut(Double_t minY, Double_t maxY){fRapMin = minY; fRapMax = maxY;}
  void SetMassCut(Double_t minM, Double_t maxM){fMassMin = minM; fMassMax = maxM;}
  void SetStoreQA(){kStoreQA = kTRUE;}
  void SetStoreGeneratorFunctions(){kStoreGen = kTRUE;}
  void SetHadronicInteractionModel(HadronicIntModel_t model){fHadronicIntModel = model;}
  
  void 	FinishEvent();
  void 	FinishProduction();
  
  Double_t GetBreakupProbability(const Double_t photonK, const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2);
  Double_t GetTotalFlux(const Double_t photonK);
  
protected:

  void BuildNucleusBreakupProbabilityTable();
  Double_t *NucleusBreakupProbability(const Double_t impactPar);
  Double_t GetBR(const Double_t energyPhoton, const Int_t nNeutrons); 
  void 	CreateNeutrons(const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2);

  void 	InitQAhistograms();
  
  void BuildPhotonFluxModulationTables();
  void BuildHadronicInteractionProbabilityTable();
  Double_t HadronicInteractionProbability(const Double_t impactPar);
  
  Double_t *PhotonFlux(const Double_t energyPhoton);
  Double_t PhotonDensity(const Double_t distance, const Double_t energyPhoton);
  Double_t PhotonDensity2(const Double_t distance, const Double_t energyPhoton);
  TH1D* CreateHist1D(const char* name, const char* title,Int_t nBins, Double_t xMin, Double_t xMax, const char* xLabel, const char* yLabel);
  TH1D* CreateHist1D(const char* name, const char* title,Int_t nBins, Double_t* xBins, const char* xLabel, const char* yLabel);
  TH2D* CreateHist2D(const char* name, const char* title,Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t yMin, Double_t yMax, const char* xLabel, const char* yLabel, const char* zLabel);
  Int_t FromMatrixToVector(Int_t i, Int_t j);
  void  FromVectorToMatrix(Int_t index, Int_t &row, Int_t &col);
  
  void BuildTwoPhotonFluxModulationTable();
  Double_t TwoPhotonFluxModulation(const Double_t VMmass, const Double_t VMrapidity); 

private:

  NeutronGenerator(const NeutronGenerator &o);
  NeutronGenerator &operator=(const NeutronGenerator &o);
  
  RunMode_t		fRunMode;
  HadronicIntModel_t	fHadronicIntModel;
  UInt_t 		iEvent;
  TH1			*hInputRapidity;
  TH1			*hInputMass;
  Double_t		fRapMin;
  Double_t		fRapMax;
  Double_t		fMassMin;
  Double_t		fMassMax;
  ifstream 		fInputStarlightAscii;
  ofstream 		fOutputStarlightAscii;
  TString 		lineString;
  
  const Double_t 	neutron_M = 939.565;
  const Double_t 	hbarc = 0.1973269718;
  const Double_t 	hbarcmev = 197.3269718;
  const Double_t 	pi = 3.14159;
  const Double_t 	alpha = 1/137.035999074;
  const Double_t 	nucleus_R = 6.624;
  const Int_t 		nSteps_impactPar = 120; 
  const Int_t 		nSteps_GG = 20;
  const Int_t 		nSteps_R = 60;  
  const Int_t 		nSteps_Phi = 40;
  const Int_t		nSteps_Energy = 100;
  const Int_t		maxNeutrons = 50;
  Int_t 		nFluxes;
  
  Int_t 	nucleus_Z; 
  Int_t 	nucleus_A;
  Double_t 	beamGamma;
  Double_t 	gammaTarget;
  Double_t 	neutronSepThr;
  Double_t	saturationEnergy;

  Double_t	energyGamma_Xn[625]; 
  Double_t	xSection_Xn[625];
  TGraph	*gSection_Nn;
  TH1D		*hSection_Nn;
  TGraph 	*gPhotonFluxTable;
  TGraph	*gNucleusBreakupTable;
  TH2D		*hTwoPhotonFluxModulationTable;
  TGraph	*gNuclearThickness;
  TGraph 	*gHadronicProbabilityTable;
  TH1D 		*hXsectionXn; 
  TGraph 	*gMeanNeutronN;
  TGraph 	*gWidthNeutronN;
  TF1		*fitMeanNeutronN;
  TF1		*fitWidthNeutronN;
  TH2D		*hENDF_2D;
  TH1D		*hENDF_1D;
  TFile		*fENDFFile; 
  TH2D		*hBranchingRatioMap;
  TH1D		*hEventBreakupMap;
  TGraph	*gUnitaryLeak;
  
  TTree 	*fEventTree;
  TClonesArray	*fParticles;
  TFile 	*fOutputFile;
  
  TList		*fQAhistList;
  TH2D		*hNeutronMultiplicity;
  TH1D		*hEnergyGen;
  TH1D		*hKinEnergyGen;
  TH1D 		*hMomGen;
  TH1D 		*hPhiGen;
  TH1D 		*hThetaGen;
  TH1D 		*hEnergyBoosted;
  TH1D 		*hMomBoosted;
  TH1D 		*hPhiBoosted;
  TH1D 		*hThetaBoosted;
  TH1D 		*hNeutronRapidity;
  TH1D 		*hNeutronEta;
  TH1D 		*hEnergyBin;
  TH1D 		*hEnergyForNeutronMulti;
  TH1D		*hProbabilityXn;
  TH1D		*hPhotonK;
  TH1D		*hRapidityVM;
  TH1D		*hMassVM;
  Bool_t	kStoreQA;
  Bool_t	kStoreGen;
  
  ClassDef(NeutronGenerator,1);
};

#endif
