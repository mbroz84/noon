#ifdef __CLING__
#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "NeutronGenerator.cxx+g"
#endif

void runStandaloneGenerator();
void computeModelBreakups();

void runBreakup(){

  runStandaloneGenerator();
  //computeModelBreakups();

}

void runStandaloneGenerator(){

#if defined(__CINT__)
  gROOT->LoadMacro("NeutronGenerator.cxx+g");
#endif
  
  NeutronGenerator *gen = new NeutronGenerator();
  gen->SetRapidityCut(-4.0,-2.5);
  //gen->SetHadronicInteractionModel(NeutronGenerator::kHardSphere);
  gen->SetStoreQA();
  //gen->SetStoreGeneratorFunctions();
  gen->SetRunMode(NeutronGenerator::kMassRapidity,"ExampleTheory.root","massHist","xSectionHist");
  //gen->SetRunMode(NeutronGenerator::k1n1n);
  gen->Initialize();
  gen->ReadENDF(kFALSE);
  //gen->LoadENDF(); 
  
  gen->Setup();
  gen->Run(1000);

}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  currentLegend->SetNColumns(columns);
  return;
}

void computeModelBreakups(){

#if defined(__CINT__)
  gROOT->LoadMacro("NeutronGenerator.cxx+g");
#endif

  NeutronGenerator *gen = new NeutronGenerator();
  //gen->SetHadronicInteractionModel(NeutronGenerator::kHardSphere);
  gen->SetRunMode(NeutronGenerator::kInterface);
  gen->Initialize();
  gen->ReadENDF(kFALSE);
  //gen->LoadENDF(); 
  gen->Setup();

  TFile *inputFile = new TFile("ExampleTheory.root","READ");
  TH1D *hInputMass = (TH1D*)inputFile->Get("massHist");
  TH1D *hInputRapidity = (TH1D*)inputFile->Get("xSectionHist");

  TString breakups[] = {"All","0n0n","Xn0n","XnXn"};
  TH1D *hRapidityBreakup[4];
  
  Double_t VMrapidity = 0;
  Double_t VMmass = 0;
  Double_t photonK_Low = 0;
  Double_t photonK_High = 0;
  Double_t probLow[4];
  Double_t probHigh[4];
  
  //for(Int_t j=1; j<=41; j++)hInputRapidity->SetBinContent(j,hInputRapidity->GetBinContent(j)*gen->GetTotalFlux(0.5*3.09*TMath::Exp(hInputRapidity->GetBinCenter(j))));
  
  hInputRapidity->SetLineWidth(2);
  hInputRapidity->SetLineColor(kBlack);
  hInputRapidity->SetStats(kFALSE);
  hInputRapidity->SetTitle("");
  hInputRapidity->GetXaxis()->SetTitle("y");
  hInputRapidity->GetYaxis()->SetTitle("d#sigma/dy[mb]");
  hInputRapidity->GetYaxis()->SetTitleOffset(1.5);
  for(Int_t k=0; k<4; k++){
    TString breakupName = breakups[k].Data();
    hRapidityBreakup[k] = (TH1D*)hInputRapidity->Clone(breakupName.Data());
    hRapidityBreakup[k]->SetLineWidth(2);
    hRapidityBreakup[k]->SetLineColor(1+k);
    hRapidityBreakup[k]->SetStats(kFALSE);
    hRapidityBreakup[k]->GetXaxis()->SetTitle("y");
    hRapidityBreakup[k]->GetYaxis()->SetTitle("#sigma [mb]");
    hRapidityBreakup[k]->GetYaxis()->SetTitleOffset(1.5);
    }
  Int_t nBinsInput = hInputRapidity->GetNbinsX()+1;
  for(Int_t j=1; j<=nBinsInput/2; j++){
    VMrapidity = hInputRapidity->GetBinCenter(j);

    for(Int_t k=0; k<4; k++){
      probLow[k] = 0;
      probHigh[k] = 0;
      for(Int_t iEvent = 0; iEvent<1000; iEvent++){
      VMmass = hInputMass->GetRandom();
      photonK_Low = 0.5*VMmass*TMath::Exp(VMrapidity);
      photonK_High = 0.5*VMmass*TMath::Exp(-1*VMrapidity);   
      if(k == 0){
  	probLow[k] += 1.0;
        probHigh[k] += 1.0;
        }
      if(k == 1){
  	probLow[k] += gen->GetBreakupProbability(photonK_Low, 0, 0); 
        probHigh[k] += gen->GetBreakupProbability(photonK_High, 0, 0);
  	}
      if(k == 2){
  	probLow[k] += gen->GetBreakupProbability(photonK_Low, -1, 0); 
        probHigh[k] += gen->GetBreakupProbability(photonK_High, -1, 0);
  	}
      if(k == 3){
  	probLow[k] += gen->GetBreakupProbability(photonK_Low, -1, -1); 
        probHigh[k] += gen->GetBreakupProbability(photonK_High, -1, -1);
  	}
      }
      probLow[k] /= 1000;
      probHigh[k] /= 1000; 
      hRapidityBreakup[k]->SetBinContent(j,hInputRapidity->GetBinContent(j)*probLow[k]+hInputRapidity->GetBinContent(nBinsInput-j)*probHigh[k]);
      hRapidityBreakup[k]->SetBinContent(nBinsInput-j,hRapidityBreakup[k]->GetBinContent(j));
      }
    }
 
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,800);
  TCanvas *c2 = new TCanvas("c2","c2",0,0,800,800);
  TCanvas *c3 = new TCanvas("c3","c3",0,0,800,800);
 
  c1->cd();  
  hRapidityBreakup[0]->GetYaxis()->SetRangeUser(0,700);
  hRapidityBreakup[0]->DrawCopy();
  for(Int_t k=1; k<4; k++) hRapidityBreakup[k]->DrawCopy("same");
  gPad->SetGridy();gPad->SetGridx();
  c2->cd();
  for(Int_t k=1; k<4; k++){ 
    hRapidityBreakup[k]->Divide(hRapidityBreakup[0]);
    hRapidityBreakup[k]->GetYaxis()->SetRangeUser(0.001,0.9);
    if(k == 1)hRapidityBreakup[k]->DrawCopy();
    hRapidityBreakup[k]->DrawCopy("same");
    }
   gPad->SetGridy();gPad->SetGridx();
   c3->cd();
   gPad->SetLogy();
   for(Int_t k=0; k<4; k++){ 
    hRapidityBreakup[k]->GetYaxis()->SetRangeUser(0.01,1);
    if(k == 0)hRapidityBreakup[k]->DrawCopy();
    hRapidityBreakup[k]->DrawCopy("same");
    }
  
  TLegend *myLegend1 = new TLegend(0.42,0.29,0.69,0.48);
  myLegendSetUp(myLegend1,0.04,1);
  myLegend1->AddEntry(hRapidityBreakup[0],"All","l");
  myLegend1->AddEntry(hRapidityBreakup[1],"0n0n","l");
  myLegend1->AddEntry(hRapidityBreakup[2],"0nXn","l");
  myLegend1->AddEntry(hRapidityBreakup[3],"XnXn","l");
  c1->cd(); myLegend1->Draw();
  
  TLatex *noontitle = new TLatex(0.45,0.83,"Hot-spot model + #bf{n_{O}^{O}n}");
  noontitle->SetNDC();
  noontitle->SetTextFont(42);
  noontitle->SetTextSize(0.04);
  //noontitle->Draw();
  
  TLegend *myLegend2 = new TLegend(0.42,0.29,0.69,0.48);
  myLegendSetUp(myLegend2,0.04,1);
  myLegend2->AddEntry(hRapidityBreakup[1],"0n0n/All","l");
  myLegend2->AddEntry(hRapidityBreakup[2],"0nXn/All","l");
  myLegend2->AddEntry(hRapidityBreakup[3],"XnXn/All","l");
  c2->cd(); myLegend2->Draw();
  c3->cd(); myLegend2->Draw();

}

