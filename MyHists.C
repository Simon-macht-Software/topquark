#include <math.h>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>
#include "MyHists.h" 
#include <TError.h> 

using namespace std;

MyHists::MyHists(){

  gErrorIgnoreLevel = 2002;
  // Set up all the histograms
  h_NMuon = new TH1F("NMuon", "Number of muons", 5, 0, 5);
  h_NMuon->SetXTitle("Number of muons");
  h_NMuon->Sumw2();

  h_MuonPt = new TH1F("MuonPt", "Muon p_{T}", 40, 0, 400);
  h_MuonPt -> SetXTitle("Muon p_{T} [GeV]");
  h_MuonPt -> Sumw2();  

  h_MuonEta = new TH1F("MuonEta", "Muon Pseudorapidity", 20, -3.5, 3.5);
  h_MuonEta -> SetXTitle("Muon #eta");
  h_MuonEta -> Sumw2(); 

  h_MuonPhi = new TH1F("MuonPhi", "Muon #phi", 20, -3.5, 3.5);
  h_MuonPhi -> SetXTitle("Muon #phi");
  h_MuonPhi -> Sumw2(); 

  h_MET = new TH1F("MET", "MET", 40, 0, 400);
  h_MET -> SetXTitle("Missing transverse energy [GeV]");
  h_MET -> Sumw2();

  h_NJets = new TH1F("NJets", "Number of jets", 12, 0, 12);
  h_NJets -> SetXTitle("Number of jets");
  h_NJets -> Sumw2(); 

  h_Jet1_pt = new TH1F("Jet1_Pt", "leading jet pT", 60, 0, 600);
  h_Jet1_pt -> SetXTitle("Leading jet p_{T} [GeV]");
  h_Jet1_pt -> Sumw2();

  h_Jet1_Eta = new TH1F("Jet1_Eta", "leading jet #eta", 20, -3.5, 3.5);
  h_Jet1_Eta -> SetXTitle("Leading jet #eta");
  h_Jet1_Eta -> Sumw2();

  h_Jet1_Phi = new TH1F("Jet1_Phi", "leading jet #phi", 20, -3.5, 3.5);
  h_Jet1_Phi -> SetXTitle("Leading jet #phi");
  h_Jet1_Phi -> Sumw2();

  h_Jet2_pt = new TH1F("Jet2_Pt", "second jet pT", 45, 0, 450);
  h_Jet2_pt -> SetXTitle("Second jet p_{T} [GeV]");
  h_Jet2_pt -> Sumw2();

  h_Jet2_Eta = new TH1F("Jet2_Eta", "second jet #eta", 20, -3.5, 3.5);
  h_Jet2_Eta -> SetXTitle("Second jet #eta");
  h_Jet2_Eta -> Sumw2();

  h_Jet2_Phi = new TH1F("Jet2_Phi", "second jet #phi", 20, -3.5, 3.5);
  h_Jet2_Phi -> SetXTitle("Second jet #phi");
  h_Jet2_Phi -> Sumw2();

  h_NbJets = new TH1F("NbJets", "Number of b-tagged jets", 5, 0,5);
  h_NbJets -> SetXTitle("Number of b-tagged jets");
  h_NbJets -> Sumw2(); 

  h_bJet1_pt = new TH1F("bJet1_Pt", "leading b-tagged jet pT", 40, 0, 400);
  h_bJet1_pt -> SetXTitle("Leading b-tagged jet p_{T} [GeV]");
  h_bJet1_pt -> Sumw2();

  h_bJet1_Eta = new TH1F("bJet1_Eta", "leading b-tagged jet #eta", 20, -3.5, 3.5);
  h_bJet1_Eta -> SetXTitle("Leading b-tagged jet #eta");
  h_bJet1_Eta -> Sumw2();

  h_bJet1_Phi = new TH1F("bJet1_Phi", "leading b-tagged jet #phi", 20, -3.5, 3.5);
  h_bJet1_Phi -> SetXTitle("Leading b-tagged jet #phi");
  h_bJet1_Phi -> Sumw2();

  h_bJet2_pt = new TH1F("bJet2_Pt", "second b-tagged jet pT", 25, 0, 250);
  h_bJet2_pt -> SetXTitle("Second b-tagged jet p_{T} [GeV]");
  h_bJet2_pt -> Sumw2();

  h_bJet2_Eta = new TH1F("bJet2_Eta", "second b-tagged jet #eta", 20, -3.5, 3.5);
  h_bJet2_Eta -> SetXTitle("Second b-tagged jet #eta");
  h_bJet2_Eta -> Sumw2();

  h_bJet2_Phi = new TH1F("bJet2_Phi", "second b-tagged jet #phi", 20, -3.5, 3.5);
  h_bJet2_Phi -> SetXTitle("Second b-tagged jet #phi");
  h_bJet2_Phi -> Sumw2();
  
  h_mtop_rec = new TH1F("mtop_rec", "reconstructed top quark mass", 35, 10, 500);
  h_mtop_rec -> SetXTitle("Reconstructed top quark mass [GeV]");
  h_mtop_rec -> Sumw2();
  
}

MyHists::~MyHists(){
  delete h_NMuon;
  delete h_MuonPt;
  delete h_MuonEta;
  delete h_MuonPhi;
  delete h_MET;
  delete h_NJets;
  delete h_Jet1_pt;
  delete h_Jet1_Eta;
  delete h_Jet1_Phi;
  delete h_Jet2_pt;
  delete h_Jet2_Eta;
  delete h_Jet2_Phi;
  delete h_NbJets;
  delete h_bJet1_pt;
  delete h_bJet1_Eta;
  delete h_bJet1_Phi;
  delete h_bJet2_pt;
  delete h_bJet2_Eta;
  delete h_bJet2_Phi;
  delete h_mtop_rec;
} 

std::vector<TH1F*> MyHists::get_histvec(){
  gErrorIgnoreLevel = 2002;
  std::vector<TH1F*> vec;
  
  TH1F* NMuon = (TH1F*)h_NMuon->Clone();
  vec.emplace_back(NMuon);
  
  TH1F* MuonPt = (TH1F*)h_MuonPt->Clone();
  vec.emplace_back(MuonPt);

  TH1F* MuonEta = (TH1F*)h_MuonEta->Clone();
  vec.emplace_back(MuonEta);

  TH1F* MuonPhi = (TH1F*)h_MuonPhi->Clone();
  vec.emplace_back(MuonPhi);

  TH1F* MET = (TH1F*)h_MET->Clone();
  vec.emplace_back(MET);

  TH1F* NJets = (TH1F*)h_NJets->Clone();
  vec.emplace_back(NJets);

  TH1F* Jet1_Pt = (TH1F*)h_Jet1_pt->Clone();
  vec.emplace_back(Jet1_Pt);

  TH1F* Jet1_Eta = (TH1F*)h_Jet1_Eta->Clone();
  vec.emplace_back(Jet1_Eta);

  TH1F* Jet1_Phi = (TH1F*)h_Jet1_Phi->Clone();
  vec.emplace_back(Jet1_Phi);  

  TH1F* Jet2_Pt = (TH1F*)h_Jet2_pt->Clone();
  vec.emplace_back(Jet2_Pt);  

  TH1F* Jet2_Eta = (TH1F*)h_Jet2_Eta->Clone();
  vec.emplace_back(Jet2_Eta);

  TH1F* Jet2_Phi = (TH1F*)h_Jet2_Phi->Clone();
  vec.emplace_back(Jet2_Phi);

  TH1F* NbJets = (TH1F*)h_NbJets->Clone();
  vec.emplace_back(NbJets);

  TH1F* bJet1_Pt = (TH1F*)h_bJet1_pt->Clone();
  vec.emplace_back(bJet1_Pt);
  
  TH1F* bJet1_Eta = (TH1F*)h_bJet1_Eta->Clone();
  vec.emplace_back(bJet1_Eta);

  TH1F* bJet1_Phi = (TH1F*)h_bJet1_Phi->Clone();
  vec.emplace_back(bJet1_Phi);

  TH1F* bJet2_Pt = (TH1F*)h_bJet2_pt->Clone();
  vec.emplace_back(bJet2_Pt);  

  TH1F* bJet2_Eta = (TH1F*)h_bJet2_Eta->Clone();
  vec.emplace_back(bJet2_Eta);

  TH1F* bJet2_Phi = (TH1F*)h_bJet2_Phi->Clone();
  vec.emplace_back(bJet2_Phi);
  
  TH1F* mtop_rec = (TH1F*)h_mtop_rec->Clone();
  vec.emplace_back(mtop_rec);
  
  return vec;
}
