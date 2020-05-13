#pragma once

#include <TH1F.h>

class MyHists{
public:
  MyHists();
  virtual ~MyHists();

  //std::unique_ptr<TH1F> h_Mmumu, h_NMuon;
  TH1F *h_NMuon, *h_MuonPt, *h_MuonEta, *h_MuonPhi, *h_NJets, *h_Jet1_pt, *h_Jet1_Eta, *h_Jet1_Phi, *h_Jet2_pt, *h_Jet2_Eta, *h_Jet2_Phi, *h_NbJets, *h_bJet1_pt, *h_bJet1_Eta, *h_bJet1_Phi, *h_bJet2_pt, *h_bJet2_Eta, *h_bJet2_Phi, *h_MET, *h_mtop_rec;

  std::vector<TH1F*> get_histvec();

};
