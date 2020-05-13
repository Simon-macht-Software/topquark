#define MyAnalysis_cxx
// The class definition in MyAnalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("MyAnalysis.C")
// Root > T->Process("MyAnalysis.C","some options")
// Root > T->Process("MyAnalysis.C+")
//

#include "MyAnalysis.h"
#include <iostream>
#include <TH1F.h>
#include <TLatex.h>
#include "MyHists.h"

using namespace std;

void MyAnalysis::BuildEvent() {

  Muons.clear();
  for (int i = 0; i < NMuon; ++i) {
    MyMuon muon(Muon_Px[i], Muon_Py[i], Muon_Pz[i], Muon_E[i]);
    muon.SetIsolation(Muon_Iso[i]);
    muon.SetCharge(Muon_Charge[i]);
    Muons.push_back(muon);
  }
   
  Electrons.clear();
  for (int i = 0; i < NElectron; ++i) {
    MyElectron electron(Electron_Px[i], Electron_Py[i], Electron_Pz[i], Electron_E[i]);
    electron.SetIsolation(Electron_Iso[i]);
    electron.SetCharge(Electron_Charge[i]);
    Electrons.push_back(electron);
  }
   
  Photons.clear();
  for (int i = 0; i < NPhoton; ++i) {
    MyPhoton photon(Photon_Px[i], Photon_Py[i], Photon_Pz[i], Photon_E[i]);
    photon.SetIsolation(Photon_Iso[i]);
    Photons.push_back(photon);
  }
   
  Jets.clear();
  for (int i = 0; i < NJet; ++i) {
    MyJet jet(Jet_Px[i], Jet_Py[i], Jet_Pz[i], Jet_E[i]);
    jet.SetBTagDiscriminator(Jet_btag[i]);
    jet.SetJetID(Jet_ID[i]);
    Jets.push_back(jet);
  }
   
  hadB.SetXYZM(MChadronicBottom_px, MChadronicBottom_py, MChadronicBottom_pz, 4.8);
  lepB.SetXYZM(MCleptonicBottom_px, MCleptonicBottom_py, MCleptonicBottom_pz, 4.8);
  hadWq.SetXYZM(MChadronicWDecayQuark_px, MChadronicWDecayQuark_py, MChadronicWDecayQuark_pz, 0.0);
  hadWqb.SetXYZM(MChadronicWDecayQuarkBar_px, MChadronicWDecayQuarkBar_py, MChadronicWDecayQuarkBar_pz, 0.0);
  lepWl.SetXYZM(MClepton_px, MClepton_py, MClepton_pz, 0.0);
  lepWn.SetXYZM(MCneutrino_px, MCneutrino_py, MCneutrino_pz, 0.0);
  met.SetXYZM(MET_px, MET_py, 0., 0.);
   
  EventWeight *= weight_factor;
}

void MyAnalysis::JEC(TString option){
 Jets.clear();
  for (int i = 0; i < NJet; ++i) {
    MyJet jet(Jet_Px[i], Jet_Py[i], Jet_Pz[i], Jet_E[i]);

    double corr_factor = -1;

    if(option == "nominal"){ corr_factor = 1; }
    else if(option == "up") {corr_factor = 1.05;}
    else if(option == "down") {corr_factor = 0.95;}
    else{ throw std::runtime_error("no JEC version declared");}
   
    TLorentzVector jet_v4_JEC = jet * corr_factor;
    MyJet jet_JEC(jet_v4_JEC.Px(),jet_v4_JEC.Py(),jet_v4_JEC.Pz(),jet_v4_JEC.E());

    jet_JEC.SetBTagDiscriminator(Jet_btag[i]);
    jet_JEC.SetJetID(Jet_ID[i]);
    Jets.push_back(jet_JEC);
  }
}

void MyAnalysis::Begin(TTree * /*tree*/) {
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   TString option = GetOption();
   
}

void MyAnalysis::SlaveBegin(TTree * /*tree*/) {
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave 
   // The tree argument is deprecated (on PROOF 0 is passed).
  
   TString option = GetOption();

   BuildHistmap();
}

//This function fills the defined histogramms for you every time its called
void MyAnalysis::FillHistos(MyHists & h){

  h.h_NMuon->Fill(N_IsoMuon, EventWeight);
  if(N_IsoMuon>0){
    h.h_MuonPt->Fill(muon1->Pt(), EventWeight);
    h.h_MuonEta->Fill(muon1->Eta(), EventWeight);
    h.h_MuonPhi->Fill(muon1->Phi(), EventWeight);
  }

  h.h_MET->Fill(met_E, EventWeight);
 
  h.h_NJets->Fill(NJet, EventWeight);
  if(N_Jets >0){  
    h.h_Jet1_pt->Fill(jet1->Pt(), EventWeight);
    h.h_Jet1_Eta->Fill(jet1->Eta(), EventWeight);
    h.h_Jet1_Phi->Fill(jet1->Phi(), EventWeight);
  }
  if(N_Jets >1){
    h.h_Jet2_pt->Fill(jet2->Pt(), EventWeight);
    h.h_Jet2_Eta->Fill(jet2->Eta(), EventWeight);
    h.h_Jet2_Phi->Fill(jet2->Phi(), EventWeight);
  }

  h.h_NbJets->Fill(N_bJets, EventWeight);
  if(N_bJets >0){
     h.h_bJet1_pt->Fill(b_jet1->Pt(), EventWeight);
     h.h_bJet1_Eta->Fill(b_jet1->Eta(), EventWeight);
     h.h_bJet1_Phi->Fill(b_jet1->Phi(), EventWeight);
  }
  if(N_bJets >1){
     h.h_bJet2_pt->Fill(b_jet2->Pt(), EventWeight);
     h.h_bJet2_Eta->Fill(b_jet2->Eta(), EventWeight);
     h.h_bJet2_Phi->Fill(b_jet2->Phi(), EventWeight);
  }
  
  if(m_top_avg >= 10.) h.h_mtop_rec->Fill(m_top_avg, EventWeight);
}

//This function fills the defined histogramms for you every time its called
vector<TLorentzVector> MyAnalysis::ReconstructNeutrino(){

  if(N_IsoMuon < 1) throw runtime_error("In MyAnalysis::ReconstructNeutrino(): Please require at least one isolated muon before calling this function. Abort.");

  vector<TLorentzVector> solutions;
  
  double mW = 80.;
  
  double dPhi = fabs(muon1->Phi() - met.Phi());
  double mu = pow(mW,2)/2 + muon1->Pt()*met.Pt()*cos(dPhi);

  double A = mu * muon1->Pz() / pow(muon1->Pt(),2);
  double B = ((pow(muon1->E(),2) * pow(met.Pt(),2)) - pow(mu,2)) / (pow(muon1->Pt(),2));

  //argument of sqrt will be: sqrt(A^2 - B). Therefore now check, if A^2 >= B.
  bool isreal = pow(A,2) >= B;

  //if the sqrt would be imaginary, instead just take the real part of the solution
  if(!isreal || pow(A,2) - B == 0){
    double pz = A;
    TLorentzVector sol;
    sol.SetXYZM(met.Px(), met.Py(), pz, 0.);
    solutions.emplace_back(sol);
    return solutions;
  }
  else{
    TLorentzVector sol;
    
    double pz = A+sqrt(pow(A,2) - B);
    sol.SetXYZM(met.Px(), met.Py(), pz, 0.);
    solutions.emplace_back(sol);
    
    pz = A-sqrt(pow(A,2) - B);
    sol.SetXYZM(met.Px(), met.Py(), pz, 0.);
    solutions.emplace_back(sol);
    return solutions;
  }
  			      
}

vector<pair<TLorentzVector,TLorentzVector>> MyAnalysis::ReconstructTTbar(int min_njets){
  if(N_IsoMuon < 1) throw runtime_error("In MyAnalysis::ReconstructTTbar(): Please require at least one isolated muon before calling this function. Abort.");
  if(N_Jets < 3) throw runtime_error("In MyAnalysis::ReconstructTTbar(): Please require at least three jets before calling this function. Abort.");
  if(min_njets < 2 || min_njets > 4) throw runtime_error("In MyAnalysis::ReconstructTTbar(): Please choose a minimum number of jets of 2, 3, or 4. Abort.");
  if(min_njets > N_Jets - 1) throw runtime_error("In MyAnalysis::ReconstructTTbar(): The number of jets in this event is not sufficient to reconstruct the ttbar system with the required number of jets in the hadronic top quark candidate. Please require a higher number of jets before running this function or reduce the minimum number of jets that shall be used for the hadronic top quark hypothesis.");
  
  vector<pair<TLorentzVector,TLorentzVector>> ttbarhypos;
  
  vector<bool> jetused;
  for(unsigned int i=0; i<Jets.size(); i++) jetused.emplace_back(false);
  
  
  // 0) Define needed variables
  TLorentzVector muon;
  muon.SetXYZM(muon1->Px(), muon1->Py(), muon1->Pz(), muon1->M());

  // 1) Get the reconstructed neutrino (1 or 2 solutions)
  vector<TLorentzVector> neutrinosols = ReconstructNeutrino();

  // 2) loop over the solutions and repeat the following for each solution.
  for(unsigned int i=0; i<neutrinosols.size(); i++){
    TLorentzVector neutrino = neutrinosols[i];

    // 3) Combine (pT-leading muon) and neutrino to W_lep
    TLorentzVector W_lep = neutrino + muon;

    // 4) assign one jet to the leptonic top quark. Loop over all possible permutations.
    for(unsigned int i=0; i<Jets.size(); i++){
      TLorentzVector jet_lep;
      jet_lep.SetXYZM(Jets[i].Px(), Jets[i].Py(), Jets[i].Pz(), Jets[i].M());

      TLorentzVector t_lep = W_lep + jet_lep;
      jetused[i] = true;

      // 5) build all possible permutations of at least 2, at most 4 jets to build hadronic top quark candidates
      
      // 5.1) first jet loop
      for(unsigned int j=0; j<Jets.size(); j++){
	if(jetused[j]) continue;

	TLorentzVector myjet1(Jets[j].Px(), Jets[j].Py(), Jets[j].Pz(), Jets[j].E());
	//jetused[j] = true;

	//5.2) second jet loop
	for(unsigned k=0; k<Jets.size(); k++){
	  if(jetused[k]) continue;
	  if(k<=j) continue;

	  TLorentzVector myjet2(Jets[k].Px(), Jets[k].Py(), Jets[k].Pz(), Jets[k].E());
	  //jetused[k] = true;

	  if(N_bJets > 0){
	    if(!(Jets[i].IsBTagged()) && !(Jets[j].IsBTagged() || Jets[k].IsBTagged())) continue;

	    if(N_bJets > 1){
	      if(!(Jets[i].IsBTagged()) || !(Jets[j].IsBTagged() || Jets[k].IsBTagged())) continue;
	    }
	  }
	  
	  // here we can already fill hypotheses with only 2 jets.
	  pair<TLorentzVector, TLorentzVector> hypo2jets = make_pair(t_lep, myjet1+myjet2);
	  if(min_njets <= 2) ttbarhypos.emplace_back(hypo2jets);
	  
	  // 5.3) third jet loop
	  if(Jets.size() < 4) continue;
	  for(unsigned int l=0; l<Jets.size(); l++){
	    if(jetused[l]) continue;
	    if(l<=k) continue;
	    TLorentzVector myjet3(Jets[l].Px(), Jets[l].Py(), Jets[l].Pz(), Jets[l].E());
	    //jetused[l] = true;

	    if(N_bJets > 0){
	      if(!(Jets[i].IsBTagged()) && !(Jets[j].IsBTagged() || Jets[k].IsBTagged() || Jets[l].IsBTagged())) continue;

	      if(N_bJets > 1){
		if(!(Jets[i].IsBTagged()) || !(Jets[j].IsBTagged() || Jets[k].IsBTagged() || Jets[l].IsBTagged())) continue;
	      }
	    }
	    
	    // here we fill hypotheses with exactly 3 jets now
	    pair<TLorentzVector, TLorentzVector> hypo3jets = make_pair(t_lep, myjet1+myjet2+myjet3);
	    if(min_njets <= 3) ttbarhypos.emplace_back(hypo3jets);


	    // 5.4) fourth and final jet loop
	    if(Jets.size() < 5) continue;
	    for(unsigned int m=0; m<Jets.size(); m++){
	      if(jetused[m]) continue;
	      if(m<=l) continue;
	      TLorentzVector myjet4(Jets[m].Px(), Jets[m].Py(), Jets[m].Pz(), Jets[m].E());

	      if(N_bJets > 0){
		if(!(Jets[i].IsBTagged()) && !(Jets[j].IsBTagged() || Jets[k].IsBTagged() || Jets[l].IsBTagged() || Jets[m].IsBTagged())) continue;

		if(N_bJets > 1){
		  if(!(Jets[i].IsBTagged()) || !(Jets[j].IsBTagged() || Jets[k].IsBTagged() || Jets[l].IsBTagged()|| Jets[m].IsBTagged())) continue;
		}
	      }
	      
	      // fill hypotheses with exactly 4 jets
	      pair<TLorentzVector, TLorentzVector> hypo4jets = make_pair(t_lep, myjet1+myjet2+myjet3+myjet4);
	      if(min_njets <= 4) ttbarhypos.emplace_back(hypo4jets);
	      
	    } // End of fourth loop over jets for the hadronic top
	  } // End of third loop over jets for the hadronic top
	} // End of second loop over jets for the hadronic top
      } // End of first loop over jets for the hadronic top

      // Reset 'jetused'
      jetused[i] = false;
      
    } // End of loop over jets being assigned to leptonic top
  } // End of loop over neutrino solutions
  
  return ttbarhypos;
}

pair<TLorentzVector, TLorentzVector> MyAnalysis::SelectBestTTbarHypothesis(vector<pair<TLorentzVector, TLorentzVector>> hypos, double max_dM){

  //Make sure that all entries consist of one leptonic and one hadronic top quark!

  //Loop over all hypotheses and search for the one with the smallest mass difference between the leptonic and the hadronic cadidate
  double dM_min = 99999;
 
  pair<TLorentzVector, TLorentzVector> besthypo;
  for(unsigned int i=0; i<hypos.size(); i++){
    double dM = fabs(hypos[i].first.M() - hypos[i].second.M());
    if(dM < dM_min){
      dM_min = dM;
      besthypo = hypos[i];
    }   
  }

  if(dM_min <= max_dM) return besthypo;
  else{
    TLorentzVector dummy(0.,0.,0.,0.);
    return make_pair(dummy,dummy);
  }
}


void MyAnalysis::PrintModule(Long64_t entry) {

//   ++++++++++++++++++++++++++++++++++++++++
//   Print Module
//   ++++++++++++++++++++++++++++++++++++++++
     cout<< "New Event: "<<endl;
     cout << "Jets: " << endl;
     for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
     cout << "pt, eta, phi, btag, id: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", " << it->IsBTagged() << ", " << it->GetJetID()
     << endl;
     }
     cout << "Muons: " << endl;
     for (vector<MyMuon>::iterator it = Muons.begin(); it != Muons.end(); ++it) {
     cout << "pt, eta, phi, iso, charge: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", "
     << it->GetIsolation() << ", " << it->GetCharge() << endl;
     }
     cout<<""<<endl;
}



Bool_t MyAnalysis::Process(Long64_t entry) {
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either MyAnalysis::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   ++TotalEvents;
   
   GetEntry(entry);
   
   if (TotalEvents % 10000 == 0)
     cout << "Next event -----> " << TotalEvents << endl;
   
   
   BuildEvent();

   //Jet Energy Correction: Change four-vector 
   JEC("nominal"); //nominal, down, up
    
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //
   //   Loop over different objects. Save number of objects and leading, second leading objects for plotting
   //  
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   double MuonRelIsoCut = 0.10;
   N_IsoMuon = 0;
   N_Jets = 0;
   N_bJets = 0;
   m_top_avg = -1.;

   // Loop over all muons. The number of isolated muons 'N_IsoMuon' is summed up. 
   // The leading muon and the second leading muon get defined by pointing to the 
   // respective iterator (muon object).
   for (vector<MyMuon>::iterator jt = Muons.begin(); jt != Muons.end(); ++jt) {
      if (jt->IsIsolated(MuonRelIsoCut)) {
         ++N_IsoMuon;
         if (N_IsoMuon == 1) muon1 = &(*jt);
         if (N_IsoMuon == 2) muon2 = &(*jt);
      }
   }

   //Loop over all jets. The number of jets N_Jets is summed up. 
   // The leading jet and the second leading jet get defined by pointing to the 
   // respective iterator (jet object).
   for(vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it){
      ++N_Jets;
     
      if (N_Jets == 1) jet1 = &(*it);
      if (N_Jets == 2) jet2 = &(*it);
          
      //The number of b-tagged jets N_bJets is summed up. 
      // The leading b-tagged jet and the second leading b-tagged jet get defined by pointing to the 
      // respective iterator (jet object).

      if(it->IsBTagged()){
	++N_bJets;
	if (N_bJets == 1) b_jet1 = &(*it);
	if (N_bJets == 2) b_jet2 = &(*it);
      }
   }

   met_E  = met.E();
   
   ///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   // +++++++++++++++++++++++++++++++++++++++++++
   // Start of Exercises: 
   // Code you have to edit
   // +++++++++++++++++++++++++++++++++++++++++++






   
   // +++++++++++++++++++++++++++++++++++++++++++
   // Exercise 1: MC / data comparisons
   // +++++++++++++++++++++++++++++++++++++++++++
   
   //Define Cut values
   //  double MuonPtCut = 30.;
	EreignisseN+=EventWeight;
   //Require at least one isolated muon
   if (N_IsoMuon >= 1 && triggerIsoMu24) {
	
	if (N_Jets >= 1 && N_Jets < 9){

		if (muon1 ->Pt()>20&&muon1 ->Pt()<240){
			if (N_IsoMuon < 3) {
				// if (jet1->Eta() < 2 && jet1->Eta()>=-2){
					if (met_Pt < 320){
	
			FillHistos(new_hist);
					 }
				// }
			}
			
		}
	}

     //This function fills for you a set of histogramms every time it is called. Don't forget to define it in MyAnalysis.h
     FillHistos(hists_nocuts);

     
     
     // Have a look at your histograms and compare the different background samples.
     // Try to enrich the fraction of ttbar by cutting on any of the distributions 
     // Plot all variables after every cut you introduce. Therefore define a new set of Histogramms in Analysis.h  e.g. hist_muoncut     
     
     /*
     if(muon1->Pt()>MuonPtCut){
   
     }
     */








     

     // +++++++++++++++++++++++++++++++++++++++++++
     // Exercise 2: Measurement of the cross section
     // +++++++++++++++++++++++++++++++++++++++++++
     // Once you have optimized your event selection, define a float (in Analysis.h) for the weighted number of selected events. 
     // Also, define a float for the weighted number of generated events. With these numbers, 
     // the selection efficiency can be calculated. Afterwards, move to the 'example.C' 
     // file to determine the efficiency and the ttbar cross section.










     
     // ++++++++++++++++++++++++++++++++++++++++++++
     // Exercise 3: Reconstruction of the top quark mass
     // ++++++++++++++++++++++++++++++++++++++++++++

     // Start off with a new event selection here (still inside the loop requiring at least one isolated muon and the trigger, but nothing from exercise 2).


     
     //The 3 lines below are working - un-comment them! Just make sure to apply your new event selection before ;)
     /*
       vector<pair<TLorentzVector,TLorentzVector>> ttbarreco = ReconstructTTbar(2);
       pair<TLorentzVector, TLorentzVector> tt_best = SelectBestTTbarHypothesis(ttbarreco);
       m_top_avg = (tt_best.first.M() + tt_best.second.M())/2;
     */


     
     // You should implement a new set of histograms to be filled here, otherwise you can't analyze anything... 












     // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // +++++++++++++++++++++++++++++++++++++++++
     // Do not edit the code below
     // +++++++++++++++++++++++++++++++++++++++++     
   }
   return kTRUE;
}




void MyAnalysis::SlaveTerminate() {
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   
}

void MyAnalysis::Terminate() {
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
}
