#include "MyAnalysis.h"
#include "Plotter.h"
#include <iostream>
#include <TChain.h>
#include <TF1.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <string>
#include <memory>

int main() {
   
   // Integrated luminosity (this value is needed for the calculation of the cross section)
   float lumi = 50.;

   vector<pair<TString, TString>> processes (9);
   
   processes[0] = make_pair("Data", "data");
   processes[1] = make_pair("TTbar", "ttbar");
   processes[2] = make_pair("Wjets", "wjets");
   processes[3] = make_pair("DY", "dy");
   processes[4] = make_pair("WW", "ww");
   processes[5] = make_pair("WZ", "wz");
   processes[6] = make_pair("ZZ", "zz");
   processes[7] = make_pair("QCD", "qcd");
   processes[8] = make_pair("single Top", "single_top");

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //
   // Set up plotters
   //
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   std::vector<Plotter> Plotters, Plotters_MC;  
   MyAnalysis* dummy = new MyAnalysis();
   dummy->BuildHistmap();
   for(std::map<TString,MyHists*>::iterator it=dummy->GetHistmap()->begin(); it!=dummy->GetHistmap()->end(); ++it){
     std::cout << "Setting up the following set of histograms: " << it->first <<  endl;
     TString outname = it->first;
     Plotter p, p_mc;
     p.SetOutname(outname);
     p_mc.SetOutname(outname+"_MC");
     Plotters.emplace_back(p);
     Plotters_MC.emplace_back(p_mc);
   }
   delete dummy;
   
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //
   // Iterate over elements of map to chain all processes
   //
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //Here you can access variables, defined in MyAnalysis 
   // Tip: use a vector to save your defined variables for all processes

   // Example: 
   // vector<float> 'vec_var' ;
   // ...


   for(unsigned int i = 0; i < processes.size(); i++){
   
     MyAnalysis* A = new MyAnalysis();
     std::unique_ptr<TChain> ch;
     ch.reset(new TChain("events"));
     TString filename = "files/" + processes[i].second + ".root";
     ch->Add(filename);
     cout << endl << "Now starting to process the sample '" << processes[i].first << "'" << endl; 
     ch->Process(A);

	cout << "gewichtete Ereignisse " << A -> EreignisseN << endl;
     // Example: 
     // 'vec_var' = A -> 'var';
     // ...

     // Add this process to the plotter
     if(processes[i].first == "Data"){
       int idx = 0;
       for(std::map<TString,MyHists*>::iterator it2=A->GetHistmap()->begin(); it2!=A->GetHistmap()->end(); ++it2){
	 Plotters[idx].SetData(A->GetHists(it2->first)->get_histvec(), string(processes[i].first));
	 idx++;
       }
     }
     else{
       int idx = 0;
       for(std::map<TString,MyHists*>::iterator it2=A->GetHistmap()->begin(); it2!=A->GetHistmap()->end(); ++it2){
	 Plotters[idx].AddBg(A->GetHists(it2->first)->get_histvec(), std::string(processes[i].first));
	 Plotters_MC[idx].AddBg(A->GetHists(it2->first)->get_histvec(), std::string(processes[i].first));
	 idx++;
       }
     }

     // In any case: save output histograms
     TString outname = "output_" + processes[i].second + ".root";
     TFile* outfile = new TFile(outname, "RECREATE");
     for(std::map<TString,MyHists*>::iterator it2=A->GetHistmap()->begin(); it2!=A->GetHistmap()->end(); ++it2){
       outfile->mkdir(it2->first);
       outfile->cd(it2->first);
       for(unsigned int j=0; j<A->GetHists(it2->first)->get_histvec().size(); j++) A->GetHists(it2->first)->get_histvec()[j]->Write();
     }
     outfile->Close();
     delete outfile;
     delete A;
   }
   
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // Plot everything
   for(unsigned int i=0; i<Plotters.size(); i++){
     Plotters[i].Plot();
     Plotters_MC[i].Plot();
   }
   cout << endl << endl << "Done processing and plotting!" << endl;


   //++++++++++++++++++++++++++++++++++++++++++++++++++++
   // Exercise 2: Measurement of the cross section
   //++++++++++++++++++++++++++++++++++++++++++++++++++++

   // Save number of events for every process in a new vector<float> before and after your cuts. 
   // Get the efficiency of your event selection on ttbar. Calculate the number of (weighted) selected events divided
   // by the number of (weighted) generated events.
   // double eff = ...;

   // Subtract the non-ttbar background from data
   // ...

   // Calculate the cross section
   // ...

   // Propagate a statistic uncertainty on the cross section
   // ...

   //investigate a systematic uncertanty introduced by JEC on the cross section
   // ...
   
   //Print/Save your results
   // ...




   //++++++++++++++++++++++++++++++++++++++++++++++++++++
   // Exercise 3: Reconstruction of the top quark mass
   //++++++++++++++++++++++++++++++++++++++++++++++++++++

   // Set the bool to 'true' in order to fit the reconstructed top quark mass and obtain a result! Only modify the four variables given
   
   bool run_this_part = false;
   TString foldername = "topreco";
   double fit_min = 130.;
   double fit_max = 210.;





   //++++++++++++++++++++++++++++++++++++++++++++++++++++
   // Do not edit the part below
   
   if(run_this_part){
     // Find the data histogram:
     TString filename_base = "output";
     TString filename_data = filename_base + "_data.root";
     TFile* infile_data = new TFile(filename_data, "READ");
     TString histname = foldername + "/mtop_rec";
     TH1D* h_mtop_data = (TH1D*)infile_data->Get(histname);

     for(unsigned int i = 0; i < processes.size(); i++){
     
       // Open all root files containing the histograms
       TString myfilename = filename_base + "_" + processes[i].second + ".root";
       TFile* infile = new TFile(myfilename, "READ");

       // Open the relevant histogram
       TH1D* myhist = (TH1D*)infile->Get(histname);

       //Subtract non-ttbar from data
       if(processes[i].first != "Data" && processes[i].first != "TTbar") h_mtop_data->Add(myhist, -1);
       for(int j=1; j<h_mtop_data->GetNbinsX()+1; j++){
	 if(h_mtop_data->GetBinContent(j) < 0.){
	   h_mtop_data->SetBinContent(j,0.);
	   h_mtop_data->SetBinError(j,0.);
	   h_mtop_data->SetMinimum(.0);
	 }
       }
     
       delete myhist;
       delete infile;
     }

     TCanvas* c = new TCanvas();
     h_mtop_data->Draw();

     // Fit The distribution to extract the reconstructed top quark mass
     //TF1* fitfunc = new TF1("myfit", "[0]", 150, 190);
     h_mtop_data->Fit("gaus", "", "", fit_min, fit_max);
     TF1* fit = h_mtop_data->GetFunction("gaus");
     double mean = fit->GetParameter(1);
     double unc = fit->GetParError(1);

     cout << endl << endl << endl << "----------------------------------------------------------" << endl;
     cout << "Fitted top quark mass: " << mean << " +- " << unc << " GeV" << endl << endl;


   
     c->SaveAs("ReconstructedTopMass.pdf");
     delete c;
     delete h_mtop_data;
     delete infile_data;
   }
}





