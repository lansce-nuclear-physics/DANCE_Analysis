#include <iostream>
#include <sstream>
#include <fstream>

void hCrID_Diagnostics() {
  
  //Location of summed rootfiles
  stringstream fpath;
  fpath.str();
  fpath << "/home/cprokop/CJP/DANCE_Analysis/stage0_root/";
  
  int start_run = 69275;
   int end_run = 71297;
   
  // int start_run = 63015;
 //  int end_run = 65265;

  // int start_run = 67339;
 //  int end_run = 68900;


  int counter_10=0;
  int counter_89=0;
  
  const int N_Runs = end_run-start_run+1;
  
  cout<<N_Runs<<endl;

  TFile *fin[N_Runs];
  TH1D *hID;
  
  double Counts[N_Runs][162];
  int Run_Num[N_Runs];

  int counter=0;
  
  TH2D *hCounts = new TH2D("hCounts","hCounts",N_Runs,start_run,end_run,162,0,162);

    ifstream testrun;
    bool fileexists=false;



  for(int eye=start_run; eye<=end_run; eye++) {

    stringstream teststring;
    teststring.str();
    teststring << fpath.str() << "Stage0_Histograms_Run_"<< eye<<"_500ns_CW_0ns_CBT_0ns_DEBT.root";
    cout<<eye<<"  "<<teststring.str()<<endl;
    
    testrun.open(teststring.str().c_str());
    fileexists=false;
    
    if(testrun.is_open()) {
      fileexists=true;
      cout<<"Run: "<<eye<<" exists"<<endl;
      testrun.close();
    }
    
    
    if(fileexists) {  
      
      fin[counter] = new TFile(Form("%s/Stage0_Histograms_Run_%d_500ns_CW_0ns_CBT_0ns_DEBT.root",fpath.str().c_str(),eye));  
      Run_Num[counter] = eye;
      hID = (TH1D*)fin[counter]->FindObjectAny("hID");
      
      cout<<hID<<endl;
      
      for(int jay=0; jay<162; jay++) {
	Counts[counter][jay] = hID->GetBinContent(jay+1);
	hCounts->Fill(eye,jay,Counts[counter][jay]);
      }
      // hID[counter]->Delete();
      fin[counter]->Close();
    }
    else {
      for(int jay=0; jay<162; jay++) {
	Counts[counter][jay] = -1;
      }
    }

    if(Counts[counter][10]>0) {
      counter_10++;
    }
    if(Counts[counter][89]>0) {
      counter_89++;
    }
    
    counter++;
  }
  
  hCounts->Draw("colz");

  cout<<"Number of runs with ch 10: "<<counter_10<<endl;
  cout<<"Number of runs with ch 89: "<<counter_89<<endl;
  
}

/*
  //Summed up root files from FARE
  TFile *fin[5];
  fin[0] = new TFile(Form("%s/Sum_Au197_2mm_60452_60510_and_60534_60571.root",fpath.str().c_str()));
  fin[1] = new TFile(Form("%s/Sum_Au197_4mm_60410_60451_and_60572_60584.root",fpath.str().c_str()));
  fin[2] = new TFile(Form("%s/Sum_Au197_7mm_60511_60533.root",fpath.str().c_str()));
  fin[3] = new TFile(Form("%s/Sum_Au197_10mm_60384_60408.root",fpath.str().c_str()));
  fin[4] = new TFile(Form("%s/Sum_Au197_15mm_60375_60381.root",fpath.str().c_str()));
  
  //Histograms
  TH3D *h3D[5];
  TH2D *h3D_yxproj[5];
  TH3D *hQgated[5][3];
  TH1D *hQgated_xproj[5][3];
  TH1D *h1D_bkgsub[5];

  //Canvasas
  TCanvas *cQGated[5];
  TCanvas *cBeamMon[5];
  TCanvas *c2D[5];
  
  //Beam Monitors
  TH1D *h_ph_U235[5];
  TH1D *h_ph_Li6[5];
  TH1D *h_ph_He3[5];
  
  //Multiplicities
  int mcl_low = 3;
  int mcl_high = 4;
  

  
  for(int eye=0; eye<5; eye++) {
    // if(fin[eye]) {
    //  if(eye==4 || eye ==3) {
      fin[eye]->cd();

      TDirectory *histodir = (TDirectory*)fin[eye]->Get("Standard/AJC");
      TDirectory *stddir = (TDirectory*)fin[eye]->Get("Standard");
      TDirectory *bmdir = (TDirectory*)fin[eye]->Get("Standard/BM_dir");
      cout<<"TDirectories"<<endl;

      cQGated[eye] = new TCanvas(Form("cQGated_%fmm",Diameter[eye]),Form("cQGated_%fmm",Diameter[eye]));
      //  c2D[eye] = new TCanvas(Form("c2D_%fmm",Diameter[eye]),Form("c2D_%fmm",Diameter[eye]));
      cBeamMon[eye] = new TCanvas(Form("cBeamMon_%fmm",Diameter[eye]),Form("cBeamMon_%fmm",Diameter[eye]));
      cBeamMon[eye]->Divide(2,2);
      h_ph_U235[eye] = ((TH1D*)bmdir->FindObjectAny("h_ph_U235"));
      cBeamMon[eye]->cd(1);
      h_ph_U235[eye]->Draw();
      U235_Counts[eye] = h_ph_U235[eye]->Integral(0,5000);
      U235_Counts_Error[eye] = TMath::Power(U235_Counts[eye],0.5);

      h_ph_Li6[eye] = ((TH1D*)bmdir->FindObjectAny("h_ph_Li6"));
      cBeamMon[eye]->cd(2);
      h_ph_Li6[eye]->Draw();
      Li6_Counts[eye] = h_ph_Li6[eye]->Integral(0,5000);
      Li6_Counts_Error[eye] = TMath::Power(Li6_Counts[eye],0.5);

      h_ph_He3[eye] = ((TH1D*)bmdir->FindObjectAny("h_ph_He3"));
      cBeamMon[eye]->cd(3);
      h_ph_He3[eye]->Draw();
      He3_Counts[eye] = h_ph_He3[eye]->Integral(0,5000);
      U235_Counts_Error[eye] = TMath::Power(U235_Counts[eye],0.5);

      cout<<eye<<"  "<<h_ph_U235[eye]<<"  "<<U235_Counts[eye]<<endl;
      cout<<eye<<"  "<<h_ph_Li6[eye]<<"  "<<Li6_Counts[eye]<<endl;
      cout<<eye<<"  "<<h_ph_He3[eye]<<"  "<<He3_Counts[eye]<<endl;

      cBeamMon[eye]->cd(4);
      //   c2D[eye]->cd();
      h3D[eye] = ((TH3D*)stddir->FindObjectAny("En_Etot_Mcl"));
      cout<<h3D[eye]<<endl;
      h3D[eye]->GetZaxis()->SetRangeUser(mcl_low,mcl_high);
      h3D[eye]->GetYaxis()->SetRangeUser(0,10);
      
      h3D_yxproj[eye] = (TH2D*)h3D[eye]->Project3D("yx");      
      h3D_yxproj[eye]->Draw("colz");
      
      cQGated[eye]->cd();
      
      Double_t temp_total=0;
      Double_t temp_bkg=0;
      
      for(int jay=0; jay<3; jay++) {
	hQgated[eye][jay] = ((TH3D*)histodir->FindObjectAny(Form("En_Ecl_Mcl_Qgated_%d",jay)));
	//	cout<<hQgated[eye][jay]<<endl;
	hQgated[eye][jay]->GetZaxis()->SetRangeUser(mcl_low,mcl_high);
	hQgated[eye][jay]->GetYaxis()->SetRangeUser(0,10);
	
	hQgated_xproj[eye][jay] = (TH1D*)hQgated[eye][jay]->Project3D("x");
	
	
	if(jay==0) {
	  h1D_bkgsub[eye] = (TH1D*)hQgated_xproj[eye][jay]->Clone();
	  // cout<<"clone: "<<h1D_bkgsub[eye]<<endl;
	  temp_total = hQgated_xproj[eye][jay]->Integral(1480,1680);
	  hQgated_xproj[eye][jay]->SetLineColor(4);
	  hQgated_xproj[eye][jay]->Draw();
	}
	if(jay==1) {
	  hQgated_xproj[eye][jay]->SetLineColor(3);
	  hQgated_xproj[eye][jay]->Draw("same");
	}
	if(jay==2) {
	  hQgated_xproj[eye][jay]->SetLineColor(2);
	  hQgated_xproj[eye][jay]->Draw("same");
	  temp_bkg = hQgated_xproj[eye][jay]->Integral(1480,1680);
	  
	  double scale = -1 - temp_bkg/temp_total;
	  
	  h1D_bkgsub[eye]->Add(hQgated_xproj[eye][jay],scale);
	  h1D_bkgsub[eye]->SetLineColor(1);
	  h1D_bkgsub[eye]->Draw("same");
	}
      } 

      //Determine 4.9 eV Resonance counts
      
      Peak_Integral[eye] = h1D_bkgsub[eye]->Integral(1187,1204);
      Peak_Error[eye] = TMath::Power(Peak_Integral[eye],0.5);
      cout<<eye<<"  "<<Peak_Integral[eye]<<endl;
      
      
      Phi_Disk[eye] = 4.0/(TMath::Pi() * TMath::Power(Diameter[eye],2.0))*(Peak_Integral[eye]/Li6_Counts[eye]);
      Phi_Disk[eye] /= Densities[eye];
      
      double ErrorFactor = TMath::Power((1/Li6_Counts[eye]+1/Peak_Integral[eye]),0.5);
      
      Phi_Disk_Error[eye] = Phi_Disk[eye]*ErrorFactor;

      cout<<eye<<" Diameter: "<<Diameter[eye]<<" Phi Disk: "<<Phi_Disk[eye]<<" +/- "<<Phi_Disk_Error[eye]<<endl;
      
      
      //  } //end check over fin[eye]
    
  }
    
  double Correction = 1.0/Phi_Disk[0];
  
  for(int eye=0; eye<5; eye++) {
    Phi_Disk[eye] *= Correction; //Normalize to 1 for smallest diameter
    Phi_Disk_Error[eye] *= Correction;
    cout<<eye<<" Diameter: "<<Diameter[eye]<<" Phi Disk: "<<Phi_Disk[eye]<<" +/- "<<Phi_Disk_Error[eye]<<endl;
 
  }

  
  //Mertz report for 4.9eV
  Double_t Mertz_Flux[5] = {1.000, 1.017, 0.959, 0.841, 0.520};
  Double_t Mertz_Error[5] = {0.009, 0.009, 0.005, 0.005, 0.002};
  
  
  TGraphErrors *gMertz = new TGraphErrors(5,Radius,Mertz_Flux,0,Mertz_Error);
  gMertz->SetMarkerStyle(20);

  TGraphErrors *gProfile = new TGraphErrors(5,Radius,Phi_Disk,0,Phi_Disk_Error);
  gProfile->SetMarkerStyle(21);
  gProfile->SetMarkerColor(2);
  gProfile->SetMarkerSize(3);
  gProfile->GetXaxis()->SetTitle("Radius (mm)");
  gProfile->GetYaxis()->SetTitle("Normalized Flux (Arb. Units)");
  

  TCanvas *cProfile = new TCanvas("cProfile","cProfile");
  gProfile->Draw("AP");
  gMertz->Draw("PSame");
  gMertz->SetMarkerSize(2.5);


}
*/
