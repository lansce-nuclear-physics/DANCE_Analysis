



void TOF_Corrections() {

  
  //Read in the root file with the dT/T vs En for a 10m flight path from Mocko

  TFile *fin = new TFile("ModeratorFunction_10MeV.root");

  TH2D *hModFunc = (TH2D*)fin->FindObjectAny("ModeratorFunc");
  cout<<hModFunc<<endl;

  TCanvas *c1 = new TCanvas("ModeratorFunction_10m","ModeratorFunction_10m");
  
  hModFunc->Draw("colz");

  
  //dT is fixed for a given neutron energy and the time of flight is proportional to the flight path length so the ratio (T2 is TOF of non 10m flight path) dt/T2 = (dt/T1)*(T1/T2) which is like (dt/T1)*(L1/L2) (T1 is TOF for 10m flight path and L1 is 10m and L2 is some other distance)

  double DANCE_FP_Length = 20.28437; //meters
  double U235_FP_Length = 22.8; //meters
  double BF3_FP_Length = 22.74; //meters
  double Li6_FP_Length = 22.607; //meters

  
  int nbinsx = hModFunc->GetNbinsX();
  const int N = nbinsx+1;  //The plus 1 is to set a high data point 
  cout<<nbinsx<<endl;
  
  TH1D *hProjy[N];
  double centroid[N];
  double Energy[N];
  double Median[N];

  //Median dt/T corrected for FP length
  double DANCE_TOF[N];
  double U235_TOF[N];
  double BF3_TOF[N];
  double Li6_TOF[N];

  //Time of flight from neutron energy
  double DANCE_Median[N];
  double U235_Median[N];
  double BF3_Median[N];
  double Li6_Median[N];

  //Moderation time
  double DANCE_TOF_Moderation[N];
  double U235_TOF_Moderation[N];
  double BF3_TOF_Moderation[N];
  double Li6_TOF_Moderation[N];

  //Time of flight plus moderation time
  double DANCE_TOF_Measured[N];
  double U235_TOF_Measured[N];
  double BF3_TOF_Measured[N];
  double Li6_TOF_Measured[N];

  //Useful constants
  double MnCsquared = 939.56563; // MeV/c^2
  double speedoflight = 299792458.0; //m/s
  
  //Temp variables for calculating stuff
  double temp_cent=0;
  double temp_bin=0;
  double total_counts=0;
  double temp_counts=0;


  DANCE_TOF_Measured[0] = DANCE_TOF[0] = 100000000.0;  //1/10 of a second (beam is 20Hz so never should have TOF > 0.05s
  U235_TOF_Measured[0] = U235_TOF[0] = 100000000.0;
  BF3_TOF_Measured[0] = BF3_TOF[0] = 100000000.0;
  Li6_TOF_Measured[0] = Li6_TOF[0] = 100000000.0;
  

  for(int eye=1; eye<N; eye++) {
    hProjy[eye] = hModFunc->ProjectionY(Form("hProjy_%d",eye),eye+1,eye+1); 
    Energy[eye] = hModFunc->GetXaxis()->GetBinCenter(eye+1);

    double kinematic_factor = TMath::Sqrt(TMath::Power(MnCsquared,2)+2*MnCsquared*Energy[eye]+TMath::Power(Energy[eye],2))/(speedoflight*TMath::Sqrt(Energy[eye])*TMath::Sqrt(2*MnCsquared+Energy[eye]));
    
    DANCE_TOF[eye] = kinematic_factor*DANCE_FP_Length*1e9; //in ns
    U235_TOF[eye] = kinematic_factor*U235_FP_Length*1e9; //in ns
    BF3_TOF[eye] = kinematic_factor*BF3_FP_Length*1e9; //in ns
    Li6_TOF[eye] = kinematic_factor*Li6_FP_Length*1e9; //in ns
    
    temp_cent=0;
    temp_bin=0;
    total_counts=0;
    temp_counts=0;
    
    for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) {
      // cout<<hProjy[eye]->GetBinContent(jay+1)<<"  "<<hProjy[eye]->GetBinCenter(jay+1)<<endl;
      temp_cent += hProjy[eye]->GetBinContent(jay+1)*hProjy[eye]->GetBinCenter(jay+1);
      temp_bin += hProjy[eye]->GetBinContent(jay+1);
      total_counts += hProjy[eye]->GetBinContent(jay+1);
    }
    
    cout<<total_counts<<"  "<<temp_bin<<endl;
    
    for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) { 
      temp_counts += hProjy[eye]->GetBinContent(jay+1);
      if(temp_counts >= 0.5*total_counts) {
	Median[eye] = hProjy[eye]->GetBinCenter(jay+1);
	break;
      }
    }
    
    centroid[eye] = temp_cent/temp_bin;
    cout<<eye<<"  "<<Energy[eye]<<"  "<<centroid[eye]<<"  "<<Median[eye]<<endl;
    
    DANCE_Median[eye] = Median[eye]*10.0/DANCE_FP_Length;
    U235_Median[eye] = Median[eye]*10.0/U235_FP_Length;
    BF3_Median[eye] = Median[eye]*10.0/BF3_FP_Length;
    Li6_Median[eye] = Median[eye]*10.0/Li6_FP_Length;
      
    DANCE_TOF_Measured[eye] = (1+DANCE_Median[eye])*DANCE_TOF[eye];
    U235_TOF_Measured[eye] = (1+U235_Median[eye])*U235_TOF[eye];
    BF3_TOF_Measured[eye] = (1+BF3_Median[eye])*BF3_TOF[eye];
    Li6_TOF_Measured[eye] = (1+Li6_Median[eye])*Li6_TOF[eye];
        
  }

  TGraph *gr = new TGraph(N,Energy,Median);
  TGraph *gr2 = new TGraph(N,Energy,centroid);
  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(2);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(1);
  c1->cd();
  gr->Draw("Psame");
  gr2->Draw("Psame");

  TGraph *gDANCE = new TGraph(N-1,DANCE_TOF_Measured,DANCE_TOF);
  TGraph *gU235 = new TGraph(N-1,U235_TOF_Measured,U235_TOF);
  TGraph *gBF3 = new TGraph(N-1,BF3_TOF_Measured,BF3_TOF);
  TGraph *gLi6 = new TGraph(N-1,Li6_TOF_Measured,Li6_TOF);
  gDANCE->SetMarkerStyle(22);
  
  TCanvas *c2 = new TCanvas();
  gDANCE->Draw("AP");

  
  for(int eye=0; eye<N-1; eye++) {
    // cout<<Energy[eye]<<"  "<<DANCE_TOF[eye]<<"  "<<DANCE_TOF_Measured[eye]<<"  "<<DANCE_TOF_Measured[eye]-DANCE_TOF[eye]<<endl;    
    cout<<Energy[eye]<<"  "<<U235_TOF[eye]<<"  "<<U235_TOF_Measured[eye]<<"  "<<U235_TOF_Measured[eye]-U235_TOF[eye]<<endl;    
    if(eye<N-1) {
      // if(DANCE_TOF_Measured[eye+1]>DANCE_TOF_Measured[eye]) {
      if(U235_TOF_Measured[eye+1]>U235_TOF_Measured[eye]) {
	//non monotonic and thus a problem
	cout<<"PROBLEM: "<<eye<<endl;
      }
    }
  }
  
  
  cout<<gDANCE->Eval(2000000)<<endl;

  
  TFile *fout = new TFile("TOF_Corrections.root","RECREATE");
  
  fout->cd();
  gDANCE->Write();
  gU235->Write();
  gBF3->Write();
  gLi6->Write();

  fout->Write();
  fout->Close();


  ofstream out;
  out.open("TOF_Corrections.txt");
  
  out<<N-1<<"\n";
  for(int eye=0; eye<N-1; eye++) {
    out<<DANCE_TOF[eye]<<"  "<<DANCE_TOF_Measured[eye]<<"   "<<U235_TOF[eye]<<"  "<<U235_TOF_Measured[eye]<<"   "<<Li6_TOF[eye]<<"  "<<Li6_TOF_Measured[eye]<<"  "<<BF3_TOF[eye]<<"  "<<BF3_TOF_Measured[eye]<<"\n";
  }
  out<<"DANCE  DANCE_Measured    U235  U235_Measured    Li6   Li6_Measured    BF3   BF3_Measured \n";

  out.close();



}





