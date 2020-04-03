void TOF_Corrections() {

  int sigma = 1;
  
  gStyle->SetPalette(54);

  //Read in the root file with the dT/T vs En for a 10m flight path from Mocko

  TFile *fin = new TFile("ModeratorFunction_10MeV.root");

  TH2D *hModFunc = (TH2D*)fin->FindObjectAny("ModeratorFunc");
  cout<<hModFunc<<endl;

  TCanvas *c1 = new TCanvas("ModeratorFunction_10m","ModeratorFunction_10m");
  
  hModFunc->Draw("colz");

  
  //dT is fixed for a given neutron energy and the time of flight is proportional to the flight path length so the ratio (T2 is TOF of non 10m flight path) dt/T2 = (dt/T1)*(T1/T2) which is like (dt/T1)*(L1/L2) (T1 is TOF for 10m flight path and L1 is 10m and L2 is some other distance)

  double DANCE_FP_Length = 20.2572; //meters
  double U235_FP_Length = 22.824; //meters
  double He3_FP_Length = 22.74; //meters
  double Li6_FP_Length = 22.6121; //meters

  
  int nbinsx = hModFunc->GetNbinsX();
  const int N = nbinsx+1;  //The plus 1 is to set a high data point 
  cout<<nbinsx<<endl;
  
  TH1D *hProjy[N];
  double centroid[N];
  double Energy[N];
  double Median[N];
  double Low[N];
  double High[N];

  //Median dt/T corrected for FP length
  double DANCE_TOF[N];
  double U235_TOF[N];
  double He3_TOF[N];
  double Li6_TOF[N];

  //Time of flight from neutron energy
  double DANCE_Median[N];
  double U235_Median[N];
  double He3_Median[N];
  double Li6_Median[N];

  double DANCE_Low[N];
  double U235_Low[N];
  double He3_Low[N];
  double Li6_Low[N];

  double DANCE_High[N];
  double U235_High[N];
  double He3_High[N];
  double Li6_High[N];


  //Moderation time
  double DANCE_TOF_Moderation[N];
  double U235_TOF_Moderation[N];
  double He3_TOF_Moderation[N];
  double Li6_TOF_Moderation[N];

  //Time of flight plus moderation time
  double DANCE_TOF_Measured[N];
  double U235_TOF_Measured[N];
  double He3_TOF_Measured[N];
  double Li6_TOF_Measured[N];

  //Time of flight plus moderation time plus + error
  double DANCE_TOF_Measured_Plus[N];
  double U235_TOF_Measured_Plus[N];
  double He3_TOF_Measured_Plus[N];
  double Li6_TOF_Measured_Plus[N];

  //Time of flight plus moderation time plus - error
  double DANCE_TOF_Measured_Minus[N];
  double U235_TOF_Measured_Minus[N];
  double He3_TOF_Measured_Minus[N];
  double Li6_TOF_Measured_Minus[N];


  //Useful constants
  double MnCsquared = 939.56563; // MeV/c^2
  double speedoflight = 299792458.0; //m/s
  
  //Temp variables for calculating stuff
  double temp_cent=0;
  double temp_bin=0;
  double total_counts=0;
  double temp_counts=0;

  //Total Errors
  double DANCE_TOF_Error_Plus[N];
  double DANCE_TOF_Error_Minus[N];
  

  DANCE_TOF_Measured[0] = DANCE_TOF[0] = 100000000.0;  //1/10 of a second (beam is 20Hz so never should have TOF > 0.05s
  U235_TOF_Measured[0] = U235_TOF[0] = 100000000.0;
  He3_TOF_Measured[0] = He3_TOF[0] = 100000000.0;
  Li6_TOF_Measured[0] = Li6_TOF[0] = 100000000.0;
  
  DANCE_TOF_Measured_Plus[0] = 100000000.0+0.0014*100000000.0;  //1/10 of a second (beam is 20Hz so never should have TOF > 0.05s
  DANCE_TOF_Measured_Minus[0] = 100000000.0-0.0004*100000000.0;  //1/10 of a second (beam is 20Hz so never should have TOF > 0.05s



  for(int eye=1; eye<N; eye++) {
    hProjy[eye] = hModFunc->ProjectionY(Form("hProjy_%d",eye),eye+1,eye+1); 
    Energy[eye] = hModFunc->GetXaxis()->GetBinCenter(eye+1);

    //   double kinematic_factor = TMath::Sqrt(TMath::Power(MnCsquared,2)+2*MnCsquared*Energy[eye]+TMath::Power(Energy[eye],2))/(speedoflight*TMath::Sqrt(Energy[eye])*TMath::Sqrt(2*MnCsquared+Energy[eye]));
    
    double kinematic_factor = TMath::Sqrt( (MnCsquared) / (2 * Energy[eye] * TMath::Power(speedoflight,2)));
    
    DANCE_TOF[eye] = kinematic_factor*DANCE_FP_Length*1e9; //in ns
    U235_TOF[eye] = kinematic_factor*U235_FP_Length*1e9; //in ns
    He3_TOF[eye] = kinematic_factor*He3_FP_Length*1e9; //in ns
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
    
    //Determine the median
    for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) { 
      temp_counts += hProjy[eye]->GetBinContent(jay+1);
      if(temp_counts >= 0.5*total_counts) {
	Median[eye] = hProjy[eye]->GetBinCenter(jay+1);
	break;
      }
    }

    temp_counts =0;
    
    if(sigma==1) {
      //Determine the lower bound 1 sigma
      for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) { 
	temp_counts += hProjy[eye]->GetBinContent(jay+1);
	if(temp_counts >= 0.159*total_counts) {
	Low[eye] = hProjy[eye]->GetBinCenter(jay+1);
	break;
	}
      }
      
      temp_counts =0;
      
      //Determine the upper bound 1 sigma
      for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) { 
	temp_counts += hProjy[eye]->GetBinContent(jay+1);
	if(temp_counts >= 0.841*total_counts) {
	  High[eye] = hProjy[eye]->GetBinCenter(jay+1);
	  break;
	}
      }
    }
    if(sigma==2) {
       //Determine the lower bound 2 sigma
      for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) { 
	temp_counts += hProjy[eye]->GetBinContent(jay+1);
	if(temp_counts >= 0.025*total_counts) {
	Low[eye] = hProjy[eye]->GetBinCenter(jay+1);
	break;
	}
      }
      
      temp_counts =0;
      
      //Determine the upper bound 2 sigma
      for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) { 
	temp_counts += hProjy[eye]->GetBinContent(jay+1);
	if(temp_counts >= 0.975*total_counts) {
	  High[eye] = hProjy[eye]->GetBinCenter(jay+1);
	  break;
	}
      }
    }
    if(sigma==3) {
      //Determine the lower bound 3 sigma
      for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) { 
	temp_counts += hProjy[eye]->GetBinContent(jay+1);
	if(temp_counts >= 0.0015*total_counts) {
	  Low[eye] = hProjy[eye]->GetBinCenter(jay+1);
	  break;
	}
      }
      
      temp_counts =0;
      
      //Determine the upper bound 3 sigma
      for(int jay=0; jay<hProjy[eye]->GetNbinsX(); jay++) { 
	temp_counts += hProjy[eye]->GetBinContent(jay+1);
	if(temp_counts >= 0.9985*total_counts) {
	  High[eye] = hProjy[eye]->GetBinCenter(jay+1);
	  break;
	}
      }
    }
    
    centroid[eye] = temp_cent/temp_bin;
    
    //Scale the dTOF/TOF by the FP lengths 
    DANCE_Median[eye] = Median[eye]*10.0/DANCE_FP_Length;
    U235_Median[eye] = Median[eye]*10.0/U235_FP_Length;
    He3_Median[eye] = Median[eye]*10.0/He3_FP_Length;
    Li6_Median[eye] = Median[eye]*10.0/Li6_FP_Length;
      
    //Scale the error bars too
    DANCE_Low[eye] = Low[eye]*10.0/DANCE_FP_Length;
    U235_Low[eye] = Low[eye]*10.0/U235_FP_Length;
    He3_Low[eye] = Low[eye]*10.0/He3_FP_Length;
    Li6_Low[eye] = Low[eye]*10.0/Li6_FP_Length;
    
    DANCE_High[eye] = High[eye]*10.0/DANCE_FP_Length;
    U235_High[eye] = High[eye]*10.0/U235_FP_Length;
    He3_High[eye] = High[eye]*10.0/He3_FP_Length;
    Li6_High[eye] = High[eye]*10.0/Li6_FP_Length;

    //Calculate a measured TOF using the dt/TOF and theoretical TOF
    DANCE_TOF_Measured[eye] = (1+DANCE_Median[eye])*DANCE_TOF[eye];
    U235_TOF_Measured[eye] = (1+U235_Median[eye])*U235_TOF[eye];
    He3_TOF_Measured[eye] = (1+He3_Median[eye])*He3_TOF[eye];
    Li6_TOF_Measured[eye] = (1+Li6_Median[eye])*Li6_TOF[eye];
    
    //Calculate a measured TOF using the dt/TOF and theoretical TOF
    DANCE_TOF_Measured_Plus[eye] = (1+DANCE_High[eye])*DANCE_TOF[eye];
    U235_TOF_Measured_Plus[eye] = (1+U235_High[eye])*U235_TOF[eye];
    He3_TOF_Measured_Plus[eye] = (1+He3_High[eye])*He3_TOF[eye];
    Li6_TOF_Measured_Plus[eye] = (1+Li6_High[eye])*Li6_TOF[eye];

    //Calculate a measured TOF using the dt/TOF and theoretical TOF
    DANCE_TOF_Measured_Minus[eye] = (1+DANCE_Low[eye])*DANCE_TOF[eye];
    U235_TOF_Measured_Minus[eye] = (1+U235_Low[eye])*U235_TOF[eye];
    He3_TOF_Measured_Minus[eye] = (1+He3_Low[eye])*He3_TOF[eye];
    Li6_TOF_Measured_Minus[eye] = (1+Li6_Low[eye])*Li6_TOF[eye];

    cout<<eye<<"  "<<Energy[eye]<<"  "<<centroid[eye]<<"  "<<Median[eye]<<" Low: "<<Low[eye]<<" High: "<<High[eye]<<"  "<< DANCE_TOF_Measured[eye]<<"  DANCE Low "<< DANCE_TOF_Measured_Minus[eye]<<"  DANCE High "<<DANCE_TOF_Measured_Plus[eye] <<endl;
    
  }

  TGraph *gr = new TGraph(N,Energy,Median);
  TGraph *gr2 = new TGraph(N,Energy,centroid);
  TGraph *gr3 = new TGraph(N,Energy,Low);
  TGraph *gr4 = new TGraph(N,Energy,High);

  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(2);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(1);
  gr3->SetMarkerStyle(21);
  gr3->SetMarkerColor(7);
  gr4->SetMarkerStyle(21);
  gr4->SetMarkerColor(7);
  c1->cd();
  gr->Draw("Psame");
  gr2->Draw("Psame");
  gr3->Draw("Psame");
  gr4->Draw("Psame");

  TGraph *gDANCE = new TGraph(N-1,DANCE_TOF_Measured,DANCE_TOF);
  TGraph *gU235 = new TGraph(N-1,U235_TOF_Measured,U235_TOF);
  TGraph *gHe3 = new TGraph(N-1,He3_TOF_Measured,He3_TOF);
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
  gHe3->Write();
  gLi6->Write();

  fout->Write();
  fout->Close();


  ofstream out;
  out.open("TOF_Corrections.txt");
  
  out<<N-1<<"\n";
  for(int eye=0; eye<N-1; eye++) {
    out<<DANCE_TOF[eye]<<"  "<<DANCE_TOF_Measured[eye]<<"   "<<U235_TOF[eye]<<"  "<<U235_TOF_Measured[eye]<<"   "<<Li6_TOF[eye]<<"  "<<Li6_TOF_Measured[eye]<<"  "<<He3_TOF[eye]<<"  "<<He3_TOF_Measured[eye]<<"\n";
  }
  out<<"DANCE  DANCE_Measured    U235  U235_Measured    Li6   Li6_Measured    He3   He3_Measured \n";

  out.close();


  //Deal with errors in TOF
 
  double DANCE_Response_Error = 0;
  double Proton_Pulse_Error = 0;

  if(sigma==1) {
    DANCE_Response_Error = 0.587; //sigma in ns from time deviation plots
    // Solve(integrate(120-x,x)==4915.364342,x)
    Proton_Pulse_Error = 52.4036; //width in ns of the central part of the triangle where 68.2% lies (FWHM of triangle is 120 ns)
  }
  if(sigma==2) {
    DANCE_Response_Error = 0.587*2.0; //sigma in ns from time deviation plots
     // Solve(integrate(120-x,x)==6872.398099,x)
   Proton_Pulse_Error = 94.4031; //width in ns of the central part of the triangle where 68.2% lies (FWHM of triangle is 120 ns)
  }
  if(sigma==3) {
    DANCE_Response_Error = 0.587*3.0; //sigma in ns from time deviation plots
     // Solve(integrate(120-x,x)==7180.561462,x)
    Proton_Pulse_Error = 113.765; //width in ns of the central part of the triangle where 68.2% lies (FWHM of triangle is 120 ns)
  }


  //Set Energy 0
  Energy[0] = 0.5*MnCsquared/TMath::Power(speedoflight,2)*TMath::Power((DANCE_FP_Length/(DANCE_TOF[0]/1e9)),2);

  for(int eye=0; eye<N; eye++) {
    //Calculate DANCE Error
    
    DANCE_TOF_Error_Plus[eye] = TMath::Sqrt( TMath::Power(DANCE_Response_Error,2) +
					     TMath::Power(Proton_Pulse_Error,2) + 
					     TMath::Power((DANCE_TOF_Measured_Plus[eye]-DANCE_TOF_Measured[eye]),2));
    
    DANCE_TOF_Error_Minus[eye] = TMath::Sqrt( TMath::Power(DANCE_Response_Error,2) +
					      TMath::Power(Proton_Pulse_Error,2) + 
					      TMath::Power((DANCE_TOF_Measured_Minus[eye]-DANCE_TOF_Measured[eye]),2));
    
    cout<<eye<<" Energy: "<<Energy[eye]<<"  TOF: "<<DANCE_TOF_Measured[eye]<<" + "<<DANCE_TOF_Error_Plus[eye]<<" ("<<100*DANCE_TOF_Error_Plus[eye]/DANCE_TOF_Measured[eye]<<" %) /- "<<DANCE_TOF_Error_Minus[eye]<<" ("<<100*DANCE_TOF_Error_Minus[eye]/DANCE_TOF_Measured[eye]<<" %)"<<endl; 
    
  }
  
  TGraph *gDANCE_Plus = new TGraph(N-1,DANCE_TOF_Measured,DANCE_TOF_Error_Plus);
  TGraph *gDANCE_Minus = new TGraph(N-1,DANCE_TOF_Measured,DANCE_TOF_Error_Minus);
  TGraph *gDANCE_Value_vs_Energy = new TGraph(N-1,Energy,DANCE_TOF_Measured);

  TCanvas *c2 = new TCanvas("c2","c2");
  gDANCE_Plus->Draw("A*");
  gDANCE_Minus->Draw("*SAME");

  TCanvas *c3 = new TCanvas("c3","c3");
  gDANCE_Value_vs_Energy->Draw("A*");

  // TCanvas *cTOF_vs_E= new TCanvas("TOF_vs_E","TOF_vs_E");
  //  gDANCE_Value_vs_Energy->Draw("A*");

  //  return 0;

  //Physics Histograms
  double x[5000];
  int NEbins=0;

  double NeutronE_From = 0.002e-6;
  double NeutronE_To = 5.0;
  double NeutronE_BinsPerDecade = 500.0;
    
  for(double lx=log10(NeutronE_From);lx<log10(NeutronE_To);lx=lx+(1./NeutronE_BinsPerDecade)){
    x[NEbins]=pow(10,lx);
    NEbins++;
  }
  NEbins--;

  
  const int NBINS = NEbins;

  double dE_dTOF[NBINS];
  double TOF_error_plus[NBINS];
  double TOF_error_minus[NBINS];
  double energy_val[NBINS];
  double TOF_val[NBINS];
  double sig_E_plus[NBINS];
  double sig_E_minus[NBINS];
  double sig_E_plus_percent[NBINS];
  double sig_E_minus_percent[NBINS];


  //Calculate Error in energy
  
  TH1D *hErrorPlus = new TH1D(Form("ErrorPlus_%dSigma",sigma),Form("ErrorPlus_%dSigma",sigma),NEbins,x);
  TH1D *hErrorMinus = new TH1D(Form("ErrorMinus_%dSigma",sigma),Form("ErrorMinus_%dSigma",sigma),NEbins,x);

  for(int eye=0; eye<NEbins; eye++) {

    //See what energy we want (in MeV)
    energy_val[eye] = hErrorPlus->GetXaxis()->GetBinCenter(eye+1);
    
    //Get the TOF from energy
    TOF_val[eye] = gDANCE_Value_vs_Energy->Eval(energy_val[eye]);

    
    //calculate derrivative of E WRT TOF
    dE_dTOF[eye] =  (MnCsquared / TMath::Power((speedoflight * 1e-9),2)) * TMath::Power(DANCE_FP_Length,2) / TMath::Power(TOF_val[eye],3);
    
    //Get the errors for this TOF val
    TOF_error_plus[eye] = gDANCE_Plus->Eval(TOF_val[eye]);
    TOF_error_minus[eye] = gDANCE_Minus->Eval(TOF_val[eye]);

    cout<<energy_val[eye]<<"  "<<TOF_val[eye]<<"  "<<TOF_error_plus[eye]<<" ("<<100*TOF_error_plus[eye]/TOF_val[eye]<<"%) "<<TOF_error_minus[eye]<<" ("<<100*TOF_error_minus[eye]/TOF_val[eye]<<"%) "<<"  "<<dE_dTOF[eye]<<endl;

    //Calculate Error in energy
    sig_E_plus[eye] = TMath::Sqrt( TMath::Power(TOF_error_plus[eye],2) * TMath::Power(dE_dTOF[eye],2));
    sig_E_minus[eye] = -1*TMath::Sqrt( TMath::Power(TOF_error_minus[eye],2) * TMath::Power(dE_dTOF[eye],2));



    sig_E_plus_percent[eye] = 100*sig_E_plus[eye]/energy_val[eye];
    sig_E_minus_percent[eye] = 100*sig_E_minus[eye]/energy_val[eye];

    hErrorPlus->Fill(energy_val[eye],sig_E_plus_percent[eye]);
    hErrorMinus->Fill(energy_val[eye],sig_E_minus_percent[eye]);


  } 

  TCanvas *c1 = new TCanvas("c1","c1");

  hErrorPlus->Draw();
  hErrorMinus->Draw("same");

  TFile *fout_error = new TFile(Form("Energy_Error_%dSigma.root",sigma),"RECREATE");

  hErrorPlus->Write();
  hErrorMinus->Write();
  fout_error->Write();
  fout_error->Close();


  
}





