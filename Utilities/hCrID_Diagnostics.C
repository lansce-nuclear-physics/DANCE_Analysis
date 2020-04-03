#include <iostream>
#include <sstream>
#include <fstream>

void hCrID_Diagnostics() {

  gStyle->SetOptStat(0);
  
  //Location of stage0 rootfiles
  stringstream fpath;
  fpath.str();
  fpath << "/home/cfry/DANCE_Analysis/stage0_root/";

 //Location of stage1 rootfiles
  stringstream fpath1;
  fpath1.str();
  fpath1 << "/home/cfry/DANCE_Analysis/stage1_root/";

  //stage0
  int CW = 500;
  int CBT = 0;
  int DEBT = 0;
  
  //stage1
  int CW1 = 10;
  int CBT1 = 2000;
  int DEBT1 = 0;

  double fraction = 0.75;  //If the value is less than the average by more than this factor exclude

  double average_T0=0;
  double average_Li6=0;
  double average_U235=0;
  double average_He3=0;

  vector<int> exclude_T0;
  vector<int> exclude_Li6;
  vector<int> exclude_U235;
  vector<int> exclude_DD;
  
  //These are the Au197 runs

  //  int start_run = 60572;
  //  int end_run = 60584;

  //   int start_run = 60409;
  //   int end_run = 60451;

  //  int start_run = 63590;
  //  int end_run = 63649;

  //  int start_run = 68104;
  //  int end_run = 68126;
  
  ////////////////////////////
  
  //These are the Cu65 runs
  
  //  int start_run = 60834;
  //  int end_run = 61251;

  //int start_run = 63015;
  //int end_run = 63589;
   
  // int start_run = 64138;
  // int end_run = 64692;

  // int start_run = 64989;
  //  int end_run = 65265;

  ////////////////////////


  //These are the Pb208 runs
  
  //  int start_run = 60622;
  //  int end_run = 60811;

  //  int start_run = 62695;
  // int end_run = 62735;
   
  //  int start_run = 63650;
     // int end_run = 64136;

  //   int start_run = 64834;
  //  int end_run = 64987;

  //   int start_run = 68040;
  //  int end_run = 68103;

  ////////////////////////


  //These are the Cu63 runs
  
    int start_run = 64693;
    int end_run = 64833;
  
  ////////////////////////


  //Output file
  TFile *fout = new TFile(Form("Histograms_Runs_%d_%d.root",start_run, end_run),"RECREATE");
  
  
  int last_drop = start_run; 
  vector<int> time_between_drops;
  int drop_counter[50];
  int total_drop_counter = 0;
  int ch_drop_counter[50][162];
  int run_counter[50];
  int total_run_counter =0;
  int group_counter;
    
  int N_Runs_per_group = 50;
  
  const int N_Runs = end_run-start_run+1;
  
  cout<<N_Runs<<endl;

  TFile *fin[N_Runs];
  TFile *fin1[N_Runs];
  
  TH1D *hID;
  TH1D *hScalers;
  TH1D *hLi6_PH;
  TH1D *hU235_PH;
  TH3D *hEn_Etot_Mcl;

  TH1D *hLi6_En;
  TH1D *hU235_En;
  TH1D *hHe3_En;

  double Counts[N_Runs][162];
  int Run_Num[N_Runs];
  double charge=0;
  
  int counter=0;
  
  TH2D *hCounts = new TH2D("hCounts","hCounts",N_Runs,start_run,end_run+1,162,0,162);


  TH1D *hT0 = new TH1D("T0","T0",N_Runs,start_run,end_run+1);
  TH1D *hLi6 = new TH1D("Li6","Li6",N_Runs,start_run,end_run+1);
  TH1D *hU235 = new TH1D("U235","U235",N_Runs,start_run, end_run+1);
  TH1D *hDANCE = new TH1D("DANCE","DANCE",N_Runs,start_run,end_run+1);
  TH1D *hHe3 = new TH1D("He3","He3",N_Runs,start_run, end_run+1);

  hT0->GetXaxis()->SetTitle("Run Number");
  hLi6->GetXaxis()->SetTitle("Run Number"); 
  hU235->GetXaxis()->SetTitle("Run Number"); 
  hDANCE->GetXaxis()->SetTitle("Run Number");
  hHe3->GetXaxis()->SetTitle("Run Number");

  hT0->GetYaxis()->SetTitle("Number of T0s");
  hLi6->GetYaxis()->SetTitle("Li6 Counts"); 
  hU235->GetYaxis()->SetTitle("U235 Counts"); 
  hDANCE->GetYaxis()->SetTitle("DANCE Counts");
  hHe3->GetYaxis()->SetTitle("He3 Counts");

  TH1D *hLi6_Gated_1to100eV = new TH1D("Li6_Gated_1to100eV ","Li6_Gated_1to100eV ",N_Runs,start_run,end_run+1);
  TH1D *hU235_Gated_1to100eV = new TH1D("U235_Gated_1to100eV ","U235_Gated_1to100eV ",N_Runs,start_run, end_run+1);
  TH1D *hHe3_Gated_1to100eV = new TH1D("He3_Gated_1to100eV ","He3_Gated_1to100eV ",N_Runs,start_run, end_run+1);

  hLi6_Gated_1to100eV->GetXaxis()->SetTitle("Run Number"); 
  hU235_Gated_1to100eV->GetXaxis()->SetTitle("Run Number"); 
  hHe3_Gated_1to100eV->GetXaxis()->SetTitle("Run Number");

  hLi6_Gated_1to100eV->GetYaxis()->SetTitle("Li6 Counts"); 
  hU235_Gated_1to100eV->GetYaxis()->SetTitle("U235 Counts"); 
  hHe3_Gated_1to100eV->GetYaxis()->SetTitle("He3 Counts");



  TH1D *hDANCE_per_T0 = new TH1D("DANCE_per_T0","DANCE_per_T0",N_Runs,start_run,end_run+1);
  TH1D *hLi6_per_T0 = new TH1D("Li6_per_T0","Li6_per_T0",N_Runs,start_run,end_run+1);
  TH1D *hU235_per_T0 = new TH1D("U235_per_T0","U235_per_T0",N_Runs,start_run,end_run+1);
  TH1D *hHe3_per_T0 = new TH1D("He3_per_T0","He3_per_T0",N_Runs,start_run,end_run+1);

  hLi6_per_T0->GetXaxis()->SetTitle("Run Number"); 
  hU235_per_T0->GetXaxis()->SetTitle("Run Number"); 
  hDANCE_per_T0->GetXaxis()->SetTitle("Run Number");
  hHe3_per_T0->GetXaxis()->SetTitle("Run Number");

  hLi6_per_T0->GetYaxis()->SetTitle("Li6 Counts per T0"); 
  hU235_per_T0->GetYaxis()->SetTitle("U235 Counts per T0"); 
  hDANCE_per_T0->GetYaxis()->SetTitle("DANCE Counts per T0");
  hHe3_per_T0->GetYaxis()->SetTitle("He3 Counts per T0");


  TH1D *hLi6_per_U235 = new TH1D("Li6_per_U235","Li6_per_U235",N_Runs,start_run, end_run+1);
  TH1D *hHe3_per_U235 = new TH1D("He3_per_U235","He3_per_U235",N_Runs,start_run, end_run+1);

  hLi6_per_U235->GetXaxis()->SetTitle("Run Number"); 
  hHe3_per_U235->GetXaxis()->SetTitle("Run Number");

  hLi6_per_U235->GetYaxis()->SetTitle("Li6 Counts per U235"); 
  hHe3_per_U235->GetYaxis()->SetTitle("He3 Counts per U235");


  TH1D *hCharge = new TH1D("Charge","Charge",N_Runs,start_run, end_run+1);
  TH1D *hCharge_per_T0 = new TH1D("Charge_per_T0","Charge_per_T0",N_Runs,start_run,end_run+1);

  hCharge->GetXaxis()->SetTitle("Run Number");
  hCharge_per_T0->GetXaxis()->SetTitle("Run Number");

  hCharge->GetYaxis()->SetTitle("Charge Scaler Counts");
  hCharge_per_T0->SetTitle("Charge per T0"); 

  TH1D *hDANCE_Mult_2to6_per_Charge = new TH1D("DANCE_Mult_2to6_per_Charge","DANCE_Mult_2to6_per_Charge",N_Runs,start_run,end_run+1);
  hDANCE_Mult_2to6_per_Charge->GetXaxis()->SetTitle("Run Number");
  hDANCE_Mult_2to6_per_Charge->GetYaxis()->SetTitle("DANCE_Mult_2to6 Counts per Charge");

  TH1D *hDANCE_per_Charge = new TH1D("DANCE_per_Charge","DANCE_per_Charge",N_Runs,start_run,end_run+1);
  TH1D *hLi6_per_Charge = new TH1D("Li6_per_Charge","Li6_per_Charge",N_Runs,start_run,end_run+1);
  TH1D *hU235_per_Charge = new TH1D("U235_per_Charge","U235_per_Charge",N_Runs,start_run,end_run+1);
  TH1D *hHe3_per_Charge = new TH1D("He3_per_Charge","He3_per_Charge",N_Runs,start_run,end_run+1);

  hLi6_per_Charge->GetXaxis()->SetTitle("Run Number"); 
  hU235_per_Charge->GetXaxis()->SetTitle("Run Number"); 
  hDANCE_per_Charge->GetXaxis()->SetTitle("Run Number");
  hHe3_per_Charge->GetXaxis()->SetTitle("Run Number");

  hLi6_per_Charge->GetYaxis()->SetTitle("Li6 Counts per Charge"); 
  hU235_per_Charge->GetYaxis()->SetTitle("U235 Counts per Charge"); 
  hDANCE_per_Charge->GetYaxis()->SetTitle("DANCE Counts per Charge");
  hHe3_per_Charge->GetYaxis()->SetTitle("He3 Counts per Charge");


  TH1D *hLi6_Gated_1to100eV_per_Charge = new TH1D("Li6_Gated_1to100eV_per_Charge","Li6_Gated_1to100eV_per_Charge",N_Runs,start_run,end_run+1);
  TH1D *hU235_Gated_1to100eV_per_Charge = new TH1D("U235_Gated_1to100eV_per_Charge","U235_Gated_1to100eV_per_Charge",N_Runs,start_run,end_run+1);
  TH1D *hHe3_Gated_1to100eV_per_Charge = new TH1D("He3_Gated_1to100eV_per_Charge","He3_Gated_1to100eV_per_Charge",N_Runs,start_run,end_run+1);

  hLi6_Gated_1to100eV_per_Charge->GetXaxis()->SetTitle("Run Number"); 
  hU235_Gated_1to100eV_per_Charge->GetXaxis()->SetTitle("Run Number"); 
  hHe3_Gated_1to100eV_per_Charge->GetXaxis()->SetTitle("Run Number");

  hLi6_Gated_1to100eV_per_Charge->GetYaxis()->SetTitle("Li6 Counts Gated 1 to 100 eV per Charge"); 
  hU235_Gated_1to100eV_per_Charge->GetYaxis()->SetTitle("U235 Counts Gated 1 to 100 eV per Charge"); 
  hHe3_Gated_1to100eV_per_Charge->GetYaxis()->SetTitle("He3 Counts Gated 1 to 100 eV per Charge");


  TH1D *hDANCE_Mult_2to6_per_Charge_Proj = new TH1D("DANCE_Mult_2to6_per_Charge_Proj","DANCE_Mult_2to6_per_Charge_Proj",500,0,0.005);
  hDANCE_Mult_2to6_per_Charge_Proj->GetXaxis()->SetTitle("DANCE_Mult_2to6 Counts per Charge");
  hDANCE_Mult_2to6_per_Charge_Proj->GetYaxis()->SetTitle("Number of Runs");

  TH1D *hDANCE_per_Charge_Proj = new TH1D("DANCE_per_Charge_Proj","DANCE_per_Charge_Proj",500,0,5);
  TH1D *hLi6_per_Charge_Proj = new TH1D("Li6_per_Charge_Proj","Li6_per_Charge_Proj",250,0,0.05);
  TH1D *hU235_per_Charge_Proj = new TH1D("U235_per_Charge_Proj","U235_per_Charge_Proj",250,0,0.05);
  TH1D *hHe3_per_Charge_Proj = new TH1D("He3_per_Charge_Proj","He3_per_Charge_Proj",250,0,0.05);

  hLi6_per_Charge_Proj->GetXaxis()->SetTitle("Li6 Counts per Charge"); 
  hU235_per_Charge_Proj->GetXaxis()->SetTitle("U235 Counts per Charge"); 
  hDANCE_per_Charge_Proj->GetXaxis()->SetTitle("DANCE Counts per Charge");
  hHe3_per_Charge_Proj->GetXaxis()->SetTitle("He3 Counts per Charge");

  hLi6_per_Charge_Proj->GetYaxis()->SetTitle("Number of Runs"); 
  hU235_per_Charge_Proj->GetYaxis()->SetTitle("Number of Runs"); 
  hDANCE_per_Charge_Proj->GetYaxis()->SetTitle("Number of Runs");
  hHe3_per_Charge_Proj->GetYaxis()->SetTitle("Number of Runs");

  TH2D *hDANCE_Mult_vs_Run = new TH2D("DANCE_Mult_vs_Run","DANCE_Mult_vs_Run",N_Runs,start_run,end_run+1,8,0,8);



  TH2D *hLi6_PH_vs_Run = new TH2D("Li6_PH_vs_Run","Li6_PH_vs_Run",N_Runs,start_run,end_run+1,1500,0,75000);
  TH2D *hU235_PH_vs_Run = new TH2D("U235_PH_vs_Run","U235_PH_vs_Run",N_Runs,start_run,end_run+1,1500,0,75000);
  TH2D *hHe3_PH_vs_Run = new TH2D("He3_PH_vs_Run","He3_PH_vs_Run",N_Runs,start_run,end_run+1,1500,0,75000);

  hLi6_PH_vs_Run->GetXaxis()->SetTitle("Run Number"); 
  hU235_PH_vs_Run->GetXaxis()->SetTitle("Run Number"); 
  hHe3_PH_vs_Run->GetXaxis()->SetTitle("Run Number"); 
  
  hLi6_PH_vs_Run->GetYaxis()->SetTitle("Li6 Pulse Height"); 
  hU235_PH_vs_Run->GetYaxis()->SetTitle("U235 Pulse Height"); 
  hHe3_PH_vs_Run->GetYaxis()->SetTitle("He3 Pulse Height"); 

  
  ifstream testrun;
  bool fileexists=false;
  bool fileexists1=false;

  for(int eye=0; eye<50; eye++) {
    drop_counter[eye] = 0;
    run_counter[eye] = 0;
    for(int jay=0; jay<162; jay++) {
      ch_drop_counter[eye][jay] = 0;
    }
  }

  for(int eye=start_run; eye<=end_run; eye++) {




   stringstream teststring1;
    teststring1.str();
    teststring1 << fpath1.str() << "Stage1_Histograms_Run_"<< eye<<"_"<<CW1<<"ns_CW_"<<CBT1<<"ns_CBT_"<<DEBT1<<"ns_DEBT.root";
    cout<<eye<<"  "<<teststring1.str()<<endl;
    
    testrun.open(teststring1.str().c_str());
    fileexists=false;
    
    if(testrun.is_open()) {
      fileexists1=true;
      cout<<"Run: "<<eye<<" stage1 exists"<<endl;
      testrun.close();
    }
    
    if(fileexists1) { 

      fin1[counter] = new TFile(Form("%s/Stage1_Histograms_Run_%d_%dns_CW_%dns_CBT_%dns_DEBT.root",fpath1.str().c_str(),eye,CW1,CBT1,DEBT1)); 

      Run_Num[counter] = eye;

      //Get the En,Etot,Mcl
      hEn_Etot_Mcl = (TH3D*)fin1[counter]->FindObjectAny("En_Etot_Mcl");
      hEn_Etot_Mcl->GetXaxis()->SetRangeUser(3,6);  //set the range to cover the big resonance in gold
      TH1D *hEn_Etot_Mcl_zy = (TH1D*)hEn_Etot_Mcl->Project3D("z");
      for(int jay=0; jay<7; jay++) {
	hDANCE_Mult_vs_Run->Fill(eye,jay+1,hEn_Etot_Mcl_zy->GetBinContent(jay+1));
      }
      
      hLi6_En = (TH1D*)fin1[counter]->FindObjectAny("Li6_En_Corrected");
      hLi6_En->GetXaxis()->SetRangeUser(1,100);
      hLi6_Gated_1to100eV->Fill(eye,hLi6_En->Integral()); 

      hU235_En = (TH1D*)fin1[counter]->FindObjectAny("U235_En_Corrected");
      hU235_En->GetXaxis()->SetRangeUser(1,100);
      hU235_Gated_1to100eV->Fill(eye,hU235_En->Integral()); 

      hHe3_En = (TH1D*)fin1[counter]->FindObjectAny("He3_En_Corrected");
      hHe3_En->GetXaxis()->SetRangeUser(1,100);
      hHe3_Gated_1to100eV->Fill(eye,hHe3_En->Integral()); 



    }

    stringstream teststring;
    teststring.str();
    teststring << fpath.str() << "Stage0_Histograms_Run_"<< eye<<"_"<<CW<<"ns_CW_"<<CBT<<"ns_CBT_"<<DEBT<<"ns_DEBT.root";
    cout<<eye<<"  "<<teststring.str()<<endl;
    
    testrun.open(teststring.str().c_str());
    fileexists=false;
    
    if(testrun.is_open()) {
      fileexists=true;
      cout<<"Run: "<<eye<<" stage0 exists"<<endl;
      testrun.close();
    }
    
 



    if(fileexists) {  
      
      fin[counter] = new TFile(Form("%s/Stage0_Histograms_Run_%d_%dns_CW_%dns_CBT_%dns_DEBT.root",fpath.str().c_str(),eye,CW,CBT,DEBT)); 

      Run_Num[counter] = eye;

      //Get the scalers
      hScalers = (TH1D*)fin[counter]->FindObjectAny("Scalers");
      charge = hScalers->GetBinContent(6);

      hID = (TH1D*)fin[counter]->FindObjectAny("hID");
      cout<<hID<<endl;

      hLi6_PH = (TH1D*)fin[counter]->FindObjectAny("Li6_PulseHeight");
      cout<<hLi6_PH<<endl;
      for(int jay=0; jay<hLi6_PH->GetNbinsX(); jay++) {
	hLi6_PH_vs_Run->Fill(eye,hLi6_PH->GetXaxis()->GetBinCenter(jay+1),hLi6_PH->GetBinContent(jay+1));
      }

      hU235_PH = (TH1D*)fin[counter]->FindObjectAny("U235_PulseHeight");
      cout<<hU235_PH<<endl;
      for(int jay=0; jay<hU235_PH->GetNbinsX(); jay++) {
	hU235_PH_vs_Run->Fill(eye,hU235_PH->GetXaxis()->GetBinCenter(jay+1),hU235_PH->GetBinContent(jay+1));
      }

      hHe3_PH = (TH1D*)fin[counter]->FindObjectAny("He3_PulseHeight");
      cout<<hHe3_PH<<endl;
      for(int jay=0; jay<hHe3_PH->GetNbinsX(); jay++) {
	hHe3_PH_vs_Run->Fill(eye,hHe3_PH->GetXaxis()->GetBinCenter(jay+1),hHe3_PH->GetBinContent(jay+1));
      }

      run_counter[group_counter]++;
      total_run_counter++;
      bool missing_channels = false;
      
      for(int jay=0; jay<162; jay++) {
	Counts[counter][jay] = hID->GetBinContent(jay+1);
	if(jay<160) {
	  hDANCE->Fill(eye,Counts[counter][jay]);
	}
	hCounts->Fill(eye,jay,Counts[counter][jay]);
	
	if(Counts[counter][jay] == 0 && jay != 76 && jay !=86) {
	  ch_drop_counter[group_counter][jay] ++;
	  missing_channels = true;
	}
      }
      
      hCharge->Fill(eye,charge);

      //T0's and beam monitors
      hT0->Fill(eye,hID->GetBinContent(201));  // ID 200 bin 201
      average_T0 += hID->GetBinContent(201);
      hLi6->Fill(eye,hID->GetBinContent(245));  // ID 244 bin 245
      average_Li6 += hID->GetBinContent(245);
      hU235->Fill(eye,hID->GetBinContent(244)); // ID 243 bin 244
      average_U235 += hID->GetBinContent(244);
      hHe3->Fill(eye,hID->GetBinContent(242)); // ID 241 bin 242
      average_He3 += hID->GetBinContent(242);
      
      
      if(missing_channels) {
	exclude_DD.push_back(eye);
	time_between_drops.push_back(eye-last_drop);
	last_drop = eye;
	drop_counter[group_counter]++;
	total_drop_counter ++;
      }

      if(run_counter[group_counter] >=N_Runs_per_group) {
	group_counter++;
      }

      // hID[counter]->Delete();
      fin[counter]->Close();
    }
    else {
      for(int jay=0; jay<162; jay++) {
	Counts[counter][jay] = -1;
      }
    }    
    counter++;
  }

  average_T0 /= N_Runs;
  average_Li6 /= N_Runs;
  average_U235 /= N_Runs;
  average_He3 /= N_Runs;

  
  hCounts->Draw("colz");

  TLine *line[15];

  for(int eye=0; eye<10; eye++) {
    line[eye] = new TLine(start_run,(eye+1)*16,end_run+1,(eye+1)*16);
    line[eye]->SetLineColor(2);
    line[eye]->SetLineWidth(1);
    line[eye]->Draw("same");
  }

  TH1D *hDrops = new TH1D("hDrops","hDrops",50,0,100);
  //TH1D *hDrops = new TH1D("hDrops","hDrops",25,0,1);

  cout<<"Overall Statistics:"<<endl;
  cout<<total_run_counter<<" Total Runs Analyzed"<<endl;
  cout<<total_drop_counter<<" Total Runs had Data Drops.  Fraction is: "<<1.0*total_drop_counter/(1.0*total_run_counter)<<endl;
    
  for(int eye=0; eye<group_counter; eye++) {
    cout<<"Group "<<eye<<": "<<endl;
    cout<<run_counter[eye]<<" Runs Analyzed"<<endl;
    cout<<drop_counter[eye]<<" Runs had Data Drops.  Fraction is: "<<1.0*drop_counter[eye]/(1.0*run_counter[eye])<<endl;
    //  hDrops->Fill(1.0*drop_counter[eye]/(1.0*run_counter[eye]));
  }

  for(int eye=0; eye<time_between_drops.size(); eye++) {
    hDrops->Fill(time_between_drops[eye]);
  }

  TCanvas *c2 = new TCanvas("cDrops","cDrops");
  hDrops->Draw();

  TCanvas *c3 = new TCanvas("cBM","cBM");
  c3->Divide(6,4,0.001,0.001);
 

  cout<<average_T0<<" "<<average_Li6<<" "<<average_U235<<endl;

  hLi6_Gated_1to100eV->Scale(42);
  hU235_Gated_1to100eV->Scale(42);
  hHe3_Gated_1to100eV->Scale(42);

  for(int eye=0; eye<hT0->GetNbinsX(); eye++) {
    if(hT0->GetBinContent(eye+1) < fraction * average_T0) {
      exclude_T0.push_back(hT0->GetXaxis()->GetBinCenter(eye+1));
    }
    if(hT0->GetBinContent(eye+1)>0) {
      hLi6_per_T0->SetBinContent(eye+1,(hLi6->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));
      hU235_per_T0->SetBinContent(eye+1,(hU235->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));
      hDANCE_per_T0->SetBinContent(eye+1,(hDANCE->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));
      hHe3_per_T0->SetBinContent(eye+1,(hHe3->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));
      hCharge_per_T0->SetBinContent(eye+1,(hCharge->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));
    }
    if(hCharge->GetBinContent(eye+1)>0) {
      hLi6_per_Charge->SetBinContent(eye+1,(hLi6->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));
      hU235_per_Charge->SetBinContent(eye+1,(hU235->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));
      hDANCE_per_Charge->SetBinContent(eye+1,(hDANCE->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));
      hHe3_per_Charge->SetBinContent(eye+1,(hHe3->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));

      hLi6_Gated_1to100eV_per_Charge->SetBinContent(eye+1,(hLi6_Gated_1to100eV->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));
      hU235_Gated_1to100eV_per_Charge->SetBinContent(eye+1,(hU235_Gated_1to100eV->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));
      hHe3_Gated_1to100eV_per_Charge->SetBinContent(eye+1,(hHe3_Gated_1to100eV->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));

      cout<<hLi6_per_Charge->GetXaxis()->GetBinCenter(eye+1)<<"  "<<hLi6_per_Charge->GetBinContent(eye+1)/ hLi6_Gated_1to100eV_per_Charge->GetBinContent(eye+1)<<endl;

      double mult2to6counts=0;
      for(int jay=3; jay<8; jay++) {
	mult2to6counts += hDANCE_Mult_vs_Run->GetBinContent(eye+1,jay);
      }
      hDANCE_Mult_2to6_per_Charge->SetBinContent(eye+1,(mult2to6counts/hCharge->GetBinContent(eye+1)));

      hLi6_per_Charge_Proj->Fill((hLi6->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));
      hU235_per_Charge_Proj->Fill((hU235->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));
      hDANCE_per_Charge_Proj->Fill((hDANCE->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));
      hHe3_per_Charge_Proj->Fill((hHe3->GetBinContent(eye+1)/hCharge->GetBinContent(eye+1)));

      hDANCE_Mult_2to6_per_Charge_Proj->Fill(mult2to6counts/hCharge->GetBinContent(eye+1));
      
      

    }
    if(hU235->GetBinContent(eye+1)>0) {
      hLi6_per_U235->SetBinContent(eye+1,(hLi6->GetBinContent(eye+1)/hU235->GetBinContent(eye+1)));
      hHe3_per_U235->SetBinContent(eye+1,(hHe3->GetBinContent(eye+1)/hU235->GetBinContent(eye+1)));

    }
  }

  c3->cd(1);
  hLi6->Draw();
  hLi6_Gated_1to100eV->SetLineColor(2);
  hLi6_Gated_1to100eV->Draw("same");
  c3->cd(2);
  hU235->Draw();
  hU235_Gated_1to100eV->SetLineColor(2);
  hU235_Gated_1to100eV->Draw("same");
  c3->cd(3);
  hHe3->Draw();
  hHe3_Gated_1to100eV->SetLineColor(2);
  hHe3_Gated_1to100eV->Draw("same");
  c3->cd(4);
  hDANCE->Draw();
  c3->cd(5);
  hLi6_per_U235->Draw();
  c3->cd(6);
  hHe3_per_U235->Draw();

  c3->cd(7);
  hLi6_per_T0->Draw();
  c3->cd(8); 
  hU235_per_T0->Draw();
  c3->cd(9);
  hHe3_per_T0->Draw();
  c3->cd(10);
  hDANCE_per_T0->Draw();
  c3->cd(11);
  hT0->Draw();
  c3->cd(12);
  hLi6_PH_vs_Run->GetYaxis()->SetRangeUser(0,30000);
  hLi6_PH_vs_Run->Draw("colz");

  c3->cd(13);
  hLi6_per_Charge->Draw();
  hLi6_Gated_1to100eV_per_Charge->SetLineColor(2);
  hLi6_Gated_1to100eV_per_Charge->Draw("same");
  c3->cd(14); 
  hU235_per_Charge->Draw();
  hU235_Gated_1to100eV_per_Charge->SetLineColor(2);
  hU235_Gated_1to100eV_per_Charge->Draw("same");
  c3->cd(15);
  hHe3_per_Charge->Draw();
  hHe3_Gated_1to100eV_per_Charge->SetLineColor(2);
  hHe3_Gated_1to100eV_per_Charge->Draw("same");
  c3->cd(16);
  hDANCE_per_Charge->Draw();
  c3->cd(17);
  hCharge->Draw();
  c3->cd(18);
  hU235_PH_vs_Run->GetYaxis()->SetRangeUser(0,30000);
  hU235_PH_vs_Run->Draw("colz");

  c3->cd(19);
  hLi6_per_Charge_Proj->Draw();
  c3->cd(20); 
  hU235_per_Charge_Proj->Draw();
  c3->cd(21);
  hHe3_per_Charge_Proj->Draw();
  c3->cd(22);
  hDANCE_per_Charge_Proj->Draw();
  c3->cd(23);
  hCharge_per_T0->Draw();
  c3->cd(24);
  hHe3_PH_vs_Run->GetYaxis()->SetRangeUser(0,30000);
  hHe3_PH_vs_Run->Draw("colz");

  c3->cd(23);
  hDANCE_Mult_vs_Run->Draw("colz");
  c3->cd(22);
  hDANCE_Mult_2to6_per_Charge->Draw();

  TCanvas *c4 = new TCanvas("c4","c4");
  hDANCE_Mult_2to6_per_Charge_Proj->Draw();

  
  cout<<"Exclude these for lack of data"<<endl;
  for(int eye=0; eye<exclude_T0.size(); eye++) {
    cout<<exclude_T0[eye]<<"\t1"<<endl;
  }

  cout<<"Exclude these for data drops"<<endl;
  for(int eye=0; eye<exclude_DD.size(); eye++) {
    cout<<exclude_DD[eye]<<"\t2"<<endl;
  }
  
  fout->Write();

}
