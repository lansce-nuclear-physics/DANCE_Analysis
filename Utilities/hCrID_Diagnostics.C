#include <iostream>
#include <sstream>
#include <fstream>

void hCrID_Diagnostics() {

  gStyle->SetOptStat(0);
  
  //Location of summed rootfiles
  stringstream fpath;
  fpath.str();
  fpath << "/home/cprokop/DANCE_Analysis/stage1_root/";
  
  int CW = 10;
  int CBT = 2000;
  int DEBT = 0;

  double fraction = 0.5;  //If the value is less than the average by more than this factor exclude

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

  //  int start_run = 60409;
    //  int end_run = 60451;

  //   int start_run = 63590;
  //   int end_run = 63649;

  //   int start_run = 68104;
  //   int end_run = 68126;
  
  ////////////////////////////
  
  //These are the Cu65 runs
  
  //   int start_run = 60834;
  //   int end_run = 61251;

  //   int start_run = 63015;
  //   int end_run = 63589;
   
  //  int start_run = 64138;
  //  int end_run = 64692;

  //   int start_run = 64989;
  //   int end_run = 65265;

  ////////////////////////


  //These are the Pb208 runs
  
  //  int start_run = 60622;
  //  int end_run = 60811;

  //  int start_run = 62695;
  //  int end_run = 62735;
   
  //  int start_run = 63650;
  //  int end_run = 64136;

  //  int start_run = 64834;
  //  int end_run = 64987;

  //  int start_run = 68040;
  //  int end_run = 68103;

    ////////////////////////


  //These are the Cu63 runs
  
  //  int start_run = 64693;
  //  int end_run = 64833;
  
  ////////////////////////













  
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
  TH1D *hID;
  
  double Counts[N_Runs][162];
  int Run_Num[N_Runs];

  int counter=0;
  
  TH2D *hCounts = new TH2D("hCounts","hCounts",N_Runs,start_run,end_run+1,162,0,162);
  TH1D *hT0 = new TH1D("T0","T0",N_Runs,start_run,end_run+1);
  TH1D *hLi6 = new TH1D("Li6","Li6",N_Runs,start_run,end_run+1);
  TH1D *hU235 = new TH1D("U235","U235",N_Runs,start_run, end_run+1);
  TH1D *hDANCE = new TH1D("DANCE","DANCE",N_Runs,start_run,end_run+1);
  TH1D *hHe3 = new TH1D("He3","He3",N_Runs,start_run, end_run+1);

  TH1D *hDANCE_per_T0 = new TH1D("DANCE_per_T0","DANCE_per_T0",N_Runs,start_run,end_run+1);

  TH1D *hLi6_per_T0 = new TH1D("Li6_per_T0","Li6_per_T0",N_Runs,start_run,end_run+1);
   TH1D *hU235_per_T0 = new TH1D("U235_per_T0","U235_per_T0",N_Runs,start_run,end_run+1);
 TH1D *hLi6_per_U235 = new TH1D("Li6_per_U235","Li6_per_U235",N_Runs,start_run, end_run+1);
   TH1D *hHe3_per_T0 = new TH1D("He3_per_T0","He3_per_T0",N_Runs,start_run,end_run+1);

  ifstream testrun;
  bool fileexists=false;

  for(int eye=0; eye<50; eye++) {
    drop_counter[eye] = 0;
    run_counter[eye] = 0;
    for(int jay=0; jay<162; jay++) {
      ch_drop_counter[eye][jay] = 0;
    }
  }

  for(int eye=start_run; eye<=end_run; eye++) {

    stringstream teststring;
    teststring.str();
    teststring << fpath.str() << "Stage1_Histograms_Run_"<< eye<<"_"<<CW<<"ns_CW_"<<CBT<<"ns_CBT_"<<DEBT<<"ns_DEBT.root";
    cout<<eye<<"  "<<teststring.str()<<endl;
    
    testrun.open(teststring.str().c_str());
    fileexists=false;
    
    if(testrun.is_open()) {
      fileexists=true;
      cout<<"Run: "<<eye<<" exists"<<endl;
      testrun.close();
    }
    
   
    if(fileexists) {  
      
      fin[counter] = new TFile(Form("%s/Stage1_Histograms_Run_%d_%dns_CW_%dns_CBT_%dns_DEBT.root",fpath.str().c_str(),eye,CW,CBT,DEBT));  
      Run_Num[counter] = eye;
      hID = (TH1D*)fin[counter]->FindObjectAny("hID");
      
      cout<<hID<<endl;

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
  c3->Divide(5,2);
  
  c3->cd(1);
  hLi6->Draw();

  c3->cd(2);
  hU235->Draw();
  
  c3->cd(5);
  hT0->Draw();
  
  c3->cd(3);
  hDANCE->Draw();

  c3->cd(4);
  hHe3->Draw();

  cout<<average_T0<<" "<<average_Li6<<" "<<average_U235<<endl;

  for(int eye=0; eye<hT0->GetNbinsX(); eye++) {
    if(hT0->GetBinContent(eye+1) < fraction * average_T0) {
      exclude_T0.push_back(hT0->GetXaxis()->GetBinCenter(eye+1));
    }
    if(hT0->GetBinContent(eye+1)>0) {
      hLi6_per_T0->SetBinContent(eye+1,(hLi6->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));
      hU235_per_T0->SetBinContent(eye+1,(hU235->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));
      hDANCE_per_T0->SetBinContent(eye+1,(hDANCE->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));
      hHe3_per_T0->SetBinContent(eye+1,(hHe3->GetBinContent(eye+1)/hT0->GetBinContent(eye+1)));

    }
    if(hU235->GetBinContent(eye+1)>0) {
      hLi6_per_U235->SetBinContent(eye+1,(hLi6->GetBinContent(eye+1)/hU235->GetBinContent(eye+1)));
    }
  }

    c3->cd(6);
    hLi6_per_T0->Draw();

    c3->cd(7); 
    hU235_per_T0->Draw();

    c3->cd(8);
    hDANCE_per_T0->Draw();

    c3->cd(9);
    hHe3_per_T0->Draw();

    c3->cd(10);
    hLi6_per_U235->Draw();

    
  cout<<"Exclude these for lack of data"<<endl;
  for(int eye=0; eye<exclude_T0.size(); eye++) {
    cout<<exclude_T0[eye]<<endl;
  }

  cout<<"Exclude these for data drops"<<endl;
  for(int eye=0; eye<exclude_DD.size(); eye++) {
    cout<<exclude_DD[eye]<<endl;
  }
  
}
