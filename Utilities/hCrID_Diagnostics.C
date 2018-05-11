#include <iostream>
#include <sstream>
#include <fstream>

void hCrID_Diagnostics() {

  gStyle->SetOptStat(0);
  
  //Location of summed rootfiles
  stringstream fpath;
  fpath.str();
  fpath << "/home/cprokop/CJP/DANCE_Analysis/stage1_root/";
  
  int CW = 10;
  int CBT = 2000;
  int DEBT = 0;

  //These are the Au197 runs

  //  int start_run = 60572;
  //  int end_run = 60584;

  //  int start_run = 60409;
  //  int end_run = 60451;

  //  int start_run = 63590;
  //  int end_run = 63649;

  //  int start_run = 68104;
  //  int end_run = 68126;
  
  ////////////////////////////
  
  //These are the Cu65 runs
  
  //  int start_run = 60834;
  //  int end_run = 61251;

  //  int start_run = 63015;
  //  int end_run = 63589;
   
    int start_run = 64138;
    int end_run = 64692;

  //  int start_run = 64989;
  //  int end_run = 65265;

  ////////////////////////
    
  int last_drop = start_run; 
  vector<int> time_between_drops;
  int drop_counter[50];
  int total_drop_counter = 0;
  int ch_drop_counter[50][162];
  int run_counter[50];
  int total_run_counter =0;
  int group_counter;
    
  int nruns_per_group = 50;
  

  const int N_Runs = end_run-start_run+1;
  
  cout<<N_Runs<<endl;

  TFile *fin[N_Runs];
  TH1D *hID;
  
  double Counts[N_Runs][162];
  int Run_Num[N_Runs];

  int counter=0;
  
  TH2D *hCounts = new TH2D("hCounts","hCounts",N_Runs,start_run,end_run+1,162,0,162);
  TH1D *hT0s = new TH1D("hT0s","hT0s",N_Runs,start_run,end_run);
  TH1D *hLi6 = new TH1D("Li6","Li6",N_Runs,start_run,end_run);
  TH1D *hU235 = new TH1D("U235","U235",N_Runs,start_run, end_run);

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
	hCounts->Fill(eye,jay,Counts[counter][jay]);

	if(Counts[counter][jay] == 0 && jay != 76 && jay !=86) {
	  ch_drop_counter[group_counter][jay] ++;
	  missing_channels = true;
	}
      }
      
      //T0's and beam monitors
      hT0s->Fill(eye,hID->GetBinContent(201));  // ID 200 bin 201
      hLi6->Fill(eye,hID->GetBinContent(245));  // ID 244 bin 245
      hU235->Fill(eye,hID->GetBinContent(244)); // ID 243 bin 244

      if(missing_channels) {
	time_between_drops.push_back(eye-last_drop);
	last_drop = eye;
	drop_counter[group_counter]++;
	total_drop_counter ++;
      }

      if(run_counter[group_counter] >=nruns_per_group) {
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
  c3->Divide(2,2);
  
  c3->cd(1);
  hLi6->Draw();

  c3->cd(2);
  hU235->Draw();
  
  c3->cd(3);
  hT0s->Draw();
  
}
