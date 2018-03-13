#include <iostream>
#include <sstream>
#include <fstream>

void hCrID_Diagnostics() {
  
  //Location of summed rootfiles
  stringstream fpath;
  fpath.str();
  fpath << "/home/cprokop/CJP/DANCE_Analysis/stage0_root/";
  
  // int start_run = 1550;
  //  int end_run = 1559;

  //  int start_run = 3010;
  // int end_run = 3019;
 
  // int start_run = 5010;
  //  int end_run = 6020;

  int start_run = 7000;
  int end_run = 7400;
  //    int start_run = 2000;
  //  int end_run = 2009;
   
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
  
  TH2D *hCounts = new TH2D("hCounts","hCounts",N_Runs,start_run,end_run+1,162,0,162);

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
    counter++;
  }
  
  hCounts->Draw("colz");

  TLine *line[10];

  for(int eye=0; eye<10; eye++) {
    line[eye] = new TLine(start_run,(eye+1)*16,end_run+1,(eye+1)*16);
    line[eye]->SetLineColor(2);
    line[eye]->SetLineWidth(1);
    line[eye]->Draw("same");
  }


}
