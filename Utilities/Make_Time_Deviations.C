#include <iostream>
#include <sstream>
#include <fstream>

#include "TStyle.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"

//TMatrix Things
int reftoindex1[200];
int reftoindex2[200];
int index1[200];
int index2[200];

int totalindex; //Number of TMatrix Pairs read in

int Read_TMatrix() {

  cout <<"Reading "<<"/home/cprokop/CJP/DANCE_Analysis/Config/TMatrix.txt" << endl;
  
  ifstream timemat;
  timemat.open("/home/cprokop/CJP/DANCE_Analysis/Config/TMatrix.txt");
    
  for(int i=0;i<200;i++){
    reftoindex1[i]=-1;
    reftoindex2[i]=-1;
  }
  int counter=0;

  if(timemat.is_open()) {
    while(!timemat.eof()) {
      
      //aaa and bbb are Detector IDs for neighboring crystals
      int aaa,bbb;
      timemat >> aaa >> bbb;

      index1[counter] = aaa;   
      index2[counter] = bbb;
      
      if(timemat.eof()) break;

      //for each pair of neighboring crystals there is a refereence to them.  
      //This way if there is a missing crystal there isnt a gap in the array.  
      reftoindex1[aaa]=counter;   //left entry of pair
      reftoindex2[bbb]=counter;   //right entry of pair
      
      // cout<<counter<<"  "<<aaa<<"  "<<bbb<<"  "<<reftoindex1[aaa]<<"  "<<reftoindex2[bbb]<<endl;
      //number of pairs read in
      counter++;
    }    
    timemat.close();
    cout<<"Read in TMatrix"<<endl;
  }
  else {
    cout<<"Couldn't open /home/cprokop/CJP/DANCE_Analysis/Config/TMatrix.txt"<<endl;
  }
  
  return counter;
}


void Make_Time_Deviations(int start_run, int end_run) {
  
  //Make a histogram
  TH2D *hTimeDeviations = new TH2D("TimeDeviations","TimeDeviations",end_run-start_run+1,start_run, end_run+1, 162,0,162);

  totalindex = Read_TMatrix();

  gStyle->SetOptStat(0);
  
  //Location of stage0 rootfiles
  stringstream fpath;
  fpath.str();
  // fpath << "/home/cprokop/CJP/DANCE_Analysis/stage0_root/";
  fpath << "/home/cprokop/CJP/DANCE_Analysis/stage0_root_automated/";

  //stage0
  int CW = 500;
  int CBT = 0;
  int DEBT = 0;
   
  const int N_Runs = end_run-start_run+1;
  cout<<N_Runs<<endl;

  TFile *fin[N_Runs];  
  TH2D *hTimeDev[N_Runs];

  int counter=0;

  for(int eye=start_run; eye<=end_run; eye++) {
    
    stringstream teststring;
    teststring.str();
    teststring << fpath.str() << "Stage0_Histograms_Run_"<< eye<<"_"<<CW<<"ns_CW_"<<CBT<<"ns_CBT_"<<DEBT<<"ns_DEBT.root";
    //cout<<eye<<"  "<<teststring.str()<<endl;
  
    ifstream testrun;
    testrun.open(teststring.str().c_str());
    bool fileexists=false;
  
    if(testrun.is_open()) {
      fileexists=true;
      // cout<<"Run: "<<eye<<" stage0 exists"<<endl;
      testrun.close();
    }
  
    if(fileexists) {  
    
      fin[counter] = new TFile(Form("%s/Stage0_Histograms_Run_%d_%dns_CW_%dns_CBT_%dns_DEBT.root",fpath.str().c_str(),eye,CW,CBT,DEBT)); 
      // cout<<fin[counter]<<endl;
      hTimeDev[counter] = (TH2D*)fin[counter]->Get("TimeDev");
      //  cout<<hTimeDev[counter];

      //Stringstream for the outpout file name of td_out
      stringstream outfilename;
      outfilename.str();
      outfilename << "/home/cprokop/CJP/TimeDeviations" << "/TimeDeviations_Run_" << eye << ".txt";
    
      //initialize the time deviation to 0.  Its cumulative...
      double time_deviation=0;
    
      //open the time deviations text file
      ofstream td_out;
      td_out.open(outfilename.str().c_str());
  
      //First detector has no offest
      td_out <<"0   \t 0 \n";
    
      // cout<<"Total Indices: "<<totalindex<<endl;
    
      //Fill in the rest
      for(int jay=0; jay<totalindex; jay++) {
	//	cout<<index1[jay]<<"  "<<index2[jay]<<endl;
      
	//Set the range of the Time Deviations 2D histogram to the proper bin 
	hTimeDev[counter]->GetYaxis()->SetRangeUser(index1[jay],index1[jay]+1);
	hTimeDev[counter]->GetXaxis()->SetRangeUser(-500,500);
      
	//Iteratively Close in on the proper range 
	for(int kay=3; kay<103; kay++) {
	  hTimeDev[counter]->GetXaxis()->SetRangeUser(hTimeDev[counter]->GetMean()-(103.0-kay), hTimeDev[counter]->GetMean()+(103.0-kay));
	  if(index2[jay]==16) {
	    cout<<"counter: "<<counter<<"  kay: "<<kay<<" TDev: "<< hTimeDev[counter]->GetMean()<<endl;
	  }
	}
	/*
	  hTimeDev[counter]->GetXaxis()->SetRangeUser(hTimeDev[counter]->GetMean()-50.0, hTimeDev[counter]->GetMean()+50.0);
	  
	hTimeDev[counter]->GetXaxis()->SetRangeUser(hTimeDev[counter]->GetMean()-10.0, hTimeDev[counter]->GetMean()+10.0);
	hTimeDev[counter]->GetXaxis()->SetRangeUser(hTimeDev[counter]->GetMean()-5.0, hTimeDev[counter]->GetMean()+5.0);
	hTimeDev[counter]->GetXaxis()->SetRangeUser(hTimeDev[counter]->GetMean()-5.0, hTimeDev[counter]->GetMean()+5.0);
	hTimeDev[counter]->GetXaxis()->SetRangeUser(hTimeDev[counter]->GetMean()-5.0, hTimeDev[counter]->GetMean()+5.0);
	hTimeDev[counter]->GetXaxis()->SetRangeUser(hTimeDev[counter]->GetMean()-4.0, hTimeDev[counter]->GetMean()+4.0);
	hTimeDev[counter]->GetXaxis()->SetRangeUser(hTimeDev[counter]->GetMean()-4.0, hTimeDev[counter]->GetMean()+4.0);
*/
	//	cout<<hTimeDev[counter]->GetMean()<<"  "<<time_deviation<<endl;
	time_deviation += hTimeDev[counter]->GetMean();
	td_out <<index2[jay]<<"   \t"<<time_deviation<<"\n";
	cout<<index2[jay]<<"  "<< time_deviation<<endl;

	hTimeDeviations->Fill(eye, index2[jay], time_deviation);
      }
    
      cout<<"Made Time Deviations for Run: "<<eye<<endl;
    
      counter++;
    }
  }
  
  hTimeDeviations->Draw("colz");
  
}
