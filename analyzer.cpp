#include "analyzer.h"
#include <iostream>
#include <fstream>


#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>

using namespace std;


/* HISTOGRAMS */
TH1D *hID;  //ID's present in datastream
TH2D *hCoinCAEN;  //coincidence matrix

//Time Diagnostics
TH1D *hEventLength;  //length in ns of the event (diagnostic)
TH2D *hTimeBetweenCrystals;  //time between subsequent hits of the same crystal

//Time Deviations
TH2D *hTimeDev_Rel0;  //Time deviations of all crystals to crystal 0
TH2D *hTimeDev;  //Time deviations relative to adjacent crystals

//Energy Histograms
TH2D *ADC_calib;  //2D PSD Plot 


//Energies


/* VARIABLES */


double last_timestamp[2000];
double current_timestamp[2000];

double last_t0_timestamp;
double current_t0_timestamp;


//TMatrix Things
int reftoindex1[200];
int reftoindex2[200];
int index1[200];
int index2[200];

int totalindex; //Number of TMatrix Pairs read in


//misc
int crap_counter=0;


int Read_TMatrix() {

  cout << "Reading Config/TMatrix.txt"  << endl;
  
  ifstream timemat;
  timemat.open("./Config/TMatrix.txt");
    
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
  }
  else {
    cout<<"Couldn't open Config/TMatrix.txt"<<endl;
  }
  
  return counter;
}


int Create_Analyzer_Histograms() {
  
  cout<<"Creating Histograms"<<endl;

  //Make Histograms
  hCoinCAEN = new TH2D("CoinCAEN","CoinCAEN",162,0,162,162,0,162);  //coincidence matrix
  hID = new TH1D("hID","hID",2000,0,2000);
  
  //Time Deviations
  hTimeDev_Rel0 = new TH2D("TimeDev_Rel0","TimeDev_Rel0",10000,-500,500,162,0,162);  //Time deviations relative to crystal 0
  hTimeDev = new TH2D("TimeDev","TimeDev",10000,-500,500,162,0,162);  //Time deviations

  //Diagnostics
  hEventLength = new TH1D("EventLength","EventLength",10000,0,10000);
  hTimeBetweenCrystals = new TH2D("TimeBetweenCrystals","TimeBetweenCrystals",10000,0,10000,162,0,162);
  
  //Energy Histograms
  ADC_calib=new TH2D("ADC_calib","ADC_calib",800,0.,16.,1000,0.,10.);

  // for(int eye=0; eye<162; eye++) {
    //cout<<"eye: "<<eye<<"  reftoindex1: "<< reftoindex1[eye] <<"  reftoindex2: "<< reftoindex2[eye] <<"     "<< index1[reftoindex1[eye]]<<"  "<<index2[reftoindex1[eye]]<<endl;
  //   hTimeDev[eye]= new TH1D(Form("hTimeDev_%d",eye),Form("Crystal_%d_RelativeTo_Crystal_%d",index2[reftoindex1[eye]],index1[reftoindex1[eye]]),10000,-500,500);
  // }

  return 0;

}


int Write_Analyzer_Histograms(TFile *fout) {
  
  cout<<"Writing Histograms"<<endl;
  
  fout->cd();
  hID->Write();
  hCoinCAEN->Write();
  hTimeDev_Rel0->Write();
  hTimeDev->Write();

  hEventLength->Write();
  hTimeBetweenCrystals->Write();
  
  ADC_calib->Write();

  /*
  TDirectory *cdTimeDev = fout->mkdir("TimeDeviations");
  cdTimeDev->cd();
 

 for(int eye=0; eye<162; eye++) {
    hTimeDev[eye]->Write();
  }

  double total_dev=0;
  for(int eye=0; eye<162; eye++) {
    hTimeDev[eye]->GetXaxis()->SetRangeUser(-40,40);
    total_dev += hTimeDev[eye]->GetMean();
    cout<<eye<<"  "<<hTimeDev[eye]->GetMean()<<"  "<<total_dev<<endl;
  }
  */
 
  return 0;
}


int Initialize_Analyzer() {

  cout<<"Initializing Analyzer"<<endl;

  totalindex = Read_TMatrix();
  Create_Analyzer_Histograms();
  
  for(int eye=0; eye<162; eye++) {
    last_timestamp[eye]=0;
    current_timestamp[eye]=0;

    last_t0_timestamp=0;
    current_t0_timestamp=0;
  }
  
  return 0;
}


int Analyze_Data(std::vector<DEVT_BANK_wWF> eventvector) {

  //Fill Event Length (mult 1 events just give 0 so require mult>1 to fill this)
  if(eventvector.size() > 1) {
    hEventLength->Fill(eventvector[eventvector.size()-1].TOF-eventvector[0].TOF);
  }

  //Loop over event 
  for(int eye=0; eye<eventvector.size(); eye++) {
    
    //Fill ID histogram
    hID->Fill(eventvector[eye].ID);

    //Fill ADC Calib
    ADC_calib->Fill(eventvector[eye].Eslow, eventvector[eye].Efast,1);

    

    int id_eye=eventvector[eye].ID;
    
    //Place the current time in the 
    current_timestamp[id_eye] = eventvector[eye].TOF;

    //Fill time between crystal hits
    hTimeBetweenCrystals->Fill((current_timestamp[id_eye]-last_timestamp[id_eye]),id_eye,1);
    
    //this is a dance crystal
    if(id_eye<162) {
      
      //Coincidences
      if(eventvector.size() > 1 && eye < (eventvector.size()-1)) {
	
	//start with the next one
	for(int jay=eye+1; jay<eventvector.size(); jay++) {
	  
	  int id_jay = eventvector[jay].ID;
	  
	  if(id_jay<162) {

	    //time difference between crystal jay and eye
	    double ddT = eventvector[jay].TOF - eventvector[eye].TOF;
	    
	    //Fill the coincidence matrix
	    hCoinCAEN->Fill(id_eye,id_jay,1);
	    hCoinCAEN->Fill(id_jay,id_eye,1);
	    
	    //Look at time deviations
	    //Check relative to 0 first
	    if(id_eye==0) { 
	      //  cout<<id_eye<<"  "<<current_timestamp[id_eye]<<"  "<<last_timestamp[0]<<" diff: "<<current_timestamp[id_eye]-last_timestamp[0]<<endl;
	      
	      //relative to crystal 0;
	      hTimeDev_Rel0->Fill(-ddT,id_jay,1);
	    } 
	    if(id_jay==0) {
	      hTimeDev_Rel0->Fill(ddT,id_eye,1);
	    }
	    
	    //id11 is the index of where this crystal is in TMatrix if it is the first one listed
	    int id11=reftoindex1[eventvector[eye].ID];
	    //same but for if it is the second one listed
	    int id12=reftoindex2[eventvector[eye].ID];
	    
	    //same references to the crystals in the tmatrix as before but for the second crystal 
	    int id21=reftoindex1[eventvector[jay].ID];
	    //same but for if it is the second one listed
	    int id22=reftoindex2[eventvector[jay].ID];
	    
	    //if the leftmost crystal referernce is >=0 
	    //if the left and right are neighbors (reftocrystals are equal) 
	    //eye on left and jay on right
	    if(id11>=0 && id11==id22){
	      //   cout<<"first: "<<id_eye<<" "<<id_jay<<"   "<<id11<<"  "<<id12<<"  "<<id21<<"  "<<id22<<"   "<<index1[id11]<<endl;
	      // hTimeDev[index1[id11]]->Fill(-ddT);
	      hTimeDev->Fill(-ddT,index1[id11]);
	    }
	    //if the leftmost crystal referernce is >=0 
	    //if the left and right are neighbors (reftocrystals are equal)
	    //eye on right and jay on left 
	    if(id12>=0 && id12==id21){
	      //   cout<<"second: "<<id_eye<<" "<<id_jay<<"   "<<id11<<"  "<<id12<<"  "<<id21<<"  "<<id22<<"   "<<index1[id21]<<endl;
	      // hTimeDev[index1[id21]]->Fill(ddT);
	      hTimeDev->Fill(ddT,index1[id21]);
	    }	  
	  }
	}
      }  //Done with coincidences 
      
      
      //Calibrate the BaF2 detectors
      
      





      
      
    
    }
    
    //this is t0
    if(id_eye==200) {
      

      
    }
   

    //Once done the time is the last timestamp
    last_timestamp[id_eye] = current_timestamp[id_eye];
  }







  if(eventvector.size()>10) {



    // cout<<"Analysis: "<<eventvector.size()<<endl;
  }
  return 0;
}
