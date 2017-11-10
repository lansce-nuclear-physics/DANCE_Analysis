#include "analyzer.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "TRandom.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCutG.h"

using namespace std;

//output binary file
ofstream outputbinfile;

//Structure to write data
DEVT_OUT devt_out;



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

//Alpha Spectra
TH2D *hAlpha;
TH2D *hAlphaCalib;


//3D Histograms

TH3F *En_Esum_Mcl;
TH3F *En_Esum_Mcr;

TH3F *En_Eg_Mcl; // this should have a Qgate on it
TH3F *En_Eg_Mcr; // this should have a Qgate on it

TH3F *Esum_Eg_Mcl; // Eg is Ecluster here
TH3F *Esum_Eg_Mcr; // Eg is Ecrystal here

//Keep track of how many entries processed (Diagnostics)
//uint64_t entry_counter=0;

//Cuts

TCutG *Gamma_Gate;
TCutG *Alpha_Gate;


//Dance Event
DANCE_Event devent;

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

//DMatrix Things
int DetMat1[167];		      
int DetMat2[167];
int DetMat3[167];
int DetMat4[167];
int DetMat5[167];
int DetMat6[167];
int DetMat7[167];

//energy calibrations
double slow_offset[200];
double slow_slope[200];
double slow_quad[200];
double fast_slope[200];
double fast_offset[200];


//This function makes the output binary file for stage 0 or 1 filled with time-ordered devt_bank structures
int Make_Output_Binfile(int RunNumber) {
  
  stringstream outfilename;
  outfilename.str();
  
  //stage0 
  if(READ_BINARY==0) {
    outfilename << STAGE0_BIN; 
    outfilename <<"/stage0_run_";
  }
  //stage1
  if(READ_BINARY==1) {
    outfilename << STAGE1_BIN;
    outfilename <<"/stage0_run_";
  }
  outfilename << RunNumber << ".bin";
  
  outputbinfile.open(outfilename.str().c_str(), ios::out | ios::binary);
  
  if(outputbinfile.is_open()) {
    cout<<"Succesfully created and opened output binary file: "<<outfilename.str()<<endl;
  }
  else {
    cout<<"ERROR (Unpacker): Failed to create output binary file: "<<outfilename.str()<<endl;
  }
}

int Read_Energy_Calibrations(int RunNumber) {
  
  int id=0;
  double temp[5] = {0,0,0,0,0};

  for(int eye=0; eye<200; eye++) {
    slow_offset[eye]=0;
    slow_slope[eye]=0;
    slow_quad[eye]=0;
    fast_offset[eye]=0;
    fast_slope[eye]=0;  
  }
  
  int retval=0;
  
  cout<<"Reading Energy Calibrations for Run "<<RunNumber<<endl;
  
  stringstream cal_name;
  cal_name.str();
  cal_name << "./Calibrations/";
  cal_name << "param_out_" << RunNumber << ".txt";
  
  ifstream encal;
  encal.open(cal_name.str().c_str());
  
  if(encal.is_open()) {
    cout<<"File: "<<cal_name.str()<<" found"<<endl;
    while(!encal.eof()) {
      encal>>id>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4];
      slow_offset[id]=temp[0];
      slow_slope[id]=temp[1];
      slow_quad[id]=temp[2];
      fast_offset[id]=temp[3];
      fast_slope[id]=temp[4];
      
      //  cout<<id<<"  "<<slow_offset[id]<<"  "<<slow_slope[id]<<"  "<<slow_quad[id]<<"  "<<fast_offset[id]<<"  "<<fast_slope[id]<<endl;

      if(encal.eof()) {
	break;
      }
    }   
  }
  else {
    cout<<"File: "<<cal_name.str()<<" NOT found..."<<endl;
    cout<<"Looking for calib_ideal.dat"<<endl;
    
    ifstream idealcal;
    idealcal.open("./Calibrations/calib_ideal.dat");
    
    if(idealcal.is_open()) {
      while(!idealcal.eof()) {
	idealcal>>id>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4];
	slow_offset[id]=temp[0];
	slow_slope[id]=temp[1];
	slow_quad[id]=temp[2];
	fast_offset[id]=temp[3];
	fast_slope[id]=temp[4];
	 
	//	cout<<id<<"  "<<slow_offset[id]<<"  "<<slow_slope[id]<<"  "<<slow_quad[id]<<"  "<<fast_offset[id]<<"  "<<fast_slope[id]<<endl;
	 
	if(idealcal.eof()) {
	  break;
	}
      }
      retval=1;
    }
    else {
      cout<<"Can't open either calibration file..."<<endl;
      retval=-1;
      
    }
  }
  return retval;
}

int Read_PI_Gates() {

  char name[200];
  int Nalphacut=0;
  int Ngammacut=0;
  double x_alphacut[200];
  double y_alphacut[200];
  double x_gammacut[200];
  double y_gammacut[200];

  sprintf(name,"Gates/%s",ALPHAGATE);
  ifstream alphacutin(name);
  if(alphacutin.is_open()) {
    cout << "Analyzer: Reading AlphaCut : " << name << endl;
    while(!alphacutin.eof()){
      alphacutin  >> x_alphacut[Nalphacut] >> y_alphacut[Nalphacut];
      Nalphacut++;
      if(alphacutin.eof()) break;
    }
    cout<<Nalphacut<<endl;
    alphacutin.close();
  }
  else {
    cout << " Alpha PI cut " << name << "  file not found." << endl;
  }
  
  
  sprintf(name,"Gates/%s",GAMMAGATE);
  ifstream gammacutin(name);
  if(gammacutin.is_open()) {
    cout << "Analyzer:  Reading GammaCut : " << name << endl;
    while(!gammacutin.eof()){
      gammacutin  >> x_gammacut[Ngammacut] >> y_gammacut[Ngammacut];
      Ngammacut++;
      if(gammacutin.eof()) break;
    }
    cout<<Ngammacut<<endl;
    gammacutin.close();
  }
  else {
    cout << " Gamma PI cut " << name << "  file not found." << endl;
  }

  Alpha_Gate=new TCutG("Alpha_Gate",Nalphacut,x_alphacut,y_alphacut);
  // Alpha_Gate->Print();	
  Gamma_Gate=new TCutG("Gamma_Gate",Ngammacut,x_gammacut,y_gammacut);
  //  Gamma_Gate->Print();
  
  return 0;
}


int Read_DMatrix() {

  ifstream matrix_in("Config/DetectorMatrix.txt");
  
  cout <<" Analyzer:  Reading in a Detector Matrix" << endl;
  
  for(int eye=0; eye<167; eye++){
    
    matrix_in >> DetMat1[eye];
    matrix_in >> DetMat2[eye];
    matrix_in >> DetMat3[eye];
    matrix_in >> DetMat4[eye];
    matrix_in >> DetMat5[eye];
    matrix_in >> DetMat6[eye];
    matrix_in >> DetMat7[eye];
    
    // cout<<eye<<"  "<<DetMat1[eye]<<" "<<DetMat2[eye]<<" "<<DetMat3[eye]<<" "<<DetMat4[eye]<<" "<<DetMat5[eye]<<" "<<DetMat6[eye]<<" "<<DetMat7[eye]<<endl;
  }
  matrix_in.close();
  
  return 0;
}


int Read_TMatrix() {

  cout <<" Analyzer: Reading ./Config/TMatrix.txt"  << endl;
  
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
  hID = new TH1D("hID","hID",256,0,256);
  
  //Time Deviations
  hTimeDev_Rel0 = new TH2D("TimeDev_Rel0","TimeDev_Rel0",10000,-500,500,162,0,162);  //Time deviations relative to crystal 0
  hTimeDev = new TH2D("TimeDev","TimeDev",10000,-500,500,162,0,162);  //Time deviations

  //Diagnostics
  hEventLength = new TH1D("EventLength","EventLength",10000,0,10000);
  hTimeBetweenCrystals = new TH2D("TimeBetweenCrystals","TimeBetweenCrystals",10000,0,10000,162,0,162);
  
  //Energy Histograms
  ADC_calib=new TH2D("ADC_calib","ADC_calib",800,0.,16.,1000,0.,10.);

  //Alpha Histograms
  hAlpha = new TH2D("hAlpha","hAlpha",1500,0,30000,162,0,162);
  hAlphaCalib = new TH2D("hAlphaCalib","hAlphaCalib",500,0.0,5.0,162,0,162);


  //Physics Histograms
  double x[5000];
  int NEbins=0;
  
  for(double lx=log10(NeutronE_From);lx<log10(NeutronE_To);lx=lx+(1./NeutronE_BinsPerDecade)){
    x[NEbins]=pow(10,lx);
    
    //   cout << NEbins << "	" << x[NEbins] << endl;
    NEbins++;
  }
  NEbins--;
  
  int NoOfEnergyBins=GammaE_NoOfBins;
  double *EtotBins=new double[NoOfEnergyBins+1];
  double *Mbins=new double[21];
  
  for(int i=0;i<21;i++){
    Mbins[i]=0.5+1.*i;
  };

  double DEGamma=(GammaE_To-GammaE_From)/GammaE_NoOfBins;
  for(int i=0;i<NoOfEnergyBins+1;i++){
    //EtotBins[i]=i*16./128.;
    EtotBins[i]=GammaE_From+i*DEGamma;
  };
  
  En_Esum_Mcl=new TH3F("En_Etot_Mcl","En_Etot_Mcl",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
  En_Esum_Mcr=new TH3F("En_Etot_Mcr","En_Etot_Mcl",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
  
  En_Eg_Mcl=new TH3F("En_Eg_Mcl","En_Eg_Mcl gated on Q",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
  En_Eg_Mcr=new TH3F("En_Eg_Mcr","En_Eg_Mcl gated on Q",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
  
  Esum_Eg_Mcl=new TH3F("Esum_Eg_Mcl","Esum_Eg_Mcl where Eg is Ecluster",NoOfEnergyBins,EtotBins,NoOfEnergyBins,EtotBins,20,Mbins);
  Esum_Eg_Mcr=new TH3F("Esum_Eg_Mcr","Esum_Eg_Mcr where Eg is Ecrystal",NoOfEnergyBins,EtotBins,NoOfEnergyBins,EtotBins,20,Mbins);
  

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

  hAlpha->Write();
  hAlphaCalib->Write();

  Alpha_Gate->Write();
  Gamma_Gate->Write();


  En_Esum_Mcl->Write();
  En_Esum_Mcr->Write();
  
  return 0;
}


int Initialize_Analyzer() {

  cout<<"Initializing Analyzer"<<endl;
  
  //Initialize Analysis Stuff
  Read_PI_Gates();
  totalindex = Read_TMatrix();
  Read_DMatrix();
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

  double ESum=0;
  int Crystal_Mult=0;

  //Initialize DANCE Event
  devent.Crystal_mult=0;
  devent.ESum=0;
  



  //Fill Event Length (mult 1 events just give 0 so require mult>1 to fill this)
  if(eventvector.size() > 1) {
    hEventLength->Fill(eventvector[eventvector.size()-1].TOF-eventvector[0].TOF);
  }

  //Loop over event 
  for(int eye=0; eye<eventvector.size(); eye++) {

    //Fill ID histogram
    hID->Fill(eventvector[eye].ID);
   
    int id_eye=eventvector[eye].ID;

  
        
    //this is a dance crystal
    if(id_eye<162) {

      //Apply Energy Calibration
      double temp_slow = eventvector[eye].lgate + gRandom->Uniform(0,1);
      double temp_fast = eventvector[eye].sgate + gRandom->Uniform(0,1);
      
      eventvector[eye].Eslow = 0.001*(temp_slow*temp_slow*slow_quad[eventvector[eye].ID] +
				  temp_slow*slow_slope[eventvector[eye].ID] +
				  slow_offset[eventvector[eye].ID]);
      
      eventvector[eye].Efast = 0.001*(temp_fast*fast_slope[eventvector[eye].ID] +
				  fast_offset[eventvector[eye].ID]);
      
      
      if(HAVE_Threshold==1) {
	if(eventvector[eye].Eslow < Energy_Threshold) {
	  eventvector[eye].Valid=0;
	}
      }
      
      
    
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
      

      //Place the current time in the 
      current_timestamp[id_eye] = eventvector[eye].TOF;

      //Fill time between crystal hits
      hTimeBetweenCrystals->Fill((current_timestamp[id_eye]-last_timestamp[id_eye]),id_eye,1);

      //Fill ADC Calib
      ADC_calib->Fill(eventvector[eye].Eslow, eventvector[eye].Efast,1);

      int Is_Alpha = Alpha_Gate->IsInside(eventvector[eye].Eslow, eventvector[eye].Efast);
      int Is_Gamma = Gamma_Gate->IsInside(eventvector[eye].Eslow, eventvector[eye].Efast);

      if(Is_Alpha) {
	//Fill some Calibration Spectra
	hAlpha->Fill(eventvector[eye].lgate, eventvector[eye].ID,1);
	hAlphaCalib->Fill(eventvector[eye].Eslow, eventvector[eye].ID,1);
      }
    
      if(Is_Gamma) {
	//Make a DANCE Event
	devent.Crystal_ID[Crystal_Mult] = eventvector[eye].ID;  //Crystal ID
	devent.Cluster_ID[Crystal_Mult] = Crystal_Mult;  //??????
	devent.Islow[Crystal_Mult] = eventvector[eye].lgate;  //Crystal long integral
	devent.Islow[Crystal_Mult] = eventvector[eye].lgate;  //Crystal short integral
	devent.tof[Crystal_Mult] = eventvector[eye].TOF;  //time of flight
	devent.Ecrystal[Crystal_Mult] = eventvector[eye].Eslow;   //Energy if calibrated 
	devent.tof[Crystal_Mult] = eventvector[eye].TOF; //time of flight for crystal hit
	devent.ESum += eventvector[eye].Eslow; //ESum 
	devent.Crystal_mult++;
	Crystal_Mult++;
      }

      


    
    }
    
    //this is t0
    if(id_eye==200) {
      
      current_t0_timestamp = eventvector[eye].TOF;
      
      //Do some T0 diagnostics
            
      last_t0_timestamp=current_t0_timestamp;
      
    }
   
    //beam monitors
    























    //Once done the time is the last timestamp
    last_timestamp[id_eye] = current_timestamp[id_eye];

    //Write to Binary
    if(WRITE_BINARY==1 && outputbinfile.is_open()) {
      devt_out.sgate = eventvector[eye].sgate;
      devt_out.lgate = eventvector[eye].lgate;
      devt_out.TOF = eventvector[eye].TOF;
      devt_out.ID = eventvector[eye].ID;
      outputbinfile.write(reinterpret_cast<char*>(&devt_out),sizeof(DEVT_OUT));
    }    
  }


  
  //Handle various events and do some physics
  
  if(devent.Crystal_mult>0) {

    //TOF now relative to last T0
    for(int kay=0; kay<devent.Crystal_mult; kay++) {
      devent.tof[kay] -= last_t0_timestamp;
    }
    
    //Calculate the neutron energy    
    devent.En = 0.5*939.565379e6*DANCE_FlightPath*DANCE_FlightPath/((devent.tof[0]+DANCE_Delay)/1e9)/((devent.tof[0]+DANCE_Delay)/1e9)/(2.997924589e8*2.997924589e8); 
    
    //    cout<<"Analyzer: Made DANCE Event of Size: "<<devent.Crystal_mult<<" ESum: "<<devent.ESum<<" "<<DANCE_FlightPath<<" "<<DANCE_Delay<<" En: "<<devent.En<<endl;
    
    //Only one crystal...
    if(devent.Crystal_mult==1) {
      devent.Cluster_mult=1;
      devent.Ecluster[0]=devent.Ecrystal[0];
    }
    
    //Multiple Crystals... Clusterize
    if(devent.Crystal_mult > 1) {
      
      // --------------------------------------------------------------------------------------------------------
      // here we will clusterize
      // result will be in Cluster_ID[i] where i runs through all the crystals
      // minimum of Cluster_ID is 1 so don't forget to lower it by 1
      // it is somewhat different then Jan Wouters routine but should lead to same results
      // if used stand alone, set Cluster_ID[i]=i+1 for i=0,Crystal_Mult
      // Label value (after this routine) is the cluster multiplicity
      
      int jj=0;
      int Label=0;
      int intcounter;
      int Hits;
      int Internal_ID[163];

      for(int loop=0; loop<devent.Crystal_mult; loop++){
	devent.Ecluster[loop]=0.;
	
	if(loop==0 || devent.Cluster_ID[loop]>Label){
	  
	  Label++;
	  intcounter=0;
	  Hits=1;
	  Internal_ID[0]=devent.Crystal_ID[loop];
	  
	  for(int l=0;l<Hits;l++) {
	    for(jj=0;jj<devent.Crystal_mult;jj++) {
	      if(Internal_ID[l]!=devent.Crystal_ID[jj]) { // no need to check the same crystals
			
		bool IsNeighbor=false;
		
		if(DetMat1[Internal_ID[l]]==DetMat2[devent.Crystal_ID[jj]]) IsNeighbor=true;
		else if(DetMat1[Internal_ID[l]]==DetMat3[devent.Crystal_ID[jj]]) IsNeighbor=true;
		else if(DetMat1[Internal_ID[l]]==DetMat4[devent.Crystal_ID[jj]]) IsNeighbor=true;
		else if(DetMat1[Internal_ID[l]]==DetMat5[devent.Crystal_ID[jj]]) IsNeighbor=true;
		else if(DetMat1[Internal_ID[l]]==DetMat6[devent.Crystal_ID[jj]]) IsNeighbor=true;
		else if(DetMat1[Internal_ID[l]]==DetMat7[devent.Crystal_ID[jj]]) IsNeighbor=true;
		
		if(IsNeighbor) {	
		  // if neighbor is found label it and add it to internal list (Internal_ID[]) for checking
		  // further in l loop (Hits is incremented by one)
		  // this way a new participant of the cluster is being chacked agains the other crystals
		    
		  if((1.*devent.Cluster_ID[jj])>=0 && (1.*devent.Cluster_ID[jj])>(1.*Label)) {
		    // here if the crystal is already labeled we skip
		    // so that we do not have to repeat already labeled crystals
		    // more speed to it
		    
		    devent.Cluster_ID[jj]=Label; // this one is important- Adding the
		    //	crystals to the cluster
		    
		    Hits++;
		    intcounter++;
		    Internal_ID[intcounter]=devent.Crystal_ID[jj];	
		  }
		}
	      }
	    }
	  }
	  if(Hits==1) devent.Cluster_ID[loop]=Label;
	}
      }
      
      devent.Cluster_mult=Label;
      
      // Fill the cluster energy
      // Labels go through 1,2,3 ... CrystalMult
      
      for(int ii=0;ii<devent.Crystal_mult;ii++){
	devent.Ecluster[devent.Cluster_ID[ii]-1]+=devent.Ecrystal[ii];
      };
      // end of cluster analysis
      
      if(0){
	cout << "-------------- After clusterization ---------------- " << endl;
	cout << "	Crystal M = " << devent.Crystal_mult << endl;
	cout << "	Cluster M = " << devent.Cluster_mult << endl;
	cout << "	Label Max = " << Label << endl << endl;
	for(int jj=0;jj<devent.Crystal_mult;jj++){
	  cout << devent.Ecrystal[jj] << "   \t" << devent.Crystal_ID[jj] << "   \t" << devent.Cluster_ID[jj] <<"    \t"<<devent.tof[jj]<<  endl;
	}
      }
    }  //Done checking mult > 1
    

    //Fill DANCE Histograms
    //  cout<<"En: "<<devent.En<<"  ESum: "<<devent.ESum<<"  Mcl:"<<devent.Cluster_mult<<" Mcr: "<<devent.Crystal_mult<<" last T0: "<<last_t0_timestamp<<endl;
    En_Esum_Mcl->Fill(devent.En,devent.ESum,devent.Cluster_mult);
    En_Esum_Mcr->Fill(devent.En,devent.ESum,devent.Crystal_mult);
    
    // TH3F *En_Eg_Mcl; // this should have a Qgate on it
    // TH3F *En_Eg_Mcr; // this should have a Qgate on it
    
    //  TH3F *Esum_Eg_Mcl->Fill(devent.ESum,devedevent.Cluster_mult);; // Eg is Ecluster here
    //   TH3F *Esum_Eg_Mcr; // Eg is Ecrystal here
    
    
  } //Done checking mult > 0
  
  





  //Done with Event Processing

  return 0;
}
