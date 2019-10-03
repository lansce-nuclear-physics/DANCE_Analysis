//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  eventbuilder.cpp       *// 
//*  Last Edit: 09/04/19    *//  
//***************************//

//File includes
#include "global.h"
#include "analyzer.h"
#include "eventbuilder.h"
#include "validator.h"
#include "calibrator.h"
#include<iomanip>

using namespace std;

stringstream emsg;

std::vector<DEVT_BANK> DANCE_eventvector;   //Vector to store dance events for analysis
std::vector<DEVT_BANK> BM_eventvector;      //Vector to store beam monitor events for analysis
std::vector<DEVT_BANK> T0_eventvector;      //Vector to store T0 monitor events for analysis

std::ofstream outputbinfile;                //Ouput binary file
DEVT_STAGE1 devt_out;                       //Ouput struct

//Graphs of TOF Corrections
TGraph *gr_DANCE_TOF_Corr;
TGraph *gr_U235_TOF_Corr;
TGraph *gr_Li6_TOF_Corr;
TGraph *gr_He3_TOF_Corr;
  
//Limits of TOF Corrections
double DANCE_TOF_Corr_Limit[2];   //[0] is lower [1] is upper
double U235_TOF_Corr_Limit[2];   //[0] is lower [1] is upper
double Li6_TOF_Corr_Limit[2];   //[0] is lower [1] is upper
double He3_TOF_Corr_Limit[2];   //[0] is lower [1] is upper

//Histograms 
#ifdef Histogram_DetectorLoad
TH1F *hDetectorLoad; //Average detector load vs TOF
TH1F *hDetectorLoad_perT0; //Average detector load vs TOF
TH1F *hDetectorLoad_En; //Average detector load vs TOF
TH1F *hDetectorLoad_En_perT0; //Average detector load vs TOF

double Detector_Load[100000000];
  
TH1F *hEn_BinWidth; //width in ns for each bin in En

#endif

//Reason for invalid
//1 Alpha
//2 Not a Gamma or Alpha
//4 ULD
//8 Threshold
//16 Retrigger
//32 Crystal Blocking
//There can be sums of these if more than one reason
TH1I *hInvalid_Reason;                    

//IDs
TH1I *hID;
TH1I *hID_Raw;
TH1I *hID_Invalid; //fill this when detectors are not valid after event processing
TH1I *hID_gamma; //fill this when detectors have gammas
TH1I *hID_alpha; //fill this when detectors have alphas
TH1I *hID_gamma_Invalid; //fill this when detectors have gammas that are invalid
TH1I *hID_alpha_Invalid; //fill this when detectors have alphas that are invalid
TH1I *hID_invalid_Invalid;


//Raw Slow Intergral - Fast Integral
TH2F *Energy_raw_ID; 

//PSD Histograms
TH2F *ADC_calib;  //2D PSD Plot (calibrated)
TH2F *ADC_raw;    //2D PSD Plot (uncalibrated)

TH3F *ADC_raw_ID;    //3D Plot of 2D PSD plot (uncalibrated) vs DANCE Crystal ID (z)
TH3F *ADC_calib_ID;  //3D Plot of 2D PSD plot (calibrated) vs DANCE Crystal ID (z)

TH2F *ADC_calib_Invalid;  //2D PSD Plot of invalid events (calibrated)
TH2F *ADC_calib_Pileup;  //2D PSD Plot of invalid events (calibrated)
TH2F *ADC_calib_Pileup_Removed;  //2D PSD Plot of invalid events (calibrated)

//Alpha Spectra
TH2F *ADC_alpha;
TH2F *hAlpha;
TH2F *hAlphaCalib;

//Gamma Spectra
TH2F *ADC_gamma;
TH2F *hGamma;
TH2F *hGammaCalib;
TH2F *hGammaCalib_PU;

//Crystal Diagnostics
TH2F *hTimeBetweenCrystals;  //time between subsequent hits of the same crystal (ns)
TH2F *hTimeBetweenCrystals_EnergyRatio; //ratio of the present and last amplitudes vs time difference between that same cyrstal (ns)
TH2F *hTimeBetweenCrystals_FastEnergyRatio; //ratio of the present and last amplitudes vs time difference between that same cyrstal (ns)
TH2F *hTimeBetweenCrystals_LongShortRatio; //ratio of the present long/short integrals vs time difference between hits in the same cyrstal (ns)
  
TH2F *hFastSlowRatio_ID;   //ID vs Ratio of Fast to Slow component



int Initialize_Eventbuilder(Input_Parameters input_params) {
  
  DANCE_Init("Eventbuilder","Initializing");
 
  int func_ret=0;

  //Make histograms
  func_ret = Create_Eventbuilder_Histograms(input_params);
 
  //Moderator Function
  func_ret = Read_Moderation_Time_Graphs();
  
  //clear the event vectors
  DANCE_eventvector.clear();
  BM_eventvector.clear();
  T0_eventvector.clear();

#ifdef Histogram_DetectorLoad
  for(int eye=0; eye<100000000; eye++) {
    Detector_Load[eye]=0;
  }
#endif

  //If we want to write the output to binary
  if(input_params.Write_Binary==1) {
    
    DANCE_Info("Eventbuilder","Output binary requested.  Creating Output binary file");

    stringstream outfilename;
    outfilename.str();
    
    //stage0 
    if(input_params.Analysis_Stage==0) {
      outfilename << STAGE0_BIN; 
      outfilename <<"/stage0_run_";
    }
    //stage1
    if(input_params.Analysis_Stage==1) {
      outfilename << STAGE1_BIN;
      outfilename <<"/stage1_run_";
    }
    outfilename << input_params.RunNumber;
    if (input_params.SingleSubrun){
      outfilename << "_" << input_params.SubRunNumber; 
    }
    outfilename << ".bin";
    
    outputbinfile.open(outfilename.str().c_str(), ios::out | ios::binary);
    
    if(outputbinfile.is_open()) {

      emsg.str("");
      emsg<<"Succesfully created and opened output binary file: "<<outfilename.str();
      
      DANCE_Success("Eventbuilder",emsg.str());
    }
    else {
      emsg.str("");
      emsg<<"Failed to create output binary file: "<<outfilename.str();
      DANCE_Error("Eventbuilder",emsg.str());
      
      func_ret += -1;
    }
  }
  
  if(func_ret==0) {
    DANCE_Success("Eventbuilder", "Initizialized");
  }
  else {
    DANCE_Error("Eventbuilder","Initialization Failed. Exiting!");
  }
  return func_ret;
}

int Build_Events(deque<DEVT_BANK> &datadeque, Input_Parameters input_params, Analysis_Parameters *analysis_params) {
  
#ifdef Eventbuilder_Verbose
  cout<<"Eventbuilder: About to event build deque size: "<<datadeque.size()<<endl;
#endif

  //Eventbuild
  while(true) {
    		  
    //check to see if the buffer is longer than the length specificed in global.h
    if((datadeque[datadeque.size()-1].timestamp - datadeque[0].timestamp) >= (double)1000000000.0*input_params.Buffer_Depth) {
	 
      analysis_params->entries_processed++;
      
      //we have started to build
      analysis_params->event_building_active=true;

      //First write the data to the output binary file if needed
      if(input_params.Write_Binary==1 && outputbinfile.is_open()) {
	devt_out.Ifast = datadeque[0].Ifast;
	devt_out.Islow = datadeque[0].Islow;
	devt_out.timestamp = datadeque[0].timestamp;
	devt_out.ID = datadeque[0].ID;
	outputbinfile.write(reinterpret_cast<char*>(&devt_out),sizeof(DEVT_STAGE1));
	analysis_params->entries_written_to_binary++;
      }    
      
      //ID
      hID_Raw->Fill(datadeque[0].channel+(datadeque[0].board*16));  //Channel + (Board *16)
      hID->Fill(datadeque[0].ID,1);

      Energy_raw_ID->Fill(datadeque[0].Islow,datadeque[0].ID,1);

      //Calculate TOF now (difference between this timestamp and the last T0
      datadeque[0].TOF = datadeque[0].timestamp - analysis_params->last_timestamp[T0_ID];

      
      //Calculate corrected DANCE TOF
      if(datadeque[0].ID < 162) {
	//Correct the TOF for moderation time between ~0 and ~10 MeV
	if( datadeque[0].TOF >= DANCE_TOF_Corr_Limit[0] && datadeque[0].TOF <= DANCE_TOF_Corr_Limit[1]) {
	  datadeque[0].TOF_Corr = gr_DANCE_TOF_Corr->Eval(datadeque[0].TOF);
	}
	else {
	  datadeque[0].TOF_Corr = -1;
	}
      }

      //Calculate corrected U235 TOF
      if(datadeque[0].ID == U235_ID) {
	//Correct the TOF for moderation time between ~0 and ~10 MeV
	if( datadeque[0].TOF >= U235_TOF_Corr_Limit[0] && datadeque[0].TOF <= U235_TOF_Corr_Limit[1]) {
	  datadeque[0].TOF_Corr = gr_U235_TOF_Corr->Eval(datadeque[0].TOF);
	}
	else {
	  datadeque[0].TOF_Corr = -1;
	}
      }

      //Calculate corrected Li6 TOF
      if(datadeque[0].ID == Li6_ID) {
	//Correct the TOF for moderation time between ~0 and ~10 MeV
	if( datadeque[0].TOF >= Li6_TOF_Corr_Limit[0] && datadeque[0].TOF <= Li6_TOF_Corr_Limit[1]) {
	  datadeque[0].TOF_Corr = gr_Li6_TOF_Corr->Eval(datadeque[0].TOF);
	}
	else {
	  datadeque[0].TOF_Corr = -1;
	}
      }

      //Calculate corrected He3 TOF
      if(datadeque[0].ID == He3_ID) {
	//Correct the TOF for moderation time between ~0 and ~10 MeV
	if( datadeque[0].TOF >= He3_TOF_Corr_Limit[0] && datadeque[0].TOF <= He3_TOF_Corr_Limit[1]) {
	  datadeque[0].TOF_Corr = gr_He3_TOF_Corr->Eval(datadeque[0].TOF);
	}
	else {
	  datadeque[0].TOF_Corr = -1;
	}
      }
      
      //Using the same corrections as DANCE (maybe this is good maybe not. IDK)
      //Calculate corrected Bkg TOF
      if(datadeque[0].ID == Bkg_ID) {
	//Correct the TOF for moderation time between ~0 and ~10 MeV
	if( datadeque[0].TOF >= DANCE_TOF_Corr_Limit[0] && datadeque[0].TOF <= DANCE_TOF_Corr_Limit[1]) {
	  datadeque[0].TOF_Corr = gr_DANCE_TOF_Corr->Eval(datadeque[0].TOF);
	}
	else {
	  datadeque[0].TOF_Corr = -1;
	}
      }
      

      //Do the calibrations and validity checks for DANCE
      if(datadeque[0].ID < 162) {
	
	//Calibrate the energy
	Calibrate_DANCE(&datadeque[0]);

	//Check to see crystal blocking time
	Check_Crystal_Blocking(&datadeque[0],analysis_params,input_params);
	
	//The crystal blocking time has to be checked first since it effects the "effective" detector load for the deadtime code. 
	//This mimics having a longer integration window for the long charge integral

	//Detector load	
#ifdef Histogram_DetectorLoad
	if(input_params.Read_Simulation==0) {
	  if(datadeque[0].Valid) {
	    if(analysis_params->last_timestamp[T0_ID] > 0) {
	      uint32_t temptof = (uint32_t) gr_DANCE_TOF_Corr->Eval((datadeque[0].timestamp-analysis_params->last_timestamp[T0_ID]));
	      if(temptof>0) {
		if(input_params.Crystal_Blocking_Time > input_params.Long_Gate) {
		  for(int el=(int)temptof; el<((int)(temptof+input_params.Crystal_Blocking_Time)); el++) {
		    if(el >=0 && el < 100000000) {
		      Detector_Load[el]+=1.0;
		    }
		  }
		}
		else {
		  for(int el=(int)temptof; el<((int)(temptof+input_params.Long_Gate)); el++) {
		    if(el >=0 && el < 100000000) {
		      Detector_Load[el]+=1.0;
		    }
		  }
		}
	      } //end check on temptof
	    } //end of check on last T0
	  } //end of check on valid
	} //end check on simulations
#endif

	  
	//Check Upper Level Discriminator
	Check_ULD(&datadeque[0]);
	  
	//Check Threshold 
	Check_Threshold(&datadeque[0],input_params);

	//Check to see if it is a Retrigger
	Check_Retrigger(&datadeque[0],analysis_params);

	//If still Valid
	if(datadeque[0].Valid == 1) {
	  
	  if(datadeque[0].Islow>0) {
	    //Fill the fast to slow ratio plots now
	    hFastSlowRatio_ID->Fill((1.0*datadeque[0].Ifast)/(1.0*datadeque[0].Islow),datadeque[0].ID,1);
	  }
	  //Fill 2D ADC Calib
	  ADC_calib->Fill(datadeque[0].Eslow, datadeque[0].Efast,1);
	  ADC_raw->Fill(datadeque[0].Islow, datadeque[0].Ifast,1);
	  if(datadeque[0].pileup_detected==1) {
	    ADC_calib_Pileup->Fill(datadeque[0].Eslow, datadeque[0].Efast,1);
	  }
	  else {
	    ADC_calib_Pileup_Removed->Fill(datadeque[0].Eslow, datadeque[0].Efast,1);
	  }

	  //Fill 3D ADC Calib vs Detector
	  ADC_calib_ID->Fill(datadeque[0].Eslow, datadeque[0].Efast, datadeque[0].ID,1);
	  ADC_raw_ID->Fill(datadeque[0].Islow, datadeque[0].Ifast,datadeque[0].ID,1);
	    
	  //Check to see if it is an Alpha
	  Check_Alpha(&datadeque[0]);
	    
	  //If it is not an alpha check to see if it is a Gamma
	  if(!datadeque[0].IsAlpha) {
	    Check_Gamma(&datadeque[0]);
	  } //End of check on Alpha
	   
	  if(datadeque[0].IsAlpha) {
	    hAlpha->Fill(datadeque[0].Islow, datadeque[0].ID,1);
	    hAlphaCalib->Fill(datadeque[0].Eslow, datadeque[0].ID,1);
	    ADC_alpha->Fill(datadeque[0].Eslow, datadeque[0].Efast,1);	// JU diagnostic histogram
	    hID_alpha->Fill(datadeque[0].ID,1);
	  }
	  if(datadeque[0].IsGamma) {
	    hGamma->Fill(datadeque[0].Islow, datadeque[0].ID,1);
	    if(datadeque[0].pileup_detected==1) {
	      hGammaCalib_PU->Fill(datadeque[0].Eslow, datadeque[0].ID,1);
	    }
	    else {
	      hGammaCalib->Fill(datadeque[0].Eslow, datadeque[0].ID,1);
	    }
	    ADC_gamma->Fill(datadeque[0].Eslow, datadeque[0].Efast,1);	// JU diagnostic histogram
	    hID_gamma->Fill(datadeque[0].ID,1);
	  }
	    
	}
        else {
          ADC_calib_Invalid->Fill(datadeque[0].Eslow, datadeque[0].Efast,1);
        }


	//Things that are not gammas are not in DANCE events
	if(!datadeque[0].IsGamma) {
	  datadeque[0].Valid = 0;
	  //Its either an alpha
	  if(datadeque[0].IsAlpha) {
	    datadeque[0].InvalidReason += 1;
	  }
	  //or not an alpha or a gamma
	  else { 
            datadeque[0].InvalidReason += 2;
	  }
	}

	//Fill some DANCE histograms
	//Time between DANCE crystals
	hTimeBetweenCrystals->Fill((datadeque[0].timestamp-analysis_params->last_timestamp[datadeque[0].ID]),datadeque[0].ID,1);

	//Ratio of Efast of n/(n-1) hits vs time between n and n-1 hit
	if(analysis_params->last_Efast[datadeque[0].ID] > 0) {
	  hTimeBetweenCrystals_FastEnergyRatio->Fill((datadeque[0].timestamp-analysis_params->last_timestamp[datadeque[0].ID]),
						 (datadeque[0].Efast/analysis_params->last_Efast[datadeque[0].ID]),1);
	}
        
        //Ratio of Energy of n/(n-1) hits vs time between n and n-1 hit
	if(analysis_params->last_Eslow[datadeque[0].ID] > 0) {
	  hTimeBetweenCrystals_EnergyRatio->Fill((datadeque[0].timestamp-analysis_params->last_timestamp[datadeque[0].ID]),
						 (datadeque[0].Eslow/analysis_params->last_Eslow[datadeque[0].ID]),1);
	}
	
	//Long/Short ratio vs Time between crystals
	if(datadeque[0].Ifast > 0) {
	  hTimeBetweenCrystals_LongShortRatio->Fill((datadeque[0].timestamp-analysis_params->last_timestamp[datadeque[0].ID]),
						    datadeque[0].Islow/datadeque[0].Ifast,1);
	}
	
	
      } //End of check on DANCE Ball

      //Check on time between T0s to remove any "junk" T0s.  
      //Nominal is 20 Hz (50e6 ns) so if they are not at least 1e6 ns apart then they are useless...
      if(datadeque[0].ID == T0_ID) {
	if(((datadeque[0].timestamp - analysis_params->last_timestamp[T0_ID]) < 1e6) && (analysis_params->last_timestamp[T0_ID])>0 ) {
	  datadeque[0].Valid=0;      
	}
      }
#ifdef InvalidDetails
      if (datadeque[0].Valid !=1) {
	  
         if (analysis_params->last_Alpha[datadeque[0].ID])
            hID_alpha_Invalid->Fill(datadeque[0].ID,1);
         if (analysis_params->last_Gamma[datadeque[0].ID])
            hID_gamma_Invalid->Fill(datadeque[0].ID,1);
         if (analysis_params->last_InvalidReason[datadeque[0].ID]>1){
           hID_invalid_Invalid->Fill(datadeque[0].ID,1);
         }
      }
#endif

      //Update analysis params
      analysis_params->last_timestamp[datadeque[0].ID] = datadeque[0].timestamp;
      analysis_params->last_Islow[datadeque[0].ID] = datadeque[0].Islow;
      analysis_params->last_Eslow[datadeque[0].ID] = datadeque[0].Eslow;
      analysis_params->last_Efast[datadeque[0].ID] = datadeque[0].Efast;
      analysis_params->last_Alpha[datadeque[0].ID] = datadeque[0].IsAlpha;
      analysis_params->last_Gamma[datadeque[0].ID] = datadeque[0].IsGamma;
      analysis_params->last_InvalidReason[datadeque[0].ID] = datadeque[0].InvalidReason;

      //If not valid
      if(datadeque[0].Valid != 1) {
	analysis_params->entries_invalid++;
	hID_Invalid->Fill(datadeque[0].ID);

	hInvalid_Reason->Fill(datadeque[0].InvalidReason);
#ifdef Eventbuilder_Verbose
	cout<<RED<<"Eventbuilder: throwing away ID "<<datadeque[0].ID<<RESET<<endl;
#endif
	datadeque.pop_front();  //remove the first entry in the deque

      }       
      //Eventbuild the valid entries
      else {

	//If valid update the analysis params
	analysis_params->last_valid_timestamp[datadeque[0].ID] = datadeque[0].timestamp;
	analysis_params->last_valid_Islow[datadeque[0].ID] = datadeque[0].Islow;
	analysis_params->last_valid_Eslow[datadeque[0].ID] = datadeque[0].Eslow;

	//Handle the DANCE Ball
	if(datadeque[0].ID < 162) {
	
	  //first thing just goes
	  if(DANCE_eventvector.size() == 0) {
	    DANCE_eventvector.push_back(datadeque[0]); //put the first event in the events vector
	    datadeque.pop_front();  //remove the first entry in the deque
	    analysis_params->entries_built++;
	  }
	  //subsequent things are subject to coincidence windows
	  else{
	    //In the window
	    if(datadeque[0].timestamp-DANCE_eventvector[0].timestamp < input_params.Coincidence_Window) {
	      DANCE_eventvector.push_back(datadeque[0]); //put the first event in the events vector
	      datadeque.pop_front();  //remove the first entry in the deque     
	      analysis_params->entries_built++;
         if (DANCE_eventvector.size()>160) {cout <<setprecision(14)<< datadeque[0].timestamp-DANCE_eventvector[0].timestamp<<" " << DANCE_eventvector.size() << endl;} 
	    }
	    //Out of the window
	    else {
	      //Analyze
#ifdef Eventbuilder_Verbose
	      cout<<"Eventbuilder: Processing DANCE Event with Size: "<<DANCE_eventvector.size()<<"  " <<datadeque.size()<<" Entries in the deque"<<endl;


#endif
	      //Send it to the analyzer
	      Analyze_Data(DANCE_eventvector, input_params, analysis_params);

	      //Clear
	      DANCE_eventvector.clear();
	    
	      //Put the entry at the start of the vector
	      DANCE_eventvector.push_back(datadeque[0]); //put the first event in the events vector
	      datadeque.pop_front();  //remove the first entry in the deque     
	      analysis_params->events_built++;
	      analysis_params->entries_built++;
          
  	    }
	  } 	
	}
	
	else if(datadeque[0].ID == Li6_ID || datadeque[0].ID == He3_ID ||  datadeque[0].ID == U235_ID ||  datadeque[0].ID == Bkg_ID) {
	
	  BM_eventvector.push_back(datadeque[0]); //put the first event in the events vector
	  datadeque.pop_front();  //remove the first entry in the deque
	  analysis_params->events_built++;
	  analysis_params->entries_built++;

#ifdef Eventbuilder_Verbose
	  cout<<"Eventbuilder: Processing BM Event with Size: "<<BM_eventvector.size()<<"  " <<datadeque.size()<<" Entries in the deque"<<endl;
#endif
	  //Send it to the analyzer
	   Analyze_Data(BM_eventvector, input_params, analysis_params);  
	
	  //Clear
	  BM_eventvector.clear();
	}
       	
	else if(datadeque[0].ID == T0_ID) {
	  T0_eventvector.push_back(datadeque[0]); //put the first event in the events vector
	  datadeque.pop_front();  //remove the first entry in the deque
	  analysis_params->events_built++;
	  analysis_params->entries_built++;

#ifdef Eventbuilder_Verbose
	  cout<<"Eventbuilder: Processing T0 Event with Size: "<<T0_eventvector.size()<<"  " <<datadeque.size()<<" Entries in the deque"<<endl;
#endif
	  //Send it to the analyzer
	  Analyze_Data(T0_eventvector, input_params, analysis_params);   
	  
	  //Clear
	  T0_eventvector.clear();
	}

	else {
#ifdef Eventbuilder_Verbose
	  cout<<RED<<"Eventbuilder: throwing away ID "<<datadeque[0].ID<<RESET<<endl;
#endif
	  analysis_params->Unknown_entries++;
	  datadeque.pop_front();  //remove the first entry in the deque     
	}
      } //End of check on valid event

      if(datadeque.size()==0) {
	//	event_build=false;
	break;
      }
      
    }
    else {
#ifdef Eventbuilder_Verbose
      cout<<"Eventbuilder: event build complete: "<<datadeque.size()<<endl;
#endif
      break;
    }
  }


  return 0;

}

//This function fetches the TOF Correction plots for the moderation time 
int Read_Moderation_Time_Graphs() {
 
  DANCE_Info("Eventbuilder","Reading TOF Corrections");

  ifstream tof_corr;
  tof_corr.open("./TOF_Corrections/TOF_Corrections.txt");

  if(tof_corr.is_open()) {
    //Find out how many points there are 
    int npoints;
    tof_corr>>npoints;
    const int N=npoints;
    
    //Time of flight from neutron energy
    double DANCE_TOF[N];
    double U235_TOF[N];
    double He3_TOF[N];
    double Li6_TOF[N];

    //Time of flight plus moderation time
    double DANCE_TOF_Measured[N];
    double U235_TOF_Measured[N];
    double He3_TOF_Measured[N];
    double Li6_TOF_Measured[N];

    for(int eye=0; eye<N; eye++) {
      tof_corr>>DANCE_TOF[eye]>>DANCE_TOF_Measured[eye]>>U235_TOF[eye]>>U235_TOF_Measured[eye]>>Li6_TOF[eye]>>Li6_TOF_Measured[eye]>>He3_TOF[eye]>>He3_TOF_Measured[eye];
    }
    
    //Graphs of TOF Corrections
    gr_DANCE_TOF_Corr = new TGraph(N,DANCE_TOF_Measured,DANCE_TOF);
    gr_U235_TOF_Corr = new TGraph(N,U235_TOF_Measured,U235_TOF);
    gr_He3_TOF_Corr = new TGraph(N,He3_TOF_Measured,He3_TOF);
    gr_Li6_TOF_Corr = new TGraph(N,Li6_TOF_Measured,Li6_TOF);
    
    DANCE_TOF_Corr_Limit[0] = DANCE_TOF_Measured[N-1];
    DANCE_TOF_Corr_Limit[1] = DANCE_TOF_Measured[0];
    U235_TOF_Corr_Limit[0] = U235_TOF_Measured[N-1];
    U235_TOF_Corr_Limit[1] = U235_TOF_Measured[0];
    Li6_TOF_Corr_Limit[0] = Li6_TOF_Measured[N-1];
    Li6_TOF_Corr_Limit[1] = Li6_TOF_Measured[0];
    He3_TOF_Corr_Limit[0] = He3_TOF_Measured[N-1];
    He3_TOF_Corr_Limit[1] = He3_TOF_Measured[0];
       
    DANCE_Success("Eventbuilder","Read TOF Corrections");
    return 0;
  }
  else {
    DANCE_Error("Eventbuilder","Faild to Read TOF Corrections");
    return -1;
  }  
}


int Create_Eventbuilder_Histograms(Input_Parameters input_params) {

  DANCE_Info("Eventbuilder","Creating Histograms");

  hInvalid_Reason = new TH1I("hInvalid_Reason","hInvalid_Reason",64,0,64);

  //ID
  hID = new TH1I("hID","hID",256,0,256);
  hID_Raw = new TH1I("hID_Raw","hID_Raw",256,0,256);
  hID_Invalid = new TH1I("hID_Invalid","hID_Invalid",256,0,256);
  hID_gamma = new TH1I("hID_Gamma","hID_Gamma",256,0,256);
  hID_alpha = new TH1I("hID_Alpha","hID_Alpha",256,0,256);
#ifdef InvalidDetails
  hID_gamma_Invalid = new TH1I("hID_Gamma_Invalid","hID_Gamma_Invalid",256,0,256);
  hID_alpha_Invalid = new TH1I("hID_Alpha_Invalid","hID_Alpha_Invalid",256,0,256);
  hID_invalid_Invalid = new TH1I("hID_Invalid_Invalid","hID_Invalid_Invalid",256,0,256);
#endif

  //Raw Energy
  Energy_raw_ID =  new TH2F("Energy_raw_ID","Energy_raw_ID",3500,0.0,70000,162,0,162);
  
  //PSD Histograms
  ADC_raw = new TH2F("ADC_raw","ADC_raw",1800,0.,72000.,180,0.,7200);
  ADC_calib = new TH2F("ADC_calib","ADC_calib",2400,0.,24.,1000,0.,10.);
  ADC_raw_ID = new TH3F("ADC_raw_ID","ADC_raw_ID",1800,0.0,72000.0,180,0.0,7200,162,0,162);
  ADC_calib_ID = new TH3F("ADC_calib_ID","ADC_calib_ID",600,0,24,250,0,10,162,0,162);
  ADC_calib_Invalid = new TH2F("ADC_calib_Invalid","ADC_calib_Invalid",2400,0,24,1000,0,10);
    ADC_calib_Pileup = new TH2F("ADC_calib_Pileup","ADC_calib_Pileup",2400,0,24,1000,0,10);
    ADC_calib_Pileup_Removed = new TH2F("ADC_calib_Pileup_Removed","ADC_calib_Pileup_Removed",2400,0,24,1000,0,10);

  //Gamma Histograms
  ADC_gamma = new TH2F("ADC_gamma","ADC_gamma",2400,0,24,1000,0,10);	// JU
  hGamma = new TH2F("hGamma","hGamma",3500,0,70000,162,0,162);
  hGammaCalib = new TH2F("hGammaCalib","hGammaCalib",2400,0,24,162,0,162);
  hGammaCalib_PU = new TH2F("hGammaCalib_PU","hGammaCalib_PU",2400,0,24,162,0,162);

  //Alpha Histograms
  ADC_alpha = new TH2F("ADC_alpha","ADC_alpha",2400,0.0,24.0,1000,0.0,10.0);	// JU
  hAlpha = new TH2F("hAlpha","hAlpha",1500,0,30000,162,0,162);
  hAlphaCalib = new TH2F("hAlphaCalib","hAlphaCalib",500,0.0,5.0,162,0,162);

  //Crystal Diagnostics
  hTimeBetweenCrystals = new TH2F("TimeBetweenCrystals","TimeBetweenCrystals",10000,0,10000,162,0,162);
  hTimeBetweenCrystals_EnergyRatio = new TH2F("TimeBetweenCrystals_EnergyRatio","TimeBetweenCrystals_EnergyRatio",1250,0,10000,1000,0,20);
  hTimeBetweenCrystals_FastEnergyRatio = new TH2F("TimeBetweenCrystals_FastEnergyRatio","TimeBetweenCrystals_EnergyRatio",1250,0,10000,1000,0,20);
  hTimeBetweenCrystals_LongShortRatio = new TH2F("TimeBetweenCrystals_LongShortRatio","TimeBetweenCrystals_LongShortRatio",1250,0,10000,1000,0,100);

  hFastSlowRatio_ID = new TH2F("FastSlowRatio_ID","FastSlowRatio_ID",1000,0,1,162,0,162);


  //Detector Load
#ifdef Histogram_DetectorLoad
  double x[5000];
  int NEbins=0;
  
  for(double lx=log10(NeutronE_From);lx<log10(NeutronE_To);lx=lx+(1./NeutronE_BinsPerDecade)){
    x[NEbins]=pow(10,lx);
    NEbins++;
  }
  NEbins--;
  
  if(input_params.Read_Simulation==0) {
    hDetectorLoad = new TH1F("DetectorLoad","DetectorLoad",10000000,0,10000000);
    hDetectorLoad_perT0 = new TH1F("DetectorLoad_perT0","DetectorLoad_perT0",10000000,0,10000000);
    hDetectorLoad_En = new TH1F("DetectorLoad_En","DetectorLoad_En",NEbins,x);
    hDetectorLoad_En_perT0 = new TH1F("DetectorLoad_En_perT0","DetectorLoad_En_perT0",NEbins,x);
  }

  hEn_BinWidth = new TH1F("En_BinWidth","En_BinWidth",NEbins,x);


#endif
  DANCE_Success("Eventbuilder","Created Histograms");
  return 0;
}

int Write_Eventbuilder_Histograms(TFile *fout,Input_Parameters input_params, Analysis_Parameters *analysis_params) {

  DANCE_Info("Eventbuilder","Writing Histograms");
  
  fout->cd();

  hInvalid_Reason->Write();
  
  hID->Write();
  hID_Raw->Write();
  hID_Invalid->Write();
  hID_gamma->Write();
  hID_alpha->Write();
#ifdef InvalidDetails
  hID_gamma_Invalid->Write();
  hID_alpha_Invalid->Write();
  hID_invalid_Invalid->Write();
#endif

  Energy_raw_ID->Write();

  ADC_calib->Write();
  ADC_raw->Write();
  ADC_calib_Invalid->Write();
      ADC_calib_Pileup->Write();
      ADC_calib_Pileup_Removed->Write();

  ADC_raw_ID->Write();
  ADC_calib_ID->Write();
    
  ADC_alpha -> Write();
  hAlpha->Write();
  hAlphaCalib->Write();
    
  ADC_gamma -> Write();
  hGamma->Write();
  hGammaCalib->Write();
  hGammaCalib_PU->Write();

  hTimeBetweenCrystals->Write();
  hTimeBetweenCrystals_EnergyRatio->Write();
  hTimeBetweenCrystals_FastEnergyRatio->Write();
  hTimeBetweenCrystals_LongShortRatio->Write();

  hFastSlowRatio_ID->Write();

#ifdef Histogram_DetectorLoad

  if(input_params.Read_Simulation==0) { 
    
    //Detector Load
    for(int eye=0; eye<10000000; eye++) {
      hDetectorLoad_perT0->SetBinContent(eye+1,Detector_Load[eye]/(160.0*analysis_params->T0_events_analyzed));
      hDetectorLoad->SetBinContent(eye+1,Detector_Load[eye]/(160.0));
    }

    //make the deadtime relative to En
    for(int eye=0; eye< hDetectorLoad_En->GetNbinsX(); eye++) {
      
      double temptof = 1.0e9*DANCE_FlightPath / TMath::Sqrt(2*(hDetectorLoad_En->GetXaxis()->GetBinCenter(eye+1))*speedoflight*speedoflight/neutronmass);
      //                        meters                          eV                                         m/s m/s eV/c^2
      
      hDetectorLoad_En->SetBinContent( eye+1, hDetectorLoad->GetBinContent((int)temptof));
      hDetectorLoad_En_perT0->SetBinContent( eye+1, hDetectorLoad_perT0->GetBinContent((int)temptof));
    }
    
    hDetectorLoad_perT0->Write();
    hDetectorLoad->Write();
    hDetectorLoad_En_perT0->Write();
    hDetectorLoad_En->Write();

    //Bin Width
    for(int eye=0; eye<hEn_BinWidth->GetNbinsX(); eye++) {
      double temptoflow = 1.0e9*DANCE_FlightPath / TMath::Sqrt(2*(hEn_BinWidth->GetXaxis()->GetBinLowEdge(eye+1))*speedoflight*speedoflight/neutronmass);
      double temptofhigh = 1.0e9*DANCE_FlightPath / TMath::Sqrt(2*(hEn_BinWidth->GetXaxis()->GetBinLowEdge(eye+1)+hEn_BinWidth->GetXaxis()->GetBinWidth(eye+1))*speedoflight*speedoflight/neutronmass);
      hEn_BinWidth->SetBinContent(eye+1,temptoflow-temptofhigh);
    }
    hEn_BinWidth->Write();
  }
#endif
  
  DANCE_Success("Eventbuilder","Wrote Histograms");
  return 0;
}

int Reset_Eventbuilder_Histograms(TFile *fout,Input_Parameters input_params, Analysis_Parameters *analysis_params) {

  fout->cd();

  hInvalid_Reason->Reset("ICES");
  
  hID->Reset("ICES");
  hID_Raw->Reset("ICES");
  hID_Invalid->Reset("ICES");
  hID_gamma->Reset("ICES");
  hID_alpha->Reset("ICES");
#ifdef InvalidDetails
  hID_gamma_Invalid->Reset("ICES");
  hID_alpha_Invalid->Reset("ICES");
  hID_invalid_Invalid->Reset("ICES");
#endif
  Energy_raw_ID->Reset("ICES");

  ADC_calib->Reset("ICES");
  ADC_raw->Reset("ICES");
  ADC_calib_Invalid->Reset("ICES");
      ADC_calib_Pileup->Reset("ICES");
      ADC_calib_Pileup_Removed->Reset("ICES");

  ADC_raw_ID->Reset("ICES");
  ADC_calib_ID->Reset("ICES");
    
  ADC_alpha -> Reset("ICES");
  hAlpha->Reset("ICES");
  hAlphaCalib->Reset("ICES");
    
  ADC_gamma -> Reset("ICES");
  hGamma->Reset("ICES");
  hGammaCalib->Reset("ICES");
  hGammaCalib_PU->Reset("ICES");

  hTimeBetweenCrystals->Reset("ICES");
  hTimeBetweenCrystals_EnergyRatio->Reset("ICES");
  hTimeBetweenCrystals_FastEnergyRatio->Reset("ICES");
  hTimeBetweenCrystals_LongShortRatio->Reset("ICES");

  hFastSlowRatio_ID->Reset("ICES");

#ifdef Histogram_DetectorLoad

  if(input_params.Read_Simulation==0) { 
    
    hDetectorLoad_perT0->Reset("ICES");
    hDetectorLoad->Reset("ICES");
    hDetectorLoad_En_perT0->Reset("ICES");
    hDetectorLoad_En->Reset("ICES");
    hEn_BinWidth->Reset("ICES");
  }
#endif
  return 0;
}
