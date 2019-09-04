//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  eventbuilder.h         *// 
//*  Last Edit: 09/04/19    *//  
//***************************//

#ifndef EVENTBUILDER_H
#define EVENTBUILDER_H

//File includes
#include "structures.h"

//C/C++ includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <deque>
#include <vector>
#include <stdlib.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TGraph.h"

int Open_Binary(Input_Parameters input_params);
bool Close_Binary();

int Initialize_Eventbuilder(Input_Parameters input_params);
  
int Build_Events(std::deque<DEVT_BANK> &datadeque, Input_Parameters input_params,Analysis_Parameters *analysis_params);

int Create_Eventbuilder_Histograms(Input_Parameters input_params);
int Write_Eventbuilder_Histograms(TFile *fout,Input_Parameters input_params, Analysis_Parameters *analysis_params);
int Reset_Eventbuilder_Histograms(TFile *fout,Input_Parameters input_params, Analysis_Parameters *analysis_params);
int Read_Moderation_Time_Graphs();
  
 
#endif
