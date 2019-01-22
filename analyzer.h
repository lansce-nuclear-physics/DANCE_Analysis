//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  analyzer.h             *// 
//*  Last Edit: 05/08/18    *//  
//***************************//

#ifndef ANALYZER_H
#define ANALYZER_H

//File Includes
#include "structures.h"

//C/C++ Includes 
#include <vector>

//ROOT Includes
#include "TFile.h"
#include "TCutG.h"
#include "TH1.h"

//Function Prototypes
int Read_TMatrix();
int Read_DMatrix();
int Initialize_Analyzer(Input_Parameters input_params);
int Create_Analyzer_Histograms(Input_Parameters input_params);

int Analyze_Data(std::vector<DEVT_BANK> eventvector,Input_Parameters input_params);
int Write_Analyzer_Histograms(TFile *fout, Input_Parameters input_params);
int Read_PI_Gates();
int Read_Energy_Calibrations( Input_Parameters input_params);
int Make_Output_Binfile( Input_Parameters input_params);
int Read_Moderation_Time_Graphs();
int Make_Time_Deviations(int RunNumber);

#endif
