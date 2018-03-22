//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  analyzer.h             *// 
//*  Last Edit: 03/22/18    *//  
//***************************//

#ifndef ANALYZER_H
#define ANALYZER_H

//File Includes
#include "structures.h"

//C/C++ Includes 
#include <vector>

//ROOT Includes
#include "TFile.h"

//Function Prototypes
int Read_TMatrix();
int Read_DMatrix();
int Initialize_Analyzer(bool read_binary, bool write_binary,int NQGates, double QGates[]);
int Analyze_Data(std::vector<DEVT_BANK> eventvector, bool read_binary, bool write_binary, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold, int NQGates, double QGates[]);
int Write_Analyzer_Histograms(TFile *fout, bool read_binary,int NQGates, double QGates[]);
int Create_Analyzer_Histograms(bool read_binary, int NQGates, double QGates[]);
int Read_PI_Gates();
int Read_Energy_Calibrations(int RunNumber, bool read_binary);
int Make_Output_Binfile(int RunNumber, bool read_binary);
int Read_Moderation_Time_Graphs();
int Make_Time_Deviations(int RunNumber);

#endif
