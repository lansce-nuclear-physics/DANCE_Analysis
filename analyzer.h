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


//Function Prototypes
int Read_TMatrix();
int Read_DMatrix();
int Initialize_Analyzer(bool read_binary, bool write_binary, bool read_simulation, int NQGates, double QGates[]);
int Analyze_Data(std::vector<DEVT_BANK> eventvector, bool read_binary, bool write_binary, bool read_simulation, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold, int NQGates, double QGates[]);
int Write_Analyzer_Histograms(TFile *fout, bool read_binary, bool read_simulation, int NQGates, double QGates[]);
int Create_Analyzer_Histograms(bool read_binary, bool read_simulation, int NQGates, double QGates[]);
int Read_PI_Gates();
int Read_Energy_Calibrations(int RunNumber, bool read_binary, bool read_simulation);
int Make_Output_Binfile(int RunNumber, bool read_binary);
int Read_Moderation_Time_Graphs();
int Make_Time_Deviations(int RunNumber);

#endif
