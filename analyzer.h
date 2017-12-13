#ifndef ANALYZER_H
#define ANALYZER_H

#include "TFile.h"

#include "structures.h"

#include <vector>


//Functions
int Read_TMatrix();
int Read_DMatrix();
int Initialize_Analyzer(bool read_binary, bool write_binary);
int Analyze_Data(std::vector<DEVT_BANK_wWF> eventvector, bool read_binary, bool write_binary, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold);
int Write_Analyzer_Histograms(TFile *fout, bool read_binary);
int Create_Analyzer_Histograms(bool read_binary);
int Read_PI_Gates();
int Read_Energy_Calibrations(int RunNumber, bool read_binary);
int Make_Output_Binfile(int RunNumber, bool read_binary);
int Read_Moderation_Time_Graphs();




//Histograms



//Diagnostics







#endif
