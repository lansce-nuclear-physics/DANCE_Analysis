#ifndef ANALYZER_H
#define ANALYZER_H

#include "TFile.h"

#include "structures.h"

#include <vector>


//Functions
int Read_TMatrix();
int Initialize_Analyzer();
int Analyze_Data(std::vector<DEVT_BANK_wWF> eventvector);
int Write_Analyzer_Histograms(TFile *fout);
int Read_Energy_Calibrations(int RunNumber);
int Create_Analyzer_Histograms();

//Histograms



//Diagnostics







#endif
