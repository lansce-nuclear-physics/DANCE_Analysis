#ifndef ANALYZER_H
#define ANALYZER_H

#include "TFile.h"

#include "structures.h"

#include <vector>


//Functions
int Read_TMatrix();
int Initialize_Analyzer();
int Analyze_Data(std::vector<DEVT_BANK_wWF> eventvector);
int WriteHistograms(TFile *fout);
int Read_BaF2_Calibrations();
int CreateHistograms();

//Histograms



//Diagnostics







#endif
