#ifndef ANALYZER_H
#define ANALYZER_H

#include "TFile.h"

#include "structures.h"

#include <vector>


//Functions
int Read_TMatrix();
int Read_DMatrix();
int Initialize_Analyzer();
int Analyze_Data(std::vector<DEVT_BANK_wWF> eventvector);
int Write_Analyzer_Histograms(TFile *fout);
int Create_Analyzer_Histograms();
int Read_PI_Gates();

//Histograms



//Diagnostics







#endif
