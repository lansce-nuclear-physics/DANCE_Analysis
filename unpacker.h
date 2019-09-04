//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  unpacker.h             *// 
//*  Last Edit: 09/04/19    *//  
//***************************//

#ifndef UNPACKER_H
#define UNPACKER_H

//C/C++ includes
#include <zlib.h>
#include <bzlib.h>
#include <string.h>
#include <iostream>
#include <stdint.h>
#include <queue>

//ROOT Includes
#include "TFile.h"


//File Includes
#include "structures.h"

using namespace std;

//Function prototypes
int Unpack_Data(queue<gzFile> &gz_queue, double begin, Input_Parameters input_params, Analysis_Parameters *analysis_params);
int Make_DANCE_Map();
int Read_TimeDeviations(Input_Parameters input_params);
int Make_Output_Diagnostics_File(int RunNumber);
int Read_DetectorLoad_Histogram(Input_Parameters input_params);

int Create_Unpacker_Histograms(Input_Parameters input_params);
int Write_Unpacker_Histograms(TFile *fout, Input_Parameters input_params);
int Reset_Unpacker_Histograms(TFile *fout, Input_Parameters input_params);
int Write_Root_File(Input_Parameters input_params, Analysis_Parameters *analysis_params);
double Calculate_Fractional_Time(uint16_t waveform[], uint32_t Ns, uint8_t dual_trace, uint16_t model, Analysis_Parameters *analysis_params);
int Make_Output_Binfile(Input_Parameters input_params);
int Initialize_Unpacker(Input_Parameters input_params);

#endif


