//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  unpacker.h             *// 
//*  Last Edit: 05/08/18    *//  
//***************************//

#ifndef UNPACKER_H
#define UNPACKER_H

//C/C++ includes
#include <zlib.h>
#include <string.h>
#include <iostream>

//ROOT Includes
#include "TFile.h"

using namespace std;

//Function prototypes
int Unpack_Data(gzFile &gz_in, double begin, int runnum, bool read_binary, bool write_binary, bool read_simulation, double CoincidenceWindow, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold, bool FitTimeDev,string DataFormat,int NQGates, double QGates[]);
int Make_DANCE_Map();
int Read_TimeDeviations(int runnum, bool FitTimeDev, bool read_simulation);
int Make_Output_Diagnostics_File(int RunNumber);

int Create_Unpacker_Histograms(bool read_binary);
int Write_Unpacker_Histograms(TFile *fout, bool read_binary);

#endif


