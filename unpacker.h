//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  unpacker.h             *// 
//*  Last Edit: 01/23/18    *//  
//***************************//

#ifndef UNPACKER_H
#define UNPACKER_H

//C/C++ includes
#include <zlib.h>
#include <string.h>
#include <iostream>

using namespace std;

//Function prototypes
int Unpack_Data(gzFile &gz_in, double begin, int runnum, bool read_binary, bool write_binary, double CoincidenceWindow, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold, bool FitTimeDev,string DataFormat);
int Make_DANCE_Map();
int Create_Unpacker_Histograms();
int Write_Unpacker_Histograms();
int Read_TimeDeviations(int runnum, bool FitTimeDev);


#endif


