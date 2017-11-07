#ifndef UNPACKER_H
#define UNPACKER_H

#include <zlib.h>

int Unpack_Data(gzFile gz_in, double begin, int runnum);
int Make_DANCE_Map();
int Create_Unpacker_Histograms();
int Write_Unpacker_Histograms();
int Read_Energy_Calibrations(int RunNumber);

#endif


