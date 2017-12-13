#ifndef UNPACKER_H
#define UNPACKER_H

#include <zlib.h>

int Unpack_Data(gzFile gz_in, double begin, int runnum, bool read_binary, bool write_binary, double CoincidenceWindow, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold);
int Make_DANCE_Map();
int Create_Unpacker_Histograms();
int Write_Unpacker_Histograms();


#endif


