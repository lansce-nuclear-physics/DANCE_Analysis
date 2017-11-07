#ifndef UNPACKER_H
#define UNPACKER_H

#include <zlib.h>

int Unpack_Data(gzFile gz_in, double begin, int runnum);
int Make_DANCE_Map();

#endif


