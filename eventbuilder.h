//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  eventbuilder.h         *// 
//*  Last Edit: 11/06/18    *//  
//***************************//

#ifndef EVENTBUILDER_H
#define EVENTBUILDER_H

//File includes
#include "structures.h"

//C/C++ includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <deque>
#include <vector>
#include <stdlib.h>

using namespace std;

int build_events(deque<DEVT_BANK> &datadeque, bool &event_building_active, Input_Parameters input_params, uint32_t &Events_Sent);



#endif
