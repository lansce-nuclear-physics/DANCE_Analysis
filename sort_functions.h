//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  sort_functions.h       *// 
//*  Last Edit: 01/23/18    *//  
//***************************//

#ifndef SORT_FUNCTIONS_H
#define SORT_FUNCTIONS_H

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

//Function Protoypes
void heapify(DEVT_BANK arr[], int n, int i);
void heapSort(DEVT_BANK arr[], int n);
void printArray(DEVT_BANK arr[], int n);

int sort_array(DEVT_BANK db_arr[], deque<DEVT_BANK> &datadeque, double smallest_timestamp, uint32_t EVTS, Input_Parameters input_params, bool &first_sort, bool event_building_active);

#endif
