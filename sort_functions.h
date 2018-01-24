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

//Function Protoypes
void heapify(DEVT_BANK arr[], int n, int i);
void heapSort(DEVT_BANK arr[], int n);
void printArray(DEVT_BANK arr[], int n);

#endif
