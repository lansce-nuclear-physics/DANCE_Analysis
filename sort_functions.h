
////////////////////////////////////////////////////////////////////////
//                                                                    //
//   Software Name: DANCE Data Acquisition and Analysis Package       //
//     Subpackage: DANCE_Analysis                                     //
//   Identifying Number: C18105                                       // 
//                                                                    //
////////////////////////////////////////////////////////////////////////
//                                                                    //
//                                                                    //
// Copyright 2019.                                                    //
// Triad National Security, LLC. All rights reserved.                 //
//                                                                    //
//                                                                    //
//                                                                    //
// This program was produced under U.S. Government contract           //
// 89233218CNA000001 for Los Alamos National Laboratory               //
// (LANL), which is operated by Triad National Security, LLC          //
// for the U.S. Department of Energy/National Nuclear Security        //
// Administration. All rights in the program are reserved by          //
// Triad National Security, LLC, and the U.S. Department of           //
// Energy/National Nuclear Security Administration. The Government    //
// is granted for itself and others acting on its behalf a            //
// nonexclusive, paid-up, irrevocable worldwide license in this       //
// material to reproduce, prepare derivative works, distribute        //
// copies to the public, perform publicly and display publicly,       //
// and to permit others to do so.                                     //
//                                                                    //
// This is open source software; you can redistribute it and/or       //
// modify it under the terms of the GPLv2 License. If software        //
// is modified to produce derivative works, such modified             //
// software should be clearly marked, so as not to confuse it         //
// with the version available from LANL. Full text of the GPLv2       //
// License can be found in the License file of the repository         //
// (GPLv2.0_License.txt).                                             //
//                                                                    //
////////////////////////////////////////////////////////////////////////



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

int sort_array(DEVT_BANK db_arr[], deque<DEVT_BANK> &datadeque, uint32_t EVTS, Input_Parameters input_params, Analysis_Parameters *analysis_params);

#endif
