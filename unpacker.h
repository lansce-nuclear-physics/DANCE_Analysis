
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
//*  Cathleen E. Fry        *//
//*  cfry@lanl.gov          *//
//*  unpacker.h             *// 
//*  Last Edit: 01/22/20    *//  
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
int Write_Root_File(Input_Parameters input_params, Analysis_Parameters *analysis_params);
double Calculate_Fractional_Time(uint16_t waveform[], uint32_t Ns, uint8_t dual_trace, uint16_t model, Analysis_Parameters *analysis_params);
int Make_Output_Binfile(Input_Parameters input_params);
int Initialize_Unpacker(Input_Parameters input_params);

#endif


