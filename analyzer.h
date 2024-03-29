
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
//*  analyzer.h             *// 
//*  Last Edit: 01/22/20    *//  
//***************************//

#ifndef ANALYZER_H
#define ANALYZER_H

//File Includes
#include "structures.h"

//C/C++ Includes 
#include <vector>

//ROOT Includes
#include "TFile.h"
#include "TCutG.h"
#include "TH1.h"

//Function Prototypes
int Read_TMatrix();
int Read_DMatrix();
int Initialize_Analyzer(Input_Parameters input_params);
int Create_Analyzer_Histograms(Input_Parameters input_params);

int Analyze_Data(std::vector<DEVT_BANK> eventvector,Input_Parameters input_params,Analysis_Parameters *analysis_params);
int Write_Analyzer_Histograms(TFile *fout, Input_Parameters input_params);
int Reset_Analyzer_Histograms(TFile *fout, Input_Parameters input_params);
int Make_Time_Deviations(int RunNumber);

#endif
