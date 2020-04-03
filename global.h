
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
//*  global.h               *// 
//*  Last Edit: 01/24/20    *//  
//***************************//

#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <string>

//Switches
#define Histogram_DetectorLoad     //Thes turns on the detector load histograms for the dead time and pileup analysis
//#define Make_Removed_Spectra       //This turns on the gamma removed spectra for the dead time and pileup analysis
//#define CheckTheDeque              //This turns on time integrity checks of the Buffer
#define Histogram_Waveforms        //This turns on histogramming of waveforms
//#define Histogram_Digital_Probes   //This turns on histogramming of digital probes (Not applicable to CAEN2015 format)
//#define CheckBufferDepth           //This turns on a check of how much of the buffer is being used
//#define MakeTimeStampHistogram     // make histogram of timestamps (very big 1D, only use for debugging)
//#define InvalidDetails             // histograms that are separated by event type before the invalid event
//#define HighRateDebug		    // histograms useful for debugging in weird conditions, normally not so useful

//Verbosity
//#define Calibrator_Verbose       //This turns on the messages from the calibrator
//#define Validator_Verbose        //This turns on the messages from the validator
//#define EventSort_Verbose        //This turns on the messages from the eventsort
//#define Eventbuilder_Verbose     //This turns on the messages from the eventbuilder
//#define Unpacker_Verbose         //This turns on the messages from the unpacker 
//#define Scaler_Verbose           //This turns on the scaler specific messages from the unpacker (stage 0 only)
//#define Diagnostic_Verbose       //This turns on the diagnostics specific messages from the unpacker (stage0 caen2018 only) 
//#define Removed_Verbose          //This turns on the messages from gamma-ray removal (Make_Removed_Spectra must be enabled)


//Status indicators
#define EventLimit 4294967295  //Event limit to shut off the unpacker (2^32 -1)
//#define EventLimit 50000000  //Event limit to shut off the unpacker (2^32 -1)
//#define EventLimit 1000000  //Event limit to shut off the unpacker (2^32 -1)
#define ProgressInterval 1000000  //Progress bar incriments

//Some Global Unpacker Variables (DONT CHANGE UNLESS NEEDED)
//size of the block buffer (this has implications for unpacking speed.  Too many sorts and there will be too much overhead.  Not enough and NlogN is too big.  N*log(N) vs k*n*log(n) where k*n=N)
#define BlockBufferSize 250000 
//size of the DEVT array in the unpacker
#define MaxDEVTArrSize 300000  //this should be a number bigger than the block buffer size but too much bigger or else the RAM load will be high. Enable CheckBufferDepth to see how it is behaving
//Number of Channel aggregates that can be in one channel aggregate unpack
#define Max_ChAgg_Size 65535  //This is the maximum number of words the channel aggregate can be for the read to work properly
#define Max_Gamma_Removed 12 //This is how far to make the gamma removed spectra

//File Inputs
#define TMatrixFile "./Config/TMatrix.txt" 
#define	DanceMapFile "./Config/DanceMap.txt" 
#define DMatrixFile "./Config/DetectorMatrix.txt"

//File Input/Output
#define STAGE0_ROOT "./stage0_root"
#define STAGE1_ROOT "./stage1_root"
#define STAGE0_BIN "./stage0_bin"
#define STAGE1_BIN "./stage1_bin"
#define STAGE0_SIM "./stage0_simulated"
#define DIAGNOSTICS "./diagnostics"
#define TIMEDEV_DIR "../TimeDeviations"
#define CALIB_DIR "../Calibrations"

//IDs
#define T0_ID 200
#define He3_ID 241
#define Li6_ID 244
#define U235_ID 243
#define Bkg_ID 242

//PI Gates
#define GAMMAGATE "Gamma.dat"
#define ALPHAGATE "Alpha.dat"
#define RETRIGGERGATE "Retrigger.dat"

//retrigger waveform ratio gates
#define wf_ratio_low 0.055  //2019
#define wf_ratio_high 0.095 //2019
//#define wf_ratio_low 0.06  //2018
//#define wf_ratio_high 0.12 //2018
//#define wf_ratio_low 0.01  //Tl 2019
//#define wf_ratio_high 0.03 //Tl 2019

//Physics Stuff
#define neutronmass 939.565379e6  //Mass of the neutron in eV/c^2
#define speedoflight 2.997924589e8  //Speed of light in m/s

#define	DANCE_FlightPath 20.234 // was 20.28437 meters and 20.2407 and 20.2572 and 20.2352
#define U235_FlightPath 22.895  // was 22.8 meters and 22.8299 and 22.882
#define He3_FlightPath 22.74 //meters
#define Li6_FlightPath 22.6557 // was 22.607 meters, 22.6121, and 22.7015
#define	Bkg_FlightPath 20.234 // was 20.28437 meters and 20.2407 and 20.2572 and 20.2352

#define DANCE_Delay 230.8239  // ns (Was 523.0 in FARE, and 205.0, 282.26, and 293.82, 204.892 in last version)
#define U235_Delay 412.721 //ns (Was 322.4 in FARE and 553.2 and 578.2201,397.571 in previous version)
#define He3_Delay 359.3 //ns (Was -207 in previous version)
#define Li6_Delay 195.9155 //ns (Was 376 in FARE, and 502.9, 203.336 in previous version)
#define Bkg_Delay 0 //ns 

//Histogramming
#define	NeutronE_From 0.005 //Neutron energy from [eV]
#define	NeutronE_To 5e6 //Neutron energy to [eV]
#define	NeutronE_BinsPerDecade 200.0 // Number of neutron energy bins per decade

#define	GammaE_From 0.0 //Gamma energy [MeV] - low limit
#define	GammaE_To 20.0 //Gamma energy [MeV] - upper limit
#define	GammaE_NoOfBins 220.0 //Number of bins





#endif
