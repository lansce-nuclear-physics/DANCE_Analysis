#ifndef LENZ_UNPACKER_H
#define LENZ_UNPACKER_H

#include <stdint.h>
#include <stdlib.h>     /* atoi */
#include <iomanip>
#include <zlib.h>
#include <vector>
#include <math.h>
#include <deque>

#include <sys/time.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#define MAXNBOARDS 8
#define CHPERBOARD 16

#define NBITSCLOCK 48		// Value 10/1/2015 (New PSD firmware)

//#define VERBOSE

int RunNum;  //run number

//Root Objects

//Histograms

//Histograms
TH2D *hEnergy2D[MAXNBOARDS];
TH1D *hEnergy[MAXNBOARDS][CHPERBOARD];
TH1D *hSWLongGate[MAXNBOARDS][CHPERBOARD];
TH2D *hWaveform[MAXNBOARDS][CHPERBOARD];
TH2D *hScaledWaveform[MAXNBOARDS][CHPERBOARD];
TH2D *hScaledWaveform2[MAXNBOARDS][CHPERBOARD];
TH2D *hCFD[MAXNBOARDS][CHPERBOARD];
TH2D *hCFDTimestamp[MAXNBOARDS];
TH2D *hTDiff[MAXNBOARDS];
TH1D *hCoinc[MAXNBOARDS];
TH2D *hTimeStampDiff[MAXNBOARDS];
TH2D *hTimeStampDiff2[MAXNBOARDS];
TH2D *hBoardMult_vs_Channel[MAXNBOARDS];

TH2D *hHit[MAXNBOARDS];

TH1D *h1m0;
TH1D *h2m1;

TH1D *hRing_m_Wedge;

//Files
TFile *fout;


//vectors
std::vector<double> CFDOut;
std::vector<UShort_t> WaveForm;
std::vector<double> Processed_WaveForm;

//CFD Parameters;
double CFD_scale;
int CFD_delay;
int CFD_threshold;
int CFD_length;
int CFD_gap;

#endif



