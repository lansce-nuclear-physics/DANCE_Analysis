
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
//*  unpack_vx735_vx730.h   *// 
//*  Last Edit: 07/23/18    *//  
//***************************//

#ifndef UNPACK_VX725_VX730_H
#define UNPACK_VX725_VX730_H

#include <stdint.h>

//This is the maximum size of the probe array
#define MAX_PROBE_LENGTH 65535

//Vx725_Vx730 board header struct
struct V1730_Header_t {
  uint32_t dataword_1;
  uint32_t dataword_2;
  uint32_t dataword_3;
  uint32_t dataword_4;
};

//Vx725_Vx730 board header masks
//WORD 1
#define Vx725_Vx730_HEADER_MASK             0xF0000000  //bits 28 to 31 inclusive
#define Vx725_Vx730_BOARDAGGSIZE_MASK       0x0FFFFFFF  //bits 0 to 27 inclusive
//WORD 2
#define Vx725_Vx730_BOARDID_MASK            0xF8000000  //bits 27 to 31 inclusive
#define Vx725_Vx730_PATTERN_MASK            0x004FFF00  //bits 8 to 22 inclusive 
#define Vx725_Vx730_CHANNEL_MASK            0x000000FF  //bits 0 to 7 inclussive
//WORD 3
#define Vx725_Vx730_BOARDAGGCOUNTER_MASK    0x004FFFFF  //bits 0 to 22 inclusive
//WORD 4
#define Vx725_Vx730_BOARDAGGTIME_MASK       0xFFFFFFFF  //bits 0 to 31 inclusive


//Vx725_Vx730 board header data
struct Vx725_Vx730_Board_Data_t {
  //WORD 1
  uint8_t header;                                       //bits 28 to 31 inclusive 
  uint32_t boardaggsize;                                //bits 0 to 27 inclusive
  //WORD 2
  uint8_t boardid;                                      //bits 27 to 31 inclusive
  uint16_t pattern;                                     //bits 8 to 22 inclusive
  uint8_t channelmask;                                  //bits 0 to 7 inclusive
  //WORD 3
  uint32_t boardaggcounter;                             //bits 0 to 22 inclusive
  //WORD 4 
  uint32_t boardaggtime;                                //bits 0 to 31 inclusive
};

int unpack_vx725_vx730_board_data(V1730_Header_t *v1730_header, Vx725_Vx730_Board_Data_t *vx725_vx730_board_data);


//Channel Aggregate Header
struct V1730_ChAgg_Header_t {
  uint32_t dataword_1;
  uint32_t dataword_2;
};

//Vx725_Vx730 PSD masks
//WORD 1
#define Vx725_Vx730_PSD_CHAGGHEADER_MASK    0x80000000  //bit 31
#define Vx725_Vx730_PSD_CHAGGSIZE_MASK      0x003FFFFF  //bits 0 to 21 inclusive 
//WORD 2
#define Vx725_Vx730_PSD_DT_MASK             0x80000000  //bit 31
#define Vx725_Vx730_PSD_EQ_MASK             0x40000000  //bit 30
#define Vx725_Vx730_PSD_ET_MASK             0x20000000  //bit 29
#define Vx725_Vx730_PSD_EE_MASK             0x10000000  //bit 28
#define Vx725_Vx730_PSD_ES_MASK             0x08000000  //bit 27
#define Vx725_Vx730_PSD_EX_MASK             0x07000000  //bits 24 to 26 inclusive
#define Vx725_Vx730_PSD_AP_MASK             0x00C00000  //bits 22 to 23 inclusive
#define Vx725_Vx730_PSD_DP2_MASK            0x00380000  //bits 19 to 21 inclusive
#define Vx725_Vx730_PSD_DP1_MASK            0x00070000  //bits 16 to 18 inclusive
#define Vx725_Vx730_PSD_NSDB8_MASK          0x0000FFFF  //bits 0 to 15 inclusive
//WORD 3
#define Vx725_Vx730_PSD_CH_MASK             0x80000000  //bit 31
#define Vx725_Vx730_PSD_TTT_MASK            0x7FFFFFFF  //bits 0 to 30 inclusive
//WORDS 4 to N-2
#define Vx725_Vx730_PSD_DP2_1_MASK          0x80080000  //bit 31
#define Vx725_Vx730_PSD_DP1_1_MASK          0x40000000  //bit 30
#define Vx725_Vx730_PSD_AP_1_MASK           0x3FFF0000  //bits 16 to 29 inclusive
#define Vx725_Vx730_PSD_DP2_0_MASK          0x00008000  //bit 15
#define Vx725_Vx730_PSD_DP1_0_MASK          0x00004000  //bit 14
#define Vx725_Vx730_PSD_AP_0_MASK           0x00003FFF  //bits 0 to 13 inclusive
//WORD N-1
#define Vx725_Vx730_PSD_EXTRAS_MASK         0xFFFFFFFF  //bits 0 to 31 inclusive
//WORD N
#define Vx725_Vx730_PSD_QLONG_MASK          0xFFFF0000  //bits 16 to 31 inclusive
#define Vx725_Vx730_PSD_PUR_MASK            0x00008000  //bit 15
#define Vx725_Vx730_PSD_QSHORT_MASK         0x00007FFF  //bits 0 to 14 inclusive

//DPP-PSD Vx725_Vx730
struct Vx725_Vx730_PSD_Data_t {
  //WORD 1 
  uint8_t chagg_header;                                 //bit 31
  uint32_t chagg_size;                                  //bits 0 to 21 inclusive 
  //WORD 2
  uint8_t dual_trace;                                   //bit 31
  uint8_t charge_enabled;                               //bit 30
  uint8_t time_enabled;                                 //bit 29
  uint8_t extras_enabled;                               //bit 28
  uint8_t waveform_enabled;                             //bit 27
  uint8_t extras_option;                                //bits 24 to 26 inclusive
  uint8_t ap;                                           //bits 22 to 23 inclusive
  uint8_t dp2;                                          //bits 19 to 21 inclusive
  uint8_t dp1;                                          //bits 16 to 18 inclusive
  uint16_t nsdb8;                                       //bits 0 to 15 inclusive
  //WORD 3
  uint8_t channel;                                      //bit 31
  uint32_t trigger_time_tag;                            //bits 0 to 30 inclusive
  //WORDS 4 to N-2
  uint16_t analog_probe1[MAX_PROBE_LENGTH];             //bits 0 to 13 inclusive
  uint16_t analog_probe2[MAX_PROBE_LENGTH];             //bits 16 to 29 inclusive
  uint8_t digital_probe1[MAX_PROBE_LENGTH];             //bit 31,15
  uint8_t digital_probe2[MAX_PROBE_LENGTH];             //bit 30,14
  //WORD N-1
  uint32_t extras;                                      //bits 0 to 31 inclusive
  uint16_t extended_time_stamp;
  uint16_t baseline_timesfour;
  uint16_t flags;
  uint16_t fine_time_stamp;
  uint16_t lost_trigger_counter;
  uint16_t total_trigger_counter;
  uint16_t cfd_sazc;
  uint16_t cfd_sbzc;
  //WORD N
  uint16_t qlong;                                       //bits 16 to 31 inclusive
  uint8_t pur;                                          //bit 15
  uint16_t qshort;                                      //bits 0 to 14 inclusive

  uint16_t individual_chagg_size;                        //size of the individual channels in the aggregate (in 32-bit words)
};

//Vx725_Vx730 PHA masks
//WORD 1
#define Vx725_Vx730_PHA_CHAGGHEADER_MASK    0x80000000  //bit 31
#define Vx725_Vx730_PHA_CHAGGSIZE_MASK      0x7FFFFFFF  //bits 0 to 30 inclusive 
//WORD 2
#define Vx725_Vx730_PHA_DT_MASK             0x80000000  //bit 31
#define Vx725_Vx730_PHA_EE_MASK             0x40000000  //bit 30
#define Vx725_Vx730_PHA_ET_MASK             0x20000000  //bit 29
#define Vx725_Vx730_PHA_E2_MASK             0x10000000  //bit 28
#define Vx725_Vx730_PHA_ES_MASK             0x08000000  //bit 27
#define Vx725_Vx730_PHA_EX_MASK             0x07000000  //bits 24 to 26 inclusive
#define Vx725_Vx730_PHA_AP1_MASK            0x00C00000  //bits 22 to 23 inclusive
#define Vx725_Vx730_PHA_AP2_MASK            0x00300000  //bits 20 to 21 inclusive
#define Vx725_Vx730_PHA_DP_MASK             0x000F0000  //bits 16 to 19 inclusive
#define Vx725_Vx730_PHA_NSDB8_MASK          0x0000FFFF  //bits 0 to 15 inclusive
//WORD 3
#define Vx725_Vx730_PHA_CH_MASK             0x80000000  //bit 31
#define Vx725_Vx730_PHA_TTT_MASK            0x7FFFFFFF  //bits 0 to 30 inclusive
//WORDS 4 to N-2
#define Vx725_Vx730_PHA_T_1_MASK            0x80080000  //bit 31
#define Vx725_Vx730_PHA_DP_1_MASK           0x40000000  //bit 30
#define Vx725_Vx730_PHA_AP_1_MASK           0x3FFF0000  //bits 16 to 29 inclusive
#define Vx725_Vx730_PHA_T_0_MASK            0x00008000  //bit 15
#define Vx725_Vx730_PHA_DP_0_MASK           0x00004000  //bit 14
#define Vx725_Vx730_PHA_AP_0_MASK           0x00003FFF  //bits 0 to 13 inclusive
//WORD N-1
#define Vx725_Vx730_PHA_EXTRAS2_MASK        0xFFFFFFFF  //bits 0 to 31 inclusive
//WORD N
#define Vx725_Vx730_PHA_EXTRAS_MASK         0x001F0000  //bits 16 to 20 inclusive
#define Vx725_Vx730_PHA_PUR_MASK            0x00008000  //bit 15
#define Vx725_Vx730_PHA_ENERGY_MASK         0x00007FFF  //bits 0 to 14 inclusive


//DPP-PHA Vx725_Vx730
struct Vx725_Vx730_PHA_Data_t {
  //WORD 1 
  uint8_t chagg_header;                                //bit 31
  uint32_t chagg_size;                                 //bits 0 to 31 inclusive 
  //WORD 2
  uint8_t dual_trace;                                  //bit 31
  uint8_t energy_enabled;                              //bit 30
  uint8_t time_enabled;                                //bit 29
  uint8_t extras2_enabled;                             //bit 28
  uint8_t waveform_enabled;                            //bit 27
  uint8_t extras2_option;                              //bits 24 to 26 inclusive
  uint8_t ap1;                                         //bits 22 to 23 inclusive
  uint8_t ap2;                                         //bits 20 to 21 inclusive
  uint8_t dp;                                          //bits 16 to 19 inclusive
  uint16_t nsdb8;                                      //bits 0 to 15 inclusive
  //WORD 3
  uint8_t channel; //bit 31
  uint32_t trigger_time_tag;                           //bits 0 to 30 inclusive
  //WORDS 4 to N-2
  uint16_t analog_probe1[MAX_PROBE_LENGTH];            //bits 0 to 13 inclusive
  uint16_t analog_probe2[MAX_PROBE_LENGTH];            //bits 16 to 29 inclusive
  uint8_t digital_probe1[MAX_PROBE_LENGTH];            //bit 31,15
  uint8_t digital_probe2[MAX_PROBE_LENGTH];            //bit 30,14
  //WORD N-1
  uint32_t extras2;                                    //bits 0 to 31 inclusive
  uint16_t extended_time_stamp;
  uint16_t baseline_timesfour;
  uint16_t flags;
  uint16_t fine_time_stamp;
  uint16_t lost_trigger_counter;
  uint16_t total_trigger_counter;
  uint16_t sample_before_zc;
  uint16_t sample_after_zc;
  //WORD N
  uint8_t extras;                                      //bits 16 to 20 inclusive
  uint8_t pur;                                         //bit 15
  uint16_t energy;                                     //bits 0 to 14 inclusive

  uint16_t individual_chagg_size;                        //size of the individual channels in the aggregate (in 32-bit words)
};


//PSD unpacking
int unpack_vx725_vx730_psd_chagg_header(V1730_ChAgg_Header_t *v1730_chagg_header, Vx725_Vx730_PSD_Data_t *vx725_vx730_psd_data);
int unpack_vx725_vx730_psd_chagg(uint32_t v1730_chagg[], Vx725_Vx730_PSD_Data_t *vx725_vx730_psd_data);

//PHA unpacking
int unpack_vx725_vx730_pha_chagg_header(V1730_ChAgg_Header_t *v1730_chagg_header, Vx725_Vx730_PHA_Data_t *vx725_vx730_pha_data);
int unpack_vx725_vx730_pha_chagg(uint32_t v1730_chagg[], Vx725_Vx730_PHA_Data_t *vx725_vx730_pha_data);





#endif
