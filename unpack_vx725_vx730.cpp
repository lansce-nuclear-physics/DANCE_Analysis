
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
//*  unpack_vx735_vx730.cpp *// 
//*  Last Edit: 07/23/18    *//  
//***************************//

#include "unpack_vx725_vx730.h"


int unpack_vx725_vx730_board_data(V1730_Header_t *v1730_header, Vx725_Vx730_Board_Data_t *vx725_vx730_board_data) {
  
  //WORD 1
  vx725_vx730_board_data->header = (v1730_header->dataword_1 & Vx725_Vx730_HEADER_MASK) >> 28;
  vx725_vx730_board_data->boardaggsize = (v1730_header->dataword_1 & Vx725_Vx730_BOARDAGGSIZE_MASK);
  
  //WORD 2
  vx725_vx730_board_data->boardid = (v1730_header->dataword_2 & Vx725_Vx730_BOARDID_MASK) >> 27;
  vx725_vx730_board_data->pattern = (v1730_header->dataword_2 & Vx725_Vx730_PATTERN_MASK) >> 8;
  vx725_vx730_board_data->channelmask = (v1730_header->dataword_2 & Vx725_Vx730_CHANNEL_MASK);
  
  //WORD 3
  vx725_vx730_board_data->boardaggcounter = (v1730_header->dataword_3 & Vx725_Vx730_BOARDAGGCOUNTER_MASK);
  
  //WORD 4
  vx725_vx730_board_data->boardaggtime = (v1730_header->dataword_4 & Vx725_Vx730_BOARDAGGTIME_MASK);
  
  return 0;
}


int unpack_vx725_vx730_psd_chagg_header(V1730_ChAgg_Header_t *v1730_chagg_header, Vx725_Vx730_PSD_Data_t *vx725_vx730_psd_data) {
  
  //WORD 1
  vx725_vx730_psd_data->chagg_header = (v1730_chagg_header->dataword_1 & Vx725_Vx730_PSD_CHAGGHEADER_MASK) >> 31;
  vx725_vx730_psd_data->chagg_size = (v1730_chagg_header->dataword_1 & Vx725_Vx730_PSD_CHAGGSIZE_MASK);

  //WORD 2
  vx725_vx730_psd_data->dual_trace = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_DT_MASK) >> 31;
  vx725_vx730_psd_data->charge_enabled = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_EQ_MASK) >> 30;
  vx725_vx730_psd_data->time_enabled = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_ET_MASK) >> 29;
  vx725_vx730_psd_data->extras_enabled = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_EE_MASK) >> 28;
  vx725_vx730_psd_data->waveform_enabled = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_ES_MASK) >> 27;
  vx725_vx730_psd_data->extras_option = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_EX_MASK) >> 24;
  vx725_vx730_psd_data->ap = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_AP_MASK) >> 22;
  vx725_vx730_psd_data->dp2 = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_DP2_MASK) >> 19;
  vx725_vx730_psd_data->dp1 = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_DP1_MASK) >> 16;
  vx725_vx730_psd_data->nsdb8 = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PSD_NSDB8_MASK);

  vx725_vx730_psd_data->individual_chagg_size = 4* vx725_vx730_psd_data->nsdb8 +  vx725_vx730_psd_data->extras_enabled + 2;
  return 0;
}


int unpack_vx725_vx730_pha_chagg_header(V1730_ChAgg_Header_t *v1730_chagg_header, Vx725_Vx730_PHA_Data_t *vx725_vx730_pha_data) {

  //WORD 1
  vx725_vx730_pha_data->chagg_header = (v1730_chagg_header->dataword_1 & Vx725_Vx730_PHA_CHAGGHEADER_MASK) >> 31;
  vx725_vx730_pha_data->chagg_size = (v1730_chagg_header->dataword_1 & Vx725_Vx730_PHA_CHAGGSIZE_MASK);
  
  //WORD 2
  vx725_vx730_pha_data->dual_trace = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_DT_MASK) >> 31;
  vx725_vx730_pha_data->energy_enabled = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_EE_MASK) >> 30;
  vx725_vx730_pha_data->time_enabled = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_ET_MASK) >> 29;
  vx725_vx730_pha_data->extras2_enabled = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_E2_MASK) >> 28;
  vx725_vx730_pha_data->waveform_enabled = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_ES_MASK) >> 27;
  vx725_vx730_pha_data->extras2_option = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_EX_MASK) >> 24;
  vx725_vx730_pha_data->ap1 = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_AP1_MASK) >> 22;
  vx725_vx730_pha_data->ap2 = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_AP2_MASK) >> 20;
  vx725_vx730_pha_data->dp = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_DP_MASK) >> 16;
  vx725_vx730_pha_data->nsdb8 = (v1730_chagg_header->dataword_2 & Vx725_Vx730_PHA_NSDB8_MASK);
  
  vx725_vx730_pha_data->individual_chagg_size = 4* vx725_vx730_pha_data->nsdb8 +  vx725_vx730_pha_data->extras2_enabled + 2;
  return 0;
}


int unpack_vx725_vx730_psd_chagg(uint32_t *v1730_chagg_data, Vx725_Vx730_PSD_Data_t *vx725_vx730_psd_data) {
  
  int word_counter=0;
  
  //WORD 1 (TTT)
  vx725_vx730_psd_data->channel = (v1730_chagg_data[word_counter] & Vx725_Vx730_PSD_CH_MASK) >> 31;
  vx725_vx730_psd_data->trigger_time_tag = (v1730_chagg_data[word_counter] & Vx725_Vx730_PSD_TTT_MASK);
  word_counter++;
  
  //WORD 2 to N-2 (or N-1 if no extras)  (Analog and Digital Probes)
  if(vx725_vx730_psd_data->waveform_enabled) {
    
    //Determine how many 32-bit words are there
    uint32_t probe_words = 4.0*vx725_vx730_psd_data->nsdb8; 
    
    //Dual trace off
    if(vx725_vx730_psd_data->dual_trace == 0) {
      for(uint32_t kay=0; kay<probe_words; kay++) {
	uint32_t probe_data = v1730_chagg_data[word_counter];
	word_counter++;
	
	//analog probes
	vx725_vx730_psd_data->analog_probe1[(int)(2*kay)] = (probe_data & Vx725_Vx730_PSD_AP_0_MASK);
	vx725_vx730_psd_data->analog_probe2[(int)(2*kay)] = 0;
	vx725_vx730_psd_data->analog_probe1[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PSD_AP_1_MASK) >> 16;
	vx725_vx730_psd_data->analog_probe2[(int)(2*kay+1)] = 0;
	
	//digital probes
	vx725_vx730_psd_data->digital_probe1[(int)(2*kay)] = (probe_data & Vx725_Vx730_PSD_DP1_0_MASK) >> 14;
	vx725_vx730_psd_data->digital_probe2[(int)(2*kay)] = (probe_data & Vx725_Vx730_PSD_DP2_0_MASK) >> 15;
	vx725_vx730_psd_data->digital_probe1[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PSD_DP1_1_MASK) >> 30;
	vx725_vx730_psd_data->digital_probe2[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PSD_DP2_1_MASK) >> 31;
	
      } //End loop over probe words
    } //End dual trace off
    
    //Dual trace on
    else if(vx725_vx730_psd_data->dual_trace == 1) {
      for(uint32_t kay=0; kay<probe_words; kay++) {
	uint32_t probe_data = v1730_chagg_data[word_counter];
	word_counter++;
	
	//analog probes
	vx725_vx730_psd_data->analog_probe1[kay] = (probe_data & Vx725_Vx730_PSD_AP_0_MASK);
	vx725_vx730_psd_data->analog_probe2[kay] = (probe_data & Vx725_Vx730_PSD_AP_1_MASK) >> 16;
	
	//digital probes
	vx725_vx730_psd_data->digital_probe1[(int)(2*kay)] = (probe_data & Vx725_Vx730_PSD_DP1_0_MASK) >> 14;
	vx725_vx730_psd_data->digital_probe2[(int)(2*kay)] = (probe_data & Vx725_Vx730_PSD_DP2_0_MASK) >> 15;
	vx725_vx730_psd_data->digital_probe1[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PSD_DP1_1_MASK) >> 30;
	vx725_vx730_psd_data->digital_probe2[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PSD_DP2_1_MASK) >> 31;
	
      } //End loop over probe words
    } //End dual trace on
  
    //Dual trace is messed up and need to exit
    else {
      return -1;
    }
  } //End of check on waveforms
  
  //WORD N-1 (or N if no extras) (Extras)
  if(vx725_vx730_psd_data->extras_enabled) {
      
    vx725_vx730_psd_data->extras = v1730_chagg_data[word_counter];
    word_counter++;
      
    //Decode the extras
    switch(vx725_vx730_psd_data->extras_option) {
    case 0:
      vx725_vx730_psd_data->extended_time_stamp = (vx725_vx730_psd_data->extras & 0xFFFF0000) >> 16;
      vx725_vx730_psd_data->baseline_timesfour = (vx725_vx730_psd_data->extras & 0x0000FFFF);
      vx725_vx730_psd_data->flags = 0;
      vx725_vx730_psd_data->fine_time_stamp = 0;
      vx725_vx730_psd_data->lost_trigger_counter = 0;
      vx725_vx730_psd_data->total_trigger_counter = 0;
      vx725_vx730_psd_data->cfd_sazc = 0;
      vx725_vx730_psd_data->cfd_sbzc = 0;	
      break;
    case 1:
      vx725_vx730_psd_data->extended_time_stamp = (vx725_vx730_psd_data->extras & 0xFFFF0000) >> 16;
      vx725_vx730_psd_data->baseline_timesfour = 0;
      vx725_vx730_psd_data->flags = (vx725_vx730_psd_data->extras & 0x0000FFFF);
      vx725_vx730_psd_data->fine_time_stamp =0;
      vx725_vx730_psd_data->lost_trigger_counter = 0;
      vx725_vx730_psd_data->total_trigger_counter = 0;
      vx725_vx730_psd_data->cfd_sazc = 0;
      vx725_vx730_psd_data->cfd_sbzc = 0;
      break;
    case 2:
      vx725_vx730_psd_data->extended_time_stamp = (vx725_vx730_psd_data->extras & 0xFFFF0000) >> 16;
      vx725_vx730_psd_data->baseline_timesfour = 0;
      vx725_vx730_psd_data->flags = (vx725_vx730_psd_data->extras & 0x0000FC00) >> 10;
      vx725_vx730_psd_data->fine_time_stamp = (vx725_vx730_psd_data->extras & 0x000003FF);
      vx725_vx730_psd_data->lost_trigger_counter = 0;
      vx725_vx730_psd_data->total_trigger_counter = 0;
      vx725_vx730_psd_data->cfd_sazc = 0;
      vx725_vx730_psd_data->cfd_sbzc = 0;
      break;
    case 4:
      vx725_vx730_psd_data->extended_time_stamp = 0;
      vx725_vx730_psd_data->baseline_timesfour = 0;
      vx725_vx730_psd_data->flags = 0;
      vx725_vx730_psd_data->fine_time_stamp = 0;
      vx725_vx730_psd_data->lost_trigger_counter = (vx725_vx730_psd_data->extras & 0xFFFF0000) >> 16;
      vx725_vx730_psd_data->total_trigger_counter = (vx725_vx730_psd_data->extras & 0x0000FFFF);
      vx725_vx730_psd_data->cfd_sazc = 0;
      vx725_vx730_psd_data->cfd_sbzc = 0;
    case 5:
      vx725_vx730_psd_data->extended_time_stamp = 0;
      vx725_vx730_psd_data->baseline_timesfour = 0;
      vx725_vx730_psd_data->flags = 0;
      vx725_vx730_psd_data->fine_time_stamp = 0;
      vx725_vx730_psd_data->lost_trigger_counter = 0;
      vx725_vx730_psd_data->total_trigger_counter = 0;
      vx725_vx730_psd_data->cfd_sazc = (vx725_vx730_psd_data->extras & 0xFFFF0000) >> 16;
      vx725_vx730_psd_data->cfd_sbzc = (vx725_vx730_psd_data->extras & 0x0000FFFF);
    default:
      vx725_vx730_psd_data->extended_time_stamp = 0;
      vx725_vx730_psd_data->baseline_timesfour = 0;
      vx725_vx730_psd_data->flags = 0;
      vx725_vx730_psd_data->lost_trigger_counter = 0;
      vx725_vx730_psd_data->total_trigger_counter = 0;
      vx725_vx730_psd_data->cfd_sazc = 0;
      vx725_vx730_psd_data->cfd_sbzc = 0;
      vx725_vx730_psd_data->fine_time_stamp = 0;
      break;   
    }
  } //End of check on extras enabled
    
  //WORD N (PSD islow,pur,ifast)
  vx725_vx730_psd_data->qlong = ( v1730_chagg_data[word_counter] & Vx725_Vx730_PSD_QLONG_MASK) >> 16;
  vx725_vx730_psd_data->pur = ( v1730_chagg_data[word_counter] & Vx725_Vx730_PSD_PUR_MASK) >> 15;
  vx725_vx730_psd_data->qshort = ( v1730_chagg_data[word_counter] & Vx725_Vx730_PSD_QSHORT_MASK);
  word_counter++;
    
  return word_counter;
}


int unpack_vx725_vx730_pha_chagg(uint32_t *v1730_chagg_data, Vx725_Vx730_PHA_Data_t *vx725_vx730_pha_data) {
  
  int word_counter=0;
  
  //WORD 1 (TTT)
  vx725_vx730_pha_data->channel = (v1730_chagg_data[word_counter] & Vx725_Vx730_PHA_CH_MASK) >> 31;
  vx725_vx730_pha_data->trigger_time_tag = (v1730_chagg_data[word_counter] & Vx725_Vx730_PHA_TTT_MASK);
  word_counter++;
  
  //WORD 2 to N-2 (or N-1 if no extras)  (Analog and Digital Probes)
  if(vx725_vx730_pha_data->waveform_enabled) {
    
    //Determine how many 32-bit words are there
    uint32_t probe_words = 4.0*vx725_vx730_pha_data->nsdb8; 
    
    //Dual trace off
    if(vx725_vx730_pha_data->dual_trace == 0) {
      for(uint32_t kay=0; kay<probe_words; kay++) {
	uint32_t probe_data = v1730_chagg_data[word_counter];
	word_counter++;
	
	//analog probes
	vx725_vx730_pha_data->analog_probe1[(int)(2*kay)] = (probe_data & Vx725_Vx730_PHA_AP_0_MASK);
	vx725_vx730_pha_data->analog_probe2[(int)(2*kay)] = 0;
	vx725_vx730_pha_data->analog_probe1[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PHA_AP_1_MASK) >> 16;
	vx725_vx730_pha_data->analog_probe2[(int)(2*kay+1)] = 0;
	
	//digital probes
	vx725_vx730_pha_data->digital_probe1[(int)(2*kay)] = (probe_data & Vx725_Vx730_PHA_T_0_MASK) >> 14;
	vx725_vx730_pha_data->digital_probe2[(int)(2*kay)] = (probe_data & Vx725_Vx730_PHA_DP_0_MASK) >> 15;
	vx725_vx730_pha_data->digital_probe1[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PHA_T_1_MASK) >> 30;
	vx725_vx730_pha_data->digital_probe2[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PHA_DP_1_MASK) >> 31;
	
      } //End loop over probe words
    } //End dual trace off
    
    //Dual trace on
    else if(vx725_vx730_pha_data->dual_trace == 1) {
      for(uint32_t kay=0; kay<probe_words; kay++) {
	uint32_t probe_data = v1730_chagg_data[word_counter];
	word_counter++;
	
	//analog probes
	vx725_vx730_pha_data->analog_probe1[kay] = (probe_data & Vx725_Vx730_PHA_AP_0_MASK);
	vx725_vx730_pha_data->analog_probe2[kay] = (probe_data & Vx725_Vx730_PHA_AP_1_MASK) >> 16;
	
	//digital probes
	vx725_vx730_pha_data->digital_probe1[(int)(2*kay)] = (probe_data & Vx725_Vx730_PHA_T_0_MASK) >> 14;
	vx725_vx730_pha_data->digital_probe2[(int)(2*kay)] = (probe_data & Vx725_Vx730_PHA_DP_0_MASK) >> 15;
	vx725_vx730_pha_data->digital_probe1[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PHA_T_1_MASK) >> 30;
	vx725_vx730_pha_data->digital_probe2[(int)(2*kay+1)] = (probe_data & Vx725_Vx730_PHA_DP_1_MASK) >> 31;
	
      } //End loop over probe words
    } //End dual trace on
  
    //Dual trace is messed up and need to exit
    else {
      return -1;
    }
  } //End of check on waveforms
  
  //WORD N-1 (or N if no extras) (Extras)
  if(vx725_vx730_pha_data->extras2_enabled) {
      
    vx725_vx730_pha_data->extras2 = v1730_chagg_data[word_counter];
    word_counter++;
      
    //Decode the extras
    switch(vx725_vx730_pha_data->extras2_option) {
      
    case 0:
      vx725_vx730_pha_data->extended_time_stamp = (vx725_vx730_pha_data->extras2 & 0xFFFF0000) >> 16;
      vx725_vx730_pha_data->baseline_timesfour = (vx725_vx730_pha_data->extras2 & 0x0000FFFF);
      vx725_vx730_pha_data->lost_trigger_counter = 0;
      vx725_vx730_pha_data->total_trigger_counter = 0;
      vx725_vx730_pha_data->fine_time_stamp = 0;
      vx725_vx730_pha_data->sample_before_zc = 0;
      vx725_vx730_pha_data->sample_after_zc = 0;
      break;
    case 2:
      vx725_vx730_pha_data->extended_time_stamp = (vx725_vx730_pha_data->extras2 & 0xFFFF0000) >> 16;
      vx725_vx730_pha_data->baseline_timesfour = 0;
      vx725_vx730_pha_data->lost_trigger_counter = 0;
      vx725_vx730_pha_data->total_trigger_counter = 0;
      vx725_vx730_pha_data->fine_time_stamp = (vx725_vx730_pha_data->extras2 & 0x0000FFFF);
      vx725_vx730_pha_data->sample_before_zc = 0;
      vx725_vx730_pha_data->sample_after_zc = 0;
      break;
    case 4:
      vx725_vx730_pha_data->extended_time_stamp = 0;
      vx725_vx730_pha_data->baseline_timesfour = 0;
      vx725_vx730_pha_data->lost_trigger_counter = (vx725_vx730_pha_data->extras2 & 0xFFFF0000) >> 16;
      vx725_vx730_pha_data->total_trigger_counter = (vx725_vx730_pha_data->extras2 & 0x0000FFFF);
      vx725_vx730_pha_data->fine_time_stamp = 0;
      vx725_vx730_pha_data->sample_before_zc = 0;
      vx725_vx730_pha_data->sample_after_zc = 0;
    case 5:
      vx725_vx730_pha_data->extended_time_stamp = 0;
      vx725_vx730_pha_data->baseline_timesfour = 0;
      vx725_vx730_pha_data->lost_trigger_counter = 0;
      vx725_vx730_pha_data->total_trigger_counter = 0;
      vx725_vx730_pha_data->fine_time_stamp = 0;
      vx725_vx730_pha_data->sample_before_zc = (vx725_vx730_pha_data->extras2 & 0xFFFF0000) >> 16;
      vx725_vx730_pha_data->sample_after_zc = (vx725_vx730_pha_data->extras2 & 0x0000FFFF);
    default:
      vx725_vx730_pha_data->extended_time_stamp = 0;
      vx725_vx730_pha_data->baseline_timesfour = 0;
      vx725_vx730_pha_data->lost_trigger_counter = 0;
      vx725_vx730_pha_data->total_trigger_counter = 0;
      vx725_vx730_pha_data->fine_time_stamp = 0;
      vx725_vx730_pha_data->sample_before_zc = 0;
      vx725_vx730_pha_data->sample_after_zc = 0;
      break;
    }
  } //End of check on extras enabled
    
  //WORD N (PHA extras,pur,energy)	
  vx725_vx730_pha_data->extras = (v1730_chagg_data[word_counter] & Vx725_Vx730_PHA_EXTRAS_MASK) >> 16;
  vx725_vx730_pha_data->pur = (v1730_chagg_data[word_counter] & Vx725_Vx730_PHA_PUR_MASK) >> 15;
  vx725_vx730_pha_data->energy = (v1730_chagg_data[word_counter] & Vx725_Vx730_PHA_ENERGY_MASK);
  word_counter++;
  
  return word_counter;
}
