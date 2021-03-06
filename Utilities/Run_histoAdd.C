//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  Run_histoAdd.C         *// 
//*  Last Edit: 01/25/18    *//  
//***************************//

//ROOT Includes
#include "TMath.h"

//C/C++ Includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <stdio.h>

void Run_histoAdd(string Isotope, int Start, int End) {

  bool force = true;

  int Coincidence_Window = 10;
  int Crystal_Blocking_Time = 2000;
  int DEvent_Blocking_Time = 500;
  
  stringstream text;
  text.str();
  
  text <<"histoAdd ";
  text << "Sum_Stage1_Histograms_"<<Isotope<<"_"<<Start<<"_"<<End<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
  
  for(int eye=Start; eye<End+1; eye++) {
    text << "Stage1_Histograms_Run_"<<eye<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
  }
  
  if(force==true) {
    text <<" -f";
  }
  
  cout<<text.str()<<endl;

  system(text.str().c_str());
}

void Script_histoAdd() {
  
  //starting time
  struct timeval tv;
  double begin, end, time_elapsed;  // start,stop, elapsed time

  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  /*
  //Fe57
  Run_histoAdd("Fe57",59236,59714);
  Run_histoAdd("Fe57",59862,60363);
  */
  
 
  //Pb208
  Run_histoAdd("Pb208",59715,59842);
  Run_histoAdd("Pb208",60622,60811);
  Run_histoAdd("Pb208",62695,62735);
  // Run_histoAdd("Pb208",63650,64136);
  Run_histoAdd("Pb208",64834,64987);
  Run_histoAdd("Pb208",68040,68103);
 
  
  /*
  //Au197 4mm
  Run_histoAdd("Au197",60409,60451);
  Run_histoAdd("Au197",60572,60584);
  Run_histoAdd("Au197",63590,63649);
  Run_histoAdd("Au197",68104,68126);
  */
 
  //Cu65
  Run_histoAdd("Cu65",60834,61251);
  Run_histoAdd("Cu65",63015,63589);
  Run_histoAdd("Cu65",64137,64692);
  Run_histoAdd("Cu65",64989,65265);
  
  //Cu63
  Run_histoAdd("Cu63",64693,64833);
  
  /*
  //Co59
  Run_histoAdd("Co59_9p9mg",67836,68039);
  Run_histoAdd("Co59_141p0mg",67339,67358);
  Run_histoAdd("Co59_329p4mg",67359,67835);
  Run_histoAdd("Co59_329p4mg",68128,68900);
  */ 

  cout<<"Done!"<<endl;
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin);
  cout<<"Elapsed Time (s): "<<time_elapsed<<endl;

}

int Merge_histoAdd(string Isotope, int Start, int End, int N) {

  //starting time
  struct timeval tv;
  double begin, end, time_elapsed;  // start,stop, elapsed time

  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  bool force = true;

  int Coincidence_Window = 10;
  int Crystal_Blocking_Time = 2000;
  int DEvent_Blocking_Time = 500;
  
  int nruns = End-Start+1;

  if(N<=0) {
    cout<<"Its hard to accomplish a task when no resources are given..."<<endl;
    N=1;
    cout<<"N is now: "<<N<<endl;
  }

  if(N > 5) {
    cout<<"There will be more base level adds than cores..."<<endl;
    N=5;
    cout<<"N is now: "<<N<<endl;
  }
  else {
    //Check to see if the requested histoadd depth is reasonable
    int base_level_adds = TMath::Power(2,N);
    int runs_per_bladd = (int)nruns/base_level_adds;

    //Loop until reasonable
    bool check_add_levels=true;
    
    while(check_add_levels) {
      if(runs_per_bladd < 2.0) {
	cout<<"Warning there are less than two runs per first level histoadd..."<<endl;
	cout<<"Reducing add level by 1"<<endl;
	N--;
	cout<<"N is now: "<<N<<endl;
	base_level_adds = TMath::Power(2,N);
	runs_per_bladd = (int)nruns/base_level_adds;
      }
      else {
	check_add_levels=false;
      }
    }

    cout<<"Runs per BL Add: "<<runs_per_bladd<<endl;

    //Number of base level adds needing to be done
    int nbladds = (int)TMath::Power(2,N);
     
    //Runs per Base Level Add
    runs_per_bladd = (int)nruns/nbladds;
    
    stringstream command_text[100];
    stringstream outfilenames[100];
    
    int first_run[100];
    int last_run[100];

    stringstream add_status[100];
    stringstream add_status_fname[100];

    int add_counter=0;
    int previous_add_counter=0;
    int ret=-1;

    int nblcounter=0;
    //deal with the first level first
    for(int eye=0; eye<nbladds; eye++) {
      first_run[eye] = Start + nblcounter*runs_per_bladd+1;
      last_run[eye] = Start + (nblcounter+1)*runs_per_bladd;
      nblcounter++;
      
    }
    //make sure the first and last run is correct
    first_run[0]=Start;
    last_run[nbladds-1]=End;

    //Do the base level add
    for(int eye=0; eye<nbladds; eye++) {
      outfilenames[eye].str();
      outfilenames[eye] << "Sum_Stage1_Histograms_"<<Isotope<<"_"<<first_run[eye]<<"_"<<last_run[eye]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
      
      command_text[eye].str();
      command_text[eye] << " xterm -e histoAdd ";
      command_text[eye] << outfilenames[eye].str()<<" ";

      add_status_fname[eye].str();
      add_status_fname[eye] <<"/tmp/add_status_"<<first_run[eye]<<"_"<<last_run[eye];

      add_status[eye] << "ps aux | grep ";
      add_status[eye] << outfilenames[eye].str();
      add_status[eye] << "> "<<add_status_fname[eye].str();


      for(int jay=first_run[eye]; jay<last_run[eye]+1; jay++) {
	command_text[eye] << "Stage1_Histograms_Run_"<<jay<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
      }      
      if(force==true) {
	command_text[eye] <<" -f";
      }
      command_text[eye] <<" &";
      
      ret = system(command_text[eye].str().c_str());
      add_counter++;
    }    
    
    int nrunning=0;
    std::ifstream fstatus;
    while(true) {
      nrunning = 0;
      for(int eye=0; eye<add_counter; eye++) {

	system(add_status[eye].str().c_str());
	fstatus.open(add_status_fname[eye].str().c_str());
	
	int numLines=0;
	string unused;
	while ( std::getline(fstatus, unused) )
	  ++numLines;
	
	if(numLines > 2) {
	  nrunning++;
	}
	fstatus.close();
      }
      
      if(nrunning==0) {
	break;
      }
      sleep(1);
    }
    
    cout<<"Base histoAdds have finished"<<endl;

    gettimeofday(&tv,NULL);  
    end=tv.tv_sec+(tv.tv_usec/1000000.0);
    time_elapsed = (double)(end-begin);
    cout<<"Elapsed Time (s): "<<time_elapsed<<endl;

    //where we leave off at each level
    previous_add_counter=add_counter;

    int level_counter=1;
    for(int eye=N-1; eye>=0; eye--) {
      cout<<TMath::Power(2,eye)<<"  Level Counter: "<<level_counter<<endl;
      
      for(int jay=0; jay<TMath::Power(2,eye); jay++) {	
	first_run[add_counter] = first_run[2*jay*level_counter];
	last_run[add_counter] = last_run[2*(jay+1)*level_counter-1];
	
	outfilenames[add_counter].str();
	outfilenames[add_counter] << "Sum_Stage1_Histograms_"<<Isotope<<"_"<<first_run[add_counter]<<"_"<<last_run[add_counter]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
	
	command_text[add_counter].str();
	command_text[add_counter] << "xterm -e histoAdd ";
	//command_text[add_counter] << "xterm -e -hold echo ";
	command_text[add_counter] << outfilenames[add_counter].str()<<" ";
	
	add_status_fname[add_counter].str();
	add_status_fname[add_counter] <<"/tmp/add_status_"<<first_run[add_counter]<<"_"<<last_run[add_counter];
	
	add_status[add_counter] << "ps aux | grep ";
	add_status[add_counter] << outfilenames[add_counter].str();
	add_status[add_counter] << "> "<<add_status_fname[add_counter].str();
	
	command_text[add_counter] << "Sum_Stage1_Histograms_"<<Isotope<<"_"<<first_run[2*jay*level_counter]<<"_"<<last_run[2*jay*level_counter+level_counter-1]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
	
	command_text[add_counter] << "Sum_Stage1_Histograms_"<<Isotope<<"_"<<first_run[2*jay*level_counter+level_counter]<<"_"<<last_run[2*(jay+1)*level_counter-1]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
	
     
	if(force==true) {
	  command_text[add_counter] <<" -f";
	}
	command_text[add_counter] <<" &";
	
       	ret = system(command_text[add_counter].str().c_str());
	

	add_counter++;
      }

      //check to see whats running before proceeding
     
      nrunning=0;
      while(true) {
	nrunning = 0;
	for(int kay=previous_add_counter; kay<add_counter; kay++) {
	  
	  system(add_status[kay].str().c_str());
	  fstatus.open(add_status_fname[kay].str().c_str());
	  
	  int numLines=0;
	  string unused;
	  while ( std::getline(fstatus, unused) )
	    ++numLines;
	  
	  if(numLines > 2) {
	    nrunning++;
	  }
	  fstatus.close();
	}
       
	if(nrunning==0) {
	  break;
	}
	sleep(1);
      }
      
      cout<<"This Level histoAdds have finished"<<endl;
      gettimeofday(&tv,NULL);  
      end=tv.tv_sec+(tv.tv_usec/1000000.0);
      time_elapsed = (double)(end-begin);
      cout<<"Elapsed Time (s): "<<time_elapsed<<endl;
    
      level_counter+=level_counter;

      //where we leave off at each level
      previous_add_counter=add_counter;
    }
  
    cout<<"Cleaning Up:"<<endl;
    
    stringstream rm_tmp[100];
    stringstream rm_root[100];
    
    
    for(int eye=0; eye<add_counter; eye++) {
      rm_tmp[eye].str();
      rm_tmp[eye] << "rm -f "<<add_status_fname[eye].str();
      system(rm_tmp[eye].str().c_str());
      
      if(eye < (add_counter-1)) {
	rm_root[eye].str();
	rm_root[eye] << "rm -f "<<outfilenames[eye].str();
	system(rm_root[eye].str().c_str());
	cout<<"Removing: "<<rm_root[eye].str()<<endl;
      }
    }
  }
  
  cout<<"Done!"<<endl;
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin);
  cout<<"Elapsed Time (s): "<<time_elapsed<<endl;
  
  return 0;

}


void Launch_Merge_histoAdd() {

  //Pb208
  // Merge_histoAdd("Pb208",59715,59842,4);
  // Merge_histoAdd("Pb208",60622,60811,4);
  // Merge_histoAdd("Pb208",62695,62735,4);
  Merge_histoAdd("Pb208",63650,64136,4);
  // Merge_histoAdd("Pb208",64834,64987,4);
  // Merge_histoAdd("Pb208",68040,68103,4);


  //Au197
  //  Merge_histoAdd("Au197",60409,60451,2);
  //  Merge_histoAdd("Au197",60572,60584,2);
  //  Merge_histoAdd("Au197",63590,63649,2);
  //  Merge_histoAdd("Au197",68104,68126,2);

  //Cu65
  //Merge_histoAdd("Cu65",60834,61251,4);
  //Merge_histoAdd("Cu65",63015,63589,4);
  //Merge_histoAdd("Cu65",64137,64692,4);
  //Merge_histoAdd("Cu65",64989,65265,3);

  //Cu63
  //Merge_histoAdd("Cu63",64693,64833,3);
}
