//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  Run_hadd.C             *// 
//*  Last Edit: 02/08/18    *//  
//***************************//

//ROOT Includes
#include "TMath.h"
#include "TRandom.h"

//C/C++ Includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <stdio.h>
#include <vector>

int Run_hadd_Bootstrap() {

  int n_bootstraps = 20;
  int n_runs_per = 25;

  vector<int> block_start;
  vector<int> block_end;
  
  int temp_start;
  int temp_end;
  
  //Read in the list of blocks of runs
  ifstream blocks;
  blocks.open("Bootstrap_Blocks.txt");
  
  if(blocks.is_open()) {
    while(!blocks.eof()) {
      blocks >> temp_start >> temp_end;
      block_start.push_back(temp_start);
      block_end.push_back(temp_end);
    }
  }
  else {
    cout<<"Can't open the blocks file"<<endl;
    return -1;
  }
  
  //Read in the exlcuded files list
  std::vector<int> Exclude_Run;
  std::vector<int> Exclude_Reason;
  
  int exc_run;
  int exc_reason;
  
  //Read in the exlcuded runs file
  ifstream exclude;
  exclude.open("exclude.txt");
  
  if(exclude.is_open()){ 
    while(true) {
      exclude >> exc_run >> exc_reason;
      if(exc_run < 0 || exc_reason < 0) {
	break;
      }
      Exclude_Run.push_back(exc_run);
      Exclude_Reason.push_back(exc_reason);
    }
    
    cout<<"There are a total of "<<Exclude_Run.size()<<" Runs in all of the data to exclude from adding"<<endl;
  }


  vector<int> run_list;
  
  for(int eye=0; eye<block_start.size(); eye++) {
    for(int jay=block_start[eye]; jay<block_end[eye]+1; jay++) {
      bool exclude_this = false;
      for(int kay=0; kay<(int)Exclude_Run.size(); kay++) {
	if(jay==Exclude_Run[kay]){
	  exclude_this = true;
	}
      }
      if(exclude_this) {
	cout<<"Run "<<jay<<" Excluded"<<endl;
      }
      else {
	run_list.push_back(jay);
      }
    }    
  }
  
  cout<<"There are "<<run_list.size()<<" Runs to choose from"<<endl;

  for(int eye=0; eye<run_list.size(); eye++) {
    cout<<eye<<"  "<<run_list[eye]<<endl;
  }
  
  
  int Coincidence_Window = 10;
  int Crystal_Blocking_Time = 2000;
  int DEvent_Blocking_Time = 0;

  string Isotope = "Au197";
  
  ofstream bootstrapped_runs;

  stringstream text;
  stringstream fname;
  
  for(int eye=0; eye<n_bootstraps; eye++) {
    
    //name the output file
    fname.str("");
    fname << "bootstrapped_runlist_"<<eye<<".txt";
    
    //output a list of runs in the rootfile
    bootstrapped_runs.open(fname.str().c_str());
    
    text.str("");
    text <<" xterm -e hadd ";
    //text <<" xterm -hold -e echo ";
    
    //Name the output rootfile
    text << "../stage1_root/Bootsrapped_Sum_Stage1_Histograms_"<<Isotope<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT_"<<eye<<".root ";
    
    //loop over the number of random runs needed
    for(int jay=0; jay<n_runs_per; jay++) {
      
      //Get a random run to add
      int this_run = gRandom->Uniform(0,run_list.size());
      
      //Output the run to the list in the Summed file
      bootstrapped_runs << run_list[this_run]<<"\n";
      
      //Add the random run to the command 
      text << "../stage1_root/Stage1_Histograms_Run_"<<run_list[this_run]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
    }
    
    //run in the background
    text << "&";
    
    //output the strin to screen
    cout<<text.str()<<endl;
    
    //close the output runlist file
    boostrapped_runs.close();
    
    //send the command
    system(text.str().c_str());
  }
  
}

/*

  
  bool force = false;

  int Coincidence_Window = 10;
  int Crystal_Blocking_Time = 2000;
  int DEvent_Blocking_Time = 0;
  
  stringstream text;
  text.str();
  
  text <<"hadd ";
  text << "../stage1_root/Sum_Stage1_Histograms_"<<Isotope<<"_"<<Start<<"_"<<End<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
  
  for(int eye=Start; eye<End+1; eye++) {
    text << "../stage1_root/Stage1_Histograms_Run_"<<eye<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
  }
  
  if(force==true) {
    text <<" -f";
  }
  
  cout<<text.str()<<endl;

  system(text.str().c_str());
}

void Script_hadd() {
  
  //starting time
  struct timeval tv;
  double begin, end, time_elapsed;  // start,stop, elapsed time

  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  
  //Fe57
  run_hadd("Fe57",59236,59714);
  run_hadd("Fe57",59862,60363);
  
  
 
  //Pb208
  run_hadd("Pb208",59715,59842);
  run_hadd("Pb208",60622,60811);
  run_hadd("Pb208",62695,62735);
  // run_hadd("Pb208",63650,64136);
  run_hadd("Pb208",64834,64987);
  run_hadd("Pb208",68040,68103);
 
  
  
  //Au197 4mm
  run_hadd("Au197",60409,60451);
  run_hadd("Au197",60572,60584);
  run_hadd("Au197",63590,63649);
  run_hadd("Au197",68104,68126);
  
 
  //Cu65
  run_hadd("Cu65",60834,61251);
  run_hadd("Cu65",63015,63589);
  run_hadd("Cu65",64137,64692);
  run_hadd("Cu65",64989,65265);
  
  //Cu63
  run_hadd("Cu63",64693,64833);
  
  
  Co59
  run_hadd("Co59_9p9mg",67836,68039);
  run_hadd("Co59_141p0mg",67339,67358);
  run_hadd("Co59_329p4mg",67359,67835);
  run_hadd("Co59_329p4mg",68128,68900);
   

  cout<<"Done!"<<endl;
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin);
  cout<<"Elapsed Time (s): "<<time_elapsed<<endl;

}

int Merge_hadd(string Isotope, int Start, int End, int N) {

  //starting time
  struct timeval tv;
  double begin, end, time_elapsed;  // start,stop, elapsed time

  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  bool force = false;

  int Coincidence_Window = 10;
  int Crystal_Blocking_Time = 2000;
  int DEvent_Blocking_Time = 0;
  
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
      outfilenames[eye] << "../stage1_root/Sum_Stage1_Histograms_"<<Isotope<<"_"<<first_run[eye]<<"_"<<last_run[eye]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT_Excluded.root ";
      
      // command_text[eye].str();
      command_text[eye] << " xterm -e hadd ";
      // command_text[add_counter] << "/usr/bin/xterm -e -hold echo ";
      //  command_text[add_counter] << " echo ";
      //  command_text[eye] << " hadd ";

      command_text[eye] << outfilenames[eye].str()<<" ";

      add_status_fname[eye].str();
      add_status_fname[eye] <<"/tmp/add_status_"<<first_run[eye]<<"_"<<last_run[eye];
      
      add_status[eye] << "ps aux | grep ";
      add_status[eye] << outfilenames[eye].str();
      add_status[eye] << "> "<<add_status_fname[eye].str();


      for(int jay=first_run[eye]; jay<last_run[eye]+1; jay++) {
	bool exclude_this = false;
	for(int kay=0; kay<(int)Exclude_Run.size(); kay++) {
	  if(jay==Exclude_Run[eye]) {
	    exclude_this = true;
	  }
	}
	if(exclude_this) {
	  cout<<"Run "<<jay<<" Excluded"<<endl;
	}
	else {
	  command_text[eye] << "../stage1_root/Stage1_Histograms_Run_"<<jay<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT.root ";
	} 
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
    
    cout<<"Base hadds have finished"<<endl;

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
	outfilenames[add_counter] << "../stage1_root/Sum_Stage1_Histograms_"<<Isotope<<"_"<<first_run[add_counter]<<"_"<<last_run[add_counter]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT_Excluded.root ";
	
	command_text[add_counter].str();
	//	command_text[add_counter] << "echo ";
       	command_text[add_counter] << "xterm -e hadd ";
       	//command_text[add_counter] << "echo ";
	command_text[add_counter] << outfilenames[add_counter].str()<<" ";

	
	add_status_fname[add_counter].str();
	add_status_fname[add_counter] <<"/tmp/add_status_"<<first_run[add_counter]<<"_"<<last_run[add_counter];
	
	add_status[add_counter] << "ps aux | grep ";
	add_status[add_counter] << outfilenames[add_counter].str();
	add_status[add_counter] << "> "<<add_status_fname[add_counter].str();
	
	command_text[add_counter] << "../stage1_root/Sum_Stage1_Histograms_"<<Isotope<<"_"<<first_run[2*jay*level_counter]<<"_"<<last_run[2*jay*level_counter+level_counter-1]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT_Excluded.root ";
	
	command_text[add_counter] << "../stage1_root/Sum_Stage1_Histograms_"<<Isotope<<"_"<<first_run[2*jay*level_counter+level_counter]<<"_"<<last_run[2*(jay+1)*level_counter-1]<<"_"<<Coincidence_Window<<"ns_CW_"<<Crystal_Blocking_Time<<"ns_CBT_"<<DEvent_Blocking_Time<<"ns_DEBT_Excluded.root ";
	
     
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
      
      cout<<"This Level hadds have finished"<<endl;
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


void Launch_Merge_hadd() {

  int exc_run;
  int exc_reason;
  
  //Read in the exlcuded runs file
  ifstream exclude;
  exclude.open("exclude.txt");
  
  if(exclude.is_open()){ 
    while(true) {
      exclude >> exc_run >> exc_reason;
      if(exc_run < 0 || exc_reason < 0) {
	break;
      }
      Exclude_Run.push_back(exc_run);
      Exclude_Reason.push_back(exc_reason);
    }
    
    cout<<"There are a total of "<<Exclude_Run.size()<<" Runs in all of the data to exclude from adding"<<endl;
  }
  
  //Pb208
  //  Merge_hadd("Pb208",59715,59842,4);
  // Merge_hadd("Pb208",60622,60811,4);
  // Merge_hadd("Pb208",62695,62735,4);
  // Merge_hadd("Pb208",63650,64136,4);
  // Merge_hadd("Pb208",64834,64987,4);
  // Merge_hadd("Pb208",68040,68103,4);

  //Au197
  Merge_hadd("Au197",60409,60451,2);
  // Merge_hadd("Au197",60572,60584,2);
  // Merge_hadd("Au197",63590,63649,2);
  // Merge_hadd("Au197",68104,68126,2);

  //Cu65
  // Merge_hadd("Cu65",60834,61251,4);
  // Merge_hadd("Cu65",63015,63589,4);
  // Merge_hadd("Cu65",64137,64692,4);
  // Merge_hadd("Cu65",64989,65265,3);

  //Cu63
  // Merge_hadd("Cu63",64693,64833,3);
}
*/
