//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  message.cpp            *// 
//*  Last Edit: 04/19/19    *//  
//***************************//

#include "message.h"

using namespace std;

void DANCE_Error(string program, string message) {
  cout<<RED<<program<<" [ERROR] "<<message<<RESET<<endl;
}

void DANCE_Info(string program, string message) {
  cout<<program<<" [INFO] "<<message<<endl;
}

void DANCE_Success(string program, string message) {
  cout<<GREEN<<program<<" [SUCCESS] "<<message<<RESET<<endl;
}

void DANCE_Init(string program, string message) {
  cout<<BLUE<<program<<" [INIT] "<<message<<RESET<<endl;
}
