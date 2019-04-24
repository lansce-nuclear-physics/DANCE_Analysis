//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  message.h              *// 
//*  Last Edit: 04/19/19    *//  
//***************************//

#ifndef MESSAGE_H
#define MESSAGE_H

#include <iostream>
#include <string>

//COLORS
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
#define CLEAR "\033[2J"  // clear screen escape code 


void DANCE_Error(std::string program, std::string message);
void DANCE_Info(std::string program, std::string message);
void DANCE_Success(std::string program, std::string message);
void DANCE_Init(std::string program, std::string message);


#endif
