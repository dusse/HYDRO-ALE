#include "Logger.hpp"

#include <iostream>
#include <string>
using namespace std;



Logger::Logger()
{
}

void Logger::writeMsg(const char* input, int level)
{
    if (level <= MINIMAL_LEVEL){
        cout << "[" << level << "] " << input <<  "." << endl;
    }
}

void Logger::writeMsg(const char* input)
{
    writeMsg(input, DEBUG);
}



Logger::~Logger()
{
}


