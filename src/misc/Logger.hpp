//
//  Logger.hpp

#ifndef Logger_hpp
#define Logger_hpp

#include <stdio.h>

#define DEBUG    2
#define INFO     1
#define CRITICAL 0

#define MINIMAL_LEVEL 1


class Logger
{
public:
    Logger();
    void writeMsg(const char* input);
    void writeMsg(const char* input, int level);
    ~Logger();
};
#endif
