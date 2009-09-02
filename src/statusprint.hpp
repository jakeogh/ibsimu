/*! \file statusprint.hpp
 *  \brief Header file for printing status.
 */


#ifndef STATUSPRINT_HPP
#define STATUSPRINT_HPP 1


#include <time.h>


class StatusPrint {
    
    size_t _width;
    time_t _time;

public:

    StatusPrint();

    ~StatusPrint();

    void print( const std::string &str );
    
};

#endif













