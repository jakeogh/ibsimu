#include <iostream>
#include <iomanip>
#include "statusprint.hpp"


StatusPrint::StatusPrint()
    : _width(0)
{
    _time = time( NULL );
}


StatusPrint::~StatusPrint()
{
    for( size_t i = 0; i < _width; i++ )
	std::cout << "\b";
}


void StatusPrint::print( const std::string &str )
{
    time_t t = time( NULL );
    if( t != _time ) {
	_width = str.length();
	_time = t;
	for( size_t i = 0; i < _width; i++ )
	    std::cout << "\b";
	std::cout << str << std::flush;
    }
}















