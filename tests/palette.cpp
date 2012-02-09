/*! \file palette.cpp
 *  \brief Test palette.
 *
 *  \test Test palette.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "palette.hpp"
#include "error.hpp"


using namespace std;


void check( const Color &c, double r, double g, double b)
{
    if( fabs( c[0] - r ) > 0.01 ) {
	std::cout << "Error with color\n";
	exit( 1 );
    }
    if( fabs( c[1] - g ) > 0.01 ) {
	std::cout << "Error with color\n";
	exit( 1 );
    }
    if( fabs( c[2] - b ) > 0.01 ) {
	std::cout << "Error with color\n";
	exit( 1 );
    }
}


void test( int argc, char **argv )
{
    std::vector<Palette::Entry> entries;
    entries.push_back( Palette::Entry( Color(1,0,0), 2 ) );
    entries.push_back( Palette::Entry( Color(0,1,0), 3 ) );
    entries.push_back( Palette::Entry( Color(1,0,1), 5 ) );
    Palette palette( entries );
    //palette.debug_print();

    Color c;
    c = palette(0.0);
    check( c, 1, 0, 0 );

    c = palette(0.3333333);
    check( c, 0, 1, 0 );

    c = palette(1.0);
    check( c, 1, 0, 1 );

    palette.clear();
    palette.push_back( Color(1,0,0), 0 );
    palette.push_back( Color(0,0,1), 10 );
    palette.normalize();
    //palette.debug_print();

    c = palette(0.0);
    check( c, 1, 0, 0 );

    c = palette(0.5);
    check( c, 0.5, 0, 0.5 );

    c = palette(1.0);
    check( c, 0, 0, 1 );
}


