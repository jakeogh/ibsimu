/*! \file plasmacyl.cpp
 *  \brief Test quasi random number generator.
 *
 *  \test Test quasi random number generator.
 */


#include <iostream>
#include <fstream>
#include <iomanip>
#include "qrandom.hpp"
#include "error.hpp"


using namespace std;


void test( int argc, char **argv )
{
    const int N = 10000;
    double sd = 0.0;
    double x[3];
    bool gaussian[3] = {false, false, true};
    QRandom qrng( 3 );

    ofstream fstr( "qrandom.dat" );
    for( int a = 0; a < N; a++ ) { 
	qrng.get_part_gaussian( gaussian, x );
	sd += x[2]*x[2];
	fstr << setw(12) << x[0] << " " 
	     << setw(12) << x[1] << " "
	     << setw(12) << x[2] << "\n";
    }
    sd = sqrt(sd/N);
    fstr << "# sd[3] = " << sd << "\n";
}


