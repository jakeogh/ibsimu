/*! \file plasmacyl.cpp
 *  \brief Test quasi random number generator.
 *
 *  \test Test quasi random number generator.
 */


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "random.hpp"
#include "error.hpp"
#include "ibsimutest.hpp"


using namespace std;


void test( int argc, char **argv )
{
    const int N = 10000;
    std::vector<double> x_vec[5];
    double x_ave[5] = {0.0, 0.0, 0.0};
    double x_sd[5] = {0.0, 0.0, 0.0};
    double x[5];

    QRandom qrng( 5 );
    qrng.set_transformation( 2, Gaussian_Transformation() );
    qrng.set_transformation( 3, Cosine_Transformation() );
    qrng.set_transformation( 4, Gamma_Transformation( 9, 0.5 ) );

    // Produce data and calculate average
    for( int a = 0; a < N; a++ ) { 

	qrng.get( x );

	x_vec[0].push_back( x[0] );
	x_vec[1].push_back( x[1] );
	x_vec[2].push_back( x[2] );
	x_vec[3].push_back( x[3] );
	x_vec[4].push_back( x[4] );

	x_ave[0] += x[0];
	x_ave[1] += x[1];
	x_ave[2] += x[2];
	x_ave[3] += x[3];
	x_ave[4] += x[4];
    }

    x_ave[0] = x_ave[0]/N;
    x_ave[1] = x_ave[1]/N;
    x_ave[2] = x_ave[2]/N;
    x_ave[3] = x_ave[3]/N;
    x_ave[4] = x_ave[4]/N;

    // Calculate standard deviation
    for( int a = 0; a < N; a++ ) { 

	double t;
	t = x_vec[0][a] - x_ave[0];
	x_sd[0] += t*t;
	t = x_vec[1][a] - x_ave[1];
	x_sd[1] += t*t;
	t = x_vec[2][a] - x_ave[2];
	x_sd[2] += t*t;
	t = x_vec[3][a] - x_ave[3];
	x_sd[3] += t*t;
	t = x_vec[4][a] - x_ave[4];
	x_sd[4] += t*t;
    }

    x_sd[0] = sqrt( x_sd[0]/N );
    x_sd[1] = sqrt( x_sd[1]/N );
    x_sd[2] = sqrt( x_sd[2]/N );
    x_sd[3] = sqrt( x_sd[3]/N );
    x_sd[4] = sqrt( x_sd[4]/N );

    // Checks
    if( fabs( x_ave[0] - 0.5 ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Linear RNG 1 failed: average = " + to_string(x_ave[0]) ) );
    if( fabs( x_ave[1] - 0.5 ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Linear RNG 2 failed: average = " + to_string(x_ave[1]) ) );
    if( fabs( x_ave[2] - 0.0 ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Gaussian RNG failed: average = " + to_string(x_ave[2]) ) );
    if( fabs( x_ave[3] - 0.0 ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Cosine RNG failed: average = " + to_string(x_ave[3]) ) );
    if( fabs( x_ave[4] - 9.0*0.5 ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Gamma(9,0.5) RNG failed: average = " + to_string(x_ave[4]) ) );

    if( fabs( x_sd[0] - sqrt(1.0/12.0) ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Linear RNG 1 failed: standard deviation = " + to_string(x_sd[0]) ) );
    if( fabs( x_sd[1] - sqrt(1.0/12.0) ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Linear RNG 2 failed: standard deviation = " + to_string(x_sd[1]) ) );
    if( fabs( x_sd[2] - 1.0 ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Gaussian RNG failed: standard deviation = " + to_string(x_sd[2]) ) );
    if( fabs( x_sd[3] - sqrt(1.0/3.0-2.0/M_PI/M_PI) ) > 1e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Cosine RNG failed: standard deviation = " + to_string(x_sd[3]) ) );
    if( fabs( x_sd[4] - sqrt(9.0)*0.5 ) > 2e-3 )
	throw( ErrorTest( ERROR_LOCATION, "Gamma(9,0.5) RNG failed: standard deviation = " + to_string(x_sd[4]) ) );

    // Output results
    ofstream fstr( "qrandom.dat" );
    fstr << "# Column 1: Linear [0,1], ave = " << x_ave[0] << ", sd = " << x_sd[0] << "\n";
    fstr << "# Column 2: Linear [0,1], ave = " << x_ave[1] << ", sd = " << x_sd[1] << "\n";
    fstr << "# Column 3: Gaussian, ave = " << x_ave[2] << ", sd = " << x_sd[2] << "\n";
    fstr << "# Column 4: Cosine, ave = " << x_ave[3] << ", sd = " << x_sd[3] << "\n";
    fstr << "# Column 5: Gamma(9,0.5), ave = " << x_ave[4] << ", sd = " << x_sd[4] << "\n";
    for( int a = 0; a < N; a++ ) { 
	fstr << setw(12) << x_vec[0][a] << " " 
	     << setw(12) << x_vec[1][a] << " "
	     << setw(12) << x_vec[2][a] << " "
	     << setw(12) << x_vec[3][a] << " "
	     << setw(12) << x_vec[4][a] << "\n";
    }
    fstr.close();
}


