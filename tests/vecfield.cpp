/*! \file vecfield.cpp 
 *  \test Test VectorField class
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "vectorfield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


using namespace std;


void test( void )
{
    bool fout[3] = {false, false, true};
    Int3D size(2,2,1);
    Vec3D origo(0,0,0);
    VectorField f( MODE_2D, fout, size, origo, 0.01 );
    f.set(0,0,Vec3D(0,0,0));
    f.set(1,0,Vec3D(0,0,1.5));
    f.set(0,1,Vec3D(0,0,1));
    f.set(1,1,Vec3D(0,0,-1.5));

    ofstream ostr( "vecfield.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " "
	 << setw(14) << "field (au)" << "\n";
    for( double y = 0.0; y <= 0.01; y += 0.002 ) {
	for( double x = 0.0; x <= 0.01; x += 0.002 ) {
	    Vec3D v(x,y,0);
	    ostr << setw(12) << x << " "
		 << setw(12) << y << " "
		 << setw(12) << f(v) << "\n";
	}
	ostr << "\n";
    }
}


int main( void )
{
    try {
	ibsimu.set_verbose_output( 0 );
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}



