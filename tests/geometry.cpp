/* geometry.cpp -- Test geometry class
 *
 *
 */



#include <fstream>
#include <iomanip>
#include "geometry.hpp"
#include "func_solid.hpp"
#include "error.hpp"


using namespace std;


bool s1( double x, double y, double z )
{
    /*
    cout << "Testing at ("
	 << x << ", "
	 << y << ", "
	 << z << ")\n";
    */

    return( x*x + y*y < 0.1*0.1 );
}


void test1( void )
{
    Geometry g( MODE_2D, Int3D(20,20,1), Vec3D(0,0,0), 0.01 );
    Solid *s = new FuncSolid( s1 );
    g.set_solid( 7, s );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET, 20.0) );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET, 20.0) );
    g.build_mesh();
    //g.debug_print();

    ofstream fstr( "geometry.dat" );
    g.save( fstr );
    fstr.close();
}


void test2( void )
{
    ifstream fstr( "geometry.dat" );
    Geometry g( fstr );
    //g.debug_print();
    int stat = 0;
    try {
	g.inside( Vec3D(0.05,0.05,0.05) );
    } catch( Error ) {
	stat = 1;
    }
    if( stat == 0 )
	exit( 1 );

    /*
    for( int j = 0; j < g.size(1); j++ ) {
	for( int i = 0; i < g.size(0); i++ )
	    cout << setw(3) << (int)g.mesh( i, j );
	cout << "\n";
    }
    */
}


int main( void )
{
    try {
	test1();
	test2();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

