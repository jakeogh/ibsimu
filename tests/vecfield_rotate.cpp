/*! \file vecfield_rotate.cpp 
 *  \test Test VectorField rotate functions
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "vectorfield.hpp"
#include "error.hpp"


using namespace std;


#define ERROR()							    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << "\n"; \
	exit( 1 );						    \
    }


void test( void )
{
    bool fout[3] = {true, true, false};
    Int3D size(2,3,1);
    Vec3D origo(2.2,1.1,0);
    VectorField f( MODE_2D, fout, size, origo, 1.0 );
    f.set( 0, 0, Vec3D(  0.1,  0.4, 0 ) );
    f.set( 1, 0, Vec3D(  0.2,  1.5, 0 ) );
    f.set( 0, 1, Vec3D(  0.3,  1.0, 0 ) );
    f.set( 1, 1, Vec3D( -0.4, -1.5, 0 ) );
    f.set( 0, 2, Vec3D(  0.5,  4.0, 0 ) );
    f.set( 1, 2, Vec3D(  0.6,  3.0, 0 ) );

    /*
    for( int32_t a = 0; a < 3; a++ ) {
	cout << "Old " << a << ":\n";
	for( int32_t j = f.size(1)-1; j >= 0; j-- ) {
	    for( int32_t i = 0; i < f.size(0); i++ )
		cout << setw(12) << f(i,j)(a) << " ";
	    cout << "\n";
	}
	cout << "\n";
    }
    f.debug_print();

    cout << "\n";
    cout << "Rotate-----------------------\n";
    cout << "\n";
    f.rotate_z( 90 );

    for( int32_t a = 0; a < 3; a++ ) {
	cout << "New " << a << ":\n";
	for( int32_t j = f.size(1)-1; j >= 0; j-- ) {
	    for( int32_t i = 0; i < f.size(0); i++ )
		cout << setw(12) << f(i,j)(a) << " ";
	    cout << "\n";
	}
	cout << "\n";
    }
    f.debug_print();
    */

    if( norm2( Vec3D(0.1,0.4,0)-f(Vec3D(2.2,1.1,0)) ) > 0.01 )
	ERROR();

    f.rotate_z( 90 );

    if( norm2( Vec3D(-0.4,0.1,0)-f(Vec3D(-1.1,2.2,0)) ) > 0.01 )
	ERROR();

    f.rotate_z( -90 );

    if( norm2( Vec3D(0.1,0.4,0)-f(Vec3D(2.2,1.1,0)) ) > 0.01 )
	ERROR();

    f.rotate_z( 180 );

    if( norm2( Vec3D(-0.1,-0.4,0)-f(Vec3D(-2.2,-1.1,0)) ) > 0.01 )
	ERROR();

    f.scale( 0.5 );

    if( norm2( Vec3D(-0.1,-0.4,0)-f(Vec3D(-1.1,-0.55,0)) ) > 0.01 )
	ERROR();

    f.scale( -2.0 );

    if( norm2( Vec3D(0.1,0.4,0)-f(Vec3D(2.2,1.1,0)) ) > 0.01 )
	ERROR();
}


int main( void )
{
    try {
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

