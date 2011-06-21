/*! \file matvec3dtest.cpp 
 *  \brief Basic tests for 3D vectors and matrices
 *
 *  \test Basic tests for 3D vectors and matrices
 *
 */


#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "vec3d.hpp"
#include "mat3d.hpp"
#include "ibsimutest.hpp"


using namespace std;


void test_equal( void )
{
    Vec3D x(1,0,3);
    Vec3D y(1,0,3*(1.0+1.0e-7));      // Same as x with small relative error and small absolute error
    Vec3D z(1,1.0e-7,3);              // Same as x with infinite relative error and small absolute error

    Vec3D x2(1,0,30);
    Vec3D y2(1,0,30*(1.0+1.0e-7));    // Same as x2 with small relative error and large absolute error
    Vec3D z2(3,0,3);

    // Test appox() and == operators
    if( !(x == x) )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( x == z )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( !(x.approx(x)) )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( !(x.approx(y)) )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( !(x.approx(z)) )
	throw( ErrorTest( ERROR_LOCATION ) );

    if( !(x2.approx(x2)) )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( !(x2.approx(y2)) )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( x.approx(x2) )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( x.approx(z2) )
	throw( ErrorTest( ERROR_LOCATION ) );

    // Test != operator
    if( x != x )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( !(x != z) )
	throw( ErrorTest( ERROR_LOCATION ) );
}


void test( int argc, char **argv )
{
    test_equal();

    /* Vector ssqr, norm2, normalize
     */
    Vec3D a(2,-1,2);
    if( fabs(a.ssqr()-9) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( fabs(a.norm2()-3) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( fabs(norm2(a)-3) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );
    a.normalize();
    if( norm2(a-Vec3D(2.0/3.0,-1.0/3.0,2.0/3.0)) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );

    /* Matrix determinant, matrix inversion,
     * vector-matrix multiplication
     */
    Vec3D x( 1, 2, 3 );
    Mat3D m(  1, -2,  3, 
	      1,  2,  2,
	     -1,  4, -2 );
    double mdet = m.determinant();
#ifdef PRINT
    cout << "x = \n" << x << "\n";
    cout << "m = \n" << m << "\n";
    cout << "mdet = " << mdet << "\n";
#endif
    if( fabs(mdet-6.0) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );
    Vec3D b = m*x;
#ifdef PRINT
    cout << "b = m*x = \n" << b << "\n";
#endif
    if( norm2(b-Vec3D(6,11,1)) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );
    Mat3D minv = m.inverse();
    Vec3D xnew = minv*b;
#ifdef PRINT
    cout << "minv = \n" << minv << "\n";
    cout << "xnew = minv*b = \n" << xnew << "\n";
#endif
    if( norm2(xnew-x) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );

    /* Dot product
     */
    double dot = x*b;
#ifdef PRINT
    cout << "x*b = " << dot << "\n";
#endif
    if( fabs(dot-31.0) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );

    /* Cross product
     */
    Vec3D y = cross(x,b);
#ifdef PRINT
    cout << "y = cross(x,b) = \n" << y << "\n";
#endif
    if( norm2(y-Vec3D(-31,17,-1)) > 1e-6 )
	throw( ErrorTest( ERROR_LOCATION ) );
}
