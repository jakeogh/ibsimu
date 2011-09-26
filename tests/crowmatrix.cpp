/*! \file crowmatrix.cpp 
 *  \brief Test compressed row matrix.
 *
 *  \test Test compressed row matrix.
 */


#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#define SPM_RANGE_CHECK 1
#include "crowmatrix.hpp"
#include "mvector.hpp"
#include "ibsimutest.hpp"


using namespace std;


void test( int argc, char **argv )
{
    /* Constructors */
    CRowMatrix A;
    if( A.columns() != 0 || A.rows() != 0 || A.nz_elements() != 0 )
	throw( ErrorTest( ERROR_LOCATION ) );
    CRowMatrix B(5,5);
    if( B.columns() != 5 || B.rows() != 5 || B.nz_elements() != 0 )
	throw( ErrorTest( ERROR_LOCATION ) );
    B.set(0,0) = 2;
    B.set(0,1) = 3;
    B.set(0,4) = 4;
    CRowMatrix C = B;
    if( C.get(0,0) != 2 || C.get(0,1) != 3 || C.get(0,2) != 0 || 
	C.get(0,3) != 0 || C.get(0,4) != 4 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( C.nz_elements() != 3 )
	throw( ErrorTest( ERROR_LOCATION ) );
    CRowMatrix D(5,6);
    if( D.columns() != 6 || D.rows() != 5 || D.nz_elements() != 0 )
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Resize */
    int rows, cols;
    D.size( rows, cols );
    if( cols != 6 || rows != 5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    D.resize( 5, 7 );
    if( D.columns() != 7 || D.rows() != 5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    A = D;
    if( A.columns() != 7 || A.rows() != 5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    A.reserve( 10 );
    if( A.capacity() < 10 )
	throw( ErrorTest( ERROR_LOCATION ) );
    A.resize( 5, 5 );
    D.resize( 5, 5 );


    /* Set row */
    A.clear();
    int col[3] = {   0,   1,   3};
    double val[3] = {0.25, 0.5, 1.5};
    A.set_row( 0, 3, col, val );
    if( A.get(0,0) != 0.25 || A.get(0,1) != 0.5 || A.get(0,2) != 0.0 || 
	A.get(0,3) != 1.5 || A.get(0,4) != 0.0 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( A.nz_elements() != 3 )
	throw( ErrorTest( ERROR_LOCATION ) );
    B = A;
    A.clear();
    if( A.nz_elements() != 0 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( B.nz_elements() != 3 )
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Range checking */
    int stat = 0;
    try {
	A.set(5,0) = 2;
    } catch( ErrorRange ) {
	stat = 1;
    }
    if( stat == 0 )
	throw( ErrorTest( ERROR_LOCATION ) );

    /* Vector-Matrix Multiplication */
    A.resize(6,5);
    A.reserve(15);
    A.set(0,0) = -2;
    A.set(1,1) = -2;
    A.set(2,2) = -2;
    A.set(3,3) = -2;
    A.set(4,4) = -2;
    A.set(0,1) = 1;
    A.set(1,2) = 1;
    A.set(2,3) = 1;
    A.set(3,4) = 1;
    A.set(1,0) = 1;
    A.set(2,1) = 1;
    A.set(3,2) = 1;
    A.set(4,3) = 1;
    A.set(5,4) = 1;
    double X_data[5] = {0.25,1,-2,-3,0.5};
    Vector X(5,X_data);
    Vector Y = A*X;
    Y = A*X;
    if( Y.size() != 6 || Y[0] != 0.5 || Y[1] != -3.75 || Y[2] != 2.0 || 
	Y[3] != 4.5 || Y[4] != -4.0 || Y[5] != 0.5 )
	throw( ErrorTest( ERROR_LOCATION ) );

    /* Merge, reserve, capacity */
    if( A.nz_elements() != 14 )
	throw( ErrorTest( ERROR_LOCATION ) );
    B.merge( A );
    if( B.nz_elements() != 14 || A.nz_elements() != 0 )
	throw( ErrorTest( ERROR_LOCATION ) );
    B.reserve( 20 );
    if( B.capacity() < 20 )
	throw( ErrorTest( ERROR_LOCATION ) );

    /* Clear element */
    int x = B.nz_elements();
    B.clear( 2, 2 );
    if( x - B.nz_elements() != 1 || B.get(2,2) != 0.0 )
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Matrix order check */
    A.resize(3,3);
    A.set(0,0) = 1;
    A.set(0,1) = 2;
    A.set(1,0) = 3;
    A.set(1,1) = 4;
    A.set(1,2) = 5;
    A.set(2,0) = 6;
    A.set(2,2) = 7;
    if( !A.check_ascending() )
	throw( ErrorTest( ERROR_LOCATION ) );

    A.clear();
    A.set(0,0) = 1;
    A.set(0,1) = 2;
    A.set(1,2) = 5;
    A.set(1,1) = 4;
    A.set(1,0) = 3;
    A.set(2,0) = 6;
    A.set(2,2) = 7;
    if( A.check_ascending() )
	throw( ErrorTest( ERROR_LOCATION ) );
}

