/*! \file ccolmatrix.cpp 
 *  \brief Test compressed column matrix.
 *
 *  \test Test compressed column matrix.
 */


#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#define SPM_RANGE_CHECK 1
#include "ccolmatrix.hpp"
#include "mvector.hpp"


using namespace std;


#define ERROR()							    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << "\n"; \
	exit( 1 );						    \
    }


void test( int argc, char **argv )
{
    /* Constructors */
    CColMatrix A;
    if( A.columns() != 0 || A.rows() != 0 || A.nz_elements() != 0 )
	ERROR();
    CColMatrix B(5,5);
    if( B.columns() != 5 || B.rows() != 5 || B.nz_elements() != 0 )
	ERROR();
    B.set(0,0) = 2;
    B.set(0,1) = 3;
    B.set(0,4) = 4;
    CColMatrix C = B;
    if( C.get(0,0) != 2 || C.get(0,1) != 3 || C.get(0,2) != 0 || 
	C.get(0,3) != 0 || C.get(0,4) != 4 )
	ERROR();
    if( C.nz_elements() != 3 )
	ERROR();
    CColMatrix D(5,6);
    if( D.columns() != 6 || D.rows() != 5 || D.nz_elements() != 0 )
	ERROR();

    /* Resize */
    int rows, cols;
    D.size( rows, cols );
    if( cols != 6 || rows != 5 )
	ERROR();
    D.resize( 5, 7 );
    if( D.columns() != 7 || D.rows() != 5 )
	ERROR();
    A = D;
    if( A.columns() != 7 || A.rows() != 5 )
	ERROR();
    A.reserve( 10 );
    if( A.capacity() < 10 )
	ERROR();
    A.resize( 5, 5 );
    D.resize( 5, 5 );

    /* Set column */
    A.clear();
    int row[3] = {   0,   1,   3};
    double val[3] = {0.25, 0.5, 1.5};
    A.set_column( 0, 3, row, val );
    if( A.get(0,0) != 0.25 || A.get(1,0) != 0.5 || A.get(2,0) != 0.0 || 
	A.get(3,0) != 1.5 || A.get(4,0) != 0.0 )
	ERROR();
    if( A.nz_elements() != 3 )
	ERROR();
    B = A;
    A.clear();
    if( A.nz_elements() != 0 )
	ERROR();
    if( B.nz_elements() != 3 )
	ERROR();

    /* Range checking */
    int stat = 0;
    try {
	A.set(5,0) = 2;
    } catch( ErrorRange ) {
	stat = 1;
    }
    if( stat == 0 )
	ERROR();

    /* Vector-Matrix Multiplication */
    A.resize(6,5);
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
    if( Y.size() != 6 || Y[0] != 0.5 || Y[1] != -3.75 || Y[2] != 2.0 || 
	Y[3] != 4.5 || Y[4] != -4.0 || Y[5] != 0.5 ) {
	cout << "X = \n" << X << "\n";
	cout << "A = \n" << A << "\n";
	cout << "Y = \n" << Y << "\n";
	ERROR();
    }

    /* Merge, reserve, capacity */
    if( A.nz_elements() != 14 )
	ERROR();
    B.merge( A );
    if( B.nz_elements() != 14 || A.nz_elements() != 0 )
	ERROR();
    B.reserve( 20 );
    if( B.capacity() < 20 )
	ERROR();


    /* Clear element */
    int x = B.nz_elements();
    B.clear( 2, 2 );
    if( x - B.nz_elements() != 1 || B.get(2,2) != 0.0 )
	ERROR();

    /* Matrix order check */
    A.resize(3,3);
    A.set(0,0) = 1;
    A.set(1,0) = 2;
    A.set(0,1) = 3;
    A.set(1,1) = 4;
    A.set(2,1) = 5;
    A.set(0,2) = 6;
    A.set(2,2) = 7;
    if( !A.check_ascending() )
	ERROR();

    A.clear();
    A.set(0,0) = 1;
    A.set(1,0) = 2;
    A.set(2,1) = 5;
    A.set(1,1) = 4;
    A.set(0,1) = 3;
    A.set(0,2) = 6;
    A.set(2,2) = 7;
    if( A.check_ascending() )
	ERROR();
}


