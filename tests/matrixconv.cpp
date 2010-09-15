/* matrixconv.cpp -- Test matrix-to-matrix conversion
 *
 *
 */


#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#define SPM_RANGE_CHECK 1
#include "matrix.hpp"
#include "crowmatrix.hpp"
#include "ccolmatrix.hpp"
#include "coordmatrix.hpp"


using namespace std;


#define ERROR()							    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << "\n"; \
	exit( 1 );						    \
    }



#define INTERNAL_ERROR()					    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << ": Test not valid\n"; \
	exit( 1 );						    \
    }



bool compare( Matrix &A, Matrix &B )
{
    if( A.columns() != B.columns() || A.rows() != B.rows() ) {
	cout << "Different matrix sizes.\n";
	return( false );
    }
    for( int i = 0; i < A.rows(); i++ ) {
	for( int j = 0; j < A.columns(); j++ ) {
	    if( A.get(i,j) != B.get(i,j) ) {
		cout << "Different matrix content.\n";
		return( false );
	    }
	}
    }
    return( true );
}


void test( void )
{
    CRowMatrix Row(6,5);
    Row.set(4,1) = -91;
    Row.set(4,0) = -49;
    Row.set(3,0) = -45;
    Row.set(1,0) = 11;
    Row.set(0,0) = -25;
    Row.set(5,1) = -33;
    Row.set(5,4) = -18;
    Row.set(0,1) = -1;
    Row.set(0,4) = -9;
    Row.set(2,2) = 31;
    Row.set(1,2) = 32;
    Row.set(4,4) = -22;
    Row.set(1,1) = -11;
    Row.set(2,3) = -9;
    Row.set(2,1) = -4;
    Row.set(3,3) = 32;
    Row.set(3,2) = 25;
    Row.set(4,3) = 5;

    /* Constructor converters */

    CColMatrix Col = Row;
    if( !compare( Col, Row ) )
	ERROR();

    CoordMatrix Coord = Row;
    if( !compare( Coord, Row ) )
	ERROR();

    CoordMatrix Coord2 = Col;
    if( !compare( Coord2, Row ) )
	ERROR();

    CColMatrix Col2 = Coord;
    if( !compare( Col2, Row ) )
	ERROR();

    CRowMatrix Row2 = Col;
    if( !compare( Row2, Row ) )
	ERROR();

    CRowMatrix Row3 = Coord;
    if( !compare( Row3, Row ) )
	ERROR();

    /* Copy converters */

    Col.clear();
    Coord.clear();
    Coord2.clear();
    Col2.clear();
    Row2.clear();
    Row3.clear();

    Col = Row;
    if( !compare( Col, Row ) )
	ERROR();

    Coord = Row;
    if( !compare( Coord, Row ) )
	ERROR();

    Coord2 = Col;
    if( !compare( Coord2, Row ) )
	ERROR();

    Col2 = Coord;
    if( !compare( Col2, Row ) )
	ERROR();

    Row2 = Col;
    if( !compare( Row2, Row ) )
	ERROR();

    Row3 = Coord;
    if( !compare( Row3, Row ) )
	ERROR();

    /* Conversion ordering check */
    if( Row.check_ascending() )
	INTERNAL_ERROR();

    CColMatrix Col5 = Row;
    if( !Col5.check_ascending() )
	ERROR();
    
    CColMatrix Col6(6,5);
    Col6.set(4,1) = -91;
    Col6.set(4,0) = -49;
    Col6.set(3,0) = -45;
    Col6.set(1,0) = 11;
    Col6.set(0,0) = -25;
    Col6.set(5,1) = -33;
    Col6.set(5,4) = -18;
    Col6.set(0,1) = -1;
    Col6.set(0,4) = -9;
    Col6.set(2,2) = 31;
    Col6.set(1,2) = 32;
    Col6.set(4,4) = -22;
    Col6.set(1,1) = -11;
    Col6.set(2,3) = -9;
    Col6.set(2,1) = -4;
    Col6.set(3,3) = 32;
    Col6.set(3,2) = 25;
    Col6.set(4,3) = 5;

    if( Col6.check_ascending() )
	INTERNAL_ERROR();

    CRowMatrix Row6 = Col6;
    if( !Row6.check_ascending() )
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



