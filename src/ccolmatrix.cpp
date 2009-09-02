#include "config.hpp"
#include "ccolmatrix.hpp"
#include "crowmatrix.hpp"
#include "coordmatrix.hpp"
#include "mvector.hpp"
#include "sort.hpp"
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iomanip>


inline void CColMatrix::allocate( void )
{
    if( !(_ptr = (uint32_t *)malloc( (_m+1)*sizeof(uint32_t) )) ) {
	_ptr = _row = NULL;
	_val = NULL;
	_n = _m = _nz = _asize = 0;	
	throw( ErrorNoMem( ERROR_LOCATION ) );
    }

    if( _asize == 0 ) {
	_row = NULL;
	_val = NULL;
    } else {
	if( !(_row = (uint32_t *)malloc( _asize*sizeof(uint32_t) )) ) {
	    free( _ptr );
	    _ptr = _row = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION) );
	}
	
	if( !(_val = (double *)malloc( _asize*sizeof(double) )) ) {
	    free( _ptr );
	    free( _row );
	    _ptr = _row = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION) );
	}
    }
}


inline void CColMatrix::reallocate( void )
{
    uint32_t *tmp;
    double *tmp2;

    if( !(tmp = (uint32_t *)realloc( _ptr, (_m+1)*sizeof(uint32_t) )) ) {
	free( _ptr );
	free( _row );
	free( _val );
	_ptr = _row = NULL;
	_val = NULL;
	_n = _m = _nz = _asize = 0;
	throw( ErrorNoMem( ERROR_LOCATION) );
    }
    _ptr = tmp;

    if( _asize == 0 ) {
	_row = NULL;
	_val = NULL;
    } else {
	if( !(tmp = (uint32_t *)realloc( _row, _asize*sizeof(uint32_t) )) ) {
	    free( _ptr );
	    free( _row );
	    free( _val );
	    _ptr = _row = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION) );
	}
	_row = tmp;
	if( !(tmp2 = (double *)realloc( _val, _asize*sizeof(double) )) ) {
	    free( _ptr );
	    free( _row );
	    free( _val );
	    _ptr = _row = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION) );
	}
	_val = tmp2;
    }
}


CColMatrix::CColMatrix()
    : _n(0), _m(0), _nz(0), _asize(0), _row(NULL), _val(NULL)
{
    allocate();
}


CColMatrix::CColMatrix( uint32_t n, uint32_t m )
{
    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    allocate();
    memset( _ptr, 0, (_m+1)*sizeof(uint32_t) );
}


CColMatrix::CColMatrix( uint32_t n, uint32_t m, uint32_t nz, 
			uint32_t *ptr, uint32_t *row, double *val )
{
    _n     = n;
    _m     = m;
    _nz    = nz;
    _asize = nz;
    _ptr   = ptr;
    _row   = row;
    _val   = val;
}


CColMatrix::CColMatrix( const CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();

    memcpy( _ptr, mat._ptr, (_m+1)*sizeof(uint32_t) );
    memcpy( _row, mat._row, _nz*sizeof(uint32_t) );
    memcpy( _val, mat._val, _nz*sizeof(double) );
}


CColMatrix &CColMatrix::operator=( const CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    memcpy( _ptr, mat._ptr, (_m+1)*sizeof(uint32_t) );
    memcpy( _row, mat._row, _nz*sizeof(uint32_t) );
    memcpy( _val, mat._val, _nz*sizeof(double) );

    return( *this );
}



CColMatrix::CColMatrix( const class CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();

    uint32_t *c = new uint32_t[_m];
    uint32_t i, j, start, cl;

    /* Count number of entries in each column to c. */
    memset( c, 0, _m*sizeof(uint32_t) );
    for( i = 0; i < _nz; i++ )
	c[mat._col[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _m; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a row pointer. */
    for( i = 0, j = 0; i < _n; i++ ) {
	for( ; j < mat._ptr[i+1]; j++ ) {
	    cl = mat._col[j];
	    start = _ptr[cl];
	    c[cl]--;
	    _row[start+c[cl]] = i;
	    _val[start+c[cl]] = mat._val[j];
	}
    }

    delete [] c;
}


CColMatrix &CColMatrix::operator=( const class CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    uint32_t *c = new uint32_t[_m];
    uint32_t i, j, start, cl;

    /* Count number of entries in each column to c. */
    memset( c, 0, _m*sizeof(uint32_t) );
    for( i = 0; i < _nz; i++ )
	c[mat._col[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _m; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a row pointer. */
    for( i = 0, j = 0; i < _n; i++ ) {
	for( ; j < mat._ptr[i+1]; j++ ) {
	    cl = mat._col[j];
	    start = _ptr[cl];
	    c[cl]--;
	    _row[start+c[cl]] = i;
	    _val[start+c[cl]] = mat._val[j];
	}
    }

    delete [] c;

    return( *this );
}


CColMatrix::CColMatrix( const class CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();

    uint32_t *c = new uint32_t[_m];
    uint32_t i, start, cl;

    /* Count number of entries in each column to c. */
    memset( c, 0, _m*sizeof(uint32_t) );
    for( i = 0; i < _nz; i++ )
	c[mat._col[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _m; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a row pointer. */
    for( i = 0; i < _nz; i++ ) {
	cl = mat._col[i];
	start = _ptr[cl];
	c[cl]--;
	_row[start+c[cl]] = mat._row[i];
	_val[start+c[cl]] = mat._val[i];
    }

    delete [] c;
}


CColMatrix &CColMatrix::operator=( const class CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    uint32_t *c = new uint32_t[_m];
    uint32_t i, start, cl;

    /* Count number of entries in each column to c. */
    memset( c, 0, _m*sizeof(uint32_t) );
    for( i = 0; i < _nz; i++ )
	c[mat._col[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _m; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a row pointer. */
    for( i = 0; i < _nz; i++ ) {
	cl = mat._col[i];
	start = _ptr[cl];
	c[cl]--;
	_row[start+c[cl]] = mat._row[i];
	_val[start+c[cl]] = mat._val[i];
    }

    delete [] c;

    return( *this );
}


CColMatrix::~CColMatrix()
{
    free( _ptr );
    free( _row );
    free( _val );
}


void CColMatrix::resize( uint32_t n, uint32_t m )
{
    free( _row );
    free( _val );
    _row = NULL;
    _val = NULL;

    if( _m != m ) {
	uint32_t *tmp;
	if( !(tmp = (uint32_t *)realloc( _ptr, (m+1)*sizeof(uint32_t) )) ) {
	    free( _ptr );
	    _ptr = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_ptr = tmp;
    }

    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    memset( _ptr, 0, (_m+1)*sizeof(uint32_t) );
}


void CColMatrix::clear( void )
{
    _nz    = 0;
    _asize = 0;
    memset( _ptr, 0, (_m+1)*sizeof(uint32_t) );

    free( _row );
    free( _val );
    _row = NULL;
    _val = NULL;
}


inline void CColMatrix::clear_no_check( uint32_t i, uint32_t j )
{
    uint32_t a;
    /* Search for the element */
    for( a = _ptr[i]; a < _ptr[i+1]; a++ )
	if( _row[a] == j )
	    break;
    /* Do nothing if the element does not exist */
    if( a == _ptr[i+1] )
	return;

    /* Move data */
    int movesize = _nz-a-1;
    memmove( &_val[a], &_val[a+1], movesize*sizeof(double) );
    memmove( &_row[a], &_row[a+1], movesize*sizeof(uint32_t) );

    /* Update pointers */
    _nz--;
    for( uint32_t a = i+1; a < _m+1; a++ )
	_ptr[a]--;

}


void CColMatrix::clear_check( uint32_t i, uint32_t j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    clear_no_check( i, j );
}


void CColMatrix::reserve( uint32_t size )
{
    if( size > _asize ) {
	_asize = size;
	reallocate();
    }
}


void CColMatrix::set_nz( uint32_t nz )
{
    if( nz > _asize ) {
	_asize = nz;
	reallocate();
    }
    _nz = nz;
}


void CColMatrix::merge( CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._asize;
    _ptr   = mat._ptr;
    _row   = mat._row;
    _val   = mat._val;

    mat._n     = 0;
    mat._m     = 0;
    mat._nz    = 0;
    mat._asize = 0;
    mat.allocate();
}


void CColMatrix::order_ascending( void )
{
    /* Sort each column. */
    for( uint32_t i = 0; i < _m; i++ )
	sort_iv( _row, _val, _ptr[i], _ptr[i+1] );
}


inline double CColMatrix::get_no_check( uint32_t i, uint32_t j ) const
{
    for( uint32_t a = _ptr[j]; a < _ptr[j+1]; a++ )
	if( _row[a] == i )
	    return( _val[a] );
    return( 0.0 );
}


inline double &CColMatrix::set_no_check( uint32_t i, uint32_t j )
{
    /* Use existing element if it exists */
    for( uint32_t a = _ptr[j]; a < _ptr[j+1]; a++ )
	if( _row[a] == i )
	    return( _val[a] );

    /* Reserve new space if necessary */
    if( _nz+1 > _asize )
	reserve( _asize+_n );

    /* Move existing data */
    int movesize = _nz-_ptr[j+1];
    memmove( &_val[_ptr[j+1]+1], &_val[_ptr[j+1]], movesize*sizeof(double) );
    memmove( &_row[_ptr[j+1]+1], &_row[_ptr[j+1]], movesize*sizeof(uint32_t) );

    /* Set new data */
    _val[_ptr[j+1]] = 0.0;
    _row[_ptr[j+1]] = i;

    /* Update pointers */
    _nz++;
    for( uint32_t a = j+1; a < _m+1; a++ )
	_ptr[a]++;

    return( _val[_ptr[j+1]-1] );
}


double CColMatrix::get_check( uint32_t i, uint32_t j ) const
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( get_no_check( i, j ) );
}


double &CColMatrix::set_check( uint32_t i, uint32_t j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( set_no_check( i, j ) );
}


void CColMatrix::set_column( uint32_t j, uint32_t N, const uint32_t *row, const double *val )
{
    if( j >= _m )
	throw( ErrorRange( ERROR_LOCATION, j, _m ) );
    else if( N > _n )
	throw( Error( ERROR_LOCATION, "too many elements" ) );

    /* Reserve new space if necessary */
    uint32_t oldsize = _ptr[j+1] - _ptr[j];
    if( _nz+N-oldsize > _asize )
	reserve( _asize+_n );

    /* Move existing data */
    int offset = N-oldsize;
    int movesize = _nz-_ptr[j+1];
    memmove( &_val[_ptr[j+1]+offset], &_val[_ptr[j+1]], movesize*sizeof(double) );
    memmove( &_row[_ptr[j+1]+offset], &_row[_ptr[j+1]], movesize*sizeof(uint32_t) );

    /* Set new data */
    int err = 0;
    for( uint32_t a = 0; a < N; a++ ) {
	for( uint32_t b = a+1; b < N; b++ ) {
	    if( row[a] == row[b] )
		err = 2;
	}
	_val[_ptr[j]+a] = val[a];
	if( (_row[_ptr[j]+a] = row[a]) >= _n )
	    err = 1;
    }

    /* Update pointers */
    _nz += offset;
    for( uint32_t a = j+1; a < _m+1; a++ )
	_ptr[a] += offset;

    if( err == 1 )
	throw( Error( ERROR_LOCATION, "row index out of range" ) );
    else if( err == 2 )
	throw( Error( ERROR_LOCATION, "repeated row index" ) );
}


void CColMatrix::construct_add( uint32_t i, uint32_t j, double val )
{
    /* Reserve new space if necessary */
    if( _nz+1 > _asize )
	reserve( _asize+_n );

    /* Set new data */
    _val[_ptr[j+1]] = val;
    _row[_ptr[j+1]] = i;

    /* Update the number of nonzeroes and only the next pointer */
    _nz++;
    _ptr[j+1]++;
    if( j+1 < _m )
	_ptr[j+2] = _ptr[j+1];

}


void CColMatrix::debug_print( void ) const
{
    std::cout << "n     = " << _n << "\n";
    std::cout << "m     = " << _m << "\n";
    std::cout << "nz    = " << _nz << "\n";
    std::cout << "asize = " << _asize << "\n";

    std::cout << "ptr[] = {";
    for( uint32_t i = 0; i < _m; i++ )
	std::cout << _ptr[i] << ", ";
    std::cout << _ptr[_m] << "}\n";

    std::cout << "row[] = {";
    if( _nz <= 0 ) {
	std::cout << "}\n";
    } else {
	for( uint32_t i = 0; i < _nz-1; i++ )
	    std::cout << _row[i] << ", ";
	std::cout << _row[_nz-1] << "}\n";
    }

    std::cout << "val[] = {";
    if( _nz <= 0 ) {
	std::cout << "}\n";
    } else {
	for( uint32_t i = 0; i < _nz-1; i++ )
	    std::cout << _val[i] << ", ";
	std::cout << _val[_nz-1] << "}\n";
    }
}



/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */


void CColMatrix::multiply_by_vector( Vector &x, const Vector &b ) const
{
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x.resize( _n );
    x.clear();

    double mult;
    for( uint32_t i = 0; i < _n; i++ ) {
	mult = b._val[i];
	for( uint32_t a = _ptr[i]; a < _ptr[i+1]; a++ )
	    x._val[_row[a]] += _val[a] * mult;
    }
}


void CColMatrix::lower_unit_solve( Vector &x, const Vector &b ) const
{
    // Make checks
    if( _n != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix not squrare" ) );
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x = b;

    double mult;
    for( uint32_t i = 0; i < _m; i++ ) {
	mult = x[i];
	for( uint32_t j = _ptr[i]; j < _ptr[i+1]; j++ )
	    x[_row[j]] -= _val[j] * mult;
    }
}


void CColMatrix::upper_diag_solve( Vector &x, const Vector &b ) const
{
    // Make checks
    if( _n != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x = b;

    double mult;
    int32_t i, j;
    for( i = _m-1; i >= 0; i-- ) {
	j = _ptr[i+1]-1;
	x[i] /= _val[j];
	mult = x[i];
	for( j--; j >= (int32_t)_ptr[i]; j-- )
	    x[_row[j]] -= _val[j] * mult;
	
    }
}

















