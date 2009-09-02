#include "config.hpp"
#include "crowmatrix.hpp"
#include "ccolmatrix.hpp"
#include "coordmatrix.hpp"
#include "mvector.hpp"
#include "sort.hpp"
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iomanip>


inline void CRowMatrix::allocate( void )
{
    if( !(_ptr = (uint32_t *)malloc( (_n+1)*sizeof(uint32_t) )) ) {
	_ptr = _col = NULL;
	_val = NULL;
	_n = _m = _nz = _asize = 0;	
	throw( ErrorNoMem( ERROR_LOCATION ) );
    }

    if( _asize == 0 ) {
	_col = NULL;
	_val = NULL;
    } else {
	if( !(_col = (uint32_t *)malloc( _asize*sizeof(uint32_t) )) ) {
	    free( _ptr );
	    _ptr = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	
	if( !(_val = (double *)malloc( _asize*sizeof(double) )) ) {
	    free( _ptr );
	    free( _col );
	    _ptr = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
    }
}


inline void CRowMatrix::reallocate( void )
{
    uint32_t *tmp;
    double *tmp2;

    if( !(tmp = (uint32_t *)realloc( _ptr, (_n+1)*sizeof(uint32_t) )) ) {
	free( _ptr );
	free( _col );
	free( _val );
	_ptr = _col = NULL;
	_val = NULL;
	_n = _m = _nz = _asize = 0;
	throw( ErrorNoMem( ERROR_LOCATION ) );
    }
    _ptr = tmp;

    if( _asize == 0 ) {
	_col = NULL;
	_val = NULL;
    } else {
	if( !(tmp = (uint32_t *)realloc( _col, _asize*sizeof(uint32_t) )) ) {
	    free( _ptr );
	    free( _col );
	    free( _val );
	    _ptr = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_col = tmp;
	if( !(tmp2 = (double *)realloc( _val, _asize*sizeof(double) )) ) {
	    free( _ptr );
	    free( _col );
	    free( _val );
	    _ptr = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_val = tmp2;
    }
}


CRowMatrix::CRowMatrix()
    : _n(0), _m(0), _nz(0), _asize(0), _col(NULL), _val(NULL)
{
    allocate();
}


CRowMatrix::CRowMatrix( uint32_t n, uint32_t m )
{
    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    allocate();
    memset( _ptr, 0, (_n+1)*sizeof(uint32_t) );
}


CRowMatrix::CRowMatrix( uint32_t n, uint32_t m, uint32_t nz, 
			uint32_t *ptr, uint32_t *col, double *val )
{
    _n     = n;
    _m     = m;
    _nz    = nz;
    _asize = nz;
    _ptr   = ptr;
    _col   = col;
    _val   = val;
}


CRowMatrix::CRowMatrix( const CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();

    memcpy( _ptr, mat._ptr, (_n+1)*sizeof(uint32_t) );
    memcpy( _col, mat._col, _nz*sizeof(uint32_t) );
    memcpy( _val, mat._val, _nz*sizeof(double) );
}


CRowMatrix &CRowMatrix::operator=( const CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    memcpy( _ptr, mat._ptr, (_n+1)*sizeof(uint32_t) );
    memcpy( _col, mat._col, _nz*sizeof(uint32_t) );
    memcpy( _val, mat._val, _nz*sizeof(double) );

    return( *this );
}


CRowMatrix::CRowMatrix( const class CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();
    
    uint32_t *c = new uint32_t[_n];
    uint32_t i, j, start, r;

    /* Count number of entries in each row to c. */
    memset( c, 0, _n*sizeof(uint32_t) );
    for( i = 0; i < _nz; i++ )
	c[mat._row[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _n; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a column pointer. */
    for( i = 0, j = 0; i < _m; i++ ) {
	for( ; j < mat._ptr[i+1]; j++ ) {
	    r = mat._row[j];
	    start = _ptr[r];
	    c[r]--;
	    _col[start+c[r]] = i;
	    _val[start+c[r]] = mat._val[j];
	}
    }

    delete [] c;
}


CRowMatrix &CRowMatrix::operator=( const class CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    uint32_t *c = new uint32_t[_n];
    uint32_t i, j, start, r;

    /* Count number of entries in each row to c. */
    memset( c, 0, _n*sizeof(uint32_t) );
    for( i = 0; i < _nz; i++ )
	c[mat._row[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _n; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a column pointer. */
    for( i = 0, j = 0; i < _m; i++ ) {
	for( ; j < mat._ptr[i+1]; j++ ) {
	    r = mat._row[j];
	    start = _ptr[r];
	    c[r]--;
	    _col[start+c[r]] = i;
	    _val[start+c[r]] = mat._val[j];
	}
    }

    delete [] c;

    return( *this );
}


CRowMatrix::CRowMatrix( const class CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();

    uint32_t *c = new uint32_t[_n];
    uint32_t i, start, r;

    /* Count number of entries in each row to c. */
    memset( c, 0, _n*sizeof(uint32_t) );
    for( i = 0; i < _nz; i++ )
	c[mat._row[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _n; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a column pointer. */
    for( i = 0; i < _nz; i++ ) {
	r = mat._row[i];
	start = _ptr[r];
	c[r]--;
	_col[start+c[r]] = mat._col[i];
	_val[start+c[r]] = mat._val[i];
    }

    delete [] c;
}


CRowMatrix &CRowMatrix::operator=( const class CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    uint32_t *c = new uint32_t[_n];
    uint32_t i, start, r;

    /* Count number of entries in each row to c. */
    memset( c, 0, _n*sizeof(uint32_t) );
    for( i = 0; i < _nz; i++ )
	c[mat._row[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _n; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a column pointer. */
    for( i = 0; i < _nz; i++ ) {
	r = mat._row[i];
	start = _ptr[r];
	c[r]--;
	_col[start+c[r]] = mat._col[i];
	_val[start+c[r]] = mat._val[i];
    }

    delete [] c;

    return( *this );
}



CRowMatrix::~CRowMatrix()
{
    free( _ptr );
    free( _col );
    free( _val );
}


void CRowMatrix::resize( uint32_t n, uint32_t m )
{
    free( _col );
    free( _val );
    _col = NULL;
    _val = NULL;

    if( _n != n ) {
	uint32_t *tmp;
	if( !(tmp = (uint32_t *)realloc( _ptr, (n+1)*sizeof(uint32_t) )) ) {
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
    memset( _ptr, 0, (_n+1)*sizeof(uint32_t) );
}


void CRowMatrix::clear( void )
{
    _nz    = 0;
    _asize = 0;
    memset( _ptr, 0, (_n+1)*sizeof(uint32_t) );

    free( _col );
    free( _val );
    _col = NULL;
    _val = NULL;
}


inline void CRowMatrix::clear_no_check( uint32_t i, uint32_t j )
{
    uint32_t a;

    /* Search for the element */
    for( a = _ptr[i]; a < _ptr[i+1]; a++ )
	if( _col[a] == j )
	    break;
    /* Do nothing if the element does not exist */
    if( a == _ptr[i+1] )
	return;

    /* Move data */
    int movesize = _nz-a-1;
    memmove( &_val[a], &_val[a+1], movesize*sizeof(double) );
    memmove( &_col[a], &_col[a+1], movesize*sizeof(uint32_t) );

    /* Update pointers */
    _nz--;
    for( a = i+1; a < _n+1; a++ )
	_ptr[a]--;
}


void CRowMatrix::clear_check( uint32_t i, uint32_t j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    clear_no_check( i, j );
}


void CRowMatrix::reserve( uint32_t size )
{
    if( size > _asize ) {
	_asize = size;
	reallocate();
    }
}


void CRowMatrix::set_nz( uint32_t nz )
{
    if( nz > _asize ) {
	_asize = nz;
	reallocate();
    }
    _nz = nz;
}


void CRowMatrix::merge( CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._asize;
    _ptr   = mat._ptr;
    _col   = mat._col;
    _val   = mat._val;

    mat._n     = 0;
    mat._m     = 0;
    mat._nz    = 0;
    mat._asize = 0;
    mat.allocate();
}


void CRowMatrix::order_ascending( void )
{
    /* Sort each row. */
    for( uint32_t i = 0; i < _n; i++ )
	sort_iv( _col, _val, _ptr[i], _ptr[i+1] );
}


inline double CRowMatrix::get_no_check( uint32_t i, uint32_t j ) const
{
    for( uint32_t a = _ptr[i]; a < _ptr[i+1]; a++ )
	if( _col[a] == j )
	    return( _val[a] );
    return( 0.0 );
}


inline double &CRowMatrix::set_no_check( uint32_t i, uint32_t j )
{
    /* Use existing element if it exists */
    for( uint32_t a = _ptr[i]; a < _ptr[i+1]; a++ )
    if( _col[a] == j )
        return( _val[a] );

    /* Reserve new space if necessary */
    if( _nz+1 > _asize )
	reserve( _asize+_m );

    /* Move existing data */
    int movesize = _nz-_ptr[i+1];
    memmove( &_val[_ptr[i+1]+1], &_val[_ptr[i+1]], movesize*sizeof(double) );
    memmove( &_col[_ptr[i+1]+1], &_col[_ptr[i+1]], movesize*sizeof(uint32_t) );

    /* Set new data */
    _val[_ptr[i+1]] = 0.0;
    _col[_ptr[i+1]] = j;

    /* Update pointers */
    _nz++;
    for( uint32_t a = i+1; a < _n+1; a++ )
	_ptr[a]++;

    return( _val[_ptr[i+1]-1] );
}


double CRowMatrix::get_check( uint32_t i, uint32_t j ) const
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( get_no_check( i, j ) );
}


double &CRowMatrix::set_check( uint32_t i, uint32_t j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( set_no_check( i, j ) );
}


void CRowMatrix::set_row( uint32_t i, uint32_t N, const uint32_t *col, const double *val )
{
    if( i >= _n )
	throw( ErrorRange( ERROR_LOCATION, i, _n ) );
    else if( N > _m )
	throw( Error( ERROR_LOCATION, "too many elements" ) );

    /* Reserve new space if necessary */
    uint32_t oldsize = _ptr[i+1] - _ptr[i];
    if( _nz+N-oldsize > _asize )
	reserve( _asize+_m );

    /* Move existing data */
    int offset = N-oldsize;
    int movesize = _nz-_ptr[i+1];
    memmove( &_val[_ptr[i+1]+offset], &_val[_ptr[i+1]], movesize*sizeof(double) );
    memmove( &_col[_ptr[i+1]+offset], &_col[_ptr[i+1]], movesize*sizeof(uint32_t) );

    /* Set new data */
    int err = 0;
    for( uint32_t a = 0; a < N; a++ ) {
	for( uint32_t b = a+1; b < N; b++ ) {
	    if( col[a] == col[b] )
		err = 2;
	}
	_val[_ptr[i]+a] = val[a];
	if( (_col[_ptr[i]+a] = col[a]) >= _m )
	    err = 1;
    }

    /* Update pointers */
    _nz += offset;
    for( uint32_t a = i+1; a < _n+1; a++ )
	_ptr[a] += offset;

    if( err == 1 )
	throw( Error( ERROR_LOCATION, "column index out of range" ) );
    else if( err == 2 )
	throw( Error( ERROR_LOCATION, "repeated column index" ) );
}


void CRowMatrix::construct_add( uint32_t i, uint32_t j, double val )
{
    /* Reserve new space if necessary */
    if( _nz+1 > _asize )
	reserve( _asize+_m );

    /* Set new data */
    _val[_ptr[i+1]] = val;
    _col[_ptr[i+1]] = j;

    /* Update the number of nonzeroes and only the next pointer */
    _nz++;
    _ptr[i+1]++;
    if( i+1 < _n )
	_ptr[i+2] = _ptr[i+1];
}


void CRowMatrix::debug_print( void ) const
{
    std::cout << "n     = " << _n << "\n";
    std::cout << "m     = " << _m << "\n";
    std::cout << "nz    = " << _nz << "\n";
    std::cout << "asize = " << _asize << "\n";

    std::cout << "ptr[] = {";
    for( uint32_t i = 0; i < _n; i++ )
	std::cout << _ptr[i] << ", ";
    std::cout << _ptr[_n] << "}\n";

    std::cout << "col[] = {";
    if( _nz <= 0 ) {
	std::cout << "}\n";
    } else {
	for( uint32_t i = 0; i < _nz-1; i++ )
	    std::cout << _col[i] << ", ";
	std::cout << _col[_nz-1] << "}\n";
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


void CRowMatrix::multiply_by_vector( Vector &x, const Vector &b ) const
{
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x.resize( _n );
    x.clear();

    double sum;
    for( uint32_t i = 0; i < _n; i++ ) {
	sum = 0;
	for( uint32_t a = _ptr[i]; a < _ptr[i+1]; a++ ) {
	    sum += _val[a] * b._val[_col[a]];
	    //std::cout << ""
	}
	x._val[i] = sum;
    }
}


void CRowMatrix::lower_unit_solve( Vector &x, const Vector &b ) const
{
    // Make checks
    if( _n != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix not squrare" ) );
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x = b;

    for( uint32_t i = 0; i < _n; i++ ) {
	for( uint32_t j = _ptr[i]; j < _ptr[i+1]; j++ )
	    x[i] -= _val[j]*x[_col[j]];
    }
}


void CRowMatrix::upper_diag_solve( Vector &x, const Vector &b ) const
{
    // Make checks
    if( _n != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x = b;

    int32_t i, j;
    for( i = _n-1; i >= 0; i-- ) {
	j = _ptr[i+1]-1;
	for( ; j > (int32_t)_ptr[i]; j-- )
	    x[i] -= _val[j] * x[_col[j]];
	x[i] /= _val[j];
    }
}

















