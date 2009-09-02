#include "config.hpp"
#include "coordmatrix.hpp"
#include "ccolmatrix.hpp"
#include "crowmatrix.hpp"
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iomanip>



inline void CoordMatrix::allocate( void )
{
    if( _asize == 0 ) {
	_row = _col = NULL;
	_val = NULL;
    } else {
	if( !(_row = (uint32_t *)malloc( _asize*sizeof(uint32_t) )) ) {
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	if( !(_col = (uint32_t *)malloc( _asize*sizeof(uint32_t) )) ) {
	    free( _row );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	if( !(_val = (double *)malloc( _asize*sizeof(double) )) ) {
	    free( _row );
	    free( _col );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
    }
}


inline void CoordMatrix::reallocate( void )
{
    if( _asize == 0 ) {
	free( _row );	
	free( _col );
	free( _val );
	_row = _col = NULL;
	_val = NULL;
    } else {
	uint32_t *tmp;
	if( !(tmp = (uint32_t *)realloc( _row, _asize*sizeof(uint32_t) )) ) {
	    free( _row );
	    free( _col );
	    free( _val );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_row = tmp;
	if( !(tmp = (uint32_t *)realloc( _col, _asize*sizeof(uint32_t) )) ) {
	    free( _row );
	    free( _col );
	    free( _val );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_col = tmp;
	double *tmp2;
	if( !(tmp2 = (double *)realloc( _val, _asize*sizeof(double) )) ) {
	    free( _row );
	    free( _col );
	    free( _val );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_val = tmp2;
    }
}


CoordMatrix::CoordMatrix( uint32_t n, uint32_t m )
{
    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    _row   = NULL;
    _col   = NULL;
    _val   = NULL;    
}


CoordMatrix::CoordMatrix( uint32_t n, uint32_t m, uint32_t nz, 
			  const uint32_t *row, const uint32_t *col, const uint32_t *val )
{
    _n     = n;
    _m     = m;
    _nz    = nz;
    _asize = nz;
    allocate();

    memcpy( _row, row, _nz*sizeof(uint32_t) );
    memcpy( _col, col, _nz*sizeof(uint32_t) );
    memcpy( _val, val, _nz*sizeof(double) );
}


CoordMatrix::CoordMatrix( const CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();

    memcpy( _row, mat._row, _nz*sizeof(uint32_t) );
    memcpy( _col, mat._col, _nz*sizeof(uint32_t) );
    memcpy( _val, mat._val, _nz*sizeof(double) );
}


CoordMatrix &CoordMatrix::operator=( const CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    memcpy( _row, mat._row, _nz*sizeof(uint32_t) );
    memcpy( _col, mat._col, _nz*sizeof(uint32_t) );
    memcpy( _val, mat._val, _nz*sizeof(double) );

    return( *this );
}


CoordMatrix::CoordMatrix( const class CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();

    uint32_t i, j;
    for( i = j = 0; i < _n; i++ ) {
	uint32_t e = mat._ptr[i+1];
	for( ; j < e; j++ ) {
	    _row[j] = i;
	    _col[j] = mat._col[j];
	    _val[j] = mat._val[j];
	}
    }
}


CoordMatrix &CoordMatrix::operator=( const CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    uint32_t i, j;
    for( i = j = 0; i < _n; i++ ) {
	uint32_t e = mat._ptr[i+1];
	for( ; j < e; j++ ) {
	    _row[j] = i;
	    _col[j] = mat._col[j];
	    _val[j] = mat._val[j];
	}
    }

    return( *this );
}


CoordMatrix::CoordMatrix( const class CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    allocate();

    uint32_t i, j;
    for( i = j = 0; j < _m; j++ ) {
	uint32_t e = mat._ptr[j+1];
	for( ; i < e; i++ ) {
	    _row[i] = mat._row[i];
	    _col[i] = j;
	    _val[i] = mat._val[i];
	}
    }
}


CoordMatrix &CoordMatrix::operator=( const CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    uint32_t i, j;
    for( i = j = 0; j < _m; j++ ) {
	uint32_t e = mat._ptr[j+1];
	for( ; i < e; i++ ) {
	    _row[i] = mat._row[i];
	    _col[i] = j;
	    _val[i] = mat._val[i];
	}
    }

    return( *this );
}


CoordMatrix::~CoordMatrix()
{
    free( _row );
    free( _col );
    free( _val );
}


void CoordMatrix::resize( uint32_t n, uint32_t m )
{
    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    free( _row );
    free( _col );
    free( _val );
    _row = NULL;
    _col = NULL;
    _val = NULL;
}


void CoordMatrix::merge( CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._asize;
    free( _row );
    free( _col );
    free( _val );
    _row = mat._row;
    _col = mat._col;
    _val = mat._val;

    mat._n     = 0;
    mat._m     = 0;
    mat._nz    = 0;
    mat._asize = 0;
    mat._row   = NULL;
    mat._col   = NULL;
    mat._val   = NULL;
}


void CoordMatrix::clear( void )
{
    _nz    = 0;
    _asize = 0;
    free( _row );
    free( _col );
    free( _val );
    _row = NULL;
    _col = NULL;
    _val = NULL;
}


inline void CoordMatrix::clear_no_check( uint32_t i, uint32_t j )
{

}


void CoordMatrix::clear_check( uint32_t i, uint32_t j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    clear_no_check( i, j );
}


void CoordMatrix::reserve( uint32_t size )
{
    if( size > _asize ) {
	_asize = size;
	reallocate();
    }
}


void CoordMatrix::set_nz( uint32_t nz )
{
    if( nz > _asize ) {
	_asize = nz;
	reallocate();
    }
    _nz = nz;
}


void CoordMatrix::order_ascending_row_column( void )
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void CoordMatrix::order_ascending_column_row( void )
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void CoordMatrix::debug_print( void ) const
{
    std::cout << "n     = " << _n << "\n";
    std::cout << "m     = " << _m << "\n";
    std::cout << "nz    = " << _nz << "\n";
    std::cout << "asize = " << _asize << "\n";

    std::cout << "row[] = {";
    if( _nz <= 0 ) {
	std::cout << "}\n";
    } else {
	for( uint32_t i = 0; i < _nz-1; i++ )
	    std::cout << _row[i] << ", ";
	std::cout << _row[_nz-1] << "}\n";
    }

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


inline double CoordMatrix::get_no_check( uint32_t i, uint32_t j ) const
{
    for( uint32_t a = 0; a < _nz; a++ )
	if( _row[a] == i && _col[a] == j  )
	    return( _val[a] );
    return( 0.0 );
}


inline double &CoordMatrix::set_no_check( uint32_t i, uint32_t j )
{
    /* Use existing element if it exists */
    for( uint32_t a = 0; a < _nz; a++ )
	if( _row[a] == i && _col[a] == j  )
	    return( _val[a] );

    /* Reserve new space if necessary */
    if( _nz >= _asize )
	reserve( _asize+_n );

    /* Set new data */
    _row[_nz] = i;
    _col[_nz] = j;
    _val[_nz] = 0.0;
    _nz += 1;

    return( _val[_nz-1] );
}


double CoordMatrix::get_check( uint32_t i, uint32_t j ) const
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( get_no_check( i, j ) );
}


double &CoordMatrix::set_check( uint32_t i, uint32_t j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( set_no_check( i, j ) );
}


void CoordMatrix::set_no_duplicate_check( uint32_t i, uint32_t j, double val )
{
    /* Reserve new space if necessary */
    if( _nz >= _asize )
	reserve( _asize+_n );

    /* Set new data */
    _row[_nz] = i;
    _col[_nz] = j;
    _val[_nz] = val;
    _nz += 1;
}



/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */


void CoordMatrix::multiply_by_vector( Vector &x, const Vector &b ) const
{
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "Matrix dimension does not match vector" ) );

    x.resize( _n );
    x.clear();

    for( uint32_t i = 0; i < _nz; i++ )
	x._val[_row[i]] += _val[i] * b._val[_col[i]];
}


void CoordMatrix::lower_unit_solve( Vector &x, const Vector &b ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void CoordMatrix::upper_diag_solve( Vector &x, const Vector &b ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}

















