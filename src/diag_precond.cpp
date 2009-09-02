#include "diag_precond.hpp"
#include "error.hpp"


Diag_Precond::Diag_Precond( const Matrix &mat )
{
    // Make checks
    if( mat.columns() != mat.rows() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not squrare" ) );

    diag.resize( mat.rows() );
    for( uint32_t i = 0; i < mat.rows(); i++ ) {
	if( mat.get(i,i) == 0 )
	    throw( Error( ERROR_LOCATION, "Zero element on diagonal at i=" + to_string(i) ) );
	diag[i] = 1.0/mat.get(i,i);
    }
}


void Diag_Precond::solve( Vector &x, const Vector &b ) const
{
    if( b.size() != diag.size() )
	throw( ErrorDim( ERROR_LOCATION, "Matrix dimension does not match vector" ) );
    x.resize( diag.size() );

    for( uint32_t i = 0; i < diag.size(); i++ )
	x[i] = diag[i]*b[i];
}
















