#include <iomanip>
#include "matrix.hpp"


MatrixMulVec Matrix::operator*( const class Vector &vec ) const
{
    if( columns() != vec.size() )
	throw( ErrorDim( ERROR_LOCATION ) );    
    MatrixMulVec res( *this, vec );
    return( res );
}


std::ostream &operator<<( std::ostream &os, const Matrix &mat )
{
    for( int i = 0; i < mat.rows(); i++ ) {
	for( int j = 0; j < mat.columns(); j++ ) {
	    os << std::setw(6) << to_string(mat.get(i,j)).substr(0,6) << " ";
	}
	if( i != mat.rows()-1 )
	    os << "\n";
    }
    return( os );
}















