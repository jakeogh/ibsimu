/*! \file color.hpp
 *  \brief Header file for plot color.
 */


#ifndef COLOR_HPP
#define COLOR_HPP 1


#include <iostream>
#include "error.hpp"


/*! \brief %Color class for plotting.
 */
class Color {

    double _comp[4]; /*!< \brief Color components rgba. */

public:

    /*! \brief Default color constructor (black).
     */
    Color() {
	_comp[0] = 0.0;
	_comp[1] = 0.0;
	_comp[2] = 0.0;
	_comp[3] = 1.0;
    }

    Color( double r, double g, double b ) {
	_comp[0] = r;
	_comp[1] = g;
	_comp[2] = b;
	_comp[3] = 1.0;
    }
    
    Color( double r, double g, double b, double a ) {
	_comp[0] = r;
	_comp[1] = g;
	_comp[2] = b;
	_comp[3] = a;
    }
    
    const double &operator[]( int i ) const {
	return( _comp[i] );
    }

    double &operator[]( int i ) {
	return( _comp[i] );
    }

    Color operator-( const Color &c ) const {
	return( Color( _comp[0]-c[0],
		       _comp[1]-c[1],
		       _comp[2]-c[2],
		       _comp[3]-c[3] ) );
    }

    Color operator+( const Color &c ) const {
	return( Color( _comp[0]+c[0],
		       _comp[1]+c[1],
		       _comp[2]+c[2],
		       _comp[3]+c[3] ) );
    }
    
    Color operator*( double x ) const {
	return( Color( x*_comp[0],
		       x*_comp[1],
		       x*_comp[2],
		       x*_comp[3] ) );
    }

    friend Color operator*( double x, const Color &c );
    friend std::ostream &operator<<( std::ostream &os, const Color &c );
};


inline Color operator*( double x, const Color &c )
{
    return( Color( x*c[0], x*c[1], x*c[2], x*c[3] ) );
}


inline std::ostream &operator<<( std::ostream &os, const Color &c ) 
{
    os << to_string(c[0]) << " "
       << to_string(c[1]) << " "
       << to_string(c[2]) << " "
       << to_string(c[3]);
    return( os );
}


#endif















