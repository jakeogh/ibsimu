/*! \file color.hpp
 *  \brief Header file for color.hpp
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * tvkalvas@cc.jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
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

















