/*! \file error.hpp
 *  \brief Header file for error.hpp
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

#ifndef ERROR_HPP
#define ERROR_HPP 1


#include <string>
#include <cstring>
#include <sstream>
#include <errno.h>


template <class T>
inline std::string to_string( const T& t )
{
    std::stringstream ss;
    ss << t;
    return( ss.str() );
}


/*! \brief Macro for setting error location when throwing errors.
 */
#define ERROR_LOCATION ErrorLocation( __FILE__, __LINE__, __func__ )


/*! \brief %Error location class.
 *
 *  Container to store the location (source file name, line number and
 *  function name) where the error happened. Used for debugging
 *  purposes.
 */
struct ErrorLocation {
    ErrorLocation() {}
    ErrorLocation( const char *file, int line, const char *func ) :
	_file(file), _line(line), _func(func) {}
    const char *_file;
    int         _line;
    const char *_func;
};


/*! \brief Base error class.
 */
struct Error {
    ErrorLocation  _loc;
    std::string    _error_str;
    Error() {}
    Error( ErrorLocation loc, const std::string &str ) 
	: _loc(loc), _error_str(str) {}
};


/*! \brief %Error class for memory allocation errors.
 */
struct ErrorNoMem : public Error {
    ErrorNoMem( ErrorLocation loc ) 
	: Error( loc, "memory allocation error" ) {}
};


/*! \brief %Error class for C-style errno errors.
 */
struct ErrorErrno : public Error {
    int _errno;
    ErrorErrno( ErrorLocation loc ) : _errno(errno) {
	_loc = loc;
	char buf[128];
	_error_str = strerror_r( _errno, buf, 128 );
    }
};


/*! \brief %Error class to use if requested feature is unimplemented.
 */
struct ErrorUnimplemented : public Error {
    ErrorUnimplemented( ErrorLocation loc ) 
	: Error( loc, "feature unimplemented" ) {}
    ErrorUnimplemented( ErrorLocation loc, const std::string &str ) 
	: Error( loc, str ) {}
};


/*! \brief %Error class for dimension mismatch errors.
 */
struct ErrorDim : public Error {
    ErrorDim( ErrorLocation loc )
	: Error( loc, "dimension mismatch" ) {}
    ErrorDim( ErrorLocation loc, const std::string &str )
	: Error( loc, str ) {}
};


/*! \brief %Error class for index range checking errors.
 */
struct ErrorRange : public Error {
    ErrorRange( ErrorLocation loc, uint32_t i, uint32_t n, uint32_t j, uint32_t m ) {
	std::ostringstream ss;
	ss << "index out of range ( " << i << " >= " << n << " || ";
	ss << j << " >= " << m << " )";
	_error_str = ss.str();
	_loc = loc;
    }
    ErrorRange( ErrorLocation loc, uint32_t i, uint32_t n ) {
	std::ostringstream ss;
	ss << "index out of range ( " << i << " >= " << n << " )";
	_error_str = ss.str();
	_loc = loc;
    }
};


#endif














