/*! \file error.hpp
 *  \brief Header file error handling.
 */


#include <string>
#include <cstring>
#include <sstream>
#include <errno.h>


#ifndef ERROR_HPP
#define ERROR_HPP 1


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














