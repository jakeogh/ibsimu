/*! \file error.hpp
 *  \brief %Error classes and handling
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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
 * taneli.kalvas@jyu.fi.
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
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <errno.h>

#include <memory.h>
#include <unistd.h>
#include <ucontext.h>
#include <dlfcn.h>
#ifndef NO_CPP_DEMANGLE
#include <cxxabi.h>
#ifdef __cplusplus
using __cxxabiv1::__cxa_demangle;
#endif
#endif

#define sigsegv_outp(x, ...)    fprintf(stderr, x "\n", ##__VA_ARGS__)

#if defined(REG_RIP)
# define SIGSEGV_STACK_IA64
# define REGFORMAT "%016lx"
#elif defined(REG_EIP)
# define SIGSEGV_STACK_X86
# define REGFORMAT "%08x"
#else
# define SIGSEGV_STACK_GENERIC
# define REGFORMAT "%x"
#endif


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
    ErrorLocation() 
	: _file(NULL), _line(0), _func(0) {}
    ErrorLocation( const char *file, int line, const char *func )
	: _file(file), _line(line), _func(func) {}
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
    Error( const std::string &str ) 
	: _error_str(str) {}
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
    int _ierrno;
    ErrorErrno( ErrorLocation loc ) : _ierrno(errno) {
	_loc = loc;
#if defined(WIN32) || defined(__MINGW32__)
	_error_str = "errno " + to_string(_ierrno);
#else
	char buf[1024];
	strerror_r( _ierrno, buf, 1024 );
	_error_str = buf;
#endif
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


/*! \brief Signal handler
 */
struct SignalHandler {

    static void signal_handler_SIGSEGV(int signum, siginfo_t* info, void*ptr) {
	static const char *si_codes[3] = {"", "SEGV_MAPERR", "SEGV_ACCERR"};

	int i, f = 0;
	ucontext_t *ucontext = (ucontext_t*)ptr;
	Dl_info dlinfo;
	void **bp = 0;
	void *ip = 0;

	sigsegv_outp("Segmentation Fault!");
	sigsegv_outp("info.si_signo = %d", signum);
	sigsegv_outp("info.si_errno = %d", info->si_errno);
	sigsegv_outp("info.si_code  = %d (%s)", info->si_code, si_codes[info->si_code]);
	sigsegv_outp("info.si_addr  = %p", info->si_addr);
	for(i = 0; i < NGREG; i++)
		sigsegv_outp("reg[%02d]       = 0x" REGFORMAT, i, ucontext->uc_mcontext.gregs[i]);

#ifndef SIGSEGV_NOSTACK
#if defined(SIGSEGV_STACK_IA64) || defined(SIGSEGV_STACK_X86)
#if defined(SIGSEGV_STACK_IA64)
	ip = (void*)ucontext->uc_mcontext.gregs[REG_RIP];
	bp = (void**)ucontext->uc_mcontext.gregs[REG_RBP];
#elif defined(SIGSEGV_STACK_X86)
	ip = (void*)ucontext->uc_mcontext.gregs[REG_EIP];
	bp = (void**)ucontext->uc_mcontext.gregs[REG_EBP];
#endif

	sigsegv_outp( "Stack trace:" );
	while( bp && ip ) {
	    if( !dladdr( ip, &dlinfo ) )
		break;
	    
	    const char *symname = dlinfo.dli_sname;
	    
#ifndef NO_CPP_DEMANGLE
	    int status;
	    char *tmp = __cxa_demangle( symname, NULL, 0, &status );

	    if( status == 0 && tmp )
		symname = tmp;
#endif

	    sigsegv_outp( "% 2d: %p <%s+%lu> (%s)",
			  ++f,
			  ip,
			  symname,
			  (unsigned long)ip - (unsigned long)dlinfo.dli_saddr,
			  dlinfo.dli_fname );
	    
#ifndef NO_CPP_DEMANGLE
	    if( tmp )
		free( tmp );
#endif

	    if( dlinfo.dli_sname && !strcmp( dlinfo.dli_sname, "main" ) )
		break;

	    ip = bp[1];
	    bp = (void**)bp[0];
	}
#else
	sigsegv_outp( "Stack trace (non-dedicated):" );
	sz = backtrace( bt, 20 );
	strings = backtrace_symbols( bt, sz );
	for( i = 0; i < sz; ++i )
	    sigsegv_outp( "%s", strings[i] );
#endif
	sigsegv_outp( "End of stack trace." );
#else
	sigsegv_outp( "Not printing stack strace." );
#endif
	_exit( -1 );
    }

};


#endif
