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


#include "config.h"
#include <iostream>
#include <iomanip>
#include "error.hpp"

#ifdef _GNU_SOURCE
#include <execinfo.h>
#endif

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




ExceptionTracer::ExceptionTracer() 
{
#ifdef _GNU_SOURCE
    _tracecount = backtrace( _traceaddress, 25 );
#else
    _tracecount = 0;
#endif
}


#ifdef _GNU_SOURCE
bool print_trace_addr( std::ostream &os, int i, void *addr )
{
    Dl_info dlinfo;
    if( !dladdr( addr, &dlinfo ) )
	return( true );

    const char *symname = dlinfo.dli_sname;
    if( !dlinfo.dli_sname )
	return( true );

#ifndef NO_CPP_DEMANGLE
    int status;
    char *tmp = __cxa_demangle( symname, NULL, 0, &status );

    if( status == 0 && tmp )
	symname = tmp;
#endif

    os << std::setw(2) << i;
    os << " 0x" << std::hex << addr;
    os << " <" << symname << "+0x" << std::hex 
       << (unsigned long)addr - (unsigned long)dlinfo.dli_saddr << ">";
    os << " (" << dlinfo.dli_fname << ")\n";

    if( !strcmp( dlinfo.dli_sname, "main" ) || !strcmp( symname, "main" ) )
	return( true );

    return( false );
}
#endif


void ExceptionTracer::print_trace( std::ostream &os ) 
{
#ifdef _GNU_SOURCE
    os << "Stack trace:\n";
    for( int i = 0; i < _tracecount; i++ ) {
	if( print_trace_addr( os, i, _traceaddress[i] ) )
	    break;
    }
    os << "End of stack trace.\n";
#else
    os << "No backtrace capability\n";
#endif
}


#ifdef _GNU_SOURCE
void SignalHandler::signal_handler_SIGSEGV( int signum, siginfo_t *info, void *ptr )
{
    static const char *si_codes[3] = {"", "SEGV_MAPERR", "SEGV_ACCERR"};
    ucontext_t *ucontext = (ucontext_t*)ptr;

    std::cerr << "Segmentation Fault!\n";
    std::cerr << "info.si_signo = " << signum << "\n";
    std::cerr << "info.si_errno = " << info->si_errno << "\n";
    std::cerr << "info.si_code  = " << info->si_code << "(" << si_codes[info->si_code] << ")\n";
    std::cerr << "info.si_addr  = 0x" << std::hex << info->si_addr << "\n";
    for( int i = 0; i < NGREG; i++ )
	std::cerr << "reg[" << std::setw(2) << i << "]       = 0x" 
		  << std::hex << ucontext->uc_mcontext.gregs[i] << "\n";

#ifndef SIGSEGV_NOSTACK
#if defined(SIGSEGV_STACK_IA64) || defined(SIGSEGV_STACK_X86)

#if defined(SIGSEGV_STACK_IA64)
    void  *ip = (void*)ucontext->uc_mcontext.gregs[REG_RIP];
    void **bp = (void**)ucontext->uc_mcontext.gregs[REG_RBP];
#elif defined(SIGSEGV_STACK_X86)
    void  *ip = (void*)ucontext->uc_mcontext.gregs[REG_EIP];
    void **bp = (void**)ucontext->uc_mcontext.gregs[REG_EBP];
#endif

    std::cerr << "Stack trace:\n";
    int i = 0;
    while( bp && ip ) {

	if( print_trace_addr( std::cerr, i, ip ) )
	    break;

	ip = bp[1];
	bp = (void**)bp[0];
	i++;
    }
#else
    std::cerr << "Stack trace (non-dedicated):\n";
    sz = backtrace( bt, 20 );
    strings = backtrace_symbols( bt, sz );
    for( int i = 0; i < sz; ++i )
	std::cerr << strings[i] << "\n";
#endif
    std::cerr << "End of stack trace.\n";
#else
    std::cerr << "Not printing stack trace.\n";
#endif
    _exit( -1 );
}
#endif



