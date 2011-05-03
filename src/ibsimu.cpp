/*! \file ibsimu.cpp
 *  \brief Ion Beam Simulator global settings
 */

/* Copyright (c) 2010-2011 Taneli Kalvas. All rights reserved.
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
#include <signal.h>
#include "ibsimu.hpp"
#include "timer.hpp"
#include "error.hpp"



IBSimu::IBSimu()
    : _hello(false), _verbose_output(0), _threadcount(1), _is_cout(true), _vout(&std::cout)
{
#ifdef _GNU_SOURCE
    // Set a catch for segmentation fault
    struct sigaction act_sigsegv;
    act_sigsegv.sa_sigaction = SignalHandler::signal_handler_SIGSEGV;
    sigemptyset( &act_sigsegv.sa_mask );
    act_sigsegv.sa_flags = SA_SIGINFO;
    sigaction( SIGSEGV, &act_sigsegv, NULL );
#endif

    // Set a catch for terminate/kill/int
    struct sigaction act_sigterm;
    act_sigterm.sa_sigaction = SignalHandler::signal_handler_SIGTERM;
    sigemptyset( &act_sigterm.sa_mask );
    act_sigterm.sa_flags = SA_SIGINFO;
    sigaction( SIGTERM, &act_sigterm, NULL );
    sigaction( SIGQUIT, &act_sigterm, NULL );
    sigaction( SIGINT, &act_sigterm, NULL );

    // Start timer for whole simulation
    _t = new Timer;
}


IBSimu::~IBSimu()
{
    // End timer
    _t->stop();

    if( _verbose_output ) {
	vout() << "Ending simulation\n";
	vout() << "  time used = " << *_t << "\n";
	vout() << std::flush;
    }

    // Close the verbose output file if it is open
    if( _fout.is_open() )
	_fout.close();

    delete _t;
}


void IBSimu::halt( void )
{
    // End timer
    _t->stop();

    if( _verbose_output ) {
	vout() << "Ending simulation\n";
	vout() << "  time used = " << *_t << "\n";
	vout() << std::flush;
    }

    // Close the verbose output file if it is open
    if( _fout.is_open() )
	_fout.close();
}


std::ostream &IBSimu::set_vout( std::ostream &vout )
{
    if( vout == std::cout )
	_is_cout = true;
    std::ostream &vout_old = *_vout;
    _vout = &vout;
    return( vout_old );
}


std::ostream &IBSimu::set_vout( const std::string &filename )
{
    if( _fout.is_open() )
	throw( Error( ERROR_LOCATION, "Trying to open verbose output file with previous file still open" ) );
	
    std::ostream &vout_old = *_vout;
    _fout.open( filename.c_str() );
    if( !_fout.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    _vout = &_fout;
    _is_cout = false;
    return( vout_old );
}


std::ostream &IBSimu::vout( void ) 
{ 
    return( *_vout ); 
}


bool IBSimu::vout_is_cout()
{
    return( _is_cout ); 
}


void IBSimu::set_verbose_output( int level )
{
    // Set verbosity level
    _verbose_output = level;

    // If setting to verbose mode and no greeting has yet been shown
    if( level > 0 && !_hello ) {
	_hello = true;
	vout() << "Ion Beam Simulator " << VERSION << " (compiled " __DATE__ ")\n";
    }
}


void IBSimu::set_thread_count( int threadcount ) 
{
    if( threadcount <= 0 )
	throw( Error( ERROR_LOCATION, "invalid parameter" ) );
    _threadcount = threadcount;
}


/* Global instance */
IBSimu ibsimu;

