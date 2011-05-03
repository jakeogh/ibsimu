/*! \file ibsimu.hpp
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


#ifndef IBSIMU_HPP
#define IBSIMU_HPP 1


#include <string>
#include <iostream>
#include <fstream>


class Timer;


/*! \brief Main class for %IBSimu.
 *
 *  Used to store global settings. One instance of the class is
 *  initialized globally with the name \a ibsimu.
 */
class IBSimu 
{
    Timer        *_t;

    bool          _hello;
    int           _verbose_output;
    int           _threadcount;

    bool          _is_cout;         // True if vout is std::cout
    std::ostream *_vout;            // Verbose output stream

    std::ofstream _fout;            // Verbose output file

    IBSimu( const IBSimu &ibs ) : _vout(ibs._vout) {}

    const IBSimu &operator=( const IBSimu &ibs ) { return( *this ); }

public:

    /*! \brief Default constructor.
     */
    IBSimu();

    /*! \brief Default destructor.
     */
    ~IBSimu();

    /*! \brief Set verbose output to stream \a vout
     *
     *  Returns a reference to the old output stream.
     */
    std::ostream &set_vout( std::ostream &vout );

    /*! \brief Set verbose output to file \a filename.
     *
     *  Returns a reference to the old output stream. If the output
     *  stream is redefined, the file is kept open in the
     *  background. IBSimu can only have one output stream opened at
     *  time with this function. The file is closed when the IBSimu
     *  object is destructed.
     */
    std::ostream &set_vout( const std::string &filename );

    /*! \brief Get a reference to verbose output stream.
     */
    std::ostream &vout( void );

    /*! \brief Return if verbose output stream is std::cout.
     */
    bool vout_is_cout();

    /*! \brief Set verbosity level.
     */
    void set_verbose_output( int level );

    /*! \brief Get verbosity level.
     */
    int get_verbose_output( void ) { return( _verbose_output ); }

    /*! \brief Set the number of threads used for calculation.
     */
    void set_thread_count( int threadcount );

    /*! \brief Get the number of threads used for calculation.
     */
    int get_thread_count( void ) { return( _threadcount ); }

    /*! \brief Halt execution
     *
     *  This function is called by the error handler in case of SIGTERM.
     */
    void halt( void );
};


/*! \brief Global instance of class %IBSimu.
 */
extern IBSimu ibsimu;


#endif
