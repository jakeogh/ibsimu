/*! \file trajectorydiagnostics.hpp
 *  \brief Header file for trajectorydiagnostics.hpp
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

#ifndef TRAJECTORYDIAGNOSTICS_HPP
#define TRAJECTORYDIAGNOSTICS_HPP 1


#include <vector>
#include "histogram.hpp"
#include "types.hpp"


/*! \brief Class for trajectory diagnostic data column.
 */
class TrajectoryDiagnosticColumn 
{

    trajectory_diagnostic_e _diag; /*!< \brief Type of diagnostic in data. */
    std::vector<double>     _data; /*!< \brief Vector of diagnostic data. */

public:

    TrajectoryDiagnosticColumn( trajectory_diagnostic_e diag ) 
	: _diag(diag) {}

    /*! \brief Add mirrored trajectory diagnostic data to the column.
     *
     *  Mirror data along plane \a axis = \a level. The mirrored data
     *  is added to the data column. This function is used to get a
     *  complete data set of a particle beam, of which only half (or
     *  quarter) has been simulated.
     */
    void mirror( coordinate_axis_e axis, double level );

    void add_data( double x ) {
	_data.push_back( x );
    }

    std::vector<double> &data( void ) { return( _data ); }

    const std::vector<double> &data( void ) const { return( _data ); }

    size_t size( void ) const { return( _data.size() ); }

    trajectory_diagnostic_e diagnostic( void ) const { return( _diag ); }

    const double &operator()( size_t i ) const { return( _data[i] ); }

    double &operator()( size_t i ) { return( _data[i] ); }

    const double &operator[]( size_t i ) const { return( _data[i] ); }

    double &operator[]( size_t i ) { return( _data[i] ); }

};


/*! \brief Class for trajectory diagnostic data.
 */
class TrajectoryDiagnosticData 
{

    std::vector<TrajectoryDiagnosticColumn> _column; /*!< \brief Vector of diagnostic data columns. */

public:

    TrajectoryDiagnosticData() {}

    TrajectoryDiagnosticData( std::vector<trajectory_diagnostic_e> diag ) {
	for( size_t a = 0; a < diag.size(); a++ )
	    _column.push_back( TrajectoryDiagnosticColumn( diag[a] ) );
    }

    void mirror( coordinate_axis_e axis, double level ) {
	for( size_t a = 0; a < _column.size(); a++ )
	    _column[a].mirror( axis, level );
    }

    void clear() {
	_column.clear();
    }
    
    void add_data_column( trajectory_diagnostic_e diag ) {
	_column.push_back( TrajectoryDiagnosticColumn( diag ) );
    }

    size_t diag_size() const {
	return( _column.size() );
    }

    size_t traj_size() const {
	if( _column.size() > 0 )
	    return( _column[0].size() );
	return( 0 );
    }

    trajectory_diagnostic_e diagnostic( size_t diag ) const {
	return( _column[diag].diagnostic() ); 
    }

    const TrajectoryDiagnosticColumn &operator()( size_t diag ) const {
	return( _column[diag] );
    }

    TrajectoryDiagnosticColumn &operator()( size_t diag ) {
	return( _column[diag] );
    }

    const double &operator()( size_t traj, size_t diag ) const {
	return( _column[diag](traj) );
    }

    double &operator()( size_t traj, size_t diag ) {
	return( _column[diag](traj) );
    }

    void add_data( size_t diag, double x ) {
	_column[diag].add_data( x );
    }
};


/*! \brief Class for emittance statistics.
 *
 *  %Emittance class does a statistical analysis on the particle
 *  distribution and it calculates averages \f$ <x> \f$ and 
 *  \f$ <x'> \f$ and the expectation values \f$ <x^2> \f$, 
 *  \f$ <x'^2> \f$ and \f$ <x x'> \f$. From these it calculates 
 *  the rms-emittance
 *  \f[ \epsilon = \sqrt{ <x^2><x'^2> - <x x'>^2 } \f]
 *  and the Twiss parameters
 *  \f[ \alpha = \frac{-<x x'>}{\epsilon}, \beta = \frac{<x^2>}{\epsilon}, \gamma = \frac{<x'^2>}{\epsilon} \f]
 *  In addition to these physical values, the class calculates the angle of the ellipse
 *  \f[ \theta = \frac{1}{2} \arctan{\left( \frac{2\alpha}{\beta - \gamma} \right)} \f]
 *  and the half-axis lengths
 *  \f[ r_1 = \sqrt{\frac{\epsilon}{2}} ( \sqrt{H+1} + \sqrt{H-1} ) \f]
 *  \f[ r_2 = \sqrt{\frac{\epsilon}{2}} ( \sqrt{H+1} - \sqrt{H-1} ), \f]
 *  where
 *  \f[ H = \frac{\beta + \gamma}{2} \f]
 *  
 */
class Emittance
{
protected:

    double _Isum;

    double _xave;
    double _xpave;

    double _x2;
    double _xp2;
    double _xxp;

    double _alpha;
    double _beta;
    double _gamma;
    double _epsilon;

    double _angle;
    double _rmajor;
    double _rminor;

public:

    /*! \brief Default constructor for emittance statistics
     */
    Emittance();

    /*! \brief Constructor for emittance statistics from trajectory
     *  diagnostic data columns \a x, \a xp and current \a I.
     */
    Emittance( const std::vector<double> &x,
	       const std::vector<double> &xp,
	       const std::vector<double> &I );

    /*! \brief Constructor for emittance statistics from trajectory
     *  diagnostic data columns \a x, \a xp, assuming even weights.
     */
    Emittance( const std::vector<double> &x,
	       const std::vector<double> &xp );

    double xave( void ) { return( _xave ); }
    double xpave( void ) { return( _xpave ); }

    double alpha( void ) { return( _alpha ); }
    double beta( void ) { return( _beta ); }
    double gamma( void ) { return( _gamma ); }
    double epsilon( void ) { return( _epsilon ); }

    double angle( void ) { return( _angle ); }
    double rmajor( void ) { return( _rmajor ); }
    double rminor( void ) { return( _rminor ); } 
};



/*! \brief Class for emittance conversion from (r,r') to (x,x')
 *
 *  The emittance converted takes the 
 */
class EmittanceConv : public Emittance
{

    Histogram2D *_grid;

public:

    /*! \brief Constructor for \a (x,x') emittance data and statistics
     *  from \a (r,r') data.
     *
     *  Reads particle diagnostic data arrays for \a r (radius), \a rp
     *  (radial angle), \a ap (skew angle) and \a I (current) and
     *  builds \a (x,x') data in a grid array of size \a n by \a
     *  m. Here the skew angle is \f$ \frac{r\omega}{v_z} \f$, where
     *  \f$ v_z \f$ is the velocity to the direction of beam
     *  propagation. The conversion is based on rotating each
     *  trajectory diagnostic points around the axis in 100 steps (to
     *  be made adjustable?).
     *
     *  The emittance statistics is built using original data and not
     *  the gridded data for maximized precision.
     */
    EmittanceConv( int n, int m,
		   const std::vector<double> &r,
		   const std::vector<double> &rp,
		   const std::vector<double> &ap,
		   const std::vector<double> &I );

    /*! \brief Destructor for emittance converter.
     */
    ~EmittanceConv();

    /*! \brief Get a const reference to histogram built.
     */
    const Histogram2D &histogram( void ) const { return( *_grid ); }

    /*! \brief Free emittance histogram.
     */
    void free_histogram( void ) { delete _grid; _grid = NULL; }
};


#endif


















