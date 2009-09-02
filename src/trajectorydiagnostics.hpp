/*! \file trajectorydiagnostics.hpp
 *  \brief Header file for trajectory diagnostics
 */


#ifndef TRAJECTORYDIAGNOSTICS_HPP
#define TRAJECTORYDIAGNOSTICS_HPP 1


#include <vector>
#include "types.hpp"


/*! \brief Class for trajectory diagnostic data column.
 */
class TrajectoryDiagnosticColumn {

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
class TrajectoryDiagnosticData {

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
 */
class Emittance
{
    double _wsum;

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

    /*! \brief Constructor for emittance statistics from trajectory
     *  diagnostic data columns \a x, \a xp and weight \a w.
     */
    Emittance( const std::vector<double> &x,
	       const std::vector<double> &xp,
	       const std::vector<double> &w );

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


#endif













