/*! \file particlediagplot.hpp
 *  \brief Header file for particlediagplot.hpp
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

#ifndef PARTICLEDIAGPLOT_HPP
#define PARTICLEDIAGPLOT_HPP 1


#include "frame.hpp"
#include "geometry.hpp"
#include "particledatabase.hpp"
#include "types.hpp"
#include "histogram.hpp"
#include "trajectorydiagnostics.hpp"

#include "xygraph.hpp"
#include "colormap.hpp"



enum particle_diag_plot_type_e {
    PARTICLE_DIAG_PLOT_NONE = 0,
    PARTICLE_DIAG_PLOT_SCATTER,
    PARTICLE_DIAG_PLOT_HISTO2D,
    PARTICLE_DIAG_PLOT_HISTO1D
};


/*! \brief %Particle diagnostic plot.
 *
 *  Two dimensional histograms have particle trajectory currents
 *  always taken in account.  Similarly profile plot
 *  (PARTICLE_DIAG_PLOT_HISTO1D) always takes in account the
 *  trajectory current. In cylindrical symmetry cases the output is
 *  scaled to have constant area per histogram bin.  One dimensional
 *  Emittance plots in (r,r') space are scaled to have constant area
 *  per histogram bin.
 *  
 */
class ParticleDiagPlot {
    
    Frame                     *_frame;

    const Geometry            *_geom;
    const ParticleDataBase    *_pdb;

    coordinate_axis_e          _axis;
    double                     _level;

    particle_diag_plot_type_e  _type;
    trajectory_diagnostic_e    _diagx;
    trajectory_diagnostic_e    _diagy;
    trajectory_diagnostic_e    _diagz;

    int                        _pdb_it_no;
    bool                       _update;
    TrajectoryDiagnosticData  *_tdata;  /*!< \brief Trajectory data (scatter) */
    Histogram                 *_histo;  /*!< \brief Histogram data */
    Emittance                 *_emit;   /*!< \brief Emittance data */

    XYGraph                   *_scatter;

    XYGraph                   *_ellipse;
    bool                       _ellipse_enable;

    Colormap                  *_colormap;
    std::vector<double>        _zdata;

    XYGraph                   *_profile;

    size_t                     _histogram_n;
    size_t                     _histogram_m;
    interpolation_e            _interpolation;
    double                     _dot_size;
    
    void build_data( void );
    void merge_bbox( double bbox[4], const double bb[4] );
    
public:

    ParticleDiagPlot( Frame *frame, const Geometry *geom, const ParticleDataBase *pdb, 
		      coordinate_axis_e axis, double level, 
		      particle_diag_plot_type_e type,
		      trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy = DIAG_NONE );

    ~ParticleDiagPlot();

    void set_emittance_ellipse( bool enable ) {
	_ellipse_enable = enable;
    }

    bool get_emittance_ellipse( void ) {
	return( _ellipse_enable );
    }

    void set_view( coordinate_axis_e axis, double level ) {	
	_update = true;
	_axis = axis;
	_level = level;
    }

    void get_view( coordinate_axis_e &axis, double &level ) {
	axis = _axis;
	level = _level;
    }

    void set_type( particle_diag_plot_type_e type ) {
	_update = true;
	_type = type;
    }

    particle_diag_plot_type_e get_type( void ) {
	return( _type );
    }

    void set_plot( particle_diag_plot_type_e type,
		   trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy ) {
	_update = true;
	_type = type;
	_diagx = diagx;
	_diagy = diagy;
    }

    void get_plot( particle_diag_plot_type_e &type,
		   trajectory_diagnostic_e &diagx, trajectory_diagnostic_e &diagy ) {
	type = _type;
	diagx = _diagx;
	diagy = _diagy;
    }

    void set_histogram_n( size_t n ) {
	_update = true;
	_histogram_n = n;
    }

    size_t get_histogram_n( void ) {
	return( _histogram_n );
    }

    void set_histogram_m( size_t m ) {
	_update = true;
	_histogram_m = m;
    }

    size_t get_histogram_m( void ) {
	return( _histogram_m );
    }

    void set_colormap_interpolation( interpolation_e interpolation ) {
	_interpolation = interpolation;
	if( _colormap )
	    _colormap->set_interpolation( interpolation );
    }

    interpolation_e get_colormap_interpolation( void ) {
	return( _interpolation );
    }

    void set_dot_size( double size ) {
	_dot_size = size;
	if( _scatter )
	    _scatter->set_point_style( XYGRAPH_POINT_CIRCLE, true, _dot_size );
    }

    double get_dot_size( void ) {
	return( _dot_size );
    }

    /*! \brief Calculate Emittance fit.
     */
    const Emittance &calculate_emittance( void );

    /*! \brief Export plotted data as ASCII.
     */
    void export_data( const std::string &filename );

    /*! \brief Rebuild plot.
     */
    void build_plot( void );
};


#endif


