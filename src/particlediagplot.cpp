/*! \file particlediagplot.cpp
 *  \brief Source code for particlediagplot.cpp
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

#include <limits>

#include "particlediagplot.hpp"
#include "histogram.hpp"


ParticleDiagPlot::ParticleDiagPlot( Frame *frame, const Geometry *geom, const ParticleDataBase *pdb, 
				    coordinate_axis_e axis, double level, 
				    particle_diag_plot_type_e type,
				    trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy )
    : _frame(frame), _geom(geom), _pdb(pdb), _axis(axis), _level(level), 
      _type(type), _diagx(diagx), _diagy(diagy),
      _scatter(NULL), _ellipse(NULL), _ellipse_enable(true), _colormap(NULL), _profile(NULL), 
      _histogram_n(50), _histogram_m(50), _interpolation(INTERPOLATION_CLOSEST), _dot_size(1.0)
{

}


ParticleDiagPlot::~ParticleDiagPlot()
{
    if( _scatter )
	delete _scatter;
    if( _colormap )
	delete _colormap;
    if( _ellipse )
	delete _ellipse;
    if( _profile )
	delete _profile;
}


void ParticleDiagPlot::build_plot( void )
{
    // Set axis labels
    if( _type == PARTICLE_DIAG_PLOT_HISTO1D ) {
	_frame->set_axis_label( PLOT_AXIS_X1, trajectory_diagnostic_string_with_unit[_diagx] );
	_frame->set_axis_label( PLOT_AXIS_Y1, "Intensity (a.u.)" );
    } else {
	_frame->set_axis_label( PLOT_AXIS_X1, trajectory_diagnostic_string_with_unit[_diagx] );
	_frame->set_axis_label( PLOT_AXIS_Y1, trajectory_diagnostic_string_with_unit[_diagy] );
    }

    // Clear old graphs
    _frame->clear_graphs();
    if( _scatter ) {
	delete _scatter;
	_scatter = NULL;
    }
    if( _colormap ) {
	delete _colormap; 
	_colormap = NULL;
    }
    if( _ellipse ) {
	delete _ellipse;
	_ellipse = NULL;
    }
    if( _profile ) {
	delete _profile;
	_profile = NULL;
    }

    // Get diagnostic data
    TrajectoryDiagnosticData tdata;
    std::vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( _diagx );
    if( _type != PARTICLE_DIAG_PLOT_HISTO1D )
	diagnostics.push_back( _diagy );
    diagnostics.push_back( DIAG_CURR );
    _pdb->trajectories_at_plane( tdata, _axis, _level, diagnostics );

    // Do data mirroring (limited to only one mirroring per
    // axis-direction, lower end dominates if both edges have
    // mirroring enabled)
    bool mirror[6];
    _pdb->get_mirror( mirror );
    // Mirror in x-direction
    if( mirror[0] ) {
	tdata.mirror( AXIS_X, _geom->origo(0) );
    } else if( mirror[1] ) {
	tdata.mirror( AXIS_X, _geom->max(0) );
    }
    // Mirror in y-direction
    if( mirror[2] ) {
	tdata.mirror( AXIS_Y, _geom->origo(1) );
    } else if( mirror[3] ) {
	tdata.mirror( AXIS_Y, _geom->max(1) );
    }
    // Mirror in z-direction
    if( mirror[4] ) {
	tdata.mirror( AXIS_Z, _geom->origo(2) );
    } else if( mirror[5] ) {
	tdata.mirror( AXIS_Z, _geom->max(2) );
    }

    if( _type == PARTICLE_DIAG_PLOT_SCATTER ) {

	// Make scatter
	_scatter = new XYGraph( tdata(0).data(), tdata(1).data() );
	_scatter->set_line_style( XYGRAPH_LINE_DISABLE );
	_scatter->set_point_style( XYGRAPH_POINT_CIRCLE, true, _dot_size );
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _scatter );

    } else if( _type == PARTICLE_DIAG_PLOT_HISTO2D ) {

	// Make colormap
	Histogram2D histogram( _histogram_n, _histogram_m, tdata(0).data(), tdata(1).data(), tdata(2).data() );
	double range[4];
	histogram.get_range( range );
	double zmin, zmax;
	histogram.get_bin_range( zmin, zmax );

	Palette palette;
	palette.clear();
	if( zmin == 0.0 && zmax >= 0.0 ) {
	    // Palette for positive beam
	    palette.push_back( Color(1,1,1), 0 );
	    palette.push_back( Color(1,1,0), 1 );
	    palette.push_back( Color(1,0,0), 2 );
	    palette.push_back( Color(0,0,0), 3 );
	} else if( zmax == 0.0 && zmin <= 0.0 ) {
	    // Palette for negative beam
	    palette.push_back( Color(1,1,1), 3 );
	    palette.push_back( Color(1,1,0), 2 );
	    palette.push_back( Color(1,0,0), 1 );
	    palette.push_back( Color(0,0,0), 0 );
	} else {
	    // Palette for positive and negative beam
	    palette.push_back( Color(0,0,0), zmin );
	    palette.push_back( Color(0,0,1), 0.67*zmin );
	    palette.push_back( Color(0,1,1), 0.33*zmin );
	    palette.push_back( Color(1,1,1), 0 );
	    palette.push_back( Color(1,1,0), 0.33*zmax );
	    palette.push_back( Color(1,0,0), 0.67*zmax );
	    palette.push_back( Color(0,0,0), zmax );
	}
	palette.norm();

	if( _geom->geom_mode() == MODE_CYL && _diagx == DIAG_R && _diagy == DIAG_RP) {
	    // Scale emittance plot to have constant area per square.
	    double dr = histogram.nstep();
	    for( size_t i = 0; i < histogram.n(); i++ ) {
		double r = fabs( histogram.icoord( i ) );
		double w = M_PI*((r+0.5*dr)*(r+0.5*dr) - (r-0.5*dr)*(r-0.5*dr));
		if( r == 0.0 )
		    w = M_PI*dr*dr;
		for( size_t j = 0; j < histogram.m(); j++ )
		    histogram(i,j) /= w;
	    }
	}
	_colormap = new Colormap( range, histogram.n(), histogram.m(), histogram.get_data() );
	_colormap->set_palette( palette );
	_colormap->set_interpolation( _interpolation );
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _colormap );

    } else if( _type == PARTICLE_DIAG_PLOT_HISTO1D ) {

	// Make XYGraph profile plot
	Histogram histogram( _histogram_n, tdata(0).data(), tdata(1).data() );
	if( _geom->geom_mode() == MODE_CYL && _diagx == DIAG_R ) {
	    std::cout << "scaling profile histogram\n";
	    // Scale profile plot to have constant area per bin.
	    double dr = histogram.step();
	    for( size_t i = 0; i < histogram.n(); i++ ) {
		double r = fabs( histogram.coord(i) );
		double w = M_PI*((r+0.5*dr)*(r+0.5*dr) - (r-0.5*dr)*(r-0.5*dr));
		if( r == 0.0 )
		    w = M_PI*dr*dr;
		histogram(i) /= w;
	    }
	}
	std::vector<double> xdata;
	xdata.reserve( histogram.n() );
	for( size_t a = 0; a < histogram.n(); a++ )
	    xdata.push_back( histogram.coord(a) );

	_profile = new XYGraph( xdata, histogram.get_data() );
	_profile->set_color( Color(1,0,0) );
	_profile->set_line_style( XYGRAPH_LINE_SOLID );
	_profile->set_point_style( XYGRAPH_POINT_DISABLE );
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _profile );
    }

    if( (_type == PARTICLE_DIAG_PLOT_HISTO2D ||
	 _type == PARTICLE_DIAG_PLOT_SCATTER) &&
	((_diagx == DIAG_X && _diagy == DIAG_XP)  ||
	 (_diagx == DIAG_Y && _diagy == DIAG_YP)  ||
	 (_diagx == DIAG_R && _diagy == DIAG_RP)  ||
	 (_diagx == DIAG_Z && _diagy == DIAG_ZP)) && _ellipse_enable ) {

	// Make emittance calculation
	Emittance em( tdata(0).data(), tdata(1).data(), tdata(2).data() );
	double a = em.rmajor();
	double b = em.rminor();
	double p = em.angle();
	double sinp = sin(p);	
	double cosp = cos(p);
	std::vector<double> xd, yd;
	for( size_t i = 0; i < 100; i++ ) {
	    double t = 2.0*M_PI*i/99.0;
	    xd.push_back( em.xave()  + a*cos(t)*sinp + b*sin(t)*cosp );
	    yd.push_back( em.xpave() + a*cos(t)*cosp - b*sin(t)*sinp );
	}

	// Add ellipse xy graph
	_ellipse = new XYGraph( xd, yd );
	_ellipse->set_color( Color(0,0,0) );
	_ellipse->set_line_style( XYGRAPH_LINE_SOLID );
	_ellipse->set_point_style( XYGRAPH_POINT_DISABLE );
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _ellipse );

	// Add emittance numbers to title
	std::stringstream ss;
	ss << "Emittance plot at " << coordinate_axis_string[_axis] << " = " << _level << " m\n"
	   << "\\alpha  = "   << em.alpha()   << ", "
	   << "\\beta  = "    << em.beta()    << " m/rad, "
	   << "\\gamma  = "   << em.gamma()   << " rad/m, "
	   << "\\epsilon  = " << em.epsilon() << " \\pi \\cdot m\\cdot rad";
	_frame->set_title( ss.str().c_str() );

    } else if( _type == PARTICLE_DIAG_PLOT_HISTO1D && 
	       (_diagx == DIAG_X || _diagx == DIAG_Y ||
		_diagx == DIAG_R || _diagx == DIAG_Z) ) {

	// Make title for profile
	std::stringstream ss;
	ss << "Profile plot at " << coordinate_axis_string[_axis] << " = " << _level << " m";
	_frame->set_title( ss.str().c_str() );

    } else {

	// Remove title
	_frame->set_title( "" );

    }
}


