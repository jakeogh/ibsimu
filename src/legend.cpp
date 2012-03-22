/*! \file legend.cpp
 *  \brief Plot legends
 */

/* Copyright (c) 2005-2009,2011-2012 Taneli Kalvas. All rights reserved.
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

#include "legend.hpp"


// Legend sizes are relative to font size
#define LEGEND_SAMPLE_WIDTH 3.0
#define LEGEND_SAMPLE_HEIGHT 1.0
#define LEGEND_SAMPLE_LINESKIP 1.5
#define LEGEND_SAMPLE_SEPARATION 0.3
#define LEGEND_MARGIN 0.5


#define COLORMAP_LEGEND_WIDTH 1.2


void LegendEntry::plot( cairo_t *cairo, double x, double y )
{
    double scale = _label.get_font_size();
    _graph.plot_sample( cairo, x, y, scale*LEGEND_SAMPLE_WIDTH, scale*LEGEND_SAMPLE_HEIGHT );
    _label.set_location( x + scale*(LEGEND_SAMPLE_WIDTH + LEGEND_SAMPLE_SEPARATION), 
			 y );
    _label.draw( cairo );
}


void LegendEntry::get_size( cairo_t *cairo, double &width, double &height ) 
{
    double bbox[4];
    double scale = _label.get_font_size();
    _label.get_bbox( cairo, bbox );
    width = scale*(LEGEND_SAMPLE_WIDTH + LEGEND_SAMPLE_SEPARATION)
	+ (bbox[2]-bbox[0]);
    height = scale*LEGEND_SAMPLE_HEIGHT;
}


void LegendEntry::set_font_size( double fontsize )
{
    _label.set_font_size( fontsize );
}


/* ********************************************************************* */



MultiEntryLegend::MultiEntryLegend()
    : _fontsize(12.0)
{

}


void MultiEntryLegend::add_entry( LegendEntry *entry )
{
    if( entry )
	entry->set_font_size( _fontsize );
    _entry.push_back( entry );
}


void MultiEntryLegend::clear_entries( void )
{
    _entry.clear();
}


void MultiEntryLegend::plot( cairo_t *cairo, double x, double y )
{
    uint32_t k = 0;
    for( int32_t i = _entry.size()-1; i >= 0; i-- ) {
	if( _entry[i] ) {
	    _entry[i]->plot( cairo, 
			     x + _fontsize*LEGEND_MARGIN, 
			     y - _fontsize*(LEGEND_SAMPLE_LINESKIP*k+LEGEND_MARGIN) );
	    k++;
	}
    }
}


void MultiEntryLegend::get_size( cairo_t *cairo, double &width, double &height )
{
    width = height = 0.0;
    for( uint32_t i = 0; i < _entry.size(); i++ ) {
	double w, h;
	if( _entry[i] ) {
	    _entry[i]->get_size( cairo, w, h );
	    if( w > width )
		width = w;
	    if( i == _entry.size()-1 )
		height += _fontsize*LEGEND_SAMPLE_HEIGHT;
	    else
		height += _fontsize*LEGEND_SAMPLE_LINESKIP;
	}
    }

    height += _fontsize*2.0*LEGEND_MARGIN;
    width += _fontsize*2.0*LEGEND_MARGIN;
}


void MultiEntryLegend::set_font_size( double fontsize )
{
    _fontsize = fontsize;
    for( uint32_t i = 0; i < _entry.size(); i++ ) {
	if( _entry[i] )
	    _entry[i]->set_font_size( fontsize );
    }
}


/* ********************************************************************* */



ColormapLegend::ColormapLegend( Colormap &colormap ) 
    : _height(0.0), _fontsize(12.0), _color((Color(0,0,0))),
      _ticlen_in(5.0), _ticlen_out(5.0),
      _ticspace(5.0), _colormap(colormap) 
{

}


void ColormapLegend::build_legend( double x, double y )
{
    _colormap.get_zrange( _range[0], _range[1] );

    // Make 4 tics
    _tic.clear();
    size_t N = 4;
    for( size_t a = 0; a < N; a++ ) {

	// zval span should depend on zscale settings
	double zval = _range[0] + (_range[1]-_range[0])*a/(N-1.0);
	double xx = COLORMAP_LEGEND_WIDTH*_fontsize + _ticlen_out*2 + _ticspace;
	double yy = y+_height*a/(N-1.0);

	char str[128];
	snprintf( str, 128, "%g", zval );

	_tic.push_back( Tic(yy,str) );
	_tic[a]._label.set_location( xx, yy );
	_tic[a]._label.set_alignment( 0.0, 0.5, true );
    }
}


void ColormapLegend::plot( cairo_t *cairo, double x, double y )
{
    build_legend( x, y );

    // Set cairo parameters
    cairo_save( cairo );
    cairo_set_source_rgba( cairo, _color[0], _color[1], _color[2], _color[3] );
    cairo_set_line_width( cairo, 1.0 );

    // Draw box and tics
    cairo_rectangle( cairo, x+_ticlen_out, y, 
		     COLORMAP_LEGEND_WIDTH*_fontsize, _height );
    for( size_t a = 0; a < _tic.size(); a++ ) {
	cairo_move_to( cairo, x,_tic[a]._loc );
	cairo_line_to( cairo, x+_ticlen_out+_ticlen_in,_tic[a]._loc );
    }
    cairo_stroke( cairo );
    for( size_t a = 0; a < _tic.size(); a++ ) {
	_tic[a]._label.draw( cairo );
    }

    cairo_restore( cairo );
}


void ColormapLegend::get_size( cairo_t *cairo, double &width, double &height )
{
    build_legend( 0.0, 0.0 );

    double maxwidth = COLORMAP_LEGEND_WIDTH*_fontsize + _ticlen_out*2;
    for( size_t a = 0; a < _tic.size(); a++ ) {
	double bbox[4];
	_tic[a]._label.get_bbox( cairo, bbox );
	if( bbox[2] > maxwidth )
	    maxwidth = bbox[2];
    }

    width = maxwidth;
    height = _height;
}


void ColormapLegend::set_font_size( double fontsize )
{
    _fontsize = fontsize;
    _ticlen_in  = 5.0*fontsize/12.0;
    _ticlen_out = 5.0*fontsize/12.0;
    _ticspace   = 5.0*fontsize/12.0;
}


void ColormapLegend::set_height( double height )
{
    _height = height;
}

