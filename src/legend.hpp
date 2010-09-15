/*! \file legend.hpp
 *  \brief Header file for legend.hpp
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

#ifndef LEGEND_HPP
#define LEGEND_HPP 1


#include <vector>
#include <string>
#include "graph.hpp"
#include "colormap.hpp"


/*! \brief Class for legend entry.
 *
 *  The legend entries contain a reference to the drawable so that if
 *  the style used in the drawable is changed, the sample in legend is
 *  automatically changed.
 */
class LegendEntry {
    
    Graph     &_drawable;
    std::string   _label;

public:
    
    /*! \brief Contructor for legend entry.
     */
    LegendEntry( Graph &drawable, std::string &label ) 
	: _drawable(drawable), _label(label) {}

    /*! \brief Copy constructor.
     */
    LegendEntry( const LegendEntry &le ) 
	: _drawable(le._drawable), _label(le._label) {}

    /*! \brief Destructor.
     */
    ~LegendEntry() {}

    /*! \brief Assignment operator.
     */
    LegendEntry &operator=( const LegendEntry &le ) {
	_drawable = le._drawable;
	_label = le._label;
	return( *this );
    }

    /*! \brief Plot legend entry.
     */
    void plot( cairo_t *cairo, double x, double y ) const;

    /*! \brief Get size of legend entry.
     */
    void get_size( double &width, double &height ) const;

};


/*! \brief Base class for legend definition.
 *
 *  Legend is an object that contains a key to the plot styles used in
 *  Graphs. The key contains a sample of the plot style used and a
 *  corresponding text label. 
 *
 *  The Colormap legend is a special case because in addition to the
 *  plot style, the plot z-range is shown in the legend.
 */
class Legend {

public:

    /*! \brief Default constructor for legend.
     */
    Legend() {}

    /*! \brief Virtual destructor.
     */
    virtual ~Legend() {}    

    /*! \brief Plot legend at (x,y).
     */
    virtual void plot( cairo_t *cairo, double x, double y ) const = 0;

    /*! \brief Get size of legend.
     */
    virtual void get_size( double &width, double &height ) const = 0;
};


/*! \brief Legend for presenting plot styles.
 */
class MultiEntryLegend {

    std::vector<LegendEntry> _entry;    /*!< \brief Legend entries. */

public:

    /*! \brief Default constructor for legend.
     */
    MultiEntryLegend() {}
    
    /*! \brief Virtual destructor.
     */
    virtual ~MultiEntryLegend() {}    

    /*! \brief Plot legend at (x,y).
     */
    virtual void plot( cairo_t *cairo, double x, double y ) const;

    /*! \brief Get size of legend.
     */
    virtual void get_size( double &width, double &height ) const;

    /*! \brief Add entry to legend.
     */
    void add_entry( const LegendEntry &entry );
};


/*! \brief Legend for presenting colormap key.
 */
class ColormapLegend {

    double      _height;
    Colormap   &_colormap;


public:

    /*! \brief Default constructor for legend.
     */
    ColormapLegend( Colormap &colormap ) : _height(0.0), _colormap(colormap) {}

    /*! \brief Virtual destructor.
     */
    virtual ~ColormapLegend() {}

    /*! \brief Plot legend at (x,y).
     */
    virtual void plot( cairo_t *cairo, double x, double y ) const;

    /*! \brief Get size of legend.
     */
    virtual void get_size( double &width, double &height ) const;

    /*! \brief Set height of legend.
     */
    void set_height( double height );
};


#endif















