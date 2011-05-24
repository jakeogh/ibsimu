/*! \file graph3d.hpp
 *  \brief Base for three dimensional plottable graphs
 */

/* Copyright (c) 2005-2009, 2011 Taneli Kalvas. All rights reserved.
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

#ifndef GRAPH3D_HPP
#define GRAPH3D_HPP 1


#include "graph.hpp"
#include "error.hpp"


/*! \brief View types.
 */
enum view_e {
    VIEW_XY = 0,
    VIEW_XZ,
    VIEW_YX,
    VIEW_YZ,
    VIEW_ZX,
    VIEW_ZY
};


/*! \brief Abstract base class for geometry slice plots.
 *
 *  Provides functionality to select a slice of geometry along axes at
 *  any mesh level in any direction. The direction is selected with \a
 *  view and level in integer (mesh count) \a level. Provides the
 *  child class with permutated coordinate indexes. The first
 *  coordinate (index \a vb[0]) will be used for the horizontal axis
 *  of the plot, the second coordinate (index \a vb[1]) for the
 *  vertical axis and the third coordinate (index \a vb[2]) for the
 *  level.  In three dimensions all the views are possible. In two
 *  dimensions, only \a VIEW_XY and \a VIEW_YX are possible. Here
 *  there is no separate redial index, Y is used instead.
 */
class Graph3D : public Graph {

protected:

    view_e     _view;    /*!< \brief Geometry view direction. */
    int        _vb[3];   /*!< \brief Coordinate index for first, second and third axes. */
    int        _level;   /*!< \brief Level of slice in mesh units. */

public:

    /*! \brief Constructor.
     */
    Graph3D() {
	_view  = VIEW_XY;
	_vb[0] = 0;
	_vb[1] = 1;
	_vb[2] = 2;
	_level = 0;
    }

    /*! \brief Virtual destructor.
     */
    virtual ~Graph3D() {}

    /*! \brief Plot drawable with cairo.
     *
     *  Plot the drawable using \a cairo and coordinate mapper \a
     *  cm. The visible range of plot is given in array \a range in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] ) = 0;

    /*! \brief Plot sample for legend.
     *
     *  Plot graph sample for legend at cairo coordinates \a x.
     */
    virtual void plot_sample( cairo_t *cairo, double x[2] ) = 0;

    /*! \brief Get bounding box of drawable.
     *
     *  Returns the bounding box of the drawable in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] ) = 0;

    /*! \brief Set the view of 3D drawable.
     *
     *  Sets view direction to \a view and the view level to \a level.
     */
    void set_view( view_e view, int level ) {
	switch( view ) {
	case VIEW_XY:
	    _vb[0] = 0;
	    _vb[1] = 1;
	    _vb[2] = 2;
	    break;
	case VIEW_XZ:
	    _vb[0] = 0;
	    _vb[1] = 2;
	    _vb[2] = 1;
	    break;
	case VIEW_YX:
	    _vb[0] = 1;
	    _vb[1] = 0;
	    _vb[2] = 2;
	    break;
	case VIEW_YZ:
	    _vb[0] = 1;
	    _vb[1] = 2;
	    _vb[2] = 0;
	    break;
	case VIEW_ZX:
	    _vb[0] = 2;
	    _vb[1] = 0;
	    _vb[2] = 1;
	    break;
	case VIEW_ZY:
	    _vb[0] = 2;
	    _vb[1] = 1;
	    _vb[2] = 0;
	    break;
	default:
	    throw( ErrorUnimplemented( ERROR_LOCATION ) );
	    break;
	}
	_view = view;
	_level = level;
    }
};


#endif
