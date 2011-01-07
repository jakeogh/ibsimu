/*! \file fieldgraph.hpp
 *  \brief %Graph for plotting fields
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
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

#ifndef FIELDGRAPH_HPP
#define FIELDGRAPH_HPP 1


#include <vector>
#include "graph3d.hpp"
#include "scalarfield.hpp"
#include "colormap.hpp"
#include "geometry.hpp"


/*! \brief Class for drawing fields with colormap
 *
 *  Implementation of %Graph3D.
 */
class FieldGraph : public Graph3D {

    const ScalarField      *_scalarfield;     /*!< \brief Scalarfield for plotting. */
    Colormap               *_colormap;        /*!< \brief Colormap for field plot. */

    view_e                  _oview;
    double                  _olevel;

    bool                    _enabled;         /*!< \brief Is mesh plotting enabled */

    void build_scalarfield_plot( void );

public:

    /*! \brief Constructor for plotting ScalarField.
     */
    FieldGraph( const ScalarField *field );

    /*! \brief Destructor.
     */
    virtual ~FieldGraph();

    /*! \brief Enable/disable plot.
     */
    void enable( bool enable );

    /*! \brief Plot drawable with cairo.

     *  Plot drawable using \a cairo and coordinate mapper \a cm. The
     *  visible range of plot is given in array \a range in order
     *  xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] );

    /*! \brief Get bounding box of drawable.
     *
     *  Returns the bounding box of the drawable in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );
};


#endif


















