/*! \file fieldgraph.hpp
 *  \brief %Graph for plotting fields
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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
#include "vectorfield.hpp"
#include "types.hpp"


/*! \brief Class for drawing fields with colormap
 *
 *  Implementation of %Graph3D.
 */
class FieldGraph : public Graph3D {

    field_type_e            _field_type;      /*!< \brief Field type used. */
    const Geometry         *_geom;            /*!< \brief Geometry. */
    const ScalarField      *_scalarfield;     /*!< \brief Scalarfield for plotting. */
    const VectorField      *_vectorfield;     /*!< \brief Vectorfield for plotting. */    
    Colormap               *_colormap;        /*!< \brief Colormap for field plot. */

    view_e                  _oview;
    double                  _olevel;

    bool                    _enabled;         /*!< \brief Is plotting enabled */
    bool                    _logscale;        /*!< \brief Logarithmic scaling */

    void build_scalarfield_plot( void );
    void build_vectorfield_plot( void );

public:

    /*! \brief Constructor for plotting ScalarField.
     */
    FieldGraph( const ScalarField *field );

    /*! \brief Constructor for plotting a VectorField \a field in geometry \a geom.
     */
    FieldGraph( const Geometry *geom, const VectorField *field, field_type_e field_type );

    /*! \brief Destructor.
     */
    virtual ~FieldGraph();

    /*! \brief Enable/disable plot.
     */
    void enable( bool enable );

    /*! \brief Set logarithmic scale.
     */
    void set_logscale( bool enable ) {
	_logscale = enable;
    }

    /*! \brief Plot drawable with cairo.

     *  Plot drawable using \a cairo and coordinate mapper \a cm. The
     *  visible range of plot is given in array \a range in order
     *  xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] );

    /*! \brief Plot sample for legend.
     *
     *  Plot graph sample for legend at cairo coordinates \a x.
     */
    virtual void plot_sample( cairo_t *cairo, double x, double y, double width, double height );

    /*! \brief Get bounding box of drawable.
     *
     *  Returns the bounding box of the drawable in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );
};


#endif




















