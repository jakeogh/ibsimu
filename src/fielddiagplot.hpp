/*! \file fielddiagplot.hpp
 *  \brief Header file for fielddiagplot.hpp
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

#ifndef FIELDDIAGPLOT_HPP
#define FIELDDIAGPLOT_HPP 1


#include "frame.hpp"
#include "xygraph.hpp"
#include "vec3d.hpp"
#include "geometry.hpp"
#include "scalarfield.hpp"
#include "efield.hpp"
#include "vectorfield.hpp"


enum field_diag_type_e {
    FIELDD_DIAG_NONE = 0,
    FIELDD_DIAG_EPOT,
    FIELDD_DIAG_EFIELD,
    FIELDD_DIAG_EFIELD_X,
    FIELDD_DIAG_EFIELD_Y,
    FIELDD_DIAG_EFIELD_Z,
    FIELDD_DIAG_SCHARGE,
    FIELDD_DIAG_BFIELD,
    FIELDD_DIAG_BFIELD_X,
    FIELDD_DIAG_BFIELD_Y,
    FIELDD_DIAG_BFIELD_Z
};


enum field_loc_type_e {
    FIELDD_LOC_NONE = 0,
    FIELDD_LOC_X,
    FIELDD_LOC_Y,
    FIELDD_LOC_Z,
    FIELDD_LOC_DIST
};


/*! \brief Field diagnostics plot.
 *
 *  A class for building xy-plots of various fields in simulations.
 */
class FieldDiagPlot {

    Frame              *_frame;

    const Geometry     *_geom;
    const ScalarField  *_epot;
    const Efield       *_efield;
    const ScalarField  *_scharge;
    const VectorField  *_bfield;

    size_t              _N;
    Vec3D               _x1;
    Vec3D               _x2;

    field_diag_type_e   _diag[2];
    field_loc_type_e    _loc[2];

    XYGraph            *_graph[2];


    void build_data( std::vector<double> coord[4], 
		     std::vector<double> fielddata[2] ) const;
    std::string diagnostic_label( field_diag_type_e diag ) const;

public:

    /*! \brief Constructor for field diagnostics plot.
     */
    FieldDiagPlot( Frame *frame, const Geometry *geom );

    /*! \brief Destructor for field diagnostics plot.
     */
    ~FieldDiagPlot();

    /*! \brief Add pointer to electric potential.
     */
    void set_epot( const ScalarField *epot ) {
	_epot = epot;
    }

    /*! \brief Add pointer to electric field.
     */
    void set_efield( const Efield *efield ) {
	_efield = efield;
    }

    /*! \brief Add pointer to space charge density map.
     */
    void set_scharge( const ScalarField *scharge ) {
	_scharge = scharge;
    }

    /*! \brief Add pointer to magnetic field.
     */
    void set_bfield( const VectorField *bfield ) {
	_bfield = bfield;
    }

    /*! \brief Set coordinates for field diagnostics.
     *
     *  The fields to be plotted are evaluated at \a N steps from \a
     *  x1 to \a x2. The first point of the plot is exactly at \a x1
     *  and the last at \a x2.
     */
    void set_coordinates( size_t N, const Vec3D &x1, const Vec3D &x2 ) {
	_N = N;
	_x1 = x1;
	_x2 = x2;
    }

    /*! \brief Set field and location plot types
     *
     *  The plot can have two x-axes and two y-axes. The diagnostic
     *  for y1-axis is set by \a diag[0] and the diagnostic for
     *  y2-axis by \a diag[1]. The x1-axis is defined by \a loc[0] and
     *  x2-axis by \a loc[1].
     */
    void set_diagnostic( const field_diag_type_e diag[2], const field_loc_type_e loc[2] ) {
	_diag[0] = diag[0];
	_diag[1] = diag[1];
	_loc[0] = loc[0];
	_loc[1] = loc[1];
    }

    /*! \brief Export plotted data as ASCII.
     */
    void export_data( const std::string &filename ) const;

    /*! \brief Rebuild plot.
     */
    void build_plot( void );

};



#endif
