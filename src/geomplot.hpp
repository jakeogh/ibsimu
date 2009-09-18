/*! \file geomplot.hpp
 *  \brief Header file for geomplot.hpp
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

#ifndef GEOMPLOT_HPP
#define GEOMPLOT_HPP 1


#include "frame.hpp"
#include "geometry.hpp"
#include "scalarfield.hpp"
#include "particledatabase.hpp"

#include "solidgraph.hpp"
#include "eqpotgraph.hpp"
#include "particlegraph.hpp"
#include "meshgraph.hpp"
#include "fieldgraph.hpp"


class GeomPlot {

    Frame                   *_frame;

    const Geometry          *_geom;
    const ScalarField       *_epot;
    const ScalarField       *_scharge;
    const ParticleDataBase  *_pdb;

    SolidGraph              *_solidgraph;
    FieldGraph              *_fieldgraph;
    EqPotGraph              *_eqpotgraph;
    ParticleGraph           *_particlegraph;
    MeshGraph               *_meshgraph;

    view_e                  _view;
    int                     _level;
    int                     _vb[3];

    size_t                  _eqlines_auto;
    std::vector<double>     _eqlines_manual;
    size_t                  _particle_div;
    bool                    _scharge_field;
    bool                    _qm_discretation;
    bool                    _mesh;

    bool                    _cache;

    void reset_graphs( void );

public:

    GeomPlot( Frame *frame, const Geometry *geom );

    ~GeomPlot();

    void disable_cache( void );

    void set_epot( const ScalarField *epot );

    void set_eqlines_manual( const std::vector<double> &pot );

    std::vector<double> get_eqlines_manual( void ) {
	return( _eqlines_manual );
    }

    void set_eqlines_auto( size_t N );

    size_t get_eqlines_auto( void ) {
	return( _eqlines_auto );
    }

    void set_scharge( const ScalarField *scharge );

    void set_scharge_field( bool enable );
    
    bool get_scharge_field( void ) {
	return( _scharge_field );
    }
    
    void set_particle_database( const ParticleDataBase *pdb );

    void set_particle_div( size_t particle_div );

    size_t get_particle_div( void ) {
	return( _particle_div );
    }

    void set_qm_discretation( bool enable );

    bool get_qm_discretation( void ) {
	return( _qm_discretation );
    }

    void set_mesh( bool enable );

    bool get_mesh( void ) {
	return( _mesh );
    }

    void set_view( view_e view, int level );

    view_e get_view( void ) const {
	return( _view );
    }

    int get_level( void ) const {
	return( _level );
    }

    double get_level_si( void ) const {
	return( _geom->origo(_vb[2])+_level*_geom->h() );
    }

    int vb( int i ) const {
	return( _vb[i] );
    }

    void get_vb( int vb[3] ) const {
	vb[0] = _vb[0];
	vb[1] = _vb[1];
	vb[2] = _vb[2];
    }
    
};


#endif
