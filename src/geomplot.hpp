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
