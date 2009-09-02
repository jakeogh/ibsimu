#ifndef PARTICLEDIAGPLOT_HPP
#define PARTICLEDIAGPLOT_HPP 1


#include "frame.hpp"
#include "geometry.hpp"
#include "particledatabase.hpp"
#include "types.hpp"

#include "xygraph.hpp"
#include "colormap.hpp"



enum particle_diag_plot_type_e {
    PARTICLE_DIAG_PLOT_SCATTER = 0,
    PARTICLE_DIAG_PLOT_HISTO2D,
    PARTICLE_DIAG_PLOT_HISTO1D
};


/*! \brief Particle diagnostic plot.
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
    
    void merge_bbox( double bbox[4], const double bb[4] );
    
public:

    ParticleDiagPlot( Frame *frame, const Geometry *geom, const ParticleDataBase *pdb, 
		      coordinate_axis_e axis, double level, 
		      particle_diag_plot_type_e type,
		      trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy );

    ~ParticleDiagPlot();

    void set_emittance_ellipse( bool enable ) {
	_ellipse_enable = enable;
    }

    bool get_emittance_ellipse( void ) {
	return( _ellipse_enable );
    }

    void set_view( coordinate_axis_e axis, double level ) {	
	_axis = axis;
	_level = level;
    }

    void get_view( coordinate_axis_e &axis, double &level ) {
	axis = _axis;
	level = _level;
    }

    void set_type( particle_diag_plot_type_e type ) {
	_type = type;
    }

    particle_diag_plot_type_e get_type( void ) {
	return( _type );
    }

    void set_plot( particle_diag_plot_type_e type,
		   trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy ) {
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
	_histogram_n = n;
    }

    size_t get_histogram_n( void ) {
	return( _histogram_n );
    }

    void set_histogram_m( size_t m ) {
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
    }

    double get_dot_size( void ) {
	return( _dot_size );
    }

    /*! \brief Rebuild plot.
     */
    void build_plot( void );
};


#endif
