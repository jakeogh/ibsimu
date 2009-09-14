/*! \file particlegraph.cpp
 *  \brief Particle trajectory plotter.
 */

#include <string.h>
#include <iostream>
#include "particlegraph.hpp"
#include "verbose.hpp"


ParticleGraph::ParticleGraph( const Geometry &g, const ParticleDataBase &pdb, 
			    int particlediv, bool qm_discr )
    : _g(g), _pdb(pdb), _particlediv(particlediv), 
      _coordsize(3), _qm_discr(qm_discr)
{
    // Check input
    if( _particlediv <= 0 )
	_particlediv = 0;

    // Allocate work space
    _coord = new double[2*_coordsize];

    // Add default colors
    _color.push_back( Color( 1.0, 0.2, 0.2 ) );
    _color.push_back( Color( 1.0, 1.0, 0.2 ) );
    _color.push_back( Color( 1.0, 0.2, 1.0 ) );
    _color.push_back( Color( 0.2, 1.0, 1.0 ) );
}


ParticleGraph::~ParticleGraph()
{
    delete [] _coord;
}


void ParticleGraph::set_particle_div( size_t particle_div )
{
    _particlediv = particle_div;
    if( _particlediv <= 0 )
	_particlediv = 0;
}


void ParticleGraph::set_qm_discretation( bool qm_discr )
{
    _qm_discr = qm_discr;
}


/*  Subroutine for drawing linear interpolation line between two
 *  particle trajectory data points.  Here x = (x,y,vx,vy,t) is a
 *  point in the trajectory. Parameter first indicates first point of
 *  the trajectory.
 */
void ParticleGraph::draw_linear( const Coordmapper *cm, LineClip &lc, 
				double x[5], bool first ) const
{
    double xout[2];
    cm->transform( xout, x );
    if( first )
	lc.move_to( xout[0], xout[1] );
    else
	lc.line_to( xout[0], xout[1] );
}


void ParticleGraph::get_point( const Coordmapper *cm, double *coord, double s, 
			      double Ax, double Bx, double Cx, double Dx, 
			      double Ay, double By, double Cy, double Dy ) const
{
    double x[2] = { ((Ax*s + Bx)*s + Cx)*s + Dx, 
		    ((Ay*s + By)*s + Cy)*s + Dy };
    cm->transform( coord, x );
}


/*  Subroutine for drawing curved line (consisting of several straight
 *  lines) between two particle trajectory data points.  Here x =
 *  (x,y,vx,vy,t) is a point in the trajectory. Parameter first
 *  indicates first point of the trajectory.
 */
void ParticleGraph::draw_curve( const Coordmapper *cm, LineClip &lc, 
			       double x[5], bool first )
{
    if( first ) {
	// Save first point
	memcpy( _ox, x, 5*sizeof(double) );
	return;
    }

    // Calculate polynomial coefficients
    double dt = x[4]-_ox[4];

    double Ax = (x[2]+_ox[2])*dt + 2.0*(_ox[0]-x[0]);
    double Bx = 3.0*(x[0]-_ox[0]) - (2.0*_ox[2]+x[2])*dt;
    double Cx = _ox[2]*dt;
    double Dx = _ox[0];

    double Ay = (x[3]+_ox[3])*dt + 2.0*(_ox[1]-x[1]);
    double By = 3.0*(x[1]-_ox[1]) - (2.0*_ox[3]+x[3])*dt;
    double Cy = _ox[3]*dt;
    double Dy = _ox[1];
    
    // Fill database with three points
    size_t coord_N = 3;
    cm->transform( &_coord[0], &_ox[0] );
    get_point( cm, &_coord[2], 0.5, Ax, Bx, Cx, Dx, Ay, By, Cy, Dy );
    cm->transform( &_coord[4], &x[0] );
    
    // Bisect curve until desired accuracy is achieved
    size_t a;
    double maxerr, err, xt, yt;
    while( 1 ) {

	//for( a = 0; a < coord_N; a++ ) {
	//    std::cout << std::setw(8) << _coord[2*a+0]
	//	      << std::setw(8) << _coord[2*a+1] << "\n";
	//}

	// Calculate maximum error between the most recent points and
	// the linear interpolation of points from the last round.
	maxerr = 0.0;
	for( a = 1; a < coord_N; a += 2 ) {
	    xt = 0.5*(_coord[2*(a-1)] + _coord[2*(a+1)]);
	    yt = 0.5*(_coord[2*(a-1)+1] + _coord[2*(a+1)+1]);
	    xt -= _coord[2*a];
	    yt -= _coord[2*a+1];
	    err = sqrt( xt*xt + yt*yt );
	    if( err > maxerr )
		maxerr = err;
	}
	
	//std::cout << "maxerr = " << std::setw(8) << maxerr << "\n";

	// If maximum error is less than 2.5 pixels, it is good enough
	if( maxerr < 2.5 )
	    break;

	// Add more points
	coord_N += coord_N-1;

	// Allocate more space if necessary
	if( _coordsize < coord_N ) {
	    double *ct = new double[2*coord_N];
	    memcpy( ct, _coord, _coordsize*2*sizeof(double) );
	    delete [] _coord;
	    _coordsize = coord_N;
	    _coord = ct;
	}

	// Fill coordinates
	for( a = coord_N-1; a > 0; a-- ) {
	    
	    // Even
	    _coord[2*a]   = _coord[a];
	    _coord[2*a+1] = _coord[a+1];
	    a--;
      
	    // Odd
	    get_point( cm, &_coord[2*a], (double)a/(coord_N-1), 
		       Ax, Bx, Cx, Dx, Ay, By, Cy, Dy );
	}
    }
    
    // Draw lines
    lc.move_to( _coord[0], _coord[1] );
    for( a = 1; a < coord_N; a++ )
	lc.line_to( _coord[2*a+0], _coord[2*a+1] );

    memcpy( _ox, x, 5*sizeof(double) );
}










void ParticleGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    if( verbose_output )
	std::cout << "  Plotting particles\n";

    // No plotting
    if( _particlediv == 0 )
	return;

    // Q/M discriminator set
    std::vector<double> qm_set;
    if( !_qm_discr )
	cairo_set_source_rgb( cairo, _color[0][0], _color[0][1], _color[0][2] );
    cairo_set_line_width( cairo, 1.0 );

    // Set clipping ranges
    double clip[4];
    cm->transform( &clip[0], &range[0] );
    cm->transform( &clip[2], &range[2] );
    LineClip lc( cairo );
    lc.set( clip[0], clip[1], clip[2], clip[3] );

    // Loop through all particles
    for( size_t a = 0; a < _pdb.size(); a += _particlediv ) {

	// No plotting if one or less trajectory points
	if( _pdb.traj_size( a ) <= 1 )
	    continue;

	// Select color for particle
	if( _qm_discr ) {
	    size_t c;
	    const ParticleBase &p = _pdb.particle(a);
	    for( c = 0; c < qm_set.size(); c++ ) {
		if( p.qm() == qm_set[c] )
		    break;
	    }
	    if( c == qm_set.size() )
		qm_set.push_back( p.qm() ); // New q/m

	    // Set color
	    size_t s = c%_color.size();
	    cairo_set_source_rgb( cairo, _color[s][0], _color[s][1], _color[s][2] );
	}

	// Loop through all particle trajectory points
	for( size_t b = 0; b < _pdb.traj_size( a ); b++ ) {

	    double t;
	    Vec3D loc, vel;
	    _pdb.trajectory_point( t, loc, vel, a, b );
	    double x[5] = { loc(_vb[0]), loc(_vb[1]), 
			    vel(_vb[0]), vel(_vb[1]), t };
	    //draw_linear( cm, lc, x, b == 0 );
	    draw_curve( cm, lc, x, b == 0 );
	}
	cairo_stroke( cairo );
    }
}


void ParticleGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _g.origo( _vb[0] );
    bbox[1] = _g.origo( _vb[1] );
    bbox[2] = _g.max( _vb[0] );
    bbox[3] = _g.max( _vb[1] );
}


void ParticleGraph::add_color( const Color &color )
{
    _color.push_back( color );
}


void ParticleGraph::clear_colors( void )
{
    _color.clear();
}













