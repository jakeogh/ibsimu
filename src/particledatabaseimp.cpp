/*! \file particledatabaseimp.cpp
 *  \brief Particle database implementation
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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


#include <fstream>
#include "particledatabaseimp.hpp"
#include "compmath.hpp"
#include "random.hpp"
#include "constants.hpp"
#include "mat3d.hpp"



/* ******************************************************************************************* *
 * ParticleDataBaseImp                                                                         *
 * ******************************************************************************************* */


ParticleDataBaseImp::ParticleDataBaseImp( ParticleDataBase *pdb )
    : _epsabs(1e-6), _epsrel(1e-6), _polyint(true), _maxsteps(1000), 
      _maxt(1e-3), _save_points(false), _trajdiv(1), _rhosum(0.0), _iteration(0), 
      _relativistic(false), _bsup_cb(NULL), _thand_cb(NULL), _tend_cb(NULL), 
      _pdb(pdb)
{
    for( size_t a = 0; a < 6; a++ )
	_mirror[a] = false;
}


ParticleDataBaseImp::ParticleDataBaseImp( ParticleDataBase *pdb, std::istream &s )
    : _pdb(pdb)
{
    ibsimu.message( 1 ) << "Constructing ParticleDataBase from stream\n";
    ibsimu.inc_indent();

    _epsabs = read_double( s );
    _epsrel = read_double( s );
    _polyint = read_int8( s );
    _maxsteps = read_int32( s );
    _maxt = read_double( s );
    _save_points = read_int8( s );
    _trajdiv = read_int32( s );
    for( uint32_t a = 0; a < 6; a++ )
	_mirror[a] = read_int8( s );
    _rhosum = read_double( s );
    _stat = ParticleStatistics( s );
    _iteration = read_int32( s );
}


ParticleDataBaseImp::ParticleDataBaseImp( const ParticleDataBaseImp &pdb )
{

}


const ParticleDataBaseImp &ParticleDataBaseImp::operator=( const ParticleDataBaseImp &pdb )
{
    return( *this );
}


ParticleDataBaseImp::~ParticleDataBaseImp()
{

}


double ParticleDataBaseImp::energy_to_velocity( double E, double m )
{
    if( E < 5.0e-9*m*SPEED_C2 ) {
	return( sqrt(2.0*E/m) );
    } else {
	double x = E/(m*SPEED_C2)+1.0;
	double beta = sqrt( 1.0 - 1.0/(x*x) );
	return( beta*SPEED_C );
    }
}


void ParticleDataBaseImp::set_accuracy( double epsabs, double epsrel )
{
    _epsabs = epsabs;
    _epsrel = epsrel;
}


void ParticleDataBaseImp::set_bfield_suppression( const CallbackFunctorD_V *bsup_cb )
{
    _bsup_cb = bsup_cb;
}


void ParticleDataBaseImp::set_trajectory_handler_callback( const TrajectoryHandlerCallback *thand_cb )
{
    _thand_cb = thand_cb;
}


void ParticleDataBaseImp::set_trajectory_end_callback( const TrajectoryEndCallback *tend_cb )
{
    _tend_cb = tend_cb;
}


void ParticleDataBaseImp::set_relativistic( bool enable )
{
    _relativistic = enable;
}


void ParticleDataBaseImp::set_polyint( bool polyint )
{
    _polyint = polyint;
}

    
bool ParticleDataBaseImp::get_polyint( void ) const
{
    return( _polyint );
}

    
void ParticleDataBaseImp::set_max_steps( uint32_t maxsteps )
{
    if( maxsteps <= 0 )
        throw( Error( ERROR_LOCATION, "invalid parameter" ) );
    _maxsteps = maxsteps;
}


void ParticleDataBaseImp::set_max_time( double maxt )
{
    if( maxt <= 0.0 )
        throw( Error( ERROR_LOCATION, "invalid parameter" ) );
    _maxt = maxt;
}


void ParticleDataBaseImp::set_save_all_points( bool save_points )
{
    _save_points = save_points;
}


void ParticleDataBaseImp::set_save_trajectories( uint32_t div )
{
    _trajdiv = div;
}


uint32_t ParticleDataBaseImp::get_save_trajectories( void ) const
{
    return( _trajdiv );
}


void ParticleDataBaseImp::set_mirror( const bool mirror[6] )
{
    _mirror[0] = mirror[0];
    _mirror[1] = mirror[1];
    _mirror[2] = mirror[2];
    _mirror[3] = mirror[3];
    if( geom_mode() == MODE_3D ) {
        _mirror[4] = mirror[4];
        _mirror[5] = mirror[5];
    } else {
        _mirror[4] = false;
        _mirror[5] = false;
    }
}


void ParticleDataBaseImp::get_mirror( bool mirror[6] ) const
{
    mirror[0] = _mirror[0];
    mirror[1] = _mirror[1];
    mirror[2] = _mirror[2];
    mirror[3] = _mirror[3];
    mirror[4] = _mirror[4];
    mirror[5] = _mirror[5];
}


int ParticleDataBaseImp::get_iteration_number( void ) const
{
    return( _iteration );
}


double ParticleDataBaseImp::get_rhosum( void ) const
{
    return( _rhosum );
}


void ParticleDataBaseImp::set_rhosum( double rhosum )
{
    _rhosum = rhosum;
}


const ParticleStatistics &ParticleDataBaseImp::get_statistics( void ) const
{
    return( _stat );
}


void ParticleDataBaseImp::save( std::ostream &os ) const
{
    write_double( os, _epsabs );
    write_double( os, _epsrel );
    write_int8( os, _polyint );
    write_int32( os, _maxsteps );
    write_double( os, _maxt );
    write_int8( os, _save_points );
    write_int32( os, _trajdiv );
    for( uint32_t a = 0; a < 6; a++ )
	write_int8( os, _mirror[a] );
    write_double( os, _rhosum );
    _stat.save( os );
    write_int32( os, _iteration );
}


void ParticleDataBaseImp::debug_print( std::ostream &os ) const 
{
    os << "epsabs = "      << _epsabs << "\n";
    os << "epsrel = "      << _epsrel << "\n";
    os << "polyint = "     << _polyint << "\n";
    os << "maxsteps = "    << _maxsteps << "\n";
    os << "maxt = "        << _maxt << "\n";
    os << "save_points = " << _save_points << "\n";
    os << "trajdiv = "     << _trajdiv << "\n";
    os << "mirror = (";
    for( uint32_t a = 0; a < 5; a++ )
	os << _mirror[a] << ", ";
    os << _mirror[5] << ")\n";
    os << "rhosum = "      << _rhosum << "\n";
    os << "iteration = "   << _iteration << "\n";
}


/* ******************************************************************************************* *
 * ParticleDataBase2DImp                                                                       *
 * ******************************************************************************************* */


ParticleDataBase2DImp::ParticleDataBase2DImp( ParticleDataBase *pdb )
    : ParticleDataBasePPImp<ParticleP2D>(pdb)
{
    
}


ParticleDataBase2DImp::ParticleDataBase2DImp( const ParticleDataBase2DImp &pdb )
    : ParticleDataBasePPImp<ParticleP2D>(pdb)
{

}


ParticleDataBase2DImp::ParticleDataBase2DImp( ParticleDataBase *pdb, std::istream &s )
    : ParticleDataBasePPImp<ParticleP2D>(pdb,s)
{
}


ParticleDataBase2DImp::~ParticleDataBase2DImp()
{

}


const ParticleDataBase2DImp &ParticleDataBase2DImp::operator=( const ParticleDataBase2DImp &pdb )
{
    return( *this );
}


void ParticleDataBase2DImp::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
						       double v, double dvp, double dvt, 
						       double x1, double y1, double x2, double y2 )
{
    ibsimu.message( 1 ) << "Defining a 2d beam\n";

    _particles.reserve( _particles.size()+N );

    QRandom qrng( 2 );
    qrng.set_transformation( 0, Gaussian_Transformation() );
    qrng.set_transformation( 1, Gaussian_Transformation() );

    m *= MASS_U;
    q *= CHARGE_E;
    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = J*s/N; // IQ is in units A/m for 2D
    double vt[2], dv[2];
    ParticleP2D x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	qrng.get( vt );
	dv[0] = transverse[0]*dvt*vt[0] + parallel[0]*dvp*vt[1];
	dv[1] = transverse[1]*dvt*vt[0] + parallel[1]*dvp*vt[1];
	x[2] = parallel[0]*v + dv[0];
	x[4] = parallel[1]*v + dv[1];

	add_particle( Particle2D( IQ, q, m, x ) );
	Isum += IQ;
    }

    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A/m\n";
}


void ParticleDataBase2DImp::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						     double E, double Tp, double Tt, 
						     double x1, double y1, double x2, double y2 )
{
    //add_2d_beam_with_velocity( N, J, q, m, energy_to_velocity(E*CHARGE_E,m*MASS_U), 
    //sqrt(Tp*CHARGE_E/(m*MASS_U)), 
    //sqrt(Tt*CHARGE_E/(m*MASS_U)), 
    //x1, y1, x2, y2 );

    ibsimu.message( 1 ) << "Defining a 2d beam\n";

    _particles.reserve( _particles.size()+N );

    m *= MASS_U;
    q *= CHARGE_E;
    E *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    double v = energy_to_velocity(E,m);
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    QRandom qrng( 2 );
    qrng.set_transformation( 0, Gaussian_Transformation() );
    qrng.set_transformation( 1, Gaussian_Transformation() );

    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = J*s/N; // IQ is in units A/m for 2D
    double vt[2];
    ParticleP2D x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	qrng.get( vt );
	double pveld = dvp*vt[1];
	double pvel = sqrt( 2.0*E/m + comp_sign(pveld)*pveld*pveld );
	x[2] = transverse[0]*dvt*vt[0] + parallel[0]*pvel;
	x[4] = transverse[1]*dvt*vt[0] + parallel[1]*pvel;

	add_particle( Particle2D( IQ, q, m, x ) );
	Isum += IQ;
    }

    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A/m\n";
}


void ParticleDataBase2DImp::add_2d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
							   double a, double b, double e,
							   double Ex, double x0, double y0 )
{
    ibsimu.message( 1 ) << "Defining a 2d beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    QRandom qrng( 2 );

    e *= 4.0; // KV-distributed emittance limited inside 4*e_rms ellipse

    m *= MASS_U;
    q *= CHARGE_E;

    double g = (1.0 + a*a)/b;
    double ymax = sqrt( b*e );
    double ypmax = sqrt( g*e );
    double rn[2];
    double IQ = I/N;

    ParticleP2D x;
    x[0] = 0.0;
    x[1] = x0;
    x[2] = sqrt(2.0*Ex*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	qrng.get( rn );

	// Randomize (y,y')
	double y = -ymax + 2.0*ymax*rn[0];
	double yp = -ypmax + 2.0*ypmax*rn[1];

	// Check if inside ellipse
	double yp1, yp2;
	int nroots = solve_quadratic( b, 2.0*a*y, g*y*y-e, &yp1, &yp2 );
	if( nroots != 2 || yp < yp1 || yp > yp2 )
	    continue;
	
	// Set up particle
	x[3] = y0 + y;
	x[4] = x[2]*yp;

	add_particle( Particle2D( IQ, q, m, x ) );
	n++;
    }
}


void ParticleDataBase2DImp::add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
								 double a, double b, double e,
								 double Ex, double x0, double y0 )
{
    ibsimu.message( 1 ) << "Defining a 2d beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    QRandom qrng( 2 );
    qrng.set_transformation( 0, Gaussian_Transformation() );
    qrng.set_transformation( 1, Gaussian_Transformation() );

    m *= MASS_U;
    q *= CHARGE_E;

    double g = (1.0 + a*a)/b;
    double h = 0.5*(b+g);
    double rmaj = sqrt(0.5*e)*(sqrt(h+1.0)+sqrt(h-1.0));
    double rmin = sqrt(0.5*e)*(sqrt(h+1.0)-sqrt(h-1.0));
    double theta = 0.5*atan2( -2.0*a, b-g );
    double w[2], rn[2];
    double IQ = I/N;

    ParticleP2D x;
    x[0] = 0.0;
    x[1] = x0;
    x[2] = sqrt(2.0*Ex*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	// Randomize point from gaussian distribution
	qrng.get( rn );
	w[0] = rmaj*rn[0];
	w[1] = rmin*rn[1];

	// Rotate to correct angle
	double y  = w[0]*cos(theta) - w[1]*sin(theta);
	double yp = w[0]*sin(theta) + w[1]*cos(theta);

	// Set up particle
	x[3] = y0 + y;
	x[4] = x[2]*yp;

	add_particle( Particle2D( IQ, q, m, x ) );
	n++;
    }
}


void ParticleDataBase2DImp::add_tdens_from_segment( MeshScalarField &tdens, double IQ,
						    ParticleP2D &x1, ParticleP2D &x2 ) const
{
    double dx = 0.0;
    double x[2];
    double t[2];
    int32_t i[2];

    for( size_t a = 0; a < 2; a++ ) {
	double xx = ( x2[2*a+1] - x1[2*a+1] );
	dx  += xx*xx;
	x[a] = 0.5*( x1[2*a+1] + x2[2*a+1] );
	i[a] = (int32_t)floor( ( x[a]-tdens.origo(a) ) * tdens.div_h() );
	t[a] = ( x[a]-(i[a]*tdens.h()+tdens.origo(a)) ) * tdens.div_h();
	
	if( i[a] < 0 || i[a] >= (int32_t)tdens.size(a)-1 )
	    continue;
    }
    dx = sqrt( dx );

    double J = IQ*dx; // A
    int p = tdens.size(0)*i[1] + i[0];
    tdens( p )                 += (1.0-t[0])*(1.0-t[1])*J;
    tdens( p+tdens.size(0) )   += (1.0-t[0])*t[1]*J;
    tdens( p+1 )               += t[0]*(1.0-t[1])*J;
    tdens( p+1+tdens.size(0) ) += t[0]*t[1]*J;
}


void ParticleDataBase2DImp::build_trajectory_density_field( MeshScalarField &tdens ) const
{
    tdens.clear();

    // Go through all particle trajectories
    std::vector< ColData<ParticleP2D> > coldata;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticleP2D x1 = _particles[a]->traj(0);
	ParticleP2D xlast = x1;
	double IQ = _particles[a]->IQ();
	for( size_t b = 1; b < N; b++ ) {
	    ParticleP2D x2 = _particles[a]->traj(b);

	    // Mesh collision points from x1 to x2
	    ColData<ParticleP2D>::build_coldata_linear( coldata, tdens, x1, x2 );
	    if( coldata.size() == 0 )
		continue;

	    // Process trajectory mesh collision points
	    for( size_t c = 0; c < coldata.size(); c++ ) {

		// Add contribution from segment from xlast to coldata[c]
		add_tdens_from_segment( tdens, IQ, xlast, coldata[c]._x );	    
		xlast = coldata[c]._x;
	    }

	    // Next trajectory line
	    x1 = x2;
	}
    }

    // Normalize to trajectory density, A/m2
    tdens /= (tdens.h()*tdens.h());

    // Fix boundaries, which only get half of the contribution they should
    // This should depend on mirroring properties!
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	tdens( i, 0 ) *= 2.0;
	tdens( i, tdens.size(1)-1 ) *= 2.0;
    }
    for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	tdens( 0, j ) *= 2.0;
	tdens( tdens.size(0)-1, j ) *= 2.0;
    }
}


void ParticleDataBase2DImp::save( std::ostream &os ) const
{
    ParticleDataBasePPImp<ParticleP2D>::save( os );
}


void ParticleDataBase2DImp::debug_print( std::ostream &os ) const 
{
    os << "**ParticleDataBase2D\n";
    ParticleDataBasePPImp<ParticleP2D>::debug_print( os );
}


/* ******************************************************************************************* *
 * ParticleDataBaseCylImp                                                                      *
 * ******************************************************************************************* */


ParticleDataBaseCylImp::ParticleDataBaseCylImp( ParticleDataBase *pdb )
    : ParticleDataBasePPImp<ParticlePCyl>(pdb)
{
    
}


ParticleDataBaseCylImp::ParticleDataBaseCylImp( const ParticleDataBaseCylImp &pdb )
    : ParticleDataBasePPImp<ParticlePCyl>(pdb)
{

}


ParticleDataBaseCylImp::ParticleDataBaseCylImp( ParticleDataBase *pdb, std::istream &s )
    : ParticleDataBasePPImp<ParticlePCyl>(pdb,s)
{
}


ParticleDataBaseCylImp::~ParticleDataBaseCylImp()
{

}


const ParticleDataBaseCylImp &ParticleDataBaseCylImp::operator=( const ParticleDataBaseCylImp &pdb )
{
    return( *this );
}


void ParticleDataBaseCylImp::add_tdens_from_segment( MeshScalarField &tdens, double IQ,
						     ParticlePCyl &x1, ParticlePCyl &x2 ) const
{
    double dx = 0.0;
    double x[2];
    double t[2];
    int i[2];

    // x-direction
    double xx = ( x2[1] - x1[1] );
    xx = xx*xx;
    dx += xx;
    x[0] = 0.5*( x1[1] + x2[1] );
    i[0] = (int)floor( ( x[0]-tdens.origo(0) ) * tdens.div_h() );
    t[0] = ( x[0]-(i[0]*tdens.h()+tdens.origo(0)) ) * tdens.div_h();

    // r-direction
    xx = ( x2[3] - x1[3] );
    xx = xx*xx;
    dx += xx;
    x[1] = 0.5*( x1[3] + x2[3] );
    i[1] = (int)floor( ( x[1]-tdens.origo(1) ) * tdens.div_h() );
    double rj1 = i[1]*tdens.h()+tdens.origo(1);
    double rj2 = rj1+tdens.h();
    rj1 = rj1*rj1;
    rj2 = rj2*rj2;
    t[1] = (x[1]*x[1]-rj1) / (rj2-rj1);
    
    dx = sqrt( dx );

    for( size_t a = 0; a < 2; a++ ) {
	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= (int32_t)tdens.size(a)-1 ) {
	    i[a] = tdens.size(a)-2;
	    t[a] = 1.0;
	}
    }

    double J = IQ*dx; // A*m
    int p = tdens.size(0)*i[1] + i[0];
    tdens( p )                 += (1.0-t[0])*(1.0-t[1])*J;
    tdens( p+tdens.size(0) )   += (1.0-t[0])*t[1]*J;
    tdens( p+1 )               += t[0]*(1.0-t[1])*J;
    tdens( p+1+tdens.size(0) ) += t[0]*t[1]*J;
}


void ParticleDataBaseCylImp::build_trajectory_density_field( MeshScalarField &tdens ) const
{
    tdens.clear();

    // Go through all particle trajectories
    std::vector< ColData<ParticlePCyl> > coldata;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticlePCyl x1 = _particles[a]->traj(0);
	ParticlePCyl xlast = x1;
	double IQ = _particles[a]->IQ();
	for( size_t b = 1; b < N; b++ ) {
	    ParticlePCyl x2 = _particles[a]->traj(b);

	    // Mesh collision points from x1 to x2
	    ColData<ParticlePCyl>::build_coldata_linear( coldata, tdens, x1, x2 );
	    if( coldata.size() == 0 )
		continue;

	    // Process trajectory mesh collision points
	    for( size_t c = 0; c < coldata.size(); c++ ) {

		// Add contribution from segment from xlast to coldata[c]
		add_tdens_from_segment( tdens, IQ, xlast, coldata[c]._x );	    
		xlast = coldata[c]._x;
	    }

	    // Next trajectory line
	    x1 = x2;
	}
    }    

    // Normalize to trajectory density, A/m2
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	    if( j == 0 ) {
		double rj2 = tdens.h()+tdens.origo(1);
		tdens( i, j ) /= (M_PI*tdens.h()*(rj2*rj2));
	    } else {
		double rj1 = (j-0.5)*tdens.h()+tdens.origo(1);
		double rj2 = (j+0.5)*tdens.h()+tdens.origo(1);
		tdens( i, j ) /= (M_PI*tdens.h()*(rj2*rj2-rj1*rj1));
	    }
	}
    }

    // Fix boundaries, which only get half of the contribution they should
    // This should depend on mirroring properties!
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	tdens( i, 0 ) *= 2.0;
	tdens( i, tdens.size(1)-1 ) *= 2.0;
    }
    for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	tdens( 0, j ) *= 2.0;
	tdens( tdens.size(0)-1, j ) *= 2.0;
    }
}

void ParticleDataBaseCylImp::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
							double v, double dvp, double dvt, 
							double x1, double y1, double x2, double y2 )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // 0: temperature parallel to defining line (transverse to beam)
    // 1: temperature perpendicular to line (parallel to beam)
    // 2: temperature in angular direction (skew velocity)
    QRandom qrng( 3 ); 
    qrng.set_transformation( 0, Gaussian_Transformation() );
    qrng.set_transformation( 1, Gaussian_Transformation() );
    qrng.set_transformation( 2, Gaussian_Transformation() ); 

    m *= MASS_U;
    q *= CHARGE_E;
    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = 2.0*M_PI*J*s/N;
    double vt[3], dv[3];
    ParticlePCyl x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	qrng.get( vt );
	dv[0] = transverse[0]*dvt*vt[0] + parallel[0]*dvp*vt[1];
	dv[1] = transverse[1]*dvt*vt[0] + parallel[1]*dvp*vt[1];
	dv[2] = dvt*vt[2];
	x[2] = parallel[0]*fabs(v) + dv[0];
	x[4] = parallel[1]*fabs(v) + dv[1];
	if( x[3] == 0.0 )
	    x[5] = 0.0;
	else
	    x[5] = dv[2]/x[3];

	add_particle( ParticleCyl( IQ*x[3], q, m, x ) );
	Isum += IQ*x[3];
    }

    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBaseCylImp::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						      double E, double Tp, double Tt, 
						      double x1, double y1, double x2, double y2 )
{
    add_2d_beam_with_velocity( N, J, q, m, energy_to_velocity(E*CHARGE_E,m*MASS_U), 
			       sqrt(Tp*CHARGE_E/(m*MASS_U)), 
			       sqrt(Tt*CHARGE_E/(m*MASS_U)), 
			       x1, y1, x2, y2 );
}


void ParticleDataBaseCylImp::add_2d_full_gaussian_beam( uint32_t N, double I, double q, double m,
							double Ex, double Tp, double Tt, 
							double x0, double dr )
{
    ibsimu.message( 1 ) << "Defining a cylindrical full gaussian beam.\n";

    _particles.reserve( _particles.size()+N );

    m *= MASS_U;
    q *= CHARGE_E;

    QRandom qrng( 5 ); // vx, y, vy, z, vz
    qrng.set_transformation( 0, Gaussian_Transformation() );
    qrng.set_transformation( 1, Gaussian_Transformation() );
    qrng.set_transformation( 2, Gaussian_Transformation() );
    qrng.set_transformation( 3, Gaussian_Transformation() );
    qrng.set_transformation( 4, Gaussian_Transformation() );
    double rn[5];

    double IQ = I/N;
    double Isum = 0.0;

    double vx0 = sqrt(2.0*Ex*CHARGE_E/m);
    double dp = sqrt(Tp*CHARGE_E/m);
    double dt = sqrt(Tt*CHARGE_E/m);

    ParticlePCyl x;
    x[0] = 0.0;
    x[1] = x0;

    uint32_t n = 0;
    while( n < N ) {

	qrng.get( rn );
	if( rn[1] < 0 || rn[3] < 0 )
	    continue;

	double vx = vx0 + dp*rn[0];
	double y = dr*rn[1];
	double vy = dt*rn[2];
	double z = dr*rn[3];
	double vz = dt*rn[4];

	// Convert to cylindrical coordinates
	double r  = sqrt( y*y + z*z );
	if( r == 0.0 )
	    continue;
	double alpha = atan2( z, y );
	double sina = sin(alpha);
	double cosa = cos(alpha);
	double vr = vy*cosa + vz*sina;
	double va = -vy*sina + vz*cosa;

	// Set up particle
	x[2] = vx;
	x[3] = r;
	x[4] = vr;
	x[5] = va/r;

	add_particle( ParticleCyl( IQ, q, m, x ) );
	Isum += IQ;
	n++;
    }

    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBaseCylImp::add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
								  double a, double b, double e,
								  double Ex, double x0 )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    m *= MASS_U;
    q *= CHARGE_E;

    QRandom qrng( 4 );
    qrng.set_transformation( 0, Gaussian_Transformation() );
    qrng.set_transformation( 1, Gaussian_Transformation() );
    qrng.set_transformation( 2, Gaussian_Transformation() );
    qrng.set_transformation( 3, Gaussian_Transformation() );

    double w[4], rn[4];

    double g = (1.0 + a*a)/b;
    double h = 0.5*(b+g);
    double rmaj = sqrt(0.5*e)*(sqrt(h+1.0)+sqrt(h-1.0));
    double rmin = sqrt(0.5*e)*(sqrt(h+1.0)-sqrt(h-1.0));
    double theta = 0.5*atan2( -2.0*a, b-g );

    double IQ = I/N;

    ParticlePCyl x;
    x[0] = 0.0;
    x[1] = x0;
    x[2] = sqrt(2.0*Ex*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	// Randomize point from gaussian distribution
	qrng.get( rn );
	w[0] = rmaj*rn[0];
	w[1] = rmin*rn[1];
	w[2] = rmaj*rn[2];
	w[3] = rmin*rn[3];

	// Rotate to correct angle
	double sint = sin(theta);
	double cost = cos(theta);
	double y  = w[0]*cost - w[1]*sint;
	double yp = w[0]*sint + w[1]*cost;
	double z  = w[2]*cost - w[3]*sint;
	double zp = w[2]*sint + w[3]*cost;

	// Convert to cylindrical coordinates
	double r  = sqrt( y*y + z*z );
	if( r == 0.0 ) // reject center point
	    continue;
	double alpha = atan2( z, y );
	double sina = sin(alpha);
	double cosa = cos(alpha);
	double rp = yp*cosa + zp*sina;
	double ap = -yp*sina + zp*cosa;

	// Set up particle
	x[3] = r;
	x[4] = x[2]*rp;
	x[5] = x[2]*ap/r;

	add_particle( ParticleCyl( IQ, q, m, x ) );
	n++;
    }
}


void ParticleDataBaseCylImp::save( std::ostream &os ) const
{
    ParticleDataBasePPImp<ParticlePCyl>::save( os );
}


void ParticleDataBaseCylImp::debug_print( std::ostream &os ) const 
{
    os << "**ParticleDataBaseCyl\n";
    ParticleDataBasePPImp<ParticlePCyl>::debug_print( os );
}


/* ******************************************************************************************* *
 * ParticleDataBase3DImp                                                                       *
 * ******************************************************************************************* */


ParticleDataBase3DImp::ParticleDataBase3DImp( ParticleDataBase *pdb )
    : ParticleDataBasePPImp<ParticleP3D>(pdb)
{
    
}


ParticleDataBase3DImp::ParticleDataBase3DImp( const ParticleDataBase3DImp &pdb )
    : ParticleDataBasePPImp<ParticleP3D>(pdb)
{

}


ParticleDataBase3DImp::ParticleDataBase3DImp( ParticleDataBase *pdb, std::istream &s )
    : ParticleDataBasePPImp<ParticleP3D>(pdb,s)
{
}


ParticleDataBase3DImp::~ParticleDataBase3DImp()
{

}


const ParticleDataBase3DImp &ParticleDataBase3DImp::operator=( const ParticleDataBase3DImp &pdb )
{
    return( *this );
}


void ParticleDataBase3DImp::add_cylindrical_beam_with_velocity( uint32_t N, double J, double q, double m, 
								double v, double dvp, double dvt, Vec3D c, 
								Vec3D dir1, Vec3D dir2, double r )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // Random number generator for two positions and three velocities (gaussian)
    QRandom qrng( 5 );
    qrng.set_transformation( 2, Gaussian_Transformation() );
    qrng.set_transformation( 3, Gaussian_Transformation() );
    qrng.set_transformation( 4, Gaussian_Transformation() );
    double qx[5];
    double px[6];

    m *= MASS_U;
    q *= CHARGE_E;
    _rhosum += J/v;
    double IQ = M_PI*r*r*J/N;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( dir1[0] != dir1[0] || dir2[0] != dir2[0] || dir3[0] != dir3[0] ) {
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );
    }

    // Prepare particle
    ParticleP3D x;
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    while( a < N ) {

	qrng.get( qx );

	// Calculate in natural (dir1,dir2,dir3) coordinates
	px[0] = -r + 2.0*r*qx[0];
	px[1] = -r + 2.0*r*qx[1];
	px[2] = 0.0;
	px[3] = dvt*qx[2];
	px[4] = dvt*qx[3];
	px[5] = dvp*qx[4] + v;

	// Reject points outside unit radius circle
	if( px[0]*px[0] + px[1]*px[1] > r*r )
	    continue;

	// Map to world coordinates
	x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0];
	x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
	x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
	x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
	x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
	x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

	add_particle( Particle3D( IQ, q, m, x ) );
	Isum += IQ;
	a++;
    }

    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_cylindrical_beam_with_energy( uint32_t N, double J, double q, double m, 
							      double E, double Tp, double Tt, Vec3D c,
							      Vec3D dir1, Vec3D dir2, double r )
{
    add_cylindrical_beam_with_velocity( N, J, q, m, 
					energy_to_velocity(E*CHARGE_E,m*MASS_U), 
					sqrt(Tp*CHARGE_E/(m*MASS_U)),
					sqrt(Tt*CHARGE_E/(m*MASS_U)),
					c, dir1, dir2, r );
}

void ParticleDataBase3DImp::add_rectangular_beam_with_velocity( uint32_t N, double J, double q, double m, 
								double v, double dvp, double dvt, Vec3D c, 
								Vec3D dir1, Vec3D dir2, double size1, double size2 )
{
    ibsimu.message( 1 ) << "Defining a rectangular beam\n";

    _particles.reserve( _particles.size()+N );

    // Random number generator for two positions and three velocities (gaussian)
    QRandom qrng( 5 );
    qrng.set_transformation( 2, Gaussian_Transformation() );
    qrng.set_transformation( 3, Gaussian_Transformation() );
    qrng.set_transformation( 4, Gaussian_Transformation() );
    double qx[5];
    double px[6];

    m *= MASS_U;
    q *= CHARGE_E;
    _rhosum += J/v;
    double IQ = 4.0*size1*size2*J/N;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( dir1[0] != dir1[0] || dir2[0] != dir2[0] || dir3[0] != dir3[0] ) {
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );
    }

    // Prepare particle
    ParticleP3D x;
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    while( a < N ) {

	qrng.get( qx );

	// Calculate in natural (dir1,dir2,dir3) coordinates
	px[0] = size1*(2.0*qx[0]-1.0);
	px[1] = size2*(2.0*qx[1]-1.0);
	px[2] = 0.0;
	px[3] = dvt*qx[2];
	px[4] = dvt*qx[3];
	px[5] = dvp*qx[4] + v;

	// Map to world coordinates
	x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0];
	x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
	x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
	x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
	x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
	x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

	add_particle( Particle3D( IQ, q, m, x ) );
	Isum += IQ;
	a++;
    }

    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_rectangular_beam_with_energy( uint32_t N, double J, double q, double m, 
							      double E, double Tp, double Tt, Vec3D c, 
							      Vec3D dir1, Vec3D dir2, double size1, double size2 )
{
    add_rectangular_beam_with_velocity( N, J, q, m, 
					energy_to_velocity(E*CHARGE_E,m*MASS_U), 
					sqrt(Tp*CHARGE_E/(m*MASS_U)),
					sqrt(Tt*CHARGE_E/(m*MASS_U)),
					c, dir1, dir2, size1, size2 );
}


void ParticleDataBase3DImp::add_3d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
							   double ay, double by, double ey,
							   double az, double bz, double ez,
							   double Ex, double x0, double y0, double z0 )
{
    ibsimu.message( 1 ) << "Defining a 3d beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    ey *= 4.0; // KV-distributed emittance limited inside 4*e_rms ellipse
    ez *= 4.0; // KV-distributed emittance limited inside 4*e_rms ellipse

    m *= MASS_U;
    q *= CHARGE_E;

    QRandom qrng( 4 );
    double rn[4];

    double gy = (1.0 + ay*ay)/by;
    double gz = (1.0 + az*az)/bz;
    
    double ymax = sqrt( by*ey );
    double ypmax = sqrt( gy*ey );
    double zmax = sqrt( bz*ez );
    double zpmax = sqrt( gz*ez );

    double IQ = I/N;

    ParticleP3D x;
    x[0] = 0.0;
    x[1] = x0;
    x[2] = sqrt(2.0*Ex*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	qrng.get( rn );

	// Randomize (y,y',z,z')
	double y = -ymax + 2.0*ymax*rn[0];
	double yp = -ypmax + 2.0*ypmax*rn[1];
	double z = -zmax + 2.0*zmax*rn[2];
	double zp = -zpmax + 2.0*zpmax*rn[3];

	// Check if inside ellipse
	double yp1, yp2;
	int nroots = solve_quadratic( by, 2.0*ay*y, gy*y*y-ey, &yp1, &yp2 );
	if( nroots != 2 || yp < yp1 || yp > yp2 )
	    continue;
	double zp1, zp2;
	nroots = solve_quadratic( bz, 2.0*az*z, gz*z*z-ez, &zp1, &zp2 );
	if( nroots != 2 || zp < zp1 || zp > zp2 )
	    continue;
	
	// Set up particle
	x[3] = y0 + y;
	x[4] = x[2]*yp;
	x[5] = z0 + z;
	x[6] = x[2]*zp;

	add_particle( Particle3D( IQ, q, m, x ) );
	n++;
    }
}


void ParticleDataBase3DImp::add_3d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
								 double E0, 
								 double a1, double b1, double e1,
								 double a2, double b2, double e2,
								 Vec3D c, Vec3D dir1, Vec3D dir2 )
{
    ibsimu.message( 1 ) << "Defining a 3d beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    m *= MASS_U;
    q *= CHARGE_E;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( dir1[0] != dir1[0] || dir2[0] != dir2[0] || dir3[0] != dir3[0] ) {
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );
    }

    QRandom qrng( 4 );
    qrng.set_transformation( 0, Gaussian_Transformation() );
    qrng.set_transformation( 1, Gaussian_Transformation() );
    qrng.set_transformation( 2, Gaussian_Transformation() );
    qrng.set_transformation( 3, Gaussian_Transformation() );
    double w[4], rn[4];

    double g1 = (1.0 + a1*a1)/b1;
    double h1 = 0.5*(b1+g1);
    double rmaj1 = sqrt(0.5*e1)*(sqrt(h1+1.0)+sqrt(h1-1.0));
    double rmin1 = sqrt(0.5*e1)*(sqrt(h1+1.0)-sqrt(h1-1.0));
    double theta1 = 0.5*atan2( -2.0*a1, b1-g1 );

    double g2 = (1.0 + a2*a2)/b2;
    double h2 = 0.5*(b2+g2);
    double rmaj2 = sqrt(0.5*e2)*(sqrt(h2+1.0)+sqrt(h2-1.0));
    double rmin2 = sqrt(0.5*e2)*(sqrt(h2+1.0)-sqrt(h2-1.0));
    double theta2 = 0.5*atan2( -2.0*a2, b2-g2 );

    double IQ = I/N;

    ParticleP3D x, px;
    px[0] = x[0] = 0.0;
    px[5] = 0.0;
    px[6] = sqrt(2.0*E0*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	// Randomize point from gaussian distribution
	qrng.get( rn );
	w[0] = rmaj1*rn[0];
	w[1] = rmin1*rn[1];
	w[2] = rmaj2*rn[2];
	w[3] = rmin2*rn[3];

	// Rotate to correct angle
	px[1] = w[0]*cos(theta1) - w[1]*sin(theta1);
	px[2] = px[6]*(w[0]*sin(theta1) + w[1]*cos(theta1));
	px[3] = w[2]*cos(theta2) - w[3]*sin(theta2);
	px[4] = px[6]*(w[2]*sin(theta2) + w[3]*cos(theta2));

	// Map to world coordinates
	x[1] = dir1[0]*px[1] + dir2[0]*px[3] + dir3[0]*px[5] + c[0];
	x[2] = dir1[0]*px[2] + dir2[0]*px[4] + dir3[0]*px[6];
	x[3] = dir1[1]*px[1] + dir2[1]*px[3] + dir3[1]*px[5] + c[1];
	x[4] = dir1[1]*px[2] + dir2[1]*px[4] + dir3[1]*px[6];
	x[5] = dir1[2]*px[1] + dir2[2]*px[3] + dir3[2]*px[5] + c[2];
	x[6] = dir1[2]*px[2] + dir2[2]*px[4] + dir3[2]*px[6];

	add_particle( Particle3D( IQ, q, m, x ) );
	n++;
    }
}


void ParticleDataBase3DImp::trajectories_at_free_plane( TrajectoryDiagnosticData &tdata, 
							Vec3D c, Vec3D o, Vec3D p,
							const std::vector<trajectory_diagnostic_e> &diagnostics ) const
{
    // Check input
    for( size_t a = 0; a < diagnostics.size(); a++ ) {
	if( diagnostics[a] != DIAG_T && 
	    diagnostics[a] != DIAG_X && diagnostics[a] != DIAG_VX && 
	    diagnostics[a] != DIAG_Y && diagnostics[a] != DIAG_VY && 
	    diagnostics[a] != DIAG_Z && diagnostics[a] != DIAG_VZ && 
	    diagnostics[a] != DIAG_O && diagnostics[a] != DIAG_VO && 
	    diagnostics[a] != DIAG_P && diagnostics[a] != DIAG_VP && 
	    diagnostics[a] != DIAG_Q && diagnostics[a] != DIAG_VQ && 
	    diagnostics[a] != DIAG_OP && diagnostics[a] != DIAG_PP && 
	    diagnostics[a] != DIAG_CURR && diagnostics[a] != DIAG_EK && 
	    diagnostics[a] != DIAG_QM && diagnostics[a] != DIAG_CHARGE && 
	    diagnostics[a] != DIAG_MASS )
	    throw( Error( ERROR_LOCATION, "invalid diagnostic type " + to_string(diagnostics[a]) ) );
    }

    // Forward vector (normal to plane)
    Vec3D q = cross(o,p);
    if( q.norm2() == 0.0 )
	throw( Error( ERROR_LOCATION, "invalid vectors defining plane" ) );

    // Orthogonalize vector p
    p = cross(q,o);

    // Normalize vectors
    o.normalize();
    p.normalize();
    q.normalize();

    ibsimu.message( 1 ) << "Making trajectory diagnostics at plane\n" 
			<< "  c = " << c << "\n"
			<< "  o = " << o << "\n"
			<< "  p = " << p << "\n"
			<< "  q = " << q << "\n";

    // Prepare output vector
    tdata.clear();
    for( size_t a = 0; a < diagnostics.size(); a++ ) {
	tdata.add_data_column( diagnostics[a] );
    }

    // Scan through particle trajectory points
    double Isum = 0.0;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticleP3D x1 = _particles[a]->traj(0);
	for( size_t b = 1; b < N; b++ ) {
	    ParticleP3D x2 = _particles[a]->traj(b);

	    // Solve trajectory crossing point using linear interpolation of position
	    // K[0]: parametric distance from x1 to x2
 	    // K[1]: distance in direction o
	    // K[2]: distance in direction p
	    Vec3D K;
	    try {
		Mat3D m( x2[1]-x1[1], -o[0], -p[0],
			 x2[3]-x1[3], -o[1], -p[1],
			 x2[5]-x1[5], -o[2], -p[2] );
		Mat3D minv = m.inverse();
		Vec3D off( c[0]-x1[1], c[1]-x1[3], c[2]-x1[5] );
		K = minv*off;
	    } catch(...) {
		continue;
	    }
	    
	    if( K[0] >= 0.0 && K[0] < 1.0 ) {
		// Crossing found between x1 and x2, accumulate total current
		Isum += _particles[a]->IQ();

		// Position and velocity in xyz coordinates
		double t = K[0];
		double nt = 1.0-K[0];
		Vec3D pos( nt*x1[1]+t*x2[1], nt*x1[3]+t*x2[3], nt*x1[5]+t*x2[5] );
		Vec3D vel( nt*x1[2]+t*x2[2], nt*x1[4]+t*x2[4], nt*x1[6]+t*x2[6] );

		// Position and velocity in opq coordinates
		Vec3D pos_opq( K[1], K[2], 0.0 );
		Vec3D vel_opq( vel*o, vel*p, vel*q );

		// Fill diagnostic data
		for( size_t c = 0; c < tdata.diag_size(); c++ ) {

		    switch( tdata.diagnostic( c ) ) {
		    case DIAG_T:
			tdata.add_data( c, nt*x1[0]+t*x2[0] );
			break;

		    case DIAG_X:
			tdata.add_data( c, pos[0] );
			break;
		    case DIAG_VX:
			tdata.add_data( c, vel[0] );
			break;
		    case DIAG_Y:
			tdata.add_data( c, pos[1] );
			break;
		    case DIAG_VY:
			tdata.add_data( c, vel[1] );
			break;
		    case DIAG_Z:
			tdata.add_data( c, pos[2] );
			break;
		    case DIAG_VZ:
			tdata.add_data( c, vel[2] );
			break;

		    case DIAG_O:
			tdata.add_data( c, pos_opq[0] );
			break;
		    case DIAG_VO:
			tdata.add_data( c, vel_opq[0] );
			break;
		    case DIAG_P:
			tdata.add_data( c, pos_opq[1] );
			break;
		    case DIAG_VP:
			tdata.add_data( c, vel_opq[1] );
			break;
		    case DIAG_Q:
			tdata.add_data( c, pos_opq[2] );
			break;
		    case DIAG_VQ:
			tdata.add_data( c, vel_opq[2] );
			break;

		    case DIAG_OP:
			tdata.add_data( c, vel_opq[0]/vel_opq[2] );
			break;
		    case DIAG_PP:
			tdata.add_data( c, vel_opq[1]/vel_opq[2] );
			break;

		    case DIAG_CURR:
			tdata.add_data( c, _particles[a]->IQ() );
			break;
		    case DIAG_EK:
			if( vel.ssqr() < 1.0e-4*SPEED_C )
			    tdata.add_data( c, 0.5*_particles[a]->m()*vel.ssqr()/CHARGE_E );
			else {
			    double beta = vel.ssqr()/SPEED_C;
			    double gamma = 1.0 / sqrt( 1.0 - beta*beta );
			    double Ek = _particles[a]->m()*SPEED_C2*( gamma - 1.0 );
			    tdata.add_data( c, Ek/CHARGE_E );
			}
			break;
		    case DIAG_QM:
			tdata.add_data( c, (_particles[a]->q()/CHARGE_E) / (_particles[a]->m()/MASS_U) );
			break;
		    case DIAG_CHARGE:
			tdata.add_data( c, _particles[a]->q()/CHARGE_E );
			break;
		    case DIAG_MASS:
			tdata.add_data( c, _particles[a]->m()/MASS_U );
			break;
		    case DIAG_NO:
			tdata.add_data( c, a );
			break;
		    default:
			throw( ErrorUnimplemented( ERROR_LOCATION ) );
			break;
		    }
		}
	    }

	    // Next trajectory step
	    x1 = x2;
	}
    }

    ibsimu.message( 1 ) << "  number of trajectories = " << tdata.traj_size() << "\n";
    ibsimu.message( 1 ) << "  total current = " << Isum << " A\n";
}


void ParticleDataBase3DImp::export_path_manager_data( std::string filename, 
						      double ref_E, double ref_q, double ref_m, 
						      Vec3D c, Vec3D o, Vec3D p ) const
{
    // Forward vector (normal to plane)
    Vec3D q = cross(o,p);
    if( q.norm2() == 0.0 )
	throw( Error( ERROR_LOCATION, "invalid vectors defining plane" ) );

    // Orthogonalize vector p
    p = cross(q,o);

    // Normalize vectors
    o.normalize();
    p.normalize();
    q.normalize();

    // Calculate reference particle
    ref_m *= MASS_U;
    double ref_v = sqrt(2.0*ref_E*CHARGE_E/ref_m);
    double ref_p = ref_m*ref_v;
    const double m_to_gevc2 = SPEED_C2/CHARGE_E/1.0e9;
    const double p_to_gevc = SPEED_C/CHARGE_E/1.0e9;

    ibsimu.message( 1 ) << "Making trajectory diagnostics at plane for Path Manager output\n" 
			<< "  ref_E = " << ref_E << " eV\n"
			<< "  ref_q = " << ref_q << " e\n"
			<< "  ref_m = " << ref_m/MASS_U << " u\n"
			<< "  c = " << c << "\n"
			<< "  o = " << o << "\n"
			<< "  p = " << p << "\n"
			<< "  q = " << q << "\n";
    ibsimu.message( 1 ) << "  Outputting to \'" << filename << "\'\n";

    std::ofstream ofile( filename.c_str() );
    if( !ofile.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\'" ) );

    // Scan through particle trajectory points and build data
    std::vector<double> data[7]; // x, x', y, y', dp/p, q, m
    bool first = true;
    double IQ_first = 0.0;
    bool IQ_warning = false;
    double Isum = 0.0;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticleP3D x1 = _particles[a]->traj(0);
	for( size_t b = 1; b < N; b++ ) {
	    ParticleP3D x2 = _particles[a]->traj(b);

	    // Solve trajectory crossing point using linear interpolation of position
	    // K[0]: parametric distance from x1 to x2
 	    // K[1]: distance in direction o
	    // K[2]: distance in direction p
	    Vec3D K;
	    try {
		Mat3D m( x2[1]-x1[1], -o[0], -p[0],
			 x2[3]-x1[3], -o[1], -p[1],
			 x2[5]-x1[5], -o[2], -p[2] );
		Mat3D minv = m.inverse();
		Vec3D off( c[0]-x1[1], c[1]-x1[3], c[2]-x1[5] );
		K = minv*off;
	    } catch(...) {
		continue;
	    }
	    
	    if( K[0] >= 0.0 && K[0] < 1.0 ) {
		// Crossing found between x1 and x2, accumulate total current
		Isum += _particles[a]->IQ();

		// Position and velocity in xyz coordinates
		double t = K[0];
		double nt = 1.0-K[0];
		Vec3D pos( nt*x1[1]+t*x2[1], nt*x1[3]+t*x2[3], nt*x1[5]+t*x2[5] );
		Vec3D vel( nt*x1[2]+t*x2[2], nt*x1[4]+t*x2[4], nt*x1[6]+t*x2[6] );

		// Position and velocity in opq coordinates
		Vec3D pos_opq( K[1], K[2], 0.0 );
		Vec3D vel_opq( vel*o, vel*p, vel*q );

		// Check IQ fixed
		if( first ) {
		    IQ_first = _particles[a]->IQ();
		} else if( IQ_first == 0.0 ) {
		    if( !IQ_warning && _particles[a]->IQ() != 0.0 ) {
			ibsimu.message( 1 ) << "  WARNING: Trajectory current is not constant.\n";
			IQ_warning = true;
		    }
		} else {
		    if( !IQ_warning && fabs( ( IQ_first - _particles[a]->IQ() ) / IQ_first ) > 1.0e-6 ) {
			ibsimu.message( 1 ) << "  WARNING: Trajectory current is not constant.\n";
			IQ_warning = true;
		    }
		}

		// Save point
		double pq = _particles[a]->q();
		double pm = _particles[a]->m();

		data[0].push_back( pos_opq[0] );                  // x (m)
		data[1].push_back( vel_opq[0]/vel_opq[2] );       // x' (rad)
		data[2].push_back( pos_opq[1] );                  // y (m)
		data[3].push_back( vel_opq[1]/vel_opq[2] );       // y' (rad)
		data[4].push_back( (pm*vel_opq[2]-ref_p)/ref_p ); // (p_z-p_ref)/p_ref
		data[5].push_back( pq/CHARGE_E );                 // q (in e)
		data[6].push_back( pm*m_to_gevc2 );               // m (in GeV/c2)

		first = false;
	    }

	    // Next trajectory step
	    x1 = x2;
	}
    }

    // Write header
    time_t tt = time(NULL);
    char *timebuf = ctime( &tt );
    ofile << "Beam data from IBSimu, " << timebuf; // contains newline
    ofile << "Total beam: " << Isum << " A\n";
    ofile << std::setw(17) << ref_p*p_to_gevc  << " !REFERENCE MOMENTUM [GeV/c]\n";
    ofile << std::setw(17) << 0.0              << " !REFERENCE PHASE [rad]\n";
    ofile << std::setw(17) << 0.0              << " !FREQUENCY [Hz]\n";
    ofile << std::setw(17) << ref_m*m_to_gevc2 << " !REFERENCE MASS [GeV/c2]\n";
    ofile << std::setw(17) << ref_q            << " !REFERENCE CHARGE STATE\n";
    ofile << std::setw(17) << data[0].size()   << " !NUMBER OF PARTICLES\n";

    for( size_t a = 0; a < data[0].size(); a++ ) {	
	ofile << std::setw(17) << a << " ";       // particle number
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[0][a] << " "; // x
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[1][a] << " "; // x'
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[2][a] << " "; // y
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[3][a] << " "; // y'
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << -M_PI + 2.0*M_PI*rand() / RAND_MAX << " "; // phase
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[4][a] << " "; // dp/p
	ofile << std::setw(17) << 0 << " "; // flag (alive=0)
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[5][a] << " "; // q
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[6][a] << " "; // m
	ofile << "\n";
    }

    ibsimu.message( 1 ) << "  number of trajectories = " << data[0].size() << "\n";
    ibsimu.message( 1 ) << "  total current = " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_tdens_from_segment( MeshScalarField &tdens, double IQ,
						    ParticleP3D &x1, ParticleP3D &x2 ) const
    
{
    double dx = 0.0;
    double x[3];
    double t[3];
    int32_t i[3];

    for( size_t a = 0; a < 3; a++ ) {
	double xx = ( x2[2*a+1] - x1[2*a+1] );
	dx  += xx*xx;
	x[a] = 0.5*( x1[2*a+1] + x2[2*a+1] );
	i[a] = (int)floor( ( x[a]-tdens.origo(a) ) * tdens.div_h() );
	t[a] = ( x[a]-(i[a]*tdens.h()+tdens.origo(a)) ) * tdens.div_h();
	
	if( i[a] < 0 || i[a] >= (int32_t)(tdens.size(a)-1) )
	    continue;
    }
    dx = sqrt( dx );

    double J = IQ*dx; // A*m
    int p = tdens.size(0)*tdens.size(1)*i[2] + tdens.size(0)*i[1] + i[0];
    tdens( p )                 += (1.0-t[0])*(1.0-t[1])*(1.0-t[2])*J;
    tdens( p+tdens.size(0) )   += (1.0-t[0])*t[1]*(1.0-t[2])*J;
    tdens( p+1 )               += t[0]*(1.0-t[1])*(1.0-t[2])*J;
    tdens( p+1+tdens.size(0) ) += t[0]*t[1]*(1.0-t[2])*J;
    
    p += tdens.size(0)*tdens.size(1);
    tdens( p )                 += (1.0-t[0])*(1.0-t[1])*t[2]*J;
    tdens( p+tdens.size(0) )   += (1.0-t[0])*t[1]*t[2]*J;
    tdens( p+1 )               += t[0]*(1.0-t[1])*t[2]*J;
    tdens( p+1+tdens.size(0) ) += t[0]*t[1]*t[2]*J;
}


void ParticleDataBase3DImp::build_trajectory_density_field( MeshScalarField &tdens ) const
{
    tdens.clear();

    // Go through all particle trajectories
    std::vector< ColData<ParticleP3D> > coldata;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticleP3D x1 = _particles[a]->traj(0);
	ParticleP3D xlast = x1;
	double IQ = _particles[a]->IQ();
	for( size_t b = 1; b < N; b++ ) {
	    ParticleP3D x2 = _particles[a]->traj(b);

	    // Mesh collision points from x1 to x2
	    ColData<ParticleP3D>::build_coldata_linear( coldata, tdens, x1, x2 );
	    if( coldata.size() == 0 )
		continue;

	    // Process trajectory mesh collision points
	    for( size_t c = 0; c < coldata.size(); c++ ) {

		// Add contribution from segment from xlast to coldata[c]
		add_tdens_from_segment( tdens, IQ, xlast, coldata[c]._x );	    
		xlast = coldata[c]._x;
	    }

	    // Next trajectory line
	    x1 = x2;
	}
    }    

    // Normalize to trajectory density, A/m2
    tdens /= (tdens.h()*tdens.h()*tdens.h());

    // Fix boundaries, which only get half of the contribution they should
    // This should depend on mirroring properties!
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	    tdens( i, j, 0 ) *= 2.0;
	    tdens( i, j, tdens.size(2)-1 ) *= 2.0;
	}
    }
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	for( uint32_t k = 0; k < tdens.size(2); k++ ) {
	    tdens( i, 0, k ) *= 2.0;
	    tdens( i, tdens.size(1)-1, k ) *= 2.0;
	}
    }	
    for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	for( uint32_t k = 0; k < tdens.size(2); k++ ) {
	    tdens( 0, j, k ) *= 2.0;
	    tdens( tdens.size(0)-1, j, k ) *= 2.0;
	}
    }
}


void ParticleDataBase3DImp::save( std::ostream &os ) const
{
    ParticleDataBasePPImp<ParticleP3D>::save( os );
}


void ParticleDataBase3DImp::debug_print( std::ostream &os ) const 
{
    os << "**ParticleDataBase3D\n";
    ParticleDataBasePPImp<ParticleP3D>::debug_print( os );
}



