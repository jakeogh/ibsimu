/*! \file particledatabase.cpp
 *  \brief %Particle databases
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


#include "particledatabase.hpp"
#include "particledatabaseimp.hpp"



ParticleDataBase::ParticleDataBase()
    : _imp(0)
{
}


ParticleDataBase::ParticleDataBase( const ParticleDataBase &pdb )
    : _imp(0)
{
}


ParticleDataBase::~ParticleDataBase()
{
}


const ParticleDataBase &ParticleDataBase::operator=( const ParticleDataBase &pdb )
{
    return( *this );
}


void ParticleDataBase::set_implementation_pointer( class ParticleDataBaseImp *imp )
{
    _imp = imp;
}


void ParticleDataBase::set_accuracy( double epsabs, double epsrel )
{
    _imp->set_accuracy( epsabs, epsrel );
}


void ParticleDataBase::set_bfield_suppression( const CallbackFunctorD_V *functor )
{
    _imp->set_bfield_suppression( functor );
}


void ParticleDataBase::set_trajectory_handler_callback( const TrajectoryHandlerCallback *trajhand )
{
    _imp->set_trajectory_handler_callback( trajhand );
}


void ParticleDataBase::set_polyint( bool polyint ) 
{
    _imp->set_polyint( polyint );
}


bool ParticleDataBase::get_polyint( void ) const 
{
    return( _imp->get_polyint() );
}


void ParticleDataBase::set_max_steps( uint32_t maxsteps ) 
{
    _imp->set_max_steps( maxsteps );
}

void ParticleDataBase::set_max_time( double maxt ) 
{
    _imp->set_max_time( maxt );
}


void ParticleDataBase::set_save_trajectories( uint32_t div ) 
{
    _imp->set_save_trajectories( div );
}


uint32_t ParticleDataBase::get_save_trajectories( void ) const
{
    return( _imp->get_save_trajectories() );
}


void ParticleDataBase::set_mirror( const bool mirror[6] )
{
    _imp->set_mirror( mirror );
}


void ParticleDataBase::get_mirror( bool mirror[6] ) const 
{
    _imp->get_mirror( mirror );
}


int ParticleDataBase::get_iteration_number( void ) const 
{
    return( _imp->get_iteration_number() );
}


double ParticleDataBase::get_rhosum( void ) const
{
    return( _imp->get_rhosum() );
}


const ParticleStatistics &ParticleDataBase::get_statistics( void ) const
{
    return( _imp->get_statistics() );
}


geom_mode_e ParticleDataBase::geom_mode() const
{
    return( _imp->geom_mode() );
}


size_t ParticleDataBase::size( void ) const
{
    return( _imp->size() );
}


size_t ParticleDataBase::traj_size( uint32_t i ) const
{
    return( _imp->traj_size( i ) );
}


void ParticleDataBase::trajectory_point( double &t, Vec3D &loc, Vec3D &vel, uint32_t i, uint32_t j ) const
{
    _imp->trajectory_point( t, loc, vel, i, j );
}


void ParticleDataBase::trajectories_at_plane( TrajectoryDiagnosticData &tdata, 
					      coordinate_axis_e axis,
					      double val,
					      const std::vector<trajectory_diagnostic_e> &diagnostics ) const
{
    _imp->trajectories_at_plane( tdata, axis, val, diagnostics );
}


void ParticleDataBase::clear( void )
{
    _imp->clear();
}


void ParticleDataBase::clear_trajectories( void )
{
    _imp->clear_trajectories();
}


void ParticleDataBase::reserve( size_t size )
{
    _imp->reserve( size );
}


void ParticleDataBase::build_trajectory_density_field( ScalarField &tdens ) const
{
    _imp->build_trajectory_density_field( tdens );
}


void ParticleDataBase::iterate_trajectories( ScalarField &scharge, const VectorField &efield, 
					     const VectorField &bfield, const Geometry &g )
{
    _imp->iterate_trajectories( scharge, efield, bfield, g );
}


void ParticleDataBase::step_particles( ScalarField &scharge, const VectorField &efield, 
				       const VectorField &bfield, const Geometry &g, double dt )
{
    _imp->step_particles( scharge, efield, bfield, g, dt );
}



/* ******************************************************************************************* *
 * ParticleDataBase2D                                                                          *
 * ******************************************************************************************* */


ParticleDataBase2D::ParticleDataBase2D()
{
    _imp = new ParticleDataBase2DImp;
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase2D::ParticleDataBase2D( const ParticleDataBase2D &pdb )
{
    _imp = new ParticleDataBase2DImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase2D::~ParticleDataBase2D()
{
    delete _imp;
}


const ParticleDataBase2D &ParticleDataBase2D::operator=( const ParticleDataBase2D &pdb )
{
    delete _imp;
    _imp = new ParticleDataBase2DImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );

    return( *this );
}


Particle2D &ParticleDataBase2D::particle( uint32_t i )
{
    return( _imp->particle( i ) );
}


const Particle2D &ParticleDataBase2D::particle( uint32_t i ) const
{
    return( _imp->particle( i ) ) ;
}


const ParticleP2D &ParticleDataBase2D::trajectory_point( uint32_t i, uint32_t j ) const
{
    return( _imp->trajectory_point( i, j ) );
}


void ParticleDataBase2D::add_particle( double IQ, double q, double m, const ParticleP2D &x )
{
    _imp->add_particle( IQ, q, m, x );
}


void ParticleDataBase2D::add_particle( const Particle2D &p )
{
    _imp->add_particle( p );
}


void ParticleDataBase2D::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
						    double v, double dvp, double dvt, 
						    double x1, double y1, double x2, double y2 )
{
    _imp->add_2d_beam_with_velocity( N, J, q, m, v, dvp, dvt, x1, y1, x2, y2 );
}


void ParticleDataBase2D::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						  double E, double Tp, double Tt, 
						  double x1, double y1, double x2, double y2 )
{
    _imp->add_2d_beam_with_energy( N, J, q, m, E, Tp, Tt, x1, y1, x2, y2 );
}


void ParticleDataBase2D::add_2d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
							double a, double b, double e,
							double Ex, double x0, double y0 )
{
    _imp->add_2d_KV_beam_with_emittance( N, I, q, m, a, b, e, Ex, x0, y0 );
}


void ParticleDataBase2D::add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
							      double a, double b, double e,
							      double Ex, double x0, double y0 )
{
    _imp->add_2d_gaussian_beam_with_emittance( N, I, q, m, a, b, e, Ex, x0, y0 );
}


void ParticleDataBase2D::debug_print( std::ostream &os ) const
{
    _imp->debug_print( os );
}


/* ******************************************************************************************* *
 * ParticleDataBaseCyl                                                                         *
 * ******************************************************************************************* */


ParticleDataBaseCyl::ParticleDataBaseCyl()
{
    _imp = new ParticleDataBaseCylImp;
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBaseCyl::ParticleDataBaseCyl( const ParticleDataBaseCyl &pdb )
{
    _imp = new ParticleDataBaseCylImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBaseCyl::~ParticleDataBaseCyl()
{
    delete _imp;
}


const ParticleDataBaseCyl &ParticleDataBaseCyl::operator=( const ParticleDataBaseCyl &pdb )
{
    delete _imp;
    _imp = new ParticleDataBaseCylImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );

    return( *this );
}


ParticleCyl &ParticleDataBaseCyl::particle( uint32_t i )
{
    return( _imp->particle( i ) );
}


const ParticleCyl &ParticleDataBaseCyl::particle( uint32_t i ) const
{
    return( _imp->particle( i ) ) ;
}


const ParticlePCyl &ParticleDataBaseCyl::trajectory_point( uint32_t i, uint32_t j ) const
{
    return( _imp->trajectory_point( i, j ) );
}


void ParticleDataBaseCyl::add_particle( double IQ, double q, double m, const ParticlePCyl &x )
{
    _imp->add_particle( IQ, q, m, x );
}


void ParticleDataBaseCyl::add_particle( const ParticleCyl &p )
{
    _imp->add_particle( p );
}


void ParticleDataBaseCyl::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
						     double v, double dvp, double dvt, 
						     double x1, double y1, double x2, double y2 )
{

}


void ParticleDataBaseCyl::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						   double E, double Tp, double Tt, 
						   double x1, double y1, double x2, double y2 )
{

}


void ParticleDataBaseCyl::add_2d_full_gaussian_beam( uint32_t N, double I, double q, double m,
						     double Ex, double Tp, double Tt, 
						     double x0, double dr )
{

}


void ParticleDataBaseCyl::add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
							       double a, double b, double e,
							       double Ex, double x0 )
{

}


void ParticleDataBaseCyl::debug_print( std::ostream &os ) const
{
    _imp->debug_print( os );
}


/* ******************************************************************************************* *
 * ParticleDataBase3D                                                                          *
 * ******************************************************************************************* */


ParticleDataBase3D::ParticleDataBase3D()
{
    _imp = new ParticleDataBase3DImp;
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase3D::ParticleDataBase3D( const ParticleDataBase3D &pdb )
{
    _imp = new ParticleDataBase3DImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase3D::~ParticleDataBase3D()
{
    delete _imp;
}


const ParticleDataBase3D &ParticleDataBase3D::operator=( const ParticleDataBase3D &pdb )
{
    delete _imp;
    _imp = new ParticleDataBase3DImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );

    return( *this );
}


Particle3D &ParticleDataBase3D::particle( uint32_t i )
{
    return( _imp->particle( i ) );
}


const Particle3D &ParticleDataBase3D::particle( uint32_t i ) const
{
    return( _imp->particle( i ) ) ;
}


const ParticleP3D &ParticleDataBase3D::trajectory_point( uint32_t i, uint32_t j ) const
{
    return( _imp->trajectory_point( i, j ) );
}


void ParticleDataBase3D::add_particle( double IQ, double q, double m, const ParticleP3D &x )
{
    _imp->add_particle( IQ, q, m, x );
}


void ParticleDataBase3D::add_particle( const Particle3D &p )
{
    _imp->add_particle( p );
}


void ParticleDataBase3D::add_cylindrical_beam_with_velocity( uint32_t N, double J, double q, double m, 
							     double v, double dvp, double dvt, Vec3D c, 
							     Vec3D dir1, Vec3D dir2, double r )
{
    _imp->add_cylindrical_beam_with_velocity( N, J, q, m, v, dvp, dvt, c, dir1, dir2, r );
}


void ParticleDataBase3D::add_cylindrical_beam_with_energy( uint32_t N, double J, double q, double m, 
							   double E, double Tp, double Tt, Vec3D c,
							   Vec3D dir1, Vec3D dir2, double r )
{
    _imp->add_cylindrical_beam_with_energy( N, J, q, m, E, Tp, Tt, c, dir1, dir2, r );
}


void ParticleDataBase3D::add_rectangular_beam_with_velocity( uint32_t N, double J, double q, double m, 
							     double v, double dvp, double dvt, Vec3D c, 
							     Vec3D dir1, Vec3D dir2, double size1, double size2 )
{
    _imp->add_rectangular_beam_with_velocity( N, J, q, m, v, dvp, dvt, c, dir1, dir2, size1, size2 );
}


void ParticleDataBase3D::add_rectangular_beam_with_energy( uint32_t N, double J, double q, double m, 
							   double E, double Tp, double Tt, Vec3D c, 
							   Vec3D dir1, Vec3D dir2, double size1, double size2 )
{
    _imp->add_rectangular_beam_with_energy( N, J, q, m, E, Tp, Tt, c, dir1, dir2, size1, size2 );
}


void ParticleDataBase3D::add_3d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
							double ay, double by, double ey,
							double az, double bz, double ez,
							double Ex, double x0, double y0, double z0 )
{
    _imp->add_3d_KV_beam_with_emittance( N, I, q, m, ay, by, ey, az, bz, ez, Ex, x0, y0, z0 );
}


void ParticleDataBase3D::add_3d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
							      double ay, double by, double ey,
							      double az, double bz, double ez,
							      double Ex, double x0, double y0, double z0 )
{
    _imp->add_3d_gaussian_beam_with_emittance( N, I, q, m, ay, by, ey, az, bz, ez, Ex, x0, y0, z0 );
}


void ParticleDataBase3D::trajectories_at_free_plane( TrajectoryDiagnosticData &tdata, 
						     Vec3D c, Vec3D o, Vec3D p,
						     const std::vector<trajectory_diagnostic_e> &diagnostics ) const
{
    _imp->trajectories_at_free_plane( tdata, c, o, p, diagnostics );
}


void ParticleDataBase3D::export_path_manager_data( std::string filename, 
						   double ref_E, double ref_q, double ref_m, 
						   Vec3D c, Vec3D o, Vec3D p ) const
{
    _imp->export_path_manager_data( filename, ref_E, ref_q, ref_m, c, o, p );
}


void ParticleDataBase3D::debug_print( std::ostream &os ) const
{
    _imp->debug_print( os );
}

