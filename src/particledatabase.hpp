/*! \file particledatabase.hpp
 *  \brief Header file for particledatabase.hpp
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

#ifndef PARTICLEDATABASE_HPP
#define PARTICLEDATABASE_HPP 1


#include <vector>
#include "timer.hpp"
#include "ibsimu.hpp"
#include "trajectory.hpp"
#include "particles.hpp"
#include "particleiterator.hpp"
#include "trajectorydiagnostics.hpp"



/* ******************************************************************************************* *
 * ParticleDataBase classes                                                                    *
 * ******************************************************************************************* */


/*! \brief %Particle database base class.
 *
 * %Particle database base class holds the definitions of particle
 * iteration parameters. Base class also provides a possibility for
 * general pointer to particle database and virtual functions for
 * accessing particles.
 */
class ParticleDataBase {

protected:

    uint32_t       _threadcount; /*!< \brief Number of threads used for calculation. */
    double         _epsabs;      /*!< \brief Absolute error limit for calculation. */
    double         _epsrel;      /*!< \brief Relative error limit for calculation. */
    bool           _polyint;     /*!< \brief Use polynomial(true)/linear(false) interpolation. */
    uint32_t       _maxsteps;    /*!< \brief Maximum number of steps to calculate. */
    double         _maxt;        /*!< \brief Maximum particle time in simulation. */
    uint32_t       _trajdiv;     /*!< \brief Divisor for saved trajectories,
				  * if 3, every third trajectory is saved. */
    bool           _mirror[6];   /*!< \brief Boundary particle mirroring. */

    double         _rhosum;      /*!< \brief Sum of space charge density in defined beams (C/m3). */

    uint32_t       _end_time;    /*!< \brief Number of time limited particle iterations. */
    uint32_t       _end_step;    /*!< \brief Number of step count limited particle iterations. */
    uint32_t       _end_out;     /*!< \brief Number of particle iterations ended with particle 
				  *   out of geometry. */
    uint32_t       _end_coll;    /*!< \brief Number of particle iterations ended with particle 
				  *   collided to an electrode. */
    uint32_t       _end_baddef;  /*!< \brief Number of bad particle definitions. */
    uint32_t       _sum_steps;   /*!< \brief Total number of steps taken. */

    int            _iteration;   /*!< \brief Iteration number. */

    /*! \brief Constructor.
     */
    ParticleDataBase()
	: _threadcount(1), _epsabs(1e-6), _epsrel(1e-6), _polyint(true), _maxsteps(1000), 
	  _maxt(1e-3), _trajdiv(1), _rhosum(0.0), _end_time(0), _end_step(0), _end_out(0), 
	  _end_coll(0), _end_baddef(0), _sum_steps(0), _iteration(-1) {
	_mirror[0] = false;
	_mirror[1] = false;
	_mirror[2] = false;
	_mirror[3] = false;
	_mirror[4] = false;
	_mirror[5] = false;
    }

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Virtual destructor.
     */
    virtual ~ParticleDataBase() {}

/* ****************************************** *
 * Particle iteration settings and statictics *
 * ****************************************** */

    /*! \brief Set the number of threads used for calculation.
     */
    void set_thread_count( uint32_t threadcount ) {
	if( threadcount <= 0 )
	    throw( Error( ERROR_LOCATION, "invalid parameter" ) );
	_threadcount = threadcount;
    }

    /*! \brief Set the accuracy requirement for calculation.
     *
     *  Accuracy requirements default to \a epsabs = 1.0e-6 and \a
     *  epsrel = 1.0e-6.
     */
    void set_accuracy( double epsabs, double epsrel ) {
	_epsabs = epsabs;
	_epsrel = epsrel;
    }

    /*! \brief Set the interpolation type to polynomial(true) or linear(false).
     *
     *  Polynomial interpolation is the default.
     */
    void set_polyint( bool polyint ) {
	_polyint = polyint;
    }
    
    /*! \brief Set maximum number of steps to iterate.
     *
     *  One thousand (1000) steps is the default
     */
    void set_max_steps( uint32_t maxsteps ) {
	if( maxsteps <= 0 )
	    throw( Error( ERROR_LOCATION, "invalid parameter" ) );
	_maxsteps = maxsteps;
    }

    /*! \brief Set maximum lifetime of particle in simulation.
     *
     *  One millisecond (1e-3 sec) is the default
     */
    void set_max_time( double maxt ) {
	if( maxt <= 0.0 )
	    throw( Error( ERROR_LOCATION, "invalid parameter" ) );
	_maxt = maxt;
    }

    /*! \brief Set trajectory saving.
     *
     *  If \a div is zero, no trajectories are saved.
     *  If \a div is one, every trajectory is saved.
     *  If \a div N>1, every Nth trajectory is saved.
     */
    void set_save_trajectories( uint32_t div ) {
	_trajdiv = div;
    }

    /*! \brief Set particle mirroring on boundaries.
     *
     *  Mirroring is set for (xmin,xmax,ymin,ymax,zmin,zmax) borders.
     */
    void set_mirror( const bool mirror[6] ) {
	_mirror[0] = mirror[0];
	_mirror[1] = mirror[1];
	_mirror[2] = mirror[2];
	_mirror[3] = mirror[3];
	_mirror[4] = mirror[4];
	_mirror[5] = mirror[5];
    }

    /*! \brief Get particle mirroring on boundaries.
     *
     *  Mirroring is read for (xmin,xmax,ymin,ymax,zmin,zmax) borders.
     */
    void get_mirror( bool mirror[6] ) const {
	mirror[0] = _mirror[0];
	mirror[1] = _mirror[1];
	mirror[2] = _mirror[2];
	mirror[3] = _mirror[3];
	mirror[4] = _mirror[4];
	mirror[5] = _mirror[5];
    }

    int get_iteration_number( void ) const {
	return( _iteration );
    }

    /*! \brief Return sum of defined beam space charge density.
     *
     *  Returns the summed beam space charge density for all beams
     *  defined with "add_beam" functions. Does not work with
     *  individually added particle trajectories. \a rhosum is cleared
     *  with particle database clearing function clear().
     *
     *  Can be used to program plasma electron density with
     *  EpotProblem::set_pexp_plasma() for example. Please note that
     *  it gives to accumulated charge density which might be
     *  incorrect for multi-beam extraction simulation defined with
     *  several calls to "add_beam" functions.
     */
    double get_rhosum( void ) {
	return( _rhosum );
    }

/* ************************************** *
 * Information and queries                *
 * ************************************** */

    /*! \brief Returns particle count.
     */
    virtual size_t size( void ) const = 0;

    /*! \brief Returns a reference to particle \a i.
     */
    virtual ParticleBase &particle( uint32_t i ) = 0;

    /*! \brief Returns a const reference to particle \a i.
     */
    virtual const ParticleBase &particle( uint32_t i ) const = 0;

    /*! \brief Returns number of trajectory points for particle \a i.
     */
    virtual size_t traj_size( uint32_t i ) const = 0;
    
    /*! \brief Gets the particle \a i trajectory point \a j into \a vel, \a loc and \a t.
     */
    virtual void trajectory_point( double &t, Vec3D &loc, Vec3D &vel, uint32_t i, uint32_t j ) const = 0;

    /*! \brief Gets trajectory diagnostic \a diagnostics at plane \a
     *  axis = \a val in trajectory diagnostic data object \a tdata.
     */
    virtual void trajectories_at_plane( TrajectoryDiagnosticData &tdata, 
					coordinate_axis_e axis,
					double val,
					const std::vector<trajectory_diagnostic_e> &diagnostics ) const = 0;

/* ************************************** *
 * Particle and trajectory clearing       *
 * ************************************** */

    /*! \brief Clears the particle database of all particles.
     *
     *  Clears the database of particles. Also clears beam space
     *  charge sum.
     */
    virtual void clear( void ) = 0;

    /*! \brief Clears the particle trajectory database.
     *
     *  The particle definitions are conserved, but existing
     *  trajectories are cleared.
     */
    virtual void clear_trajectories( void ) = 0;

/* ************************************** *
 * Particle definition                    *
 * ************************************** */

    /*! \brief Reserve memory for \a size particles.
     */
    virtual void reserve( size_t size ) = 0;

/* ************************************** *
 * Debugging, plotting and saving         *
 * ************************************** */

    virtual void debug_print( void ) const = 0;

};


/*! \brief Templated particle database base class for different
    dimensions.
 *
 *  %Particle database holds the definitions of particles, particle
 *  iteration parameters and possibly the trajectories of particles if
 *  the particle iterator has saved them. %Particle database also
 *  provides a variety of convenience functions for defining particle
 *  beams.
 *
 *  Specialized derived classes ParticleDataBase2D,
 *  ParticleDataBaseCyl and ParticleDataBase3D exist for use with
 *  different dimensionalities. The templated %ParticleDataBasePP
 *  class provides all the common functions.
 */
template<class PP> class ParticleDataBasePP : public ParticleDataBase {

    /*! \brief Add requested diagnostics to \a tdata from particle point \a x.
     */
    static void add_diagnostics( TrajectoryDiagnosticData &tdata, const PP &x, 
				 const Particle<PP> &p, int crd ) {
	//std::cout << "add_diagnostics():\n";
	for( size_t a = 0; a < tdata.diag_size(); a++ ) {
	    //std::cout << "  diagnostic[" << a << "] = " << tdata.diagnostic(a) << "\n";
	    
	    double data = 0.0;
	    switch( tdata.diagnostic( a ) ) {
	    case DIAG_NONE:
		data = 0.0;
		break;
	    case DIAG_T:
		data = x[0];
		break;
	    case DIAG_X:
		data = x[1];
		break;
	    case DIAG_VX:
		data = x[2];
		break;
	    case DIAG_Y:
	    case DIAG_R:
		data = x[3];
		break;
	    case DIAG_VY:
	    case DIAG_VR:
		data = x[4];
		break;
	    case DIAG_Z:
		data = x[5];
		break;
	    case DIAG_VZ:
		data = x[6];
		break;
	    case DIAG_W:
		data = x[5];
		break;
	    case DIAG_VTHETA:
		data = x[5]*x[3];
		break;
	    case DIAG_XP:
		data = x[2]/x[2*crd+2];
		break;
	    case DIAG_YP:
	    case DIAG_RP:
		data = x[4]/x[2*crd+2];
		break;
	    case DIAG_AP:
		data = x[3]*x[5]/x[2*crd+2];
		break;
	    case DIAG_ZP:
		data = x[6]/x[2*crd+2];
		break;
	    case DIAG_CURR:
		data = p.IQ();
		break;
	    case DIAG_QM:
		data = p.qm();
		break;
	    case DIAG_EK:
		Vec3D velocity = x.velocity();
		data = velocity.norm2();
		break;
	    }
	    //std::cout << "  adding data = " << data << "\n";
	    tdata.add_data( a, data );
	}
    }

protected:

    std::vector<Particle<PP> > _particles;      /*!< \brief Vector of type PP particles */

    /*! \brief Constructor.
     */
    ParticleDataBasePP() {}

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Destructor.
     */
    ~ParticleDataBasePP() {}

/* ************************************** *
 * Information and queries                *
 * ************************************** */

    /*! \brief Returns particle count.
     */
    virtual size_t size( void ) const { return( _particles.size() ); }

    /*! \brief Returns a reference to particle \a i.
     */
    virtual Particle<PP> &particle( uint32_t i ) { return( _particles[i] ); }

    /*! \brief Returns a const reference to particle \a i.
     */
    virtual const Particle<PP> &particle( uint32_t i ) const { return( _particles[i] ); }
    
    /*! \brief Returns number of trajectory points for particle \a i.
     */
    virtual size_t traj_size( uint32_t i ) const { return( _particles[i].traj_size() ); }
    
    /*! \brief Gets the particle \a i trajectory point \a j into \a vel, \a loc and \a t.
     */
    virtual void trajectory_point( double &t, Vec3D &loc, Vec3D &vel, uint32_t i, uint32_t j ) const {
	PP x = _particles[i].traj(j); 
	t = x[0];
	loc = x.location();
	vel = x.velocity();
    }

    /*! \brief Gets trajectory diagnostic \a diagnostics at plane \a
     *  axis = \a val in trajectory diagnostic data object \a tdata.
     */
    virtual void trajectories_at_plane( TrajectoryDiagnosticData &tdata, 
					coordinate_axis_e axis,
					double val,
					const std::vector<trajectory_diagnostic_e> &diagnostics ) const {

	if( ibsimu.get_verbose_output() )
	    std::cout << "Making trajectory diagnostics at " 
		      << coordinate_axis_string[axis] << " = " << val << "\n";

	// Check query
	switch( PP::geom_mode() ) {
	case MODE_1D:
	    throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	    break;
	case MODE_2D:
	    if( axis == AXIS_R || axis == AXIS_Z )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	case MODE_CYL:
	    if( axis == AXIS_Y || axis == AXIS_Z )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	case MODE_3D:
	    if( axis == AXIS_R )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	}

	// Check diagnostics query validity
	for( size_t a = 0; a < diagnostics.size(); a++ ) {
	    if( diagnostics[a] == DIAG_NONE )
		throw( Error( ERROR_LOCATION, "invalid diagnostics query \'DIAG_NONE\'" ) );
	    else if( PP::geom_mode() != MODE_CYL && (diagnostics[a] == DIAG_R ||
						     diagnostics[a] == DIAG_VR ||
						     diagnostics[a] == DIAG_RP ||
						     diagnostics[a] == DIAG_W ||
						     diagnostics[a] == DIAG_VTHETA ||
						     diagnostics[a] == DIAG_AP) )
		throw( Error( ERROR_LOCATION, "invalid diagnostics query for geometry type" ) );
	}

	// Prepare output vector
	tdata.clear();
	for( size_t a = 0; a < diagnostics.size(); a++ ) {
	    tdata.add_data_column( diagnostics[a] );
	}
	
	// Set coordinate index
	int crd;
	switch( axis ) {
	case AXIS_X:
	    crd = 0;
	    break;
	case AXIS_Y:
	case AXIS_R:
	    crd = 1;
	    break;
	case AXIS_Z:
	    crd = 2;
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unsupported axis" ) );
	}

	// Scan through particle trajectory points
	double Isum = 0.0;
	std::vector<PP> intsc;
	for( size_t a = 0; a < _particles.size(); a++ ) {
	    size_t N = _particles[a].traj_size();
	    if( N < 2 )
		continue;
	    PP x1 = _particles[a].traj(0);
	    size_t nintsum = 0;
	    for( size_t b = 1; b < N; b++ ) {
		PP x2 = _particles[a].traj(b);
		intsc.clear();
		size_t nintsc = PP::trajectory_intersections_at_plane( intsc, crd, val, x1, x2 );
		nintsum += nintsc;
		for( size_t c = 0; c < nintsc; c++ ) {
		    Isum += _particles[a].IQ();
		    add_diagnostics( tdata, intsc[c], _particles[a], crd );
		}

		x1 = x2;
	    }
	}

	if( ibsimu.get_verbose_output() ) {
	    std::cout << "  number of trajectories = " << tdata.traj_size() << "\n";
	    if( PP::geom_mode() == MODE_2D )
		std::cout << "  total current = " << Isum << " A/m\n";
	    else
		std::cout << "  total current = " << Isum << " A\n";
	}
    }

/* ************************************** *
 * Particle and trajectory clearing       *
 * ************************************** */

    /*! \brief Clears the particle database of all particles.
     */
    virtual void clear( void ) { 
	_particles.clear();
	_rhosum = 0.0;
    }

    /*! \brief Clears the particle trajectory database.
     *
     *  The particle definitions are conserved, but existing
     *  trajectories are cleared.
     */
    virtual void clear_trajectories( void ) {
	for( uint32_t a = 0; a < _particles.size(); a++ )
	    _particles[a].clear_trajectory();
    }

/* ************************************** *
 * Particle definition                    *
 * ************************************** */

    /*! \brief Reserve memory for \a size particles.
     */
    virtual void reserve( size_t size ) { _particles.reserve( size ); }

    /*! \brief Add one particle.
     *
     *  Adds one particle to particle database. Particle properties
     *  are: \a IQ is the current (A) in time-independent or charge
     *  (C) in time-dependent simulations carried by the particle
     *  cloud that the simulated particle represents, \a q is the
     *  charge state of the microscopic particle (in multiples of e),
     *  \a m is the mass of the microscopic particle (u) and \a x
     *  contains the time, position (m) and velocity (m/s) of the
     *  particle.
     */
    void add_particle( double IQ, double q, double m, const PP &x ) {
	_particles.push_back( Particle<PP>( IQ, CHARGE_E*q, MASS_U*m, x ) );
    }

/* ************************************** *
 * Particle iterators                     *
 * ************************************** */

    /*! \brief Iterate particles through the geometry.
     *
     *  The particles defined in particle database \a pdb are iterated
     *  through electric field \a efield and magnetic field \a bfield 
     *  in geometry \a g. Space charge density field \a scharge is set
     *  from the particle trajectories.
     */
    virtual void iterate_trajectories( ScalarField &scharge, const Efield &efield, 
				       const VectorField &bfield, const Geometry &g ) {

	ScalarField                         *schmap[_threadcount];
	std::vector<ParticleIterator<PP> *>  iterators;

	Timer t;
	if( ibsimu.get_verbose_output() )
	    std::cout << "Calculating particle trajectories\n";
	_iteration++;

	// Clear space charge
	scharge.clear();

	// Clear statistics
	_end_time   = 0;
	_end_step   = 0;
	_end_out    = 0;
	_end_coll   = 0;
	_end_baddef = 0;
	_sum_steps  = 0;

	// Check number of particles
	if( _particles.size() == 0 ) {
	    std::cout << "  no particles to calculate\n";
	    return;
	}

	// Make separate space charge maps for all threads and build iterators
	for( uint32_t a = 0; a < _threadcount; a++ ) {
	    if( a == 0 ) schmap[a] = &scharge;
	    else schmap[a] = new ScalarField( g );
	    iterators.push_back( new ParticleIterator<PP>( _epsabs, _epsrel, _polyint, _maxsteps, 
							   _maxt, _trajdiv, _mirror, schmap[a], 
							   &efield, &bfield, &g, &_particles[0] ) );
	}

	// Make Scheduler
	Scheduler<ParticleIterator<PP>,Particle<PP>,Error> scheduler( iterators );

	// Add problems
	for( size_t a = 0; a < _particles.size(); a++ )
	    scheduler.add_problem( &_particles[a] );

	// Wait for completition
	scheduler.run();
	scheduler.finish();

	if( scheduler.is_error() ) {
	    // Throw the error
	    std::vector<Error> err;
	    std::vector<Particle<PP> *> part;
	    scheduler.get_errors( err, part );
	    throw( err[0] );
	}

	// Combine separate space charge maps and collect
	// statistics. Free all allocated memory.
	for( uint32_t a = 0; a < _threadcount; a++ ) {
	    if( a != 0 ) {
		scharge += *schmap[a];
		delete schmap[a];
	    }
	    uint32_t end_time, end_step, end_out, end_coll, end_baddef, sum_steps;
	    iterators[a]->get_statistics( end_time, end_step, end_out, 
					  end_coll, end_baddef, sum_steps );
	    _end_time   += end_time;
	    _end_step   += end_step;
	    _end_out    += end_out;
	    _end_coll   += end_coll;
	    _end_baddef += end_baddef;
	    _sum_steps  += sum_steps;
	    delete iterators[a];
	}

	scharge_finalize( scharge );
	
	t.stop();
	if( ibsimu.get_verbose_output() ) {
	    std::cout << "  Particle histories (" << _particles.size() << " total):\n";
	    std::cout << "    time limited = " << _end_time << "\n";
	    std::cout << "    step count limited = " << _end_step << "\n";
	    std::cout << "    out of geometry = " << _end_out << "\n";
	    std::cout << "    collided = " << _end_coll << "\n";
	    std::cout << "    bad definitions = " << _end_baddef << "\n";
	    std::cout << "    total steps = " << _sum_steps << "\n";
	    std::cout << "    steps per particle (ave) = " << _sum_steps/(double)_particles.size() << "\n";
	    std::cout << "  time used = " << t << "\n";
	}
    }

    /*! \brief Step particles forward by time step dt.
     *
     *  The particles defined in particle database \a pdb are stepped
     *  forward one time step in electric field \a efield and geometry
     *  \a g.
     */
    virtual void step_particles( const Efield &efield, const Geometry &g, double dt ) {

    }

/* ************************************** *
 * Debugging, plotting and saving         *
 * ************************************** */

    /*! \brief Prints internal data to std::cout.
     */
    virtual void debug_print( void ) const {
	std::cout << "threadcount = " << _threadcount << "\n";
	std::cout << "epsabs = "      << _epsabs << "\n";
	std::cout << "epsrel = "      << _epsrel << "\n";
	std::cout << "maxsteps = "    << _maxsteps << "\n";
	std::cout << "maxt = "        << _maxt << "\n";
	std::cout << "trajdiv = "     << _trajdiv << "\n";
	std::cout << "mirror = (";
	for( uint32_t a = 0; a < 5; a++ )
	    std::cout << _mirror[a] << ", ";
	std::cout << _mirror[5] << ")\n";
	
	for( uint32_t a = 0; a < _particles.size(); a++ ) {
	    std::cout << "Particle " << a << ":\n";
	    _particles[a].debug_print();
	}
    }

};



/*! \brief %Particle database class for two dimensions.
 *
 *  %Particle database holds the definitions of particles and possibly
 *  the trajectories of particles it the particle iterator has saved
 *  them. %ParticleDataBase2D provides a variety of convenience
 *  functions for defining particle beams.
 *
 *  Particles are always stored in the database in the order they are
 *  defined. When reading back the simulation results, the order can
 *  be used to identify the particles.
 */
class ParticleDataBase2D : public ParticleDataBasePP<ParticleP2D> {


public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Constructor.
     */
    ParticleDataBase2D() {}

    /*! \brief Destructor.
     */
    ~ParticleDataBase2D() {}

/* ************************************** *
 * Particle beam definition               *
 * ************************************** */

    /*! \brief Add a 2d beam with energies.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean energy \a E (eV). The beam has parallel
     *  temperature \a Tp (eV) and transverse temperature \a Tt (eV).
     *
     *  The particle speeds of the beam in direction \a i are sampled
     *  from a gaussian distribution with standard deviation dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
				  double E, double Tp, double Tt, 
				  double x1, double y1, double x2, double y2 );

    /*! \brief Add a 2d beam with velocities.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean velocity \a v (m/s). The beam has parallel
     *  gaussian velocity distribution with standard deviation \a dvp
     *  (m/s) and transverse gaussian velocity distribution with
     *  standard deviation \a dvt (m/s).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
				    double v, double dvp, double dvt, 
				    double x1, double y1, double x2, double y2 );


    /*! \brief Add a 2d beam with defined KV emittance.
     *
     */
    void add_2d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
					double a, double b, double e,
					double Ex, double x0, double y0 );

    /*! \brief Add a 2d beam with defined gaussian emittance.
     *
     */
    void add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double a, double b, double e,
					      double Ex, double x0, double y0 );
};




/*! \brief %Particle database class for cylindrical systems.
 *
 *  %Particle database holds the definitions of particles and possibly
 *  the trajectories of particles it the particle iterator has saved
 *  them. %ParticleDataBaseCyl provides a variety of convenience
 *  functions for defining particle beams.
 *
 *  Particles are always stored in the database in the order they are
 *  defined. When reading back the simulation results, the order can
 *  be used to identify the particles.
 */
class ParticleDataBaseCyl : public ParticleDataBasePP<ParticlePCyl> {


public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Constructor.
     */
    ParticleDataBaseCyl() {}

    /*! \brief Destructor.
     */
    ~ParticleDataBaseCyl () {}

/* ************************************** *
 * Particle beam definition               *
 * ************************************** */

    /*! \brief Add a 2d beam with energies.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean energy \a E (eV). The beam has parallel
     *  temperature \a Tp (eV) and transverse temperature \a Tt (eV).
     *
     *  The particle speeds of the beam in direction \a i are sampled
     *  from a gaussian distribution with standard deviation dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
				  double E, double Tp, double Tt, 
				  double x1, double y1, double x2, double y2 );

    /*! \brief Add a 2d beam with velocities.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean velocity \a v (m/s). The beam has parallel
     *  gaussian velocity distribution with standard deviation \a dvp
     *  (m/s) and transverse gaussian velocity distribution with
     *  standard deviation \a dvt (m/s).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
				    double v, double dvp, double dvt, 
				    double x1, double y1, double x2, double y2 );

    /*! \brief Add a 2d beam with defined gaussian emittance.
     *
     *  The beam is made to match Twiss parameters \f$ \alpha \f$ (a), \f$ \beta \f$ (b), 
     *  \f$ \epsilon \f$ (e) in projectional directions (y,y') and (z,z'). The particles 
     *  are defined to have equal currents. This is different from other beam definitions 
     *  for cylindrical coordinates. This might change soon.
     *
     *  \todo Redo particle distribution definition from emittance for cylindrical 
     *  symmetry systems.
     *
     */
    void add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double a, double b, double e,
					      double Ex, double x0 );
};



/*! \brief %Particle database class for three dimensions.
 *
 *  %Particle database holds the definitions of particles and possibly
 *  the trajectories of particles it the particle iterator has saved
 *  them. %ParticleDataBase3D provides a variety of convenience
 *  functions for defining particle beams.
 *
 *  Particles are always stored in the database in the order they are
 *  defined. When reading back the simulation results, the order can
 *  be used to identify the particles.
 */
class ParticleDataBase3D : public ParticleDataBasePP<ParticleP3D> {


public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Constructor.
     */
    ParticleDataBase3D() {}

    /*! \brief Destructor.
     */
    ~ParticleDataBase3D() {}

/* ************************************** *
 * Particle beam definition               *
 * ************************************** */

    /*! \brief Add a cylindrical beam with energies.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particles is \a q (in
     *  multiples of e), mass is \a m (u). The beam starting surface
     *  is a disc of radius \a r centered at \a c. The normal
     *  direction of the disc is \a dir3 = \a dir1 x \a dir2. The
     *  first tangent direction is \a dir1 and the second is \a dir1 x
     *  \a dir3. If you want beam to go to positive x-direction, \a
     *  dir1 could be (0,1,0) and \a dir2 (0,0,1) for example.  The
     *  beam energy \a E (eV) is defined in the normal
     *  direction. Temperatures are defined in normal (parallel)
     *  direction and transverse direction as \a Tp (eV) and \a Tt
     *  (eV), respectively.
     *
     *  The particle speeds of the beam in direction \a i are sampled
     *  from a gaussian distribution with standard deviation dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_cylindrical_beam_with_energy( uint32_t N, double J, double q, double m, 
					   double E, double Tp, double Tt, Vec3D c, 
					   Vec3D dir1, Vec3D dir2, double r );

    /*! \brief Add a 3d beam with velocities.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particles is \a q (in
     *  multiples of e), mass is \a m (u). The beam starting surface
     *  is a disc of radius \a r centered at \a c. The normal
     *  direction of the disc is \a dir3 = \a dir1 x \a dir2. The
     *  first tangent direction is \a dir1 and the second is \a dir1 x
     *  \a dir3. The beam velocity \a v (m/s) is defined in the normal
     *  direction. The velocity distribution standard deviation are
     *  defined in normal (parallel) direction and transverse
     *  directions as \a dvp (m/s) and \a dvt (m/s), respectively
     *
     *  The particle velocities of the beam in direction \a i are
     *  sampled from a gaussian distribution with standard deviation
     *  dv_i, which is related to beam temperature by dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_cylindrical_beam_with_velocity( uint32_t N, double J, double q, double m, 
					     double v, double dvp, double dvt, Vec3D c, 
					     Vec3D dir1, Vec3D dir2, double r );


    /*! \brief Add a 3d beam with defined KV emittance.
     *
     */
    void add_3d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
					double ay, double by, double ey,
					double az, double bz, double ez,
					double Ex, double x0, double y0, double z0 );

    /*! \brief Add a 3d beam with defined gaussian emittance.
     *
     */
    void add_3d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double ay, double by, double ey,
					      double az, double bz, double ez,
					      double Ex, double x0, double y0, double z0 );
};






#endif















