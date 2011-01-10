/*! \file particleiterator.hpp
 *  \brief %Particle iterator
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

#ifndef PARTICLEITERATOR_HPP
#define PARTICLEITERATOR_HPP 1


#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_poly.h>
#include "geometry.hpp"
#include "particles.hpp"
#include "efield.hpp"
#include "scalarfield.hpp"
#include "scharge.hpp"
#include "scheduler.hpp"
#include "polysolver.hpp"


//#define DEBUG_PARTICLE_ITERATOR 1


/*! \brief %Particle iterator type.
 */
enum particle_iterator_type_e {
    PARTICLE_ITERATOR_ADAPTIVE = 0,
    PARTICLE_ITERATOR_FIXED_STEP_LEN
};


/*! \brief %Mesh intersection (collision) coordinate data
 *
 * Contains one coordinate data and the direction of particle for one
 * intersection between particle trajectory and mesh plane. Templated
 * for particle point type (see ParticlePBase).
 *
 * Templated for particle point type (see ParticlePBase).
 */
template <class PP> class ColData {
public:
    PP                _x;         /*!< \brief %Mesh intersection coordinates. */
    int               _dir;       /*!< \brief Direction of particle at intersection.
				   *  i: -1/+1, j: -2/+2, k: -3:/+3. */
    
    /*! \brief Constructor for collision at \a x into direction \a dir.
     */
    ColData( PP x, int dir ) : _x(x), _dir(dir) {}
    
    /*! \brief Compare coldata entry times.
     *
     *  Used for sorting coldata entries.
     */
    bool operator<( const ColData &cd ) const {
	return( _x[0] < cd._x[0] );
    }

    /*! \brief Find mesh intersections of linearly interpolated
     *  particle trajectory segment.
     *
     *  Makes a linear interpolation between points \a x1 and \a x2
     *  and searches intersection points of this line and \a
     *  mesh. Intersection points are saved to vector \a coldata in
     *  increasing time order.
     */
    static void build_coldata_linear( std::vector<ColData> &coldata, const Mesh &mesh,
				      const PP &x1, const PP &x2 ) {
	
	coldata.clear();

	for( size_t a = 0; a < PP::dim(); a++ ) {
	    
            int a1 = (int)floor( (x1[2*a+1]-mesh.origo(a))/mesh.h() );
            int a2 = (int)floor( (x2[2*a+1]-mesh.origo(a))/mesh.h() );
            if( a1 > a2 ) {
                int a = a2;
                a2 = a1;
                a1 = a;
            }
	    
            for( int b = a1+1; b <= a2; b++ ) {
        
                // Save intersection coordinates
                double K = (b*mesh.h() + mesh.origo(a) - x1[2*a+1]) / 
                    (x2[2*a+1] - x1[2*a+1]);
                if( K < 0.0 ) K = 0.0;
                else if( K > 1.0 ) K = 1.0;
                //std::cout << "Found valid root: " << K << "\n";

                if( x2[2*a+1] > x1[2*a+1] )
                    coldata.push_back( ColData( x1 + (x2-x1)*K, a+1 ) );
                else
                    coldata.push_back( ColData( x1 + (x2-x1)*K, -a-1 ) );
            }
        }

	// Sort intersections in increasing time order
	sort( coldata.begin(), coldata.end() );
    }

    /*! \brief Find mesh intersections of polynomially interpolated
        particle trajectory segment.
     *
     *  Makes a polynomial interpolation between points \ə x1 and \a
     *  x2 and searches intersection points of this line and \a
     *  mesh. Intersection points are saved to vector \a coldata in
     *  increasing time order.
     */
    static void build_coldata_poly( std::vector<ColData> &coldata, const Mesh &mesh,
				    const PP &x1, const PP &x2 ) {
	
#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Building coldata using polynomial interpolation\n";
#endif

	coldata.clear();

	// Construct trajectory representation
	TrajectoryRep1D traj[PP::dim()];
	for( size_t a = 0; a < PP::dim(); a++ ) {
	    traj[a].construct( x2[0]-x1[0], 
			       x1[2*a+1], x1[2*a+2], 
			       x2[2*a+1], x2[2*a+2] );
	}

	// Solve trajectory intersections
	for( size_t a = 0; a < PP::dim(); a++ ) {

	    // Mesh number of x1 (start point)
	    int i = (int)floor( (x1[2*a+1]-mesh.origo(a))/mesh.h() );
	    
	    // Search to negative (dj = -1) and positive (dj = +1) mesh directions
	    for( int dj = -1; dj <= 1; dj += 2 ) {
		int j = i;
		if( dj == +1 )
		    j = i+1;
		int Kcount;  // Solution counter
		double K[3]; // Solution array
		while( 1 ) {

		    // Intersection point
		    double val = mesh.origo(a) + mesh.h() * j;
		    if( val < mesh.origo(a) )
			break;
		    else if( val > mesh.max(a) )
			break;

#ifdef DEBUG_PARTICLE_ITERATOR
		    std::cout << "  Searching intersections at coord(" << a << ") = " << val << "\n";
#endif
		    Kcount = traj[a].solve( K, val );
		    if( Kcount == 0 )
			break; // No valid roots

#ifdef DEBUG_PARTICLE_ITERATOR
		    std::cout << "  Found " << Kcount << " valid roots: ";
		    for( int p = 0; p < Kcount; p++ )
			std::cout << K[p] << " ";
		    std::cout << "\n";
#endif

		    // Save roots to coldata
		    for( int b = 0; b < Kcount; b++ ) {
			PP xcol;
			double x, v;
			xcol(0) = x1[0] + K[b]*(x2[0]-x1[0]);
			for( size_t c = 0; c < PP::dim(); c++ ) {
			    traj[c].coord( x, v, K[b] );
			    if( a == c )
				xcol[2*c+1] = val; // limit numerical inaccuracy
			    else
				xcol[2*c+1] = x;
			    xcol[2*c+2] = v;
			}
			if( mesh.geom_mode() == MODE_CYL )
			    xcol[5] = x1[5] + K[b]*(x2[5]-x1[5]);
			if( xcol[2*a+2] >= 0.0 )
			    coldata.push_back( ColData( xcol, a+1 ) );
			else
			    coldata.push_back( ColData( xcol, -a-1 ) );
		    }

		    j += dj;
		}
	    }
	}

	// Sort intersections in increasing time order
	sort( coldata.begin(), coldata.end() );

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "  Coldata built\n";
#endif
    }

};


/*! \brief %Particle iterator class for continuous Vlasov-type iteration.
 *
 * Templated for particle point type (see ParticlePBase).
 *
 *  \todo Detailed documentation needed.
 *  \todo PIC style iterator needed.
 */
template <class PP> class ParticleIterator {

    gsl_odeiv_system      _system;    /**< \brief GSL ODE integrator system. */
    gsl_odeiv_step       *_step;      /**< \brief GSL ODE integrator stepper. */
    gsl_odeiv_control    *_control;   /**< \brief GSL ODE integrator constrol. */
    gsl_odeiv_evolve     *_evolve;    /**< \brief GSL ODE integrator integrator. */

    particle_iterator_type_e _type;   /**< \brief Iteratory type. */

    bool                  _polyint;   /*!< \brief Interpolation type to use. True means use polynomial */
    double                _epsabs;    /*!< \brief Absolute error limit. */
    double                _epsrel;    /*!< \brief Relative error limit. */
    uint32_t              _maxsteps;  /*!< \brief Maximum number of simulation steps for particle. */
    double                _maxt;      /*!< \brief Maximum particle lifetime. */
    uint32_t              _trajdiv;   /*!< \brief Divisor for saved trajectories,
				       * if 3, every third trajectory is saved. */
    bool                  _mirror[6]; /*!< \brief Is particle mirrored on boundary? */

    Particle<PP>         *_first;     /*!< \brief Pointer to first particle of the database. */
    ParticleIteratorData  _pidata;    /*!< \brief User data provided to PP::get_derivatives(). */

    PP                    _xi;        /*!< \brief Previous mesh intersection coordinates 
				       *   or starting point. */
    std::vector<PP>       _traj;      /*!< \brief %Particle trajectory data for current trajectory. */
    std::vector<ColData<PP> >  _coldata;   /*!< \brief Mesh intersection coordinate data. */

    uint32_t              _end_time;  /*!< \brief Number of time limited particle iterations. */
    uint32_t              _end_step;  /*!< \brief Number of step count limited particle iterations. */
    uint32_t              _end_out;   /*!< \brief Number of particle iterations ended with particle 
				       *   out of geometry. */
    uint32_t              _end_coll;  /*!< \brief Number of particle iterations ended with particle 
				       *   collided to an electrode. */
    uint32_t              _end_baddef;/*!< \brief Number of bad particle definitions. */
    uint32_t              _sum_steps; /*!< \brief Total number of steps taken. */


    /*! \brief Check for particle collision with solid
     *
     *  Particle propagates from x1 to x2, where x1 is in
     *  vacuum. Checks if x2 is inside solid and if it is, the
     *  collision point is bracketed between x1 and x2. Particle
     *  coordinates at status_x are set to collision coordinates.
     *  Returns false if particle collided.
     * */
    bool check_collision( Particle<PP> &particle, const PP &x1, const PP &x2, PP &status_x ) {

	size_t a;
	double K;
	Vec3D v1, v2, vc;

	// Convert PP to Vec3D
	for( a = 0; a < (PP::size()-1)/2; a++ ) {
	    v1[a] = x1[2*a+1];
	    v2[a] = x2[2*a+1];
	}

	// If inside solid, bracket for collision point
	if( (a = _pidata._g->inside( v2 )) >= 7 ) {
	    K = _pidata._g->bracket_surface( a, v2, v1, vc );
	} else {
	    return( true ); // No collision happened.
	}

	// Convert Vec3D to PP
	for( a = 0; a < PP::size(); a++ )
	    status_x[a] = x2[a] + K*(x1[a]-x2[a]);

	// Remove all points from _traj after time status_x[0].
	for( a = _traj.size()-1; a > 0; a-- ) {
	    if( _traj[a][0] > status_x[0] )
		_traj.pop_back();
	    else
		break;
	}

	// Save last trajectory point and update status
	_traj.push_back( status_x );
	particle.set_status( PARTICLE_COLL );
	_end_coll++;

	return( false ); // Collision happened.
    }


    /*! \brief Mirror trajectory.
     *
     *  Trajectory is mirrored at \a _coldata[c] on axis \a at \a
     *  border, where -1 is the negative side and +1 is the positive
     *  side.. Already saved trajectory points are checked back to
     *  _xi.
     */
    void handle_mirror( size_t c, int i[3], size_t a, int border, PP &x2 ) {

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "    handle_mirror( c = " << c 
		  << ", i = (" << i[0] << ", " << i[1] << ", " << i[2]
		  << "), a = " << a << ", border = " << border 
		  << ")\n";
#endif

	double xmirror;
	if( border < 0 ) {
	    xmirror = _pidata._g->origo(a);
	    i[a] = -i[a]-1;
	} else {
	    xmirror = _pidata._g->max(a);
	    i[a] = 2*_pidata._g->size(a)-i[a]-3;
	}

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "    xmirror = " << xmirror << "\n";
	std::cout << "    i = (" << i[0] << ", " << i[1] << ", " << i[2] << ")\n";
	std::cout << "    xi = " << _xi << "\n";
#endif
	
	// Check if found edge at first encounter
	bool caught_at_boundary = false;
	if( _coldata[c]._dir == border*((int)a+1) && 
	    ( i[a] == 0 || i[a] == (int)_pidata._g->size(a)-2 ) ) {
	    caught_at_boundary = true;
#ifdef DEBUG_PARTICLE_ITERATOR
	    std::cout << "   caught_at_boundary\n";
#endif
	}

	// Mirror traj back to _xi
	if( caught_at_boundary ) {
	    _traj.push_back( _coldata[c]._x );
	} else {
	    for( int b = _traj.size()-1; b > 0; b-- ) {
		if( _traj[b][0] >= _xi[0] ) {
		    
#ifdef DEBUG_PARTICLE_ITERATOR
		    std::cout << "    mirroring traj[" << b << "] = " << _traj[b] << "\n";
#endif
		    _traj[b][2*a+1] = 2.0*xmirror - _traj[b][2*a+1];
		    _traj[b][2*a+2] *= -1.0;
		} else
		    break;
	    }
	}

	// Mirror rest of the coldata
	for( size_t b = c; b < _coldata.size(); b++ ) {
	    if( (size_t)abs(_coldata[b]._dir) == a+1 )
		_coldata[b]._dir *= -1;
	    _coldata[b]._x[2*a+1] = 2.0*xmirror - _coldata[b]._x[2*a+1];
	    _coldata[b]._x[2*a+2] *= -1.0;
	}

	if( caught_at_boundary )
	    _traj.push_back( _coldata[c]._x );

	// Mirror calculation point
	x2[2*a+1] = 2.0*xmirror - x2[2*a+1];
	x2[2*a+2] *= -1.0;
	
	// Coordinates changed, reset integrator
	gsl_odeiv_step_reset( _step );
	gsl_odeiv_evolve_reset( _evolve );
    }


    void handle_collision( Particle<PP> &particle, size_t c, PP &status_x ) {

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "    handle_collision()\n";
#endif

	_traj.push_back( _coldata[c]._x );
	status_x = _coldata[c]._x;
	particle.set_status( PARTICLE_OUT );
	_end_out++;
    }


    /*! \brief Handle particle mesh intersection.
     *
     *  Particle mesh coordinates \a i are advanced through
     *  intersection described in \a _coldata[c]. Collisions with
     *  solids and particle mirroring on boundaries are handled. Final
     *  particle coordinates \a x2 are updated accordingly. Returns
     *  false if particle collided.
     */
    bool handle_trajectory_advance( Particle<PP> &particle, size_t c, int i[3], PP &x2 ) {

	// Check for collisions with solids and advance coordinates i.
	if( PP::dim() == 2 ) {
	    if( _coldata[c]._dir == -1 ) {
		if( ( abs(_pidata._g->mesh(i[0],  i[1]  )) >= 7 || 
		      abs(_pidata._g->mesh(i[0],  i[1]+1)) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[0]--;
	    } else if( _coldata[c]._dir == +1 ) {
		if( ( abs(_pidata._g->mesh(i[0]+1,i[1]  )) >= 7 || 
		      abs(_pidata._g->mesh(i[0]+1,i[1]+1)) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[0]++;
	    } else if( _coldata[c]._dir == -2 ) {
		if( ( abs(_pidata._g->mesh(i[0],  i[1]  )) >= 7 || 
		      abs(_pidata._g->mesh(i[0]+1,i[1]  )) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[1]--;
	    } else {
		if( ( abs(_pidata._g->mesh(i[0],  i[1]+1)) >= 7 || 
		      abs(_pidata._g->mesh(i[0]+1,i[1]+1)) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[1]++;
	    }
	} else if( PP::dim() == 3 ) {
	    if( _coldata[c]._dir == -1 ) {
		if( ( abs(_pidata._g->mesh(i[0],  i[1],  i[2]  )) >= 7 || 
		      abs(_pidata._g->mesh(i[0],  i[1]+1,i[2]  )) >= 7 ||
		      abs(_pidata._g->mesh(i[0],  i[1],  i[2]+1)) >= 7 ||
		      abs(_pidata._g->mesh(i[0],  i[1]+1,i[2]+1)) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[0]--;
	    } else if( _coldata[c]._dir == +1 ) {
		if( ( abs(_pidata._g->mesh(i[0]+1,i[1],  i[2]  )) >= 7 || 
		      abs(_pidata._g->mesh(i[0]+1,i[1]+1,i[2]  )) >= 7 ||
		      abs(_pidata._g->mesh(i[0]+1,i[1],  i[2]+1)) >= 7 ||
		      abs(_pidata._g->mesh(i[0]+1,i[1]+1,i[2]+1)) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[0]++;
	    } else if( _coldata[c]._dir == -2 ) {
		if( ( abs(_pidata._g->mesh(i[0],  i[1],i[2]  )) >= 7 || 
		      abs(_pidata._g->mesh(i[0]+1,i[1],i[2]  )) >= 7 ||
		      abs(_pidata._g->mesh(i[0],  i[1],i[2]+1)) >= 7 ||
		      abs(_pidata._g->mesh(i[0]+1,i[1],i[2]+1)) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[1]--;
	    } else if( _coldata[c]._dir == +2 ) {
		if( ( abs(_pidata._g->mesh(i[0],  i[1]+1,i[2]  )) >= 7 || 
		      abs(_pidata._g->mesh(i[0]+1,i[1]+1,i[2]  )) >= 7 ||
		      abs(_pidata._g->mesh(i[0],  i[1]+1,i[2]+1)) >= 7 ||
		      abs(_pidata._g->mesh(i[0]+1,i[1]+1,i[2]+1)) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[1]++;
	    } else if( _coldata[c]._dir == -3 ) {
		if( ( abs(_pidata._g->mesh(i[0],  i[1],  i[2]  )) >= 7 || 
		      abs(_pidata._g->mesh(i[0]+1,i[1],  i[2]  )) >= 7 ||
		      abs(_pidata._g->mesh(i[0],  i[1]+1,i[2])) >= 7 ||
		      abs(_pidata._g->mesh(i[0]+1,i[1]+1,i[2])) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[2]--;
	    } else {
		if( ( abs(_pidata._g->mesh(i[0],  i[1],  i[2]+1)) >= 7 || 
		      abs(_pidata._g->mesh(i[0]+1,i[1],  i[2]+1)) >= 7 ||
		      abs(_pidata._g->mesh(i[0],  i[1]+1,i[2]+1)) >= 7 ||
		      abs(_pidata._g->mesh(i[0]+1,i[1]+1,i[2]+1)) >= 7 ) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2 ) )
		    return( false );
		i[2]++;
	    }
	} else {
	    throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	}
	
	// Check for collisions/mirroring with simulation boundary. Here
	// coordinates i are already advanced to next mesh.
	for( size_t a = 0; a < PP::dim(); a++ ) {

	    if( i[a] < 0 ) {
		if( _mirror[2*a] )
		    handle_mirror( c, i, a, -1, x2 );
		else {
		    handle_collision( particle, c, x2 );
		    return( false );
		}
	    } else if( i[a] >= (_pidata._g->size(a)-1) ) {
		if( _mirror[2*a+1] )
		    handle_mirror( c, i, a, +1, x2 );
		else {
		    handle_collision( particle, c, x2 );
		    return( false );
		}
	    }
	}

	return( true );
    }

    /*! \brief Limit trajectory advance from \a x1 to \a x2 to double
     *  the simulation box.
     *
     *  Return true if limitation is done, false if data is left untouched.
     */
    bool limit_trajectory_advance( const PP &x1, PP &x2 ) {

	bool touched = false;

	for( size_t a = 0; a < PP::dim(); a++ ) {

	    double lim1 = _pidata._g->origo(a) - 
		(_pidata._g->size(a)-1)*_pidata._g->h();
	    double lim2 = _pidata._g->origo(a) + 
		2*(_pidata._g->size(a)-1)*_pidata._g->h();

	    if( x2[2*a+1] < lim1 ) {
		
		double K = (lim1 - x1[2*a+1]) / (x2[2*a+1] - x1[2*a+1]);
		x2 = x1 + K*(x2-x1);
		touched = true;
#ifdef DEBUG_PARTICLE_ITERATOR
		std::cout << "Limiting step to:\n";
		std::cout << "  x2: " << x2 << "\n";
#endif
	    } else if(x2[2*a+1] > lim2 ) {

		double K = (lim2 - x1[2*a+1]) / (x2[2*a+1] - x1[2*a+1]);
		x2 = x1 + K*(x2-x1);
		touched = true;
#ifdef DEBUG_PARTICLE_ITERATOR
		std::cout << "Limiting step to:\n";
		std::cout << "  x2: " << x2 << "\n";
#endif
	    }
	}

	return( touched );
    }

    /*! \brief Handle particle iteration step from coordinates \a x1 to \a x2.
     *
     *  Searches mesh intersections between points \a x1 and \a x2 and
     *  builds ColData. Checks for collisions with solids and
     *  boundaries and sets space charge. If particle collides with
     *  mirroring boundary, the \a x2 coordinates are changed and GSL
     *  ODE integrator is resetted. Space charge is deposited at each
     *  mesh line crossing.
     *
     *  If \a force_linear is true, linear interpolation of trajectory
     *  is used regardless of interpolation settings.
     *
     *  Return true if particle status is PARTICLE_OK after trajectory
     *  step, false otherwise.
     */
    bool handle_trajectory( Particle<PP> &particle, const PP &x1, PP &x2, 
			    bool force_linear=false ) {

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Handle trajectory from x1 to x2:\n";
	std::cout << "  x1: " << x1 << "\n";
	std::cout << "  x2: " << x2 << "\n";
#endif

	// Limit trajectory advance to double the simulation box
	// If limitation done, force to linear interpolation
	if( limit_trajectory_advance( x1, x2 ) )
	    force_linear = true;

	// Make coldata
	if( _polyint && !force_linear )
	    ColData<PP>::build_coldata_poly( _coldata, *_pidata._g, x1, x2 );
	else
	    ColData<PP>::build_coldata_linear( _coldata, *_pidata._g, x1, x2 );

	// No intersections, nothing to do
	if( _coldata.size() == 0 ) {
#ifdef DEBUG_PARTICLE_ITERATOR
	    std::cout << "No coldata\n";
#endif
	    return( true );
	}

	// Starting mesh index
	int i[3] = {0, 0, 0};
	for( size_t cdir = 0; cdir < PP::dim(); cdir++ )
	    i[cdir] = (int)floor( (x1[2*cdir+1]-_pidata._g->origo(cdir))/_pidata._g->h() );

	// Process intersection points
#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Process coldata points:\n";
#endif
	for( size_t a = 0; a < _coldata.size(); a++ ) {

#ifdef DEBUG_PARTICLE_ITERATOR
	    std::cout << "  Coldata " << std::setw(4) << a << ": " 
		      << _coldata[a]._x << ", " 
		      << std::setw(3) << i[0] << " "
		      << std::setw(3) << i[1] << " "
		      << std::setw(3) << i[2] << " "
	    	      << std::setw(3) << _coldata[a]._dir << "\n";
#endif

	    // Advance particle in mesh, check for possible collisions and
	    // do mirroring.
	    handle_trajectory_advance( particle, a, i, x2 );

	    // Update space charge for one mesh.
	    if( _pidata._scharge )
		scharge_add_from_trajectory( *_pidata._scharge, particle.IQ(), 
					     _xi, _coldata[a]._x );

#ifdef DEBUG_PARTICLE_ITERATOR
	    if( particle.get_status() == PARTICLE_OUT ) {
		std::cout << "  Particle out\n";
		std::cout << "  x = " << x2 << "\n";
	    } else if( particle.get_status() == PARTICLE_COLL ) {
		std::cout << "  Particle collided\n";
		std::cout << "  x = " << x2 << "\n";
	    }
#endif
	    // Clear coldata and exit if particle collided.
	    if( particle.get_status() != PARTICLE_OK ) {
		_coldata.clear();
		return( false );
	    }

	    // Update last intersection point xi.
	    _xi = _coldata[a]._x;
	}

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Coldata done\n";
#endif
	_coldata.clear();
	return( true );
    }


    /*! \brief Is particle mirroring required at axis in cylindrical symmetry?
     */
     bool axis_mirror_required( const PP &x2 ) {
	 return( _pidata._g->geom_mode() == MODE_CYL && 
		 x2[4] < 0.0 && 
		 x2[3] <= 0.01*_pidata._g->h() &&
		 x2[3]*fabs(x2[5]) <= 1.0e-9*fabs(x2[4]) );
		 
     }


    /*! \brief Handle particle mirrored at axis in cylindrical symmetry.
     *
     *  Return true if particle status is PARTICLE_OK after trajectory
     *  step, false otherwise.
     */
    bool handle_axis_mirror_step( Particle<PP> &particle, const PP &x1, PP &x2 ) {

	// Get acceleration at x2
	double dxdt[5];
	PP::get_derivatives( x2[0], &x2[1], dxdt, (void *)&_pidata );

	// Calculate crossover point assuming zero acceleration in
	// r-direction and constant acceleration in x-direction
	double dt = -x2[3]/x2[4];
	PP xc;
	xc[0] = x2[0]+dt;
	xc[1] = x2[1]+(x2[2]+0.5*dxdt[1]*dt)*dt;
	xc[2] = x2[2];
	xc[3] = x2[3]+x2[4]*dt;
	xc[4] = x2[4];
	xc[5] = x2[5];

	// Mirror x2 to x3
	PP x3 = 2*xc - x2;
	x3[3] *= -1.0;
	x3[4] *= -1.0;
	x3[5] *= -1.0;

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Particle mirror:\n";
	std::cout << "  x1: " << x1 << "\n";
	std::cout << "  x2: " << x2 << "\n";
	std::cout << "  xc: " << xc << "\n";
	std::cout << "  x3: " << x3 << "\n";
#endif
	
	// Handle step with linear interpolation to avoid going to r<=0
	if( !handle_trajectory( particle, x2, x3, true ) )
	    return( false ); // Particle done

	// Save trajectory calculation points
	_traj.push_back( x2 );
	_traj.push_back( xc );
	xc[4] *= -1.0;
	xc[5] *= -1.0;
	_traj.push_back( xc );
	
	// Next step not a continuation of previous one, reset
	// integrator
	gsl_odeiv_step_reset( _step );
	gsl_odeiv_evolve_reset( _evolve );

	// Continue iteration at mirrored point
	x2 = x3;
	return( true );
    }
    
    /*! \brief Check particle definition.
     *
     *  False is returned if particle is defined out of simulation
     *  area or inside a solid.  False is also returned if a particle
     *  is defined on the simulation border going out of simulation
     *  with the corresponding border not having mirroring enabled. If
     *  the mirroring is enabled, the particle is mirrored and true is
     *  returned. Otherwise the particle definition is ok and true is
     *  returned.
     */
    bool check_particle_definition( PP &x ) {

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Particle defined at:\n";
	std::cout << "  x = " << x << "\n";
#endif

	// Check if inside solids of outside geometry.
	if( _pidata._g->inside( x.location() ) )
	    return( false );

	// Check if particle on simulation geometry border and directed outwards
	/*
	for( size_t a = 0; a < PP::dim(); a++ ) {
	    if( x[2*a+1] == _pidata._g->origo(a) && x[2*a+2] < 0.0 ) {
		if( _mirror[2*a] ) {
		    x[2*a+2] *= -1.0;
#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Mirroring to:\n";
	std::cout << "  x = " << x << "\n";
#endif
		} else {
		    return( false );
		}

	    } else if( x[2*a+1] == _pidata._g->max(a) & x[2*a+2] > 0.0 ) {
		if( _mirror[2*a+1] ) {
		    x[2*a+2] *= -1.0;
#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Mirroring to:\n";
	std::cout << "  x = " << x << "\n";
#endif
		} else {
		    return( false );
		}
	    }
	}

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "Definition ok\n";
#endif

	*/
	return( true );
    }
    
    double calculate_dt( const PP &x, const double *dxdt ) {

	double spd = 0.0, acc = 0.0;

	for( size_t a = 0; a < (PP::size()-1)/2; a++ ) {
	    //std::cout << "spd += " << dxdt[2*a]*dxdt[2*a] << "\n";
	    spd += dxdt[2*a]*dxdt[2*a];
	    //std::cout << "acc += " << dxdt[2*a+1]*dxdt[2*a+1] << "\n";
	    acc += dxdt[2*a+1]*dxdt[2*a+1];
	}
	if( _pidata._g->geom_mode() == MODE_CYL ) {
	    //std::cout << "MODE_CYL\n";
	    //std::cout << "spd += " << x[3]*x[3]*x[5]*x[5] << "\n";
	    spd += x[3]*x[3]*x[5]*x[5];
	    //std::cout << "acc += " << x[3]*x[3]*dxdt[4]*dxdt[4] << "\n";
	    acc += x[3]*x[3]*dxdt[4]*dxdt[4];
	}
	//std::cout << "spd = " << sqrt(spd) << "\n";
	//std::cout << "acc = " << sqrt(acc) << "\n";
	spd = _pidata._g->h() / sqrt(spd);
	acc = sqrt( 2.0*_pidata._g->h() / sqrt(acc) );

	return( spd < acc ? spd : acc );
    }

public:

    /*! \brief Constructor for new particle iterator.
     *
     *  New particle iterator is initialized with given settings.
     *
     *  \param type Particle iterator type used
     *  \param epsabs Absolute error limit in iteration
     *  \param epsrel Relative error limit in iteration
     *  \param polyint Interpolation type to use. True means use polynomial
     *  interpolation, false means use linear interpolation
     *  \param maxsteps Maximum number of steps to take before particle is killed
     *  \param maxt Maximum flight time for a particle
     *  \param trajdiv Trajectory saving divisor. Only every trajdiv:th particle 
     *  trajectory saved.
     *  \param mirror %Particle mirroring on surfaces
     *  \param scharge Space charge field to save to
     *  \param efield Electric field in the geometry
     *  \param bfield Magnetic field in the geometry
     *  \param g %Geometry definition
     *  \param first Pointer to first particle of the database.
     *
     *  The particle iterator is given the settings for calculation
     *  and geometry, electric field and space charge map to
     *  build. Pointer to first particle in the particle database
     *  vector is used to calculate the particle number from the
     *  particle memory location.
     */
    ParticleIterator( particle_iterator_type_e type, double epsabs, double epsrel, 
		      bool polyint, uint32_t maxsteps, double maxt, 
		      uint32_t trajdiv, bool mirror[6], ScalarField *scharge, 
		      const Efield *efield, const VectorField *bfield, 
		      const Geometry *g, Particle<PP> *first )
	: _type(type), _polyint(polyint), _epsabs(epsabs), _epsrel(epsrel), _maxsteps(maxsteps), _maxt(maxt), 
	  _trajdiv(trajdiv), _first(first), _pidata(scharge,efield,bfield,g), _end_time(0), 
	  _end_step(0), _end_out(0), _end_coll(0), _end_baddef(0), _sum_steps(0) {

	// Initialize mirroring
	_mirror[0] = mirror[0];
	_mirror[1] = mirror[1];
	_mirror[2] = mirror[2];
	_mirror[3] = mirror[3];
	_mirror[4] = mirror[4];
	_mirror[5] = mirror[5];

	// Initialize system of ordinary differential equations (ODE)
	_system.jacobian  = NULL;
	_system.params    = (void *)&_pidata;
	_system.function  = PP::get_derivatives;
	_system.dimension = PP::size()-1; // Time is not part of differential equation dimensions

	// Make scale
	// 2D:  x vx y vy
	// Cyl: x vx r vr omega
	// 3D:  x vx y vy z vz
	double scale_abs[PP::size()-1];
	for( uint32_t a = 0; a < (uint32_t)PP::size()-2; a+=2 ) {
	    scale_abs[a+0] = 1.0;
	    scale_abs[a+1] = 1.0e6;
	}
	if( _pidata._g->geom_mode() == MODE_CYL )
	    scale_abs[4] = 1.0;

	// Initialize ODE solver
	_step    = gsl_odeiv_step_alloc( gsl_odeiv_step_rkck, _system.dimension );
	//_control = gsl_odeiv_control_standard_new( _epsabs, _epsrel, 1.0, 1.0 );
	_control = gsl_odeiv_control_scaled_new( _epsabs, _epsrel, 1.0, 1.0, scale_abs, PP::size()-1 );
	_evolve  = gsl_odeiv_evolve_alloc( _system.dimension );
    }


    /*! \brief Destructor.
     */
    ~ParticleIterator() {
	gsl_odeiv_evolve_free( _evolve );
	gsl_odeiv_control_free( _control );
	gsl_odeiv_step_free( _step );
    }


    /*! \brief Enable plasma threshold model for NSIMP plasma model.
     */
    void enable_nsimp_plasma_threshold( const ScalarField *epot, double phi_plasma ) {
	_pidata._epot = epot;
	_pidata._phi_plasma = phi_plasma;
    }


    /*! \brief Get particle iterator statistics.
     *
     *  Returns statistics in variables \a end_time, \a end_step, \a
     *  end_out, \a end_coll, \a end_baddef and \a sum_steps. See
     *  ParticleDatabase for more details about the particle iterator
     *  statistics.
     */
    void get_statistics( uint32_t &end_time, uint32_t &end_step, uint32_t &end_out,
			 uint32_t &end_coll, uint32_t &end_baddef, uint32_t &sum_steps ) {
	end_time   = _end_time;
	end_step   = _end_step;
	end_out    = _end_out;
	end_coll   = _end_coll;
	end_baddef = _end_baddef;
	sum_steps  = _sum_steps;
    }

    /*! \brief Iterate a particle from start to end.
     *
     *  Iterate particle \a particle from start to end. This function
     *  is called by the Scheduler \a scheduler, which provides
     *  particles to be solved. Reference to \a scheduler is provided
     *  for the possibility to add secondary particles to particle
     *  database.
     */
    void operator()( Particle<PP> *particle,
		     Scheduler<ParticleIterator<PP>,Particle<PP>,Error> &scheduler ) {

	// Copy starting point to x and 
	PP x = particle->x();

	// Check particle definition
	if( !check_particle_definition( x ) ) {
	    particle->set_status( PARTICLE_BADDEF );
	    _end_baddef++;
	    return;
	}
	particle->x() = x;

	// Reset trajectory and save first trajectory point.
	_traj.clear();
	_traj.push_back( x );
#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << x[0] << " " 
		  << x[1] << " " 
		  << x[2] << " " 
		  << x[3] << " " 
		  << x[4] << "\n";
#endif
	_pidata._qm = particle->qm();
	_xi = x;

	// Reset integrator
	gsl_odeiv_step_reset( _step );
	gsl_odeiv_evolve_reset( _evolve );
	
	// Make initial guess for step size
	//std::cout << "Guess dt ------------------------------------------------\n";
	double dxdt[PP::size()-1];
	PP::get_derivatives( 0.0, &x[1], dxdt, (void *)&_pidata );
	double dt = calculate_dt( x, dxdt );

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "dxdt = ";
	for( size_t a = 0; a < PP::size()-1; a++ )
	    std::cout  << dxdt[a] << " ";
	std::cout << "\n";
	std::cout << "dt = " << dt << "\n";
	std::cout << "*** Starting iteration\n";
#endif
	
	// Iterate ODEs until maximum steps are done, time is used 
	// or particle collides.
	PP x2;
	size_t nstp = 0;
	while( nstp < _maxsteps && x[0] < _maxt ) {

#ifdef DEBUG_PARTICLE_ITERATOR
	    std::cout << "\n*** Step ***\n";
	    std::cout << "  x  = " << x2 << "\n";
	    std::cout << "  dt = " << dt << " (proposed)\n";
#endif
	    
	    // Take a step.
	    x2 = x;

	    while( nstp < _maxsteps ) {
		int retval = gsl_odeiv_evolve_apply( _evolve, _control, _step, &_system, 
						     &x2[0], _maxt, &dt, &x2[1] );
		if( retval == IBSIMU_DERIV_ERROR ) {
#ifdef DEBUG_PARTICLE_ITERATOR
		    std::cout << "Step rejected\n";
		    std::cout << "  x2 = " << x2 << "\n";
		    std::cout << "  dt = " << dt << "\n";
#endif
		    x2[0] = x[0]; // Reset time (this shouldn't be necessary - there 
		                  // is a bug in GSL-1.12, report has been sent)
		    dt *= 0.5;
		    nstp++;
		    continue;
		} else if( retval == GSL_SUCCESS ) {
		    break;
		} else {
		    throw( Error( ERROR_LOCATION, "gsl_odeiv_evolve_apply failed" ) );
		}
	    }
	    
	    // Check step count number and step size validity
	    if( nstp >= _maxsteps )
		break;
	    if( x2[0] == x[0] )
		throw( Error( ERROR_LOCATION, "too small step size" ) );

#ifdef DEBUG_PARTICLE_ITERATOR
	    std::cout << "Step accepted from x1 to x2:\n";
	    std::cout << "  dt = " << dt << " (taken)\n";
	    std::cout << "  x1 = " << x << "\n";
	    std::cout << "  x2 = " << x2 << "\n";
#endif

	    // Handle collisions and space charge of step.
	    if( !handle_trajectory( *particle, x, x2 ) ) {
		x = x2;
		break; // Particle done
	    }

	    // Check if particle mirroring is required to avoid 
	    // singularity at symmetry axis.
	    if( axis_mirror_required( x2 ) ) {
		if( !handle_axis_mirror_step( *particle, x, x2 ) )
		    break; // Particle done
	    }

	    // Propagate coordinates
	    x = x2;

	    // Save trajectory point
	    _traj.push_back( x2 );
	    
	    // Increase step count.
	    nstp++;
	}

#ifdef DEBUG_PARTICLE_ITERATOR
	std::cout << "\n*** Done stepping ***\n";
#endif

	// Check if step count or time limited 
	if( nstp == _maxsteps ) {
	    particle->set_status( PARTICLE_NSTP );
	    _end_step++;
	} else if( x[0] >= _maxt ) {
	    particle->set_status( PARTICLE_TIME );
	    _end_time++;
	}

	// Save step count
	_sum_steps += nstp;

	// Save trajectory of current particle
	if( _trajdiv != 0 && (particle-_first) % _trajdiv == 0 )
	    particle->copy_trajectory( _traj );

	// Save last particle location
	particle->x() = x;
    }

/*
    std::cout << "Kala\n";

    std::cout << _coldata.capacity() << "\n";
    std::cout << _traj.capacity() << "\n";
*/
};


#endif
