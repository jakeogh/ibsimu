/*! \file epot_problem.cpp
 *  \brief Source code for epot_problem.cpp
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

#define SPM_RANGE_CHECK 1
#include "mvector.hpp"
#include "crowmatrix.hpp"
#include "epot_problem.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


// From NIST, official value of Apr 7 2008
#define EPSILON0 8.854187817e-12


/* ************************************** *
 * Node2DoF                               *
 * ************************************** */

void EpotProblem::Node2DoF::debug_print( void ) const
{
    int a;
    std::cout << "**Node2DoF\n";
    std::cout << "size = ("
       	<< _size[0] << ", "
	<< _size[1] << ", "
	<< _size[2] << ")\n";
    std::cout << "n2d = (";
    for( a = 0; a < _size[0]*_size[1]*_size[2]-1; a++ )
	std::cout << _n2d[a] << ", ";
    if( a < _size[0]*_size[1]*_size[2] )
	std::cout << _n2d[a] << ")\n";
}


/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */


EpotProblem::EpotProblem() 
    : _nodecount(0), _dof(0), _fd_mat(0), _fd_vec(0), _neumann_order(2), 
      _smooth_solid(true), _plasma(PLASMA_NONE), _rhoe(0.0), _Te(0.0), _Up(0.0), 
      _force_pot(0.0), _force_pot_func(0), _init_plasma_func(0), _solver(0)
{
    
}


EpotProblem::EpotProblem( std::istream &s )
{
    
}


EpotProblem::~EpotProblem()
{
    if( _fd_mat )
	delete _fd_mat;
    if( _fd_vec )
	delete _fd_vec;
}


void EpotProblem::clear_problem( void )
{
    if( _solver )
	_solver->reset();
    if( _fd_mat )
	delete _fd_mat;
    if( _fd_vec )
	delete _fd_vec;
    _fd_mat = 0;
    _fd_vec = 0;
    _dof = 0;
}


/* ************************************** *
 * EpotProblem constructing                   *
 * ************************************** */


void EpotProblem::set_neumann_order( int32_t order )
{
    if( order < 1 || order > 2 )
	throw( Error( ERROR_LOCATION, "illegal neumann order" ) );
    _neumann_order = order;
    clear_problem();
}


void EpotProblem::enable_smooth_solids( bool enable )
{
    _smooth_solid = enable;
    clear_problem();
}


void EpotProblem::set_forced_potential_volume( double force_pot, 
					       bool (*force_pot_func)(double,double,double) )
{
    _force_pot = force_pot;
    _force_pot_func = force_pot_func;
}


void EpotProblem::set_initial_plasma( double Up, 
				      bool (*plasma_func)(double,double,double) )
{
    _plasma     = PLASMA_INITIAL;
    _Up         = Up;
    if( !plasma_func )
	throw( Error( ERROR_LOCATION, "NULL initial plasma function" ) );
    _init_plasma_func = plasma_func;
    clear_problem();
}


void EpotProblem::set_pexp_plasma( double rhoe, double Te, double Up )
{
    _plasma     = PLASMA_PEXP;
    _rhoe       = -fabs(rhoe); // Ensure correct sign of charge density
    _Te         = Te;
    _Up         = Up;
    clear_problem();
}



void EpotProblem::set_nsimp_initial_plasma( bool (*plasma_func)(double,double,double) )
{
    _plasma     = PLASMA_INITIAL;
    _Up         = 0.0;
    if( !plasma_func )
	throw( Error( ERROR_LOCATION, "NULL initial plasma function" ) );
    _init_plasma_func = plasma_func;
    clear_problem();
}


void EpotProblem::set_nsimp_plasma( double rhop, double Ep, 
				    std::vector<double> rhoi, std::vector<double> Ei )
{
    _plasma = PLASMA_NSIMP;
    _rhoi.clear();
    _Ei.clear();

    _rhoi.push_back( fabs(rhop) ); // Ensure correct sign of charge density
    _Ei.push_back( Ep );

    size_t size = rhoi.size() < Ei.size() ? rhoi.size() : Ei.size();
    for( size_t a = 0; a < size; a++ ) {
	_rhoi.push_back( fabs(rhoi[a]) ); // Ensure correct sign of charge density
	_Ei.push_back( Ei[a] );
    }

    clear_problem();
}


/*! \brief Set a link in the linear system A*X=B.
 * 
 * Makes the node a dependent on node b. If b is positive this means
 * that element (a,b) of matrix A is set to value val. If b is
 * negative, the node has a fixed potential and the dependence should
 * be added to the vector side of the system of equations on row
 * a. The value of b indicates the boundary value number of the node.
 */
void EpotProblem::set_link( CRowMatrix &A, Vector &B, 
			    int32_t a, int32_t b, double val )
{
    if( b < 0 )
	B[a] += -val*_g->get_boundary(-b).val;
    else
	A.construct_add( a, b, val );
}


/*! \brief Adds an initial plasma node (i,j,k) to the linear system.
 *
 *  Initial guess for plasma node is a fixed potential value at plasma potential.
 */
void EpotProblem::add_initial_plasma( int32_t i, int32_t j, int32_t k, 
				      CRowMatrix &A, Vector &B, Node2DoF &n2d )
{
    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
    B[n2d(i,j,k)] += _Up;
}


/*! \brief Adds an forced potential node (i,j,k) to the linear system.
 *
 *  The node is forced to \a _force_pot.
 */
void EpotProblem::add_forced_pot( int32_t i, int32_t j, int32_t k, 
				  CRowMatrix &A, Vector &B, Node2DoF &n2d )
{
    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
    B[n2d(i,j,k)] += _force_pot;
}


/*! \brief Adds a vacuum node (i,j,k) to the linear system.
 */
void EpotProblem::add_vacuum_node( int32_t i, int32_t j, int32_t k, 
				   CRowMatrix &A, Vector &B, Node2DoF &n2d )
{
    switch( _g->geom_mode() ) {
    case MODE_1D:
	set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j,k), -2.0 );
	break;
    case MODE_2D:
	set_link( A, B, n2d(i,j,k), n2d(i,j-1,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j,k), -4.0 );
	set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j+1,k), 1.0 );
	break;
    case MODE_CYL:
	set_link( A, B, n2d(i,j,k), n2d(i,j-1,k), 1.0-0.5/j );
	set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j,k), -4.0 );
	set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j+1,k), 1.0+0.5/j );
	break;
    default:
	set_link( A, B, n2d(i,j,k), n2d(i,j,k-1), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j-1,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j,k), -6.0 );
	set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j+1,k), 1.0 );
	set_link( A, B, n2d(i,j,k), n2d(i,j,k+1), 1.0 );
	break;
    }

    // Space charge is added when problem is being solved
}


/* Adds a Neumann boundary node (i,j,k) to the linear system.
 *
 * Second order approximation is used if enabled and the approximation
 * does not reference any solid nodes (edges are ok if smooth edges enabled)
 */
void EpotProblem::add_neumann_node( signed char a, int32_t i, int32_t j, int32_t k, 
				    CRowMatrix &A, Vector &B, Node2DoF &n2d )
{
    double coef = 0.0;



    switch( a ) {
    case -1:
	if( _neumann_order == 2 && 
	    ( (_smooth_solid && 
	       _g->mesh_check(i+1,j,k) <= 0 && _g->mesh_check(i+2,j,k) <= 0) || 
	      (_g->mesh_check(i+1,j,k) == 0 && _g->mesh_check(i+2,j,k) == 0) ) ) {
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 3.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), -4.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i+2,j,k), 1.0 );
	    coef = -2.0;
	} else {
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), -1.0 );
	    coef = -1.0;
	}
	break;
    case -2:
	if( _neumann_order == 2 &&
	    ( (_smooth_solid && 
	       _g->mesh_check(i-1,j,k) <= 0 && _g->mesh_check(i-2,j,k) <= 0) || 
	      (_g->mesh_check(i-1,j,k) == 0 && _g->mesh_check(i-2,j,k) == 0) ) ) {
	    set_link( A, B, n2d(i,j,k), n2d(i-2,j,k), 1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), -4.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 3.0 );
	    coef = -2.0;
	} else {
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), -1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
	    coef = -1.0;
	}
	break;
    case -3:
	if( _g->geom_mode() == MODE_CYL ) {
	    // Axis boundary, exception, not really neumann
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), 1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), -6.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), 1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j+1,k), 4.0 );
	    // Space charge is added when problem is being solved
	    return;
	} else {
	    if( _neumann_order == 2 &&
		( (_smooth_solid && 
		   _g->mesh_check(i,j+1,k) <= 0 && _g->mesh_check(i,j+2,k) <= 0) || 
		  (_g->mesh_check(i,j+1,k) == 0 && _g->mesh_check(i,j+2,k) == 0) ) ) {
		set_link( A, B, n2d(i,j,k), n2d(i,j,k), 3.0 );
		set_link( A, B, n2d(i,j,k), n2d(i,j+1,k), -4.0 );
		set_link( A, B, n2d(i,j,k), n2d(i,j+2,k), 1.0 );
		coef = -2.0;
	    } else {
		set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
		set_link( A, B, n2d(i,j,k), n2d(i,j+1,k), -1.0 );
		coef = -1.0;
	    }
	}
	break;
    case -4:
	if( _neumann_order == 2 &&
	    ( (_smooth_solid && 
	       _g->mesh_check(i,j-1,k) <= 0 && _g->mesh_check(i,j-2,k) <= 0) || 
	      (_g->mesh_check(i,j-1,k) == 0 && _g->mesh_check(i,j-2,k) == 0) ) ) {
	    set_link( A, B, n2d(i,j,k), n2d(i,j-2,k), 1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j-1,k), -4.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 3.0 );
	    coef = -2.0;
	} else {
	    set_link( A, B, n2d(i,j,k), n2d(i,j-1,k), -1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
	    coef = -1.0;
	}
	break;
    case -5:
	if( _neumann_order == 2 &&
	    ( (_smooth_solid && 
	       _g->mesh_check(i,j,k+1) <= 0 && _g->mesh_check(i,j,k+2) <= 0) || 
	      (_g->mesh_check(i,j,k+1) == 0 && _g->mesh_check(i,j,k+2) == 0) ) ) {
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 3.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k+1), -4.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k+2), 1.0 );
	    coef = -2.0;
	} else {
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k+1), -1.0 );
	    coef = -1.0;
	}
	break;
    default:
	if( _neumann_order == 2 &&
	    ( (_smooth_solid && 
	       _g->mesh_check(i,j,k-1) <= 0 && _g->mesh_check(i,j,k-2) <= 0) || 
	      (_g->mesh_check(i,j,k-1) == 0 && _g->mesh_check(i,j,k-2) == 0) ) ) {
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k-2), 1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k-1), -4.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 3.0 );
	    coef = -2.0;
	} else {
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k-1), -1.0 );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
	    coef = -1.0;
	}
	break;
    }

    // Add right hand side coefficient
    B[n2d(i,j,k)] = coef * _g->h() * _g->get_boundary(-_g->mesh(i,j,k)).val;
}


bool vac_or_edg( const Geometry *g, int32_t i, int32_t j, int32_t k ) 
{
    signed char a = g->mesh_check( i, j, k ); 
    
    return( a >= -6 && a <= 6 ); // Vacuum or Neumann or Dirichlet
}


bool sol_or_dir( const Geometry *g, int32_t i, int32_t j, int32_t k ) 
{
    signed char a = g->mesh_check( i, j, k ); 
    
    return( a > 0 || a <= -7 ); // Solid or Dirichlet (or outside simulation geometry)
}


/*! \brief Adds a solid edge node (i,j,k) to the linear system.
 */
void EpotProblem::add_solid_edge_node( signed char a, int32_t i, int32_t j, int32_t k, 
				       CRowMatrix &A, Vector &B, Node2DoF &n2d )
{
    Vec3D xsurf;

    switch( _g->geom_mode() ) {
    case MODE_1D:

	// Edge (i+1)
	if( vac_or_edg(_g,i+1,j,k) &&
	    sol_or_dir(_g,i-1,j,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k), 
					    _g->origo()+_g->h()*Vec3D(i+1,j,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0-x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Edge (i-1)
	else if( vac_or_edg(_g,i-1,j,k) &&
		 sol_or_dir(_g,i+1,j,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k), 
					    _g->origo()+_g->h()*Vec3D(i-1,j,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), x );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0-x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Otherwise like normal solid node
	else {
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}
	break;

    case MODE_2D:
    case MODE_CYL:  // *******************************************************************

	// Face (i+1,j)
	if( vac_or_edg(_g,i+1,j,  k) && 
	    sol_or_dir(_g,i-1,j,  k) &&
	    sol_or_dir(_g,i,  j+1,k) &&
	    sol_or_dir(_g,i,  j-1,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0-x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	} 

	// Face (i-1,j)
	else if( sol_or_dir(_g,i+1,j,  k) &&
		 vac_or_edg(_g,i-1,j,  k) &&
		 sol_or_dir(_g,i,  j+1,k) &&
		 sol_or_dir(_g,i,  j-1,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), x );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0-x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}
	
	// Face (i,j+1)
	else if( sol_or_dir(_g,i+1,j,  k) &&
		 sol_or_dir(_g,i-1,j,  k) &&
		 vac_or_edg(_g,i,  j+1,k) &&
		 sol_or_dir(_g,i,  j-1,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j+1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0-x );
	    set_link( A, B, n2d(i,j,k), n2d(i,j+1,k), x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}
	
	// Face (i,j-1) 
	else if( sol_or_dir(_g,i+1,j,  k) &&
		 sol_or_dir(_g,i-1,j,  k) &&
		 sol_or_dir(_g,i,  j+1,k) &&
		 vac_or_edg(_g,i,  j-1,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j-1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,j-1,k), x );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0-x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i+1,j+1)
	else if( vac_or_edg(_g,i+1,j,  k) &&
		 sol_or_dir(_g,i-1,j,  k) &&
		 vac_or_edg(_g,i,  j+1,k) &&
		 sol_or_dir(_g,i,  j-1,k) &&
		 vac_or_edg(_g,i+1,j+1,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j+1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j+1,k), x*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	} 

	// Corner (i+1,j-1)
	else if( vac_or_edg(_g,i+1,j,  k) &&
		 sol_or_dir(_g,i-1,j,  k) &&
		 sol_or_dir(_g,i,  j+1,k) &&
		 vac_or_edg(_g,i,  j-1,k) &&
		 vac_or_edg(_g,i+1,j-1,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j-1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j-1,k), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k), (1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	} 

	// Corner (i-1,j-1)
	else if( sol_or_dir(_g,i+1,j,  k) &&
		 vac_or_edg(_g,i-1,j,  k) &&
		 sol_or_dir(_g,i,  j+1,k) &&
		 vac_or_edg(_g,i,  j-1,k) &&
		 vac_or_edg(_g,i-1,j-1,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j-1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j-1,k), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k), (1.0-x)*(1.0-x) );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	} 

	// Corner (i-1,j+1)
	else if( sol_or_dir(_g,i+1,j,  k) &&
		 vac_or_edg(_g,i-1,j,  k) &&
		 vac_or_edg(_g,i,  j+1,k) &&
		 sol_or_dir(_g,i,  j-1,k) &&
		 vac_or_edg(_g,i-1,j+1,k) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j+1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j+1,k), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k), (1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Otherwise like normal solid node
	else {

	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	break;

    case MODE_3D: // *******************************************************************

	// Face (i+1,j,k)
	if( vac_or_edg(_g,i+1,j,  k  ) && 
	    sol_or_dir(_g,i-1,j,  k  ) &&
	    sol_or_dir(_g,i,  j+1,k  ) &&
	    sol_or_dir(_g,i,  j-1,k  ) &&
	    sol_or_dir(_g,i,  j,  k+1) &&
	    sol_or_dir(_g,i,  j,  k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,k), 1.0-x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,k), x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Face (i-1,j,k)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) && 
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,k), x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,k), 1.0-x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Face (i,j+1,k)
	else if( sol_or_dir(_g,i+1,j,  k  ) && 
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j+1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,  k), 1.0-x );
	    set_link( A, B, n2d(i,j,k), n2d(i,j+1,k), x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Face (i,j-1,k)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j-1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,  k), 1.0-x );
	    set_link( A, B, n2d(i,j,k), n2d(i,j-1,k), x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Face (i,j,k+1)
	else if( sol_or_dir(_g,i+1,j,  k  ) && 
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0-x );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k+1), x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Face (i,j,k-1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k-1), x );
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0-x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// XY-Edge (i+1,j+1,k)
	else if( vac_or_edg(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i+1,j+1,k  ) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j+1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j+1,k), x*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// XY-Edge (i-1,j+1,k)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i-1,j+1,k  ) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j+1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j+1,k), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k), (1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// XY-Edge (i+1,j-1,k)
	else if( vac_or_edg(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i+1,j-1,k  ) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j-1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j-1,k), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k), (1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// XY-Edge (i-1,j-1,k)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i-1,j-1,k  ) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j-1,k), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j-1,k), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k), (1.0-x)*(1.0-x) );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// XZ-Edge (i+1,j,k+1)
	else if( vac_or_edg(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i+1,j,  k+1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k  ), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k+1), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k+1), x*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// XZ-Edge (i-1,j,k+1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i-1,j,  k+1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k  ), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k+1), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k+1), (1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// XZ-Edge (i+1,j,k-1)
	else if( vac_or_edg(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i+1,j,  k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k-1), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k-1), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k  ), (1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// XZ-Edge (i-1,j,k-1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i-1,j,  k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k-1), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k-1), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k  ), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x) );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// YZ-Edge (i,j+1,k+1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i,  j+1,k+1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j+1,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k  ), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k+1), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k+1), x*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// YZ-Edge (i,j-1,k+1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i,  j-1,k+1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j-1,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k  ), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k+1), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k+1), (1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// YZ-Edge (i,j+1,k-1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i,  j+1,k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j+1,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k-1), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k-1), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k  ), (1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// YZ-Edge (i,j-1,k-1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i,  j-1,k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i,j-1,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k-1), x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k-1), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k  ), (1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x) );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i+1,j+1,k+1)
	else if( vac_or_edg(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i+1,j+1,k+1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j+1,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j+1,k  ), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k+1), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k+1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k+1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j+1,k+1), x*x*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i-1,j+1,k+1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i-1,j+1,k+1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j+1,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j+1,k  ), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k+1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k+1), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j+1,k+1), x*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k+1), (1.0-x)*x*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i+1,j-1,k+1)
	else if( vac_or_edg(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i+1,j-1,k+1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j-1,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j-1,k  ), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k+1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j-1,k+1), x*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k+1), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k+1), (1.0-x)*x*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i-1,j-1,k+1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 vac_or_edg(_g,i,  j,  k+1) &&
		 sol_or_dir(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i-1,j-1,k+1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j-1,k+1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j-1,k  ), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j-1,k+1), x*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k+1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k+1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k+1), (1.0-x)*(1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i+1,j+1,k-1)
	else if( vac_or_edg(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i+1,j+1,k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j+1,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k-1), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k-1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k-1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j+1,k-1), x*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j+1,k  ), (1.0-x)*x*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i-1,j+1,k-1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) &&
		 vac_or_edg(_g,i,  j+1,k  ) &&
		 sol_or_dir(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i-1,j+1,k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j+1,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k-1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k-1), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j+1,k-1), x*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k-1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j+1,k  ), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j+1,k  ), (1.0-x)*(1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i+1,j-1,k-1)
	else if( vac_or_edg(_g,i+1,j,  k  ) &&
		 sol_or_dir(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i+1,j-1,k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i+1,j-1,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k-1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j-1,k-1), x*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k-1), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k-1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j-1,k  ), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x)*(1.0-x) );
	    set_link( A, B, n2d(i,j,k), n2d(i+1,j,  k  ), (1.0-x)*(1.0-x)*x );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Corner (i-1,j-1,k-1)
	else if( sol_or_dir(_g,i+1,j,  k  ) &&
		 vac_or_edg(_g,i-1,j,  k  ) &&
		 sol_or_dir(_g,i,  j+1,k  ) &&
		 vac_or_edg(_g,i,  j-1,k  ) &&
		 sol_or_dir(_g,i,  j,  k+1) &&
		 vac_or_edg(_g,i,  j,  k-1) &&
		 vac_or_edg(_g,i-1,j-1,k-1) ) {
	    double x = _g->bracket_surface( -a, _g->origo()+_g->h()*Vec3D(i,j,k),
					    _g->origo()+_g->h()*Vec3D(i-1,j-1,k-1), xsurf );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j-1,k-1), x*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k-1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k-1), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k-1), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j-1,k  ), (1.0-x)*x*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j-1,k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i-1,j,  k  ), (1.0-x)*(1.0-x)*x );
	    set_link( A, B, n2d(i,j,k), n2d(i,  j,  k  ), (1.0-x)*(1.0-x)*(1.0-x) );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}

	// Otherwise like normal solid node
	else {
	    set_link( A, B, n2d(i,j,k), n2d(i,j,k), 1.0 );
	    B[n2d(i,j,k)] += _g->get_boundary(-a).val;
	}
	break;

    default:

	throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	break;
    }
}


void EpotProblem::construct( const Geometry &g )
{
    signed char a;
    int32_t i, j, k;
    double x, y, z;

    if( ibsimu.get_verbose_output() ) {
	if( _plasma == PLASMA_NONE )
	    std::cout << "Constructing linear electric potential problem\n";
	else if( _plasma == PLASMA_INITIAL )
	    std::cout << "Constructing linear starting point for plasma problem\n";
	else
	    std::cout << "Constructing nonlinear plasma problem\n";
    }

    if( !g.built() )
	throw( Error( ERROR_LOCATION, "geometry mesh not built" ) );

    _nodecount = g.nodecount();
    _g = &g;
    if( _plasma == PLASMA_INITIAL ) {
	if( ibsimu.get_verbose_output() )
	    std::cout << "  Using initial plasma volume at Up = " << _Up << " V\n";
    } else if( _plasma == PLASMA_PEXP ) {
	if( ibsimu.get_verbose_output() ) {
	    std::cout << "  Using exponential plasma model for positive ion extraction\n";
	    std::cout << "  Te = " << _Te << " eV, Up = " << _Up 
		      << " V, rhoe = " << _rhoe << " C/m^3\n";
	}
    } else if( _plasma == PLASMA_NSIMP ) {
	if( ibsimu.get_verbose_output() ) {
	    std::cout << "  Using negative ion extraction plasma model\n";
	    std::cout << "  Ep = " << _Ei[0] << " eV, rhop = " << _rhoi[0] << " C/m^3\n";
	    for( size_t a = 1; a < _Ei.size(); a++ )
		std::cout << "  Ei[" << a << "] = " << _Ei[a] << " eV, rhop[" 
			  << a << "] = " << _rhoi[a] << " C/m^3\n";
	}
    }

    // Delete old problem formulation
    clear_problem();

    // Build n2d array and calculate degrees of freedom
    _n2d.resize( _g->size() );
    _dof = 0;
    for( i = 0; i < _nodecount; i++ ) {
	a = _g->mesh(i);
	if( (_smooth_solid && a <= 0) || (!_smooth_solid && a <= 0 && a >= -6) ) { 
	    // Vacuum, solid edge (if smooth solid edges enabled) or Neumann boundary
	    _n2d(i) = _dof;
	    _dof++;
	} else { 
	    // Dirichlet boundary or solid interior
	    _n2d(i) = -abs(a);
	}
    }

    if( _dof == 0 ) {
	_fd_mat = 0;
	_fd_vec = 0;
	throw( Error( ERROR_LOCATION, "zero degrees of freedom" ) );
    }

    if( ibsimu.get_verbose_output() ) {
	std::cout << "  dof = " << _dof << "\n";
    }

    // Allocate problem matrix and vector
    _fd_mat = new CRowMatrix( _dof, _dof );
    _fd_vec = new Vector( _dof );

    // Build matrix and vector contents
    for( k = 0; k < _g->size(2); k++ ) {
	for( j = 0; j < _g->size(1); j++ ) {
	    for( i = 0; i < _g->size(0); i++ ) {

		a = _g->mesh(i,j,k);
		x = i*_g->h()+_g->origo(0);
		y = j*_g->h()+_g->origo(1);
		z = k*_g->h()+_g->origo(2);

		if( _smooth_solid && a < -6 ) { 
		    // Solid smooth edge
		    add_solid_edge_node( a, i, j, k, *_fd_mat, *_fd_vec, _n2d );
		} else if( _plasma == PLASMA_INITIAL && (a <= 0 && a >= -6) &&
			   _init_plasma_func( x, y, z ) ) {
		    // Initial plasma area (vacuum or neumann)
		    add_initial_plasma( i, j, k, *_fd_mat, *_fd_vec, _n2d );
		} else if( _force_pot_func != 0 && (a <= 0 && a >= -6) && 
			   _force_pot_func( x, y, z ) ) {
		    // Forced potential area (vacuum or neumann)
		    add_forced_pot( i, j, k, *_fd_mat, *_fd_vec, _n2d );
		} else if( a >= -6 && a < 0 ) {
		    // Neumann boundary
		    add_neumann_node( a, i, j, k, *_fd_mat, *_fd_vec, _n2d );
		} else if( a == 0 ) {
		    // Vacuum
		    add_vacuum_node( i, j, k, *_fd_mat, *_fd_vec, _n2d );
		}
	    }
	}
    }
}


void EpotProblem::set_solver( Solver &s )
{
    _solver = &s;
}


/* ************************************** *
 * Solver interface                       *
 * ************************************** */


void EpotProblem::solve( ScalarField &epot, const ScalarField &scharge ) const
{
    if( ibsimu.get_verbose_output() ) {
	std::cout << "Solving problem\n";
    }

    if( _solver == 0 )
	throw( Error( ERROR_LOCATION, "solver undefined" ) );
    if( _dof == 0 )
	throw( Error( ERROR_LOCATION, "problem undefined" ) );
    if( scharge.size() != _g->size() )
	throw( Error( ERROR_LOCATION, "incorrect size of space charge field" ) );

    // Resize epot if necessary
    if( epot.size() != _g->size() )
	epot.reset( _g->geom_mode(), _g->size(), _g->origo(), _g->h() );

    // Allocate temporary space for calculation
    Vector X( _dof );
    _fd_vec2 = new Vector( _dof );
    if( !linear() ) {
	_fd_vec3 = new Vector( _dof );
	_fd_mat2 = new CRowMatrix( _dof, _dof );
    }

    // 1. Load content from epot to X.  2. Construct rhs vector with
    // space charge in vacuum nodes (note the exception of cylindrical
    // axis, which is marked as Neumann boundary).
    int32_t a, b;
    double coef = _g->h()*_g->h()/EPSILON0;
    for( a = 0; a < _nodecount; a++ ) {
	if( (b = _n2d(a)) >= 0 ) {
	    X(b) = epot(a);
	    if( _g->mesh(a) == 0 || (_g->mesh(a) == -3 && _g->geom_mode() == MODE_CYL) ) {
		// Plain vacuum or cylindrical axis node
		(*_fd_vec2)(b) = (*_fd_vec)(b) - coef*scharge(a);
	    } else {
		// Neumann or edge node
		(*_fd_vec2)(b) = (*_fd_vec)(b);
	    }
	}
    }

    // Call solver
     _solver->solve( *this, X );

    // Load content from X to epot.
    for( a = 0; a < _nodecount; a++ ) {
	if( (b = _n2d(a)) < 0 )
	    epot(a) = _g->get_boundary(-b).val;
	else
	    epot(a) = X(b);
    }

    // Free working vectors and matrix
    delete _fd_vec2;
    _fd_vec2 = 0;
    if( !linear() ) {
	delete _fd_vec3;
	delete _fd_mat2;
	_fd_vec3 = 0;
	_fd_mat2 = 0;
    }
}


void EpotProblem::get_vecmat( const Matrix **A, const Vector **B ) const
{
    // fd_vec2 contains geometry defined right-hand-side plus
    // contribution from space charge.
    *A = _fd_mat;
    *B = _fd_vec2;
}


#define GSCOEF  1.12837916709551257390

void EpotProblem::get_resjac( const Matrix **J, const Vector **R, const Vector &X ) const
{
    // Precalculate coefficients
    double Q = 0.0, K = 0.0;
    if( _plasma == PLASMA_PEXP ) {
	Q = 1.0/_Te;
	K = _Up/_Te;
    }

    // Construct jacobian to _fd_mat2 from general (linear) problem matrix _fd_mat;
    // Calculate R = J0*X - B(X) and J = J0 + I*D(X)

    // Contribution from linear coefficients (linear Poisson)
    (*_fd_vec3) = (*_fd_mat) * X;
    (*_fd_mat2) = (*_fd_mat);

    Vector Diag( _dof );
    int32_t a, b;
    double coef = _g->h()*_g->h()/EPSILON0;

    for( int32_t k = 0; k < _g->size(2); k++ ) {
	for( int32_t j = 0; j < _g->size(1); j++ ) {
	    for( int32_t i = 0; i < _g->size(0); i++ ) {

		a = i + j*_g->size(0) + k*_g->size(0)*_g->size(1);
		//for( a = 0; a < _nodecount; a++ ) {

		// For each vacuum node calculate contribution of plasma model to matrix and vector
		if( (b = _n2d(a)) >= 0 ) {
		    if( _g->mesh(a) == 0 || (_g->mesh(a) == -3 && _g->geom_mode() == MODE_CYL) ) {
			
			if( _force_pot_func != 0 && (_g->mesh(a) <= 0 && _g->mesh(a) >= -6) && 
			    // Forced potential node
			    _force_pot_func(i*_g->h()+_g->origo(0),
					    j*_g->h()+_g->origo(1),
					    k*_g->h()+_g->origo(2)) ) {
			    (*_fd_vec3)(b) = 0.0;
			} else if( _plasma == PLASMA_PEXP ) {
			    double Y = exp( Q*X(b) - K );
			    (*_fd_vec3)(b) -= (*_fd_vec2)(b) - coef*_rhoe*Y;
			    _fd_mat2->set(b,b) += coef*_rhoe*Q*Y;
			} else if( _plasma == PLASMA_NSIMP ) {
			    double xx = X(b)/_Ei[0];
			    double f = _rhoi[0]*(1.0 + erf( -xx ));
			    double df = -_rhoi[0] * GSCOEF*exp(-xx*xx)/_Ei[0];
			    for( size_t a = 1; a < _Ei.size(); a++ ) {
				xx = exp( -X(b)/_Ei[a] );
				f += _rhoi[a]*xx;
				df -= _rhoi[a]*xx/_Ei[a];
			    }
			    (*_fd_vec3)(b) -= (*_fd_vec2)(b) - coef*f;
			    _fd_mat2->set(b,b) += coef*df;
			} else {
			    throw( Error( ERROR_LOCATION, "unknown plasma model type" ) );
			}
		    } else {
			// Neumann or edge node
			(*_fd_vec3)(b) -= (*_fd_vec2)(b);
		    }
		}
	
	    }
	}
    }
    
    *J = _fd_mat2;
    *R = _fd_vec3;
}


bool EpotProblem::linear( void ) const 
{
    if( _plasma == PLASMA_NONE || _plasma == PLASMA_INITIAL )
	return( true );
    else
	return( false );
}

/* ************************************** *
 * Misc                                   *
 * ************************************** */


void EpotProblem::debug_print( void ) const 
{
    std::cout << "**EpotProblem\n";
    std::cout << "dof = " << _dof << "\n";
    std::cout << "n2d:\n";
    _n2d.debug_print();
    std::cout << "fd_mat = \n" << *_fd_mat << "\n";
    std::cout << "fd_vec = \n" << *_fd_vec << "\n";
    std::cout << "neumann_order = " << _neumann_order << "\n";
    std::cout << "plasma = " << _plasma << "\n";
    std::cout << "rhoe = " << _rhoe << "\n";
    std::cout << "Tc = " << _Te << "\n";
    std::cout << "Up = " << _Up << "\n";
}


void EpotProblem::save( std::ostream &s ) const
{

}















