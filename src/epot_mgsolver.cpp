/*! \file epot_mgsolver.cpp
 *  \brief Multigrid solver for electric potential problem
 */

/* Copyright (c) 2011 Taneli Kalvas. All rights reserved.
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


#include "epot_mgsolver.hpp"
#include "ibsimu.hpp"
#include "constants.hpp"
#include "compmath.hpp"




/* *****************************************************************************
 * EpotMGSolver
 */


EpotMGSolver::EpotMGSolver( Geometry &geom )
    : EpotSolver( geom ), _geom_prepared(false), _levels(1), _npre(5), 
      _npost(5), _mgcycmax(100), _mgcyc(0), _mgeps(1.0e-6), _gamma(1), _res(0.0), 
      _eps(1.0e-12), _w(1.7), _imax(10000)
{

}


EpotMGSolver::EpotMGSolver( Geometry &geom, std::istream &s )
    : EpotSolver(geom,s)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


EpotMGSolver::~EpotMGSolver()
{
    reset_problem();
}


void EpotMGSolver::reset_problem( void )
{
    //std::cout << "Reset!\n";

    // Free geom and subsolver arrays made in prepare_mg_geom()
    for( uint32_t a = 0; a < _geomv.size(); a++ ) {
	if( a >= 1 )
	    delete _geomv[a];
	delete _epotsolverv[a];
    }

    // Free other arrays made in preprocess()
    for( uint32_t a = 0; a < _epotv.size(); a++ ) {
	delete _rhsv[a];
	delete _workv[a];
	if( _work2v[a] )
	    delete _work2v[a];
	if( a >= 1 )
	    delete _epotv[a];
    }

    _epotv.clear();
    _geomv.clear();
    _epotsolverv.clear();
    _rhsv.clear();
    _workv.clear();
    _work2v.clear();

    _res = 0.0;
    _geom_prepared = false;
}


void EpotMGSolver::set_eps( double eps )
{
    _eps = eps;
}

void EpotMGSolver::set_w( double w )
{
    _w = w;
}


void EpotMGSolver::set_imax( uint32_t imax )
{
    _imax = imax;
}

double EpotMGSolver::get_residual( void ) const
{
    return( _res );
}


uint32_t EpotMGSolver::get_mgcyc( void ) const
{
    return( _mgcyc );
}

void EpotMGSolver::set_levels( uint32_t levels )
{
    if( levels < 1 )
	throw( Error( ERROR_LOCATION, "Invalid number of multigrid levels" ) ); 
    _levels = levels;
    reset_problem();
}


void EpotMGSolver::set_mgcycmax( uint32_t mgcycmax )
{
    _mgcycmax = mgcycmax;
}


void EpotMGSolver::set_mgeps( double mgeps )
{
    _mgeps = mgeps;
}


void EpotMGSolver::set_gamma( uint32_t gamma )
{
    if( gamma < 1 )
	throw( Error( ERROR_LOCATION, "Invalid multigrid cycle coefficient" ) ); 
    _gamma = gamma;
}


void EpotMGSolver::set_npre( uint32_t npre )
{
    if( npre < 1 )
	throw( Error( ERROR_LOCATION, "Invalid number of multigrid pre-smoothing cycles" ) ); 
    _npre = npre;
}


void EpotMGSolver::set_npost( uint32_t npost )
{
    if( npost < 1 )
	throw( Error( ERROR_LOCATION, "Invalid number of multigrid post-smoothing cycles" ) ); 
    _npost = npost;
}


/* *****************************************************************************
 * Common
 */


uint32_t EpotMGSolver::number_of_dimensions( void ) const
{
    switch( _geom.geom_mode() ) {
    case MODE_1D:
	return( 1 );
	break;
    case MODE_2D:
	return( 2 );
	break;
    case MODE_CYL:
	return( 2 );
	break;
    case MODE_3D:
	return( 3 );
	break;
    }

    throw( ErrorAssert( ERROR_LOCATION ) );
}


/* Create hierarchy of geometries and solvers.
 */
void EpotMGSolver::prepare_mg_geom( void )
{
    if( ibsimu.get_verbose_output() )
	std::cout << "  Preparing multigrid geometries\n";

    geom_mode_e geom_mode = _geom.geom_mode();
    Int3D size = _geom.size();
    Vec3D origo = _geom.origo();
    double h = _geom.h();

    for( uint32_t a = 0; a < _levels; a++ ) {

	if( a == 0 ) {

	    // First (finest) level
	    _geomv.push_back( &_geom );
	    EpotMGSubSolver *mgss = new EpotMGSubSolver( *this, _geom );
	    _epotsolverv.push_back( mgss );

	} else {

	    // Build mesh density	    
	    for( uint32_t b = 0; b < number_of_dimensions(); b++ ) {
		if( size[b] % 2 == 0 )
		    throw( Error( ERROR_LOCATION, "Incorrect parity of mesh size " + 
				  to_string(size[b]) + "in direction " + to_string(b) + 
				  " at level " + to_string(a) ) );		    
		size[b] = (size[b]+1)/2;
	    }
	    h *= 2.0;

	    // Build geometry
	    Geometry *geom = new Geometry( geom_mode, size, origo, h );
	    for( uint32_t b = 7; b <= _geom.number_of_boundaries(); b++ )
		geom->set_solid( b, _geom.get_solid(b) );
	    for( uint32_t b = 1; b <= _geom.number_of_boundaries(); b++ ) {
		Bound bound = _geom.get_boundary(b);
		bound.val = 0.0;
		geom->set_boundary( b, bound );
	    }
	    geom->build_mesh();
	    _geomv.push_back( geom );

	    // Make sub-solver copying parameters from this
	    EpotMGSubSolver *mgss = new EpotMGSubSolver( *this, *_geomv.back() );
	    _epotsolverv.push_back( mgss );
	}
    }

    _geom_prepared = true;
}


void EpotMGSolver::preprocess( MeshScalarField &epot, const MeshScalarField &scharge )
{
    if( ibsimu.get_verbose_output() )
	std::cout << "  Preprocessing\n";

    // Run preprocess geometry for all EpotSolvers, allocate rhs and epot fields 
    // and one working field
    for( uint32_t a = 0; a < _levels; a++ ) {

	if( a == 0 ) {
	    _epotv.push_back( &epot );
	} else {
	    MeshScalarField *ep = new MeshScalarField( (const Mesh)(*_geomv[a]) );
	    _epotv.push_back( ep );
	}

	MeshScalarField *rhs = new MeshScalarField( (const Mesh)(*_geomv[a]) );
	_rhsv.push_back( rhs );

	MeshScalarField *work = new MeshScalarField( (const Mesh)(*_geomv[a]) );
	_workv.push_back( work );

	if( !linear() ) {
	    if( a == 0 ) {
		_work2v.push_back( NULL );
	    } else {
		MeshScalarField *work2 = new MeshScalarField( (const Mesh)(*_geomv[a]) );
		_work2v.push_back( work2 );
	    }
	}

	// Preprocess solid meshes
	_epotsolverv[a]->preprocess( *_epotv[a] );
    }

    // Build rhs for top level
    for( uint32_t b = 0; b < _geom.nodecount(); b++ ) {

	uint32_t mesh = _geomv[0]->mesh(b);
	uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	if( node_id == SMESH_NODE_ID_NEAR_SOLID ||
	    node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		
	    // Ordinary vacuum/near solid
	    (*_rhsv[0])(b) = -epot.h()*epot.h()*scharge(b)/EPSILON0;

	} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
	    
	    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
	    if( _geom.geom_mode() == MODE_CYL && boundary == 3 ) {
		
		// Symmetry axis (vacuum)
		(*_rhsv[0])(b) = -epot.h()*epot.h()*scharge(b)/EPSILON0;
		
	    } else {
		
		// Ordinary Neumann node
		if( _neumann_order == 2 )
		    (*_rhsv[0])(b) = 2.0*epot.h()*_geom.get_boundary( boundary ).val;
		else
		    (*_rhsv[0])(b) = epot.h()*_geom.get_boundary( boundary ).val;
	    }
	}
    }
}


void EpotMGSolver::postprocess( void )
{
    if( ibsimu.get_verbose_output() )
	std::cout << "  Postprocessing\n";

    for( uint32_t a = 0; a < _levels; a++ ) {

	if( a != 0 )
	    delete _epotv[a];
	delete _rhsv[a];
	delete _workv[a];
	if( !linear() && a != 0 )
	    delete _work2v[a];

	// Postprocess solid meshes
	_epotsolverv[a]->postprocess();
    }

    _epotv.clear();
    _rhsv.clear();
    _workv.clear();
    _work2v.clear();
}


void EpotMGSolver::restrict_2d( MeshScalarField *out, const MeshScalarField *in, bool defect )
{
    // Go through internal nodes of rougher level (out)
    int32_t s = out->size(0)-1;
    int32_t t = out->size(1)-1;
    for( int32_t j = 1; j < t; j++ ) {
        for( int32_t i = 1; i < s; i++ ) {

	    int32_t ii = 2*i;
	    int32_t jj = 2*j;
	    double D = 
		1.0/16.0*(*in)( ii-1, jj-1 ) + 
		1.0/8.0 *(*in)( ii-1, jj   ) + 
		1.0/16.0*(*in)( ii-1, jj+1 ) + 

		1.0/8.0 *(*in)( ii,   jj-1 ) + 
		1.0/4.0 *(*in)( ii,   jj   ) + 
		1.0/8.0 *(*in)( ii,   jj+1 ) + 

		1.0/16.0*(*in)( ii+1, jj-1 ) + 
		1.0/8.0 *(*in)( ii+1, jj   ) + 
		1.0/16.0*(*in)( ii+1, jj+1 );

	    if( defect )
		(*out)(i,j) = 4.0*D;
	    else
		(*out)(i,j) = D;
	}
    }

    double cof = 1.0;
    if( _neumann_order == 1 && defect )
	cof = 2.0;

    // i=0 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii, jj-1 ) + 
	    1.0/2.0*(*in)( ii, jj   ) + 
	    1.0/4.0*(*in)( ii, jj+1 );
	
	(*out)(i,j) = cof*D;
    }

    // i=size(0)-1 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = s;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii, jj-1 ) + 
	    1.0/2.0*(*in)( ii, jj   ) + 
	    1.0/4.0*(*in)( ii, jj+1 );
	
	(*out)(i,j) = cof*D;
    }

    // j=0 boundary nodes of level+1
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj ) + 
	    1.0/2.0*(*in)( ii,   jj ) + 
	    1.0/4.0*(*in)( ii+1, jj );
	
	(*out)(i,j) = cof*D;
    }

    // j=size(1)-1 boundary nodes of level+1
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = t;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj ) + 
	    1.0/2.0*(*in)( ii,   jj ) + 
	    1.0/4.0*(*in)( ii+1, jj );
	
	(*out)(i,j) = cof*D;
    }

    // And corners
    int32_t i = 0;
    int32_t j = 0;
    double D = (*in)( 2*i, 2*j );
    (*out)(i,j) = cof*D;

    i = s;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = cof*D;

    j = t;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = cof*D;

    i = 0;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = cof*D;
}


void EpotMGSolver::restrict_1d( MeshScalarField *out, const MeshScalarField *in, bool defect )
{
    // Go through internal nodes of rougher level (out)
    int32_t s = out->size(0)-1;
    for( int32_t i = 1; i < s; i++ ) {

	int32_t ii = 2*i;
	double D = 
	    1.0/4.0*(*in)( ii-1 ) + 
	    1.0/2.0*(*in)( ii   ) + 
	    1.0/4.0*(*in)( ii+1 );
	
	if( defect )
	    (*out)(i) = 4.0*D;
	else
	    (*out)(i) = D;
    }

    // Boundary nodes completely separated
    if( _neumann_order == 1 && defect ) {
	(*out)(0) = 2.0*(*in)(0);
	(*out)(s) = 2.0*(*in)(2*s);
    } else {
	(*out)(0) = (*in)(0);
	(*out)(s) = (*in)(2*s);
    }
}


void EpotMGSolver::restrict( MeshScalarField *out, const MeshScalarField *in, bool defect )
{
    switch( _geom.geom_mode() ) {
    case MODE_1D:
	restrict_1d( out, in, defect );
	break;
    case MODE_2D:
	restrict_2d( out, in, defect );
	break;
    default:
	break;
    }
}


void EpotMGSolver::prolong_add_2d( MeshScalarField *out, int32_t i, int32_t j, double C )
{
    if( i >= 0 && i < (int32_t)out->size(0) &&
	j >= 0 && j < (int32_t)out->size(1) )
	(*out)(i,j) += C;
}


void EpotMGSolver::prolong_2d( MeshScalarField *out, const MeshScalarField *in )
{
    // Clear output field
    out->clear();

    // Loop through all input nodes
    int32_t s = in->size(0);
    int32_t t = in->size(1);
    for( int32_t j = 0; j < t; j++ ) {
	for( int32_t i = 0; i < s; i++ ) {

	    int32_t ii = 2*i;
	    int32_t jj = 2*j;
	    double C = (*in)(i,j);
	    prolong_add_2d( out, ii-1, jj-1, 1.0/4.0*C );
	    prolong_add_2d( out, ii-1, jj,   1.0/2.0*C );
	    prolong_add_2d( out, ii-1, jj+1, 1.0/4.0*C );

	    prolong_add_2d( out, ii,   jj-1, 1.0/2.0*C );
	    prolong_add_2d( out, ii,   jj,   1.0*C );
	    prolong_add_2d( out, ii,   jj+1, 1.0/2.0*C );

	    prolong_add_2d( out, ii+1, jj-1, 1.0/4.0*C );
	    prolong_add_2d( out, ii+1, jj,   1.0/2.0*C );
	    prolong_add_2d( out, ii+1, jj+1, 1.0/4.0*C );
	}
    }    
}


void EpotMGSolver::prolong_add_1d( MeshScalarField *out, int32_t i, double C )
{
    if( i >= 0 && i < (int32_t)out->size(0) )
	(*out)(i) += C;
}


void EpotMGSolver::prolong_1d( MeshScalarField *out, const MeshScalarField *in )
{
    // Clear output field
    out->clear();

    // Loop through all input nodes
    int32_t s = in->size(0);
    for( int32_t i = 0; i < s; i++ ) {
	
	int32_t ii = 2*i;
	double C = (*in)(i);
	prolong_add_1d( out, ii-1, 1.0/2.0*C );
	prolong_add_1d( out, ii,       1.0*C );
	prolong_add_1d( out, ii+1, 1.0/2.0*C );
    }    
}


void EpotMGSolver::prolong( MeshScalarField *out, const MeshScalarField *in )
{
    switch( _geom.geom_mode() ) {
    case MODE_1D:
	prolong_1d( out, in );
	break;
    case MODE_2D:
	prolong_2d( out, in );
	break;
    default:
	break;
    }
}

// Make correction Xnew=X+V
void EpotMGSolver::correct( const Geometry *geom, MeshScalarField *sol, const MeshScalarField *corr )
{
    // Loop through all nodes, only correct non-fixed nodes
    int32_t s = geom->size(0);
    int32_t t = geom->size(1);
    for( int32_t j = 0; j < t; j++ ) {
	for( int32_t i = 0; i < s; i++ ) {

	    if( !(geom->mesh(i,j) & SMESH_NODE_FIXED) )
		(*sol)(i,j) += (*corr)(i,j);
	}
    }
}


void EpotMGSolver::mg_recurse( uint32_t level )
{
    //std::cout << "MG recursion, level " << level << "\n";

    if( level == _levels-1 ) {
 
	// Linear case: clear last field, expected result closer to zero than result 
	// from last round. Never clear top level field.
	if( linear() && level != 0 )
	    _epotv[level]->clear();

	if( ibsimu.get_verbose_output() ) {
	    for( uint32_t b = 0; b < level; b++ )
		std::cout << "  ";
	    std::cout << "  Roughest level solution\n";
	}
	
	// Last level, solve the roughest problem until convergence
	uint32_t a = 0;
	for( a = 1; a <= _imax; a++ ) {
	    _res = _res_coef * _epotsolverv[level]->mg_solve( _epotv[level], _rhsv[level], _w );
	    if( _res < _eps )
		break;
	}

	//std::cout << "epot (level = " << level << ") accurate solve:\n";
	//print_field( _epotv[level] );

	if( ibsimu.get_verbose_output() ) {
	    for( uint32_t b = 0; b < level; b++ )
		std::cout << "  ";
	    std::cout << "  " << a << " iterations done\n";
	    for( uint32_t b = 0; b < level; b++ )
		std::cout << "  ";
	    std::cout << "  " << _res << " accuracy reached\n";
	    if( a >= _imax ) {
		for( uint32_t b = 0; b < level; b++ )
		    std::cout << "  ";
		std::cout << "  maximum number of iterations done\n";
	    }
	}

	return;
    } 

    // Do gamma cycles of next level
    for( uint32_t a = 0; a < _gamma; a++ ) {

	if( ibsimu.get_verbose_output() ) {
	    for( uint32_t b = 0; b < level; b++ )
		std::cout << "  ";
	    if( _gamma > 1 )
		std::cout << "  Doing gamma cycle " << a+1 << "/" << _gamma << " at level " << level << "\n";
	    else
		std::cout << "  Doing cycle at level " << level << "\n";
	}

	//std::cout << "epot (level = " << level << "):\n";
	//print_field( _epotv[level] );

	// Pre smoothing
	for( uint32_t a = 0; a < _npre; a++ )
	    _epotsolverv[level]->mg_smooth( _epotv[level], _rhsv[level] );

	//std::cout << "epot (level = " << level << ") pre smoothed:\n";
	//print_field( _epotv[level] );

	// Defect
	//std::cout << "  Calculating defect for level " << level << "\n";
	_epotsolverv[level]->defect( _workv[level], _epotv[level], _rhsv[level] );

	//std::cout << "defect (level = " << level << "):\n";
	//print_field( _workv[level] );

	if( linear() ) {

	    //std::cout << "  Restricing defect from level " << level << " to level " << level+1 << "\n";
	    restrict( _rhsv[level+1], _workv[level], true );

	    (*_rhsv[level+1]) *= -1.0;

	    //std::cout << "defect (level = " << level+1 << "):\n";
	    //print_field( _rhsv[level+1] );

	} else {

	    // Restrict defect
	    //std::cout << "  Restricing defect from level " << level << " to level " << level+1 << "\n";
	    restrict( _workv[level+1], _workv[level], true );

	    //std::cout << "defect (level = " << level+1 << "):\n";
	    //print_field( _workv[level+1] );

	    // Restrict solution approximation
	    //std::cout << "  Restricing solution from level " << level << " to level " << level+1 << "\n";
	    restrict( _epotv[level+1], _epotv[level], false );

	    // Store solution
	    (*_work2v[level+1]) = (*_epotv[level+1]);

	    //std::cout << "solution (level = " << level+1 << "):\n";
	    //print_field( _epotv[level+1] );

	    // Calculate defect
	    //std::cout << "  Calculating right hand side for level " << level+1 << "\n";
	    _epotsolverv[level+1]->defect( _rhsv[level+1], _epotv[level+1], _workv[level+1] );

	    //std::cout << "rhs (level = " << level+1 << "):\n";
	    //print_field( _rhsv[level+1] );
	}

	// Recurse to next level to solve for correction V (epot)
	// A*V=-D
	mg_recurse( level+1 );
	
	if( linear() ) {

	    //std::cout << "correction (level = " << level+1 << "):\n";
	    //print_field( _epotv[level+1] );

	    // Prolong correction (epot) to work field
	    //std::cout << "  Prolonging correction from level " << level+1 << " to level " << level << "\n";
	    prolong( _workv[level], _epotv[level+1] );
	    
	} else {

	    //std::cout << "refined solution (level = " << level+1 << "):\n";
	    //print_field( _epotv[level+1] );

	    // Compute correction
	    (*_workv[level+1]) = (*_epotv[level+1]);
	    (*_workv[level+1]) -= (*_work2v[level+1]);

	    //std::cout << "correction (level = " << level+1 << "):\n";
	    //print_field( _workv[level+1] );

	    //std::cout << "  Prolonging correction from level " << level+1 << " to level " << level << "\n";
	    prolong( _workv[level], _workv[level+1] );
	}

	//std::cout << "correction (level = " << level << "):\n";
	//print_field( _workv[level] );

	// Make correction Xnew=X+V
	//std::cout << "  Calculating correction for level " << level << "\n";
	correct( _geomv[level], _epotv[level], _workv[level] );

	//std::cout << "epot (level = " << level << "):\n";
	//print_field( _epotv[level] );

	// Post smoothing
	for( uint32_t a = 0; a < _npost; a++ )
	    _res = _res_coef * _epotsolverv[level]->mg_smooth( _epotv[level], _rhsv[level] );

	//std::cout << "epot (level = " << level << ") post smoothed:\n";
	//print_field( _epotv[level] );

	if( ibsimu.get_verbose_output() ) {
	    for( uint32_t b = 0; b < level; b++ )
		std::cout << "  ";
	    std::cout << "  " << _res << " accuracy reached\n";
	}
     }
}


void EpotMGSolver::print_field( const MeshScalarField *F )
{
    for( size_t i = 0; i < F->size(0); i++ ) {
	std::cout << std::setw(8) << (*F)(i) << " ";
    }
    std::cout << "\n";
}


void EpotMGSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    //std::cout << "\n";

    if( ibsimu.get_verbose_output() ) {
	std::cout << "  Using Multigrid solver (" 
		  << "levels = " << _levels
		  << ", npre = " << _npre
		  << ", npost = " << _npost
		  << ", mgcycmax = " << _mgcycmax
		  << ", mgeps = " << _mgeps
		  << ", gamma = " << _gamma
		  << ", w = " << _w
		  << ", eps = " << _eps
		  << ", imax = " << _imax
		  << ")\n";
	if( linear() )
	    std::cout << "  Linear problem\n";
	else
	    std::cout << "  Nonlinear problem\n";
    }

    if( !_geom_prepared )
	prepare_mg_geom();
    preprocess( epot, scharge );

    //std::cout << "rhs (level = " << 0 << "):\n";
    //print_field( _rhsv[0] );

    // Call iterator mgcyc times
    uint32_t a;
    for( a = 0; a < _mgcycmax; a++ ) {
	if( ibsimu.get_verbose_output() )
	    std::cout << "  Cycle " << a+1 << "\n";
	mg_recurse( 0 );
	if( _res < _mgeps )
	    break;
    }
    _mgcyc = a;

    if( ibsimu.get_verbose_output() ) {
	if( a == _mgcycmax )
	    std::cout << "  Maximum number of cycles done\n";
	else
	    std::cout << "  Convergence criterion reached\n";
    }

    postprocess();
}


void EpotMGSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}




void EpotMGSolver::debug_print( std::ostream &os ) const 
{
    EpotSolver::debug_print( os );
    os << "**EpotMGSolver\n";
    os << "levels = " << _levels << "\n";
    os << "npre = " << _npre << "\n";
    os << "npost = " << _npost << "\n";
    os << "mgcyc = " << _mgcyc << "\n";
    os << "mgeps = " << _mgeps << "\n";
    os << "gamma = " << _gamma << "\n";
    os << "res = " << _res << "\n";
    os << "eps = " << _eps << "\n";
}

