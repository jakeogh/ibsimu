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
#include "statusprint.hpp"


EpotMGSolver::EpotMGSolver( Geometry &geom )
    : EpotSolver( geom ), _geom_prepared(false), _levels(1), _npre(5), 
      _npost(5), _ncyc(1), _res(0.0), _eps(1.0e-6)
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
    for( uint32_t a = 0; a < _epotv.size(); a++ ) {
	if( a >= 1 && _epotv[a] ) {
	    delete _epotv[a];
	    delete _rhsv[a];
	}
    }

    for( uint32_t a = 0; a < _geomv.size(); a++ ) {
	if( a >= 1 && _geomv[a] ) {
	    delete _geomv[a];
	    //delete _epotsolverv[a];
	}
    }

    _epotv.clear();
    _geomv.clear();
    //_epotsolverv.clear();
    _rhsv.clear();
    _res = 0.0;
    _geom_prepared = false;
}


void EpotMGSolver::set_eps( double eps )
{
    _eps = eps;
}

double EpotMGSolver::get_residual( void ) const
{
    return( _res );
}


void EpotMGSolver::set_levels( uint32_t levels )
{
    _levels = levels;
    reset_problem();
}


void EpotMGSolver::set_ncyc( uint32_t ncyc )
{
    ncyc = ncyc;
}


void EpotMGSolver::set_npre( uint32_t npre )
{
    npre = npre;
}


void EpotMGSolver::set_npost( uint32_t npost )
{
    npost = npost;
}


/* *****************************************************************************
 * 2D
 */

double EpotMGSolver::gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
					       uint32_t a, uint32_t dj ) const
{
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    double cof = 2.0/(alpha*beta);
    double epf = 2.0/(alpha+beta)*( (*_epot)(a-1)/alpha + (*_epot)(a+1)/beta );

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    cof += 2.0/(alpha*beta);
    epf += 2.0/(alpha+beta)*( (*_epot)(a-dj)/alpha + (*_epot)(a+dj)/beta );

    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( p + ( epf - cof*p - (*_rhs)(a) - rhst ) / ( cof + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( epf - cof*p - (*_rhs)(a) - rhst ) / ( cof + drhst ) );	
    }

    return( (1.0/cof) * ( epf - (*_rhs)(a) ) );
}


double EpotMGSolver::gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(a+1) + (*_epot)(a-1) 
		      + (*_epot)(a+dj) + (*_epot)(a-dj) - 4.0*p
		      - (*_rhs)(a) - rhst ) / ( 4.0 + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(a+1) + (*_epot)(a-1) 
		      + (*_epot)(a+dj) + (*_epot)(a-dj) - 4.0*p
		      - (*_rhs)(a) - rhst ) / ( 4.0 + drhst ) );
    }

    return( (1.0/4.0) * ( (*_epot)(a+1) + (*_epot)(a-1) +
			  (*_epot)(a+dj) + (*_epot)(a-dj) - (*_rhs)(a) ) );
}


double EpotMGSolver::gs_process_neumann_2d( uint32_t boundary, uint32_t a, uint32_t dj ) const
{
    switch( boundary ) {
    case 1:
	if( _neumann_order == 2 )
	    // (3*phi_i - 4*phi_i+1 + phi_i+2) / 2h = q_0
	    return( (4.0*(*_epot)(a+1) - (*_epot)(a+2) + (*_rhs)(a))/3.0 );
	else
	    // (phi_i - phi_i+1) / h = q_0
	    return( (*_epot)(a+1) + (*_rhs)(a) );
	break;
    case 2:
	if( _neumann_order == 2 )
	    // (3*phi_i - 4*phi_i-1 + phi_i-2) / 2h = q_0
	    return( (4.0*(*_epot)(a-1) - (*_epot)(a-2) + (*_rhs)(a))/3.0 );
	else
	    // (phi_i-1 - phi_i) / h = q_0
	    return( (*_epot)(a-1) + (*_rhs)(a) );
	break;
    case 3:
	if( _neumann_order == 2 )
	    // (3*phi_j - 4*phi_j+1 + phi_j+2) / 2h = q_0
	    return( (4.0*(*_epot)(a+dj) - (*_epot)(a+2*dj) + (*_rhs)(a))/3.0 );
	else
	    // (phi_j - phi_j+1) / h = q_0
	    return( (*_epot)(a+dj) + (*_rhs)(a) );
	break;
    case 4:
	if( _neumann_order == 2 )
	    // (3*phi_j - 4*phi_j-1 + phi_j-2) / 2h = q_0
	    return( (4.0*(*_epot)(a-dj) - (*_epot)(a-2*dj) + (*_rhs)(a))/3.0 );
	else
	    // (phi_j-1 - phi_j) / h = q_0
	    return( (*_epot)(a-dj) + (*_rhs)(a) );
	break;
    }

    return( 0.0 );
}


double EpotMGSolver::gs_loop_2d( void ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    const uint32_t dj = _geom.size(0);
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	uint32_t a = j*dj;
	for( uint32_t i = 0; i < _geom.size(0); i++, a++ ) {
	    
	    double Vold = (*_epot)(a);
	    double Vnew;
	    uint32_t mesh = _geom.mesh(a);
	    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		Vnew = gs_process_near_solid_2d( nearsolid_ptr, a, dj );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		Vnew = gs_process_pure_vacuum_2d( a, dj );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		Vnew = gs_process_neumann_2d( boundary, a, dj );
	    } else {
		// Dirichlet
		continue;
	    }
	    (*_epot)(a) = Vnew;
	    double err = fabs( Vnew - Vold );
	    if( err > maxerr )
		maxerr = err;
	    if( comp_isinf(err) ) {
		throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) + 
			      ", " + to_string(j) ) );
	    } else if( comp_isnan(err) ) {
		throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) + 
			      ", " + to_string(j) ) );
	    }
	}
    }
    return( maxerr );
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
    geom_mode_e geom_mode = _geom.geom_mode();
    Int3D size = _geom.size();
    Vec3D origo = _geom.origo();
    double h = _geom.h();

    // _geom is given, construct _geomv vector
    for( uint32_t a = 0; a < _levels; a++ ) {
	if( a == 0 ) {
	    _geomv.push_back( &_geom );
	    //_epotsolverv.push_back( (EpotSolver *)this );
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
	    for( uint32_t b = 7; b < _geom.number_of_boundaries(); b++ )
		geom->set_solid( b, _geom.get_solid(b) );
	    for( uint32_t b = 0; b < _geom.number_of_boundaries(); b++ ) {
		Bound bound = _geom.get_boundary(b);
		bound.val = 0.0;
		geom->set_boundary( b, bound );
	    }
	    geom->build_mesh();
	    _geomv.push_back( geom );
	}
    }
}


void EpotMGSolver::preprocess( MeshScalarField &epot, const MeshScalarField &scharge )
{
    // Run preprocess geometry for all EpotSolvers, allocate rhs and epot fields
    for( uint32_t a = 0; a < _levels; a++ ) {

	if( a == 0 ) {
	    _epotv.push_back( &epot );
	} else {
	    MeshScalarField *ep = new MeshScalarField( (const Mesh)(*_geomv[a]) );
	    _epotv.push_back( ep );
	}

	MeshScalarField *rhs = new MeshScalarField( (const Mesh)(*_geomv[a]) );
	_rhsv.push_back( rhs );

	// Preprocess solid meshes
	//_epotsolverv[a]->preprocess( _epotv[a] );
    }

    // Build rhs
    for( uint32_t a = 0; a < _geom.nodecount(); a++ ) {

	uint32_t mesh = _geom.mesh(a);
	uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	if( node_id == SMESH_NODE_ID_NEAR_SOLID ||
	    node_id == SMESH_NODE_ID_PURE_VACUUM ) {

	    // Ordinary vacuum/near solid
	    (*_rhs)(a) = -scharge(a)*_geom.h()*_geom.h()/EPSILON0;

	} else if( node_id == SMESH_NODE_ID_NEUMANN ) {

	    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
	    if( _geom.geom_mode() == MODE_CYL && boundary == 3 ) {
			
		// Symmetry axis (vacuum)
		(*_rhs)(a) = -scharge(a)*_geom.h()*_geom.h()/EPSILON0;

	    } else {
		
		// Ordinary Neumann node
		if( _neumann_order == 2 )
		    (*_rhs)(a) = 2.0*_geom.h()*_geom.get_boundary( boundary ).val;
		else
		    (*_rhs)(a) = _geom.h()*_geom.get_boundary( boundary ).val;
	    }
	}
    }
}


void EpotMGSolver::postprocess( void )
{
    delete _rhs;
    _rhs = NULL;

    EpotSolver::postprocess();
}


void EpotMGSolver::mg_smooth( void )
{
    if( _geom.geom_mode() == MODE_3D )
	_res = gs_loop_3d();
    else if( _geom.geom_mode() == MODE_2D )
	_res = gs_loop_2d();
    else if( _geom.geom_mode() == MODE_CYL )
	_res = gs_loop_cyl();
    else if( _geom.geom_mode() == MODE_1D )
	_res = gs_loop_1d();
    else
	throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void EpotMGSolver::mg_recurse( uint32_t level )
{
    // Do ncyc cycles of next level
    for( uint32_t a = 0; a < _ncyc; a++ ) {

	// Pre smoothing
	mg_smooth();

	if( level == _levels ) {
	    // Last level, solve problem
	} else {
	    // Recurse to next level
	    mg_recurse( level+1 );
	}

	// Post smoothing
	mg_smooth();
    }
}


void EpotMGSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    if( ibsimu.get_verbose_output() ) {
	std::cout << "  Using Multigrid solver (" 
		  << "levels = " << _levels
		  << ", npre = " << _npre
		  << ", npost = " << _npost
		  << ", ncyc = " << _ncyc
		  << ", eps = " << _eps
		  << ")\n";
    }

    if( !_geom_prepared )
	prepare_mg_geom();
    preprocess( epot, scharge );

    // Do first level iteration only once
    mg_recurse( 1 );

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
    os << "ncyc = " << _ncyc << "\n";
    os << "res = " << _res << "\n";
    os << "eps = " << _eps << "\n";
}

