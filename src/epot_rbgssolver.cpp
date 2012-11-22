/*! \file epot_rbgssolver.cpp
 *  \brief Red-Black Gauss-Seidel solver for electric potential problem
 */

/* Copyright (c) 2011,2012 Taneli Kalvas. All rights reserved.
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


#include "epot_rbgssolver.hpp"
#include "ibsimu.hpp"
#include "timer.hpp"
#include "constants.hpp"
#include "compmath.hpp"


EpotRBGSSolver::EpotRBGSSolver( Geometry &geom )
    : EpotSolver( geom ), _rhs(NULL), _imax(10000), _eps(1.0e-6), _w(1.0), _sp(NULL)
{
    pthread_mutex_init( &_mutex, NULL );
    pthread_cond_init( &_cond, NULL );
}



EpotRBGSSolver::EpotRBGSSolver( Geometry &geom, std::istream &s )
    : EpotSolver(geom)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );    
}


EpotRBGSSolver::~EpotRBGSSolver()
{
    if( _rhs )
	delete _rhs;
    if( _sp )
	delete _sp;
    pthread_mutex_destroy( &_mutex );
    pthread_cond_destroy( &_cond );
}


void EpotRBGSSolver::set_eps( double eps )
{
    _eps = eps;
}

double EpotRBGSSolver::get_residual( void ) const
{
    return( _res );
}


void EpotRBGSSolver::set_imax( uint32_t imax )
{
    _imax = imax;
}


void EpotRBGSSolver::set_w( double w )
{
    _w = w;
}


/* *****************************************************************************
 * 3D
 */

double EpotRBGSSolver::gs_process_near_solid_3d( const uint8_t *nearsolid_ptr, uint32_t a,
						 uint32_t dj, uint32_t dk ) const
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

    // Zmin direction
    alpha = 1.0;
    if( sflag & 0x10 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Zmax direction
    beta = 1.0;
    if( sflag & 0x20 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Z axis
    cof += 2.0/(alpha*beta);
    epf += 2.0/(alpha+beta)*( (*_epot)(a-dk)/alpha + (*_epot)(a+dk)/beta );


    return( (1.0/cof) * ( epf - (*_rhs)(a) ) );
}


double EpotRBGSSolver::gs_process_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const
{
    return( (1.0/6.0) * ( (*_epot)(a+1)  + (*_epot)(a-1) +
			  (*_epot)(a+dj) + (*_epot)(a-dj) + 
			  (*_epot)(a+dk) + (*_epot)(a-dk) - (*_rhs)(a) ) );
}


double EpotRBGSSolver::gs_process_neumann_3d( uint32_t boundary, uint32_t a,
					      uint32_t dj, uint32_t dk ) const
{
    switch( boundary ) {
    case 1:
	// (phi_i - phi_i+1) / h = q_0
	return( (*_epot)(a+1) + (*_rhs)(a) );
	break;
    case 2:
	// (phi_i-1 - phi_i) / h = q_0
	return( (*_epot)(a-1) + (*_rhs)(a) );
	break;
    case 3:
	// (phi_j - phi_j+1) / h = q_0
	return( (*_epot)(a+dj) + (*_rhs)(a) );
	break;
    case 4:
	// (phi_j-1 - phi_j) / h = q_0
	return( (*_epot)(a-dj) + (*_rhs)(a) );
	break;
    case 5:
	// (phi_k - phi_k+1) / h = q_0
	return( (*_epot)(a+dk) + (*_rhs)(a) );
	break;
    case 6:
	// (phi_k-1 - phi_k) / h = q_0
	return( (*_epot)(a-dk) + (*_rhs)(a) );
	break;
    }

    return( 0.0 );
}


double EpotRBGSSolver::iterator_loop( uint32_t thno, uint32_t thcount, uint32_t rb, uint32_t dj, uint32_t dk )
{
    const double w2 = 1.0-_w;
    double maxerr = 0.0;

    uint32_t kmax = _geom.size(2)*(thno+1)/thcount;
    for( uint32_t k = _geom.size(2)*thno/thcount; k < kmax; k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    
	    uint32_t i = (k+j) % 2;
	    if( i == rb )
		i = 0;
	    else
		i = 1;
		    
	    uint32_t a = k*dk+j*dj+i;
	    for( ; i < _geom.size(0); i += 2, a += 2 ) {

		double Vold = (*_epot)(a);
		double Vnew;
		uint32_t mesh = _geom.mesh(a);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		    Vnew = gs_process_near_solid_3d( nearsolid_ptr, a, dj, dk );
		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		    Vnew = gs_process_pure_vacuum_3d( a, dj, dk );
		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    Vnew = gs_process_neumann_3d( boundary, a, dj, dk );
		} else {
		    // Dirichlet
		    continue;
		}
		Vnew = _w*Vnew + w2*Vold;
		(*_epot)(a) = Vnew;
		double err = fabs( Vnew - Vold );
		if( err > maxerr )
		    maxerr = err;
		if( comp_isinf(err) ) {
		    throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) + 
				  ", " + to_string(j) + ", " + to_string(k) ) );
		} else if( comp_isnan(err) ) {
		    throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) + 
				  ", " + to_string(j) + ", " + to_string(k) ) );
		}
	    }
	}
    }

    return( maxerr );
}


void *EpotRBGSSolver::iterator_main( uint32_t thno, uint32_t thcount )
{
    const uint32_t dj = _geom.size(0);
    const uint32_t dk = _geom.size(0)*_geom.size(1);

    uint32_t rb = 0;
    uint32_t iter_end = _imax*2;
    uint32_t iter = 0;
    double maxerr = 0.0;
    while( 1 ) {

	try {
	    maxerr = iterator_loop( thno, thcount, rb, dj, dk );
	} catch( Error e ) {
	    pthread_mutex_lock( &_mutex );
	    _err = e;
	    _error = true;
	    _done = true;
	    pthread_mutex_unlock( &_mutex );
	}

	iter++;

	// Synchronize threads
	pthread_mutex_lock( &_mutex );

	// Store largest residual
	if( maxerr > _res )
	    _res = maxerr;

	if( rb == 0 )
	    rb = 1;
	else
	    rb = 0;

	_done_count++;

	if( _done_count == thcount ) {
	    // Last thread

	    if( iter == iter_end )
		_done = true;
	    _iter = iter/2;
	    _done_count = 0;

	    // If done both Red and Black
	    if( rb == 0 ) {
		std::stringstream ss;
		ss << "  " << std::setw(5) << _iter << " " << std::scientific << std::setw(20) << _res;
		_sp->print( ss.str() );
		if( _res < _eps )
		    _done = true;
		else
		    _res = 0.0;
	    }

	    pthread_cond_broadcast( &_cond );
	} else {
	    
	    // Wait until all threads done
	    pthread_cond_wait( &_cond, &_mutex );
	}

	// Quit iterator
	if( _done )
	    break;

	pthread_mutex_unlock( &_mutex );
    }
    
    pthread_mutex_unlock( &_mutex );
    return( NULL );
}


struct EpotIteratorData {
    EpotRBGSSolver        *_solver;
    pthread_t              _thread;
    uint32_t               _thno;
    uint32_t               _thcount;
};


void *EpotRBGSSolver::iterator_entry( void *data )
{
    EpotIteratorData *itdata = (EpotIteratorData *)data;
    EpotRBGSSolver *solver = itdata->_solver;
    return( solver->iterator_main( itdata->_thno, itdata->_thcount ) );
}


/* *****************************************************************************
 * 2D
 */

/*
double EpotRBGSSolver::gs_process_near_solid_2d( MeshScalarField &epot, const MeshScalarField &rhs,
					       const uint8_t *nearsolid_ptr, 
					       uint32_t i, uint32_t j ) const
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
    double epf = 2.0/(alpha+beta)*( epot(i-1,j)/alpha + epot(i+1,j)/beta );

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
    epf += 2.0/(alpha+beta)*( epot(i,j-1)/alpha + epot(i,j+1)/beta );

    return( (1.0/cof) * ( epf - rhs(i,j) ) );
}


double EpotRBGSSolver::gs_process_pure_vacuum_2d( MeshScalarField &epot, const MeshScalarField &rhs,
						uint32_t i, uint32_t j ) const
{
    return( (1.0/4.0) * ( epot(i+1,j) + epot(i-1,j) +
			  epot(i,j+1) + epot(i,j-1) - rhs(i,j) ) );
}


double EpotRBGSSolver::gs_process_neumann_2d( MeshScalarField &epot, const MeshScalarField &rhs,
					    uint32_t boundary, uint32_t i, uint32_t j ) const
{
    switch( boundary ) {
    case 1:
	// (phi_i - phi_i+1) / h = q_0
	return( epot(i+1,j) + rhs(i,j) );
	break;
    case 2:
	// (phi_i-1 - phi_i) / h = q_0
	return( epot(i-1,j) + rhs(i,j) );
	break;
    case 3:
	// (phi_j - phi_j+1) / h = q_0
	return( epot(i,j+1) + rhs(i,j) );
	break;
    case 4:
	// (phi_j-1 - phi_j) / h = q_0
	return( epot(i,j-1) + rhs(i,j) );
	break;
    }

    return( 0.0 );
}


double EpotRBGSSolver::gs_loop_2d( MeshScalarField &epot, const MeshScalarField &rhs ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    double w2 = 1.0-_w;
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	//double y = j*_geom.h()+_geom.origo(1);
	for( uint32_t i = 0; i < _geom.size(0); i++ ) {
	    //double x = i*_geom.h()+_geom.origo(0);
	    
	    double Vold = epot(i,j);
	    double Vnew;
	    uint32_t mesh = _geom.mesh(i,j);
	    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		Vnew = gs_process_near_solid_2d( epot, rhs, nearsolid_ptr, i, j );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		Vnew = gs_process_pure_vacuum_2d( epot, rhs, i, j );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		Vnew = gs_process_neumann_2d( epot, rhs, boundary, i, j );
	    } else {
		// Dirichlet
		continue;
	    }
	    Vnew = _w*Vnew + w2*Vold;
	    epot(i,j) = Vnew;
	    double err = fabs( Vnew - Vold );
	    if( err > maxerr )
		maxerr = err;
	}
    }
    return( maxerr );
}
*/

/* *****************************************************************************
 * 1D
 */

/*
double EpotRBGSSolver::gs_process_near_solid_1d( MeshScalarField &epot, const MeshScalarField &rhs,
					       const uint8_t *nearsolid_ptr, 
					       uint32_t i ) const
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
    double epf = 2.0/(alpha+beta)*( epot(i-1)/alpha + epot(i+1)/beta );

    return( (1.0/cof) * ( epf - rhs(i) ) );
}


double EpotRBGSSolver::gs_process_pure_vacuum_1d( MeshScalarField &epot, const MeshScalarField &rhs,
						uint32_t i ) const
{
    return( (1.0/2.0) * ( epot(i+1) + epot(i-1) - rhs(i) ) );
}


double EpotRBGSSolver::gs_process_neumann_1d( MeshScalarField &epot, const MeshScalarField &rhs,
					    uint32_t boundary, uint32_t i ) const
{
    switch( boundary ) {
    case 1:
	// (phi_i - phi_i+1) / h = q_0
	return( epot(i+1) + rhs(i) );
	break;
    case 2:
	// (phi_i-1 - phi_i) / h = q_0
	return( epot(i-1) + rhs(i) );
	break;
    }

    return( 0.0 );
}


double EpotRBGSSolver::gs_loop_1d( MeshScalarField &epot, const MeshScalarField &rhs ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    double w2 = 1.0-_w;
    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
	//double x = i*_geom.h()+_geom.origo(0);
	    
	double Vold = epot(i);
	double Vnew;
	uint32_t mesh = _geom.mesh(i);
	uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
	    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
	    Vnew = gs_process_near_solid_1d( epot, rhs, nearsolid_ptr, i );
	} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
	    Vnew = gs_process_pure_vacuum_1d( epot, rhs, i );
	} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
	    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
	    Vnew = gs_process_neumann_1d( epot, rhs, boundary, i );
	} else {
	    // Dirichlet
	    continue;
	}
	Vnew = _w*Vnew + w2*Vold;
	epot(i) = Vnew;
	double err = fabs( Vnew - Vold );
	if( err > maxerr )
	    maxerr = err;
    }

    return( maxerr );
}
*/

/* *****************************************************************************
 * Common
 */


void EpotRBGSSolver::preprocess( const MeshScalarField &scharge )
{
    // Build right-hand-side
    if( _rhs )
	delete _rhs;
    _rhs = new MeshScalarField( (const Mesh)scharge );

    // Change near solid nodes on Neumann boundaries to
    // NODE_ID_NEUMANN and store near solid indexes. Take care to
    // process in the same order as in Geometry class (x overrides y,
    // which overrides z.

    // Clear near solid indexes vector
    _nsind.clear();

    // Xmin and Xmax
    for( uint32_t bound = 1; bound <= 2; bound++ ) {
	uint32_t i = 0;
	if( bound == 2 ) i = _geom.size(0)-1;
	if( _geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
	    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
		for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		    uint32_t mesh = _geom.mesh(i,j,k);
		    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
			uint32_t index = mesh & SMESH_NEAR_SOLID_INDEX_MASK;
			_nsind.push_back( index );
			_geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
		    }
		}
	    }
	}
    }
    if( _geom.geom_mode() == MODE_2D || _geom.geom_mode() == MODE_CYL ||
	_geom.geom_mode() == MODE_3D ) {
	// Ymin and Ymax
	for( uint32_t bound = 3; bound <= 4; bound++ ) {
	    uint32_t j = 0;
	    if( bound == 4 ) j = _geom.size(1)-1;
	    if( _geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
		for( uint32_t k = 0; k < _geom.size(2); k++ ) {
		    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
			uint32_t mesh = _geom.mesh(i,j,k);
			uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
			if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
			    uint32_t index = mesh & SMESH_NEAR_SOLID_INDEX_MASK;
			    _nsind.push_back( index );
			    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
			}
		    }
		}
	    }
	}
    }
    if( _geom.geom_mode() == MODE_3D ) {
	// Zmin and Zmax
	for( uint32_t bound = 5; bound <= 6; bound++ ) {
	    uint32_t k = 0;
	    if( bound == 6 ) k = _geom.size(2)-1;
	    if( _geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
		for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
			uint32_t mesh = _geom.mesh(i,j,k);
			uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
			if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
			    uint32_t index = mesh & SMESH_NEAR_SOLID_INDEX_MASK;
			    _nsind.push_back( index );
			    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
			}
		    }
		}
	    }
	}
    }

    // Build rhs and set forced vacuum nodes and dirichlet nodes to
    // epot. Mark fixed vacuum nodes with a tag.
    Vec3D x;
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	x[2] = _geom.origo(2) + _geom.h()*k;
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    x[1] = _geom.origo(1) + _geom.h()*j;
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		x[0] = _geom.origo(0) + _geom.h()*i;

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID ||
		    node_id == SMESH_NODE_ID_PURE_VACUUM ) {

		    // Vacuum
		    if( _force_pot_func && (*_force_pot_func)( x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			(*_epot)(i,j,k) = _force_pot;

		    } else if ( _init_plasma_func && (*_init_plasma_func)( x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			(*_epot)(i,j,k) = _Up;

		    } else {

			// Ordinary vacuum/near solid
			(*_rhs)(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
		    }

		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {

		    // Neumann
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    (*_rhs)(i,j,k) = _geom.h()*_geom.get_boundary( boundary ).value( x );

		} else if( node_id == SMESH_NODE_ID_DIRICHLET ) {
		    
		    // Dirichlet
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    (*_epot)(i,j,k) = _geom.get_boundary( boundary ).value( x );

		}
	    }
	}
    }
    

}


void EpotRBGSSolver::postprocess( void )
{
    delete _rhs;
    _rhs = NULL;

    // Remove fixed vacuum tags
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID_FIX ) {
		    // Change to near solid node, keeping index pointer
		    uint32_t index = SMESH_NEAR_SOLID_INDEX_MASK & mesh;
		    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | index;
		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX ) {
		    // Change to vacuum node
		    _geom.mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM;
		}
	    }
	}
    }

    // Change Neumann boundaries next to solid nodes back to
    // NEAR_SOLID and retrieve stored solid indexes. Take care to
    // process in the same order as in Geometry class (x overrides y,
    // which overrides z.
    uint32_t near_solid_index = 0;

    // Xmin and Xmax
    for( uint32_t bound = 1; bound <= 2; bound++ ) {
	uint32_t i = 0;
	if( bound == 2 ) i = _geom.size(0)-1;
	if( _geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
	    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
		for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		    uint32_t mesh = _geom.mesh(i,j,k);
		    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		    if( node_id == SMESH_NODE_ID_NEUMANN && _geom.is_near_solid(i,j,k) ) {
			uint32_t index = _nsind[near_solid_index++];
			_geom.mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | index;
		    }
		}
	    }
	}
    }
    if( _geom.geom_mode() == MODE_2D || _geom.geom_mode() == MODE_CYL ||
	_geom.geom_mode() == MODE_3D ) {
	// Ymin and Ymax
	for( uint32_t bound = 3; bound <= 4; bound++ ) {
	    uint32_t j = 0;
	    if( bound == 4 ) j = _geom.size(1)-1;
	    if( _geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
		for( uint32_t k = 0; k < _geom.size(2); k++ ) {
		    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
			uint32_t mesh = _geom.mesh(i,j,k);
			uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
			if( node_id == SMESH_NODE_ID_NEUMANN && _geom.is_near_solid(i,j,k) ) {
			    uint32_t index = _nsind[near_solid_index++];
			    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | index;
			}
		    }
		}
	    }
	}
    }
    if( _geom.geom_mode() == MODE_3D ) {
	// Zmin and Zmax
	for( uint32_t bound = 5; bound <= 6; bound++ ) {
	    uint32_t k = 0;
	    if( bound == 6 ) k = _geom.size(2)-1;
	    if( _geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
		for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
			uint32_t mesh = _geom.mesh(i,j,k);
			uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
			if( node_id == SMESH_NODE_ID_NEUMANN && _geom.is_near_solid(i,j,k) ) {
			    uint32_t index = _nsind[near_solid_index++];
			    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | index;
			}
		    }
		}
	    }
	}
    }

    // Clear near solid indexes vector
    _nsind.clear();
}


void EpotRBGSSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    Timer t;

	ibsimu.message( 1 ) << "  Using Red-Black Gauss-Seidel solver (" 
		  << "w = " << _w
		  << ", imax = " << _imax
		  << ", eps = " << _eps
		  << ", threads = " << ibsimu.get_thread_count()
		  << ")\n";

    // Build rhs and save fields
    _epot = &epot;
    preprocess( scharge );

    // Initialize
    if( _sp )
	delete _sp;
    //_sp = new StatusPrint( ibsimu.message( 1 ) );
    _sp = new StatusPrint;
    std::stringstream ss;
    ss << "  " << std::setw(5) << 0 << " " << std::scientific << std::setw(20) << 0;
    _sp->print( ss.str() );
    _res = 0.0;
    _done = false;
    _done_count = 0;

    // Start threads
    EpotIteratorData itdata[ibsimu.get_thread_count()];
    for( uint32_t a = 0; a < ibsimu.get_thread_count(); a++ ) {

	itdata[a]._solver = this;
	itdata[a]._thno = a;
	itdata[a]._thcount = ibsimu.get_thread_count();
	pthread_create( &itdata[a]._thread, NULL, EpotRBGSSolver::iterator_entry, (void *)&itdata[a] );
    }

    // Join threads
    for( uint32_t a = 0; a < ibsimu.get_thread_count(); a++ )
	pthread_join( itdata[a]._thread, NULL );
    if( _error )
	throw( _err );

    // Postprocess
    postprocess();

    // End timer
    t.stop();

    ss.str( "" );
    ss << "  " << std::setw(5) << _iter << " " << std::scientific << std::setw(20) << _res;
    _sp->print( ss.str(), true );
    ibsimu.message( 1 ) << "\n";
    if( _iter == _imax )
	ibsimu.message( 1 ) << "  Maximum number of iteration rounds done.\n";
    ibsimu.message( 1 ) << "  residual error = " << _res << "\n";
    ibsimu.message( 1 ) << "  iterations = " << _iter << "\n";
    ibsimu.message( 1 ) << "  time used = " << t << "\n";
    ibsimu.message( 1 ) << std::flush;

    // Free rhs and sp
    delete _rhs;
    _rhs = NULL;
    delete _sp;
    _sp = NULL;
}


void EpotRBGSSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}




void EpotRBGSSolver::debug_print( std::ostream &os ) const 
{
    EpotSolver::debug_print( os );
    os << "**EpotRBGSSolver\n";
    os << "imax = " << _imax << "\n";
    os << "eps = " << _eps << "\n";
    os << "w = " << _w << "\n";

}

