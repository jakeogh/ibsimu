/*! \file epot_gssolver.cpp
 *  \brief Gauss-Seidel solver for electric potential problem
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


#include "epot_gssolver.hpp"
#include "ibsimu.hpp"
#include "constants.hpp"
#include "compmath.hpp"
#include "statusprint.hpp"


EpotGSSolver::EpotGSSolver( Geometry &geom )
    : EpotSolver( geom ), _epot(NULL), _rhs(NULL), _iter(0), _imax(10000), 
      _eps(1.0e-6), _res(0.0), _w(1.66)
{
    
}


EpotGSSolver::EpotGSSolver( Geometry &geom, std::istream &s )
    : EpotSolver(geom,s)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void EpotGSSolver::reset_problem( void )
{
    // Do nothing
}


void EpotGSSolver::set_eps( double eps )
{
    _eps = eps;
}


double EpotGSSolver::get_residual( void ) const
{
    return( _res );
}


uint32_t EpotGSSolver::get_iter( void ) const
{
    return( _iter );
}


void EpotGSSolver::set_imax( uint32_t imax )
{
    _imax = imax;
}


void EpotGSSolver::set_w( double w )
{
    _w = w;
}


/* *****************************************************************************
 * 3D
 */

double EpotGSSolver::gs_process_near_solid_3d( const uint8_t *nearsolid_ptr, uint32_t a,
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


double EpotGSSolver::gs_process_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(a+1) + (*_epot)(a-1) 
		      + (*_epot)(a+dj) + (*_epot)(a-dj) 
		      + (*_epot)(a+dk) + (*_epot)(a-dk) - 6.0*p
		      - (*_rhs)(a) - rhst ) / ( 6.0 + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(a+1) + (*_epot)(a-1) 
		      + (*_epot)(a+dj) + (*_epot)(a-dj) 
		      + (*_epot)(a+dk) + (*_epot)(a-dk) - 6.0*p
		      - (*_rhs)(a) - rhst ) / ( 6.0 + drhst ) );
    }

    return( (1.0/6.0) * ( (*_epot)(a+1)  + (*_epot)(a-1) +
			  (*_epot)(a+dj) + (*_epot)(a-dj) + 
			  (*_epot)(a+dk) + (*_epot)(a-dk) - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_process_neumann_3d( uint32_t boundary, uint32_t a,
					    uint32_t dj, uint32_t dk ) const
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
    case 5:
	if( _neumann_order == 2 )
	    // (3*phi_k - 4*phi_k+1 + phi_k+2) / 2h = q_0
	    return( (4.0*(*_epot)(a+dk) - (*_epot)(a+2*dk) + (*_rhs)(a))/3.0 );
	else
	    // (phi_k - phi_k+1) / h = q_0
	    return( (*_epot)(a+dk) + (*_rhs)(a) );
	break;
    case 6:
	if( _neumann_order == 2 )
	    // (3*phi_k - 4*phi_k-1 + phi_k-2) / 2h = q_0
	    return( (4.0*(*_epot)(a-dk) - (*_epot)(a-2*dk) + (*_rhs)(a))/3.0 );
	else
	    // (phi_k-1 - phi_k) / h = q_0
	    return( (*_epot)(a-dk) + (*_rhs)(a) );
	break;
    }

    return( 0.0 );
}


double EpotGSSolver::gs_loop_3d( void ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    const double w2 = 1.0-_w;
    const uint32_t dj = _geom.size(0);
    const uint32_t dk = _geom.size(0)*_geom.size(1);

    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    uint32_t a = k*dk+j*dj;
	    for( uint32_t i = 0; i < _geom.size(0); i++, a++ ) {

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
		    // Dirichlet or otherwise fixed
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


/* *****************************************************************************
 * CYL
 */

double EpotGSSolver::gs_process_near_solid_cyl( const uint8_t *nearsolid_ptr, 
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
    double epf = 2.0/(alpha+beta)*( (*_epot)(i-1,j)/alpha + (*_epot)(i+1,j)/beta );

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
    epf += 1.0/(alpha+beta)*( (2.0/alpha-1.0/j)*(*_epot)(i,j-1) + (2.0/beta+1.0/j)*(*_epot)(i,j+1) );

    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( p + ( epf - cof*p - (*_rhs)(i,j) - rhst ) / ( cof + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( epf - cof*p - (*_rhs)(i,j) - rhst ) / ( cof + drhst ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(i,j) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_cyl( uint32_t i, uint32_t j ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(i+1,j) + (*_epot)(i-1,j) 
		      + (1.0+0.5/j)*(*_epot)(i,j+1) 
		      + (1.0-0.5/j)*(*_epot)(i,j-1) - 4.0*p
		      - (*_rhs)(i,j) - rhst ) / ( 4.0 + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(i+1,j) + (*_epot)(i-1,j) 
		      + (1.0+0.5/j)*(*_epot)(i,j+1) 
		      + (1.0-0.5/j)*(*_epot)(i,j-1) - 4.0*p
		      - (*_rhs)(i,j) - rhst ) / ( 4.0 + drhst ) );
    }
 
    return( (1.0/4.0) * ( (*_epot)(i+1,j) + (*_epot)(i-1,j)
			  + (1.0+0.5/j)*(*_epot)(i,j+1)
			  + (1.0-0.5/j)*(*_epot)(i,j-1) 
			  - (*_rhs)(i,j) ) );
}


double EpotGSSolver::gs_process_neumann_cyl( uint32_t boundary, uint32_t a, uint32_t dj ) const
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
	// Special axis condition: (phi_{i-1,0} + phi_{i+1,j} + 4phi_{i,j+1} - 6phi_{i,j}) = rhs
	if( _plasma == PLASMA_PEXP ) {
	    double p = (*_epot)(a);
	    double rhst, drhst;
	    pexp_newton( rhst, drhst, p );
	    return( p + ( (*_epot)(a+1) + (*_epot)(a-1) 
			  + 4.0*(*_epot)(a+dj) - 6.0*p
			  - (*_rhs)(a) - rhst ) / ( 6.0 + drhst ) );
	}  else if( _plasma == PLASMA_NSIMP ) {
	    double p = (*_epot)(a);
	    double rhst, drhst;
	    nsimp_newton( rhst, drhst, p );
	    return( p + ( (*_epot)(a+1) + (*_epot)(a-1) 
			  + 4.0*(*_epot)(a+dj) - 6.0*p
			  - (*_rhs)(a) - rhst ) / ( 6.0 + drhst ) );
	} else {
	    return( (1.0/6.0) * ( (*_epot)(a-1) + (*_epot)(a+1) + 4.0*(*_epot)(a+dj) - (*_rhs)(a) ) );
	}
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


double EpotGSSolver::gs_loop_cyl( void ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    const double w2 = 1.0-_w;
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
		Vnew = gs_process_near_solid_cyl( nearsolid_ptr, i, j );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		Vnew = gs_process_pure_vacuum_cyl( i, j );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		Vnew = gs_process_neumann_cyl( boundary, a, dj );
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
 * 2D
 */

double EpotGSSolver::gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
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


double EpotGSSolver::gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const
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


double EpotGSSolver::gs_process_neumann_2d( uint32_t boundary, uint32_t a, uint32_t dj ) const
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


double EpotGSSolver::gs_loop_2d( void ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    const double w2 = 1.0-_w;
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
	    Vnew = _w*Vnew + w2*Vold;
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
 * 1D
 */

double EpotGSSolver::gs_process_near_solid_1d( const uint8_t *nearsolid_ptr, 
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
    double epf = 2.0/(alpha+beta)*( (*_epot)(i-1)/alpha + (*_epot)(i+1)/beta );

    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(i);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( p + ( epf - cof*p - (*_rhs)(i) - rhst ) / ( cof + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(i);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( epf - cof*p - (*_rhs)(i) - rhst ) / ( cof + drhst ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(i) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_1d( uint32_t i ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(i);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(i+1) + (*_epot)(i-1) 
		      - 2.0*p - (*_rhs)(i) - rhst ) / ( 2.0 + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(i);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(i+1) + (*_epot)(i-1) 
		      - 2.0*p - (*_rhs)(i) - rhst ) / ( 2.0 + drhst ) );
    }

    return( (1.0/2.0) * ( (*_epot)(i+1) + (*_epot)(i-1) - (*_rhs)(i) ) );
}


double EpotGSSolver::gs_process_neumann_1d( uint32_t boundary, uint32_t i ) const
{
    switch( boundary ) {
    case 1:
	if( _neumann_order == 2 )
	    // (3*phi_i - 4*phi_i+1 + phi_i+2) / 2h = q_0
	    return( (4.0*(*_epot)(i+1) - (*_epot)(i+2) + (*_rhs)(i))/3.0 );
	else
	    // (phi_i - phi_i+1) / h = q_0
	    return( (*_epot)(i+1) + (*_rhs)(i) );
	break;
    case 2:
	if( _neumann_order == 2 )
	    // (3*phi_i - 4*phi_i-1 + phi_i-2) / 2h = q_0
	    return( (4.0*(*_epot)(i-1) - (*_epot)(i-2) + (*_rhs)(i))/3.0 );
	else
	    // (phi_i-1 - phi_i) / h = q_0
	    return( (*_epot)(i-1) + (*_rhs)(i) );
	break;
    }

    return( 0.0 );
}


double EpotGSSolver::gs_loop_1d( void ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    double w2 = 1.0-_w;
    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
	    
	double Vold = (*_epot)(i);
	double Vnew;
	uint32_t mesh = _geom.mesh(i);
	uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
	    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
	    Vnew = gs_process_near_solid_1d( nearsolid_ptr, i );
	} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
	    Vnew = gs_process_pure_vacuum_1d( i );
	} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
	    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
	    Vnew = gs_process_neumann_1d( boundary, i );
	} else {
	    // Dirichlet
	    continue;
	}
	Vnew = _w*Vnew + w2*Vold;
	(*_epot)(i) = Vnew;
	double err = fabs( Vnew - Vold );
	if( err > maxerr )
	    maxerr = err;
	if( comp_isinf(err) ) {
	    throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) ) );
	} else if( comp_isnan(err) ) {
	    throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) ) );
	}
    }

    return( maxerr );
}


/* *****************************************************************************
 * Common
 */


void EpotGSSolver::preprocess( const MeshScalarField &scharge )
{
    EpotSolver::preprocess( *_epot );

    // Build right-hand-side
    if( _rhs )
	delete _rhs;
    _rhs = new MeshScalarField( (const Mesh)scharge );

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


void EpotGSSolver::postprocess( void )
{
    delete _rhs;
    _rhs = NULL;

    EpotSolver::postprocess();
}


void EpotGSSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    StatusPrint sp( std::cout );
    if( ibsimu.get_verbose_output() ) {
	std::cout << "  Using Gauss-Seidel solver (" 
		  << "w = " << _w
		  << ", imax = " << _imax
		  << ", eps = " << _eps
		  << ")\n";
	std::stringstream ss;
	ss << "  " << std::setw(5) << 0 << " " << std::scientific << std::setw(20) << 0;
	sp.print( ss.str() );
    }

    // Set epot pointer and preprocess
    _epot = &epot;
    preprocess( scharge );

    // Loop until converged
    _iter = 0;
    while( _iter < _imax ) {
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

	if( ibsimu.get_verbose_output() ) {
	    std::stringstream ss;
	    ss << "  " << std::setw(5) << _iter << " " << std::scientific << std::setw(20) << _res;
	    sp.print( ss.str() );
	}

	_iter++;
	if( _res < _eps )
	    break;
	if( comp_isinf(_res) || comp_isnan(_res) )
	    break;
    }

    // Postprocess
    postprocess();

    if( ibsimu.get_verbose_output() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << _iter << " " << std::scientific << std::setw(20) << _res;
	sp.print( ss.str(), true );
	std::cout << "\n";
	if( _iter == _imax )
	    std::cout << "  Maximum number of iteration rounds done.\n";
	std::cout << "  residual error = " << _res << "\n";
	std::cout << "  iterations = " << _iter << "\n";
    }
}


void EpotGSSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}




void EpotGSSolver::debug_print( std::ostream &os ) const 
{
    EpotSolver::debug_print( os );
    os << "**EpotGSSolver\n";
    os << "imax = " << _imax << "\n";
    os << "eps = " << _eps << "\n";
    os << "w = " << _w << "\n";

}

