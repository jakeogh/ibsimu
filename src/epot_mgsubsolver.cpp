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


#include "epot_mgsubsolver.hpp"
#include "ibsimu.hpp"
#include "constants.hpp"
#include "compmath.hpp"


/* *****************************************************************************
 * EpotMGSubSolver
 */

EpotMGSubSolver::EpotMGSubSolver( const EpotSolver &epsolver, Geometry &geom )
    : EpotSolver( epsolver, geom )
{
}


double EpotMGSubSolver::mg_smooth( MeshScalarField *epot, const MeshScalarField *rhs )
{
    _epot = epot;
    _rhs  = rhs;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	return( rbgs_loop_1d() );
	break;
    case MODE_2D:
	return( rbgs_loop_2d() );
	break;
    default:
	break;
    }

    throw( ErrorAssert( ERROR_LOCATION ) );
}


double EpotMGSubSolver::mg_solve( MeshScalarField *epot, const MeshScalarField *rhs, double w )
{
    _epot = epot;
    _rhs  = rhs;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	return( sor_loop_1d( w ) );
	break;
    case MODE_2D:
	return( sor_loop_2d( w ) );
	break;
    default:
	break;
    }

    throw( ErrorAssert( ERROR_LOCATION ) );
}


void EpotMGSubSolver::defect( MeshScalarField *defect, MeshScalarField *epot, const MeshScalarField *rhs )
{
    _defect = defect;
    _epot   = epot;
    _rhs    = rhs;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	defect_1d();
	break;
    case MODE_2D:
	defect_2d();
	break;
    default:
	break;
    }
}


void EpotMGSubSolver::preprocess( MeshScalarField &epot )
{
    EpotSolver::preprocess( epot );
}


void EpotMGSubSolver::postprocess( void )
{
    EpotSolver::postprocess();
}





/* *****************************************************************************
 * 1D Gauss-Seidel
 */


double EpotMGSubSolver::gs_process_near_solid_1d( const uint8_t *nearsolid_ptr, 
						  uint32_t a ) const
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


double EpotMGSubSolver::gs_process_pure_vacuum_1d( uint32_t a ) const
{
    // (phi_i-1 - 2*phi_i + phi_i+1) / h^2 = rho_i/eps_0
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(a+1) + (*_epot)(a-1) 
		      - 2.0*p - (*_rhs)(a) - rhst ) / ( 2.0 + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( (*_epot)(a+1) + (*_epot)(a-1) 
		      - 2.0*p - (*_rhs)(a) - rhst ) / ( 2.0 + drhst ) );
    }

    return( (1.0/2.0) * ( (*_epot)(a+1) + (*_epot)(a-1) - (*_rhs)(a) ) );
}


double EpotMGSubSolver::gs_process_neumann_1d( uint32_t boundary, uint32_t a ) const
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
	    // (phi_i - phi_i-1) / h = q_0
	    return( (*_epot)(a-1) + (*_rhs)(a) );
	break;
    }

    return( 0.0 );
}


double EpotMGSubSolver::rbgs_loop_1d( void ) const
{
    double maxerr = 0.0;

    // Go through internal nodes once using Red-Black ordering
    for( uint32_t rb = 0; rb < 2; rb++ ) {
	    
	uint32_t i = 0;
	if( rb == 0 )
	    i = 0;
	else
	    i = 1;

	for( ; i < _geom.size(0); i+=2 ) {
	    
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
    }

    return( maxerr );
}


double EpotMGSubSolver::sor_loop_1d( double w ) const
{
    double maxerr = 0.0;
    const double w2 = 1.0-w;

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
	Vnew = w*Vnew + w2*Vold;
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
 * 1D Defect
 */


double EpotMGSubSolver::defect_near_solid_1d( const uint8_t *nearsolid_ptr, uint32_t a ) const
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

    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    }

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}


double EpotMGSubSolver::defect_pure_vacuum_1d( uint32_t a ) const
{
    // (phi_i-1 - 2*phi_i + phi_i+1) / h^2 = rho_i/eps_0
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) - 2.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) - 2.0*p - (*_rhs)(a) - rhst );
    }

    return( (*_epot)(a+1) + (*_epot)(a-1) - 2.0*(*_epot)(a) - (*_rhs)(a) );
}


double EpotMGSubSolver::defect_neumann_1d( uint32_t boundary, uint32_t a ) const
{
    switch( boundary ) {
    case 1:
	if( _neumann_order == 2 )
	    // 3*phi_i - 4*phi_i+1 + phi_i+2 = 2h*q_0
	    return( 3.0*(*_epot)(a) - 4.0*(*_epot)(a+1) + (*_epot)(a+2) - (*_rhs)(a) );
	else
	    // phi_i - phi_i+1 = h*q_0
	    return( (*_epot)(a) - (*_epot)(a+1) - (*_rhs)(a) );
	break;
    case 2:
	if( _neumann_order == 2 )
	    // 3*phi_i - 4*phi_i-1 + phi_i-2 = 2h*q_0
	    return( 3.0*(*_epot)(a) - 4.0*(*_epot)(a-1) + (*_epot)(a-2) - (*_rhs)(a) );
	else
	    // phi_i-1 - phi_i = h*q_0
	    return( (*_epot)(a) - (*_epot)(a-1) - (*_rhs)(a) );
	break;
    }

    return( 0.0 );
}


void EpotMGSubSolver::defect_1d( void ) const
{
    // Go through all nodes
    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

	uint32_t mesh = _geom.mesh(i);
	uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	double D;
	if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
	    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
	    D = defect_near_solid_1d( nearsolid_ptr, i );
	} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
	    D = defect_pure_vacuum_1d( i );
	} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
	    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
	    D = defect_neumann_1d( boundary, i );
	} else {
	    D = 0.0;
	}
	(*_defect)(i) = D;
    }
}


/* *****************************************************************************
 * 2D Gauss-Seidel
 */


double EpotMGSubSolver::gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
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


double EpotMGSubSolver::gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const
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


double EpotMGSubSolver::gs_process_neumann_2d( uint32_t boundary, uint32_t a, uint32_t dj ) const
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


double EpotMGSubSolver::rbgs_loop_2d( void ) const
{
    double maxerr = 0.0;

    // Go through all nodes once using Red-Black ordering
    for( uint32_t rb = 0; rb < 2; rb++ ) {
	const uint32_t dj = _geom.size(0);
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    
	    uint32_t i = j % 2;
	    if( i == rb )
		i = 0;
	    else
		i = 1;

	    for( ; i < _geom.size(0); i+=2 ) {
	    
		uint32_t a = j*dj+i;
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
    }

    // Return largest change in any node
    return( maxerr );
}


double EpotMGSubSolver::sor_loop_2d( double w ) const
{
    double maxerr = 0.0;
    const double w2 = 1.0-w;
    const uint32_t dj = _geom.size(0);

    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	for( uint32_t i = 0; i < _geom.size(0); i++ ) {
	    
	    uint32_t a = j*dj+i;
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
	    Vnew = w*Vnew + w2*Vold;
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

    // Return largest change in any node
    return( maxerr );
}


/* *****************************************************************************
 * 2D Defect
 */


double EpotMGSubSolver::defect_near_solid_2d( const uint8_t *nearsolid_ptr, 
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
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    }

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}



double EpotMGSubSolver::defect_pure_vacuum_2d( uint32_t a, uint32_t dj ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) + (*_epot)(a+dj) 
		+ (*_epot)(a-dj) - 4.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) + (*_epot)(a+dj) 
		+ (*_epot)(a-dj) - 4.0*p - (*_rhs)(a) - rhst );
    }

    return( (*_epot)(a+1) + (*_epot)(a-1) +
	    (*_epot)(a+dj) + (*_epot)(a-dj) - 4.0*(*_epot)(a) - (*_rhs)(a) );
}


double EpotMGSubSolver::defect_neumann_2d( uint32_t boundary, uint32_t a, uint32_t dj ) const
{
    switch( boundary ) {
    case 1:
	if( _neumann_order == 2 )
	    // 3*phi_i - 4*phi_i+1 + phi_i+2 = 2h*q_0
	    return( 3.0*(*_epot)(a) - 4.0*(*_epot)(a+1) + (*_epot)(a+2) - (*_rhs)(a) );
	else
	    // phi_i - phi_i+1 = h*q_0
	    return( (*_epot)(a) - (*_epot)(a+1) - (*_rhs)(a) );
	break;
    case 2:
	if( _neumann_order == 2 )
	    // 3*phi_i - 4*phi_i-1 + phi_i-2 = 2h*q_0
	    return( 3.0*(*_epot)(a) - 4.0*(*_epot)(a-1) + (*_epot)(a-2) - (*_rhs)(a) );
	else
	    // phi_i-1 - phi_i = h*q_0
	    return( (*_epot)(a) - (*_epot)(a-1) - (*_rhs)(a) );
	break;
    case 3:
	if( _neumann_order == 2 )
	    // 3*phi_j - 4*phi_j+1 + phi_j+2 = 2h*q_0
	    return( 3.0*(*_epot)(a) - 4.0*(*_epot)(a+dj) + (*_epot)(a+2*dj) - (*_rhs)(a) );
	else
	    // phi_j - phi_j+1 = h*q_0
	    return( (*_epot)(a) - (*_epot)(a+dj) - (*_rhs)(a) );
	break;
    case 4:
	if( _neumann_order == 2 )
	    // 3*phi_j - 4*phi_j-1 + phi_j-2 = 2h*q_0
	    return( (*_epot)(a) - 4.0*(*_epot)(a-dj) + (*_epot)(a-2*dj) - (*_rhs)(a) );
	else
	    // phi_j-1 - phi_j = h*q_0
	    return( (*_epot)(a) - (*_epot)(a-dj) - (*_rhs)(a) );
	break;
    }

    return( 0.0 );
}


void EpotMGSubSolver::defect_2d( void ) const
{
    // Go through all nodes
    const uint32_t dj = _geom.size(0);
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	for( uint32_t i = 0; i < _geom.size(0); i++ ) {

	    uint32_t a = j*dj+i;
	    uint32_t mesh = _geom.mesh(a);
	    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	    double D;
	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		D = defect_near_solid_2d( nearsolid_ptr, a, dj );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		D = defect_pure_vacuum_2d( a, dj );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		D = defect_neumann_2d( boundary, a, dj );
	    } else {
		D = 0.0;
	    }
	    (*_defect)(a) = D;
	}
    }
}



