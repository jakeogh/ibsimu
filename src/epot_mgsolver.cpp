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




/* *****************************************************************************
 * EpotMGSubSolver
 */

EpotMGSubSolver::EpotMGSubSolver( const EpotSolver &epsolver, Geometry &geom )
    : EpotSolver( epsolver, geom )
{
}


double EpotMGSubSolver::mg_relax( MeshScalarField *epot, const MeshScalarField *rhs, double w )
{
    _epot = epot;
    _rhs  = rhs;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	if( w == 1.0 )
	    return( rbgs_loop_1d() );
	else
	    return( sor_loop_1d( w ) );
	break;
    case MODE_2D:
	if( w == 1.0 )
	    return( rbgs_loop_2d() );
	else
	    return( sor_loop_2d( w ) );
	break;
    default:
	break;
    }

    return( 0.0 );
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

    return( (1.0/cof) * ( epf - (*_rhs)(i) ) );
}


double EpotMGSubSolver::gs_process_pure_vacuum_1d( uint32_t i ) const
{
    // (phi_i-1 - 2*phi_i + phi_i+1) / h^2 = rho_i/eps_0
    //return( (1.0/2.0) * ( (*_epot)(i+1) + (*_epot)(i-1) - _epot->h()*_epot->h()*(*_rhs)(i) ) );
    return( (1.0/2.0) * ( (*_epot)(i+1) + (*_epot)(i-1) - (*_rhs)(i) ) );
}


double EpotMGSubSolver::gs_process_neumann_1d( uint32_t boundary, uint32_t i ) const
{
    switch( boundary ) {
    case 1:
	if( _neumann_order == 2 )
	    // (3*phi_i - 4*phi_i+1 + phi_i+2) / 2h = q_0
	    //return( (4.0*(*_epot)(i+1) - (*_epot)(i+2) + 2.0*_epot->h()*(*_rhs)(i))/3.0 );
	    return( (4.0*(*_epot)(i+1) - (*_epot)(i+2) + (*_rhs)(i))/3.0 );
	else
	    // (phi_i - phi_i+1) / h = q_0
	    return( (*_epot)(i+1) + (*_rhs)(i) );
	break;
    case 2:
	if( _neumann_order == 2 )
	    // (3*phi_i - 4*phi_i-1 + phi_i-2) / 2h = q_0
	    //return( (4.0*(*_epot)(i-1) - (*_epot)(i-2) + 2.0*_epot->h()*(*_rhs)(i))/3.0 );
	    return( (4.0*(*_epot)(i-1) - (*_epot)(i-2) + (*_rhs)(i))/3.0 );
	else
	    // (phi_i - phi_i-1) / h = q_0
	    return( (*_epot)(i-1) + (*_rhs)(i) );
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


// D=A*X-B
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

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}


// D=A*X-B
double EpotMGSubSolver::defect_pure_vacuum_1d( uint32_t a ) const
{
    // (phi_i-1 - 2*phi_i + phi_i+1) / h^2 = rho_i/eps_0
    //return( ((*_epot)(a+1) + (*_epot)(a-1) - 2.0*(*_epot)(a))/(_epot->h()*_epot->h()) - (*_rhs)(a) );
    return( (*_epot)(a+1) + (*_epot)(a-1) - 2.0*(*_epot)(a) - (*_rhs)(a) );
}


// D=A*X-B
double EpotMGSubSolver::defect_neumann_1d( uint32_t boundary, uint32_t a ) const
{
    switch( boundary ) {
    case 1:
	if( _neumann_order == 2 )
	    // 3*phi_i - 4*phi_i+1 + phi_i+2 = 2h*q_0
	    //return( (3.0*(*_epot)(a) - 4.0*(*_epot)(a+1) + (*_epot)(a+2))/(2.0*_epot->h()) - (*_rhs)(a) );
	    return( 3.0*(*_epot)(a) - 4.0*(*_epot)(a+1) + (*_epot)(a+2) - (*_rhs)(a) );
	else
	    // phi_i - phi_i+1 = h*q_0
	    return( (*_epot)(a) - (*_epot)(a+1) - (*_rhs)(a) );
	break;
    case 2:
	if( _neumann_order == 2 )
	    // 3*phi_i - 4*phi_i-1 + phi_i-2 = 2h*q_0
	    //return( (3.0*(*_epot)(a) - 4.0*(*_epot)(a-1) + (*_epot)(a-2))/(2.0*_epot->h()) - (*_rhs)(a) );
	    return( 3.0*(*_epot)(a) - 4.0*(*_epot)(a-1) + (*_epot)(a-2) - (*_rhs)(a) );
	else
	    // phi_i-1 - phi_i = h*q_0
	    return( (*_epot)(a) - (*_epot)(a-1) - (*_rhs)(a) );
	break;
    }

    return( 0.0 );
}


// D=A*X-B
void EpotMGSubSolver::defect_1d( void ) const
{
    // Go through all nodes
    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

	uint32_t mesh = _geom.mesh(i);
	uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	double D;
	if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
	    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
	    D = 4.0*defect_near_solid_1d( nearsolid_ptr, i );
	} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
	    D = 4.0*defect_pure_vacuum_1d( i );
	} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
	    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
	    D = 2.0*defect_neumann_1d( boundary, i );
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


// D=A*X-B
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
	return( p + ( epf - cof*p - (*_rhs)(a) - rhst ) / ( cof + drhst ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( p + ( epf - cof*p - (*_rhs)(a) - rhst ) / ( cof + drhst ) );	
    }

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}


// D=A*X-B
double EpotMGSubSolver::defect_pure_vacuum_2d( uint32_t a, uint32_t dj ) const
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

    return( (*_epot)(a+1) + (*_epot)(a-1) +
	    (*_epot)(a+dj) + (*_epot)(a-dj) - 4.0*(*_epot)(a) - (*_rhs)(a) );
}


// D=A*X-B
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


// D=A*X-B
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
		D = 4.0*defect_near_solid_2d( nearsolid_ptr, a, dj );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		D = 4.0*defect_pure_vacuum_2d( a, dj );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		D = 2.0*defect_neumann_2d( boundary, a, dj );
	    } else {
		D = 0.0;
	    }
	    (*_defect)(a) = D;
	}
    }
}


/* *****************************************************************************
 * EpotMGSolver
 */


EpotMGSolver::EpotMGSolver( Geometry &geom )
    : EpotSolver( geom ), _geom_prepared(false), _levels(1), _npre(5), 
      _npost(5), _ncyc(1), _res(0.0), _eps(1.0e-6), _imax(10000)
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

void EpotMGSolver::set_imax( uint32_t imax )
{
    _imax = imax;
}

double EpotMGSolver::get_residual( void ) const
{
    return( _res );
}


void EpotMGSolver::set_levels( uint32_t levels )
{
    if( levels < 1 )
	throw( Error( ERROR_LOCATION, "Invalid number of multigrid levels" ) ); 
    _levels = levels;
    reset_problem();
}


void EpotMGSolver::set_ncyc( uint32_t ncyc )
{
    if( ncyc < 1 )
	throw( Error( ERROR_LOCATION, "Invalid number of multigrid cycles" ) ); 
    _ncyc = ncyc;
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
    std::cout << "  Preparing multi grid geometries\n";

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

	    // Make sub-solver
	    EpotMGSubSolver *mgss = new EpotMGSubSolver( *this, *_geomv.back() );
	    _epotsolverv.push_back( mgss );
	}
    }
}


void EpotMGSolver::preprocess( MeshScalarField &epot, const MeshScalarField &scharge )
{
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
    std::cout << "  Postprocessing\n";

    for( uint32_t a = 0; a < _levels; a++ ) {

	if( a != 0 )
	    delete _epotv[a];
	delete _workv[a];
	delete _rhsv[a];

	// Postprocess solid meshes
	_epotsolverv[a]->postprocess();
    }

    _epotv.clear();
    _rhsv.clear();
    _workv.clear();
}


// Calculate defect to workfield for level, D=A*X-B
void EpotMGSolver::defect( int level )
{
    std::cout << "  Calculating defect for level " << level << "\n";
    _epotsolverv[level]->defect( _workv[level], _epotv[level], _rhsv[level] );
}


// Restrict defect from level (work field) to level+1 (rhs)
//
// Also negates defect
//
// 2D symbol is:
// 1/16   1/8   1/16
//  1/8   1/4   1/8
// 1/16   1/8   1/16
//
// Return 0.0 contribution from outside mesh
void EpotMGSolver::restrict_2d( int level )
{
    std::cout << "  Restricing defect from level " << level << " to level " << level+1 << "\n";

    // Go through internal nodes of level+1
    int32_t s = _geomv[level+1]->size(0)-1;
    int32_t t = _geomv[level+1]->size(1)-1;
    for( int32_t j = 1; j < t; j++ ) {
        for( int32_t i = 1; i < s; i++ ) {

	    int32_t ii = 2*i;
	    int32_t jj = 2*j;
	    double D = 
		1.0/16.0*(*_workv[level])( ii-1, jj-1 ) + 
		1.0/8.0 *(*_workv[level])( ii-1, jj   ) + 
		1.0/16.0*(*_workv[level])( ii-1, jj+1 ) + 

		1.0/8.0 *(*_workv[level])( ii,   jj-1 ) + 
		1.0/4.0 *(*_workv[level])( ii,   jj   ) + 
		1.0/8.0 *(*_workv[level])( ii,   jj+1 ) + 

		1.0/16.0*(*_workv[level])( ii+1, jj-1 ) + 
		1.0/8.0 *(*_workv[level])( ii+1, jj   ) + 
		1.0/16.0*(*_workv[level])( ii+1, jj+1 );

	    (*_rhsv[level+1])(i,j) = -D;
	}
    }

    double fac = -1.0;
    if( _neumann_order == 2 )
	fac = -0.5;

    // i=0 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*_workv[level])( ii, jj-1 ) + 
	    1.0/2.0*(*_workv[level])( ii, jj   ) + 
	    1.0/4.0*(*_workv[level])( ii, jj+1 );
	
	(*_rhsv[level+1])(i,j) = fac*D;
    }

    // i=size(0)-1 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = s;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*_workv[level])( ii, jj-1 ) + 
	    1.0/2.0*(*_workv[level])( ii, jj   ) + 
	    1.0/4.0*(*_workv[level])( ii, jj+1 );
	
	(*_rhsv[level+1])(i,j) = fac*D;
    }

    // j=0 boundary nodes of level+1
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*_workv[level])( ii-1, jj ) + 
	    1.0/2.0*(*_workv[level])( ii,   jj ) + 
	    1.0/4.0*(*_workv[level])( ii+1, jj );
	
	(*_rhsv[level+1])(i,j) = fac*D;
    }

    // j=size(1)-1 boundary nodes of level+1
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = t;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*_workv[level])( ii-1, jj ) + 
	    1.0/2.0*(*_workv[level])( ii,   jj ) + 
	    1.0/4.0*(*_workv[level])( ii+1, jj );
	
	(*_rhsv[level+1])(i,j) = fac*D;
    }

    // And corners
    int32_t i = 0;
    int32_t j = 0;
    double D = (*_workv[level])( 2*i, 2*j );
    (*_rhsv[level+1])(i,j) = fac*D;

    i = s;
    D = (*_workv[level])( 2*i, 2*j );
    (*_rhsv[level+1])(i,j) = fac*D;

    j = t;
    D = (*_workv[level])( 2*i, 2*j );
    (*_rhsv[level+1])(i,j) = fac*D;

    i = 0;
    D = (*_workv[level])( 2*i, 2*j );
    (*_rhsv[level+1])(i,j) = fac*D;
}


void EpotMGSolver::restrict_1d( int level )
{
    std::cout << "  Restricing defect from level " << level << " to level " << level+1 << "\n";

    // Go through internal nodes of level+1
    int32_t s = _geomv[level+1]->size(0)-1;
    for( int32_t i = 1; i < s; i++ ) {

	int32_t ii = 2*i;
	double D = 
	    1.0/4.0*(*_workv[level])( ii-1 ) + 
	    1.0/2.0*(*_workv[level])( ii   ) + 
	    1.0/4.0*(*_workv[level])( ii+1 );
	
	(*_rhsv[level+1])(i) = -D;
    }

    // Boundary nodes completely separated
    if( _neumann_order == 1 ) {
	(*_rhsv[level+1])(0) = -(*_workv[level])(0);
	(*_rhsv[level+1])(s) = -(*_workv[level])(2*s);
    } else {
	(*_rhsv[level+1])(0) = -0.5*(*_workv[level])(0);
	(*_rhsv[level+1])(s) = -0.5*(*_workv[level])(2*s);
    }
}


void EpotMGSolver::restrict( int level )
{
    switch( _geom.geom_mode() ) {
    case MODE_1D:
	restrict_1d( level );
	break;
    case MODE_2D:
	restrict_2d( level );
	break;
    default:
	break;
    }
}


void EpotMGSolver::prolong_add_2d( int level, int32_t i, int32_t j, double C )
{
    if( i >= 0 && i < (int32_t)_geomv[level]->size(0) &&
	j >= 0 && j < (int32_t)_geomv[level]->size(1) )
	(*_workv[level])(i,j) += C;
}


// Prolong correction (epot) to work field from level+1 to level
//
// 2D symbol is:
// 1/4   1/2   1/4
// 1/2    1    1/2
// 1/4   1/2   1/4
//
void EpotMGSolver::prolong_2d( int level )
{
    std::cout << "  Prolonging 2d correction from level " << level+1 << " to level " << level << "\n";

    // Clear work field
    _workv[level]->clear();

    // Loop through all level+1 nodes
    for( int32_t j = 0; j < (int32_t)_geomv[level+1]->size(1); j++ ) {
	for( int32_t i = 0; i < (int32_t)_geomv[level+1]->size(0); i++ ) {

	    int32_t ii = 2*i;
	    int32_t jj = 2*j;
	    double C = (*_epotv[level+1])(i,j);
	    prolong_add_2d( level, ii-1, jj-1, 1.0/4.0*C );
	    prolong_add_2d( level, ii-1, jj,   1.0/2.0*C );
	    prolong_add_2d( level, ii-1, jj+1, 1.0/4.0*C );

	    prolong_add_2d( level, ii,   jj-1, 1.0/2.0*C );
	    prolong_add_2d( level, ii,   jj,   1.0*C );
	    prolong_add_2d( level, ii,   jj+1, 1.0/2.0*C );

	    prolong_add_2d( level, ii+1, jj-1, 1.0/4.0*C );
	    prolong_add_2d( level, ii+1, jj,   1.0/2.0*C );
	    prolong_add_2d( level, ii+1, jj+1, 1.0/4.0*C );
	}
    }    
}


void EpotMGSolver::prolong_add_1d( int level, int32_t i, double C )
{
    if( i >= 0 && i < (int32_t)_geomv[level]->size(0) )
	(*_workv[level])(i) += C;
}


void EpotMGSolver::prolong_1d( int level )
{
    std::cout << "  Prolonging 1d correction from level " << level+1 << " to level " << level << "\n";

    // Clear work field
    _workv[level]->clear();

    // Loop through all level+1 nodes
    for( int32_t i = 0; i < (int32_t)_geomv[level+1]->size(0); i++ ) {
	
	int32_t ii = 2*i;
	double C = (*_epotv[level+1])(i);
	prolong_add_1d( level, ii-1, 1.0/2.0*C );
	prolong_add_1d( level, ii,       1.0*C );
	prolong_add_1d( level, ii+1, 1.0/2.0*C );
    }    
}


void EpotMGSolver::prolong( int level )
{
    switch( _geom.geom_mode() ) {
    case MODE_1D:
	prolong_1d( level );
	break;
    case MODE_2D:
	prolong_2d( level );
	break;
    default:
	break;
    }
}

// Make correction Xnew=X+V
void EpotMGSolver::correct( int level )
{
    std::cout << "  Calculating correction for level " << level << "\n";

    // Loop through all level nodes, only correct non-fixed nodes
    for( int32_t j = 0; j < (int32_t)_geomv[level]->size(1); j++ ) {
	for( int32_t i = 0; i < (int32_t)_geomv[level]->size(0); i++ ) {

	    if( !(_geomv[level]->mesh(i,j) & SMESH_NODE_FIXED) )
		(*_epotv[level])(i,j) += (*_workv[level])(i,j);
	}
    }
}


void EpotMGSolver::mg_recurse( uint32_t level )
{
    std::cout << "MG recursion, level " << level << "\n";

    if( level == _levels-1 ) {
 
	std::cout << "  Roughest level\n";

	// Last level, solve the roughest problem until convergence
	uint32_t a = 0;
	double res = 0.0;
	for( a = 1; a <= _imax; a++ ) {
	    res = _epotsolverv[level]->mg_relax( _epotv[level], _rhsv[level], 1.7 );
	    if( res < _eps )
		break;
	}

	//std::cout << "epot (level = " << level << ") accurate solve:\n";
	//print_field( _epotv[level] );

	std::cout << "  " << a << " iterations done\n";
	std::cout << "  " << res << " accuracy reached\n";
	if( a == _imax )
	    std::cout << "  maximum number of iterations done\n";

	return;
    } 

    // Do ncyc cycles of next level
    double res = 0.0;
    for( uint32_t a = 0; a < _ncyc; a++ ) {

	std::cout << "  Doing cycle " << a+1 << "/" << _ncyc << "\n";

	//std::cout << "epot (level = " << level << "):\n";
	//print_field( _epotv[level] );

	// Pre smoothing
	for( uint32_t a = 0; a < _npre; a++ )
	    _epotsolverv[level]->mg_relax( _epotv[level], _rhsv[level] );

	//std::cout << "epot (level = " << level << ") pre smoothed:\n";
	//print_field( _epotv[level] );

	// Calculate defect to work field, D=A*X-B
	defect( level );

	//std::cout << "defect (level = " << level << "):\n";
	//print_field( _workv[level] );

	// Restrict defect (work) to next level (rhs)
	restrict( level );

	//std::cout << "defect (level = " << level+1 << "):\n";
	//print_field( _rhsv[level+1] );

	// Recurse to next level to solve for correction V (epot)
	// A*V=-D
	mg_recurse( level+1 );

	//std::cout << "correction (level = " << level+1 << "):\n";
	//print_field( _epotv[level+1] );

	// Prolong correction (epot) to work field
	prolong( level );

	//std::cout << "correction (level = " << level << "):\n";
	//print_field( _workv[level] );

	// Make correction Xnew=X+V
	correct( level );

	//std::cout << "epot (level = " << level << "):\n";
	//print_field( _epotv[level] );

	// Post smoothing
	for( uint32_t a = 0; a < _npost; a++ )
	    res = _epotsolverv[level]->mg_relax( _epotv[level], _rhsv[level] );

	//std::cout << "epot (level = " << level << ") post smoothed:\n";
	//print_field( _epotv[level] );

	std::cout << "  " << res << " accuracy reached\n";
     }

    // Save residual from last level
    _res = res;
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
    if( ibsimu.get_verbose_output() ) {
	std::cout << "  Using Multigrid solver (" 
		  << "levels = " << _levels
		  << ", npre = " << _npre
		  << ", npost = " << _npost
		  << ", ncyc = " << _ncyc
		  << ", eps = " << _eps
		  << ", imax = " << _imax
		  << ")\n";
    }

    if( !_geom_prepared )
	prepare_mg_geom();
    preprocess( epot, scharge );

    // Do first level iteration only once
    mg_recurse( 0 );

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

