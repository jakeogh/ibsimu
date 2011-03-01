#include "epot_gssolver.hpp"
#include "ibsimu.hpp"
#include "timer.hpp"
#include "constants.hpp"
#include "compmath.hpp"
#include "statusprint.hpp"


EpotGSSolver::EpotGSSolver( Geometry &geom )
    : EpotSolver( geom ), _epot(NULL), _rhs(NULL), _imax(1000), _eps(1.0e-6), _w(1.66)
{
    
}


EpotGSSolver::EpotGSSolver( Geometry &geom, std::istream &s )
    : EpotSolver(geom)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );    
}


void EpotGSSolver::set_eps( double eps )
{
    _eps = eps;
}

double EpotGSSolver::get_residual( void ) const
{
    return( _res );
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


    return( (1.0/cof) * ( epf - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const
{
    return( (1.0/6.0) * ( (*_epot)(a+1)  + (*_epot)(a-1) +
			  (*_epot)(a+dj) + (*_epot)(a-dj) + 
			  (*_epot)(a+dk) + (*_epot)(a-dk) - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_process_neumann_3d( uint32_t boundary, uint32_t a,
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
 * 2D
 */

double EpotGSSolver::gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
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
    epf += 2.0/(alpha+beta)*( (*_epot)(i,j-1)/alpha + (*_epot)(i,j+1)/beta );

    return( (1.0/cof) * ( epf - (*_rhs)(i,j) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_2d( uint32_t i, uint32_t j ) const
{
    return( (1.0/4.0) * ( (*_epot)(i+1,j) + (*_epot)(i-1,j) +
			  (*_epot)(i,j+1) + (*_epot)(i,j-1) - (*_rhs)(i,j) ) );
}


double EpotGSSolver::gs_process_neumann_2d( uint32_t boundary, uint32_t i, uint32_t j ) const
{
    switch( boundary ) {
    case 1:
	// (phi_i - phi_i+1) / h = q_0
	return( (*_epot)(i+1,j) + (*_rhs)(i,j) );
	break;
    case 2:
	// (phi_i-1 - phi_i) / h = q_0
	return( (*_epot)(i-1,j) + (*_rhs)(i,j) );
	break;
    case 3:
	// (phi_j - phi_j+1) / h = q_0
	return( (*_epot)(i,j+1) + (*_rhs)(i,j) );
	break;
    case 4:
	// (phi_j-1 - phi_j) / h = q_0
	return( (*_epot)(i,j-1) + (*_rhs)(i,j) );
	break;
    }

    return( 0.0 );
}


double EpotGSSolver::gs_loop_2d( void ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    double w2 = 1.0-_w;
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	for( uint32_t i = 0; i < _geom.size(0); i++ ) {
	    
	    double Vold = (*_epot)(i,j);
	    double Vnew;
	    uint32_t mesh = _geom.mesh(i,j);
	    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		Vnew = gs_process_near_solid_2d( nearsolid_ptr, i, j );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		Vnew = gs_process_pure_vacuum_2d( i, j );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		Vnew = gs_process_neumann_2d( boundary, i, j );
	    } else {
		// Dirichlet
		continue;
	    }
	    Vnew = _w*Vnew + w2*Vold;
	    (*_epot)(i,j) = Vnew;
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

    return( (1.0/cof) * ( epf - (*_rhs)(i) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_1d( uint32_t i ) const
{
    return( (1.0/2.0) * ( (*_epot)(i+1) + (*_epot)(i-1) - (*_rhs)(i) ) );
}


double EpotGSSolver::gs_process_neumann_1d( uint32_t boundary, uint32_t i ) const
{
    switch( boundary ) {
    case 1:
	// (phi_i - phi_i+1) / h = q_0
	return( (*_epot)(i+1) + (*_rhs)(i) );
	break;
    case 2:
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
	if( _geom.get_boundary(bound).type == BOUND_NEUMANN ) {
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
	    if( _geom.get_boundary(bound).type == BOUND_NEUMANN ) {
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
	    if( _geom.get_boundary(bound).type == BOUND_NEUMANN ) {
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
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	double z = k*_geom.h()+_geom.origo(2);
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    double y = j*_geom.h()+_geom.origo(1);
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		double x = i*_geom.h()+_geom.origo(0);

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID ||
		    node_id == SMESH_NODE_ID_PURE_VACUUM ) {

		    // Vacuum
		    if( _force_pot_func && (*_force_pot_func)(Vec3D(x,y,z)) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			(*_epot)(i,j,k) = _force_pot;

		    } else if ( _init_plasma_func && (*_init_plasma_func)(Vec3D(x,y,z)) ) {

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
		    (*_rhs)(i,j,k) = _geom.h()*_geom.get_boundary( boundary ).val;

		} else if( node_id == SMESH_NODE_ID_DIRICHLET ) {
		    
		    // Dirichlet
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    (*_epot)(i,j,k) = _geom.get_boundary( boundary ).val;

		}
	    }
	}
    }
    

}


void EpotGSSolver::postprocess( void )
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
	if( _geom.get_boundary(bound).type == BOUND_NEUMANN ) {
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
	    if( _geom.get_boundary(bound).type == BOUND_NEUMANN ) {
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
	    if( _geom.get_boundary(bound).type == BOUND_NEUMANN ) {
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


void EpotGSSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    Timer t;

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
    uint32_t iter = 0;
    while( iter < _imax ) {
	if( _geom.geom_mode() == MODE_3D )
	    _res = gs_loop_3d();
	else if( _geom.geom_mode() == MODE_2D )
	    _res = gs_loop_2d();
	else if( _geom.geom_mode() == MODE_1D )
	    _res = gs_loop_1d();
	else
	    throw( ErrorUnimplemented( ERROR_LOCATION ) );

	if( ibsimu.get_verbose_output() ) {
	    std::stringstream ss;
	    ss << "  " << std::setw(5) << iter << " " << std::scientific << std::setw(20) << _res;
	    sp.print( ss.str() );
	}

	iter++;
	if( _res < _eps )
	    break;
	if( comp_isinf(_res) || comp_isnan(_res) )
	    break;
    }

    // Postprocess
    postprocess();

    // End timer
    t.stop();

    if( ibsimu.get_verbose_output() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << iter << " " << std::scientific << std::setw(20) << _res;
	sp.print( ss.str(), true );
	std::cout << "\n";
	if( iter == _imax )
	    std::cout << "  Maximum number of iteration rounds done.\n";
	std::cout << "  residual error = " << _res << "\n";
	std::cout << "  iterations = " << iter << "\n";
	std::cout << "  time used = " << t << "\n";
	std::cout << std::flush;
    }
}


void EpotGSSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}




void EpotGSSolver::debug_print( std::ostream &os ) const 
{
    os << "**EpotGSSolver\n";
    debug_print_base( os );
    os << "imax = " << _imax << "\n";
    os << "eps = " << _eps << "\n";
    os << "w = " << _w << "\n";

}

