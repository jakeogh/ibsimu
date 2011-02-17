#include "epot_gssolver.hpp"
#include "ibsimu.hpp"
#include "timer.hpp"
#include "constants.hpp"
#include "compmath.hpp"
#include "statusprint.hpp"


EpotGSSolver::EpotGSSolver( const Geometry &geom )
    : EpotSolver( geom ), _imax(1000), _eps(1.0e-6), _w(1.66)
{
    
}



EpotGSSolver::EpotGSSolver( const Geometry &geom, std::istream &s )
    : EpotSolver(geom)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );    
}


void EpotGSSolver::set_eps( double eps )
{
    _eps = eps;
}

double EpotGSSolver::get_err( void ) const
{
    return( _err );
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

double EpotGSSolver::gs_process_near_solid_3d( MeshScalarField &epot, const MeshScalarField &rhs,
					       const uint8_t *nearsolid_ptr, uint32_t a,
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
    double epf = 2.0/(alpha+beta)*( epot(a-1)/alpha + epot(a+1)/beta );

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
    epf += 2.0/(alpha+beta)*( epot(a-dj)/alpha + epot(a+dj)/beta );

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
    epf += 2.0/(alpha+beta)*( epot(a-dk)/alpha + epot(a+dk)/beta );


    return( (1.0/cof) * ( epf - rhs(a) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_3d( MeshScalarField &epot, const MeshScalarField &rhs,
						uint32_t a, uint32_t dj, uint32_t dk ) const
{
    return( (1.0/6.0) * ( epot(a+1)  + epot(a-1) +
			  epot(a+dj) + epot(a-dj) + 
			  epot(a+dk) + epot(a-dk) - rhs(a) ) );
}


double EpotGSSolver::gs_process_neumann_3d( MeshScalarField &epot, const MeshScalarField &rhs,
					    uint32_t boundary, uint32_t a,
					    uint32_t dj, uint32_t dk ) const
{
    switch( boundary ) {
    case 1:
	// (phi_i - phi_i+1) / h = q_0
	return( epot(a+1) + rhs(a) );
	break;
    case 2:
	// (phi_i-1 - phi_i) / h = q_0
	return( epot(a-1) + rhs(a) );
	break;
    case 3:
	// (phi_j - phi_j+1) / h = q_0
	return( epot(a+dj) + rhs(a) );
	break;
    case 4:
	// (phi_j-1 - phi_j) / h = q_0
	return( epot(a-dj) + rhs(a) );
	break;
    case 5:
	// (phi_k - phi_k+1) / h = q_0
	return( epot(a+dk) + rhs(a) );
	break;
    case 6:
	// (phi_k-1 - phi_k) / h = q_0
	return( epot(a-dk) + rhs(a) );
	break;
    }

    return( 0.0 );
}


double EpotGSSolver::gs_loop_3d( MeshScalarField &epot, const MeshScalarField &rhs ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    const double w2 = 1.0-_w;
    const uint32_t dj = _geom.size(0);
    const uint32_t dk = _geom.size(0)*_geom.size(1);

    //for( uint32_t rb = 0; rb < 2; rb++ ) {
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	//double z = k*_geom.h()+_geom.origo(2);
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    uint32_t a = k*dk+j*dj;
	    for( uint32_t i = 0; i < _geom.size(0); i++, a++ ) {
		//double x = i*_geom.h()+_geom.origo(0);

		//uint32_t s = i+j+k;
		//if( s % 2 == rb )
		//continue;

		double Vold = epot(a);
		double Vnew;
		uint32_t mesh = _geom.mesh(a);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		    Vnew = gs_process_near_solid_3d( epot, rhs, nearsolid_ptr, a, dj, dk );
		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		    Vnew = gs_process_pure_vacuum_3d( epot, rhs, a, dj, dk );
		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    Vnew = gs_process_neumann_3d( epot, rhs, boundary, a, dj, dk );
		} else {
		    // Dirichlet
		    continue;
		}
		Vnew = _w*Vnew + w2*Vold;
		epot(a) = Vnew;
		double err = fabs( Vnew - Vold );
		if( err > maxerr )
		    maxerr = err;
		if( comp_isinf(err) || comp_isnan(err) ) {
		    std::cout << "location = " << i << ", " << j << ", " << k << "\n";
		    return( err );
		}
	    }
	}
    }
    //}
    return( maxerr );
}


/* *****************************************************************************
 * 2D
 */

double EpotGSSolver::gs_process_near_solid_2d( MeshScalarField &epot, const MeshScalarField &rhs,
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


double EpotGSSolver::gs_process_pure_vacuum_2d( MeshScalarField &epot, const MeshScalarField &rhs,
						uint32_t i, uint32_t j ) const
{
    return( (1.0/4.0) * ( epot(i+1,j) + epot(i-1,j) +
			  epot(i,j+1) + epot(i,j-1) - rhs(i,j) ) );
}


double EpotGSSolver::gs_process_neumann_2d( MeshScalarField &epot, const MeshScalarField &rhs,
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


double EpotGSSolver::gs_loop_2d( MeshScalarField &epot, const MeshScalarField &rhs ) const
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


/* *****************************************************************************
 * 1D
 */

double EpotGSSolver::gs_process_near_solid_1d( MeshScalarField &epot, const MeshScalarField &rhs,
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


double EpotGSSolver::gs_process_pure_vacuum_1d( MeshScalarField &epot, const MeshScalarField &rhs,
						uint32_t i ) const
{
    return( (1.0/2.0) * ( epot(i+1) + epot(i-1) - rhs(i) ) );
}


double EpotGSSolver::gs_process_neumann_1d( MeshScalarField &epot, const MeshScalarField &rhs,
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


double EpotGSSolver::gs_loop_1d( MeshScalarField &epot, const MeshScalarField &rhs ) const
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


/* *****************************************************************************
 * Common
 */


void EpotGSSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge ) const
{
    Timer t;

    if( ibsimu.get_verbose_output() )
	std::cout << "  Using Gauss-Seidel solver (" 
		  << "w = " << _w
		  << ", imax = " << _imax
		  << ", eps = " << _eps
		  << ")\n";

    // Build rhs
    MeshScalarField *rhs = new MeshScalarField( (const Mesh)scharge );

    // Build rhs and set forced and dirichlet nodes to epot
    //std::cout << "\nPreprocess\n------------------------\n";
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	//double z = k*_geom.h()+_geom.origo(2);
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    //double y = j*_geom.h()+_geom.origo(1);
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		//double x = i*_geom.h()+_geom.origo(0);

		//std::cout << "(" 
		//	  << i << ", "
		//	  << j << ", "
		//	  << k << "): ";
		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID ||
		    node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		    //std::cout << "vacuum\n";
		    (*rhs)(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    //std::cout << "neumann\n";
		    (*rhs)(i,j,k) = _geom.h()*_geom.get_boundary( boundary ).val;
		} else if( node_id == SMESH_NODE_ID_DIRICHLET ) {
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    //std::cout << "dirichlet\n";
		    epot(i,j,k) = _geom.get_boundary( boundary ).val;
		}
	    }
	}
    }
    
    StatusPrint sp( std::cout );
    if( ibsimu.get_verbose_output() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << 0 << " " << std::scientific << std::setw(20) << 0;
	sp.print( ss.str() );
    }

    // Loop until converged
    //std::cout << "\nIterate\n------------------------\n";
    uint32_t iter = 0;
    while( iter < _imax ) {
	if( _geom.geom_mode() == MODE_3D )
	    _err = gs_loop_3d( epot, *rhs );
	else if( _geom.geom_mode() == MODE_2D )
	    _err = gs_loop_2d( epot, *rhs );
	else if( _geom.geom_mode() == MODE_1D )
	    _err = gs_loop_1d( epot, *rhs );
	else
	    throw( ErrorUnimplemented( ERROR_LOCATION ) );

	//if( ibsimu.get_verbose_output() ) 
	//std::cout << "  " << std::setw(5) << iter << " " << err << "\n";
	
	if( ibsimu.get_verbose_output() ) {
	    std::stringstream ss;
	    ss << "  " << std::setw(5) << iter << " " << std::scientific << std::setw(20) << _err;
	    sp.print( ss.str() );
	}

	iter++;
	if( _err < _eps )
	    break;
	if( comp_isinf(_err) || comp_isnan(_err) )
	    break;
    }

    // Free rhs
    delete rhs;

    // End timer
    t.stop();

    if( ibsimu.get_verbose_output() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << iter << " " << std::scientific << std::setw(20) << _err;
	sp.print( ss.str(), true );
	std::cout << "\n";
	if( iter == _imax )
	    std::cout << "  Maximum number of iteration rounds done.\n";
	std::cout << "  residual error = " << _err << "\n";
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
    debug_print_base( os );
    os << "imax = " << _imax << "\n";
    os << "eps = " << _eps << "\n";
    os << "w = " << _w << "\n";

}






