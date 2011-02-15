#include "epot_gssolver.hpp"
#include "constants.hpp"


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


void EpotGSSolver::set_imax( uint32_t imax )
{
    _imax = imax;
}

double EpotGSSolver::gs_process_near_solid( MeshScalarField &epot, const MeshScalarField &rhs,
					    const uint8_t *nearsolid_ptr, 
					    uint32_t i, uint32_t j, uint32_t k ) const
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
    double epf = 2.0/(alpha+beta)*( epot(i-1,j,k)/alpha + epot(i+1,j,k)/beta );

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
    epf += 2.0/(alpha+beta)*( epot(i,j-1,k)/alpha + epot(i,j+1,k)/beta );

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
    epf += 2.0/(alpha+beta)*( epot(i,j-1,k)/alpha + epot(i,j+1,k)/beta );


    return( (1.0/cof) * ( epf - rhs(i,j,k) ) );
}


double EpotGSSolver::gs_process_pure_vacuum( MeshScalarField &epot, const MeshScalarField &rhs,
					     uint32_t i, uint32_t j, uint32_t k ) const
{
    return( (1.0/6.0) * ( epot(i+1,j,k) + epot(i-1,j,k) +
			  epot(i,j+1,k) + epot(i,j-1,k) + 
			  epot(i,j,k+1) + epot(i,j,k-1) - rhs(i,j,k) ) );
}


double EpotGSSolver::gs_process_neumann( MeshScalarField &epot, const MeshScalarField &rhs,
					 uint32_t boundary, uint32_t i, uint32_t j, uint32_t k ) const
{
    switch( boundary ) {
    case 1:
	// (phi_i - phi_i+1) / h = q_0
	return( epot(i+1,j,k) + rhs(i,j,k) );
	break;
    case 2:
	// (phi_i-1 - phi_i) / h = q_0
	return( epot(i-1,j,k) + rhs(i,j,k) );
	break;
    case 3:
	// (phi_j - phi_j+1) / h = q_0
	return( epot(i,j+1,k) + rhs(i,j,k) );
	break;
    case 4:
	// (phi_j-1 - phi_j) / h = q_0
	return( epot(i,j-1,k) + rhs(i,j,k) );
	break;
    case 5:
	// (phi_k - phi_k+1) / h = q_0
	return( epot(i,j,k+1) + rhs(i,j,k) );
	break;
    case 6:
	// (phi_k-1 - phi_k) / h = q_0
	return( epot(i,j,k-1) + rhs(i,j,k) );
	break;
    }

    return( 0.0 );
}


double EpotGSSolver::gs_process_node( MeshScalarField &epot, const MeshScalarField &rhs,
				      uint32_t i, uint32_t j, uint32_t k ) const
{
    uint32_t mesh = _geom.mesh(i,j,k);
    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
	const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
	return( gs_process_near_solid( epot, rhs, nearsolid_ptr, i, j, k ) );
    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
	return( gs_process_pure_vacuum( epot, rhs, i, j, k ) );
    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
	uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
	return( gs_process_neumann( epot, rhs, boundary, i, j, k ) );
    }

    return( 0.0 );
}			   


double EpotGSSolver::gs_loop( MeshScalarField &epot, const MeshScalarField &rhs ) const
{
    // Go through all nodes
    double maxerr = 0.0;
    double w2 = 1.0-_w;
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	//double z = k*_geom.h()+_geom.origo(2);
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    //double y = j*_geom.h()+_geom.origo(1);
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		//double x = i*_geom.h()+_geom.origo(0);

		double Vold = epot(i,j,k);
		double Vnew = gs_process_node( epot, rhs, i, j, k );
		Vnew = _w*Vnew + w2*Vold;
		epot(i,j,k) = Vnew;
		double err = fabs( Vnew - Vold );
		if( err > maxerr )
		    maxerr = err;
	    }
	}
    }
    return( maxerr );
}


void EpotGSSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge ) const
{
    // Build rhs
    MeshScalarField *rhs = new MeshScalarField( scharge );
    (*rhs) *= -_geom.h()*_geom.h()/EPSILON0;
    
    // Loop until converged
    for( uint32_t iter = 0; iter < _imax; iter++ ) {
	double err = gs_loop( epot, *rhs );
	if( err < _eps )
	    break;
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






