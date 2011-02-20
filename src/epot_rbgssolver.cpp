#include "epot_rbgssolver.hpp"
#include "ibsimu.hpp"
#include "timer.hpp"
#include "constants.hpp"
#include "compmath.hpp"


EpotRBGSSolver::EpotRBGSSolver( const Geometry &geom )
    : EpotSolver( geom ), _rhs(NULL), _imax(10000), _eps(1.0e-6), _w(1.66), _sp(NULL)
{
    pthread_mutex_init( &_mutex, NULL );
    pthread_cond_init( &_cond, NULL );
}



EpotRBGSSolver::EpotRBGSSolver( const Geometry &geom, std::istream &s )
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

    // Debug
    //pthread_mutex_lock( &_mutex );
    //if( rb == 0 )
    //std::cout << "Red (rb=0)\n";
    //else
    //std::cout << "Black (rb=1)\n";
    //pthread_mutex_unlock( &_mutex );

    for( uint32_t k = thno; k < _geom.size(2); k += thcount ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    
	    uint32_t i = (k+j) % 2;
	    if( i == rb )
		i = 0;
	    else
		i = 1;
		    
	    uint32_t a = k*dk+j*dj+i;
	    for( ; i < _geom.size(0); i += 2, a += 2 ) {
		
		// Debug
		//pthread_mutex_lock( &_mutex );
		//std::cout << "Thread " << thno << " processing node (" 
		//<< i << ", " << j << ", " << k << ")\n";
		//pthread_mutex_unlock( &_mutex );

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
	//std::cout << "Thread " << thno << ": done_count = " << _done_count << "\n";
	//std::cout << "Thread " << thno << ": thcount = " << thcount << "\n";

	if( _done_count == thcount ) {
	    // Last thread

	    if( iter == iter_end )
		_done = true;
	    _iter = iter/2;
	    _done_count = 0;

	    // If done both Red and Black
	    if( rb == 0 ) {
		if( ibsimu.get_verbose_output() ) {
		    std::stringstream ss;
		    ss << "  " << std::setw(5) << iter << " " << std::scientific << std::setw(20) << _res;
		    _sp->print( ss.str() );
		}
		if( _res < _eps )
		    _done = true;
		else
		    _res = 0.0;
	    }

	    // Debug
	    //std::cout << "Thread " << thno << " broadcast\n";

	    pthread_cond_broadcast( &_cond );
	} else {

	    // Debug
	    //std::cout << "Thread " << thno << " waiting\n";
	    
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


void EpotRBGSSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    Timer t;

    if( ibsimu.get_verbose_output() )
	std::cout << "  Using Red-Black Gauss-Seidel solver (" 
		  << "w = " << _w
		  << ", imax = " << _imax
		  << ", eps = " << _eps
		  << ", threads = " << ibsimu.get_thread_count()
		  << ")\n";

    // Build rhs and save fields
    _epot = &epot;
    if( _rhs )
	delete _rhs;
    _rhs = new MeshScalarField( (const Mesh)scharge );

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
		    (*_rhs)(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    //std::cout << "neumann\n";
		    (*_rhs)(i,j,k) = _geom.h()*_geom.get_boundary( boundary ).val;
		} else if( node_id == SMESH_NODE_ID_DIRICHLET ) {
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    //std::cout << "dirichlet\n";
		    epot(i,j,k) = _geom.get_boundary( boundary ).val;
		}
	    }
	}
    }

    // Initialize
    if( _sp )
	delete _sp;
    _sp = new StatusPrint( std::cout );
    if( ibsimu.get_verbose_output() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << 0 << " " << std::scientific << std::setw(20) << 0;
	_sp->print( ss.str() );
    }
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

    // End timer
    t.stop();

    if( ibsimu.get_verbose_output() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << _iter << " " << std::scientific << std::setw(20) << _res;
	_sp->print( ss.str(), true );
	std::cout << "\n";
	if( _iter == _imax )
	    std::cout << "  Maximum number of iteration rounds done.\n";
	std::cout << "  residual error = " << _res << "\n";
	std::cout << "  iterations = " << _iter << "\n";
	std::cout << "  time used = " << t << "\n";
	std::cout << std::flush;
    }

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
    os << "**EpotRBGSSolver\n";
    debug_print_base( os );
    os << "imax = " << _imax << "\n";
    os << "eps = " << _eps << "\n";
    os << "w = " << _w << "\n";

}

