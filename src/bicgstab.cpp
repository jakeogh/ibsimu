#include <iostream>
#include <iomanip>
#include <sstream>

#include "bicgstab.hpp"
#include "verbose.hpp"
#include "statusprint.hpp"


bool bicgstab( const Matrix &mat, const Vector &rhs, Vector &sol,
	       const Precond &pc, uint32_t &imax, double &eps )
{
    // Checks
    if( mat.columns() != mat.rows() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );
    if( mat.rows() != rhs.size() )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    double resid, omega = 0, alpha = 0, beta, rho_1, rho_2 = 0;
    Vector p, phat, s, shat, t, v;
    double norm_rhs = norm2(rhs);
    if( sol.size() != mat.columns() ) {
	sol.resize( mat.columns() );
	sol.clear();
    }

    Vector r = mat * sol;
    r = rhs - r;
    Vector rtilde = r;

    if( norm_rhs == 0.0 )
	norm_rhs = 1;
  
    if( (resid = norm2(r) / norm_rhs) <= eps ) {
	eps = resid;
	imax = 0;
	return( true );
    }

    StatusPrint sp;
    if( verbose_output ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << 0 << " " << std::setw(20) << resid;
	sp.print( ss.str() );
    }

    uint32_t i;
    bool retval = false;
    for( i = 1; i <= imax; i++ ) {
	rho_1 = dot_prod( rtilde, r );
	if( rho_1 == 0 ) {
	    eps = norm2(r) / norm_rhs;
	    imax = i;
	    break;
	}
	if( i == 1 )
	    p = r;
	else {
	    beta = (rho_1/rho_2) * (alpha/omega);
	    p = r + beta * (p - omega * v);
	}
	pc.solve( phat, p );
	v = mat * phat;
	alpha = rho_1 / dot_prod( rtilde, v );
	s = r - alpha * v;
	if( (resid = norm2(s)/norm_rhs) < eps ) {
	    sol += alpha * phat;
	    eps = resid;
	    imax = i;
	    retval = true;
	    break;
	}
	pc.solve( shat, s );
	t = mat * shat;
	omega = dot_prod( t, s ) / dot_prod( t, t );
	sol += alpha * phat + omega * shat;
	r = s - omega * t;

	rho_2 = rho_1;
	if( (resid = norm2(r) / norm_rhs) < eps ) {
	    eps = resid;
	    imax = i;
	    retval = true;
	    break;
	}
	if( omega == 0 ) {
	    eps = norm2(r) / norm_rhs;
	    imax = i;
	    retval = false;
	    break;
	}

	if( verbose_output ) {
	    std::stringstream ss;
	    ss << "  " << std::setw(5) << i << " " << std::setw(20) << resid;
	    sp.print( ss.str() );
	}
    }

    if( i > imax ) {
	eps = resid;
	imax = i;
    }
    
    return( retval );
}
















