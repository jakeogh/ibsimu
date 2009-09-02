#include <cmath>
#include "epot_efield.hpp"


EpotEfield::EpotEfield( const Geometry &g, const ScalarField &epot )
  : _g(g), _epot(epot)
{
    _extrpl[0] = EFIELD_EXTRAPOLATE;
    _extrpl[1] = EFIELD_EXTRAPOLATE;
    _extrpl[2] = EFIELD_EXTRAPOLATE;
    _extrpl[3] = EFIELD_EXTRAPOLATE;
    _extrpl[4] = EFIELD_EXTRAPOLATE;
    _extrpl[5] = EFIELD_EXTRAPOLATE;
}    


/* Return second order polynomial interpolation of 2D field at s.
 * Field is defined by three values p1, p2 and p3 at points s=-1, s=0
 * and s=+1 respectively.
 */
inline double epot_efield_interpolate_1d
( double p1, double p2, double p3, double s )
{
    return( (0.5*(p1+p3)-p2)*s*s + 0.5*(p3-p1)*s + p2 );
}


/* Return second order polynomial interpolation of 2D field at
 * (s,t). Field is defined by nine values p1-p9 at points (-1,+1),
 * (0,+1), (+1,1), (-1,0), (0,0), (+1,0), (-1,-1), (0,-1) and (+1,-1)
 * respectively. */
inline double epot_efield_interpolate_2d
( double p1, double p2, double p3, 
  double p4, double p5, double p6, 
  double p7, double p8, double p9, 
  double s, double t )
{
    double t2;
    double A, B, C, D, E, F, G, H, I;
    
    A = 0.25*(p1-2*p2+p3-2*p4+4*p5-2*p6+p7-2*p8+p9);
    B = 0.25*(-p1+2*p2-p3+p7-2*p8+p9);
    C = 0.5*(p4+p6)-p5;
    D = 0.25*(-p1+p3+2*p4-2*p6-p7+p9);
    E = 0.25*(p1-p3-p7+p9);
    F = 0.5*(p6-p4);
    G = 0.5*(p2+p8)-p5;
    H = 0.5*(p8-p2);
    I = p5;
    t2 = t*t;

    return( (A*t2+B*t+C)*s*s + (D*t2+E*t+F)*s + G*t2+H*t+I );
}


/* Return first derivative of field based on second order polynomial
 * interpolation of field at s. Field is defined by three points p1,
 * p2, p3 at points s=-1, s=0, s=1 respectively. */
inline double epot_efield_interpolate_deriv
( double p1, double p2, double p3, double s )
{
    return( ((p1-2.0*p2+p3)*s + 0.5*(p3-p1)) );
}


Vec3D EpotEfield::operator()( Vec3D x ) const
{
    Vec3D E;
    double inv_h = 1.0/_g.h();

    switch( _g.geom_mode() ) {
    case MODE_1D:
    {
	double sign = -1.0;
	if( x[0] < _g.origo(0) ) {
	    if( _extrpl[0] == EFIELD_MIRROR ) {
		sign *= -1.0;
		x[0]  = 2.0*_g.origo(0) - x[0];
	    } else if( _extrpl[0] == EFIELD_ZERO ) {
		// return zero
		return( E );
	    }
	} else if( x[0] > _g.max(0) ) {
	    if( _extrpl[1] == EFIELD_MIRROR ) {
		sign *= -1.0;
		x[0]  = 2.0*_g.max(0) - x[0];
	    } else if( _extrpl[1] == EFIELD_ZERO ) {
		// return zero
		return( E );
	    }
	}

	int32_t i = (int32_t)floor( (x[0]-_g.origo(0))*inv_h + 0.5 );
	if( i < 1 )
	    i = 1;
	else if( i >= _g.size(0)-1 )
	    i = _g.size(0)-2;
	
	if( _g.mesh( i ) > 0 )
	    // Inside solid, return zero
	    return( E );
	else if( _g.mesh( i ) < 0 ) {
	    // At solid edge, search free direction
	    if( _g.mesh(i+1) == 0 ) i++;
	    else if( _g.mesh(i-1) == 0 ) i--;
	}
	
	double t = ( x[0]-(i*_g.h()+_g.origo(0)) )*inv_h;
	E[0] = sign*epot_efield_interpolate_deriv( 
	    _epot( i-1 ), _epot( i ), _epot( i+1 ), t )*inv_h;
	break;
    }
    case MODE_2D:
    case MODE_CYL:
    {
	double sign[2] = {-1.0, -1.0};
	for( int a = 0; a < 2; a++ ) {
	    if( x[a] < _g.origo(a) ) {
		if( _extrpl[2*a] == EFIELD_MIRROR ) {
		    sign[a] *= -1.0;
		    x[a]     = 2.0*_g.origo(a) - x[a];
		} else if( _extrpl[2*a] == EFIELD_ZERO ) {
		    // return zero
		    return( E );
		}
	    } else if( x[0] > _g.max(0) ) {
		if( _extrpl[2*a+1] == EFIELD_MIRROR ) {
		    sign[a] *= -1.0;
		    x[a]     = 2.0*_g.max(a) - x[a];
		} else if( _extrpl[2*a+1] == EFIELD_ZERO ) {
		    // return zero
		    return( E );
		}
	    }
	}
	
	int32_t i = (int32_t)floor( (x[0]-_g.origo(0))*inv_h + 0.5 );
	int32_t j = (int32_t)floor( (x[1]-_g.origo(1))*inv_h + 0.5 );

	if( i < 1 )
	    i = 1;
	else if( i >= _g.size(0)-1 )
	    i = _g.size(0)-2;
	if( j < 1 )
	    j = 1;
	else if( j >= _g.size(1)-1 )
	    j = _g.size(1)-2;

	if( _g.mesh(i,j) > 0 )
	    // Inside solid, return zero
	    return( E );
	else if( _g.mesh(i,j) < 0 ) {
	    // At solid edge, search free direction
	    if( _g.mesh(i-1,j) == 0 && _g.mesh(i+1,j) != 0 && _g.mesh(i,j+1) != 0 && _g.mesh(i,j-1) != 0 )
		i--;
	    else if( _g.mesh(i+1,j) == 0 && _g.mesh(i-1,j) != 0 && _g.mesh(i,j+1) != 0 && _g.mesh(i,j-1) != 0 )
		i++;
	    else if( _g.mesh(i,j-1) == 0 && _g.mesh(i,j+1) != 0 && _g.mesh(i+1,j) != 0 && _g.mesh(i-1,j) != 0 )
		j--;
	    else if( _g.mesh(i,j+1) == 0 && _g.mesh(i,j-1) != 0 && _g.mesh(i+1,j) != 0 && _g.mesh(i-1,j) != 0 )
		j++;
	    else if( _g.mesh(i-1,j-1) == 0 && _g.mesh(i-1,j) == 0 && _g.mesh(i,j-1) == 0 && 
		     _g.mesh(i+1,j) != 0 && _g.mesh(i,j+1) != 0 ) {
		i--; j--;
	    } else if( _g.mesh(i+1,j-1) == 0 && _g.mesh(i+1,j) == 0 && _g.mesh(i,j-1) == 0 && 
		       _g.mesh(i-1,j) != 0 && _g.mesh(i,j+1) != 0 ) {
		i++; j--;
	    } else if( _g.mesh(i-1,j+1) == 0 && _g.mesh(i-1,j) == 0 && _g.mesh(i,j+1) == 0 && 
		       _g.mesh(i+1,j) != 0 && _g.mesh(i,j-1) != 0 ) {
		i--; j++;
	    } else if( _g.mesh(i+1,j+1) == 0 && _g.mesh(i+1,j) == 0 && _g.mesh(i,j+1) == 0 && 
		       _g.mesh(i-1,j) != 0 && _g.mesh(i,j-1) != 0 ) {
		i++; j++;
	    }
	}
	
	double t = ( x[0]-(i*_g.h()+_g.origo(0)) )*inv_h;
	double u = ( x[1]-(j*_g.h()+_g.origo(1)) )*inv_h;

	E[0] = sign[0]*epot_efield_interpolate_deriv( 
	    epot_efield_interpolate_1d( _epot( i-1, j-1 ), 
					_epot( i-1, j   ), 
					_epot( i-1, j+1 ), u ),
	    epot_efield_interpolate_1d( _epot( i,   j-1 ), 
					_epot( i,   j   ), 
					_epot( i,   j+1 ), u ),
	    epot_efield_interpolate_1d( _epot( i+1, j-1 ), 
					_epot( i+1, j   ), 
					_epot( i+1, j+1 ), u ), t )*inv_h;
	E[1] = sign[1]*epot_efield_interpolate_deriv( 
	    epot_efield_interpolate_1d( _epot( i-1, j-1 ), 
					_epot( i,   j-1 ), 
					_epot( i+1, j-1 ), t ),
	    epot_efield_interpolate_1d( _epot( i-1, j   ), 
					_epot( i,   j   ), 
					_epot( i+1, j   ), t ),
	    epot_efield_interpolate_1d( _epot( i-1, j+1 ), 
					_epot( i,   j+1 ), 
					_epot( i+1, j+1 ), t ), u )*inv_h;
	break;
    }
    case MODE_3D:
    {
	double sign[3] = {-1.0, -1.0, -1.0};
	for( int a = 0; a < 3; a++ ) {
	    if( x[a] < _g.origo(a) ) {
		if( _extrpl[2*a] == EFIELD_MIRROR ) {
		    sign[a] *= -1.0;
		    x[a]     = 2.0*_g.origo(a) - x[a];
		} else if( _extrpl[2*a] == EFIELD_ZERO ) {
		    // return zero
		    return( E );
		}
	    } else if( x[0] > _g.max(0) ) {
		if( _extrpl[2*a+1] == EFIELD_MIRROR ) {
		    sign[a] *= -1.0;
		    x[a]     = 2.0*_g.max(a) - x[a];
		} else if( _extrpl[2*a+1] == EFIELD_ZERO ) {
		    // return zero
		    return( E );
		}
	    }
	}
	
	int32_t i = (int32_t)floor( (x[0]-_g.origo(0))*inv_h + 0.5 );
	int32_t j = (int32_t)floor( (x[1]-_g.origo(1))*inv_h + 0.5 );
	int32_t k = (int32_t)floor( (x[2]-_g.origo(2))*inv_h + 0.5 );

	if( i < 1 )
	    i = 1;
	else if( i >= _g.size(0)-1 )
	    i = _g.size(0)-2;
	if( j < 1 )
	    j = 1;
	else if( j >= _g.size(1)-1 )
	    j = _g.size(1)-2;
	if( k < 1 )
	    k = 1;
	else if( k >= _g.size(2)-1 )
	    k = _g.size(2)-2;

	if( _g.mesh(i,j,k) > 0 )
	    // Inside solid, return zero
	    return( E );
	else if( _g.mesh(i,j,k) < 0 ) {
	    // At solid edge, search free direction
	    // On face
	    if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		_g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) != 0 &&
		_g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) != 0 )
		i--;
	    else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		     _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) != 0 &&
		     _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) != 0 )
		i++;
	    else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) != 0 && 
		     _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 &&
		     _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) != 0 )
		j--;
	    else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) != 0 && 
		     _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 &&
		     _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) != 0 )
		j++;
	    else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) != 0 && 
		     _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) != 0 &&
		     _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 )
		k--;
	    else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) != 0 && 
		     _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) != 0 &&
		     _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 )
		k++;
	    // On edge
	    else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		     _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 && 
		     _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) != 0 &&
		     _g.mesh(i-1,j-1,k) == 0 ) {
		i--; j--;
	    } else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i+1,j-1,k) == 0 ) {
		i++; j--;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		       _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i-1,j+1,k) == 0 ) {
		i--; j++;
	    } else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i+1,j+1,k) == 0 ) {
		i++; j++;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 &&
		       _g.mesh(i-1,j,k-1) == 0 ) {
		i--; k--;
	    } else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 &&
		       _g.mesh(i+1,j,k-1) == 0 ) {
		i++; k--;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i-1,j,k+1) == 0 ) {
		i--; k++;
	    } else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i+1,j,k+1) == 0 ) {
		i++; k++;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 &&
		       _g.mesh(i,j-1,k-1) == 0 ) {
		j--; k--;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 &&
		       _g.mesh(i,j+1,k-1) == 0 ) {
		j++; k--;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 && 
		       _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i,j-1,k+1) == 0 ) {
		j--; k++;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i,j+1,k+1) == 0 ) {
		j++; k++;
	    }
	    // On corner
	    else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		     _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 && 
		     _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 &&
		     _g.mesh(i-1,j-1,k-1) == 0 ) {
		i--; j--; k--;
	    } else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 &&
		       _g.mesh(i+1,j-1,k-1) == 0 ) {
		i++; j--; k--;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		       _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 &&
		       _g.mesh(i-1,j+1,k-1) == 0 ) {
		i--; j++; k--;
	    } else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) != 0 && _g.mesh(i,j,k-1) == 0 &&
		       _g.mesh(i+1,j+1,k-1) == 0 ) {
		i++; j++; k--;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 && 
		       _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i-1,j-1,k+1) == 0 ) {
		i--; j--; k++;
	    } else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) != 0 && _g.mesh(i,j-1,k) == 0 && 
		       _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i+1,j-1,k+1) == 0 ) {
		i++; j--; k++;
	    } else if( _g.mesh(i+1,j,k) != 0 && _g.mesh(i-1,j,k) == 0 && 
		       _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i-1,j+1,k+1) == 0 ) {
		i--; j++; k++;
	    } else if( _g.mesh(i+1,j,k) == 0 && _g.mesh(i-1,j,k) != 0 && 
		       _g.mesh(i,j+1,k) == 0 && _g.mesh(i,j-1,k) != 0 && 
		       _g.mesh(i,j,k+1) == 0 && _g.mesh(i,j,k-1) != 0 &&
		       _g.mesh(i+1,j+1,k+1) == 0 ) {
		i++; j++; k++;
	    } 
	}
	
	double t = ( x[0]-(i*_g.h()+_g.origo(0)) )*inv_h;
	double u = ( x[1]-(j*_g.h()+_g.origo(1)) )*inv_h;
	double v = ( x[2]-(k*_g.h()+_g.origo(2)) )*inv_h;

	E[0] = sign[0]*epot_efield_interpolate_deriv( 
	    epot_efield_interpolate_2d( _epot( i-1, j-1, k-1 ), 
					_epot( i-1, j,   k-1 ), 
					_epot( i-1, j+1, k-1 ),
					_epot( i-1, j-1, k   ),
					_epot( i-1, j,   k   ),
					_epot( i-1, j+1, k   ),
					_epot( i-1, j-1, k+1 ),
					_epot( i-1, j,   k+1 ), 
					_epot( i-1, j+1, k+1 ), u, v ),
	    epot_efield_interpolate_2d( _epot( i,   j-1, k-1 ), 
					_epot( i,   j,   k-1 ), 
					_epot( i,   j+1, k-1 ),
					_epot( i,   j-1, k   ),
					_epot( i,   j,   k   ),
					_epot( i,   j+1, k   ),
					_epot( i,   j-1, k+1 ),
					_epot( i,   j,   k+1 ), 
					_epot( i,   j+1, k+1 ), u, v ),
	    epot_efield_interpolate_2d( _epot( i+1, j-1, k-1 ), 
					_epot( i+1, j,   k-1 ), 
					_epot( i+1, j+1, k-1 ),
					_epot( i+1, j-1, k   ),
					_epot( i+1, j,   k   ),
					_epot( i+1, j+1, k   ),
					_epot( i+1, j-1, k+1 ),
					_epot( i+1, j,   k+1 ), 
					_epot( i+1, j+1, k+1 ), u, v ), t )*inv_h;

	E[1] = sign[1]*epot_efield_interpolate_deriv( 
	    epot_efield_interpolate_2d( _epot( i-1, j-1, k-1 ), 
					_epot( i,   j-1, k-1 ), 
					_epot( i+1, j-1, k-1 ),
					_epot( i-1, j-1, k   ), 
					_epot( i,   j-1, k   ), 
					_epot( i+1, j-1, k   ),
					_epot( i-1, j-1, k+1 ), 
					_epot( i,   j-1, k+1 ), 
					_epot( i+1, j-1, k+1 ), t, v ),
	    epot_efield_interpolate_2d( _epot( i-1, j,   k-1 ), 
					_epot( i,   j,   k-1 ), 
					_epot( i+1, j,   k-1 ),
					_epot( i-1, j,   k   ), 
					_epot( i,   j,   k   ), 
					_epot( i+1, j,   k   ),
					_epot( i-1, j,   k+1 ), 
					_epot( i,   j,   k+1 ), 
					_epot( i+1, j,   k+1 ), t, v ),
	    epot_efield_interpolate_2d( _epot( i-1, j+1, k-1 ), 
					_epot( i,   j+1, k-1 ), 
					_epot( i+1, j+1, k-1 ),
					_epot( i-1, j+1, k   ), 
					_epot( i,   j+1, k   ), 
					_epot( i+1, j+1, k   ),
					_epot( i-1, j+1, k+1 ), 
					_epot( i,   j+1, k+1 ), 
					_epot( i+1, j+1, k+1 ), t, v ), u )*inv_h;
	
	E[2] = sign[2]*epot_efield_interpolate_deriv( 
	    epot_efield_interpolate_2d( _epot( i-1, j-1, k-1 ), 
					_epot( i,   j-1, k-1 ), 
					_epot( i+1, j-1, k-1 ),
					_epot( i-1, j,   k-1 ), 
					_epot( i,   j,   k-1 ), 
					_epot( i+1, j,   k-1 ),
					_epot( i-1, j+1, k-1 ), 
					_epot( i,   j+1, k-1 ), 
					_epot( i+1, j+1, k-1 ), t, u ),
	    epot_efield_interpolate_2d( _epot( i-1, j-1, k   ), 
					_epot( i,   j-1, k   ), 
					_epot( i+1, j-1, k   ),
					_epot( i-1, j,   k   ), 
					_epot( i,   j,   k   ), 
					_epot( i+1, j,   k   ),
					_epot( i-1, j+1, k   ), 
					_epot( i,   j+1, k   ), 
					_epot( i+1, j+1, k   ), t, u ),
	    epot_efield_interpolate_2d( _epot( i-1, j-1, k+1 ), 
					_epot( i,   j-1, k+1 ), 
					_epot( i+1, j-1, k+1 ),
					_epot( i-1, j,   k+1 ), 
					_epot( i,   j,   k+1 ), 
					_epot( i+1, j,   k+1 ),
					_epot( i-1, j+1, k+1 ), 
					_epot( i,   j+1, k+1 ), 
					_epot( i+1, j+1, k+1 ), t, u ), v )*inv_h;
	break;
    }
    default:
    {
	throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	break;
    }
    }
    
    return( E );
}














