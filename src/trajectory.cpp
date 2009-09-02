#include <cmath>
#include <iostream>
#include "trajectory.hpp"
#include "error.hpp"
#include "polysolver.hpp"


void TrajectoryRep1D::construct( double dt, double x1, double v1, double x2, double v2,
				 trajectory_rep_e force )
{
    double sca = 1.0/(x2-x1);

    switch( force ) {
    case TRAJ_EMPTY:
	_A = (v2+v1)*dt + 2.0*(x1-x2);
	if( fabs(_A*sca) < 1.0e-6 ) {
	    _A = x2-x1-v1*dt;
	    if( fabs(_A*sca) < 1.0e-6 ) {
		_A = x2-x1;
		_B = x1;
		_C = v2-v1;
		_D = v1;
		_rep = TRAJ_LINEAR;
		return;
	    }
	    _B = v1*dt;
	    _C = x1;
	    _D = v2-v1;
	    _E = v1;
	    _rep = TRAJ_QUADRATIC;
	    return;
	}
	_B = 3.0*(x2-x1) - (2.0*v1+v2)*dt;
	_C = v1*dt;
	_D = x1;
	_E = dt;
	_rep = TRAJ_CUBIC;
	break;
    case TRAJ_LINEAR:
	_A = x2-x1;
	_B = x1;
	_C = v2-v1;
	_D = v1;
	_rep = TRAJ_LINEAR;
	break;
    case TRAJ_QUADRATIC:
	_A = x2-x1-v1*dt;
	_B = v1*dt;
	_C = x1;
	_D = v2-v1;
	_E = v1;
	_rep = TRAJ_QUADRATIC;
	break;
    case TRAJ_CUBIC:
	_A = (v2+v1)*dt + 2.0*(x1-x2);
	_B = 3.0*(x2-x1) - (2.0*v1+v2)*dt;
	_C = v1*dt;
	_D = x1;
	_E = dt;
	_rep = TRAJ_CUBIC;
	break;
    };
}


TrajectoryRep1D::TrajectoryRep1D( double dt, double x1, double v1, double x2, double v2,
				  trajectory_rep_e force )
{
    construct( dt, x1, v1, x2, v2, force );
}


void TrajectoryRep1D::coord( double &x, double &v, double K )
{
    switch( _rep ) {
    case TRAJ_EMPTY:
	throw( Error( ERROR_LOCATION, "empty representation" ) );
	break;
    case TRAJ_LINEAR:
	x = _A*K + _B;
	v = _C*K + _D;
	break;
    case TRAJ_QUADRATIC:
	x = (_A*K + _B)*K + _C;
	v = _D*K + _E;
	break;
    case TRAJ_CUBIC:
	x = ((_A*K + _B)*K + _C)*K + _D;
	v = ((3.0*_A*K + 2.0*_B)*K + _C)/_E;
	break;
    };
}


int TrajectoryRep1D::solve( double K[3], double x )
{
    switch( _rep ) {
    case TRAJ_EMPTY:
	throw( Error( ERROR_LOCATION, "empty representation" ) );
	break;
    case TRAJ_LINEAR:
	if( _A == 0.0 )
	    return( 0 );
	K[0] = (x-_B) / _A;
	if( K[0] > 0.0 && K[0] <= 1.0 )
	    return( 1 );
	break;
    case TRAJ_QUADRATIC:
    {
	int nroots = solve_quadratic( _A, _B, _C-x, &K[0], &K[1] );
	if( nroots == 0 ) {
	    return( 0 );
	} else if( nroots == 1 ) {
	    if( K[0] > 0.0 && K[0] <= 1.0 ) {
		return( 1 );
	    } else {
		return( 0 );
	    }
	} else /* nroots = 2 */ {
	    if( K[0] > 0.0 && K[0] <= 1.0 ) {
		if( K[1] > 0.0 && K[1] <= 1.0 ) {
		    return( 2 );
		} else {
		    return( 1 );
		}
	    } else {
		if( K[1] > 0.0 && K[1] <= 1.0 ) {
		    K[0] = K[1];
		    return( 1 );
		} else {
		    return( 0 );
		}
	    }
	}
	break;
    }
    case TRAJ_CUBIC:
    {
	int nroots = solve_cubic( _A, _B, _C, _D-x, &K[0], &K[1], &K[2] );
	if( nroots == 0 ) {
	    return( 0 );
	} else if( nroots == 1 ) {
	    if( K[0] > 0.0 && K[0] <= 1.0 ) {
		return( 1 );
	    }
	} else if( nroots == 2 ) {
	    if( K[0] > 0.0 && K[0] <= 1.0 ) {
		if( K[1] > 0.0 && K[1] <= 1.0 ) {
		    return( 2 );
		} else {
		    return( 1 );
		}
	    } else {
		if( K[1] > 0.0 && K[1] <= 1.0 ) {
		    K[0] = K[1];
		    return( 1 );
		} else {
		    return( 0 );
		}
	    }
	} else /* nroots = 3 */ {
	    if( K[0] > 0.0 && K[0] <= 1.0 ) {
		if( K[1] > 0.0 && K[1] <= 1.0 ) {
		    if( K[2] > 0.0 && K[2] <= 1.0 ) {
			return( 3 );
		    } else {
			return( 2 );
		    }
		} else {
		    if( K[2] > 0.0 && K[2] <= 1.0 ) {
			K[1] = K[2];
			return( 2 );
		    } else {
			return( 1 );
		    }
		}
	    } else {
		if( K[1] > 0.0 && K[1] <= 1.0 ) {
		    if( K[2] > 0.0 && K[2] <= 1.0 ) {
			K[0] = K[1];
			K[1] = K[2];
			return( 2 );
		    } else {
			K[0] = K[1];
			return( 1 );
		    }
		} else {
		    if( K[2] > 0.0 && K[2] <= 1.0 ) {
			K[0] = K[2];
			return( 1 );
		    } else {
			return( 0 );
		    }
		}
	    }
	}
	break;
    }
    };

    return( 0 );
}


void TrajectoryRep1D::debug_print( void ) const
{
    switch( _rep ) {
    case TRAJ_EMPTY:
	std::cout << "rep = EMPTY\n";
	break;
    case TRAJ_LINEAR:
	std::cout << "rep = LINEAR\n";
	std::cout << "x   = A*K + B\n";
	std::cout << "v   = C*K + D\n";
	std::cout << "A   = " << _A << "\n";
	std::cout << "B   = " << _B << "\n";
	std::cout << "C   = " << _C << "\n";
	std::cout << "D   = " << _D << "\n";
	break;
    case TRAJ_QUADRATIC:
	std::cout << "rep = QUADRATIC\n";
	std::cout << "x   = A*K^2 + B*K + C\n";
	std::cout << "v   = D*K + E\n";
	std::cout << "A   = " << _A << "\n";
	std::cout << "B   = " << _B << "\n";
	std::cout << "C   = " << _C << "\n";
	std::cout << "D   = " << _D << "\n";
	std::cout << "E   = " << _E << "\n";
	break;
    case TRAJ_CUBIC:
	std::cout << "rep = CUBIC\n";
	std::cout << "x   = A*K^3 + B*K^2 + C*K + D\n";
	std::cout << "v   = (3*A*K^2 + 2*B*K + C)/E\n";
	std::cout << "A   = " << _A << "\n";
	std::cout << "B   = " << _B << "\n";
	std::cout << "C   = " << _C << "\n";
	std::cout << "D   = " << _D << "\n";
	std::cout << "E   = " << _E << "\n";
	break;
    }
}













