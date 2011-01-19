/*! \file convergence.cpp
 *  \brief Vlasov system convergence tester
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


#include "convergence.hpp"
#include "ibsimu.hpp"


Convergence::EmitPoint::EmitPoint()
    : _epsilon(0.0), _alpha(0.0), _beta(0.0), _xave(0.0), _xpave(0.0)
{
}


Convergence::EmitPoint::EmitPoint( double x )
    : _epsilon(x), _alpha(x), _beta(x), _xave(x), _xpave(x)
{
}


Convergence::EmitPoint::EmitPoint( const Emittance &emit )
    : _epsilon(emit.epsilon()), _alpha(emit.alpha()), _beta(emit.beta()), 
      _xave(emit.xave()), _xpave(emit.xpave())
{    
}


double &Convergence::EmitPoint::operator[]( int i )
{
    if( i == 0 ) 
	return( _epsilon );
    else if( i == 1 ) 
	return( _alpha );
    else if( i == 2 ) 
	return( _beta );
    else if( i == 3 ) 
	return( _xave );
    else if( i == 4 ) 
	return( _xpave );
    else
	throw( Error( ERROR_LOCATION, "invalid index" ) );
}


Convergence::Convergence()
    : _iter(0), _epot_old(NULL), _epot(NULL), _scharge_old(NULL), _scharge(NULL),
      _emit_pdb(NULL)
{

}


Convergence::~Convergence()
{
    if( _epot_old )
	delete _epot_old;
    if( _scharge_old )
	delete _scharge_old;
}


bool Convergence::evaluate_iteration( void )
{
    bool convergence = true;
    
    if( _iter == 0 )
	convergence = false;
    
    if( ibsimu.get_verbose_output() )
	std::cout << "Iteration round " << _iter << "\n";

    // Calculate convergence estimate for epot
    if( _epot ) {
	if( !_epot_old ) {
	    // First round 
	    _epot_hist.push_back( 0.0 );
	    _epot_old = new ScalarField( *_epot );
	} else {
	    // Evaluate error estimate
	    double error = 0.0;
	    int nc = _epot->nodecount();
	    for( int a = 0; a < nc; a++ ) {
		double abs = _epot_absf*fabs( (*_epot)(a) - (*_epot_old)(a) );
		double rel = _epot_relf*fabs( ((*_epot)(a) - (*_epot_old)(a))/(*_epot)(a) );
		if( abs < rel ) {
		    if( abs > error )
			error = abs;
		} else {
		    if( rel > error )
			error = rel;
		}
	    }
	    _epot_hist.push_back( error );
	    delete _epot_old;
	    _epot_old = new ScalarField( *_epot );

	    if( error > 1.0 )
		convergence = false;
	}

	if( ibsimu.get_verbose_output() )
	    std::cout << "  Epot error = " << _epot_hist.back() << "\n";
    }
    
    // Calculate convergence estimate for scharge
    if( _scharge ) {
	if( !_scharge_old ) {
	    // First round 
	    _scharge_hist.push_back( 0.0 );
	    _scharge_old = new ScalarField( *_scharge );
	} else {
	    // Evaluate error estimate
	    double error = 0.0;
	    int nc = _scharge->nodecount();
	    for( int a = 0; a < nc; a++ ) {
		double abs = _scharge_absf*fabs( (*_scharge)(a) - (*_scharge_old)(a) );
		double rel = _scharge_relf*fabs( ((*_scharge)(a) - (*_scharge_old)(a))/(*_scharge)(a) );
		if( abs < rel ) {
		    if( abs > error )
			error = abs;
		} else {
		    if( rel > error )
			error = rel;
		}
	    }
	    _scharge_hist.push_back( error );
	    delete _scharge_old;
	    _scharge_old = new ScalarField( *_scharge );

	    if( error > 1.0 )
		convergence = false;
	}

	if( ibsimu.get_verbose_output() )
	    std::cout << "  Space charge error = " << _scharge_hist.back() << "\n";
    }
    
    // Calculate convergence estimate for emittances
    EmitPoint emitp[2];
    int verbose_back = ibsimu.get_verbose_output();
    ibsimu.set_verbose_output( 0 );
    if( _emit_pdb ) {
	// Make trajectory diagnostics
	// NOTE: Trajectory diagnostics most probably not done like 
	// what user sees in plots. The trajectory data is not mirrored.
	// Mirroring would require mesh-type information, which is currently
	// not known in particle database.
	TrajectoryDiagnosticData tdata;
	std::vector<trajectory_diagnostic_e> diagnostics;
	if( _emit_pdb->geom_mode() == MODE_2D ) {

	    if( _emit_axis == AXIS_X ) {
		diagnostics.push_back( DIAG_Y );
		diagnostics.push_back( DIAG_YP );
		diagnostics.push_back( DIAG_CURR );
	    } else {
		diagnostics.push_back( DIAG_X );
		diagnostics.push_back( DIAG_XP );
		diagnostics.push_back( DIAG_CURR );
	    }
	    _emit_pdb->trajectories_at_plane( tdata, _emit_axis, _emit_val, diagnostics );
	    Emittance emit( tdata(0).data(), tdata(1).data(), tdata(2).data() );
	    emitp[0] = EmitPoint( emit );

	} else if( _emit_pdb->geom_mode() == MODE_CYL ) {

	    if( _emit_axis == AXIS_X ) {
		diagnostics.push_back( DIAG_R );
		diagnostics.push_back( DIAG_RP );
		diagnostics.push_back( DIAG_CURR );
	    } else {
		diagnostics.push_back( DIAG_X );
		diagnostics.push_back( DIAG_XP );
		diagnostics.push_back( DIAG_CURR );
	    }
	    _emit_pdb->trajectories_at_plane( tdata, _emit_axis, _emit_val, diagnostics );
	    Emittance emit( tdata(0).data(), tdata(1).data(), tdata(2).data() );
	    emitp[0] = EmitPoint( emit );

	} else {

	    if( _emit_axis == AXIS_X ) {
		diagnostics.push_back( DIAG_Y );
		diagnostics.push_back( DIAG_YP );
		diagnostics.push_back( DIAG_Z );
		diagnostics.push_back( DIAG_ZP );
		diagnostics.push_back( DIAG_CURR );
	    } else if( _emit_axis == AXIS_Y ) {
		diagnostics.push_back( DIAG_X );
		diagnostics.push_back( DIAG_XP );
		diagnostics.push_back( DIAG_Z );
		diagnostics.push_back( DIAG_ZP );
		diagnostics.push_back( DIAG_CURR );
	    } else {
		diagnostics.push_back( DIAG_X );
		diagnostics.push_back( DIAG_XP );
		diagnostics.push_back( DIAG_Y );
		diagnostics.push_back( DIAG_YP );
		diagnostics.push_back( DIAG_CURR );
	    }
	    _emit_pdb->trajectories_at_plane( tdata, _emit_axis, _emit_val, diagnostics );
	    Emittance emit1( tdata(0).data(), tdata(1).data(), tdata(4).data() );
	    emitp[0] = EmitPoint( emit1 );
	    Emittance emit2( tdata(2).data(), tdata(3).data(), tdata(4).data() );
	    emitp[1] = EmitPoint( emit2 );

	}
    } else {
	if( _emit[0] )
	    emitp[0] = EmitPoint( *_emit[0] );
	if( _emit[1] )
	    emitp[1] = EmitPoint( *_emit[1] );
    }
    ibsimu.set_verbose_output( verbose_back );

    for( int i = 0; i < 2; i++ ) {
	if( !_emit[i] && !_emit_pdb )
	    continue;
	if( i == 1 && _emit_pdb && _emit_pdb->geom_mode() != MODE_3D  )
	    continue;
	else if( _iter == 0 ) {
	    // First round
	    _emit_hist[i].push_back( 0.0 );
	    _emit_old[i] = emitp[i];
	} else {
	    // Evaluate error estimate
	    double error = 0.0;
	    for( int a = 0; a < 5; a++ ) {
		double abs = _emit_absf*fabs( emitp[i][a] - _emit_old[i][a] );
		double rel = _emit_relf*fabs( (emitp[i][a] - _emit_old[i][a]) / emitp[i][a] );
		if( abs < rel ) {
		    if( abs > error )
			error = abs;
		} else {
		    if( rel > error )
			error = rel;
		}
	    }
	    _emit_hist[i].push_back( error );
	    _emit_old[i] = emitp[i];

	    if( error > 1.0 )
		convergence = false;
	}

	if( ibsimu.get_verbose_output() )
	    std::cout << "  Emittance " << i << " error = " << _emit_hist[i].back() << "\n";
    }

    // Increase iteration counter
    _iter++;

    return( convergence );
}


void Convergence::print_history( std::ostream &os ) const
{
    // Otsikot!

    for( int a = 1; a < _iter; a++ ) {
	os << std::setw(4) << a << " ";
	if( _epot )
	    os << std::setw(12) << _epot_hist[a] << " ";
	if( _scharge )
	    os << std::setw(12) << _scharge_hist[a] << " ";
	if( _emit_hist[0].size() >= (unsigned int)a )
	    os << std::setw(12) << _emit_hist[0][a] << " ";
	if( _emit_hist[1].size() >= (unsigned int)a )
	    os << std::setw(12) << _emit_hist[1][a] << " ";
	os << "\n";
    }
}


void Convergence::add_epot( const ScalarField &epot, double absf, double relf, double lim )
{
    _epot = &epot;
    _epot_old = NULL;
    _epot_hist.clear();
    _epot_absf = fabs(absf);
    _epot_relf = fabs(relf);
    _epot_lim = fabs(lim);
    _iter = 0;
}


void Convergence::add_scharge( const ScalarField &scharge, double absf, double relf, double lim )
{
    _scharge = &scharge;
    _scharge_old = NULL;
    _scharge_hist.clear();
    _scharge_absf = fabs(absf);
    _scharge_relf = fabs(relf);
    _scharge_lim = fabs(lim);
    _iter = 0;
}


void Convergence::add_tdiag( const ParticleDataBase &pdb, coordinate_axis_e axis,
			     double val, double absf, double relf, double lim )
{
    _emit_pdb = &pdb;
    _emit_axis = axis;
    _emit_val = val;
    _emit_hist[0].clear();
    _emit_hist[1].clear();
    _emit_absf = fabs(absf);
    _emit_relf = fabs(relf);
    _emit_lim = fabs(lim);
    _emit[0] = NULL;
    _emit[1] = NULL;
    _iter = 0;
}


void Convergence::add_tdiag( const Emittance *emit1, const Emittance *emit2, 
			     double absf, double relf, double lim )
{
    _emit_pdb = NULL;
    _emit_hist[0].clear();
    _emit_hist[1].clear();
    _emit_absf = fabs(absf);
    _emit_relf = fabs(relf);
    _emit_lim = fabs(lim);
    _emit[0] = emit1;
    _emit[1] = emit2;
    _iter = 0;
}


void Convergence::set_emittance( const Emittance *emit1, const Emittance *emit2 )
{
    _emit[0] = emit1;
    _emit[1] = emit2;
}
