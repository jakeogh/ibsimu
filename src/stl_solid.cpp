/*! \file stl_solid.cpp
 *  \brief %Solid definition using Stereolithography CAD format
 */

/* Copyright (c) 2011,2012 Taneli Kalvas. All rights reserved.
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


#include <iostream>
#include "stlfile.hpp"
#include "stl_solid.hpp"
#include "ibsimu.hpp"


STLSolid::STLSolid( const std::string &filename )
{
    ibsimu.message( 1 ) << "Defining electrode from STL file \'" << filename << "\'\n";
    ibsimu.inc_indent();

    _stl = new STLFile( filename );

    if( _stl ) {
	Vec3D min;
	Vec3D max;
	_stl->get_bbox( min, max );
	ibsimu.message( 1 ) << "bbox min = " << min << "\n";
	ibsimu.message( 1 ) << "bbox max = " << max << "\n";
    }

    ibsimu.dec_indent();
}


STLSolid::STLSolid( std::istream &is )
{
    ibsimu.message( 1 ) << "Warning: loading of STLSolid not implemented\n";
    ibsimu.flush();
}


STLSolid::~STLSolid()
{
    if( _stl )
	delete _stl;
}


bool STLSolid::inside( const Vec3D &x ) const
{
    // Transform 3D -> 3D
    Vec3D y = _T.transform_point( x );

    return( _stl->inside( y ) );
}

STLFile *STLSolid::get_stl_file( void ) const
{
    return( _stl );
}


void STLSolid::debug_print( std::ostream &os ) const
{
    os << "**STLSolid\n";
    //_T.debug_print( os );
    if( _stl )
	_stl->debug_print( os );
    else
	os << "stl = NULL\n";
}


void STLSolid::save( std::ostream &os ) const
{
    write_int32( os, FILEID_STLSOLID );
    ibsimu.message( 1 ) << "Warning: saving of STLSolid not implemented\n";
    ibsimu.flush();
}






