/*! \file dxf_solid.cpp
 *  \brief Source code for dxf_solid.hpp
 */

/* Copyright (c) 2010 Taneli Kalvas. All rights reserved.
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
#include "dxf_solid.hpp"
#include "ibsimu.hpp"


DXFSolid::DXFSolid( MyDXFFile *dxffile, const std::string &layername )
{
    if( ibsimu.get_verbose_output() )
	std::cout << "Defining electrode \'" << layername << "\'\n";

    MyDXFEntities *ent = dxffile->get_entities();
    MyDXFEntitySelection *layer = ent->selection_layer( layername );
    if( layer->size() == 0 )
	throw( Error( ERROR_LOCATION, "No entities in layer" ) );
    MyDXFEntitySelection *loop = ent->selection_path_loop( layer );
    if( loop->size() == 0 )
	throw( Error( ERROR_LOCATION, "No loops defined in layer" ) );

    _entities = new MyDXFEntities( ent, loop );
    _selection = _entities->selection_all();

    if( ibsimu.get_verbose_output() ) {
	if( (int)layer->size()-(int)loop->size() > 0 )
	    std::cout << "  removed " << (int)layer->size()-(int)loop->size() << " entities\n";
	std::cout << "  solid defined using " << _entities->size() << " entities\n";
    }

    delete layer;
    delete loop;
}


DXFSolid::~DXFSolid()
{
    if( _entities )
	delete _entities;
    if( _selection )
	delete _selection;
}


bool DXFSolid::inside( const Vec3D &x ) const
{
    return( _entities->inside_loop( _selection, x(0), x(1) ) );
}


void DXFSolid::debug_print( void ) const
{
    for( size_t a = 0; a < _entities->size(); a++ ) {
	std::cout << _entities->get_entity( a );
    }
}


void DXFSolid::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}





