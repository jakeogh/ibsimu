/*! \file mydxftables.cpp
 *  \brief DXF Tables 
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
#include "mydxftables.hpp"
#include "error.hpp"


#define MYDXF_DEBUG 1


/* ************************************************************************** *
 * DXF Entry                                                                  *
 * ************************************************************************** */

MyDXFEntry::MyDXFEntry()
{

}


void MyDXFEntry::process_group( class MyDXFFile *dxf )
{
    if( dxf->group_get_code() == 5 || dxf->group_get_code() == 105 )
	_handle = dxf->group_get_string();
    if( dxf->group_get_code() == 330 )
	_handle_to_owner = dxf->group_get_string();
}


void MyDXFEntry::write_common( class MyDXFFile *dxf, std::ofstream &ostr )
{
    //if( dynamic_cast<MyDXFEntry_DimStyle *>( this ) )
    //dxf->write_group( 105, _handle.c_str() );
    //else
    dxf->write_group( 5, _handle.c_str() );
    //dxf->write_group( 330, _handle_to_owner.c_str() );
}


void MyDXFEntry::debug_print_common( std::ostream &os ) const
{
    os << "  handle = \'" << _handle << "\'\n";
    os << "  handle_to_owner = \'" << _handle_to_owner << "\'\n";
}


std::ostream &operator<<( std::ostream &os, const MyDXFEntry &e )
{
    e.debug_print( os );
    e.debug_print_common( os );
    return( os );
}


/* ************************************************************************** *
 * DXF Entry BlockRecord                                                      *
 * ************************************************************************** */

MyDXFEntry_BlockRecord::MyDXFEntry_BlockRecord( class MyDXFFile *dxf )
    : _units(0), _explodability(0), _scalability(0)
{
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with record

	else if( dxf->group_get_code() == 2 )
	    _name = dxf->group_get_string();
	else if( dxf->group_get_code() == 70 )
	    _units = dxf->group_get_int16();
	else if( dxf->group_get_code() == 340 )
	    _handle_to_layout = dxf->group_get_string();
	else if( dxf->group_get_code() == 280 )
	    _explodability = dxf->group_get_int8();
	else if( dxf->group_get_code() == 281 )
	    _scalability = dxf->group_get_int8();

	else
	    process_group( dxf );
    }

}


MyDXFEntry_BlockRecord::~MyDXFEntry_BlockRecord()
{

}


void MyDXFEntry_BlockRecord::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "BLOCK_RECORD" );
    write_common( dxf, ostr );

    dxf->write_group( 100, "AcDbSymbolTableRecord" );
    dxf->write_group( 100, "AcDbBlockTableRecord" );

    dxf->write_group( 2, _name.c_str() );
    //dxf->write_group( 70, _units );
    //dxf->write_group( 280, _explodability );
    //dxf->write_group( 281, _scalability );
    dxf->write_group( 340, _handle_to_layout.c_str() );
}


void MyDXFEntry_BlockRecord::debug_print( std::ostream &os ) const
{
    os << "  name = \'" << _name << "\'\n";
    os << "  units = " << _units << "\n";
    os << "  explodability = " << (int)_explodability << "\n";
    os << "  scalability = " << (int)_scalability << "\n";
    os << "  handle_to_layout = \'" << _handle_to_layout << "\'\n";
}


/* ************************************************************************** *
 * DXF Entry Layer                                                            *
 * ************************************************************************** */

MyDXFEntry_Layer::MyDXFEntry_Layer( class MyDXFFile *dxf )
    : _flags(0), _color(0), _plotting(true), _lineweight(0)
{
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with record

	else if( dxf->group_get_code() == 2 )
	    _name = dxf->group_get_string();
	else if( dxf->group_get_code() == 6 )
	    _linetype = dxf->group_get_string();
	else if( dxf->group_get_code() == 70 )
	    _flags = dxf->group_get_int16();
	else if( dxf->group_get_code() == 62 )
	    _color = dxf->group_get_int16();
	else if( dxf->group_get_code() == 290 )
	    _plotting = dxf->group_get_bool();
	else if( dxf->group_get_code() == 370 )
	    _lineweight = dxf->group_get_int8();
	else if( dxf->group_get_code() == 390 )
	    _handle_to_plot_style_name = dxf->group_get_string();
	else if( dxf->group_get_code() == 347 )
	    _handle_to_material = dxf->group_get_string();

	else
	    process_group( dxf );
    }

}


MyDXFEntry_Layer::~MyDXFEntry_Layer()
{

}


void MyDXFEntry_Layer::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "LAYER" );
    write_common( dxf, ostr );

    dxf->write_group( 100, "AcDbSymbolTableRecord" );
    dxf->write_group( 100, "AcDbLayerTableRecord" );

    dxf->write_group( 2, _name.c_str() );
    dxf->write_group( 6, _linetype.c_str() );
    dxf->write_group( 70, _flags );
    dxf->write_group( 62, _color );
    dxf->write_group( 290, _plotting );
    dxf->write_group( 370, _lineweight );
    dxf->write_group( 390, _handle_to_plot_style_name.c_str() );
    dxf->write_group( 347, _handle_to_material.c_str() );
}


void MyDXFEntry_Layer::debug_print( std::ostream &os ) const
{
    os << "  name = \'" << _name << "\'\n";
    os << "  linetype = \'" << _linetype << "\'\n";
    os << "  flags = " << _flags << "\n";
    os << "  color = " << _color << "\n";
    os << "  plotting = " << _plotting << "\n";
    os << "  lineweight = " << (int)_lineweight << "\n";
    os << "  handle_to_plot_style_name = \'" << _handle_to_plot_style_name << "\'\n";
    os << "  handle_to_layout = \'" << _handle_to_material << "\'\n";
}


/* ************************************************************************** *
 * DXF Table                                                                  *
 * ************************************************************************** */

MyDXFTable::MyDXFTable( const std::string &name, class MyDXFFile *dxf )
    : _name(name)
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading Table " << _name << "\n";
#endif

    while( dxf->group_get_code() != -1 ) {
	
	if( dxf->group_get_code() == 0 ) {
	    if( dxf->group_get_string() == "ENDSEC" || dxf->group_get_string() == "ENDTAB" )
		break; // Done with table
	    if( dxf->group_get_string() == _name ) {
		if( dxf->group_get_string() == "BLOCK_RECORD" ) 
		    _entries.push_back( new MyDXFEntry_BlockRecord( dxf ) );
		else if( dxf->group_get_string() == "LAYER" ) 
		    _entries.push_back( new MyDXFEntry_Layer( dxf ) );
		else
		    throw Error( ERROR_LOCATION, "Error reading table entry on line " + 
				 to_string(dxf->linec()) );
	    }
	} else {

	    if( dxf->group_get_code() == 5 )
		_handle = dxf->group_get_string();
	    else if( dxf->group_get_code() == 330 )
		_handle_to_owner = dxf->group_get_string();

	    dxf->read_group();
	}
    }
}



MyDXFTable::~MyDXFTable()
{
    for( uint32_t i = 0; i < _entries.size(); i++ )
	delete _entries[i];
}


void MyDXFTable::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "TABLE" );
    dxf->write_group( 2, _name.c_str() );

    dxf->write_group( 5, _handle.c_str() );

    dxf->write_group( 100, "AcDbSymbolTable" );

    //dxf->write_group( 330, _handle_to_owner.c_str() );
    dxf->write_group( 70, (int16_t)_entries.size() );

    for( uint32_t i = 0; i < _entries.size(); i++ )
	_entries[i]->write( dxf, ostr );

    dxf->write_group( 0, "ENDTAB" );
}




void MyDXFTable::debug_print( std::ostream &os ) const
{
    os << "*** Table " << _name << " ****************************************\n";

    os << "  handle = \'" << _handle << "\'\n";
    os << "  handle_to_owner = \'" << _handle_to_owner << "\'\n";
    os << "\n";

    for( uint32_t i = 0; i < _entries.size(); i++ )
	os << (*_entries[i]) << "\n";
}






/* ************************************************************************** *
 * DXF Tables                                                                 *
 * ************************************************************************** */


MyDXFTables::MyDXFTables( class MyDXFFile *dxf )
    : _blockrecord(0), _layer(0)
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading section TABLES\n";
#endif

    // Read TABLES section
    //
    // Ends in ENDSEC
    dxf->read_group();

    while( dxf->group_get_code() != -1 ) {
	
	if( dxf->group_get_code() != 0 ) {
	    dxf->read_group();
	    continue;
	}

	else if( dxf->group_get_string() == "ENDSEC" )
	    break; // Done with tables
	else if( dxf->group_get_string() == "TABLE" ) {
	    if( dxf->read_group() != 2 )
		throw Error( ERROR_LOCATION, "Error at start of table on line " + 
			     to_string(dxf->linec()) );

	    if( dxf->group_get_string() == "BLOCK_RECORD" )
		_blockrecord = new MyDXFTable( "BLOCK_RECORD", dxf );
	    else if( dxf->group_get_string() == "LAYER" )
		_layer = new MyDXFTable( "LAYER", dxf );
	    else {
		// Unknown table
#ifdef MYDXF_DEBUG
		std::cout << "Skipping unknown table " << dxf->group_get_string() << "\n";
#endif
		while( dxf->read_group() != -1 ) {
		    if( dxf->group_get_code() == 0 && dxf->group_get_string() == "ENDTAB" )
			break;
		}
		dxf->read_group();
	    }
	} else {
	    dxf->read_group();
	}
    }
}



MyDXFTables::~MyDXFTables()
{
    if( _layer ) 
	delete _layer;
    if( _blockrecord ) 
	delete _blockrecord;
}


void MyDXFTables::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "SECTION" );
    dxf->write_group( 2, "TABLES" );

    if( _blockrecord ) 
	_blockrecord->write( dxf, ostr );
    if( _layer ) 
	_layer->write( dxf, ostr );

    dxf->write_group( 0, "ENDSEC" );
}




void MyDXFTables::debug_print( std::ostream &os ) const
{
    os << "*** Section TABLES ****************************************\n";

    if( _blockrecord ) {
	_blockrecord->debug_print( os );
	os << "\n";
    }

    if( _layer ) {
	_layer->debug_print( os );
	os << "\n";
    }
}

