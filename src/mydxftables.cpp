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
    : _flags(0)
{

}


/* ************************************************************************** *
 * DXF Entry BlockRecord                                                      *
 * ************************************************************************** */

MyDXFEntry_BlockRecord::MyDXFEntry_BlockRecord( class MyDXFFile *dxf )
{
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with record

	else if( dxf->group_get_code() == 5 )
    }

}


MyDXFEntry_BlockRecord::~MyDXFEntry_BlockRecord()
{

}


void MyDXFEntry_BlockRecord::write( class MyDXFFile *dxf, std::ofstream &ostr )
{

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

}


void MyDXFTable::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "TABLE" );
    dxf->write_group( 2, _name.c_str() );



    dxf->write_group( 0, "ENDTAB" );
}




void MyDXFTable::debug_print( std::ostream &os ) const
{
    os << "*** Table " << _name << "****************************************\n";


    os << "***************************************************************\n";
}






/* ************************************************************************** *
 * DXF Tables                                                                 *
 * ************************************************************************** */


MyDXFTables::MyDXFTables( class MyDXFFile *dxf )
    : _blockrecord(0)
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading section TABLES\n";
#endif

    // Read TABLES section
    //
    // Ends in ENDSEC

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() != 0 )
	    continue;
	else if( dxf->group_get_string() == "ENDSEC" )
	    break; // Done with tables
	else if( dxf->group_get_string() == "TABLE" ) {
	    if( dxf->read_group() != 2 )
		throw Error( ERROR_LOCATION, "Error at start of table on line " + 
			     to_string(dxf->linec()) );

	    if( dxf->group_get_string() == "BLOCK_RECORD" ) {
		_blockrecord = new MyDXFTable( "BLOCK_RECORD", dxf );
	    } else {
		// Unknown table
#ifdef MYDXF_DEBUG
		std::cout << "Skipping unknown table " << dxf->group_get_string() << "\n";
#endif
		while( dxf->read_group() != -1 ) {
		    if( dxf->group_get_code() == 0 && dxf->group_get_string() == "ENDTAB" )
			break;
		}
	    }
	}
    }
}



MyDXFTables::~MyDXFTables()
{
    if( _blockrecord ) 
	delete _blockrecord;
}


void MyDXFTables::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "SECTION" );
    dxf->write_group( 2, "TABLES" );

    _blockrecord->write( dxf, ostr );

    dxf->write_group( 0, "ENDSEC" );
}




void MyDXFTables::debug_print( std::ostream &os ) const
{
    os << "*** Section TABLES ****************************************\n";


    os << "***********************************************************\n";
}

