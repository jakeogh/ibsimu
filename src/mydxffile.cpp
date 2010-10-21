#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "mydxffile.hpp"
#include "error.hpp"


MyDXFFile::MyDXFFile( const std::string &filename )
    : _header(0), _entities(0)
{
    // Open file
    _fstr.open( filename.c_str() );
    if( !_fstr.good() )
	throw Error( ERROR_LOCATION, "Couldn't open file \'" + filename + "\'" );

    // Check if binary type
    char buf[22];
    _fstr.get( buf, 22 );
    if( !strncmp( buf, "AutoCAD Binary DXF\x0d\x0a\x1a\0", 22 ) )
	_ascii = false;
    else {
	_linec = 0;
	_ascii = true;
    }
    _fstr.seekg( 0 );
    
    // Search for known section types
    while( read_group() != -1 ) {

	if( group_get_code() == 0 && group_get_string() == "SECTION" ) {

	    if( read_group() != 2 )
		throw Error( ERROR_LOCATION, "Error at string of section on line " + 
			     to_string(_linec) );	

	    if( group_get_string() == "HEADER" ) {
		_header = new MyDXFHeader( this );
	    } else if( group_get_string() == "ENTITIES" ) {
		_entities = new MyDXFEntities( this );
	    } else {
		// Unknown section
#ifdef MYDXF_DEBUG
		std::cout << "Skipping unknown section " << group_get_string() << "\n";
#endif
		while( read_group() != -1 ) {
		    if( group_get_code() == 0 && group_get_string() == "ENDSEC" )
			break;
		}
	    }
	}
    }

    // Close file
    _fstr.close();
}


MyDXFFile::~MyDXFFile()
{
    if( _header )
	delete _header;
    if( _entities )
	delete _entities;
}


int MyDXFFile::group_get_code( void ) const
{
    return( _group_code );
}


#define GROUP_TYPE_UNKNOWN 0
#define GROUP_TYPE_STRING  1
#define GROUP_TYPE_DOUBLE  2
#define GROUP_TYPE_BOOL    3
#define GROUP_TYPE_INT8    4
#define GROUP_TYPE_INT16   5
#define GROUP_TYPE_INT32   6
#define GROUP_TYPE_INT64   7


std::string MyDXFFile::group_get_string( void ) const
{
    if( _group_type != GROUP_TYPE_STRING )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_string );
}


double MyDXFFile::group_get_double( void ) const
{
    if( _group_type != GROUP_TYPE_DOUBLE )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_double );
}


int8_t MyDXFFile::group_get_int8( void ) const
{
    if( _group_type != GROUP_TYPE_INT8 )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_int8 );
}


int16_t MyDXFFile::group_get_int16( void ) const
{
    if( _group_type != GROUP_TYPE_INT16 )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_int16 );
}


int32_t MyDXFFile::group_get_int32( void ) const
{
    if( _group_type != GROUP_TYPE_INT32 )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_int32 );
}


int64_t MyDXFFile::group_get_int64( void ) const
{
    if( _group_type != GROUP_TYPE_INT64 )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_int64 );
}


int MyDXFFile::read_group( void )
{
    char buf[257];
    char *endptr;

    if( _ascii ) {

	// Read group code
	_linec++;
	_fstr.getline( buf, 256 );
	if( _fstr.eof() ) {
	    _group_code = -1;
	    return( -1 );
	}

	_group_code = strtol( buf, &endptr, 10 );
	if( endptr == buf ) {
	    throw Error( ERROR_LOCATION, "Error reading group code on line " + 
			 to_string(_linec) );
	}
	while( isspace( *endptr ) ) endptr++;
	if( *endptr != '\0' ) {
	    throw Error( ERROR_LOCATION, "Error reading group code on line " + 
			 to_string(_linec) );
	}

	// Read group value
	_linec++;
	_fstr.getline( buf, 256 );
	if( _fstr.eof() ) {
	    throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			 to_string(_linec) + ", premature end of file" );
	}
	if( (_group_code >= 0 && _group_code <= 9) ||
	    _group_code == 100 || _group_code == 102 || _group_code == 105 || 
	    (_group_code >= 300 && _group_code <= 369) ||
	    (_group_code >= 390 && _group_code <= 399) ||
	    (_group_code >= 410 && _group_code <= 419) ||
	    (_group_code >= 430 && _group_code <= 439) ||
	    (_group_code >= 470 && _group_code <= 481) ||
	    (_group_code >= 999 && _group_code <= 1009) ) {
	    // String value 
	    _group_type = GROUP_TYPE_STRING;
	    _group_string = buf;
	} else if( (_group_code >= 290 && _group_code <= 299) ) {
	    // bool
	    int value = strtol( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_BOOL;
	    _group_bool = value;
	} else if( false ) {
	    // 8-bit integer value
	    int value = strtol( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_INT8;
	    _group_int8 = value;
	} else if( (_group_code >= 60 && _group_code <= 79) ||
		   (_group_code >= 170 && _group_code <= 179) ||
		   (_group_code >= 270 && _group_code <= 289) ||
		   (_group_code >= 370 && _group_code <= 389) ||
		   (_group_code >= 400 && _group_code <= 409) ||
		   (_group_code >= 1060 && _group_code <= 1070) ) {
	    // 16-bit integer value
	    int value = strtol( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_INT16;
	    _group_int16 = value;
	} else if( (_group_code >= 90 && _group_code <= 99) ||
		   (_group_code >= 420 && _group_code <= 429) ||
		   (_group_code >= 440 && _group_code <= 459) ||
		   _group_code == 1071 ) {
	    // 32-bit integer value
	    int value = strtol( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_INT32;
	    _group_int32 = value;
	} else if( (_group_code >= 160 && _group_code <= 169) ) {
	    // 64-bit integer value
	    int64_t value = strtoll( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_INT64;
	    _group_int64 = value;
	} else if( (_group_code >= 10 && _group_code <= 59) ||
		   (_group_code >= 110 && _group_code <= 149) ||
		   (_group_code >= 210 && _group_code <= 239) ||
		   (_group_code >= 460 && _group_code <= 469) ||
		   (_group_code >= 1010 && _group_code <= 1059) ) {
	    // Double
	    double value = strtod( buf, &endptr );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_DOUBLE;
	    _group_double = value;
	} else {
	    throw Error( ERROR_LOCATION, "Unknown group code " + to_string(_group_code)
			 + " on line " + to_string(_linec) );	
	}
	
    } else { // Read binary

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }

    return( _group_code );
}
