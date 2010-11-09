#include <iostream>
#include "mydxfheader.hpp"
#include "error.hpp"


MyDXFHeader::MyDXFHeader( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading section HEADER\n";
#endif

    // Read HEADER section
    //
    // Contains pairs of groups:
    // A. Variable: code=9, data=header name (string)
    // B. Header data: code and data according to field
    //
    // Ends in ENDSEC

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 && dxf->group_get_string() == "ENDSEC" )
	    break; // Done with header
	if( dxf->group_get_code() != 9 )
	    continue; // Skip unknown input

	// Read data group
	std::string field = dxf->group_get_string();
	dxf->read_group();

	if( dxf->group_get_code() == 1 && field == "$ACADVER" ) {
	    acadver = dxf->group_get_string();
#ifdef MYDXF_DEBUG
	    std::cout << "  acadver = " << acadver << "\n";
#endif
	} else if( dxf->group_get_code() == 50 && field == "$ANGBASE" ) {
	    angdir = dxf->group_get_double();
#ifdef MYDXF_DEBUG
	    std::cout << "  angbase = " << angbase << "\n";
#endif
	} else if( dxf->group_get_code() == 70 && field == "$ANGDIR" ) {
	    angdir = dxf->group_get_int16();
#ifdef MYDXF_DEBUG
	    std::cout << "  angdir = " << angdir << "\n";
#endif
	}
    }

}



MyDXFHeader::~MyDXFHeader()
{

}



void MyDXFHeader::debug_print( std::ostream &os ) const
{
    os << "*** Section HEADER ****************************************\n";

    os << "  acadver = " << acadver << "\n";
    os << "  angbase = " << angbase << "\n";
    os << "  angdir = " << angdir << "\n";

    os << "\n";
}
