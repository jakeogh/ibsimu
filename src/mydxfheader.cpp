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
    // Contains packages of groups:
    // A. Variable name: code=9, data=header name (string)
    // B. Variable data: code(s) and data according to field(s)
    //
    // Ends in ENDSEC

    std::string field;
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 && dxf->group_get_string() == "ENDSEC" )
	    break; // Done with header

	else if( dxf->group_get_code() == 9 )
	    field = dxf->group_get_string(); // Read header variable name

	else if( dxf->group_get_code() == 1 && field == "$ACADVER" )
	    acadver = dxf->group_get_string();
	else if( dxf->group_get_code() == 50 && field == "$ANGBASE" )
	    angdir = dxf->group_get_double();
	else if( dxf->group_get_code() == 70 && field == "$ANGDIR" )
	    angdir = dxf->group_get_int16();
	else if( dxf->group_get_code() == 70 && field == "$INSUNITS" )
	    insunits = dxf->group_get_int16();
	else if( dxf->group_get_code() == 5 && field == "$HANDSEED" )
	    handseed = dxf->group_get_string();
	else if( dxf->group_get_code() == 40 && field == "$DIMASZ" )
	    dimasz = dxf->group_get_double();
	else if( dxf->group_get_code() == 40 && field == "$DIMGAP" )
	    dimgap = dxf->group_get_double();
	else if( dxf->group_get_code() == 40 && field == "$DIMEXO" )
	    dimexo = dxf->group_get_double();
	else if( dxf->group_get_code() == 40 && field == "$DIMEXE" )
	    dimexe = dxf->group_get_double();
	else if( dxf->group_get_code() == 40 && field == "$DIMTXT" )
	    dimtxt = dxf->group_get_double();

	else if( dxf->group_get_code() == 10 && field == "$PLIMMAX" )
	    plimmax[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PLIMMAX" )
	    plimmax[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 10 && field == "$PLIMMIN" )
	    plimmin[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 10 && field == "$PLIMMIN" )
	    plimmin[1] = dxf->group_get_double();
    }

}



MyDXFHeader::~MyDXFHeader()
{

}


void MyDXFHeader::write( class MyDXFFile *dxf, std::ofstream &_ostr )
{
    dxf->write_group( 0, "SECTION" );
    dxf->write_group( 2, "HEADER" );

    dxf->write_group( 9, "$ACADVER" );
    dxf->write_group( 1, acadver.c_str() );

    dxf->write_group( 9, "$ANGBASE" );
    dxf->write_group( 50, angbase );

    dxf->write_group( 9, "$ANGDIR" );
    dxf->write_group( 70, angdir );

    dxf->write_group( 9, "$HANDSEED" );
    dxf->write_group( 5, handseed.c_str() );

    dxf->write_group( 9, "$DIMASZ" );
    dxf->write_group( 40, dimasz );

    dxf->write_group( 9, "$DIMGAP" );
    dxf->write_group( 40, dimgap );

    dxf->write_group( 9, "$DIMEXO" );
    dxf->write_group( 40, dimexo );

    dxf->write_group( 9, "$DIMEXE" );
    dxf->write_group( 40, dimexe );

    dxf->write_group( 9, "$DIMTXT" );
    dxf->write_group( 40, dimtxt );

    dxf->write_group( 9, "$INSUNITS" );
    dxf->write_group( 70, insunits );

    dxf->write_group( 9, "$PLIMMAX" );
    dxf->write_group( 10, plimmax[0] );
    dxf->write_group( 20, plimmax[1] );

    dxf->write_group( 9, "$PLIMMIN" );
    dxf->write_group( 10, plimmin[0] );
    dxf->write_group( 20, plimmin[1] );

    dxf->write_group( 0, "ENDSEC" );
}



void MyDXFHeader::debug_print( std::ostream &os ) const
{
    os << "*** Section HEADER ****************************************\n";

    os << "  acadver = \'" << acadver << "\'\n";
    os << "  angbase = " << angbase << "\n";
    os << "  angdir = " << angdir << "\n";

    os << "  handseed = \'" << handseed << "\'\n";
    os << "  dimasz = " << dimasz << "\n";
    os << "  dimgap = " << dimgap << "\n";
    os << "  dimexo = " << dimexo << "\n";
    os << "  dimexe = " << dimexe << "\n";
    os << "  dimtxt = " << dimtxt << "\n";
    os << "  insunits = " << insunits << "\n";
    os << "  plimmax = " << plimmax << "\n";
    os << "  plimmin = " << plimmin << "\n";

    os << "\n";
}



