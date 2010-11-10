#include <iostream>
#include "mydxfblocks.hpp"
#include "error.hpp"



MyDXFBlock::MyDXFBlock( class MyDXFFile *dxf )
    : _entities(0)
{
#ifdef MYDXF_DEBUG
    std::cout << "*** Reading BLOCK ***\n";
#endif

    // Load default values
    _type = 0;

    while( dxf->read_group() != -1 ) {

	if( dxf->group_get_code() == 0 ) {

	    // Check if done with block
	    if( dxf->group_get_string() == "ENDBLK" ||
		dxf->group_get_string() == "ENDSEC" )
		break;

	    // Read entities
	    _entities = new MyDXFEntities( dxf, true );
	    break;
	}

	else if( dxf->group_get_code() == 1 )
	    _path = dxf->group_get_string();
	else if( dxf->group_get_code() == 2 || dxf->group_get_code() == 3 )
	    _name = dxf->group_get_string();
	else if( dxf->group_get_code() == 5 )
	    _handle = dxf->group_get_string();
	else if( dxf->group_get_code() == 8 )
	    _layer = dxf->group_get_string();
	else if( dxf->group_get_code() == 70 )
	    _type = dxf->group_get_int16();
	else if( dxf->group_get_code() == 330 )
	    _owner_handle = dxf->group_get_string();

	else if( dxf->group_get_code() == 10 )
	    _p[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _p[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _p[2] = dxf->group_get_double();
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


MyDXFBlock::~MyDXFBlock()
{
    if( _entities )
	delete _entities;
}


std::ostream &operator<<( std::ostream &os, const MyDXFBlock &blk )
{
    os << "BLOCK\n";
    os << "  path = \'" << blk._path << "\'\n";
    os << "  handle = \'" << blk._handle << "\'\n";
    os << "  layer = \'" << blk._layer << "\'\n";
    os << "  owner_handle = \'" << blk._owner_handle << "\'\n";
    os << "  name = \'" << blk._name << "\'\n";
    os << "  type = " << blk._type << "\n";
    os << "  p = " << blk._p << "\n";

    if( blk._entities ) {
	for( size_t a = 0; a < blk._entities->size(); a++ ) {
	    MyDXFEntity *e = blk._entities->get_entity(a);
	    os << *e;
	}
    }

    os << "ENDBLK\n\n";

    return( os );
}


/* ************************************************************************** *
 * DXFBlocks                                                                  *
 * ************************************************************************** */


MyDXFBlocks::MyDXFBlocks( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading section BLOCKS\n";
#endif

    // Read BLOCKS section
    // 
    // Block entities between (0,BLOCK) and (0,ENDBLK) 
    //
    // Ends in ENDSEC

    dxf->read_group();
    while( dxf->group_get_code() != -1 ) {
	
	if( dxf->group_get_code() == 0 && dxf->group_get_string() == "ENDSEC" )
	    break; // Done with blocks
	if( dxf->group_get_code() != 0 ) {
	    dxf->read_group();
	    continue; // Skip unknown input
	}

	// Check for blocks
	if( dxf->group_get_string() == "BLOCK" ) {
	    _blocks.push_back( new MyDXFBlock(dxf) );
	} else {
	    dxf->read_group();
        }
    }
}



MyDXFBlocks::~MyDXFBlocks()
{

}



void MyDXFBlocks::debug_print( std::ostream &os ) const
{
    os << "*** Section BLOCKS ****************************************\n";

    for( size_t a = 0; a < _blocks.size(); a++ ) {
	MyDXFBlock *b = _blocks[a];
	os << *b;
    }

    os << "\n";
}



