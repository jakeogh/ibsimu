#include <iostream>
#include <limits>
#include "mydxfblocks.hpp"
#include "error.hpp"



MyDXFBlock::MyDXFBlock( class MyDXFFile *dxf )
    : _type(0), _entities(0)
{
#ifdef MYDXF_DEBUG
    std::cout << "\n*** Reading BLOCK ***\n";
#endif

    while( dxf->read_group() != -1 ) {

	if( dxf->group_get_code() == 0 ) {
	    if( dxf->group_get_string() != "ENDBLK" && dxf->group_get_string() != "ENDSEC" )
		_entities = new MyDXFEntities( dxf, true );
	    break; // Done with block
	}

	else if( dxf->group_get_code() == 1 )
	    _path = dxf->group_get_string();
	else if( dxf->group_get_code() == 2 || dxf->group_get_code() == 3 )
	    _name = dxf->group_get_string();
	else if( dxf->group_get_code() == 5 )
	    _block_handle = dxf->group_get_string();
	else if( dxf->group_get_code() == 8 )
	    _block_layer = dxf->group_get_string();

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


    if( dxf->group_get_string() == "ENDBLK" ) {
	// Read end of block
	while( dxf->read_group() != -1 ) {
	    if( dxf->group_get_code() == 5 )
		_endblk_handle = dxf->group_get_string();
	    else if( dxf->group_get_code() == 8 )
		_endblk_layer = dxf->group_get_string();
	    else if( dxf->group_get_code() == 0 )
		break; // Done with end of block
	}
    }

#ifdef MYDXF_DEBUG
    std::cout << "\n*** Done with BLOCK ***\n";
    //std::cout << *this;
#endif
}


MyDXFBlock::~MyDXFBlock()
{
    if( _entities )
	delete _entities;
}


void MyDXFBlock::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "BLOCK" );

    dxf->write_group( 1, _path.c_str() );
    dxf->write_group( 2, _name.c_str() );
    dxf->write_group( 3, _name.c_str() );
    dxf->write_group( 5, _block_handle.c_str() );
    dxf->write_group( 8, _block_layer.c_str() );

    dxf->write_group( 70, _type );
    dxf->write_group( 330, _owner_handle.c_str() );

    dxf->write_group( 10, _p[0] );
    dxf->write_group( 20, _p[0] );
    dxf->write_group( 30, _p[0] );

    // Write entities within block
    if( _entities )
	_entities->write_entities( dxf, ostr );

    dxf->write_group( 0, "ENDBLK" );
    dxf->write_group( 5, _endblk_handle.c_str() );
    dxf->write_group( 8, _endblk_layer.c_str() );
}


void MyDXFBlock::explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const
{
    if( !_entities )
	return;

    Transformation t2 = *t;
    t2.translate( -_p );
    _entities->explode( ent, dxf, &t2 );
}


void MyDXFBlock::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const
{
    if( !_entities )
	return;

    Transformation t2 = *t;
    t2.translate( -_p );
    _entities->plot( NULL, dxf, cairo, &t2, range );
}


void MyDXFBlock::get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const
{
    if( !_entities ) {
	min = Vec3D( std::numeric_limits<double>::infinity(),
		     std::numeric_limits<double>::infinity(),
		     std::numeric_limits<double>::infinity() );
	max = Vec3D( -std::numeric_limits<double>::infinity(),
		     -std::numeric_limits<double>::infinity(),
		     -std::numeric_limits<double>::infinity() );
	return;
    }

    Transformation t2 = *t;
    t2.translate( -_p );
    _entities->get_bbox( NULL, min, max, dxf, &t2 );
}


void MyDXFBlock::scale( class MyDXFFile *dxf, double s )
{
    if( !_entities )
	return;

    _entities->scale( NULL, dxf, s );
}

std::ostream &operator<<( std::ostream &os, const MyDXFBlock &blk )
{
    os << "BLOCK\n";
    os << "  path = \'" << blk._path << "\'\n";
    os << "  block_handle = \'" << blk._block_handle << "\'\n";
    os << "  block_layer = \'" << blk._block_layer << "\'\n";
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

    os << "ENDBLK\n";
    os << "  endblk_handle = \'" << blk._endblk_handle << "\'\n";
    os << "  endblk_layer = \'" << blk._endblk_layer << "\'\n";
    os << "\n";

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
    for( size_t i = 0; i < _blocks.size(); i++ )
	delete _blocks[i];
}


void MyDXFBlocks::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "SECTION" );
    dxf->write_group( 2, "BLOCKS" );

    for( size_t a = 0; a < _blocks.size(); a++ )
	_blocks[a]->write( dxf, ostr );

    dxf->write_group( 0, "ENDSEC" );
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


MyDXFBlock *MyDXFBlocks::get_by_name( const std::string &name )
{
    for( int i = 0; i < (int)_blocks.size(); i++ ) {
	if( _blocks[i]->name() == name ) {
	    return( _blocks[i] );
	}
    }

    return( 0 );
}


const MyDXFBlock *MyDXFBlocks::get_by_name( const std::string &name ) const
{
    for( int i = 0; i < (int)_blocks.size(); i++ ) {
	if( _blocks[i]->name() == name ) {
	    return( _blocks[i] );
	}
    }

    return( 0 );
}


void MyDXFBlocks::clear( void )
{
    for( size_t i = 0; i < _blocks.size(); i++ )
	delete _blocks[i];
    _blocks.resize( 0 );
}
