#include "particleiterator.hpp"


ParticleStatistics::ParticleStatistics()
    : _nboundaries(0), _end_time(0), _end_step(0), _end_baddef(0), _sum_steps(0)
{
    _bound_collisions = NULL;
    _bound_current = NULL;
}


ParticleStatistics::ParticleStatistics( uint32_t nboundaries )
    : _nboundaries(nboundaries), _end_time(0), _end_step(0), _end_baddef(0), _sum_steps(0)
{
    if( _nboundaries ) {
	_bound_collisions = new uint32_t[_nboundaries];
	_bound_current = new double[_nboundaries];
	for( uint32_t a = 0; a < _nboundaries; a++ ) {
	    _bound_collisions[a] = 0;
	    _bound_current[a] = 0.0;
	}
    } else {
	_bound_collisions = NULL;
	_bound_current = NULL;
    }
}


ParticleStatistics::ParticleStatistics( const ParticleStatistics &stat )
    : _nboundaries(stat._nboundaries), _end_time(stat._end_time), 
      _end_step(stat._end_step), _end_baddef(stat._end_baddef), _sum_steps(stat._sum_steps)
{
    if( stat._bound_collisions ) {
	_bound_collisions = new uint32_t[_nboundaries];
	_bound_current = new double[_nboundaries];
	for( uint32_t a = 0; a < _nboundaries; a++ ) {
	    _bound_collisions[a] = stat._bound_collisions[a];
	    _bound_current[a] = stat._bound_current[a];
	}
    } else {
	_bound_collisions = NULL;
	_bound_current = NULL;
    }
}


ParticleStatistics::~ParticleStatistics()
{
    if( _bound_collisions )
	delete [] _bound_collisions;
    if( _bound_current )
	delete [] _bound_current;
}


const ParticleStatistics &ParticleStatistics::operator=( const ParticleStatistics &stat )
{
    _end_time = stat._end_time;
    _end_step = stat._end_step;
    _end_baddef = stat._end_baddef;
    _sum_steps = stat._sum_steps;

    if( _nboundaries == stat._nboundaries ) {
	// No need to reallocate, just copy
	for( uint32_t a = 0; a < _nboundaries; a++ ) {
	    _bound_collisions[a] = stat._bound_collisions[a];
	    _bound_current[a] = stat._bound_current[a];
	}
    } else {
	// Reallocate and copy
	if( _bound_collisions )
	    delete [] _bound_collisions;
	if( _bound_current )
	    delete [] _bound_current;

	if( stat._bound_collisions ) {
	    _bound_collisions = new uint32_t[_nboundaries];
	    _bound_current = new double[_nboundaries];
	    for( uint32_t a = 0; a < _nboundaries; a++ ) {
		_bound_collisions[a] = stat._bound_collisions[a];
		_bound_current[a] = stat._bound_current[a];
	    }
	} else {
	    _bound_collisions = NULL;
	    _bound_current = NULL;
	}
    }

    return( *this );
}


void ParticleStatistics::reset( uint32_t nboundaries )
{
    _end_time   = 0;
    _end_step   = 0;
    _end_baddef = 0;
    _sum_steps  = 0;

    if( nboundaries == _nboundaries ) {
	// No need to reallocate, just clear
	for( uint32_t a = 0; a < _nboundaries; a++ ) {
	    _bound_collisions[a] = 0;
	    _bound_current[a] = 0.0;
	}
    } else {
	// Reallocate and clear
	_nboundaries = nboundaries;
	if( _bound_collisions )
	    delete [] _bound_collisions;
	if( _bound_current )
	    delete [] _bound_current;

	if( _nboundaries ) {
	    _bound_collisions = new uint32_t[_nboundaries];
	    _bound_current = new double[_nboundaries];
	    for( uint32_t a = 0; a < _nboundaries; a++ ) {
		_bound_collisions[a] = 0;
		_bound_current[a] = 0.0;
	    }
	} else {
	    _bound_collisions = NULL;
	    _bound_current = NULL;
	}
    }
}


const ParticleStatistics &ParticleStatistics::operator+=( const ParticleStatistics &stat )
{
    if( _nboundaries != stat._nboundaries )
	throw( Error( ERROR_LOCATION, "different number of boundaries" ) );

    _end_time += stat._end_time;
    _end_step += stat._end_step;
    _end_baddef += stat._end_baddef;
    _sum_steps += stat._sum_steps;

    for( uint32_t a = 0; a < _nboundaries; a++ ) {
	_bound_collisions[a] += stat._bound_collisions[a];
	_bound_current[a] += stat._bound_current[a];
    }

    return( *this );
}


void ParticleStatistics::clear( void )
{
    _end_time   = 0;
    _end_step   = 0;
    _end_baddef = 0;
    _sum_steps  = 0;
    for( uint32_t a = 0; a < _nboundaries; a++ ) {
	_bound_collisions[a] = 0;
	_bound_current[a] = 0.0;
    }
}


uint32_t ParticleStatistics::number_of_boundaries( void ) const
{
    return( _nboundaries );
}


uint32_t ParticleStatistics::bound_collisions( uint32_t bound ) const
{
    if( bound == 0 || bound > _nboundaries )
	throw( Error( ERROR_LOCATION, "invalid boundary number" ) );
    return( _bound_collisions[bound-1] );
}


uint32_t ParticleStatistics::bound_collisions( void ) const
{
    uint32_t col = 0;
    for( uint32_t a = 0; a < _nboundaries; a++ )
	col += _bound_collisions[a];
    return( col );
}


double ParticleStatistics::bound_current( uint32_t bound ) const
{
    if( bound == 0 || bound > _nboundaries )
	throw( Error( ERROR_LOCATION, "invalid boundary number" ) );
    return( _bound_current[bound-1] );
}


double ParticleStatistics::bound_current( void ) const
{
    double current = 0;
    for( uint32_t a = 0; a < _nboundaries; a++ )
	current += _bound_current[a];
    return( current );
}


uint32_t ParticleStatistics::end_time( void ) const
{
    return( _end_time );
}


uint32_t ParticleStatistics::end_step( void ) const
{
    return( _end_step );
}


uint32_t ParticleStatistics::end_baddef( void ) const
{
    return( _end_baddef );
}


uint32_t ParticleStatistics::sum_steps( void ) const
{
    return( _sum_steps );
}


void ParticleStatistics::add_bound_collision( uint32_t bound, double IQ )
{
    if( bound > _nboundaries )
	throw( Error( ERROR_LOCATION, "invalid boundary number: (bound = " + to_string(bound)
		      + ") >= (nboundaries = " + to_string(_nboundaries + ")" ) ) );
    _bound_collisions[bound-1]++;
    _bound_current[bound-1] += IQ;
}

