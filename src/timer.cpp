#include "timer.hpp"


Timer::Timer()
{
    _cpu1 = clock();
    gettimeofday( &_time1, NULL );
}


void Timer::start( void )
{
    _cpu1 = clock();
    gettimeofday( &_time1, NULL );
}


void Timer::stop( void )
{
    _cpu2 = clock();
    gettimeofday( &_time2, NULL );
}


double Timer::get_real_time( void ) const
{
    return( (double)(_time2.tv_sec - 
                     _time1.tv_sec) + 
            ((double)(_time2.tv_usec - 
                      _time1.tv_usec))/1.0e6 );
}


double Timer::get_cpu_time( void ) const
{
    return( (_cpu2-_cpu1)/(double)CLOCKS_PER_SEC );
}


std::ostream &operator<<( std::ostream &os, const Timer &t )
{
    os << t.get_cpu_time() << " s (" << t.get_real_time() << " s realtime)";
    return( os );
}













