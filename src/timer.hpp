/*! \file timer.hpp
 *  \brief Header file defining timer.
 */


#ifndef TIMER_HPP
#define TIMER_HPP 1


#include <iostream>
#include <time.h>
#include <sys/time.h>


/*! \brief Class for measuring code runtime in cpu time and realtime.
 */
class Timer {

    clock_t         _cpu1, _cpu2;
    struct timeval  _time1, _time2;

public:

    /*! \brief Construct and start timer.
     */
    Timer();

    /*! \brief Destruct timer.
     */
    ~Timer() {}

    /*! \brief Start timer.
     */
    void start( void );

    /*! \brief Stop timer.
     */
    void stop( void );
    
    /*! \brief Get real time elapsed in seconds.
     */
    double get_real_time( void ) const;

    /*! \brief Get cpu time elapsed in seconds.
     */
    double get_cpu_time( void ) const;

    /*! \brief Print report of elapsed time to stream.
     */    
    friend std::ostream &operator<<( std::ostream &os, const Timer &t );
};

#endif













