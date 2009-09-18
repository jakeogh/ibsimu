/*! \file scheduler.hpp
 *  \brief Header file for scheduler.hpp
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
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
 * tvkalvas@cc.jyu.fi.
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

#ifndef SCHEDULER_HPP
#define SCHEDULER_HPP 1


#include <pthread.h>
#include <iostream>
#include <vector>
#include <deque>
//#include <sys/time.h>


//pthread_mutex_t cout_mutex = PTHREAD_MUTEX_INITIALIZER;


/*! \brief %Scheduler class for implementing consumer-producer
 *  threading.
 *
 *  %Scheduler uses a manager thread and a given number of working
 *  threads for solving problems. %Scheduler is templated with Solver,
 *  Problem and Error classes. Solver class has to provide an operator
 *  \code void operator()( Problem *p, Scheduler<Solver,Problem,Error> &s ) \endcode 
 *  to solve problem p. Handle to scheduler itself in
 *  this function is provided to solver so that solver can access
 *  add_problems() for hierarchial construction of problems. Problem
 *  is only used by user defined solver class and Scheduler does not
 *  pose any requirements to it.
 *
 *  Error class has to have a default constructor. Scheduler catches
 *  the errors of this type from the working threads and saves caught
 *  errors in a container. If an error is caught, all the working
 *  threads are interrupted and problem solving is finished. The
 *  scheduler does indicate the error state by returning false from
 *  finish(). Error state can also be queried with is_error(). Errors
 *  can be fetched from the internal containers with get_errors().
 *
 *  The %Scheduler can be used for static and dynamic problems. All
 *  threads can add problems using add_problems(). The functions are
 *  mutex protected. %Scheduler provides functions to diagnose the
 *  ending and status of the process.
 */
template <class Solv, class Prob, class Err>
class Scheduler {

    class Consumer {

	/*
	enum consumer_status_e {
	    CONSUMER_CREATED = 0,
	    CONSUMER_RUNNING,
	    CONSUMER_FINISHED
	};
	*/

	//pthread_mutex_t      _mutex;            //!< \brief Mutex for active check
	pthread_t            _thread;
	Solv                *_solver;
	Scheduler           *_scheduler;
	//struct timeval       _t0;
	//std::vector<struct timeval> _t;
    
	void *consumer_main( void ) {
	    Prob *p;
	    //struct timeval t;
	    
	    //pthread_mutex_lock( &_mutex );
	    //_status = CONSUMER_RUNNING;
	    //pthread_mutex_unlock( &_mutex );

	    while( (p = _scheduler->get_next_problem()) ) {
		try {
		    //gettimeofday( &t, NULL );
		    //_t.push_back( t );
		    (*_solver)( p, *_scheduler );
		    //gettimeofday( &t, NULL );
		    //_t.push_back( t );
		} catch( Err e ) {
		    //std::cout << "on_error\n";
		    // Handle error and stop solving
		    _scheduler->on_error( e, p );
		    break;
		};
		_scheduler->put_solved_problem( p );
	    }
      
	    //std::cout << "Exiting consumer\n";
	    //pthread_mutex_lock( &_mutex );
	    //_status = CONSUMER_FINISHED;
	    //pthread_mutex_unlock( &_mutex );
	    return( NULL );
	}
    
    public:

	static void *consumer_entry( void *data ) {
	    Consumer *consumer = (Consumer *)data;
	    return( consumer->consumer_main() );
	}

	Consumer( Solv *solver, Scheduler *scheduler ) : _solver(solver), _scheduler(scheduler) { 

	    //pthread_mutex_init( &_mutex, NULL );
	    //std::cout << "Start\n";
	    //gettimeofday( &_t0, NULL );
	}

	~Consumer() {
	    //pthread_mutex_lock( &cout_mutex );
	    //std::cout << "End\n";
	    //for( size_t a = 0; a < _t.size(); a++ ) {
	    //std::cout << (_t[a].tv_sec-_t0.tv_sec) + 
	    //(_t[a].tv_usec-_t0.tv_usec)/1e6 << "\n";
	    //a++;
	    //std::cout << (_t[a].tv_sec-_t0.tv_sec) + 
	    //(_t[a].tv_usec-_t0.tv_usec)/1e6 << "\n\n\n";
	    //}
	    //pthread_mutex_unlock( &cout_mutex );
	}

	void run( void ) {
	    pthread_create( &_thread, NULL, consumer_entry, (void *)this );
	}

	void join( void ) {
	    //pthread_mutex_lock( &_mutex );
	    //if( _status == CONSUMER_FINISHED ) {
	    //pthread_mutex_unlock( &_mutex );
	    //return;
	    //} else if( _status == CONSUMER_CREATED ) {
	    //
	    //}
	    //pthread_mutex_unlock( &_mutex );
	    pthread_join( _thread, NULL );
	}

    };


    pthread_mutex_t         _mutex;            //!< \brief Mutex for all shared data
    pthread_cond_t          _scheduler_cond;   //!< \brief Wake up scheduler
    pthread_cond_t          _producer_cond;    //!< \brief Wake up producer
    pthread_cond_t          _consumer_cond;    //!< \brief Wake up consumer

    size_t                  _problems_in_c;    //!< \brief Total problems in count
    size_t                  _problems_out_c;   //!< \brief Total problems out count
    size_t                  _problems_err_c;   //!< \brief Total error problems out count
    std::deque<Prob*>       _problems_in;      //!< \brief Problems to be solved
    std::deque<Prob*>       _problems_out;     //!< \brief Problems already solved

    pthread_t               _scheduler_thread; //!< \brief Scheduler main thread
    std::vector<Consumer *> _consumers;        //!< \brief Consumer objects

    bool                    _running;          //!< \brief Are we running
    bool                    _error;            //!< \brief Finish as soon as possible
    bool                    _done;             //!< \brief Exit after current problem
    bool                    _finish;           //!< \brief Finish all problems and exit
    std::vector<Err>        _err;              //!< \brief Error class
    std::vector<Prob *>     _prob;             //!< \brief Problem causing error


    /*! \brief %Error handler
     *
     *  Saves caught error \a e and problem causing the error \a p and
     *  broadcasts other threads signalling about the error condition.
     */
    void on_error( Err &e, Prob *p ) {
	pthread_mutex_lock( &_mutex );
	_err.push_back( e );
	_prob.push_back( p );
	_problems_err_c++;
	_error = true;
	pthread_cond_broadcast( &_scheduler_cond );
	pthread_mutex_unlock( &_mutex );
    }


    Prob *get_next_problem( void ) {
	Prob *ret;
	pthread_mutex_lock( &_mutex );
    
	if( _done || _error ) {
	    pthread_mutex_unlock( &_mutex );
	    return( NULL );
	}

	if( _problems_in.empty() ) {
	    // Signal producer that problems are spent
	    pthread_cond_signal( &_scheduler_cond );
	    while( _problems_in.empty() ) {
		// Wait for new problems
		pthread_cond_wait( &_consumer_cond, &_mutex );
		if( _done || _error ) {
		    pthread_mutex_unlock( &_mutex );
		    return( NULL );
		}
	    }
	}

	// Return next problem
	ret = _problems_in.front();
	_problems_in.pop_front();
	pthread_mutex_unlock( &_mutex );
	return( ret );
    }


    void put_solved_problem( Prob *p ) {
	pthread_mutex_lock( &_mutex );
	_problems_out_c++;
	//std::cout << "put_solved_problem(): " << _problems_out_c << "\n";	
	_problems_out.push_back( p );
	pthread_mutex_unlock( &_mutex );
    }


    void *scheduler_main( void ) {

	// Moved from
	for( size_t a = 0; a < _consumers.size(); a++ )
	    _consumers[a]->run();

	pthread_mutex_lock( &_mutex );

	while( 1 ) {
	    // Wait until all consumers are done with all problems or error occurs
	    while( !(_problems_in.empty() || _done || _error) ) {
		//std::cout << "scheduler_main(): scheduler_cond wait 1\n";
		pthread_cond_wait( &_scheduler_cond, &_mutex );
	    }

	    if( (_finish && _problems_in_c == _problems_out_c+_problems_err_c) || 
		_done || _error )
		break;

	    // Problems temporarily done
	    pthread_cond_wait( &_scheduler_cond, &_mutex );
	    //std::cout << "scheduler_main(): prob_in = " << _problems_in_c
	    //<< " prob_out = " << _problems_out_c << "\n";
	    //std::cout << "scheduler_main(): scheduler_cond wait 2\n";

	    // Signal consumers to wake up
	    pthread_cond_broadcast( &_consumer_cond );
	}

	// Broadcast done
	_done = true;
	pthread_cond_broadcast( &_consumer_cond );
	pthread_mutex_unlock( &_mutex );

	// Join all consumers
	//std::cout << "scheduler_main(): Scheduler waiting in join\n";
	for( size_t a = 0; a < _consumers.size(); a++ )
	    _consumers[a]->join();

	pthread_cond_broadcast( &_producer_cond );
	//std::cout << "scheduler_main(): Exiting scheduler\n";
	_running = false;
	return( NULL );
    }




    static void *scheduler_entry( void *data ) {
	Scheduler *scheduler = (Scheduler *)data;
	return( scheduler->scheduler_main() );
    }


public:


    /*! \brief Constructor
     *
     *  \param s Vector of solvers to be used by \a N problem
     *  consumers, where \a N is the size of the vector.
     */
    Scheduler( std::vector<Solv *> s )
	: _problems_in_c(0), _problems_out_c(0), _problems_err_c(0), _running(false) {

	pthread_mutex_init( &_mutex, NULL );
	pthread_cond_init( &_scheduler_cond, NULL );
	pthread_cond_init( &_consumer_cond, NULL );
	pthread_cond_init( &_producer_cond, NULL );

	// Create consumer threads
	for( size_t a = 0; a < s.size(); a++ )
	    _consumers.push_back( new Consumer( s[a], this ) );
    }


    /*! \brief Destructor
     */
    ~Scheduler() {
	finish();
	pthread_mutex_destroy( &_mutex );
	pthread_cond_destroy( &_scheduler_cond );
	pthread_cond_destroy( &_consumer_cond );
	pthread_cond_destroy( &_producer_cond );

	// Delete consumer threads
	for( size_t a = 0; a < _consumers.size(); a++ )
	    delete _consumers[a];
    }


    /*! \brief Fetch solved problems
     *
     *  \param c Container object where pointers to solved problems
     *  are appended.
     */
    template <class Cont>
    size_t get_solved_problems( Cont &c ) {
	pthread_mutex_lock( &_mutex );
	size_t r = _problems_out.size();
	while( !_problems_out.empty() ) {
	    c.push_back( _problems_out.front() );
	    _problems_out.pop_front();
	}
	pthread_mutex_unlock( &_mutex );
	return( r );
    }

    
    /*! \brief Return true on errors.
     */
    bool is_error( void ) {
	// No mutex needed for one bit read
	return( _error );
    }


    /*! \brief Return true if scheduler is running.
     */
    bool is_running( void ) {
	// No mutex needed for one bit read
	return( _running );
    }


    /*! \brief Fetch errors and corresponding problems.
     *
     *  \param e Container where errors are appended
     *  \param p Container where problems are appended
     *  \return Number of error problems
     */
    template <class Cont1, class Cont2>
    size_t get_errors( Cont1 &e, Cont2 &p ) {
	pthread_mutex_lock( &_mutex );
	size_t r = _err.size();
	for( size_t a = 0; a < _err.size(); a++ ) {
	    e.push_back( _err[a] );
	    p.push_back( _prob[a] );
	}
	_err.clear();
	_prob.clear();
	pthread_mutex_unlock( &_mutex );
	return( r );
    }    

    /*! \brief Run threads.
     *
     *  Returns immediately after creating working threads. Use
     *  finish() or destructor of class to wait for work to be
     *  completed.
     */
    void run( void ) {

	if( _running )
	    return;
	_running = true;
	_error = false;
	_done = false;
	_finish = false;
	_err.clear();
	_prob.clear();
	pthread_create( &_scheduler_thread, NULL, scheduler_entry, (void *)this );
    }


    /*! \brief Add problem to input queue.
     */
    void add_problem( Prob *p ) {

	pthread_mutex_lock( &_mutex );
	_problems_in_c++;
	_problems_in.push_back( p );
	pthread_cond_broadcast( &_scheduler_cond );	
	pthread_mutex_unlock( &_mutex );
    }


    /*! \brief Add multiple problems to input queue.
     */
    void add_problems( std::vector<Prob *> p ) {

	pthread_mutex_lock( &_mutex );
	_problems_in_c += p.size();
	_problems_in.insert( _problems_in.end(), p.begin(), p.end() );
	pthread_cond_broadcast( &_scheduler_cond );	
	pthread_mutex_unlock( &_mutex );
    }


    /*! \brief Wait for all problems to be solved.
     *
     *  \return True if finished already or finish() already called. 
     *  False if any error occured during run or scheduler not running.
     */
    bool finish( void ) {
	if( _finish )
	    return( true );
	if( !_running )
	    return( false );

	pthread_mutex_lock( &_mutex );
	_finish = true;
	//std::cout << "finish(): scheduler_cond broadcast\n";
	pthread_cond_broadcast( &_scheduler_cond );

	//std::cout << "finish(): producer_cond wait\n";
	pthread_cond_wait( &_producer_cond, &_mutex );
	pthread_mutex_unlock( &_mutex );

	if( _error )
	    return( false );
	return( true );
    }


    friend class Consumer;
};



#endif














