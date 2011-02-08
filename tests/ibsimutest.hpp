#include "error.hpp"


/*! \brief %Error class for test errors/failures.
 */
class ErrorTest : public Error {

public:

    /*! \brief Constructor for test error.
     */
    ErrorTest( const ErrorLocation &loc, const std::string &str );
};


/* Main function for tests.
 */
int main( int argc, char **argv );


/* Test function.
 */
void test( int argc, char **argv );
