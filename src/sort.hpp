/*! \brief Sort index-value pairs in ascending index order. 
 *  
 *  Sort uses insertion sort algorithm, which is rather fast with
 *  small arrays, which should be the case with most sparse matrix
 *  applications.
 */
inline void sort_iv( int *ind, double *val, int start, int end )
{
    int keyind, k, l;
    double keyval;
    
    for( k = start+1; k < end; k++ ) {
	// Take a new key
	keyind = ind[k];
	keyval = val[k];
	// Move smaller values up one position
	for( l = k; l > start && ind[l-1] > keyind; l-- ) {
	    ind[l] = ind[l-1];
	    val[l] = val[l-1];
	}
	// Insert key
	ind[l] = keyind;
	val[l] = keyval;
    }
}














