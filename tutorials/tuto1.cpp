#include "kadath.hpp"

using namespace Kadath ;
int main() {

	
	// 1d array of 4 integers
	Array<int> onedint (4) ;

	// 2d array of (8x4) doubles :
	Array<double> twoddoubles (8, 4) ;

	// Array of dimension 5 (of booleans), constructed by a Dim_array :
	Dim_array dimensions(5) ;
	// The size of the various dimensions must be initialized
	for (int i=0 ; i<5 ; i++)
       	dimensions.set(i) = i+1 ;
	Array<bool> multidbool (dimensions) ;

	// 1d array 
	for (int i=0 ; i<4 ; i++)
      		onedint.set(i) = i ;
	cout << onedint << endl ; // Print the whole Array`enter code here`

	// 2d array, same value everywhere
	twoddoubles = 2.3 ;
	// Printing one particular element
	cout << twoddoubles(0, 2) << endl ;

	// Array of dimension 5  :
	Index pos(dimensions) ;
	// By default it is affected to the first element
	// Modify the index in dimension 2
	pos.set(2) = 1 ;
	cout << pos << endl ;
	// pos now points at one particular element of the Array
	multidbool = false ; // Everything is false
	multidbool.set(pos) = true ; // but this one...

	// Sets the index to the first index
	pos.set_start() ;
	// Loop
	do {
      	if (multidbool(pos)) {
                 cout << "True value found at " << endl ;
                 cout << pos << endl ;
		}
     	}
	while (pos.inc()) ;

	return EXIT_SUCCESS ;
}
