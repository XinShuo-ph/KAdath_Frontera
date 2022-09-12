#include "kadath.hpp"

using namespace Kadath ;
int main() {

	// Opening the file in read mode.
	FILE* fin = fopen ("file.dat", "r") ;

	// Two values to be read
	int valint  ;
	double valdouble  ;

	// Read using fread_be
	fread_be (&valint, sizeof(int), 1, fin) ;
	fread_be (&valdouble, sizeof(double), 1, fin) ;

	cout << "Integer read " << valint << endl ;
	cout << "Double read " << valdouble << endl ;

	// Reading Kadath classes
	// Space first, its exact type must be known
	Space_spheric space (fin) ;
	cout << space << endl ;

	// Then the fields
	Scalar field (space, fin) ;
	Vector vecU (space, fin) ;

	// Don't forget to close the file in the end
	fclose(fin) ;

	return EXIT_SUCCESS ;
}

