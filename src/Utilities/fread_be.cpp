/*
    Copyright 2017 Philippe Grandclement

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

// C headers
#include <stdio.h>
#include <assert.h>
#include <stdexcept>
namespace Kadath {
			//-------------------------//
			//	int version 	   //
			//-------------------------//
			
int fread_be(int* aa, int size, int nb, FILE* fich) {
	if (fich == NULL) {
    throw std::runtime_error("Filepointer is null.\n");
  }

	assert(size == 4) ;
	
	// Determines whether the default storage is big endian
	//  or large endians
	
	int itest = 1 ;
	bool little_endian = ( *( reinterpret_cast<char*>(&itest) ) == 1) ;
	if (little_endian) {

		int size_tot = 4 * nb ;

		char* bytes_big = new char[size_tot] ;
		
		int nr = static_cast<int>(fread(bytes_big, 1, size_tot, fich)) ;
		
		char* pbig =  bytes_big ;
		char* plit = reinterpret_cast<char*>( aa ) ;
		
		for (int j=0; j< nb; j++) {
		
			for (int i=0; i<4; i++) {
				plit[i] = pbig[3-i] ;
			}
		
			plit += 4 ; 	// next item
			pbig += 4 ;
			
		}
		
		delete [] bytes_big ; 

		return nr / 4 ;		

	}
	else {  // Big endian case: nothing to do:
	
		return static_cast<int>(fread(aa, size, nb, fich)) ;
	}
		
}

			//-------------------------//
			//	double version 	   //
			//-------------------------//
			
int fread_be(double* aa, int size, int nb, FILE* fich) {
	if (fich == NULL) {
    throw std::runtime_error("Filepointer is null.\n");
  }

	assert(size == 8) ;
	
	// Determines whether the default storage is big endian
	//  or large endians
	
	int itest = 1 ;
	bool little_endian = ( *( reinterpret_cast<char*>(&itest) ) == 1) ;

	if (little_endian) {

		int size_tot = 8 * nb ;

		char* bytes_big = new char[size_tot] ;

		int nr = static_cast<int>(fread(bytes_big, 1, size_tot, fich)) ;
		
		char* pbig =  bytes_big ;
		char* plit = reinterpret_cast<char*>( aa ) ;
		
		for (int j=0; j< nb; j++) {
		
			for (int i=0; i<8; i++) {
				plit[i] = pbig[7-i] ;
			}
		
			plit += 8 ; 	// next item
			pbig += 8 ;
			
		}

		delete [] bytes_big ;

		return nr / 8 ;		
		
	}
	else {  // Big endian case: nothing to do:
	
		return static_cast<int>(fread(aa, size, nb, fich)) ;
	}
		
}
}
