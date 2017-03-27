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
namespace Kadath {
			//-------------------------//
			//	int version 	   //
			//-------------------------//
			

int fwrite_be(const int* aa, int size, int nb, FILE* fich) {

	assert(size == 4) ;
	
	// Determines whether the default storage is big endian
	//  or large endians
	
	int itest = 1 ;
	bool little_endian =  ( *( reinterpret_cast<char*>(&itest) ) == 1)  ;
	
	if (little_endian) {

		int size_tot = 4 * nb ;

		char* bytes_big = new char[size_tot] ;
		char* pbig =  bytes_big ;
		const char* plit =  reinterpret_cast<const char*>(aa) ;
		
		for (int j=0; j< nb; j++) {
		
			for (int i=0; i<4; i++) {
				pbig[i] = plit[3-i] ;
			}
		
			plit += 4 ; 	// next item
			pbig += 4 ;
			
		}
		
		
		int nx =  static_cast<int>(fwrite(bytes_big, 1, size_tot, fich)) / 4 ;

		delete [] bytes_big ; 
		
		return nx ; 

	}
	else {  // Big endian case: nothing to do:
	
		return static_cast<int>(fwrite(aa, size, nb, fich)) ;
	}
		
}


			//-------------------------//
			//	double version 	   //
			//-------------------------//
			

int fwrite_be(const double* aa, int size, int nb, FILE* fich) {

	assert(size == 8) ;
	
	// Determines whether the default storage is big endian
	//  or large endians
	
	int itest = 1 ;
	bool little_endian = ( *( reinterpret_cast<char*>(&itest) ) == 1) ;
	
	if (little_endian) {

		int size_tot = 8 * nb ;

		char* bytes_big = new char[size_tot] ;
		char* pbig =  bytes_big ;
		const char* plit = reinterpret_cast<const char*>(aa) ;
		
		for (int j=0; j< nb; j++) {
		
			for (int i=0; i<8; i++) {
				pbig[i] = plit[7-i] ;
			}
		
			plit += 8 ; 	// next item
			pbig += 8 ;
			
		}
		
		int nx = static_cast<int>(fwrite(bytes_big, 1, size_tot, fich)) / 8 ;		
		delete [] bytes_big ; 
		
		return nx ; 

	}
	else {  // Big endian case: nothing to do:
	
		return static_cast<int>(fwrite(aa, size, nb, fich)) ;
	}
		
}}
