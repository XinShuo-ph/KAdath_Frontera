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

#include "headcpp.hpp"
#include "space.hpp"
#include "scalar.hpp"

#include <cpgplot.h>  
namespace Kadath {
void des_equipot(float* uutab, int nx, int ny, float xmin, float xmax, 
		 float ymin, float ymax, int ncour, const char* nomx, const char* nomy, 
		 const char* title, const char* device = 0x0, int newgraph = 3, 
		 int nxpage = 1, int nypage = 1) ;

void des_equipot(float* uutab, int nx, int ny, float xmin, float xmax, 
		 float ymin, float ymax, int ncour, const char* nomx, const char* nomy, 
		 const char* title, const char* device, int newgraph, int nxpage, 
		 int nypage) {
		 
    // Search for the extremal values of the field : 
    // -------------------------------------------

    float uumin = uutab[0] ;
    float uumax = uutab[0] ;
    for (int i=1; i<nx*ny; i++) {
	uumin = (uutab[i] < uumin) ? uutab[i] : uumin ;
	uumax = (uutab[i] > uumax) ? uutab[i] : uumax ;	
    }

    cout << "  " << title << " : min, max : " << uumin << "   " << uumax 
         << endl ; 

    // Values of equipotentials
    // -------------------------
 
    float* isopot = new float [ncour] ;
    float hh = float(uumax-uumin) / float(ncour) ; 
    for (int i=0; i<ncour; i++) {
	isopot[i] = uumin + hh * float(i) ;
    }
    
    // Array defining the grid for pgcont_
    // -----------------------------------
    float hx = (xmax - xmin)/float(nx-1) ; 
    float hy = (ymax - ymin)/float(ny-1) ; 

    float tr[6] ;
    tr[0] = xmin - hx ;
    tr[1] = hx ;
    tr[2] = 0 ;
    tr[3] = ymin - hy ; 
    tr[4] = 0 ;
    tr[5] = hy ;
     
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_equipot: problem in opening PGPLOT display !" << endl ;
	}

    }

    // Taille des caracteres:
    float size = float(1.3) ;
    cpgsch(size) ;
    
    // Epaisseur des traits:
    int lepais = 1 ; 
    cpgslw(lepais) ;
    
    // Fonte axes: caracteres romains:
    cpgscf(2) ;
    
    // Cadre de la figure
    cpgenv(xmin, xmax, ymin, ymax, 1, 0 ) ; 
    cpglab(nomx,nomy,title) ;

    // On n'effectue le dessin que si la dynamique est suffisante
    
    float dynamique = float(fabs(uumax - uumin)) ; 

    if (dynamique > 1.e-14) {
    
	cpgcont(uutab, nx, ny, 1, nx, 1, ny, isopot, ncour, tr) ;
	
    }
    
    // Closing the graphical output
    // ----------------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }
    
    
    delete [] isopot ; 

}


void des_coupe (const Scalar& uu, const Point& x0, 
		int num_un, double var_un_min, double var_un_max, 
		int num_deux, double var_deux_min, double var_deux_max, 
		const char* title, const char* axis_one, const char* axis_two, int ncour, int n_un, int n_deux) {

    assert ((num_un>0) && (num_un<=uu.get_ndim())) ;
    assert ((num_deux>0) && (num_deux<=uu.get_ndim())) ;
    assert (num_un != num_deux) ;

    // Plot of isocontours
    // -------------------
    float* uutab = new float[n_un*n_deux] ; 
    
    double h_un = (var_un_max - var_un_min) / double(n_un-1) ; 
    double h_deux = (var_deux_max - var_deux_min) / double(n_deux-1) ; 

    Point points(x0) ;

    for (int j=0; j<n_deux; j++) {
	
	double var_deux = var_deux_min + h_deux * j ; 
	
	for (int i=0; i<n_un; i++) {
    
	    double var_un = var_un_min + h_un * i ; 
	    
	    points.set(num_un) = var_un ;
            points.set(num_deux) = var_deux ;
	    
	    uutab[n_deux*j+i] = float(uu.val_point(points)) ;
	}
    }

    const char* nomy = (axis_two==0x0) ?  "" : axis_two ; 
    const char* nomx = (axis_one==0x0) ? "" : axis_one ; 
    
    const char* titi = (title==0x0) ? "" : title ;

    char* device = 0x0 ; 
    int newgraph = 3 ; 
    
    des_equipot(uutab, n_un, n_deux, float(var_un_min), float(var_un_max), float(var_deux_min), float(var_deux_max), 
		ncour, nomx, nomy,
		titi, device, newgraph) ;
		
    delete [] uutab ; 
} 



void des_coupe_zeronotdef (const Scalar& uu, const Point& x0, 
		int num_un, double var_un_min, double var_un_max, 
		int num_deux, double var_deux_min, double var_deux_max, 
		const char* title, const char* axis_one, const char* axis_two, int ncour, int n_un, int n_deux) {

    assert ((num_un>0) && (num_un<=uu.get_ndim())) ;
    assert ((num_deux>0) && (num_deux<=uu.get_ndim())) ;
    assert (num_un != num_deux) ;

    // Plot of isocontours
    // -------------------
    float* uutab = new float[n_un*n_deux] ; 
    
    double h_un = (var_un_max - var_un_min) / double(n_un-1) ; 
    double h_deux = (var_deux_max - var_deux_min) / double(n_deux-1) ; 

    Point points(x0) ;

    for (int j=0; j<n_deux; j++) {
	
	double var_deux = var_deux_min + h_deux * j ; 
	
	for (int i=0; i<n_un; i++) {
    
	    double var_un = var_un_min + h_un * i ; 
	    
	    points.set(num_un) = var_un ;
            points.set(num_deux) = var_deux ;
	    
	    uutab[n_deux*j+i] = float(uu.val_point_zeronotdef(points)) ;
	}
    }

    const char* nomy = (axis_two==0x0) ?  "" : axis_two ; 
    const char* nomx = (axis_one==0x0) ? "" : axis_one ; 
    
    const char* titi = (title==0x0) ? "" : title ;

    char* device = 0x0 ; 
    int newgraph = 3 ; 
    
    des_equipot(uutab, n_un, n_deux, float(var_un_min), float(var_un_max), float(var_deux_min), float(var_deux_max), 
		ncour, nomx, nomy,
		titi, device, newgraph) ;
		
    delete [] uutab ; 
} 
void des_sphere (const Scalar& uu, const Point& x0, double rad, const char* title, int ncour, int n_theta, int n_phi) {

  // Check dim :
  if (uu.get_space().get_ndim()!=3) {
      cerr << "des_sphere only defined for 3-dimensional spaces" << endl ;
      abort() ;
  }
  
    // Plot of isocontours
    // -------------------
    float* uutab = new float[n_theta*n_phi] ; 
    
    double h_phi = 2*M_PI / double(n_phi-1) ; 
    double h_theta = M_PI / double(n_theta-1) ; 
    
    double theta = 0 ;
    double phi ;
    Point MM(3) ;

    for (int j=0 ; j<n_theta ; j++) {
      phi = 0 ;
      for (int i=0 ; i<n_phi ; i++) {


	      MM.set(1) = rad*sin(theta)*cos(phi) + x0(1) ;
	      MM.set(2) = rad*sin(theta)*sin(phi) + x0(2) ;
	      MM.set(3) = rad*cos(theta) + x0(3) ;
	      uutab[j*n_theta+i] = float(uu.val_point(MM)) ;
	      phi += h_phi ;
    }
    theta += h_theta ;
  }

   const char* nomy =   "theta" ; 
   const char* nomx =  "phi" ; 
    
    const char* titi = (title==0x0) ? "" : title ;

    char* device = 0x0 ; 
    int newgraph = 3 ; 
    
    des_equipot(uutab, n_phi, n_theta, 0, float(2*M_PI), 0, float(M_PI), ncour, nomx, nomy, titi, device, newgraph) ;
		
    delete [] uutab ; 
}}
