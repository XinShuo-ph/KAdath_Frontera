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

/*!\mainpage KADATH --- Reference manual
 *
 */

/**
 * \defgroup util Utilities.
 *
 * The class template \c Array implements mult-dimensional arrays (mainly of integers and doubles).
 * The size of the various dimensions are stored in the class \c Dim_array and the position in a given array by the class 
 *\c Index.
 *
 * The class \c Matrice is used for various matrices operations, mainly to solve linear systems.
 * 
 * \c Param is the storage class.
 *
 * All these classes are independent of the numerical method 
 * (spectral method).
 *
 */


/**
 * \defgroup domain Description of the physical space.
 *
 * Each domain class are derived from the abstract class \c Domain. The domains contain the informations 
 * on the number of collocation points and of coefficients and, more importantly define the coordinates use by 
 * implementing the relations between the numerical coordinates and the physical ones.
 * 
 * The class \c Space is an ensemble of \c Domains describing the computational space.
 *
 */

/**
 * \defgroup spectral Spectral representation.
 *
 * The fields are given, in each \c Domain, by a \c Val_domain objetcs. The basis are stored inside the 
\c Base_spectral class.
 */

/** 
 * \defgroup fields  Physical fields.
 *
 * Here are regrouped all the objects used to describe fields, without making reference to systems of equations.
 *
 */

/** 
 * \defgroup metric Metric handling.
 * 
 * They differ from the usual fields because they give access to many more methods (covariant derivatives, Ricci etc...)
 * They are also more closely linked to systems of equations.
 */



/** 
 * \defgroup systems Equations management
 * 
 * Are regrouped here the objects needed to solve systems of equations, like the \c Term_eq, the various operators, types of equations...
 */

