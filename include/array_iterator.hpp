//
// Created by sauliac on 20/03/2020.
//

#ifndef __ARRAY_ITERATOR_HPP_
#define __ARRAY_ITERATOR_HPP_

#include <vector>

#include "headcpp.hpp"
#include "dim_array.hpp"

namespace Kadath {

    /**
     * Version of the \c Index class optimized for incremental access to Array components.
     */
    class Array_iterator : public Memory_mapped
    {
    protected:
        /**
        * Sizes of the associated \c Array.
        * When used with a \c Tensor, it is the dimension, for each tensorial index.
        */
        Dim_array sizes ;
        //! The incremental index step with respect to each dimensions.
        Memory_mapped_array<int> steps;
        int position; ///< Corresponding value for 1D indexing .

    public:
        explicit Array_iterator(int ndim, int tndim) : sizes{ndim}, steps(ndim + 1), position{0}
        {steps.front()=1; for(int i=0;i<ndim;i++) {steps[i+1] = steps[i]*tndim; sizes.set(i) = tndim;}}
        /**Standard constructor.
        * All the positions are set to zero.
        * @param dim [input] Sizes in each dimensions.
        **/
        explicit Array_iterator (const Dim_array& dim) : sizes{dim}, steps(dim.get_ndim() + 1), position{0}
        {steps.front()=1; for(int i=0;i<dim.get_ndim();i++) steps[i+1] = steps[i]*sizes(i);}

        //Copy constructor and assignment operator, move constructor and assignment as well as the destructor are the
        // default ones generated by the compiler.
        /**
        * Returns the number of dimensions.
        */
        int get_ndim() const {return sizes.get_ndim() ;} ;
        //! Accessor for the \c value member.
        int get_position() const {return position;}

        //Array_iterator & set_position

        bool operator==(Array_iterator const & xx) const {return (position == xx.position && get_ndim() == xx.get_ndim());}

        Array_iterator & set_value(int _position) { position = _position; return *this;}
        /**
         * Set the position of the iterator at the same as the one of the source.
         * @param _so iterator to copy the position from.
         * @return a reference toward the current object.
         */
        Array_iterator & set(Array_iterator const &_so) {position = _so.position; return *this;}
        bool check_value();
        void set_start() { position = 0;}
        /**
        * Returns all the dimensions
        */
        Dim_array const& get_sizes() const {return sizes ;} ;
        /**
        * Increments the position of the \c Array_iterator.
        * If one reaches the last point of a dimension, then the next one is increased.
        * @param increm [input] value of the increment.
        * @param var [input] dimension to be incremented.
        * @return \c false if the result is outside the \c Array and \c true otherwise.
        */
        bool inc (int increm, int var=0) { position += increm * steps[var]; return check_value();}
        /**
         * Optimized unit increment with respect to the passed index number.
         * @param var index number
         * @return \c true while the last component index is not passed.
         */
        bool inc1(int var)  { position += steps[var]; return check_value();}
        /**
         * Optimized unit increment.
         * @return \c true while the last component index is not passed.
         */
        bool inc() {position++; return check_value(); }
        template <class> friend class Array ;
    };

    inline bool Array_iterator::check_value() {
        if(position >= steps.back()) {
            position = steps.back() - 1;
            return false;
        }
        else return true;
    }

}

#endif //__ARRAY_ITERATOR_HPP_
