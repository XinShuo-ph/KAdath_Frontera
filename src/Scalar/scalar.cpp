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
#include "scalar.hpp"
#include "tensor_impl.hpp"

namespace Kadath {

    void Scalar::save (FILE* fd) const {
        for (int i=0 ; i<get_nbr_domains() ; i++)
            val_zones[i]->save(fd) ;
    }

    Scalar & Scalar::operator= (const Scalar& so)  {
        assert (&espace==&so.espace) ;
        for (int i=0 ; i<ndom ; i++) {
            if (val_zones[i] != nullptr) delete val_zones[i] ;
            val_zones[i] = new Val_domain{*so.val_zones[i]} ;
        }
        return *this;
    }

    Scalar & Scalar::operator= (const Tensor& so)  {
        assert (&espace==&so.espace) ;
        assert (so.valence==0) ;
        for (int i=0 ; i<ndom ; i++) {
            if (val_zones[i] != nullptr) delete val_zones[i] ;
            val_zones[i] = new Val_domain{*so.cmp[0]->val_zones[i]} ;
        }
        return *this;
    }

#ifdef TENSOR_MOVE_SEMANTIC
    Scalar & Scalar::operator=(Tensor && so) noexcept
    {
        assert(so.valence==0);
        this->do_move(std::move(so),false);
        val_zones.swap(so.cmp[0]->val_zones);
        so.cmp[0] = nullptr;
        assert(cmp[0] == this);
        return *this;
    }
#endif

    double Scalar::val_point(const Point& xx, int sens) const {
        assert ((sens==+1) || (sens==-1)) ;

        bool* inside = new bool[ndom] ;
        for (int l=ndom-1 ; l>=0 ; l--)
            inside[l] = get_domain(l)->is_in(xx) ;
        // First domain in which the point is :
        int ld = -1 ;
        if (sens == -1) {
          for (int l=ndom-1 ; l>=0 ; l--)
              if ((ld==-1) && (inside[l])) ld = l ;
        }
        else {
           for (int l=0 ; l<ndom ; l++)
              if ((ld==-1) && (inside[l])) ld = l ;
        }
#ifndef REMOVE_ALL_CHECKS
        if (ld==-1) {
             cout << "Point " << xx << "not found in the computational space..." << endl ;
             abort() ;
        }
#endif
        else {
            if (val_zones[ld]->check_if_zero())
                return 0. ;
            else {
                Point num(get_domain(ld)->absol_to_num(xx)) ;

                coef() ;
                delete [] inside ;
                return val_zones[ld]->base.summation(num, *val_zones[ld]->cf) ;
            }
        }
    }

    double Scalar::val_point_zeronotdef(const Point& xx, int sens) const {
        assert ((sens==+1) || (sens==-1)) ;

        bool* inside = new bool[ndom] ;
        for (int l=ndom-1 ; l>=0 ; l--)
            inside[l] = get_domain(l)->is_in(xx) ;
        // First domain in which the point is :
        int ld = -1 ;
        if (sens == -1) {
          for (int l=ndom-1 ; l>=0 ; l--)
              if ((ld==-1) && (inside[l])) ld = l ;
        }
        else {
           for (int l=0 ; l<ndom ; l++)
              if ((ld==-1) && (inside[l])) ld = l ;
        }

        if (ld==-1) {
            return 0 ;
        }
        else {
            if (val_zones[ld]->check_if_zero())
                return 0. ;
            else {
                Point num(get_domain(ld)->absol_to_num(xx)) ;

                coef() ;
                delete [] inside ;
                return val_zones[ld]->base.summation(num, *val_zones[ld]->cf) ;
            }
        }
    }


    void Scalar::filter_phi (int dom, int ncf) {
        coef() ;
        Index pcf (get_space().get_domain(dom)->get_nbr_coefs()) ;
        int np = get_space().get_domain(dom)->get_nbr_coefs()(2) ;
        do {
            if (pcf(2)+ncf > np-1)
                set_domain(dom).set_coef(pcf) = 0 ;
        }
        while (pcf.inc()) ;
    }


    Scalar Scalar::der_var(int var) const {
        Scalar res (*this, false) ;
        for (int dom=0 ; dom<ndom ; dom++)
            res.set_domain(dom) = val_zones[dom]->der_var(var) ;
        return res ;
    }

    Scalar Scalar::der_abs(int var) const {
        Scalar res (*this, false) ;
        for (int dom=0 ; dom<ndom ; dom++)
            res.set_domain(dom) = val_zones[dom]->der_abs(var) ;

        return res ;
    }

    Scalar Scalar::der_spher(int var) const
    {
       Scalar res(*this, false) ;
       for (int dom(0) ; dom < ndom ; ++dom)
          res.set_domain(dom) = val_zones[dom]->der_spher(var) ;
       return res ;
    }

    Scalar Scalar::der_r() const {
        Scalar res (*this, false) ;
        for (int dom=0 ; dom<ndom ; dom++)
            res.set_domain(dom) = val_zones[dom]->der_r() ;

        return res ;
    }

    double Scalar::integrale() const {
      double res = 0 ;
      for (int dom=0 ; dom<ndom ; dom++)
        res += val_zones[dom]->integrale() ;

        return res ;
    }

    ostream& operator<< (ostream& o, const Scalar& so) {

        o << "Scalar" << endl ;
        for (int l=0 ; l<so.ndom ; l++) {
            o << "Domain : " << l << endl ;
        o << *so.val_zones[l] << endl ;
        }
        return o ;
    }

    double Scalar::integ_volume() const {
      double res = 0 ;
      for (int dom=0 ; dom<ndom ; dom++)
        res += val_zones[dom]->integ_volume() ;

        return res ;
    }

}
