/*
 * BICYCL Implements CryptographY in CLass groups
 * Copyright (C) 2022  Cyril Bouvier <cyril.bouvier@lirmm.fr>
 *                     Guilhem Castagnos <guilhem.castagnos@math.u-bordeaux.fr>
 *                     Laurent Imbert <laurent.imbert@lirmm.fr>
 *                     Fabien Laguillaumie <fabien.laguillaumie@lirmm.fr>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef BICYCL_QFI_HPP
#define BICYCL_QFI_HPP

#include <iostream>
#include <cstddef>
#include <stdexcept>

#include "bicycl/openssl_wrapper.hpp"
#include "bicycl/gmp_extras.hpp"

namespace BICYCL
{
  /* forward declaration */
  class ClassGroup; /* needed to declare it friend of QFI */
  class QFICompressedRepresentation;

  /**
   * Binary forms of imaginary quadratic fields.
   *
   * The form is always primitive (gcd (a,b,c) = 1), in reduced form and with a
   * positive (It implies that c is also positive).
   * All public methods assumed that the input forms satisfy these assumptions
   * and ensure that the output forms also satisfy these assumptions.
   *
   * Protected methods do not necessarily respect these assumptions.
   */
  class QFI
  {
    protected:
      Mpz a_, b_, c_;

    public:
      /** default ctor, set the form to (1,1,1) */
      QFI ();
      /** Use with care */
      QFI (const Mpz &a, const Mpz &b, const Mpz &c, bool bypass_check=false);
      QFI (const QFI & other) = default;
      QFI (const QFICompressedRepresentation &compressed_form, const Mpz &disc);
      QFI (QFI && other) = default;

      /* assignment operators */
      QFI & operator= (const QFI & other) = default;
      QFI & operator= (QFI && other) = default;

      /* comparison operators */
      bool operator== (const QFI & other) const;

      /** Getter for the coefficient \f$a\f$ of the form */
      const Mpz & a () const;
      /** Getter for the coefficient \f$b\f$ of the form */
      const Mpz & b () const;
      /** Getter for the coefficient \f$c\f$ of the form */
      const Mpz & c () const;

      /** Return the discriminant of the form */
      Mpz discriminant () const;

      /** Return \c true if the form the neutral form, \c false otherwise */
      bool is_one () const;

      /* */
      void neg ();

      /* */
      Mpz eval (const Mpz &, const Mpz &) const;

      /* */
      QFICompressedRepresentation compressed_repr () const;

      /* */
      void lift (const Mpz &);
      void lift_2exp (unsigned int);
      void to_maximal_order (const Mpz &, const Mpz &, bool);
      void to_maximal_order_2exp (unsigned int, const Mpz &, bool);

      Mpz kernel_representative (const Mpz &, const Mpz &) const;
      Mpz kernel_representative_2exp (size_t, const Mpz &) const;

      /* I/O */
      friend std::ostream & operator<< (std::ostream &, const QFI &);

      /* */
      struct OpsAuxVars
      {
        Mpz Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, m00, m01, m10, m11, t0, t1,
            m, s,  F, u, v, x, y, H, l, q, by;
      };

    protected:
      void set_c_from_disc (const Mpz &disc);

      /* methods for reduction */
      void normalize ();
      void normalize (Mpz &tmp0, Mpz &tmp1);
      void rho ();
      void rho (Mpz &tmp0, Mpz &tmp1);
      void reduction ();
      void reduction (Mpz &tmp0, Mpz &tmp1);

      /* */
      void prime_to (const Mpz &l);
      void prime_to_2exp ();

      /* compositions and exponentiations */
      static void nucomp (QFI &, const QFI &, const QFI &, const Mpz &,
                          bool negf2);
      static void nucomp (QFI &, const QFI &, const QFI &, const Mpz &,
                          bool negf2, OpsAuxVars &);
      static void nudupl (QFI &, const QFI &, const Mpz &);
      static void nudupl (QFI &, const QFI &, const Mpz &, OpsAuxVars &);
      static void nupow (QFI &, const QFI &, const Mpz &, const Mpz &) ;
      static void nupow (QFI &, const QFI &, const Mpz &, const QFI &,
                          const Mpz &, const Mpz &);
      static void nupow (QFI &, const QFI &, const Mpz &, size_t, size_t,
                          const QFI &, const QFI &, const QFI &, const Mpz &);

      /* friend class */
      friend ClassGroup;
  };

  /*
   * Ref: https://eprint.iacr.org/2020/196.pdf (algo 2 and 3)
   */
  class QFICompressedRepresentation
  {
    public:
      const Mpz ap;
      const Mpz g;
      const Mpz tp;
      const Mpz b0;
      const bool is_neg;

      QFICompressedRepresentation () = delete;
      QFICompressedRepresentation (const Mpz &, const Mpz &, const Mpz &,
                                   const Mpz &, bool);

      /* getters */
      size_t nbits () const;

      /* I/O */
      friend std::ostream & operator<< (std::ostream &,
                                        const QFICompressedRepresentation &);
  };

  /** Class groups of binary forms of imaginary quadratic fields.
   */
  class ClassGroup
  {
    protected:
      Mpz disc_;
      Mpz default_nucomp_bound_;
      mutable Mpz class_number_bound_;

    public:
      explicit ClassGroup (const Mpz &discriminant);

      /* getters */
      const Mpz & discriminant () const;
      const Mpz & default_nucomp_bound () const;

      /* create qfi */
      QFI one () const;
      QFI primeform (const Mpz &) const;
      template <size_t ngens=10, unsigned long coeff_bound=(1U << 16)>
      QFI random (RandGen &randgen) const;

      /* */
      const Mpz & class_number_bound () const;

      /* */
      void nucomp (QFI &, const QFI &, const QFI &) const;
      void nucomp (QFI &, const QFI &, const QFI &, const Mpz &) const;
      void nucompinv (QFI &, const QFI &, const QFI &) const;
      void nudupl (QFI &, const QFI &) const;
      void nudupl (QFI &, const QFI &, size_t) const;
      void nupow (QFI &, const QFI &, const Mpz &) const;
      void nupow (QFI &, const QFI &, const Mpz &, const QFI &,
                  const Mpz &) const;
      void nupow (QFI &, const QFI &, const Mpz &, size_t, size_t, const QFI &,
                  const QFI &, const QFI &) const;
  };

  #include "qfi.inl"

} /* BICYCL namespace */

#endif /* BICYCL_QFI_HPP */
