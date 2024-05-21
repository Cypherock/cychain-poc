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
#ifndef BICYCL_GMP_EXTRAS_HPP
#define BICYCL_GMP_EXTRAS_HPP

#include <gmp.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <vector>

namespace BICYCL
{
  namespace _gmp_impl
  {
#if 0
    static inline void
    mpz_invert_2exp_scratch (mpz_ptr r, mpz_srcptr a, mp_bitcnt_t n, mpz_ptr t)
    {
      mpz_set_ui (r, 1);
      for (size_t i = 1; i < n; i <<= 1)
      {
        mpz_set_ui (t, 2);
        mpz_submul (t, r, a);
        mpz_mul (t, r, t);
        mpz_fdiv_r_2exp (r, t, i << 1);
      }
    }

    /*
     * Input: a, r and n
     *
     * Output: set r to a modulo 2^n with r in ]-2^(n-1), 2^(n-1)]
     *
     * Assumption:
     *    - n >= 1
     *
     */
    static inline void
    mpz_mod_with_centered_remainder_2exp (mpz_ptr r, mpz_srcptr a,
                                          mp_bitcnt_t n)
    {
        mpz_fdiv_r_2exp (r, a, n);

        /* substract 2^n if needed */
        if (mpz_tstbit (r, n-1))
        {
          mpn_com (r->_mp_d, r->_mp_d, r->_mp_size);
          mpz_add_ui (r, r, 1);
          mpz_fdiv_r_2exp (r, r, n);
          mpz_neg_inplace (r);
        }
    }
#endif
  } /* _gmp_impl namespace */

  /* forward declaration, needed to declare it friend of Mpz */
  class RandGen;

  class Mpz
  {
    protected:
      mpz_t mpz_;

    public:
      /* constructors */
      Mpz ();
      Mpz (const Mpz &);
      Mpz (Mpz &&);
      explicit Mpz (unsigned long);
      explicit Mpz (long);
      explicit Mpz (const std::string &);
      explicit Mpz (mpf_srcptr);
      explicit Mpz (const std::vector<unsigned char> &);

      /* destructor */
      ~Mpz ();

      /* assignment */
      Mpz & operator= (const Mpz &);
      Mpz & operator= (Mpz &&);
      Mpz & operator= (unsigned long);
      Mpz & operator= (long);
      Mpz & operator= (const std::string &);
      Mpz & operator= (mpf_srcptr); /* only needed once in ClassGroup */
      Mpz & operator= (const std::vector<unsigned char> &);

      /* comparison */
      bool operator== (const Mpz &) const;
      bool operator!= (const Mpz &) const;
      bool operator< (const Mpz &) const;
      bool operator> (const Mpz &) const;
      bool operator<= (const Mpz &) const;
      bool operator>= (const Mpz &) const;
      bool operator== (unsigned long) const;
      bool operator!= (unsigned long) const;
      bool operator< (unsigned long) const;
      bool operator> (unsigned long) const;
      bool operator<= (unsigned long) const;
      bool operator>= (unsigned long) const;
      bool operator== (long) const;
      bool operator!= (long) const;
      bool operator< (long) const;
      bool operator> (long) const;
      bool operator<= (long) const;
      bool operator>= (long) const;

      /* conversion */
      explicit operator mpz_srcptr() const;
      explicit operator unsigned long() const;
      explicit operator long() const;

      /* getters */
      /* */
      size_t nbits () const;
      size_t ndigits () const;
      size_t nlimbs () const;
      int sgn () const;

      /* tests */
      bool is_zero () const;
      bool is_odd () const;
      bool is_even () const;
      bool is_one () const;
      bool is_prime (int reps=BICYCL_GMP_PRIMALITY_TESTS_ITERATION) const;
      bool is_divisible_by (const Mpz &d) const;

      /* misc */
      void neg ();
      mp_limb_t extract_bits (size_t, size_t) const;
      int tstbit (size_t) const;
      void setbit (size_t);
      unsigned long mod4 () const;
      unsigned long mod8 () const;
      size_t val2 () const;
      void nextprime ();
      int legendre (const Mpz &) const;
      int jacobi (const Mpz &) const;
      int kronecker (const Mpz &) const;

      /* */
      static void swap (Mpz &, Mpz &);

      /* arithmetic */
      static void abs (Mpz &, const Mpz &);
      static int cmpabs (const Mpz &, const Mpz &);

      static void add (Mpz &, const Mpz &, const Mpz &);
      static void add (Mpz &, const Mpz &, unsigned long);

      static char *get_str(const Mpz &op1, size_t base);
      static void set_str(Mpz &op1, const char *str, size_t base);

      static void from_bytes(Mpz &op1, uint8_t *data, size_t size);
      static void to_bytes(const Mpz &op1, uint8_t *data, size_t size);

      static void sub (Mpz &, const Mpz &, const Mpz &);
      static void sub (Mpz &, const Mpz &, unsigned long op2);

      static void mul (Mpz &, const Mpz &, const Mpz &);
      static void mul (Mpz &, const Mpz &, unsigned long);
      static void mulby2k (Mpz &, const Mpz &, mp_bitcnt_t k);
      static void mulby2k (Mpz &, unsigned long, mp_bitcnt_t k);
      static void mulby2 (Mpz &, const Mpz &);
      static void mulby4 (Mpz &, const Mpz &);

      static void addmul (Mpz &, const Mpz &, const Mpz &);
      static void submul (Mpz &, const Mpz &, const Mpz &);

      static void divby2k (Mpz &, const Mpz &, mp_bitcnt_t k);
      static void divby2 (Mpz &, const Mpz &);
      static void divby4 (Mpz &, const Mpz &);
      static void divexact (Mpz &, const Mpz &, const Mpz &);
      static void divexact (Mpz &, const Mpz &, unsigned long);
      static void cdiv_qr (Mpz &, Mpz &, const Mpz &, const Mpz &);
      static void fdiv_qr (Mpz &, Mpz &, const Mpz &, const Mpz &);

      static void mod (Mpz &, const Mpz &, const Mpz &);
      static void mod2k (Mpz &, const Mpz &, mp_bitcnt_t);
      static void mod2k_centered (Mpz &, const Mpz &, mp_bitcnt_t);
      // TODO inverse -> invert
      static void mod_inverse (Mpz &, const Mpz &, const Mpz &);
      static void mod_inverse_2k (Mpz &, const Mpz &, mp_bitcnt_t k);
      static void mod_inverse_2k (Mpz &, const Mpz &, mp_bitcnt_t k, Mpz &);
      static void pow_mod (Mpz &, const Mpz &, const Mpz &, const Mpz &);
      static void pow_mod (Mpz &, const Mpz &, const Mpz &, const Mpz &,
                           size_t, const Mpz &, const Mpz &, const Mpz &);

      static void gcd (Mpz &, const Mpz &, const Mpz &);
      static void gcdext (Mpz &, Mpz &, Mpz &, const Mpz &, const Mpz &);
      static void lcm (Mpz &, const Mpz &, const Mpz &);
      static void sqrt (Mpz &, const Mpz &);
      static void root4th (Mpz &, const Mpz &);
      static void sqrt_mod_prime (Mpz &, const Mpz &, const Mpz &);

      static size_t remove (Mpz &, const Mpz &, const Mpz &);

      /* */
      static void CRT (Mpz &, const Mpz &, const Mpz &, const Mpz &,
                       const Mpz &);

      /* */
      static void ceil_abslog_square (Mpz &, const Mpz &);

      /* */
      static void partial_euclid (Mpz &, Mpz &, Mpz &, Mpz &, Mpz &, Mpz &,
                                  mp_size_t, Mpz &, Mpz &);
      static void partial_euclid (Mpz &, Mpz &, Mpz &u10, Mpz &u11, Mpz &,
                                  Mpz &, mp_size_t);

      /* I/O */
      friend std::ostream & operator<< (std::ostream &, const Mpz &);
      friend std::istream & operator>> (std::istream &, Mpz &);

      /* exception */
      class ModInverseException;

    protected:
      static int richcmp (const Mpz &, const Mpz &);
      static int richcmp (const Mpz &, unsigned long);
      static int richcmp (const Mpz &, long);

    friend RandGen;
  }; /* Mpz class */


  class RandGen
  {
    protected:
      gmp_randstate_t rand_;

    public:
      RandGen ();
      RandGen (const RandGen &);
      explicit RandGen (const Mpz &);
      ~RandGen ();

      void set_seed (const Mpz &);

      Mpz random_mpz (const Mpz &);
      Mpz random_mpz_2exp (mp_bitcnt_t);

      unsigned char random_uchar ();
      std::vector<unsigned char> random_bytes (size_t);
      unsigned long random_ui (unsigned long);
      unsigned long random_ui_2exp (mp_bitcnt_t);

      Mpz random_negative_discriminant (mp_bitcnt_t);

      bool random_bool ();

      Mpz random_prime (mp_bitcnt_t);
  }; /* RandGen */


  /* */
  class JSF : protected std::vector<uint8_t>
  {
    public:
      using std::vector<uint8_t>::size;

      JSF (const Mpz &n0, const Mpz &n1);

      uint8_t operator[] (size_t i) const;

    protected:
      void set (size_t i, int d0, int d1);
  };

  #include "gmp_extras.inl"

} /* BICYCL namespace */

#endif /* BICYCL_GMP_EXTRAS_HPP */
