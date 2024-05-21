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
#ifndef BICYCL_PAILLIER_HPP
#define BICYCL_PAILLIER_HPP

#include <iostream>

#ifdef BICYCL_WITH_PTHREADS
#include <thread>
#endif

#include "bicycl/gmp_extras.hpp"
#include "bicycl/seclevel.hpp"

namespace BICYCL
{
  /**
   * Paillier cryptosystem
   */
  class Paillier
  {
    protected:
      size_t N_nbits_;

    public:
      /*** SecretKey ***/
      class SecretKey
      {
        protected:
          Mpz p_;
          Mpz q_;

          Mpz lambda_;
          Mpz lambda_inv_;

        public:
          /* constructors */
          SecretKey (const Paillier &, RandGen &);

          /* getters */
          const Mpz & p() const;
          const Mpz & q() const;
          const Mpz & lambda() const;
          const Mpz & lambda_inv() const;
      };

      /*** PublicKey ***/
      class PublicKey
      {
        protected:
          Mpz N_;
          Mpz NN_;

        public:
          /* constructors */
          explicit PublicKey (const SecretKey &);

          /* getters */
          const Mpz & N () const;
          const Mpz & N_square () const;
      };

      /* Forward declaration */
      class CipherText;

      /*** ClearText ***/
      class ClearText : public Mpz
      {
        public:
          /* constructors */
          ClearText (const Paillier &, const PublicKey &, const Mpz &);
          ClearText (const Paillier &, const PublicKey &, RandGen &);
          ClearText (const PublicKey &, const SecretKey &, const CipherText &);
      };

      /*** CipherText ***/
      class CipherText : public Mpz
      {
        public:
          /* constructors */
          CipherText (const PublicKey &, const ClearText &, const Mpz &);
      };

      /* ctor */
      explicit Paillier (size_t N_nbits);
      explicit Paillier (SecLevel seclevel);

      /* getters */
      size_t N_nbits () const;
      const Mpz & cleartext_bound (const PublicKey &) const;
      const Mpz & encrypt_randomness_bound (const PublicKey &) const;

      /* crypto protocol */
      SecretKey keygen (RandGen &) const;
      PublicKey keygen (const SecretKey &) const;
      CipherText encrypt (const PublicKey &, const ClearText &,
                          const Mpz&) const;
      CipherText encrypt (const PublicKey &, const ClearText &,
                          RandGen &) const;
      ClearText decrypt (const PublicKey &, const SecretKey &,
                         const CipherText &) const;
  };

  /**
   * CamenischShoup cryptosystem
   *
   * Based on description of Section 3.3 of https://eprint.iacr.org/2020/385.pdf
   */
  class CamenischShoup
  {
    protected:
      Mpz N_;
      Mpz NN_;
      Mpz Nover4_;
      Mpz gen_;

      /** Precomputation data: a positive integer */
      size_t e_;
      /** Precomputation data: gen_^(2^e_), gen_^(2^2e_), gen_^(2^3e_) */
      Mpz gen_e_precomp_;
      Mpz gen_2e_precomp_;
      Mpz gen_3e_precomp_;

    public:
      using SecretKey = _Utils::CL_HSM_SecretKey<CamenischShoup>;
      class PublicKey
      {
        protected:
          Mpz pk_;

          /** Precomputation data: a positive integer */
          size_t e_;
          /** Precomputation data: pk_^(2^e_), pk_^(2^2e_), pk_^(2^3e_) */
          Mpz pk_e_precomp_;
          Mpz pk_2e_precomp_;
          Mpz pk_3e_precomp_;


        public:
          /* constructors */
          PublicKey (const CamenischShoup &, const SecretKey &);

          void exponentiation (const CamenischShoup &C, Mpz &r, const Mpz &n)
                                                                          const;

          /* getter */
          const Mpz & elt() const;
      };
      using ClearText = _Utils::CL_HSM_ClearText<CamenischShoup>;
      class CipherText
      {
        protected:
          /** Two Mpzs */
          Mpz c1_, c2_;

        public:
          CipherText (const CamenischShoup &, const PublicKey &,
                      const ClearText &, const Mpz &);

          /* getters */
          const Mpz & c1() const;
          const Mpz & c2() const;
      };

      /* constructors */
      CamenischShoup (const Mpz &N, const Mpz &r);
      CamenischShoup (const Mpz &N, RandGen &randgen);
      CamenischShoup (size_t N_nbits, RandGen &randgen);
      CamenischShoup (SecLevel seclevel, RandGen &randgen);

      /* getters */
      const Mpz & gen() const;
      const Mpz & N() const;
      const Mpz & Nsquare() const;

      const Mpz & secretkey_bound () const;
      const Mpz & cleartext_bound () const;
      const Mpz & encrypt_randomness_bound () const;

      /* */
      void power_of_gen (Mpz &r, const Mpz &e) const;

      /* crypto protocol */
      SecretKey keygen (RandGen &) const;
      PublicKey keygen (const SecretKey &) const;
      CipherText encrypt (const PublicKey &, const ClearText &,
                          const Mpz&) const;
      CipherText encrypt (const PublicKey &, const ClearText &,
                          RandGen &) const;
      ClearText decrypt (const SecretKey &, const CipherText &) const;

    protected:
      /* utils for ctor */
      static Mpz random_N (RandGen &randgen, size_t N_nbits);
      static Mpz random_r (RandGen &randgen, const Mpz &N);
  };

  #include "Paillier.inl"

} /* BICYCL namespace */

#endif /* BICYCL_PAILLIER_HPP */
