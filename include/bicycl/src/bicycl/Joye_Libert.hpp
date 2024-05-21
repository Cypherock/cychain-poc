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
#ifndef BICYCL_JOYE_LIBERT_HPP
#define BICYCL_JOYE_LIBERT_HPP

#include <iostream>

#ifdef BICYCL_WITH_PTHREADS
#include <thread>
#endif

#include "bicycl/gmp_extras.hpp"
#include "bicycl/seclevel.hpp"

namespace BICYCL
{
  class JoyeLibert
  {
    protected:
      size_t N_nbits_;
      size_t k_;
      Mpz twotothek_;

    public:
      /*** SecretKey ***/
      class SecretKey
      {
        protected:
          Mpz r_;
          Mpz p_;
          Mpz q_;
          Mpz y_;
          Mpz D_;

        public:
          /* constructors */
          SecretKey (const JoyeLibert &, RandGen &);

          /* getters */
          const Mpz & r() const;
          const Mpz & p() const;
          const Mpz & q() const;
          const Mpz & y() const;
          const Mpz & D() const;
      };

      /*** PublicKey ***/
      class PublicKey
      {
        protected:
          Mpz N_;
          Mpz y_;

        public:
          /* constructors */
          explicit PublicKey (const SecretKey &);

          /* getters */
          const Mpz & N() const;
          const Mpz & y() const;
      };

      /* Forward declaration */
      class CipherText;

      /*** ClearText ***/
      class ClearText : public Mpz
      {
        public:
          /* constructors */
          ClearText (const JoyeLibert &, const Mpz &);
          ClearText (const JoyeLibert &, RandGen &);
          ClearText (const JoyeLibert &, const SecretKey &, const CipherText &);
      };

      /*** CipherText ***/
      class CipherText : public Mpz
      {
        public:
          /* constructors */
          CipherText (const JoyeLibert &, const PublicKey &, const ClearText &,
                      const Mpz &);
      };

      /* ctor */
      JoyeLibert (size_t N_nbits, size_t k);
      JoyeLibert (SecLevel seclevel, size_t k);

      /* getters */
      size_t k() const;
      size_t N_nbits() const;
      const Mpz & cleartext_bound () const;
      const Mpz & encrypt_randomness_bound (const PublicKey &) const;

      /* crypto protocol */
      SecretKey keygen (RandGen &) const;
      PublicKey keygen (const SecretKey &) const;
      CipherText encrypt (const PublicKey &, const ClearText &,
                          const Mpz&) const;
      CipherText encrypt (const PublicKey &, const ClearText &,
                          RandGen &) const;
      ClearText decrypt (const SecretKey &, const CipherText &) const;
  };

  #include "Joye_Libert.inl"

} /* BICYCL namespace */

#endif /* BICYCL_JOYE_LIBERT_HPP */
