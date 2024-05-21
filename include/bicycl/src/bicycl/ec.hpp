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
#ifndef BICYCL_EC_HPP
#define BICYCL_EC_HPP

#include <functional>
#include <tuple>

#include "bicycl/seclevel.hpp"
#include "bicycl/gmp_extras.hpp"
#include "bicycl/openssl_wrapper.hpp"

namespace BICYCL
{
  /*****/
  class ECDSA : public OpenSSL::ECGroup
  {
    public:
      using PublicKey = OpenSSL::ECPoint;
      using Message = std::vector<unsigned char>;

      /*** SecretKey ***/
      class SecretKey
      {
        public:
          explicit SecretKey (const ECDSA &C);

          const OpenSSL::BN & d () const;
          const OpenSSL::ECPoint & Q () const;

        private:
          OpenSSL::BN d_;
          OpenSSL::ECPoint Q_;
      };

      /*** Signature ***/
      class Signature
      {
        public:
          /* constructors */
          Signature (const ECDSA &C, const SecretKey &sk, const Message &m);

          bool verify (const ECDSA &C, const PublicKey &Q,
                                       const Message &m) const;

        private:
          OpenSSL::BN r_, s_;
      };

      /* constructors */
      explicit ECDSA (SecLevel seclevel);

      /* crypto protocol */
      SecretKey keygen () const;
      PublicKey keygen (const SecretKey &sk) const;
      Signature sign (const SecretKey &sk, const Message &m) const;
      bool verif (const Signature &s, const PublicKey &Q,
                                      const Message &m) const;

      /* utils */
      static Message random_message ();

    protected:
      OpenSSL::BN hash_message (const Message &m) const;

    private:
      mutable OpenSSL::HashAlgo H_;
  }; /* ECDSA */

  /*****/
  class ECNIZKProof
  {
    public:
      using SecretValue = OpenSSL::BN;
      using PublicValue = OpenSSL::ECPoint;

      ECNIZKProof (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                   const SecretValue &s, const PublicValue &Q);
      ECNIZKProof (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                   const SecretValue &s);
      ECNIZKProof (const OpenSSL::ECGroup &E, const ECNIZKProof &p);

      bool verify (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                   const PublicValue &Q) const;

    protected:
      static OpenSSL::BN hash_for_challenge (OpenSSL::HashAlgo &H,
                                             const OpenSSL::ECGroup &E,
                                             const OpenSSL::ECPoint &R,
                                             const OpenSSL::ECPoint &Q);
      static OpenSSL::ECPoint compute_Q_from_secret (const OpenSSL::ECGroup &E,
                                                     const SecretValue &s);

    private:
      OpenSSL::ECPoint R_;
      OpenSSL::BN z_;
  }; /* ECNIZKProof */

  /*****/
  class ECNIZKAoK
  {
    public:
      using SecretValue = OpenSSL::BN;
      using PublicValue = OpenSSL::ECPoint;

      ECNIZKAoK (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                 const OpenSSL::ECPoint &R, const SecretValue &x,
                 const SecretValue &y, const SecretValue &rho,
                 const PublicValue &V, const PublicValue &A);
      ECNIZKAoK (const OpenSSL::ECGroup &E, const ECNIZKAoK &p);

      bool verify (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                   const OpenSSL::ECPoint &R, const PublicValue &V,
                   const PublicValue &A) const;

    protected:
      static OpenSSL::BN hash_for_challenge (OpenSSL::HashAlgo &Hash,
                                             const OpenSSL::ECGroup &E,
                                             const OpenSSL::ECPoint &R,
                                             const OpenSSL::ECPoint &V,
                                             const OpenSSL::ECPoint &A,
                                             const OpenSSL::ECPoint &H);

    private:
      OpenSSL::ECPoint H_;
      OpenSSL::BN t1_;
      OpenSSL::BN t2_;
  }; /* ECNIZKAoK */


  #include "ec.inl"

} /* BICYCL namespace */

#endif /* BICYCL_EC_HPP */
