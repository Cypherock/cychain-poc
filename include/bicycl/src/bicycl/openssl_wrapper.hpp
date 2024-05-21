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
#ifndef BICYCL_OPENSSL_WRAPPER_HPP
#define BICYCL_OPENSSL_WRAPPER_HPP

#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <openssl/bn.h>
#include <openssl/ec.h>
#include <openssl/evp.h>
#include <openssl/obj_mac.h> /* for NID_* */
#include <openssl/rand.h>

#include "bicycl/gmp_extras.hpp"
#include "bicycl/seclevel.hpp"

namespace BICYCL
{
  namespace OpenSSL
  {
    /*****/
    void random_bytes (unsigned char *buf, int num);

    /*****/
    class HashAlgo
    {
      public:
        using Digest = std::vector<unsigned char>;

        static const int SHAKE128 = NID_shake128;
        static const int SHAKE256 = NID_shake256;
        static const int SHA3_224 = NID_sha3_224;
        static const int SHA3_256 = NID_sha3_256;
        static const int SHA3_384 = NID_sha3_384;
        static const int SHA3_512 = NID_sha3_512;

        /* constructors */
        explicit HashAlgo (SecLevel seclevel); /* Use SHA3 with desired
                                                * security level
                                                */
        explicit HashAlgo (int nid);
        HashAlgo (const HashAlgo &H);
        HashAlgo (HashAlgo &&H);

        /* destructor */
        ~HashAlgo ();

        /* assignment */
        HashAlgo & operator= (const HashAlgo &H);
        HashAlgo & operator= (HashAlgo &&H);

        /* */
        template <typename... Args>
        Digest operator() (const Args & ...args);

        template <typename T>
        void hash (const T &t);

        /* */
        int digest_nbytes () const;
        int digest_nbits () const;

      private:
        static EVP_MD_CTX * new_ctx_ ();

        void hash_update ();
        template <typename First, typename... Args>
        void hash_update (const First & first, const Args & ...args);
        void hash_bytes (const void *ptr, size_t n);

        const EVP_MD *md_;
        EVP_MD_CTX *mdctx_;
    };

    /*****/
    class ECGroup; /* forward declaration */

    /*****/
    class BN
    {
      friend ECGroup;

      public:
        /* constructors */
        BN ();
        BN (const BN &other);
        BN (BN &&other);
        explicit BN (const std::vector<unsigned char> &bytes);
        explicit BN (const Mpz &v);
        explicit BN (unsigned long v);

        /* destructor */
        ~BN ();

        /* assignment */
        BN & operator= (const BN &other);
        BN & operator= (BN &&other);
        BN & operator= (const Mpz &other);

        /* */
        explicit operator Mpz () const;

        /* interface with unsigned long */
        BN & operator= (unsigned long other);
        BN & operator*= (unsigned long other);

        /* comparisons */
        bool operator== (const BN &other) const;
        bool operator!= (const BN &other) const;
        bool is_zero () const;

        /* */
        int num_bytes () const;
        void neg ();
        static void add (BN &r, const BN &a, const BN &b);

        /* */
        friend std::ostream & operator<< (std::ostream &o, const BN &v);
        friend void random_BN_2exp (const BN &r, int nbits);
        friend void HashAlgo::hash<> (const BN &v);

      private:
        static std::vector<unsigned char> BIGNUM_abs_to_bytes (const BIGNUM *v);

        BIGNUM *bn_;
    }; /* BN */

    /*****/
    class ECPoint
    {
      friend ECGroup;

      public:
        /* constructors */
        explicit ECPoint (const ECGroup &E);
        ECPoint (const ECGroup &E, const ECPoint &Q);
        ECPoint (const ECGroup &E, const BN & m);
        ECPoint (const ECPoint &) = delete;
        ECPoint (ECPoint &&);

        /* assignment */
        ECPoint & operator= (const ECPoint &);
        ECPoint & operator= (ECPoint &&);

        /* destructor */
        ~ECPoint ();

      private:
        EC_POINT *P_;
    }; /* ECPoint */

    /****/
    using ECPointGroupCRefPair = std::tuple<const ECPoint &, const ECGroup &>;

    /****/
    class ECGroup
    {
      /* Constructors of ECPoint need to access ec_group_ to create EC_POINT* */
      friend ECPoint::ECPoint (const ECGroup &);
      friend ECPoint::ECPoint (const ECGroup &, const ECPoint &);

      public:
        static const int P224 = NID_secp224r1;
        static const int P256 = NID_X9_62_prime256v1;
        static const int P384 = NID_secp384r1;
        static const int P521 = NID_secp521r1;

        /* constructors */
        explicit ECGroup (SecLevel seclevel);
        ECGroup (const ECGroup &G) = delete;
        ECGroup (ECGroup &&G);

        /* destructor */
        ~ECGroup ();

        /* assignment */
        ECGroup & operator= (const ECGroup &G) = delete;
        ECGroup & operator= (ECGroup &&G);

        /* getters */
        const Mpz & order () const;
        BN a () const;
        BN b () const;
        BN p () const;

        /* */
        bool is_on_curve (const ECPoint &P) const;
        bool is_in_group (const ECPoint &P) const;
        bool is_at_infinity (const ECPoint &P) const;

        /* elliptic operations */
        void coords_of_point (BN &x, BN &y, const ECPoint &P) const;
        void x_coord_of_point (BN &x, const ECPoint &P) const;
        bool ec_point_eq (const ECPoint &P, const ECPoint &Q) const;
        void ec_add (ECPoint &R, const ECPoint &P,
                                 const ECPoint &Q) const;
        void ec_neg (ECPoint &R) const;
        void scal_mul_gen (ECPoint &R, const BN &n) const;
        void scal_mul (ECPoint &R, const BN &n,
                                   const ECPoint &P) const;
        void scal_mul (ECPoint &R, const BN &m, const BN &n,
                                                    const ECPoint &P) const;

        /* arithmetic operations modulo the group order */
        void mod_order (BN &r, const BN &a) const;
        void add_mod_order (BN &r, const BN &a, const BN &b) const;
        void sub_mod_order (BN &r, const BN &a, const BN &b) const;
        void mul_mod_order (BN &r, const BN &a, const BN &b) const;
        void mul_by_word_mod_order (BN &r, BN_ULONG w) const;
        void inverse_mod_order (BN &r, const BN &a) const;
        BN random_mod_order () const;

        /* utils */
        bool has_correct_prime_order (const ECPoint &P) const;
        bool is_positive_less_than_order (const BN &v) const;

        friend std::ostream & operator<< (std::ostream &o, const ECGroup &E);
        friend std::ostream & operator<< (std::ostream &o,
                                          const ECPointGroupCRefPair &pt);
        friend void HashAlgo::hash<> (const ECGroup &E);
        friend void HashAlgo::hash<> (const ECPointGroupCRefPair &pt);

      private:
        const BIGNUM * get_order() const;
        void coords_of_point (BN &x, BN &y, const EC_POINT *P) const;
        std::vector<unsigned char> EC_POINT_to_bytes (const EC_POINT *P) const;
        std::vector<unsigned char> ECPoint_to_bytes (const ECPoint &P) const;

        EC_GROUP *ec_group_;
        Mpz order_;
        BN_CTX *ctx_;
    }; /* ECGroup */

    #include "openssl_wrapper.inl"

  }; /* namespace OpenSSL */

} /* namespace BICYCL */

#endif /* BICYCL_OPENSSL_WRAPPER_HPP */
