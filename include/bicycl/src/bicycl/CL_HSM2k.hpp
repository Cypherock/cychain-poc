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
#ifndef BICYCL_CL_HSM2K_HPP
#define BICYCL_CL_HSM2K_HPP

#include <iostream>

#include "bicycl/gmp_extras.hpp"
#include "bicycl/qfi.hpp"
#include "bicycl/CL_HSM_utils.hpp"
#include "bicycl/seclevel.hpp"

namespace BICYCL
{
  /**
   * Class for the cryptosystem based on the hidden subgroup membership problem,
   * with subgroup a power of 2.
   *
   * Ref: ??
   */
  class CL_HSM2k
  {
    protected:
      /** a product of two odd primes. */
      Mpz N_;

      /** a positive integer */
      size_t k_;

      /** 2^k */
      Mpz M_;

      /** \f$ \ClDeltaK \f$ : the class group of the maximal order.
       * Its discriminant is equal to \f$ -8N \f$.
       */
      ClassGroup Cl_DeltaK_;

      /** \f$ \ClDelta \f$: the class group of the order of conductor
       * \f$2^{k+1}\f$.
       * Its discriminant is equal to \f$ -2^{2k+5} N \f$.
       * It contains the subgroup \f$F\f$.
       */
      ClassGroup Cl_Delta_;

      /** \c true if the compact variant is used, \c false otherwise. */
      bool compact_variant_;

      /** \c true if the large-message variant is used, \c false otherwise. */
      bool large_message_variant_;

      /** The generator of the group \f$G\f$.
       * If the compact variant is not used, the generator is an element of
       * \f$ \ClDelta \f$, else it is an element of \f$ \ClDeltaK \f$.
       */
      QFI h_;

      /** The distance parameter used to produced a almost uniform distribution.
       * Given a bound on the class number of \f$ \ClDeltaK \f$, this bound is
       * multiplied by 2^(distance_-2) to produced a random distribution that
       * is at distance 2^(distance_) of being uniform.
       */
      unsigned int distance_;

      /** Actual bound use to draw random values
       * It is equal to 2^(distance_-2) times Cl_Delta_.class_number_bound_
       */
      Mpz exponent_bound_;

      /** Precomputation data: a positive integer */
      size_t d_;
      size_t e_;
      /** Precomputation data: h_^(2^e_), h_^(2^d_), h_^(d_+e_) */
      QFI h_e_precomp_;
      QFI h_d_precomp_;
      QFI h_de_precomp_;

    public:
      /** Class used to represent a secret key of the cryptosystem */
      using SecretKey = _Utils::CL_HSM_SecretKey<CL_HSM2k>;
      /** Class used to represent a public key of the cryptosystem */
      using PublicKey = _Utils::CL_HSM_PublicKey<CL_HSM2k>;
      /** Class used to represent a cleartext for the cryptosystem */
      using ClearText = _Utils::CL_HSM_ClearText<CL_HSM2k>;
      /** Class used to represent a ciphertext for the cryptosystem */
      using CipherText = _Utils::CL_HSM_CipherText<CL_HSM2k>;

      /**
       * @name Constructors
       *
       * Setup of the cryptosystem
       *
       *@{
       */
      /**
       * Setup of the cryptosystem given @p N and @p k.
       */
      CL_HSM2k (const Mpz &N, size_t k, unsigned int distance,
                bool compact_variant);
      /**
       * Same as above, using default value `false` for @p compact_variant.
       */
      CL_HSM2k (const Mpz &N, size_t k, unsigned int distance);
      /**
       * Same as above, using default value for @p distance.
       */
      CL_HSM2k (const Mpz &N, size_t k, bool compact_variant);
      /**
       * Same as above, using default values.
       */
      CL_HSM2k (const Mpz &N, size_t k);
      /**
       * Copy constructor, only the value of compact variant can be changed.
       */
      CL_HSM2k (const CL_HSM2k &C, bool compact_variant);
      /**
       * Setup of the cryptosystem given the size of @p N and @p k.
       */
      template <class... Ts>
      CL_HSM2k (size_t N_bits, size_t k, RandGen &randgen, Ts... other_args);
      /**
       * Setup of the cryptosystem given the desired security level and @p k.
       *
       * The equivalence between security level and the size of \f$\Delta_K\f$
       * can be found in the class \ref SecLevel.
       */
      template <class... Ts>
      CL_HSM2k (SecLevel seclevel, size_t k, RandGen &randgen,
                                             Ts... other_args);
      /**@}*/

      /**
       * @name Public methods to retrieve the public parameters
       *@{
       */
      /** Return N, a product of two primes. */
      const Mpz & N () const;
      /**
       * Return k; the cardinality of the cyclic subgroup \f$F\f$ is \f$2^k\f$.
       */
      size_t k () const;
      /** Return \f$M=2^{k}\f$. */
      const Mpz & M () const;
      /** Return \f$\Delta_K = -8N\f$. */
      const Mpz & DeltaK () const;
      /** Return \f$\Delta = -2^{2k+2} \Delta_K\f$. */
      const Mpz & Delta () const;
      /**
       * Return \f$\ClDeltaK\f$: the class group of discriminant
       * \f$\Delta_K = -8N\f$.
       */
      const ClassGroup & Cl_DeltaK () const;
      /**
       * Return \f$\ClDelta\f$: the class group of discriminant
       * \f$\Delta = -2^{2k+5}N\f$.
       */
      const ClassGroup & Cl_Delta () const;
      const ClassGroup & Cl_G () const;
      /** Return \f$h\f$, the generator of the cyclic subgroup \f$H\f$ */
      const QFI & h () const;
      /** Return whether the compact variant is used or not */
      bool compact_variant () const;
      /** Return whether the large message variant is used or not */
      bool large_message_variant () const;
      /** Return the bound for secret keys: the bound on the size of \f$H\f$ */
      const Mpz & secretkey_bound () const;
      /** Return the bound for cleartexts: \f$M=2^k\f$ */
      const Mpz & cleartext_bound () const;
      /** Return the bound for random exponents: same as #secretkey_bound */
      const Mpz & encrypt_randomness_bound () const;
      /** Return the distance */
      unsigned int lambda_distance () const;
      /**@}*/

      /**
       * @name Public methods for computation in subgroups
       *@{
       */
      /** Set @p r to \f$h^e\f$, where #h is the generator of \f$H\f$. */
      void power_of_h (QFI &r, const Mpz &e) const;
      /** Return \f$f^m\f$, where `f` is the generator of \f$F\f$. */
      QFI power_of_f (const Mpz &m) const;
      /** Return the discrete logarithm of the form @p fm. */
      Mpz dlog_in_F (const QFI &) const;
      /**
       * Compute \f$\psi_{2^k}(f)\f$ to move @p f from \f$\Delta_K\f$ to
       * \f$\Delta\f$.
       */
      void from_Cl_DeltaK_to_Cl_Delta (QFI &f) const;
      /**@}*/

      /**
       * @name Public methods implementing the cryptographic functionalities
       *@{
       */
      /** Generate a random secret key */
      SecretKey keygen (RandGen &) const;
      /** Compute the public key associated to a secret key */
      PublicKey keygen (const SecretKey &) const;
      /** Encrypt @p m using public key @p pk */
      CipherText encrypt (const PublicKey &, const ClearText &,
                          const Mpz&) const;
      /** Encrypt @p m using public key @p pk and randomness @p r*/
      CipherText encrypt (const PublicKey &, const ClearText &,
                          RandGen &) const;
      /** Decrypt @p c using secret key @p sk*/
      ClearText decrypt (const SecretKey &, const CipherText &) const;
      /** Homomorphically add ciphertexts @p ca and @p cb */
      CipherText add_ciphertexts (const PublicKey &, const CipherText &,
                                  const CipherText &, RandGen &) const;
      /** Homomorphically add ciphertexts @p ca and @p cb using @p r */
      CipherText add_ciphertexts (const PublicKey &, const CipherText &,
                                  const CipherText &, const Mpz &) const;
      /** Add the two cleartexts @p ma and @p mb */
      ClearText add_cleartexts (const ClearText &, const ClearText &) const;
      /** Homomorphically compute @p s times @p c */
      CipherText scal_ciphertexts (const PublicKey &, const CipherText &,
                                   const Mpz &, RandGen &) const;
      /** Homomorphically compute @p s times @p c using @p r*/
      CipherText scal_ciphertexts (const PublicKey &, const CipherText &,
                                   const Mpz &, const Mpz &) const;
      /** Compute @p s times @p m */
      ClearText scal_cleartexts (const ClearText &, const Mpz &) const;
      /**@}*/

      /** Print the public parameters of the cryptosystem */
      friend std::ostream & operator<< (std::ostream &, const CL_HSM2k &);

    protected:
      /* utils for ctor */
      static Mpz random_N (RandGen &randgen, size_t N_nbits);
      static Mpz compute_DeltaK (const Mpz &);
      static Mpz compute_Delta (const Mpz &, size_t);
      /* utils */
      void raise_to_power_M (const ClassGroup &Cl, QFI &f) const;
      void F_kerphi_square (Mpz &, Mpz &, Mpz &, Mpz &) const;
      void F_kerphi_div (Mpz &, const Mpz &, Mpz &, Mpz &, Mpz &) const;
  };

  #include "CL_HSM2k.inl"

} /* BICYCL namespace */

#endif /* BICYCL_CL_HSM2K_HPP */
