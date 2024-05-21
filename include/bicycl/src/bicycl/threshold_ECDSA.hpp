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
#ifndef BICYCL_THRESHOLD_ECDSA_HPP
#define BICYCL_THRESHOLD_ECDSA_HPP

#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "bicycl/gmp_extras.hpp"
#include "bicycl/openssl_wrapper.hpp"
#include "bicycl/ec.hpp"
#include "bicycl/CL_HSMqk.hpp"

namespace BICYCL
{
  /****/
  class thresholdECDSA
  {
    public:
      using Message = std::vector<unsigned char>;
      using PublicKey = OpenSSL::ECPoint;
      using Commitment = OpenSSL::HashAlgo::Digest;
      using CommitmentSecret = std::vector<unsigned char>;
      using ParticipantsList = std::vector<unsigned int>;
      template <class T>
      using ParticipantsMap = std::unordered_map<unsigned int, T>;

      /* exception */
      class ProtocolAbortError;

      /* */
      class KeygenPart1
      {
        public:
          /* ctor */
          KeygenPart1 (const thresholdECDSA &C, unsigned int n, unsigned int t,
                       unsigned int i);

          /* getters */
          unsigned int n () const;
          unsigned int t () const;
          unsigned int i () const;
          const OpenSSL::ECPoint & Q_part () const;
          const OpenSSL::BN & u_part () const;
          const Commitment & commitment () const;
          const CommitmentSecret & commitment_secret () const;
          const OpenSSL::ECPoint & V (size_t k) const;
          const OpenSSL::BN & sigma (size_t j) const;

        private:
          unsigned int i_;
          OpenSSL::BN u_;
          OpenSSL::ECPoint Q_; /* Q_i = [u_i] P */
          Commitment c_;
          CommitmentSecret cs_;
          std::vector<OpenSSL::BN> a_; /* t scalars a_i,k */
          std::vector<OpenSSL::ECPoint> V_; /* V_i,k = [a_i,k] P */
          std::vector<OpenSSL::BN> sigma_; /* n scalars corresponding to
                                            * evaluations of the polynomial
                                            * u_i + sum_{k=1}^{t}{a_i,k X^k}
                                            */
      }; /* KeygenPart1 */

      /* */
      class KeygenPart2
      {
        public:
          /* ctor */
          KeygenPart2 (const thresholdECDSA &C, const KeygenPart1 &data1,
                       RandGen &randgen, const std::vector<Commitment> &Co,
                       const std::vector<OpenSSL::ECPoint> &Q,
                       const std::vector<CommitmentSecret> &CoSec,
                       const std::vector<std::vector<OpenSSL::ECPoint>> &V,
                       const std::vector<OpenSSL::BN> &Sigma);

          /* getters */
          const OpenSSL::ECPoint & Q () const;
          const OpenSSL::BN & x () const;
          const ECNIZKProof & zk_proof () const;
          const CL_HSMqk::SecretKey & CL_secret_key () const;
          const CL_HSMqk::PublicKey & CL_public_key () const;

        private:
          OpenSSL::ECPoint Q_;
          OpenSSL::BN x_;
          ECNIZKProof zk_proof_;
          CL_HSMqk::SecretKey sk_;
          CL_HSMqk::PublicKey pk_;
      }; /* KeygenPart2 */

      /* */
      class SecretKey
      {
        public:
          /* ctor */
          SecretKey (const thresholdECDSA &C, unsigned int i,
                     const KeygenPart1 &data1, const KeygenPart2 &data2,
                     const std::vector<std::vector<OpenSSL::ECPoint>> &V,
                     const std::vector<ECNIZKProof> &ZK,
                     const std::vector<CL_HSMqk::PublicKey> &PK);

          /* getters */
          const PublicKey & public_key () const;
          const CL_HSMqk::SecretKey & CL_secret_key () const;
          const CL_HSMqk::PublicKey & CL_public_key (unsigned int) const;
          const OpenSSL::ECPoint & X (unsigned int) const;
          const OpenSSL::BN & x_part () const;

        private:
          CL_HSMqk::SecretKey sk_;
          OpenSSL::BN x_;
          std::vector<CL_HSMqk::PublicKey> PK_;
          std::vector<OpenSSL::ECPoint> X_;
          OpenSSL::ECPoint Q_;
      }; /* Secret Key */

      /* */
      class SignPart1
      {
        public:
          /* ctor */
          SignPart1 (const thresholdECDSA &C, unsigned int i,
                     const ParticipantsList &S,
                     const thresholdECDSA::SecretKey &sk, RandGen &randgen);

          /* getters */
          const ParticipantsList & S () const;
          const OpenSSL::BN & gamma () const;
          const OpenSSL::ECPoint & Gamma () const;
          const ECNIZKProof & zk_gamma () const;
          const Commitment & commitment () const;
          const CommitmentSecret & commitment_secret () const;
          unsigned int i () const;
          const OpenSSL::BN & omega() const;
          const Mpz & k_part () const;
          const CL_HSMqk::CipherText & ciphertext () const;
          const CL_HSMqk_ZKAoKProof & zk_encrypt_proof () const;

        private:
          unsigned int i_;
          ParticipantsList S_;
          OpenSSL::BN omega_;
          OpenSSL::BN gamma_;
          OpenSSL::ECPoint Gamma_;
          ECNIZKProof zk_gamma_;
          Commitment co_;
          CommitmentSecret cos_;
          CL_HSMqk::ClearText k_;
          Mpz r_;
          CL_HSMqk::CipherText c_;
          CL_HSMqk_ZKAoKProof zk_encrypt_;
      }; /* SignPart1 */

      /* */
      class SignPart2
      {
        public:
          SignPart2 (const thresholdECDSA &C, const SignPart1 &data,
                     const thresholdECDSA::SecretKey &sk,
                     const ParticipantsMap<Commitment> & commitment_map,
                     const ParticipantsMap<CL_HSMqk::CipherText> &c_map,
                     const ParticipantsMap<CL_HSMqk_ZKAoKProof> &zk_map,
                     RandGen &randgen);

          /* getters */
          const Commitment & commitment (unsigned int j) const;
          const OpenSSL::BN & nu (unsigned int j) const;
          const OpenSSL::ECPoint & B (unsigned int j) const;
          const CL_HSMqk::ClearText & beta (unsigned int j) const;
          const CL_HSMqk::CipherText & c_kg (unsigned int j) const;
          const CL_HSMqk::CipherText & c_kw (unsigned int j) const;

        private:
          ParticipantsMap<Commitment> commitment_map_;
          ParticipantsMap<OpenSSL::BN> nu_map_;
          ParticipantsMap<OpenSSL::ECPoint> B_map_;
          ParticipantsMap<CL_HSMqk::ClearText> beta_map_;
          ParticipantsMap<CL_HSMqk::CipherText> c_kg_map_;
          ParticipantsMap<CL_HSMqk::CipherText> c_kw_map_;
      }; /* SignPart2 */

      /* */
      class SignPart3
      {
        public:
          SignPart3 (const thresholdECDSA &C, const SignPart1 &data1,
                     const SignPart2 &data2,
                     const thresholdECDSA::SecretKey &sk,
                     const ParticipantsMap<CL_HSMqk::CipherText> &c_kg_map,
                     const ParticipantsMap<CL_HSMqk::CipherText> &c_kw_map,
                     const ParticipantsMap<OpenSSL::ECPoint> &B_map);

          const OpenSSL::BN & delta_part () const;
          const OpenSSL::BN & sigma_part () const;

        private:
          OpenSSL::BN delta_;
          OpenSSL::BN sigma_;
      }; /* SignPart3 */

      /* */
      class SignPart4
      {
        public:
          SignPart4 (const thresholdECDSA &C, const SignPart1 &data1,
                     const ParticipantsMap<OpenSSL::BN> &delta_map);

          const OpenSSL::BN & delta () const;

        private:
          OpenSSL::BN delta_;
      }; /* SignPart4 */

      /* */
      class SignPart5
      {
        public:
          SignPart5 (const thresholdECDSA &C, const SignPart1 &data1,
                     const SignPart2 &data2, const SignPart3 &data3,
                     const SignPart4 &data4, const Message &m,
                     const ParticipantsMap<OpenSSL::ECPoint> &Gamma_map,
                     const ParticipantsMap<CommitmentSecret> &CoSec_map,
                     const ParticipantsMap<ECNIZKProof> &zk_proof_map);

          /* getters */
          const Message & m () const;
          const OpenSSL::ECPoint & R () const;
          const OpenSSL::BN & r () const;
          const OpenSSL::BN & z () const;
          const OpenSSL::BN & s_part () const;
          const OpenSSL::ECPoint & V_part () const;
          const OpenSSL::BN & ell () const;
          const OpenSSL::BN & rho () const;
          const OpenSSL::ECPoint & A_part () const;
          const Commitment & commitment () const;
          const CommitmentSecret & commitment_secret () const;

        private:
          const Message m_;
          OpenSSL::ECPoint R_;
          OpenSSL::BN r_;
          OpenSSL::BN z_;
          OpenSSL::BN s_;
          OpenSSL::BN ell_;
          OpenSSL::BN rho_;
          OpenSSL::ECPoint V_;
          OpenSSL::ECPoint A_;
          Commitment c_;
          CommitmentSecret cs_;
      }; /* SignPart5 */

      /* */
      class SignPart6
      {
        public:
          SignPart6 (const thresholdECDSA &C, const SignPart5 &data5,
                     const ParticipantsMap<Commitment> &Co_map);

          /* getters */
          const Commitment & commitment (unsigned int j) const;
          const ECNIZKAoK & aok () const;

        private:
          ParticipantsMap<Commitment> commitment_map_;
          ECNIZKAoK zk_aok_;
      }; /* SignPart6 */

      /* */
      class SignPart7
      {
        public:
          SignPart7 (const thresholdECDSA &C, const SignPart1 &data1,
                     const SignPart5 &data5, const SignPart6 &data6,
                     const thresholdECDSA::SecretKey &sk,
                     const ParticipantsMap<OpenSSL::ECPoint> &V_map,
                     const ParticipantsMap<OpenSSL::ECPoint> &A_map,
                     const ParticipantsMap<CommitmentSecret> &CoSec_map,
                     const ParticipantsMap<ECNIZKAoK> &zk_aok_map);

          /* getters */
          const Commitment & commitment () const;
          const CommitmentSecret & commitment_secret () const;
          const OpenSSL::ECPoint & U_part () const;
          const OpenSSL::ECPoint & T_part () const;

        private:
          OpenSSL::ECPoint U_;
          OpenSSL::ECPoint T_;
          Commitment c_;
          CommitmentSecret cs_;
      }; /* SignPart7 */

      /* */
      class SignPart8
      {
        public:
          SignPart8 (const thresholdECDSA &C, const SignPart1 &data1,
                     const SignPart7 &data7,
                     const ParticipantsMap<Commitment> &Co_map,
                     const ParticipantsMap<OpenSSL::ECPoint> &U_map,
                     const ParticipantsMap<OpenSSL::ECPoint> &T_map,
                     const ParticipantsMap<CommitmentSecret> &CoSec_map);
      }; /* SignPart8 */

      /* */
      class Signature
      {
        public:
          Signature (const thresholdECDSA &C, const SignPart1 &data1,
                     const SignPart5 &data5, const SecretKey &sk,
                     const ParticipantsMap<OpenSSL::BN> &s_map);

          bool verify (const thresholdECDSA &C, const PublicKey &Q,
                                                const Message &m) const;

          bool operator== (const Signature &other) const;
          bool operator!= (const Signature &other) const;

        private:
          OpenSSL::BN r_;
          OpenSSL::BN s_;
      }; /* Signature */

      /* constructors */
      thresholdECDSA (SecLevel seclevel, RandGen &randgen);

      /* getters */
      const OpenSSL::ECGroup & ec_group () const;

      /* utils */
      std::tuple<Commitment, CommitmentSecret> commit (
                                              const OpenSSL::ECPoint &Q) const;
      std::tuple<Commitment, CommitmentSecret> commit (
                                              const OpenSSL::ECPoint &Q1,
                                              const OpenSSL::ECPoint &Q2) const;
      bool open (const Commitment &c, const OpenSSL::ECPoint &Q,
                                      const CommitmentSecret &r) const;
      bool open (const Commitment &c, const OpenSSL::ECPoint &Q1,
                                      const OpenSSL::ECPoint &Q2,
                                      const CommitmentSecret &r) const;

      /* crypto */
      bool verify (const Signature &s, const PublicKey &Q,
                                       const Message &m) const;

      /* utils */
      OpenSSL::BN lagrange_at_zero (const ParticipantsList &S,
                                    unsigned int i) const;
      static Message random_message ();

      friend std::ostream & operator<< (std::ostream &o,
                                        const thresholdECDSA &C);

    private:
      /* utils */
      OpenSSL::BN sum (const std::vector<OpenSSL::BN> &Operands) const;
      OpenSSL::BN hash_message (const Message &m) const;

      const SecLevel seclevel_;
      const OpenSSL::ECGroup ec_group_;
      const CL_HSMqk CL_HSMq_;
      mutable OpenSSL::HashAlgo H_;
  };

  #include "threshold_ECDSA.inl"

} /* BICYCL namespace */

#endif /* BICYCL_THRESHOLD_ECDSA_HPP */
