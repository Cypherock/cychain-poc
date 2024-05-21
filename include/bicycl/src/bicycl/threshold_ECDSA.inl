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
#ifndef BICYCL_THRESHOLD_ECDSA_INL
#define BICYCL_THRESHOLD_ECDSA_INL

/******************************************************************************/
class thresholdECDSA::ProtocolAbortError : public std::runtime_error
{
  public:
    using runtime_error::runtime_error;
};

/******************************************************************************/
inline
thresholdECDSA::KeygenPart1::KeygenPart1 (const thresholdECDSA &C,
                                          unsigned int n, unsigned int t,
                                          unsigned int i)
  : i_(i),
    u_(C.ec_group_.random_mod_order()),
    Q_(C.ec_group_, u_),
    c_(),
    cs_(),
    a_(),
    V_(),
    sigma_(n)
{
  /* check */
  if (n < 2)
  {
    throw std::runtime_error ("n must be >= 2 in thresholdECDSA keygen");
  }
  if (t >= n)
  {
    throw std::runtime_error ("t must be < n in thresholdECDSA keygen");
  }
  if (t < 1)
  {
    throw std::runtime_error ("t must be >= 1 in thresholdECDSA keygen");
  }
  if (i >= n)
  {
    throw std::runtime_error ("i must be < n in thresholdECDSA keygen");
  }

  /* compute commitment */
  std::tie(c_, cs_) = C.commit (Q_);

  /* draw a_k, compute V_k as [a_k]P */
  a_.reserve (t);
  V_.reserve (t);
  for (unsigned int k = 0; k < t; k++)
  {
    a_.push_back (C.ec_group_.random_mod_order());
    V_.push_back (OpenSSL::ECPoint (C.ec_group_, a_[k]));
  }

  /* compute sigma_j values */
  OpenSSL::BN v;
  for (unsigned int j = 0; j < n; j++)
  {
    /* Compute sigma_[j] = u + sum_{k=1}^{t_}{V_k (j+1)^k} using
     * Horner's method.
     */
    sigma_[j] = a_[t-1];
    for (unsigned int k = t-1; k > 0; k--)
    {
      C.ec_group_.mul_by_word_mod_order (sigma_[j], j+1);
      C.ec_group_.add_mod_order (sigma_[j], sigma_[j], a_[k-1]);
    }
    C.ec_group_.mul_by_word_mod_order (sigma_[j], j+1);
    C.ec_group_.add_mod_order (sigma_[j], sigma_[j], u_);
  }
}

/* */
inline
unsigned int thresholdECDSA::KeygenPart1::n () const
{
  return sigma_.size();
}

/* */
inline
unsigned int thresholdECDSA::KeygenPart1::t () const
{
  return a_.size();
}

/* */
inline
unsigned int thresholdECDSA::KeygenPart1::i () const
{
  return i_;
}

/* */
inline
const thresholdECDSA::Commitment & thresholdECDSA::KeygenPart1::commitment () const
{
  return c_;
}

/* */
inline
const thresholdECDSA::CommitmentSecret & thresholdECDSA::KeygenPart1::commitment_secret () const
{
  return cs_;
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::KeygenPart1::Q_part () const
{
  return Q_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::KeygenPart1::u_part () const
{
  return u_;
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::KeygenPart1::V (size_t k) const
{
  return V_[k];
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::KeygenPart1::sigma (size_t j) const
{
  return sigma_[j];
}

/******************************************************************************/
inline
thresholdECDSA::KeygenPart2::KeygenPart2 (const thresholdECDSA &C,
                            const KeygenPart1 &data1, RandGen &randgen,
                            const std::vector<Commitment> &Co,
                            const std::vector<OpenSSL::ECPoint> &Q,
                            const std::vector<CommitmentSecret> &CoSec,
                            const std::vector<std::vector<OpenSSL::ECPoint>> &V,
                            const std::vector<OpenSSL::BN> &Sigma)
  : Q_(C.ec_group_, data1.Q_part()),
    x_(C.sum(Sigma)),
    zk_proof_(C.ec_group_, C.H_, x_),
    sk_(C.CL_HSMq_, randgen),
    pk_(C.CL_HSMq_, sk_)
{
  const unsigned int n = data1.n();
  const unsigned int t = data1.t();
  const unsigned int i = data1.i();

  /* check */
  if (Co.size() != n || Q.size() != n || CoSec.size() != n || V.size() != n
      || Sigma.size() != n)
  {
    throw std::runtime_error ("Vectors must have size n in thresholdECDSA keygen");
  }
  if (std::any_of (V.cbegin(), V.cend(),
      [t](const std::vector<OpenSSL::ECPoint> & Vi) { return Vi.size() != t; }))
  {
    throw std::runtime_error ("Vi must have size t in thresholdECDSA keygen");
  }

  /* open commitments */
  for (unsigned int j = 0; j < n; j++)
  {
    if (j != i)
    {
      bool b = C.open (Co[j], Q[j], CoSec[j]);
      if (!b)
      {
        throw ProtocolAbortError ("cannot verify commitment");
      }
      C.ec_group_.ec_add (Q_, Q_, Q[j]);
    }
    /* Note that Q_common_ was initialized with Qi */
  }

  /* check that [Sigma[j]]P = Q[j] + sum_{k=1}^{t}{i^k V[j][k]}*/
  OpenSSL::ECPoint T1 (C.ec_group_), T2 (C.ec_group_), T (C.ec_group_);
  OpenSSL::BN m;
  for (unsigned int j = 0; j < n; j++)
  {
    if (j != i)
    {
      C.ec_group_.scal_mul_gen (T1, Sigma[j]);
      T2 = Q[j];
      m = (i+1);
      for (unsigned int k = 0; k < t; k++)
      {
        C.ec_group_.scal_mul (T, m, V[j][k]);
        C.ec_group_.ec_add (T2, T2, T);
        m *= (i+1);
      }
      if (!C.ec_group_.ec_point_eq (T1, T2))
      {
        throw ProtocolAbortError ("cannot verify sigma_{j,i} equality");
      }
    }
  }
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::KeygenPart2::Q () const
{
  return Q_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::KeygenPart2::x () const
{
  return x_;
}

/* */
inline
const ECNIZKProof & thresholdECDSA::KeygenPart2::zk_proof () const
{
  return zk_proof_;
}

/* */
inline
const CL_HSMqk::PublicKey & thresholdECDSA::KeygenPart2::CL_public_key () const
{
  return pk_;
}

/* */
inline
const CL_HSMqk::SecretKey & thresholdECDSA::KeygenPart2::CL_secret_key () const
{
  return sk_;
}

/******************************************************************************/
inline
thresholdECDSA::SecretKey::SecretKey (const thresholdECDSA &C, unsigned int i,
                            const KeygenPart1 &data1, const KeygenPart2 &data2,
                            const std::vector<std::vector<OpenSSL::ECPoint>> &V,
                            const std::vector<ECNIZKProof> &ZK,
                            const std::vector<CL_HSMqk::PublicKey> &PK)
  : sk_(data2.CL_secret_key()),
    x_(data2.x()),
    PK_(PK),
    Q_(C.ec_group_, data2.Q())
{
  const unsigned int n = data1.n();
  const unsigned int t = data1.t();

  /* check */
  if (i >= n)
  {
    throw std::runtime_error ("index i must be < n in thresholdECDSA keygen");
  }
  if (V.size() != n || ZK.size() != n || PK.size() != n)
  {
    throw std::runtime_error ("Vectors must have size n in thresholdECDSA keygen");
  }
  if (std::any_of (V.cbegin(), V.cend(),
      [t](const std::vector<OpenSSL::ECPoint> & Vi) { return Vi.size() != t; }))
  {
    throw std::runtime_error ("Vi must have size t in thresholdECDSA keygen");
  }

  X_.reserve (n);
  OpenSSL::ECPoint Xj (C.ec_group_), T (C.ec_group_);
  OpenSSL::BN m;
  for (unsigned int j = 0; j < n; j++)
  {
    if (j != i)
    {
      /* Compute X_j as Q + sum_{l=1}^{n}{sum_{k=1}^{t}{j^k V_l,k}} */
      Xj = data2.Q();
      m = (j+1);
      for (unsigned int k = 0; k < t; k++)
      {
        /* we know that n >= 2 */
        C.ec_group_.ec_add (T, V[0][k], V[1][k]);
        for (unsigned int l = 2; l < n; l++)
        {
          C.ec_group_.ec_add (T, T, V[l][k]);
        }
        C.ec_group_.scal_mul (T, m, T);
        C.ec_group_.ec_add (Xj, Xj, T);
        m *= (j+1);
      }

      /* check zero-knowledge proof */
      bool b = ZK[j].verify (C.ec_group_, C.H_, Xj);
      if (!b)
      {
        throw ProtocolAbortError ("cannot verify ZK proof");
      }
      X_.push_back (OpenSSL::ECPoint (C.ec_group_, Xj));
    }
    else
    {
      X_.push_back (OpenSSL::ECPoint (C.ec_group_));/* XXX Do we need the actual point ? */
    }
  }
}


/* */
inline
const thresholdECDSA::PublicKey & thresholdECDSA::SecretKey::public_key () const
{
  return Q_;
}

/* */
inline
const CL_HSMqk::SecretKey & thresholdECDSA::SecretKey::CL_secret_key () const
{
  return sk_;
}

/* */
inline
const CL_HSMqk::PublicKey & thresholdECDSA::SecretKey::CL_public_key (
                                                          unsigned int i) const
{
  return PK_[i];
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::SecretKey::X (unsigned int i) const
{
  return X_[i];
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SecretKey::x_part () const
{
  return x_;
}

/******************************************************************************/
inline
thresholdECDSA::SignPart1::SignPart1 (const thresholdECDSA &C,
                                      unsigned int i, const ParticipantsList &S,
                                      const thresholdECDSA::SecretKey &sk,
                                      RandGen &randgen)
  : i_ (i),
    S_ (S),
    gamma_(C.ec_group_.random_mod_order()),
    Gamma_(C.ec_group_, gamma_),
    zk_gamma_ (C.ec_group_, C.H_, gamma_, Gamma_), co_(), cos_(),
    k_ (C.CL_HSMq_, randgen),
    r_ (randgen.random_mpz (C.CL_HSMq_.encrypt_randomness_bound())),
    c_ (C.CL_HSMq_, sk.CL_public_key (i), k_, r_),
    zk_encrypt_ (C.CL_HSMq_, C.H_, sk.CL_public_key (i), c_, k_, r_, randgen)
{
  /* check */
  if (std::find (S.begin(), S.end(), i) == S.end()) /* if i not in S */
  {
    throw std::runtime_error ("i must be in the list of participants S");
  }

  /* */
  OpenSSL::BN tmp (C.lagrange_at_zero (S, i));
  C.ec_group_.mul_mod_order (tmp, tmp, sk.x_part());
  omega_ = tmp;

  /* compute commitment */
  std::tie(co_, cos_) = C.commit (Gamma_);
}

/* */
inline
const thresholdECDSA::ParticipantsList & thresholdECDSA::SignPart1::S () const
{
  return S_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart1::gamma () const
{
  return gamma_;
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::SignPart1::Gamma () const
{
  return Gamma_;
}

/* */
inline
const ECNIZKProof & thresholdECDSA::SignPart1::zk_gamma () const
{
  return zk_gamma_;
}

/* */
inline
const thresholdECDSA::Commitment & thresholdECDSA::SignPart1::commitment () const
{
  return co_;
}

/* */
inline
const thresholdECDSA::CommitmentSecret & thresholdECDSA::SignPart1::commitment_secret () const
{
  return cos_;
}

/* */
inline
unsigned int thresholdECDSA::SignPart1::i () const
{
  return i_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart1::omega () const
{
  return omega_;
}

/* */
inline
const Mpz & thresholdECDSA::SignPart1::k_part () const
{
  return k_;
}

/* */
inline
const CL_HSMqk::CipherText & thresholdECDSA::SignPart1::ciphertext () const
{
  return c_;
}

/* */
inline
const CL_HSMqk_ZKAoKProof & thresholdECDSA::SignPart1::zk_encrypt_proof () const
{
  return zk_encrypt_;
}

/******************************************************************************/
/* */
inline
thresholdECDSA::SignPart2::SignPart2 (const thresholdECDSA &C,
                            const SignPart1 &data,
                            const thresholdECDSA::SecretKey &sk,
                            const ParticipantsMap<Commitment> & commitment_map,
                            const ParticipantsMap<CL_HSMqk::CipherText> &c_map,
                            const ParticipantsMap<CL_HSMqk_ZKAoKProof> &zk_map,
                            RandGen &randgen)
  : commitment_map_ (commitment_map)
{
  unsigned int Slen = data.S().size();
  unsigned int i = data.i();

  nu_map_.reserve (Slen);
  B_map_.reserve (Slen);
  beta_map_.reserve (Slen);
  c_kg_map_.reserve (Slen);
  c_kw_map_.reserve (Slen);

  Mpz omega_i (data.omega());

  for (unsigned int j: data.S())
  {
    if (j != i)
    {
      const CL_HSMqk_ZKAoKProof &zkj = zk_map.at (j);
      const CL_HSMqk::PublicKey &pkj = sk.CL_public_key (j);
      const CL_HSMqk::CipherText &cj = c_map.at (j);
      bool bj = zkj.verify (C.CL_HSMq_, C.H_, pkj, cj);
      if (!bj)
      {
        throw ProtocolAbortError ("cannot verify NIZKAoK");
      }

      nu_map_.emplace (j, C.ec_group_.random_mod_order());
      B_map_.emplace (j, OpenSSL::ECPoint(C.ec_group_, nu_map_.at(j)));

      CL_HSMqk::ClearText beta (C.CL_HSMq_, randgen);
      CL_HSMqk::CipherText c_beta (C.CL_HSMq_.encrypt (pkj, beta, randgen));
      Mpz::sub (beta, C.ec_group().order(), beta); /* -beta modulo q */
      beta_map_.emplace (j, beta);

      Mpz gamma (data.gamma());
      CL_HSMqk::CipherText cj_scal (C.CL_HSMq_.scal_ciphertexts (pkj, cj, gamma,
                                                                      randgen));
      c_kg_map_.emplace (j, C.CL_HSMq_.add_ciphertexts (pkj, cj_scal, c_beta,
                                                                      randgen));

      Mpz nu (nu_map_.at (j));
      CL_HSMqk::ClearText minus_nu (C.CL_HSMq_, nu);
      Mpz::sub (minus_nu, C.ec_group().order(), minus_nu); /* neg mod q */
      CL_HSMqk::CipherText c_nu (C.CL_HSMq_.encrypt (pkj, minus_nu, randgen));

      CL_HSMqk::CipherText cj_scal2 (C.CL_HSMq_.scal_ciphertexts (pkj, cj,
                                                            omega_i, randgen));
      c_kw_map_.emplace (j, C.CL_HSMq_.add_ciphertexts (pkj, cj_scal2, c_nu,
                                                                      randgen));
    }
  }
}

/* */
inline
const thresholdECDSA::Commitment & thresholdECDSA::SignPart2::commitment (unsigned int j) const
{
  return commitment_map_.at (j);
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart2::nu (unsigned int j) const
{
  return nu_map_.at (j);
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::SignPart2::B (unsigned int j) const
{
  return B_map_.at (j);
}

/* */
inline
const CL_HSMqk::ClearText & thresholdECDSA::SignPart2::beta (
                                                          unsigned int j) const
{
  return beta_map_.at (j);
}

/* */
inline
const CL_HSMqk::CipherText & thresholdECDSA::SignPart2::c_kg (
                                                          unsigned int j) const
{
  return c_kg_map_.at (j);
}
/* */
inline
const CL_HSMqk::CipherText & thresholdECDSA::SignPart2::c_kw (
                                                          unsigned int j) const
{
  return c_kw_map_.at (j);
}

/******************************************************************************/
/* */
inline
thresholdECDSA::SignPart3::SignPart3 (const thresholdECDSA &C,
                         const SignPart1 &data1, const SignPart2 &data2,
                         const thresholdECDSA::SecretKey &sk,
                         const ParticipantsMap<CL_HSMqk::CipherText> &c_kg_map,
                         const ParticipantsMap<CL_HSMqk::CipherText> &c_kw_map,
                         const ParticipantsMap<OpenSSL::ECPoint> &B_map)
  : delta_(),
    sigma_()
{
  unsigned int i = data1.i();

  OpenSSL::BN tmp;
  OpenSSL::BN ki (data1.k_part());

  C.ec_group_.mul_mod_order (delta_, ki, data1.gamma());
  C.ec_group_.mul_mod_order (sigma_, ki, data1.omega());

  for (unsigned int j: data1.S())
  {
    if (j != i)
    {
      CL_HSMqk::ClearText alpha (C.CL_HSMq_, sk.CL_secret_key(), c_kg_map.at(j));
      CL_HSMqk::ClearText mu (C.CL_HSMq_, sk.CL_secret_key(), c_kw_map.at(j));

      OpenSSL::ECPoint T (C.ec_group_, B_map.at(j));
      OpenSSL::ECPoint T2 (C.ec_group_);

      tmp = mu;
      C.ec_group_.scal_mul_gen (T2, tmp); /* mu_i,j P */
      C.ec_group_.ec_add (T, T, T2); /* mu_i,j P + B_j,i */

      C.ec_group_.add_mod_order (sigma_, sigma_, tmp); /* += mu_i,j */
      tmp = data2.nu(j);
      C.ec_group_.add_mod_order (sigma_, sigma_, tmp); /* += nu_i,j */

      tmp = C.lagrange_at_zero (data1.S(), j);
      C.ec_group_.mul_mod_order (tmp, tmp, ki);
      C.ec_group_.scal_mul (T2, tmp, sk.X(j)); /* (ki * lj(0)) Xj */
      bool b = C.ec_group_.ec_point_eq (T, T2);
      if (!b)
      {
        throw ProtocolAbortError ("cannot verify mu_i,j P+B_ij =k_i W_j");
      }

      tmp = alpha;
      C.ec_group_.add_mod_order (delta_, delta_, tmp); /* += alpha_i,j */
      tmp = data2.beta (j);
      C.ec_group_.add_mod_order (delta_, delta_, tmp); /* += beta_i,j */
    }
  }
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart3::delta_part () const
{
  return delta_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart3::sigma_part () const
{
  return sigma_;
}

/******************************************************************************/
/* */
inline
thresholdECDSA::SignPart4::SignPart4 (const thresholdECDSA &C,
                                  const SignPart1 &data1,
                                  const ParticipantsMap<OpenSSL::BN> &delta_map)
  : delta_()
{
  delta_ = 0UL;
  for (unsigned int j: data1.S())
  {
    C.ec_group().add_mod_order (delta_, delta_, delta_map.at(j));
  }
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart4::delta () const
{
  return delta_;
}

/******************************************************************************/
/* */
inline
thresholdECDSA::SignPart5::SignPart5 (const thresholdECDSA &C,
                            const SignPart1 &data1, const SignPart2 &data2,
                            const SignPart3 &data3, const SignPart4 &data4,
                            const Message &m,
                            const ParticipantsMap<OpenSSL::ECPoint> &Gamma_map,
                            const ParticipantsMap<CommitmentSecret> &CoSec_map,
                            const ParticipantsMap<ECNIZKProof> &zk_proof_map)
  : m_(m),
    R_(C.ec_group_),
    r_(),
    z_(C.hash_message (m)),
    ell_(C.ec_group_.random_mod_order()),
    rho_(C.ec_group_.random_mod_order()),
    V_(C.ec_group_),
    A_(C.ec_group_, rho_),
    c_(),
    cs_()
{
  unsigned int i = data1.i();

  OpenSSL::ECPoint Gamma (C.ec_group_, data1.Gamma());

  for (unsigned int j: data1.S())
  {
    if (j != i)
    {
      /* open commitment */
      if (!C.open (data2.commitment(j), Gamma_map.at(j), CoSec_map.at(j)))
      {
        throw ProtocolAbortError ("cannot verify commitment");
      }
      /* check zero-knowledge proof */
      bool b = zk_proof_map.at(j).verify (C.ec_group_, C.H_, Gamma_map.at(j));
      if (!b)
      {
        throw ProtocolAbortError ("cannot verify ZK proof");
      }

      C.ec_group_.ec_add (Gamma, Gamma, Gamma_map.at(j));
    }
  }

  OpenSSL::BN tmp;

  C.ec_group_.inverse_mod_order (tmp, data4.delta());
  C.ec_group_.scal_mul (R_, tmp, Gamma); /* R = delta^(-1) Gamma */
  C.ec_group_.x_coord_of_point (r_, R_);
  C.ec_group_.mod_order (r_, r_);
  /* the condition r != 0 will be checked at the end */

  /* s_i = z * k_i + r * sigma_i*/
  tmp = data1.k_part();
  C.ec_group_.mul_mod_order (s_, z_, tmp);
  C.ec_group_.mul_mod_order (tmp, r_, data3.sigma_part());
  C.ec_group_.add_mod_order (s_, tmp, s_);

  /* V_i = s_i R + l_i P = s_i R + L_i */
  C.ec_group_.scal_mul (V_, ell_, s_, R_);

  /* compute commitment */
  std::tie(c_, cs_) = C.commit (V_, A_);
}

/* */
inline
const thresholdECDSA::Message & thresholdECDSA::SignPart5::m () const
{
  return m_;
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::SignPart5::R () const
{
  return R_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart5::r () const
{
  return r_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart5::z () const
{
  return z_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart5::s_part () const
{
  return s_;
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::SignPart5::V_part () const
{
  return V_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart5::ell () const
{
  return ell_;
}

/* */
inline
const OpenSSL::BN & thresholdECDSA::SignPart5::rho () const
{
  return rho_;
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::SignPart5::A_part () const
{
  return A_;
}

/* */
inline
const thresholdECDSA::Commitment & thresholdECDSA::SignPart5::commitment () const
{
  return c_;
}

/* */
inline
const thresholdECDSA::CommitmentSecret & thresholdECDSA::SignPart5::commitment_secret () const
{
  return cs_;
}

/******************************************************************************/
/* */
inline
thresholdECDSA::SignPart6::SignPart6 (const thresholdECDSA &C,
                                      const SignPart5 &data5,
                                      const ParticipantsMap<Commitment> &Co_map)
  : commitment_map_(Co_map),
    zk_aok_(C.ec_group_, C.H_, data5.R(), data5.s_part(),
            data5.ell(), data5.rho(), data5.V_part(),
            data5.A_part())
{
}

/* */
inline
const thresholdECDSA::Commitment & thresholdECDSA::SignPart6::commitment (unsigned int j) const
{
  return commitment_map_.at (j);
}

/* */
inline
const ECNIZKAoK & thresholdECDSA::SignPart6::aok () const
{
  return zk_aok_;
}

/******************************************************************************/
/* */
inline
thresholdECDSA::SignPart7::SignPart7 (const thresholdECDSA &C,
                            const SignPart1 &data1, const SignPart5 &data5,
                            const SignPart6 &data6,
                            const thresholdECDSA::SecretKey &sk,
                            const ParticipantsMap<OpenSSL::ECPoint> &V_map,
                            const ParticipantsMap<OpenSSL::ECPoint> &A_map,
                            const ParticipantsMap<CommitmentSecret> &CoSec_map,
                            const ParticipantsMap<ECNIZKAoK> &zk_aok_map)
  : U_(C.ec_group_),
    T_(C.ec_group_)
{
  unsigned int i = data1.i();

  OpenSSL::ECPoint V (C.ec_group_);
  OpenSSL::ECPoint A (C.ec_group_, data5.A_part());

  C.ec_group_.scal_mul (V, data5.z(), data5.r(), sk.public_key());
  C.ec_group_.ec_neg (V); /* V = - z P - r Q */
  C.ec_group_.ec_add (V, V, data5.V_part());

  for (unsigned int j: data1.S())
  {
    if (j != i)
    {
      /* open commitment */
      bool b = C.open (data6.commitment(j), V_map.at(j), A_map.at(j),
                                                             CoSec_map.at(j));
      if (!b)
      {
        throw ProtocolAbortError ("cannot verify commitment");
      }

      /* check zero-knowledge */
      b = zk_aok_map.at(j).verify (C.ec_group_, C.H_, data5.R(),
                                                      V_map.at(j), A_map.at(j));
      if (!b)
      {
        throw ProtocolAbortError ("cannot verify ZK proof");
      }

      C.ec_group_.ec_add (V, V, V_map.at(j));
      C.ec_group_.ec_add (A, A, A_map.at(j));
    }
  }

  C.ec_group_.scal_mul (U_, data5.rho(), V);
  C.ec_group_.scal_mul (T_, data5.ell(), A);

  /* compute commitment */
  std::tie(c_, cs_) = C.commit (U_, T_);
}

/* */
inline
const thresholdECDSA::Commitment & thresholdECDSA::SignPart7::commitment () const
{
  return c_;
}

/* */
inline
const thresholdECDSA::CommitmentSecret & thresholdECDSA::SignPart7::commitment_secret () const
{
  return cs_;
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::SignPart7::U_part () const
{
  return U_;
}

/* */
inline
const OpenSSL::ECPoint & thresholdECDSA::SignPart7::T_part () const
{
  return T_;
}


/******************************************************************************/
/* */
inline
thresholdECDSA::SignPart8::SignPart8 (const thresholdECDSA &C,
                            const SignPart1 &data1, const SignPart7 &data7,
                            const ParticipantsMap<Commitment> &Co_map,
                            const ParticipantsMap<OpenSSL::ECPoint> &U_map,
                            const ParticipantsMap<OpenSSL::ECPoint> &T_map,
                            const ParticipantsMap<CommitmentSecret> &CoSec_map)
{
  unsigned int i = data1.i();

  OpenSSL::ECPoint T (C.ec_group_, data7.T_part());
  OpenSSL::ECPoint U (C.ec_group_, data7.U_part());

  for (unsigned int j: data1.S())
  {
    if (j != i)
    {
      /* open commitment */
      bool b = C.open (Co_map.at(j), U_map.at(j), T_map.at(j), CoSec_map.at(j));
      if (!b)
      {
        throw ProtocolAbortError ("cannot verify commitment");
      }

      C.ec_group_.ec_add (T, T, T_map.at(j));
      C.ec_group_.ec_add (U, U, U_map.at(j));
    }
  }

  if (!C.ec_group_.ec_point_eq (U, T))
  {
    throw ProtocolAbortError ("cannot verify equality U == T");
  }
}

/******************************************************************************/
/* */
inline
thresholdECDSA::Signature::Signature (const thresholdECDSA &C,
                                      const SignPart1 &data1,
                                      const SignPart5 &data5,
                                      const SecretKey &sk,
                                      const ParticipantsMap<OpenSSL::BN> &s_map)
  : r_ (data5.r()),
    s_ (data5.s_part())
{
  unsigned int i = data1.i();

  for (unsigned int j: data1.S())
  {
    if (j != i)
    {
      C.ec_group_.add_mod_order (s_, s_, s_map.at(j));
    }
  }

  bool b = verify (C, sk.public_key(), data5.m());
  if (!b)
  {
    throw ProtocolAbortError ("cannot verify signature");
  }
}

/* */
inline
bool thresholdECDSA::Signature::verify (const thresholdECDSA &C,
                                        const PublicKey &Q,
                                        const Message &m) const
{
  OpenSSL::BN z, sinv, u1, u2, x1, tmp;

  const OpenSSL::ECGroup & E = C.ec_group_;

  if (!E.has_correct_prime_order (Q)) /* check that Q as order n */
    return false;

  if (!E.is_positive_less_than_order (r_))
    return false;

  if (!E.is_positive_less_than_order (s_))
    return false;

  bool ok = true;
  OpenSSL::ECPoint T (E);
  z = C.hash_message (m);
  E.inverse_mod_order (sinv, s_);
  E.mul_mod_order (u1, sinv, z);
  E.mul_mod_order (u2, sinv, r_);

  E.scal_mul (T, u1, u2, Q); /* u1*G + u2*Q */

  if (E.is_at_infinity (T))
    ok = false;
  else
  {
    E.x_coord_of_point (tmp, T);
    E.mod_order (x1, tmp);

    ok = (x1 == r_);
  }

  return ok;
}

/* */
inline
bool thresholdECDSA::Signature::operator== (const Signature &other) const
{
  return r_ == other.r_ && s_ == other.s_;
}

/* */
inline
bool thresholdECDSA::Signature::operator!= (const Signature &other) const
{
  return !(*this == other);
}

/******************************************************************************/
/* */
inline
thresholdECDSA::thresholdECDSA (SecLevel seclevel, RandGen &randgen)
  : seclevel_(seclevel),
    ec_group_ (seclevel_),
    CL_HSMq_ (ec_group_.order(), 1, seclevel_, randgen),
    H_(seclevel_)
{
}

/* */
inline
const OpenSSL::ECGroup & thresholdECDSA::ec_group () const
{
  return ec_group_;
}

/* */
inline
std::tuple<thresholdECDSA::Commitment, thresholdECDSA::CommitmentSecret>
thresholdECDSA::commit (const OpenSSL::ECPoint &Q) const
{
  size_t nbytes = static_cast<unsigned int>(seclevel_) >> 3; /* = seclevel/8 */
  CommitmentSecret r(nbytes);
  OpenSSL::random_bytes (r.data(), nbytes);
  return std::make_tuple (H_ (r, OpenSSL::ECPointGroupCRefPair (Q, ec_group_)),
                          r);
}

/* */
inline
std::tuple<thresholdECDSA::Commitment, thresholdECDSA::CommitmentSecret>
thresholdECDSA::commit (const OpenSSL::ECPoint &Q1,
                        const OpenSSL::ECPoint &Q2) const
{
  size_t nbytes = static_cast<unsigned int>(seclevel_) >> 3; /* = seclevel/8 */
  CommitmentSecret r(nbytes);
  OpenSSL::random_bytes (r.data(), nbytes);
  return std::make_tuple (H_ (r, OpenSSL::ECPointGroupCRefPair (Q1, ec_group_),
                                 OpenSSL::ECPointGroupCRefPair (Q2, ec_group_)),
                          r);
}

/* */
inline
bool thresholdECDSA::open (const Commitment &c, const OpenSSL::ECPoint &Q,
                                                const CommitmentSecret &r) const
{
  Commitment c2 (H_ (r, OpenSSL::ECPointGroupCRefPair (Q, ec_group_)));
  return c == c2;
}

/* */
inline
bool thresholdECDSA::open (const Commitment &c, const OpenSSL::ECPoint &Q1,
                                                const OpenSSL::ECPoint &Q2,
                                                const CommitmentSecret &r) const
{
  Commitment c2 (H_ (r, OpenSSL::ECPointGroupCRefPair (Q1, ec_group_),
                        OpenSSL::ECPointGroupCRefPair (Q2, ec_group_)));
  return c == c2;
}

/* */
inline
bool thresholdECDSA::verify (const Signature &s, const PublicKey &Q,
                                                 const Message &m) const
{
  return s.verify (*this, Q, m);
}

/* */
inline
thresholdECDSA::Message thresholdECDSA::random_message ()
{
  return ECDSA::random_message ();
}

/* */
inline
std::ostream & operator<< (std::ostream &o, const thresholdECDSA &C)
{
  return o << C.CL_HSMq_ << C.ec_group_;
}

/* */
inline
OpenSSL::BN thresholdECDSA::sum (const std::vector<OpenSSL::BN> &Operands) const
{
  OpenSSL::BN r;
  r = 0UL;
  for (auto const &v: Operands)
  {
    ec_group_.add_mod_order (r, r, v);
  }
  return r;
}

/*
 * returns prod_{s \in S}{\frac{-(s+1)}{(s+1)-(i+1)}}
 * Note: in code participants ids are 0-based, in the protocol they are 1-based,
 * so we need to add 1 to all participants ids.
 */
inline
OpenSSL::BN thresholdECDSA::lagrange_at_zero (const ParticipantsList &S,
                                              unsigned int i) const
{
  OpenSSL::BN num, den, r;

  num = 1UL;
  den = 1UL;
  for (unsigned int s: S)
  {
    if (s != i)
    {
      ec_group_.mul_by_word_mod_order (num, s+1);
      if (s > i)
      {
        ec_group_.mul_by_word_mod_order (den, s-i);
      }
      else
      {
        ec_group_.mul_by_word_mod_order (den, i-s);
        den.neg();
      }
    }
  }
  ec_group_.inverse_mod_order (r, den);
  ec_group_.mul_mod_order (r, r, num);

  return r;
}

/* */
inline
OpenSSL::BN thresholdECDSA::hash_message (const Message &m) const
{
  return OpenSSL::BN (H_ (m));
}

#endif /* BICYCL_THRESHOLD_ECDSA_INL */
