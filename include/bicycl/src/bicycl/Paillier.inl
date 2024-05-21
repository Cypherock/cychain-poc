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
#ifndef BICYCL_PAILLIER_INL
#define BICYCL_PAILLIER_INL

/*
 */
inline
Paillier::Paillier (size_t N_nbits) : N_nbits_ (N_nbits)
{
}

/*
 */
inline
Paillier::Paillier (SecLevel seclevel)
  : Paillier(seclevel.RSA_modulus_bitsize())
{
}

/* */
inline
Paillier::SecretKey::SecretKey (const Paillier &C, RandGen &randgen)
  : p_(randgen.random_prime (C.N_nbits_/2)),
    q_(randgen.random_prime ((C.N_nbits_+1)/2))
{
  Mpz pm1, qm1, N;
  Mpz::mul (N, p_, q_);
  Mpz::sub (pm1, p_, 1UL);
  Mpz::sub (qm1, q_, 1UL);
  Mpz::lcm (lambda_, pm1, qm1);
  Mpz::mod_inverse (lambda_inv_, lambda_, N);
}

/* */
inline
const Mpz & Paillier::SecretKey::p () const
{
  return p_;
}

/* */
inline
const Mpz & Paillier::SecretKey::q () const
{
  return q_;
}

/* */
inline
const Mpz & Paillier::SecretKey::lambda () const
{
  return lambda_;
}

/* */
inline
const Mpz & Paillier::SecretKey::lambda_inv () const
{
  return lambda_inv_;
}

/* */
inline
Paillier::SecretKey Paillier::keygen (RandGen &randgen) const
{
  return SecretKey (*this, randgen);
}

/* */
inline
Paillier::PublicKey::PublicKey (const SecretKey &sk)
{
  Mpz::mul (N_, sk.p(), sk.q());
  Mpz::mul (NN_, N_, N_);
}

/* */
inline
const Mpz & Paillier::PublicKey::N () const
{
  return N_;
}

/* */
inline
const Mpz & Paillier::PublicKey::N_square () const
{
  return NN_;
}

/* */
inline
Paillier::PublicKey Paillier::keygen (const SecretKey &sk) const
{
  return PublicKey (sk);
}

/* */
inline
Paillier::CipherText::CipherText (const PublicKey &pk, const ClearText &m,
                                  const Mpz &r)
{
  Mpz rN;
  Mpz M(1UL);

  Mpz::pow_mod (rN, r, pk.N(), pk.N_square());  /* r^N mod N^2 */
  Mpz::addmul (M, m, pk.N());                   /* 1+m*N */
  Mpz::mul (*this, M, rN);
  Mpz::mod (*this, *this, pk.N_square());       /* c = (1+m*N)*r^N mod N^2 */
}

/* */
inline
Paillier::CipherText Paillier::encrypt (const PublicKey &pk, const ClearText &m,
                                        const Mpz &r) const
{
  return CipherText (pk, m, r);
}

/*
 * Same as above but without the randomness
 */
inline
Paillier::CipherText Paillier::encrypt (const PublicKey &pk, const ClearText &m,
                                        RandGen &randgen) const
{
  return encrypt (pk, m, randgen.random_mpz (encrypt_randomness_bound(pk)));
}


/* */
inline
Paillier::ClearText::ClearText (const Paillier &C, const PublicKey &pk,
                                const Mpz &v)
  : Mpz (v)
{
  if (!(v.sgn() >= 0 && v < C.cleartext_bound(pk)))
    throw std::range_error ("Cleartext is negative or too large");
}

/* */
inline
Paillier::ClearText::ClearText (const Paillier &C, const PublicKey &pk,
                                RandGen &randgen)
  : Mpz (randgen.random_mpz (C.cleartext_bound(pk)))
{
}

/* */
inline
Paillier::ClearText::ClearText (const PublicKey &pk, const SecretKey &sk,
                                const CipherText &c)
{
  Mpz::pow_mod (*this, c, sk.lambda(), pk.N_square());  /* c^lambda mod N^2 */
  Mpz::sub (*this, *this, 1UL);         /* c^lambda mod N^2 -1 = lambda*m*N */
  Mpz::divexact (*this, *this, pk.N()); /* lambda*m mod N */
  Mpz::mul (*this, *this, sk.lambda_inv());       /* mul by 1/lambda mod N */
  Mpz::mod (*this, *this, pk.N());
}

/* */
inline
Paillier::ClearText Paillier::decrypt (const PublicKey &pk, const SecretKey &sk,
                                       const CipherText &c) const
{
  return ClearText (pk, sk, c);
}

inline
size_t Paillier::N_nbits () const
{
  return N_nbits_;
}

inline
const Mpz & Paillier::encrypt_randomness_bound (const PublicKey & pk) const
{
  return pk.N();
}

inline
const Mpz & Paillier::cleartext_bound (const PublicKey & pk) const
{
  return pk.N();
}

/*
 */
inline
CamenischShoup::CamenischShoup (const Mpz &N, const Mpz &r)
  : N_(N)
{
  Mpz::mul (NN_, N_, N_);
  Mpz::divby4 (Nover4_, N_);
  Mpz::mul (gen_, r, r);
  Mpz::mod (gen_, gen_, NN_);
  Mpz::pow_mod (gen_, gen_, N_, NN_); /* gen = r^(2*N) mod N^2 */

  /*
   * Precomputation
   */
  e_ = (encrypt_randomness_bound().nbits () + 3)/4;
  gen_3e_precomp_ = gen_;
  for (size_t i = 0; i < 3*e_; i++)
  {
    if (i == e_)
      gen_e_precomp_ = gen_3e_precomp_;
    if (i == 2*e_)
      gen_2e_precomp_ = gen_3e_precomp_;
    Mpz::mul (gen_3e_precomp_, gen_3e_precomp_, gen_3e_precomp_);
    Mpz::mod (gen_3e_precomp_, gen_3e_precomp_, NN_);
  }
}

/*
 */
inline
CamenischShoup::CamenischShoup (const Mpz &N, RandGen &randgen)
  : CamenischShoup (N, CamenischShoup::random_r (randgen, N))
{
}

/*
 */
inline
CamenischShoup::CamenischShoup (size_t N_nbits, RandGen &randgen)
  : CamenischShoup (CamenischShoup::random_N (randgen, N_nbits), randgen)
{
}

/*
 */
inline
CamenischShoup::CamenischShoup (SecLevel seclevel, RandGen &randgen)
  : CamenischShoup (seclevel.RSA_modulus_bitsize(), randgen)
{
}

inline
void CamenischShoup::power_of_gen (Mpz &r, const Mpz &n) const
{
  Mpz::pow_mod (r, gen_, n, NN_, e_, gen_e_precomp_, gen_2e_precomp_,
                                                              gen_3e_precomp_);
}

/*
 */
inline
CamenischShoup::SecretKey CamenischShoup::keygen (RandGen &randgen) const
{
  return SecretKey (*this, randgen);
}

/*
 */
inline
CamenischShoup::PublicKey CamenischShoup::keygen (const SecretKey &sk) const
{
  return PublicKey (*this, sk);
}

/*
 */
inline
CamenischShoup::PublicKey::PublicKey (const CamenischShoup &C,
                                      const SecretKey &sk)
{
  Mpz::pow_mod (pk_, C.gen(), sk, C.Nsquare()); /* pk = gen^sk mod N^2 */

  /*
   * Precomputation
   */
  e_ = (C.encrypt_randomness_bound().nbits () + 3)/4;
  pk_3e_precomp_ = pk_;
  for (size_t i = 0; i < 3*e_; i++)
  {
    if (i == e_)
      pk_e_precomp_ = pk_3e_precomp_;
    if (i == 2*e_)
      pk_2e_precomp_ = pk_3e_precomp_;
    Mpz::mul (pk_3e_precomp_, pk_3e_precomp_, pk_3e_precomp_);
    Mpz::mod (pk_3e_precomp_, pk_3e_precomp_, C.Nsquare());
  }
}

/* */
inline
void CamenischShoup::PublicKey::exponentiation (const CamenischShoup &C, Mpz &r,
                                                const Mpz &n) const
{
  Mpz::pow_mod (r, pk_, n, C.Nsquare(), e_, pk_e_precomp_, pk_2e_precomp_,
                                                                pk_3e_precomp_);
}

/*
 */
inline
const Mpz & CamenischShoup::PublicKey::elt () const
{
  return pk_;
}

/*
 */
inline
CamenischShoup::CipherText CamenischShoup::encrypt (const PublicKey &pk,
                                                    const ClearText &m,
                                                    const Mpz &r) const
{
  return CipherText (*this, pk, m, r);
}

/*
 * Same as above but without the randomness
 */
inline
CamenischShoup::CipherText CamenischShoup::encrypt (const PublicKey &pk,
                                                    const ClearText &m,
                                                    RandGen &randgen) const
{
  return encrypt (pk, m, randgen.random_mpz (encrypt_randomness_bound()));
}

/*
 */
inline
CamenischShoup::CipherText::CipherText (const CamenischShoup &C,
                                        const PublicKey &pk,
                                        const ClearText &m, const Mpz &r)
{
#ifdef BICYCL_WITH_PTHREADS
  using FctType = void (CamenischShoup::*) (Mpz &, const Mpz &) const;
  auto fctptr = static_cast<FctType>(&CamenischShoup::power_of_gen);
  std::thread th (fctptr, C, std::ref(c1_), std::cref(r)); /* c1 = g^r */
#else
  C.power_of_gen (c1_, r); /* c1 = h^r */
#endif

  Mpz t;
  pk.exponentiation (C, c2_, r);
  Mpz::mul (t, m, C.N());
  Mpz::add (t, t, 1UL);
  Mpz::mod (t, t, C.Nsquare());
  Mpz::mul (c2_, c2_, t);
  Mpz::mod (c2_, c2_, C.Nsquare());

#ifdef BICYCL_WITH_PTHREADS
  th.join ();
#endif
}

/* */
inline
const Mpz & CamenischShoup::CipherText::c1 () const
{
  return c1_;
}

/* */
inline
const Mpz & CamenischShoup::CipherText::c2 () const
{
  return c2_;
}

/*
 */
inline
CamenischShoup::ClearText CamenischShoup::decrypt (const SecretKey &sk,
                                                   const CipherText &c) const
{
  Mpz t, m;
  Mpz::pow_mod (m, c.c1(), sk, NN_); /* c1^sk mod N^2 */
  Mpz::mod_inverse (m, m, NN_);
  Mpz::mul (m, m, c.c2());
  Mpz::sub (m, m, 1UL);
  Mpz::mod (m, m, NN_); /* (c2/c1^sk) - 1 mod N^2 */
  Mpz::divexact (m, m, N_);
  return ClearText (*this, m);
}

/*
 */
inline
const Mpz & CamenischShoup::gen () const
{
  return gen_;
}

/*
 */
inline
const Mpz & CamenischShoup::N () const
{
  return N_;
}

/*
 */
inline
const Mpz & CamenischShoup::Nsquare () const
{
  return NN_;
}

/*
 */
inline
const Mpz & CamenischShoup::encrypt_randomness_bound () const
{
  return Nover4_;
}

/*
 */
inline
const Mpz & CamenischShoup::cleartext_bound () const
{
  return N_;
}

/*
 */
inline
const Mpz & CamenischShoup::secretkey_bound () const
{
  return Nover4_;
}

/*
 */
inline
Mpz CamenischShoup::random_N (RandGen &randgen, size_t N_nbits)
{
  Mpz p, q, t;
  do
  {
    t = randgen.random_prime (N_nbits/2);
    Mpz::mulby2 (t, t);
    Mpz::add (p, t, 1UL);
  } while (p.is_prime());

  do
  {
    t = randgen.random_prime (N_nbits/2);
    Mpz::mulby2 (t, t);
    Mpz::add (q, t, 1UL);
  } while (q.is_prime());

  Mpz::mul (t, p, q);
  return t;
}

/*
 */
inline
Mpz CamenischShoup::random_r (RandGen &randgen, const Mpz &N)
{
  Mpz NN, r, g;
  Mpz::mul (NN, N, N);
  do
  {
    r = randgen.random_mpz (NN);
    Mpz::gcd (g, r, N);
  } while (g != 1UL);
  return r;
}

#endif /* BICYCL_PAILLIER_INL */
