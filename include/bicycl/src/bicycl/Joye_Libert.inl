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
#ifndef BICYCL_JOYE_LIBERT_INL
#define BICYCL_JOYE_LIBERT_INL

/*
 */
inline
JoyeLibert::JoyeLibert (size_t N_nbits, size_t k) : N_nbits_(N_nbits), k_(k)
{
  twotothek_ = 1UL; Mpz::mulby2k (twotothek_, twotothek_, k); // FIXME temporary
}

/*
 */
inline
JoyeLibert::JoyeLibert (SecLevel seclevel, size_t k)
  : JoyeLibert(seclevel.RSA_modulus_bitsize(), k)
{
  // TODO check that k < log (N)/4 - seclevel
}

/* */
inline
size_t JoyeLibert::k () const
{
  return k_;
}

/* */
inline
size_t JoyeLibert::N_nbits () const
{
  return N_nbits_;
}

/* */
inline
JoyeLibert::SecretKey::SecretKey (const JoyeLibert &C, RandGen &randgen)
{
  size_t r_nbits = C.N_nbits_/2-C.k_;
  do
  {
    r_ = randgen.random_mpz_2exp (r_nbits);
    // TODO ensure r has r_nbits ie set MSB to 1
    Mpz::mulby2k (p_, r_, C.k_);
    Mpz::add (p_, p_, 1UL);
  } while (!p_.is_prime());

  q_ = randgen.random_prime ((C.N_nbits_+1)/2);
  Mpz N;
  Mpz::mul (N, p_, q_);

  do
  {
    y_ = randgen.random_mpz (N);
  } while (y_.kronecker (p_) != -1 || y_.kronecker (q_) != -1);

  Mpz e(r_);
  e.neg();
  Mpz::pow_mod (D_, y_, e, p_);
}

/* */
inline
const Mpz & JoyeLibert::SecretKey::p () const
{
  return p_;
}

/* */
inline
const Mpz & JoyeLibert::SecretKey::r () const
{
  return r_;
}

/* */
inline
const Mpz & JoyeLibert::SecretKey::q () const
{
  return q_;
}

/* */
inline
const Mpz & JoyeLibert::SecretKey::D () const
{
  return D_;
}

/* */
inline
const Mpz & JoyeLibert::SecretKey::y () const
{
  return y_;
}

/* */
inline
JoyeLibert::SecretKey JoyeLibert::keygen (RandGen &randgen) const
{
  return SecretKey (*this, randgen);
}

/* */
inline
JoyeLibert::PublicKey::PublicKey (const SecretKey &sk) : y_(sk.y())
{
  Mpz::mul (N_, sk.p(), sk.q());
}

/* */
inline
const Mpz & JoyeLibert::PublicKey::N () const
{
  return N_;
}

/* */
inline
const Mpz & JoyeLibert::PublicKey::y () const
{
  return y_;
}

/* */
inline
JoyeLibert::PublicKey JoyeLibert::keygen (const SecretKey &sk) const
{
  return PublicKey (sk);
}

/* */
inline
JoyeLibert::CipherText::CipherText (const JoyeLibert &C, const PublicKey &pk,
                                    const ClearText &m, const Mpz &r)
{
  Mpz ym;
  Mpz x (r);

  for (size_t i = 0 ; i < C.k_; i++)
  {
    Mpz::mul (x, x, x);               /* x^2 */
    Mpz::mod (x, x, pk.N());            /* x^2 mod N */
  }

  Mpz::pow_mod (ym, pk.y(), m, pk.N());   /* y^m mod N */
  Mpz::mul (*this, ym, x);            /* y^m * x^(2^k) */
  Mpz::mod (*this, *this, pk.N());      /* y^m * x^(2^k) mod N */
}

/* */
inline
JoyeLibert::CipherText JoyeLibert::encrypt (const PublicKey &pk,
                                            const ClearText &m,
                                            const Mpz &r) const
{
  return CipherText (*this, pk, m, r);
}

/*
 * Same as above but without the randomness
 */
inline
JoyeLibert::CipherText JoyeLibert::encrypt (const PublicKey &pk,
                                            const ClearText &m,
                                            RandGen &randgen) const
{
  return encrypt (pk, m, randgen.random_mpz (encrypt_randomness_bound(pk)));
}


/* */
inline
JoyeLibert::ClearText::ClearText (const JoyeLibert &C, const Mpz &v) : Mpz (v)
{
  if (!(v.sgn() >= 0 && v < C.cleartext_bound()))
    throw std::range_error ("Cleartext is negative or too large");
}

/* */
inline
JoyeLibert::ClearText::ClearText (const JoyeLibert &C, RandGen &randgen)
    : Mpz (randgen.random_mpz (C.cleartext_bound()))
{
}

/* */
inline
JoyeLibert::ClearText::ClearText (const JoyeLibert &C, const SecretKey &sk,
                                  const CipherText &c) : Mpz (0UL)
{
  Mpz T, z, t;
  Mpz d(sk.D());
  Mpz B(1UL);

  Mpz::pow_mod (T, c, sk.r(), sk.p());      /* c^((p-1)/2^k) mod p */

  for (size_t j = 1 ; j < C.k_ ; j++)
  {
    t = 1UL;
    Mpz::mulby2k (t, t, C.k_ - j);      /* 2^(k-j) */
    Mpz::pow_mod (z, T, t, sk.p());       /* c^(2^(k-j)) mod p */
    if (z != 1UL)
    {
      Mpz::add (*this, *this, B);
      Mpz::mul (T, T, d);
      Mpz::mod (T, T, sk.p());
    }
    Mpz::add (B, B, B);
    Mpz::mul (d, d, d);
    Mpz::mod (d, d, sk.p());
  }
  if (T != 1UL)
  {
    Mpz::add (*this, *this, B);
  }
}

/* */
inline
JoyeLibert::ClearText JoyeLibert::decrypt (const SecretKey &sk,
                                           const CipherText &c) const
{
  return ClearText (*this, sk, c);
}

const Mpz & JoyeLibert::encrypt_randomness_bound (const PublicKey & pk) const
{
  return pk.N();
}

const Mpz & JoyeLibert::cleartext_bound () const
{
  return twotothek_;
}
#endif /* BICYCL_JOYE_LIBERT_INL */
