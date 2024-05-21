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
#ifndef BICYCL_CL_HSM_UTILS_INL
#define BICYCL_CL_HSM_UTILS_INL

/* */
template <class Cryptosystem>
inline
CL_HSM_SecretKey<Cryptosystem>::CL_HSM_SecretKey (const Cryptosystem &C,
                                                  const Mpz &v)
  : Mpz (v)
{
  if (!(v.sgn() >= 0 && v < C.secretkey_bound()))
    throw std::range_error ("Secret key is negative or too large");
}

/* */
template <class Cryptosystem>
inline
CL_HSM_SecretKey<Cryptosystem>::CL_HSM_SecretKey (const Cryptosystem &C,
                                                  RandGen &r)
  : Mpz (r.random_mpz (C.secretkey_bound()))
{
}

/* */
template <class Cryptosystem>
inline
CL_HSM_PublicKey<Cryptosystem>::CL_HSM_PublicKey (const Cryptosystem &C,
                                      const CL_HSM_SecretKey<Cryptosystem> &sk)
{
  C.power_of_h (pk_, sk);

  d_ = (C.encrypt_randomness_bound().nbits () + 1)/2;
  e_ = d_/2 + 1;

  pk_de_precomp_ = pk_;
  for (size_t i = 0; i < d_+e_; i++)
  {
    if (i == e_)
      pk_e_precomp_ = pk_de_precomp_;
    if (i == d_)
      pk_d_precomp_ = pk_de_precomp_;
    C.Cl_G().nudupl (pk_de_precomp_, pk_de_precomp_);
  }
}

/* */
template <class Cryptosystem>
inline
CL_HSM_PublicKey<Cryptosystem>::CL_HSM_PublicKey (const Cryptosystem &C,
                                      const QFI &pk)
{
  pk_ = pk;

  d_ = (C.encrypt_randomness_bound().nbits () + 1)/2;
  e_ = d_/2 + 1;

  pk_de_precomp_ = pk_;
  for (size_t i = 0; i < d_+e_; i++)
  {
    if (i == e_)
      pk_e_precomp_ = pk_de_precomp_;
    if (i == d_)
      pk_d_precomp_ = pk_de_precomp_;
    C.Cl_G().nudupl (pk_de_precomp_, pk_de_precomp_);
  }
}

template <class Cryptosystem>
inline
const QFI & CL_HSM_PublicKey<Cryptosystem>::elt () const
{
  return pk_;
}

/* */
template <class Cryptosystem>
inline
void CL_HSM_PublicKey<Cryptosystem>::exponentiation (const Cryptosystem &C,
                                                     QFI &r, const Mpz &n) const
{
  C.Cl_G().nupow (r, pk_, n, d_, e_, pk_e_precomp_, pk_d_precomp_,
                                                    pk_de_precomp_);
}

/* */
template <class Cryptosystem>
inline
std::ostream & operator<< (std::ostream &o,
                           const CL_HSM_PublicKey<Cryptosystem> &pk)
{
  return o << pk.pk_;
}

/* */
template <class Cryptosystem>
inline
CL_HSM_ClearText<Cryptosystem>::CL_HSM_ClearText (const Cryptosystem &C,
                                                  const Mpz &v)
  : Mpz (v)
{
  if (!(v.sgn() >= 0 && v < C.cleartext_bound()))
    throw std::range_error ("Cleartext is negative or too large");
}

/* */
template <class Cryptosystem>
inline
CL_HSM_ClearText<Cryptosystem>::CL_HSM_ClearText (const Cryptosystem &C,
                                                  RandGen &r)
  : Mpz (r.random_mpz (C.cleartext_bound()))
{
}

/* */
template <class Cryptosystem>
inline
CL_HSM_ClearText<Cryptosystem>::CL_HSM_ClearText (const Cryptosystem &C,
                                      const CL_HSM_SecretKey<Cryptosystem> &sk,
                                      const CL_HSM_CipherText<Cryptosystem> &c)
{
  QFI fm;

  C.Cl_G().nupow (fm, c.c1(), sk);
  if (C.compact_variant())
    C.from_Cl_DeltaK_to_Cl_Delta (fm);

  C.Cl_Delta().nucompinv (fm, c.c2(), fm); /* c2/c1^sk */

  Mpz::operator= (C.dlog_in_F(fm));
}

/* */
template <class Cryptosystem>
inline
CL_HSM_ClearText<Cryptosystem>::CL_HSM_ClearText (const Cryptosystem &C,
                                                  const CL_HSM_ClearText &ma,
                                                  const CL_HSM_ClearText &mb)
{
  Mpz::add (*this, ma, mb);
  Mpz::mod (*this, *this, C.cleartext_bound());
}

/* */
template <class Cryptosystem>
inline
CL_HSM_ClearText<Cryptosystem>::CL_HSM_ClearText (const Cryptosystem &C,
                                                  const CL_HSM_ClearText &m,
                                                  const Mpz &s)
{
  Mpz::mul (*this, m, s);
  Mpz::mod (*this, *this, C.cleartext_bound());
}

/* */
template <class Cryptosystem>
inline
CL_HSM_CipherText<Cryptosystem>::CL_HSM_CipherText (const Cryptosystem &C,
                                      const CL_HSM_PublicKey<Cryptosystem> &pk,
                                      const CL_HSM_ClearText<Cryptosystem> &m,
                                      const Mpz &r)
{
#ifdef BICYCL_WITH_PTHREADS
  using FctType = void (Cryptosystem::*) (QFI &, const Mpz &) const;
  auto fctptr = static_cast<FctType>(&Cryptosystem::power_of_h);
  std::thread th (fctptr, C, std::ref(c1_), std::cref(r)); /* c1 = h^r */
#else
  C.power_of_h (c1_, r); /* c1 = h^r */
#endif

  QFI fm = C.power_of_f (m); /* fm = [q^2, q, ...]^m */
  pk.exponentiation (C, c2_, r); /* pk^r */

  if (C.compact_variant())
    C.from_Cl_DeltaK_to_Cl_Delta (c2_);
  C.Cl_Delta().nucomp (c2_, c2_, fm); /* c2 = f^m*pk^r */

#ifdef BICYCL_WITH_PTHREADS
  th.join ();
#endif
}

/* */
template <class Cryptosystem>
inline
CL_HSM_CipherText<Cryptosystem>::CL_HSM_CipherText (const Cryptosystem &C,
                                      const CL_HSM_PublicKey<Cryptosystem> &pk,
                                      const CL_HSM_CipherText &ca,
                                      const CL_HSM_CipherText &cb, const Mpz &r)
{
#ifdef BICYCL_WITH_PTHREADS
  using FctType = void (Cryptosystem::*) (QFI &, const Mpz &) const;
  auto fctptr = static_cast<FctType>(&Cryptosystem::power_of_h);
  std::thread th (fctptr, C, std::ref(c1_), std::cref(r)); /* c1 = h^r */
#else
  C.power_of_h (c1_, r); /* c1 = h^r */
#endif

  pk.exponentiation (C, c2_, r); /* pk^r */

  if (C.compact_variant())
    C.from_Cl_DeltaK_to_Cl_Delta (c2_);
  C.Cl_Delta().nucomp (c2_, c2_, ca.c2_);
  C.Cl_Delta().nucomp (c2_, c2_, cb.c2_);

#ifdef BICYCL_WITH_PTHREADS
  th.join ();
#endif

  C.Cl_G().nucomp (c1_, c1_, ca.c1_);
  C.Cl_G().nucomp (c1_, c1_, cb.c1_);
}

/* */
template <class Cryptosystem>
inline
CL_HSM_CipherText<Cryptosystem>::CL_HSM_CipherText (const Cryptosystem &C,
                                      const CL_HSM_PublicKey<Cryptosystem> &pk,
                                      const CL_HSM_CipherText &c, const Mpz &s,
                                      const Mpz &r)
{
  QFI tmp;

#ifdef BICYCL_WITH_PTHREADS
  using FctType = void (Cryptosystem::*) (QFI &, const Mpz &) const;
  auto fctptr = static_cast<FctType>(&Cryptosystem::power_of_h);
  std::thread th (fctptr, C, std::ref(c1_), std::cref(r)); /* c1 = h^r */
#else
  C.power_of_h (c1_, r); /* c1 = h^r */
#endif

  pk.exponentiation (C, c2_, r); /* pk^r */
  if (C.compact_variant())
    C.from_Cl_DeltaK_to_Cl_Delta (c2_);
  C.Cl_Delta().nupow (tmp, c.c2_, s);
  C.Cl_Delta().nucomp (c2_, c2_, tmp);

#ifdef BICYCL_WITH_PTHREADS
  th.join ();
#endif

  C.Cl_G().nupow (tmp, c.c1_, s);
  C.Cl_G().nucomp (c1_, c1_, tmp);
}

/* */
template <class Cryptosystem>
inline
const QFI & CL_HSM_CipherText<Cryptosystem>::c1 () const
{
  return c1_;
}

/* */
template <class Cryptosystem>
inline
const QFI & CL_HSM_CipherText<Cryptosystem>::c2 () const
{
  return c2_;
}

#endif /* BICYCL_CL_HSM_UTILS_INL */
