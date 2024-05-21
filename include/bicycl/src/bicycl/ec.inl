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
#ifndef BICYCL_EC_INL
#define BICYCL_EC_INL

/******************************************************************************/
/* */
inline
ECDSA::SecretKey::SecretKey (const ECDSA &C)
  : d_(C.random_mod_order()),
    Q_(C, d_)
{
}

/* */
inline
const OpenSSL::BN & ECDSA::SecretKey::d () const
{
  return d_;
}

/* */
inline
const OpenSSL::ECPoint & ECDSA::SecretKey::Q () const
{
  return Q_;
}

/* */
inline
ECDSA::ECDSA (SecLevel seclevel) : ECGroup(seclevel), H_(seclevel)
{
}

/* */
inline
ECDSA::SecretKey ECDSA::keygen () const
{
  return SecretKey (*this);
}

/* */
inline
ECDSA::PublicKey ECDSA::keygen (const SecretKey &sk) const
{
  return PublicKey (*this, sk.Q());
}

/* */
inline
OpenSSL::BN ECDSA::hash_message (const Message &m) const
{
  return OpenSSL::BN (H_ (m));
}

/* */
inline
ECDSA::Signature ECDSA::sign (const SecretKey &sk, const Message &m) const
{
  return Signature (*this, sk, m);
}

/* */
inline
ECDSA::Signature::Signature (const ECDSA &C, const SecretKey &sk,
                                             const Message &m)
{
  OpenSSL::BN z, tmp;

  z = C.hash_message (m);

  do
  {
    OpenSSL::BN k (C.random_mod_order());
    if (k.is_zero())
      continue;

    OpenSSL::ECPoint K (C, k);
    C.x_coord_of_point (tmp, K);
    C.mod_order (r_, tmp); /* r = x([k] P) mod n */
    if (r_.is_zero())
      continue;

    C.mul_mod_order (s_, r_, sk.d());

    OpenSSL::BN::add (s_, s_, z);

    C.inverse_mod_order (tmp, k);
    C.mul_mod_order (s_, s_, tmp); /* s = k^(-1)*(z + r*d) */
  } while (s_.is_zero());
}

/* */
inline
bool ECDSA::Signature::verify (const ECDSA &C, const PublicKey &Q,
                                               const Message &m) const
{
  OpenSSL::BN z, sinv, u1, u2, x1, tmp;

  if (!C.has_correct_prime_order (Q)) /* check that Q as order n */
    return false;

  if (!C.is_positive_less_than_order (r_))
    return false;

  if (!C.is_positive_less_than_order (s_))
    return false;

  bool ok = true;
  OpenSSL::ECPoint T (C);
  z = C.hash_message (m);
  C.inverse_mod_order (sinv, s_);
  C.mul_mod_order (u1, sinv, z);
  C.mul_mod_order (u2, sinv, r_);

  C.scal_mul (T, u1, u2, Q); /* u1*G + u2*Q */

  if (C.is_at_infinity (T))
    ok = false;
  else
  {
    C.x_coord_of_point (tmp, T);
    C.mod_order (x1, tmp);

    ok = (x1 == r_);
  }

  return ok;
}

/* */
inline
bool ECDSA::verif (const Signature &signature, const PublicKey &Q,
                                               const Message &m) const
{
  return signature.verify (*this, Q, m);
}

/* random message of random length between 4 and UCHAR_MAX */
inline
ECDSA::Message ECDSA::random_message ()
{
  unsigned char size;
  OpenSSL::random_bytes (&size, 1 * sizeof (unsigned char));
  size = (size < 4) ? 4 : size;
  Message m (size);
  OpenSSL::random_bytes (m.data(), m.size() * sizeof (unsigned char));
  return m;
}

/******************************************************************************/
/* */
inline
ECNIZKProof::ECNIZKProof (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                          const SecretValue &s, const PublicValue &Q)
  : R_(E)
{
  OpenSSL::BN r (E.random_mod_order());
  E.scal_mul_gen (R_, r);

  OpenSSL::BN c (hash_for_challenge (H, E, R_, Q)); /* c = Hash (E, R_, Q) */

  E.mul_mod_order (z_, c, s);
  E.sub_mod_order (z_, r, z_); /* z = r - c*s */
}

/* */
inline
ECNIZKProof::ECNIZKProof (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                          const SecretValue &s)
  : ECNIZKProof (E, H, s, compute_Q_from_secret (E, s))
{
}

/* */
inline
ECNIZKProof::ECNIZKProof (const OpenSSL::ECGroup &E, const ECNIZKProof &p)
  : R_ (E, p.R_), z_(p.z_)
{
}

/* */
inline
bool ECNIZKProof::verify (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                          const PublicValue &Q) const
{
  if (!E.is_in_group (Q))
    return false;

  OpenSSL::BN c (hash_for_challenge (H, E, R_, Q)); /* c = Hash (E, R_, Q) */

  OpenSSL::ECPoint rhs (E);
  E.scal_mul (rhs, z_, c, Q); /* z*P + cQ */

  return E.ec_point_eq (R_, rhs);
}

/* */
inline
OpenSSL::BN ECNIZKProof::hash_for_challenge (OpenSSL::HashAlgo & H,
                                             const OpenSSL::ECGroup &E,
                                             const OpenSSL::ECPoint &R,
                                             const OpenSSL::ECPoint &Q)
{
  return OpenSSL::BN (H (E, OpenSSL::ECPointGroupCRefPair (R, E),
                            OpenSSL::ECPointGroupCRefPair (Q, E)));
}

/* */
inline
OpenSSL::ECPoint ECNIZKProof::compute_Q_from_secret (const OpenSSL::ECGroup &E,
                                                     const SecretValue &s)
{
  OpenSSL::ECPoint Q (E);
  E.scal_mul_gen (Q, s);
  return Q;
}

/******************************************************************************/
/* */
inline
ECNIZKAoK::ECNIZKAoK (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                      const OpenSSL::ECPoint & R, const SecretValue &x,
                      const SecretValue &y, const SecretValue &rho,
                      const PublicValue &V, const PublicValue &A)
  : H_(E)
{
  OpenSSL::BN v (E.random_mod_order());
  OpenSSL::BN u (E.random_mod_order());

  E.scal_mul (H_, v, u, R); /* H = u R + v P */

  /* c = Hash(E, R, V, A, H) */
  OpenSSL::BN c (hash_for_challenge (H, E, R, V, A, H_));

  E.mul_mod_order (t1_, c, x);
  E.add_mod_order (t1_, t1_, u); /* t1 = u + c * x */

  E.mul_mod_order (t2_, c, y);
  E.mul_mod_order (u, c, c); /* use u as temp var */
  E.mul_mod_order (u, u, rho); /* use u as temp var */
  E.add_mod_order (t2_, t2_, u);
  E.add_mod_order (t2_, t2_, v); /* t2 = v + c*y + c^2 * rho */
}

/* */
inline
ECNIZKAoK::ECNIZKAoK (const OpenSSL::ECGroup &E, const ECNIZKAoK &p)
  : H_ (E, p.H_), t1_(p.t1_), t2_(p.t2_)
{
}

/* */
inline
bool ECNIZKAoK::verify (const OpenSSL::ECGroup &E, OpenSSL::HashAlgo &H,
                        const OpenSSL::ECPoint & R, const PublicValue &V,
                        const PublicValue &A) const
{
  /* c = Hash(E, R, V, A, H) */
  OpenSSL::BN c (hash_for_challenge (H, E, R, V, A, H_));

  OpenSSL::ECPoint lhs (E);
  OpenSSL::ECPoint rhs (E);

  E.scal_mul (lhs, t2_, t1_, R); /* t1 R + t2 P */

  E.scal_mul (rhs, c, A);
  E.ec_add (rhs, V, rhs);
  E.scal_mul (rhs, c, rhs);
  E.ec_add (rhs, rhs, H_); /* c V + c^2 A + H */

  return E.ec_point_eq (lhs, rhs);
}

/* */
inline
OpenSSL::BN ECNIZKAoK::hash_for_challenge (OpenSSL::HashAlgo &Hash,
                                           const OpenSSL::ECGroup &E,
                                           const OpenSSL::ECPoint &R,
                                           const OpenSSL::ECPoint &V,
                                           const OpenSSL::ECPoint &A,
                                           const OpenSSL::ECPoint &H)
{
  return OpenSSL::BN (Hash (E, OpenSSL::ECPointGroupCRefPair (R, E),
                               OpenSSL::ECPointGroupCRefPair (V, E),
                               OpenSSL::ECPointGroupCRefPair (A, E),
                               OpenSSL::ECPointGroupCRefPair (H, E)));
}

#endif /* BICYCL_EC_INL */
