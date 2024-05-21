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
#ifndef BICYCL_OPENSSL_WRAPPER_INL
#define BICYCL_OPENSSL_WRAPPER_INL

/******************************************************************************/
inline
void random_bytes (unsigned char *buf, int num)
{
  int ret = RAND_bytes (buf, num);
  if (ret != 1)
    throw std::runtime_error ("RAND_bytes failed in random_bytes");
}

/******************************************************************************/
/* */
inline
EVP_MD_CTX * HashAlgo::new_ctx_ ()
{
  EVP_MD_CTX *r = EVP_MD_CTX_new ();
  if (r == NULL)
    throw std::runtime_error ("EVP_MD_CTX_new failed in HashAlgo");
  return r;
}

/* */
inline
HashAlgo::HashAlgo (int nid) : md_(EVP_get_digestbynid (nid)),
                               mdctx_ (new_ctx_())
{
  if (md_ == NULL)
    throw std::runtime_error ("could not set EVP from nid in HashAlgo");
}

/* */
inline
HashAlgo::HashAlgo (SecLevel seclevel) : HashAlgo (seclevel.sha3_openssl_nid())
{
}

/* */
inline
HashAlgo::HashAlgo (const HashAlgo &H) : md_ (H.md_), mdctx_ (new_ctx_())
{
  operator= (H);
}

/* */
inline
HashAlgo::HashAlgo (HashAlgo &&H) : md_ (H.md_), mdctx_ (H.mdctx_)
{
  H.mdctx_ = NULL;
}

/* */
inline
HashAlgo::~HashAlgo ()
{
  EVP_MD_CTX_free (mdctx_);
}

/* */
inline
HashAlgo & HashAlgo::operator= (const HashAlgo &H)
{
  md_ = H.md_;
  int ret = EVP_MD_CTX_copy_ex (mdctx_, H.mdctx_);
  if (ret != 1)
    throw std::runtime_error ("could not copy EVP_MD_CTX");
  return *this;
}

/* */
inline
HashAlgo & HashAlgo::operator= (HashAlgo &&H)
{
  md_ = H.md_;
  mdctx_ = H.mdctx_;
  H.mdctx_ = NULL;
  return *this;
}

/* */
template <typename... Args>
inline
HashAlgo::Digest HashAlgo::operator() (const Args & ...args)
{
  int ret = EVP_DigestInit_ex (mdctx_, md_, NULL);
  if (ret != 1)
    throw std::runtime_error ("EVP_DigestInit_ex failed in HashAlgo");

  hash_update (args...);

  Digest h (digest_nbytes ());
  ret = EVP_DigestFinal_ex (mdctx_, h.data(), NULL);
  if (ret != 1)
    throw std::runtime_error ("EVP_DigestFinal_ex failed in HashAlgo");
  return h;
}

/* */
template <typename First, typename... Args>
inline
void HashAlgo::hash_update (const First & first, const Args & ...args)
{
  hash (first);
  hash_update (args...);
}

/* */
inline
void HashAlgo::hash_update ()
{
}

/* */
inline
void HashAlgo::hash_bytes (const void *ptr, size_t n)
{
  int ret = EVP_DigestUpdate (mdctx_, ptr, n);
  if (ret != 1)
    throw std::runtime_error ("EVP_DigestUpdate failed in hash_bytes");
}

/* */
template <>
inline
void HashAlgo::hash (const int &v)
{
  hash_bytes (&v, sizeof (v));
}

template <>
inline
void HashAlgo::hash (const unsigned int &v)
{
  hash_bytes (&v, sizeof (v));
}

/* */
template <>
inline
void HashAlgo::hash (const std::vector<unsigned char> &m)
{
  hash_bytes (m.data(), m.size() * sizeof(unsigned char));
}

/* */
template <>
inline
void HashAlgo::hash (const Mpz &v)
{
  mpz_srcptr vptr = static_cast<mpz_srcptr> (v);
  hash (v.sgn());
  hash_bytes (mpz_limbs_read (vptr), v.nlimbs()*sizeof(mp_limb_t));
}

/* */
inline
int HashAlgo::digest_nbytes () const
{
  return EVP_MD_size (md_);
}

/* */
inline
int HashAlgo::digest_nbits () const
{
  return 8 * digest_nbytes();
}

/******************************************************************************/
/* */
inline
BN::BN () : bn_(BN_new())
{
  if (bn_ == NULL)
    throw std::runtime_error ("could not allocate BIGNUM");
}

/* */
inline
BN::BN (const BN &other) : bn_ (BN_dup (other.bn_))
{
  if (bn_ == NULL)
    throw std::runtime_error ("could not duplicate BIGNUM");
}

/* */
inline
BN::BN (BN &&other) : bn_(other.bn_)
{
  other.bn_ = NULL;
}

/* */
inline
BN::BN (const std::vector<unsigned char> &bytes)
  : bn_ (BN_bin2bn (bytes.data(), bytes.size(), NULL))
{
  if (bn_ == NULL)
    throw std::runtime_error ("Could not set BIGNUM from binary");
}

/* */
inline
BN::BN (const Mpz &v) : BN()
{
  *this = v;
}

/* */
inline
BN::BN (unsigned long v) : BN()
{
  *this = v;
}

/* */
inline
BN & BN::operator= (const BN &other)
{
  const BIGNUM *ret = BN_copy (bn_, other.bn_);
  if (ret == NULL)
    throw std::runtime_error ("could not copy BIGNUM");
  return *this;
}

/* */
inline
BN & BN::operator= (BN &&other)
{
  bn_ = other.bn_;
  other.bn_ = NULL;
  return *this;
}

/* */
inline
BN & BN::operator= (const Mpz &v)
{
  unsigned char *tmp = (unsigned char *) malloc ((v.nbits() + 7) / 8);
  if (tmp == NULL)
    throw std::runtime_error ("could not allocate temporary buffer");

  size_t len;
  mpz_export (tmp, &len, -1, 1, 0, 0, static_cast<mpz_srcptr>(v));
  BIGNUM *ret = BN_lebin2bn (tmp, len, bn_);
  if (ret == NULL)
    throw std::runtime_error ("BN_lebin2bn failed");
  free (tmp);
  if (v.sgn() < 0)
    neg();
  return *this;
}

/* */
inline
BN::~BN ()
{
  BN_free (bn_);
}

/* */
inline
BN::operator Mpz () const
{
  int nbytes = num_bytes();
  if (nbytes > 0)
  {
    Mpz r (BN::BIGNUM_abs_to_bytes (bn_));
    if (BN_is_negative (bn_))
      r.neg();
    return r;
  }
  else
    return Mpz(0UL);
}

/* */
inline
BN & BN::operator= (unsigned long other)
{
  BN_set_word (bn_, other);
  return *this;
}

/* */
inline
BN & BN::operator*= (unsigned long other)
{
  BN_mul_word (bn_, other);
  return *this;
}

/* */
inline
bool BN::operator== (const BN &other) const
{
  return BN_cmp (bn_, other.bn_) == 0;
}

/* */
inline
bool BN::operator!= (const BN &other) const
{
  return !(*this == other);
}

/* */
inline
bool BN::is_zero () const
{
  return BN_is_zero (bn_);
}

/* */
inline
int BN::num_bytes () const
{
  return BN_num_bytes (bn_);
}

/* */
inline
void BN::neg ()
{
  BN_set_negative (bn_, !BN_is_negative (bn_));
}

/* */
inline
void BN::add (BN &r, const BN &a, const BN &b)
{
  int ret = BN_add (r.bn_, a.bn_, b.bn_);
  if (ret != 1)
    throw std::runtime_error ("BN_add failed");
}

/* */
inline
std::ostream & operator<< (std::ostream &o, const BN &v)
{
  char *buf = BN_bn2dec (v.bn_);
  if (buf == NULL)
    throw std::runtime_error ("BN_bn2dec failed in operator<<");
  o << buf;
  OPENSSL_free (buf);
  return o;
}

/* */
inline
void random_BN_2exp (const BN &r, int nbits)
{
  int ret = BN_rand (r.bn_, nbits, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
  if (ret != 1)
  {
    throw std::runtime_error ("BN_rand failed in random_BN_2exp");
  }
}

/* */
inline
std::vector<unsigned char> BN::BIGNUM_abs_to_bytes (const BIGNUM *v)
{
  int nbytes = BN_num_bytes (v);
  std::vector<unsigned char> r (nbytes);
  int ret = BN_bn2binpad (v, r.data(), nbytes);
  if (ret < 0)
  {
    throw std::runtime_error ("BN_bn2binpad failed in BIGNUM_abs_to_bytes");
  }
  return r;
}

template <>
inline
void HashAlgo::hash (const BN &v)
{
  std::vector<unsigned char> b (BN::BIGNUM_abs_to_bytes (v.bn_));
  hash (BN_is_negative (v.bn_));
  hash (b);
}

/****************************************************************************/
/* */
inline
ECPoint::ECPoint (const ECGroup &E) : P_(NULL)
{
  P_ = EC_POINT_new (E.ec_group_);
  if (P_ == NULL)
    throw std::runtime_error ("EC_POINT_new failed in ECPoint constructor");
}

/* */
inline
ECPoint::ECPoint (const ECGroup &E, const ECPoint &Q)
  : P_(EC_POINT_dup (Q.P_, E.ec_group_))
{
  if (P_ == NULL)
    throw std::runtime_error ("EC_POINT_dup failed in ECPoint constructor");
}

/* */
inline
ECPoint::ECPoint (const ECGroup &E, const BN & m)
  : ECPoint (E)
{
  E.scal_mul_gen (*this, m);
}

/* */
inline
ECPoint::ECPoint (ECPoint &&Q) : P_(Q.P_)
{
  Q.P_ = NULL;
}

/*
 * Assumes Q can be copied into P_ (must be init with compatible ECGroup).
 */
inline
ECPoint & ECPoint::operator= (const ECPoint &Q)
{
  int ret = EC_POINT_copy (P_, Q.P_);
  if (ret != 1)
    throw std::runtime_error ("EC_POINT_copy failed in ECPoint::operator=");
  return *this;
}

/*
 * Assumes Q can be copied into P_ (must be init with compatible ECGroup).
 */
inline
ECPoint & ECPoint::operator= (ECPoint &&Q)
{
  P_ = Q.P_;
  Q.P_ = NULL;
  return *this;
}

/* */
inline
ECPoint::~ECPoint ()
{
  EC_POINT_free (P_);
}

/******************************************************************************/
/* */
inline
ECGroup::ECGroup (SecLevel seclevel) : ctx_ (BN_CTX_new())
{
  int nid = seclevel.elliptic_curve_openssl_nid(); /* openssl curve id */
  ec_group_ = EC_GROUP_new_by_curve_name (nid);
  if (ec_group_ == NULL)
    throw std::runtime_error ("could not allocate elliptic curve");

  if (ctx_ == NULL)
    throw std::runtime_error ("could not allocate BN_CTX");

  /* convert group order to Mpz */
  order_ = BN::BIGNUM_abs_to_bytes (get_order());
}

/* */
inline
ECGroup::ECGroup (ECGroup &&G) : ec_group_ (G.ec_group_),
                                 order_ (std::move(G.order_)),
                                 ctx_ (G.ctx_)
{
  G.ec_group_ = NULL;
  G.ctx_ = NULL;
}

/* */
inline
ECGroup::~ECGroup ()
{
  EC_GROUP_free (ec_group_);
  BN_CTX_free (ctx_);
}

/* */
inline
ECGroup & ECGroup::operator= (ECGroup &&G)
{
  ec_group_ = G.ec_group_;
  G.ec_group_ = NULL;
  ctx_ = G.ctx_;
  G.ctx_ = NULL;
  order_ = std::move (G.order_);
  return *this;
}

/* */
inline
const Mpz & ECGroup::order () const
{
  return order_;
}

/* */
inline
BN ECGroup::a () const
{
  BN a;
  int ret = EC_GROUP_get_curve (ec_group_, NULL, a.bn_, NULL, ctx_);
  if (ret != 1)
  {
    throw std::runtime_error ("EC_GROUP_get_curve failed for a");
  }
  return a;
}

/* */
inline
BN ECGroup::b () const
{
  BN b;
  int ret = EC_GROUP_get_curve (ec_group_, NULL, NULL, b.bn_, ctx_);
  if (ret != 1)
  {
    throw std::runtime_error ("EC_GROUP_get_curve failed for b");
  }
  return b;
}

/* */
inline
BN ECGroup::p () const
{
  BN p;
  int ret = EC_GROUP_get_curve (ec_group_, p.bn_, NULL, NULL, ctx_);
  if (ret != 1)
  {
    throw std::runtime_error ("EC_GROUP_get_curve failed for p");
  }
  return p;
}

/* */
inline
bool ECGroup::is_on_curve (const ECPoint &P) const
{
  return EC_POINT_is_on_curve (ec_group_, P.P_, ctx_);
}

/* */
inline
bool ECGroup::is_in_group (const ECPoint &P) const
{
  if (!is_on_curve (P))
    return false;
  ECPoint T (*this);
  int ret = EC_POINT_mul (ec_group_, T.P_, NULL, P.P_,
                                     EC_GROUP_get0_cofactor (ec_group_), ctx_);
  if (ret != 1)
    throw std::runtime_error ("EC_POINT_mul failed in is_on_curve");
  return !is_at_infinity (T);
}

/* */
inline
bool ECGroup::is_at_infinity (const ECPoint &P) const
{
  return EC_POINT_is_at_infinity (ec_group_, P.P_);
}

/* */
inline
void ECGroup::x_coord_of_point (BN &x, const ECPoint &P) const
{
  int ret = EC_POINT_get_affine_coordinates (ec_group_, P.P_, x.bn_, NULL,
                                                                          ctx_);
  if (ret != 1)
    throw std::runtime_error ("Could not get x coordinate");
}

inline
void ECGroup::coords_of_point (BN &x, BN &y, const ECPoint &P) const
{
  coords_of_point (x, y, P.P_);
}

/* */
inline
bool ECGroup::ec_point_eq (const ECPoint &P, const ECPoint &Q) const
{
  return EC_POINT_cmp (ec_group_, P.P_, Q.P_, ctx_) == 0;
}

/* */
inline
void ECGroup::ec_add (ECPoint &R, const ECPoint &P, const ECPoint &Q) const
{
  int ret = EC_POINT_add (ec_group_, R.P_, P.P_, Q.P_, ctx_);
  if (ret != 1)
    throw std::runtime_error ("EC_POINT_add failed in add");
}

/* */
inline
void ECGroup::ec_neg (ECPoint &R) const
{
  int ret = EC_POINT_invert (ec_group_, R.P_, ctx_);
  if (ret != 1)
    throw std::runtime_error ("EC_POINT_invert failed in ec_neg");
}

/* */
inline
void ECGroup::scal_mul_gen (ECPoint &R, const BN &n) const
{
  int ret = EC_POINT_mul (ec_group_, R.P_, n.bn_, NULL, NULL, ctx_);
  if (ret != 1)
    throw std::runtime_error ("EC_POINT_mul failed in scal_mul_gen");
}

/* */
inline
void ECGroup::scal_mul (ECPoint &R, const BN &n, const ECPoint &P) const
{
  int ret = EC_POINT_mul (ec_group_, R.P_, NULL, P.P_, n.bn_, ctx_);
  if (ret != 1)
    throw std::runtime_error ("EC_POINT_mul failed in scal_mul");
}

/* */
inline
void ECGroup::scal_mul (ECPoint &R, const BN &m, const BN &n,
                                                     const ECPoint &P) const
{
  int ret = EC_POINT_mul (ec_group_, R.P_, m.bn_, P.P_, n.bn_, ctx_);
  if (ret != 1)
    throw std::runtime_error ("EC_POINT_mul failed in scal_mul");
}

/* We assume that the order is prime (which must be the case for NIST curves) */
inline
bool ECGroup::has_correct_prime_order (const ECPoint &G) const
{
  if (is_at_infinity (G))
    return false;

  if (!is_on_curve (G))
    return false;

  ECPoint T (*this);
  int ret = EC_POINT_mul (ec_group_, T.P_, NULL, G.P_, get_order (), ctx_);
  if (ret != 1)
    throw std::runtime_error ("EC_POINT_mul failed in has_correct_prime_order");
  return is_at_infinity (T);
}

/* */
inline
void ECGroup::mod_order (BN &r, const BN &a) const
{
  int ret = BN_nnmod (r.bn_, a.bn_, get_order (), ctx_);
  if (ret != 1)
    throw std::runtime_error ("BN_nnmod failed");
}

/* */
inline
void ECGroup::add_mod_order (BN &r, const BN &a, const BN &b) const
{
  int ret = BN_mod_add (r.bn_, a.bn_, b.bn_, get_order (), ctx_);
  if (ret != 1)
    throw std::runtime_error ("BN_mod_add failed");
}

/* */
inline
void ECGroup::sub_mod_order (BN &r, const BN &a, const BN &b) const
{
  int ret = BN_mod_sub (r.bn_, a.bn_, b.bn_, get_order (), ctx_);
  if (ret != 1)
    throw std::runtime_error ("BN_mod_sub failed");
}

/* */
inline
void ECGroup::mul_by_word_mod_order (BN &r, BN_ULONG w) const
{
  int ret = BN_mul_word (r.bn_, w);
  if (ret != 1)
    throw std::runtime_error ("BN_mod_word failed");
  mod_order (r, r);
}

/* */
inline
void ECGroup::mul_mod_order (BN &r, const BN &a, const BN &b) const
{
  int ret = BN_mod_mul (r.bn_, a.bn_, b.bn_, get_order (), ctx_);
  if (ret != 1)
    throw std::runtime_error ("BN_mod_mul failed");
}

/* */
inline
void ECGroup::inverse_mod_order (BN &r, const BN &a) const
{
  const BIGNUM *ret = BN_mod_inverse (r.bn_, a.bn_, get_order (), ctx_);
  if (ret == NULL)
    throw std::runtime_error ("could not inverse modulo order");
}

/* */
inline
BN ECGroup::random_mod_order () const
{
  BN r;
  int ret = BN_priv_rand_range (r.bn_, get_order ());
  if (ret != 1)
    throw std::runtime_error ("BN_priv_rand_range failed");
  return r;
}

/* */
inline
bool ECGroup::is_positive_less_than_order (const BN &v) const
{
  return !BN_is_negative (v.bn_) && !BN_is_zero (v.bn_)
                                 && BN_cmp (v.bn_, get_order ()) < 0;
}

/* */
inline
const BIGNUM * ECGroup::get_order () const
{
  return EC_GROUP_get0_order (ec_group_);
}

/* */
inline
void ECGroup::coords_of_point (BN &x, BN &y, const EC_POINT *P) const
{
  int ret = EC_POINT_get_affine_coordinates (ec_group_, P, x.bn_, y.bn_, ctx_);
  if (ret != 1)
    throw std::runtime_error ("Could not get x, y coordinates");
}

/* */
inline
std::vector<unsigned char> ECGroup::EC_POINT_to_bytes (const EC_POINT *P) const
{
  size_t s = EC_POINT_point2oct (ec_group_, P, POINT_CONVERSION_UNCOMPRESSED,
                                 NULL, 0, ctx_);
  if (s == 0)
  {
    throw std::runtime_error ("could not get size in point_to_bytes");
  }
  std::vector<unsigned char> r(s);
  s = EC_POINT_point2oct (ec_group_, P, POINT_CONVERSION_UNCOMPRESSED,
                          r.data(), s, ctx_);
  if (s == 0)
  {
    throw std::runtime_error ("could not convert in point_to_bytes");
  }
  return r;
}

/* */
inline
std::vector<unsigned char> ECGroup::ECPoint_to_bytes (const ECPoint &P) const
{
  return EC_POINT_to_bytes (P.P_);
}

/* */
inline
std::ostream & operator<< (std::ostream &o, const ECGroup &E)
{
  BN x, y;
  E.coords_of_point (x, y, EC_GROUP_get0_generator (E.ec_group_));
  return o << "E = EllipticCurve(GF(" << E.p() << "), ["
           << E.a() << "," << E.b() << "])" << std::endl
           << "P = E(" << x << ", " << y << ")" << std::endl;
}

/* */
inline
std::ostream & operator<< (std::ostream &o, const ECPointGroupCRefPair &pt)
{
  const ECPoint & P = std::get<0> (pt);
  const ECGroup & E = std::get<1> (pt);
  BN x, y;
  E.coords_of_point (x, y, P);
  return o << "(" << x << ", " << y << ")";
}

/* */
template<>
inline
void HashAlgo::hash (const ECPointGroupCRefPair &pt)
{
  const ECPoint & P = std::get<0> (pt);
  const ECGroup & E = std::get<1> (pt);
  std::vector<unsigned char> bin (E.ECPoint_to_bytes (P));
  hash (bin);
}

/* */
template<>
inline
void HashAlgo::hash (const ECGroup &E)
{
  const EC_POINT *gen = EC_GROUP_get0_generator (E.ec_group_);
  std::vector<unsigned char> gen_bytes (E.EC_POINT_to_bytes (gen));
  hash (E.a());
  hash (E.b());
  hash (E.p());
  hash (gen_bytes);
}

#endif /* BICYCL_OPENSSL_WRAPPER_INL */
