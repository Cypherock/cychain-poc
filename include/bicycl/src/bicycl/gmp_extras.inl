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
#ifndef BICYCL_GMP_EXTRAS_INL
#define BICYCL_GMP_EXTRAS_INL

/*
 * GMP internal macros copied from gmp-impl.h.
 * Should be undef at the end of header
 */
#define ALLOC(x) ((x)->_mp_alloc)
#define PTR(x) ((x)->_mp_d)
#define SIZ(x) ((x)->_mp_size)
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define ABSIZ(x) ABS (SIZ (x))

#define GMP_NUMB_HIGHBIT  (((mp_limb_t) 1) << (GMP_NUMB_BITS-1))

#define MPN_EXTRACT_NUMB(count, xh, xl)                       \
    ((((xh) << ((count) - GMP_NAIL_BITS)) & GMP_NUMB_MASK) |  \
    ((xl) >> (GMP_LIMB_BITS - (count))))

#define MPN_NORMALIZE(DST, NLIMBS) do {   \
    while ((NLIMBS) > 0)                  \
    {                                     \
      if ((DST)[(NLIMBS) - 1] != 0)       \
        break;                            \
      (NLIMBS)--;                         \
    }                                     \
  } while (0)


/*
 * Should be undef at the end of header
 * Note: ideally we want to get the same definition as in longlong.h from GMP,
 * but there is a lot of cases depending on architectures and compilers. For
 * example, for gcc, it uses __builtin_clzll. To be generic, we need to use
 * mpz_sizeinbase.
 * Warning: x must be a variable, not a constant, as we take its address.
 */
//#define count_leading_zeros(count, x) (count) = __builtin_clzll((x))
#define count_leading_zeros(count, x) (count) = GMP_LIMB_BITS - mpn_sizeinbase (&(x), 1, 2)

#define mpn_hgcd2 __MPN(hgcd2)
#define mpn_matrix22_mul1_inverse_vector __MPN(matrix22_mul1_inverse_vector)
#define mpn_hgcd_mul_matrix1_vector __MPN (hgcd_mul_matrix1_vector)
#define mpn_binvert __MPN (binvert)
#define mpn_binvert_itch __MPN (binvert_itch)
#define mpn_redc_n __MPN (redc_n)

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
namespace
{
  /*
   * GMP internal struct copied from gmp-impl.h.
   */
  struct hgcd_matrix1
  {
    mp_limb_t u[2][2];
  };

  /*
   * GMP internal functions, defined in gmp-impl.h, implemented in mpn/ files,
   * exported in libgmp but not defined in gmp.h.
   */
  extern "C"
  {
    extern
    int mpn_hgcd2 (mp_limb_t, mp_limb_t, mp_limb_t, mp_limb_t,
                    struct hgcd_matrix1 *);
    extern
    mp_size_t mpn_matrix22_mul1_inverse_vector (const struct hgcd_matrix1 *,
                                                mp_ptr, mp_srcptr, mp_ptr,
                                                mp_size_t);
    extern
    mp_size_t mpn_hgcd_mul_matrix1_vector (const struct hgcd_matrix1 *, mp_ptr,
                                            mp_srcptr, mp_ptr, mp_size_t);
    extern
    void mpn_binvert (mp_ptr rp, mp_srcptr up, mp_size_t n, mp_ptr scratch);
    extern
    mp_size_t mpn_binvert_itch (mp_size_t n);
    extern
    void mpn_redc_n (mp_ptr rp, mp_ptr up, mp_srcptr mp, mp_size_t n,
                                                         mp_srcptr ip);
  }

  /*
   * Return the 2*GMP_NUMB_BITS most significant bits of ap and bp.
   * It assumes that one of ap[n-1] or bp[n-1] is nonzero.
   * If l = max (log2(ap), log2(bp)), then
   *    ah*2^GMP_NUMB_BITS + al = floor(ap/2^(l-2*GMP_NUMB_BITS))
   *    bh*2^GMP_NUMB_BITS + bl = floor(bp/2^(l-2*GMP_NUMB_BITS))
   *  with one of the most significant bit of ah or bh equals to 1.
   */
  static inline void
  mpn_highest_two_limbs (mp_limb_t *ah, mp_limb_t *al, mp_limb_t *bh,
                         mp_limb_t *bl, const mp_ptr ap, const mp_ptr bp,
                         mp_size_t n)
  {
    mp_limb_t mask = ap[n-1] | bp[n-1]; /* nonzero by assumption */

    if (mask & GMP_NUMB_HIGHBIT) /* if we are lucky, no shift is necessary */
    {
      *ah = ap[n-1]; *al = ap[n-2];
      *bh = bp[n-1]; *bl = bp[n-2];
    }
    else if (n == 2)
    {
      /* We use the full inputs without truncation, so we can
        safely shift left. */
      int shift;

      count_leading_zeros (shift, mask);
      *ah = MPN_EXTRACT_NUMB (shift, ap[1], ap[0]);
      *al = ap[0] << shift;
      *bh = MPN_EXTRACT_NUMB (shift, bp[1], bp[0]);
      *bl = bp[0] << shift;
    }
    else
    {
      int shift;

      count_leading_zeros (shift, mask);
      *ah = MPN_EXTRACT_NUMB (shift, ap[n-1], ap[n-2]);
      *al = MPN_EXTRACT_NUMB (shift, ap[n-2], ap[n-3]);
      *bh = MPN_EXTRACT_NUMB (shift, bp[n-1], bp[n-2]);
      *bl = MPN_EXTRACT_NUMB (shift, bp[n-2], bp[n-3]);
    }
  }

  /*
   * Input: rp of length <= n, sp of length <= n and qp of length exactly qn.
   * Output: overwrite rp with rp+qp*sp and returns its length.
   *
   * Assumption: at least of of rp or sp must be of length n, and rp should be
   * large enough to store the result.
   *
   * Scratch variable(s): tq (must be of length at least qn + n)
   */
  static inline mp_size_t
  mpn_addmul (mp_ptr rp, mp_srcptr qp, mp_size_t qn, mp_srcptr sp,
              mp_size_t n, mp_ptr tp)
  {
    mp_limb_t cy;
    if (qn == 1) /* common case: q has only 1 limb */
    {
      mp_limb_t q = qp[0];
      if (q == 1)
        cy = mpn_add_n (rp, rp, sp, n);
      else
        cy = mpn_addmul_1 (rp, sp, n, q);
    }
    else
    {
      mp_size_t spn = n;
      MPN_NORMALIZE (sp, spn);
      if (spn > 0)
      {
        if (qn > spn)
          mpn_mul (tp, qp, qn, sp, spn);
        else
          mpn_mul (tp, sp, spn, qp, qn);
        mp_size_t tpn = spn + qn;
        tpn -= (tp[tpn-1] == 0);

        if (tpn >= n)
        {
          cy = mpn_add (rp, tp, tpn, rp, n);
          n = tpn;
        }
        else /* tpn < n */
        {
          /* In this case, rp has exactly n limbs before the addition because
           * qn >= 1 and tpn < n implies spn < n.
           */
          cy = mpn_add (rp, rp, n, tp, tpn);
        }
      }
      else /* sp is zero => no mul, no add */
        cy = 0;
    }
    rp[n] = cy;
    n += (cy > 0);
    return n;
  }

  /* */
  static inline void
  redcify (mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr mp, mp_size_t n,
           mp_ptr tp, mp_ptr qp)
  {
    mpn_zero (tp, n);
    mpn_copyi (tp + n, up, un);
    mpn_tdiv_qr (qp, rp, 0L, tp, un + n, mp, n);
  }

} /* anonymous namespace */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
class Mpz::ModInverseException : public std::invalid_argument
{
  public:
    ModInverseException() : std::invalid_argument("not invertible") { }
};

/* */
inline
Mpz::Mpz ()
{
  mpz_init (mpz_);
}

/* */
inline
Mpz::Mpz (const Mpz &v)
{
  mpz_init_set (mpz_, v.mpz_);
}

/* */
inline
Mpz::Mpz (Mpz &&v)
{
  ALLOC(mpz_) = ALLOC (v.mpz_);
  SIZ(mpz_) = SIZ (v.mpz_);
  PTR(mpz_) = PTR (v.mpz_);

  ALLOC(v.mpz_) = 0;
  SIZ(v.mpz_) = 0;
  PTR(v.mpz_) = NULL;
}

#if 0
/* */
inline
Mpz::Mpz (mpz_srcptr v)
{
  mpz_init_set (mpz_, v);
}
#endif

/* */
inline
Mpz::Mpz (unsigned long v)
{
  mpz_init_set_ui (mpz_, v);
}

/* */
inline
Mpz::Mpz (long v)
{
  mpz_init_set_si (mpz_, v);
}

/* */
inline
Mpz::Mpz (const std::string &s)
{
  int ret = mpz_init_set_str (mpz_, s.c_str(), 0);
  if (ret)
    throw std::runtime_error (std::string ("could not parse '") + s
                                        + std::string ("' as a valid number"));
}

/* */
inline
Mpz::Mpz (mpf_srcptr v) : Mpz()
{
  *this = v;
}

/* */
inline
Mpz::Mpz (const std::vector<unsigned char> &data) : Mpz()
{
  *this = data;
}

/* */
inline
Mpz::~Mpz ()
{
  mpz_clear (mpz_);
}

/* */
inline
Mpz & Mpz::operator= (const Mpz &v)
{
  mpz_set (mpz_, v.mpz_);
  return *this;
}

/* */
inline
Mpz & Mpz::operator= (Mpz &&v)
{
  mpz_swap (mpz_, v.mpz_);
  return *this;
}

/* */
inline
Mpz & Mpz::operator= (mpf_srcptr v)
{
  mpz_set_f (mpz_, v);
  return *this;
}

/* */
inline
Mpz & Mpz::operator= (unsigned long v)
{
  mpz_set_ui (mpz_, v);
  return *this;
}

/* */
inline
Mpz & Mpz::operator= (long v)
{
  mpz_set_si (mpz_, v);
  return *this;
}

/* */
inline
Mpz & Mpz::operator= (const std::string &s)
{
  int ret = mpz_set_str (mpz_, s.c_str(), 0);
  if (ret)
    throw std::runtime_error (std::string ("could not parse '") + s
                                        + std::string ("' as a valid number"));
  return *this;
}

/* */
inline
Mpz & Mpz::operator= (const std::vector<unsigned char> &data)
{
  /* the binary data is interpreted most significant bit first */
  mpz_import (mpz_, data.size(), 1, 1, 0, 0, data.data());
  return *this;
}

/* */
inline
bool Mpz::operator== (const Mpz &other) const
{
  return richcmp (*this, other) == 0;
}

/* */
inline
bool Mpz::operator!= (const Mpz &other) const
{
  return richcmp (*this, other) != 0;
}

/* */
inline
bool Mpz::operator<  (const Mpz &other) const
{
  return richcmp (*this, other) <  0;
}

/* */
inline
bool Mpz::operator>  (const Mpz &other) const
{
  return richcmp (*this, other) >  0;
}

/* */
inline
bool Mpz::operator<= (const Mpz &other) const
{
  return richcmp (*this, other) <= 0;
}

/* */
inline
bool Mpz::operator>= (const Mpz &other) const
{
  return richcmp (*this, other) >= 0;
}

/* */
inline
bool Mpz::operator== (unsigned long v) const
{
  return richcmp (*this, v) == 0;
}

/* */
inline
bool Mpz::operator!= (unsigned long v) const
{
  return richcmp (*this, v) != 0;
}

/* */
inline
bool Mpz::operator<  (unsigned long v) const
{
  return richcmp (*this, v) <  0;
}

/* */
inline
bool Mpz::operator>  (unsigned long v) const
{
  return richcmp (*this, v) >  0;
}

/* */
inline
bool Mpz::operator<= (unsigned long v) const
{
  return richcmp (*this, v) <= 0;
}

/* */
inline
bool Mpz::operator>= (unsigned long v) const
{
  return richcmp (*this, v) >= 0;
}

/* */
inline
bool Mpz::operator== (long v) const
{
  return richcmp (*this, v) == 0;
}

/* */
inline
bool Mpz::operator!= (long v) const
{
  return richcmp (*this, v) != 0;
}

/* */
inline
bool Mpz::operator<  (long v) const
{
  return richcmp (*this, v) <  0;
}

/* */
inline
bool Mpz::operator>  (long v) const
{
  return richcmp (*this, v) >  0;
}

/* */
inline
bool Mpz::operator<= (long v) const
{
  return richcmp (*this, v) <= 0;
}

/* */
inline
bool Mpz::operator>= (long v) const
{
  return richcmp (*this, v) >= 0;
}

/* */
inline
Mpz::operator mpz_srcptr() const
{
  return mpz_;
}

/* */
inline
Mpz::operator unsigned long () const
{
  if (!mpz_fits_uint_p (mpz_))
    throw std::runtime_error ("mpz value could not be parsed as an unsigned long");
  else
    return mpz_get_ui (mpz_);
}

/* */
inline
Mpz::operator long () const
{
  if (!mpz_fits_sint_p (mpz_))
    throw std::runtime_error ("mpz value could not be parsed as an long");
  else
    return mpz_get_si (mpz_);
}

/* */
inline
size_t Mpz::nbits () const
{
  return mpz_sizeinbase (mpz_, 2);
}

/* */
inline
size_t Mpz::ndigits () const
{
  return mpz_sizeinbase (mpz_, 10);
}

/* */
inline
size_t Mpz::nlimbs () const
{
  return mpz_size (mpz_);
}

/* */
inline
int Mpz::sgn () const
{
  return mpz_sgn (mpz_);
}

/* */
inline
bool Mpz::is_zero () const
{
  return sgn() == 0;
}

/* */
inline
bool Mpz::is_odd () const
{
  return mpz_odd_p (mpz_);
}

/* */
inline
bool Mpz::is_even () const
{
  return mpz_even_p (mpz_);
}

/* */
inline
bool Mpz::is_one () const
{
  return mpz_cmp_ui (mpz_, 1UL) == 0;
}

/*
 * Return 1 if n is prime (or probably prime), return 0 otherwise.
 */
inline
bool Mpz::is_prime (int reps) const
{
  int r = mpz_probab_prime_p (mpz_, reps);
  return (r > 0);
}

/* */
inline
bool Mpz::is_divisible_by (const Mpz &d) const
{
  return mpz_divisible_p (mpz_, d.mpz_);
}

/* */
inline
void Mpz::neg ()
{
  SIZ(mpz_) = -SIZ(mpz_);
}

/*
 * Return integer with bits from index to index-len+1 (inclusive)
 * Assumes len is positive and less than GMP_LIMB_BITS (i.e., 32 or 64)
 * The sign is not taken into account.
 */
inline
mp_limb_t Mpz::extract_bits (size_t index, size_t len) const
{
  mp_limb_t r;
  const mp_limb_t mask = ((1UL << len)-1UL);
  const size_t limb_high = mpz_getlimbn (mpz_, index/GMP_NUMB_BITS);
  const size_t im = index % GMP_NUMB_BITS;

  if (len <= im+1) /* all bits are in limb_high */
  {
    r = limb_high >> (im-len+1);
  }
  else
  {
    const size_t limb_low = len > index+1 ? 0 :
                          mpz_getlimbn (mpz_, (index - len + 1)/GMP_NUMB_BITS);
    r = (limb_high << (len-im-1)) | (limb_low >> (im+1+GMP_NUMB_BITS-len));
  }

  return r & mask;
}

/* */
inline
int Mpz::tstbit (size_t i) const
{
  return mpz_tstbit (mpz_, i);
}

/* */
inline
void Mpz::setbit (size_t i)
{
  mpz_setbit (mpz_, i);
}

/* */
inline
unsigned long Mpz::mod4 () const
{
  unsigned long r = mpz_get_ui (mpz_);
  return (sgn() < 0 ? -r : r) & 0x3UL;
}

/* */
inline
unsigned long Mpz::mod8 () const
{
  unsigned long r = mpz_get_ui (mpz_);
  return (sgn() < 0 ? -r : r) & 0x7UL;
}

/* */
inline
size_t Mpz::val2 () const
{
  return sgn() == 0 ? SIZE_MAX : mpn_scan1 (mpz_->_mp_d, 0);
}

/* */
inline
void Mpz::nextprime ()
{
  mpz_nextprime (mpz_, mpz_);
}

/* */
inline
int Mpz::legendre (const Mpz &l) const
{
  return mpz_legendre (mpz_, l.mpz_);
}

/* */
inline
int Mpz::jacobi (const Mpz &l) const
{
  return mpz_jacobi (mpz_, l.mpz_);
}

/* */
inline
int Mpz::kronecker (const Mpz &l) const
{
  return mpz_kronecker (mpz_, l.mpz_);
}

/* */
inline
void Mpz::swap (Mpz &op1, Mpz &op2)
{
  mpz_swap (op1.mpz_, op2.mpz_);
}

/* */
inline
void Mpz::abs (Mpz &r, const Mpz &op)
{
  mpz_abs (r.mpz_, op.mpz_);
}

/* */
inline
int Mpz::cmpabs (const Mpz &a, const Mpz &b)
{
  return mpz_cmpabs (a.mpz_, b.mpz_);
}

/* */
inline
void Mpz::add (Mpz &r, const Mpz &op1, const Mpz &op2)
{
  mpz_add (r.mpz_, op1.mpz_, op2.mpz_);
}

/* */
inline
void Mpz::add (Mpz &r, const Mpz &op1, unsigned long op2)
{
  mpz_add_ui (r.mpz_, op1.mpz_, op2);
}

/* */
inline
char * Mpz::get_str(const Mpz &op1, size_t base)
{
  return mpz_get_str (NULL, base, op1.mpz_);
}

/* */
inline
void Mpz::set_str(Mpz &op1, const char *str, size_t base)
{
  int ret = mpz_set_str (op1.mpz_, str, base);
  if (ret)
    throw std::runtime_error (std::string ("could not parse '") + str
                                        + std::string ("' as a valid number"));
}

/* */
inline
void Mpz::from_bytes(Mpz &op1, uint8_t *data, size_t size)
{
  mpz_import (op1.mpz_, size, 1, sizeof(data[0]), 0, 0, data);
}

/* */
inline
void Mpz::to_bytes(const Mpz &op1, uint8_t *data, size_t size)
{
  size_t count;
  mpz_export (data, &count, 1, sizeof(data[0]), 0, 0, op1.mpz_);
  if (count != size)
    throw std::runtime_error ("unexpected size of exported data");
}

/* */
inline
void Mpz::sub (Mpz &r, const Mpz &op1, const Mpz &op2)
{
  mpz_sub (r.mpz_, op1.mpz_, op2.mpz_);
}

/* */
inline
void Mpz::sub (Mpz &r, const Mpz &op1, unsigned long op2)
{
  mpz_sub_ui (r.mpz_, op1.mpz_, op2);
}

/* */
inline
void Mpz::mul (Mpz &r, const Mpz &op1, const Mpz &op2)
{
  mpz_mul (r.mpz_, op1.mpz_, op2.mpz_);
}

/* */
inline
void Mpz::mul (Mpz &r, const Mpz &op1, unsigned long op2)
{
  mpz_mul_ui (r.mpz_, op1.mpz_, op2);
}

/* */
inline
void Mpz::mulby2k (Mpz &r, const Mpz &op, mp_bitcnt_t k)
{
  mpz_mul_2exp (r.mpz_, op.mpz_, k);
}

/* */
inline
void Mpz::mulby2k (Mpz &r, unsigned long op, mp_bitcnt_t k)
{
  r = op;
  mulby2k (r, r, k);
}

/* */
inline
void Mpz::mulby2 (Mpz &r, const Mpz &op)
{
  mulby2k (r, op, 1); /* maybe an add is faster ?? */
}

/* */
inline
void Mpz::mulby4 (Mpz &r, const Mpz &op)
{
  mulby2k (r, op, 2);
}

/* */
inline
void Mpz::addmul (Mpz &r, const Mpz &op1, const Mpz &op2)
{
  mpz_addmul (r.mpz_, op1.mpz_, op2.mpz_);
}

/* */
inline
void Mpz::submul (Mpz &r, const Mpz &op1, const Mpz &op2)
{
  mpz_submul (r.mpz_, op1.mpz_, op2.mpz_);
}

/* */
inline
void Mpz::divby2k (Mpz &r, const Mpz &op, mp_bitcnt_t k)
{
  mpz_fdiv_q_2exp (r.mpz_, op.mpz_, k);
}

/* */
inline
void Mpz::divby2 (Mpz &r, const Mpz &op)
{
  divby2k (r, op, 1);
}

/* */
inline
void Mpz::divby4 (Mpz &r, const Mpz &op)
{
  divby2k (r, op, 2);
}

/* */
inline
void Mpz::divexact (Mpz &r, const Mpz &op1, const Mpz &op2)
{
  mpz_divexact (r.mpz_, op1.mpz_, op2.mpz_);
}

/* */
inline
void Mpz::divexact (Mpz &r, const Mpz &op1, unsigned long op2)
{
  mpz_divexact_ui (r.mpz_, op1.mpz_, op2);
}

/* */
inline
void Mpz::cdiv_qr (Mpz &q, Mpz &r, const Mpz &n, const Mpz &d)
{
  mpz_cdiv_qr (q.mpz_, r.mpz_, n.mpz_, d.mpz_);
}

/* */
inline
void Mpz::fdiv_qr (Mpz &q, Mpz &r, const Mpz &n, const Mpz &d)
{
  mpz_fdiv_qr (q.mpz_, r.mpz_, n.mpz_, d.mpz_);
}

/* */
inline
void Mpz::mod (Mpz &r, const Mpz &op1, const Mpz &op2)
{
  mpz_mod (r.mpz_, op1.mpz_, op2.mpz_);
}

/* */
inline
void Mpz::mod2k (Mpz &r, const Mpz &op, mp_bitcnt_t k)
{
  mpz_fdiv_r_2exp (r.mpz_, op.mpz_, k);
}

/* */
inline
void Mpz::mod2k_centered (Mpz &r, const Mpz &op, mp_bitcnt_t k)
{
  Mpz::mod2k (r, op, k);

  /* substract 2^k if needed */
  if (r.tstbit (k-1))
  {
    mpn_com (r.mpz_->_mp_d, r.mpz_->_mp_d, r.mpz_->_mp_size);
    Mpz::add (r, r, 1UL);
    Mpz::mod2k (r, r, k);
    r.neg();
  }
}

/* */
inline
void Mpz::mod_inverse (Mpz &r, const Mpz &op1, const Mpz &op2)
{
  int ret = mpz_invert (r.mpz_, op1.mpz_, op2.mpz_);
  if (ret == 0)
    throw ModInverseException ();
}

/*
 * Input: op1 and k
 *
 * Output: set r to op1^-1 modulo 2^k
 *
 * Assumption:
 *    - op1 is odd
 *    - k >= 1
 *    - r and op1 are different variables
 *
 * Scratch variable(s): t
 */
inline
void Mpz::mod_inverse_2k (Mpz &r, const Mpz &op1, mp_bitcnt_t k, Mpz &t)
{
  if (op1.is_even())
    throw ModInverseException ();
  r = 1UL;
  for (size_t i = 1; i < k; i <<= 1)
  {
    t = 2UL;
    Mpz::submul (t, r, op1);
    Mpz::mul (t, r, t);
    Mpz::mod2k (r, t, i << 1);
  }
}

/* */
inline
void Mpz::pow_ui (Mpz &r, const Mpz &b, unsigned long e)
{
  mpz_pow_ui (r.mpz_, b.mpz_, e);
}

/* */
inline
void Mpz::pow_mod (Mpz &r, const Mpz &b, const Mpz &e, const Mpz &m)
{
  mpz_powm (r.mpz_, b.mpz_, e.mpz_, m.mpz_);
}

/*
 * Multiple exponentiation.
 *
 * Input:
 *  - f: positive integer (the base)
 *  - n: integer (the exponent)
 *  - m: positive integer (the modulus)
 *  - e: a positive integer
 *  - fe: positive integer such that fe = f^(2^e) mod m
 *  - f2e: positive integer such that fe = f^(2^2e) mod m
 *  - f3e: positive integer such that fe = f^(2^3e) mod m
 * Ouput:
 *  - r: output integer corresponding to f^n mod m.
 *
 * Assumption: the modulus m is odd
 *
 */
inline
void Mpz::pow_mod (Mpz &r, const Mpz &f, const Mpz &n, const Mpz &m, size_t e,
                   const Mpz &fe, const Mpz &f2e, const Mpz &f3e)
{
  if (n.is_zero ())
  {
    r = 1UL;
  }
  else if (n.nbits() < e)
  {
    pow_mod (r, f, n, m);
  }
  else /* n != 0: exponentiation with abs(n) and handle sign after */
  {
    /* */
    mp_srcptr mp = PTR (m.mpz_);
    mp_size_t nlimbs = ABSIZ(m.mpz_);
    mp_size_t itch = mpn_binvert_itch (nlimbs);
    if (itch < 2*nlimbs)
      itch = 2*nlimbs;

    Mpz mi, T0, T1;
    mp_ptr mip = mpz_limbs_write (mi.mpz_, nlimbs);
    mp_ptr t0 = mpz_limbs_write (T0.mpz_, itch);
    mp_ptr t1 = mpz_limbs_write (T1.mpz_, nlimbs+1);

    mpn_binvert (mip, mp, nlimbs, t0);

#define SQR_REDC(r, a) do {                           \
    mpn_sqr (t0, PTR((a).mpz_), nlimbs);              \
    mpn_redc_n (PTR((r).mpz_), t0, mp, nlimbs, mip);  \
  } while (0)
#define MUL_REDC(r, a, b) do {                            \
    mpn_mul_n (t0, PTR((a).mpz_), PTR((b).mpz_), nlimbs); \
    mpn_redc_n (PTR((r).mpz_), t0, mp, nlimbs, mip);      \
  } while (0)

    /* precomputations */
    /* tab[i] = f^b0 * fe^b1 * f2e^b2 * f3e^b3
     * where [ b0, b1, b2, b3] is (i+1) written in basis 2.
     */
    Mpz tab[15];

    for (size_t i = 0; i < 15; i++)
    {
      _mpz_realloc (tab[i].mpz_, nlimbs);
    }

    redcify (PTR(tab[0].mpz_), PTR(f.mpz_), f.nlimbs(), mp, nlimbs, t0, t1);
    redcify (PTR(tab[1].mpz_), PTR(fe.mpz_), fe.nlimbs(), mp, nlimbs, t0, t1);
    redcify (PTR(tab[3].mpz_), PTR(f2e.mpz_), f2e.nlimbs(), mp, nlimbs, t0, t1);
    redcify (PTR(tab[7].mpz_), PTR(f3e.mpz_), f3e.nlimbs(), mp, nlimbs, t0, t1);

    MUL_REDC (tab[2], tab[1], tab[0]);

    for (size_t i = 0; i < 3; i++)
      MUL_REDC (tab[4+i], tab[3], tab[i]);

    for (size_t i = 0; i < 7; i++)
      MUL_REDC (tab[8+i], tab[7], tab[i]);

    /* */
    mp_ptr rp = mpz_limbs_write (r.mpz_, nlimbs);
    rp[0] = 1UL;
    redcify (rp, rp, 1, mp, nlimbs, t0, t1);

    for (size_t j = n.nbits(); j > 4*e; j--)
    {
      int b = n.tstbit (j-1);
      SQR_REDC (r, r);
      if (b)
        MUL_REDC (r, r, tab[7]);
    }

    for (size_t j = e; j > 0; j--)
    {
      int b0 = n.tstbit (j-1);
      int b1 = n.tstbit (j-1+e);
      int b2 = n.tstbit (j-1+2*e);
      int b3 = n.tstbit (j-1+3*e);
      int idx = (b3 << 3) | (b2 << 2) | (b1 << 1) | b0;

      SQR_REDC (r, r);
      if (idx)
        MUL_REDC (r, r, tab[idx-1]);
    }

#undef SQR_REDC
#undef MUL_REDC

    mpn_copyi (t0, rp, nlimbs);
    mpn_zero (t0 + nlimbs, nlimbs);
    mpn_redc_n (rp, t0, mp, nlimbs, mip);
    SIZ(r.mpz_) = nlimbs;
    MPN_NORMALIZE (rp, SIZ(r.mpz_));
  }
}

/* */
inline
void Mpz::gcd (Mpz &g, const Mpz &a, const Mpz &b)
{
  mpz_gcd (g.mpz_, a.mpz_, b.mpz_);
}

/* */
inline
void Mpz::gcdext (Mpz &g, Mpz &u, Mpz &v, const Mpz &a, const Mpz &b)
{
  mpz_gcdext (g.mpz_, u.mpz_, v.mpz_, a.mpz_, b.mpz_);
}

/* */
inline
void Mpz::lcm (Mpz &g, const Mpz &a, const Mpz &b)
{
  mpz_lcm (g.mpz_, a.mpz_, b.mpz_);
}

/* */
inline
void Mpz::sqrt (Mpz &r, const Mpz &op)
{
  mpz_sqrt (r.mpz_, op.mpz_);
}

/* */
inline
void Mpz::root4th (Mpz &r, const Mpz &op)
{
  mpz_root (r.mpz_, op.mpz_, 4);
}

/*
 * Set r to a square root of s modulo l, i.e., r^2 = s mod l
 *
 * Implementation of Tonelli-Shanks algorithm.
 * Doc: https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
 */
inline
void Mpz::sqrt_mod_prime (Mpz &r, const Mpz &s, const Mpz &l)
{
  Mpz Q, z, c, t, b, tmp;
  size_t M;

  sub (Q, l, 1UL);
  for (M = 0; Q.is_even(); M++, divby2 (Q, Q));
  /* Now Q*2^M = l-1 */

  for (z = 1UL; z.kronecker (l) != -1; add (z, z, 1UL));
  /* z is a non-quadratic residue */

  pow_mod (c, z, Q, l);
  pow_mod (t, s, Q, l);
  add (tmp, Q, 1UL);
  divby2 (tmp, tmp);
  pow_mod (r, s, tmp, l);

  while (!t.is_zero() && t != 1UL)
  {
    size_t i;
    tmp = t;

    for (i = 0; tmp != 1UL; i++, mul (tmp, tmp, tmp), mod (tmp, tmp, l));

    tmp = 1UL;
    mulby2k (tmp, tmp, M-i-1);
    pow_mod (b, c, tmp, l);

    M = i;
    mul (c, b, b); mod (c, c, l); /* c <- b^2 mod l */
    mul (t, t, c); mod (t, t, l); /* t <- t*c mod l = t*b^2 mod l */
    mul (r, r, b); mod (r, r, l); /* r <- r*b mod l */
  }

  if (t.is_zero())
    r = 0UL;
}

/* */
inline
size_t Mpz::remove (Mpz &r, const Mpz &n, const Mpz &f)
{
  return mpz_remove (r.mpz_, n.mpz_, f.mpz_);
}

/* */
inline
void Mpz::CRT (Mpz &r, const Mpz &a1, const Mpz &m1, const Mpz &a2,
                       const Mpz &m2)
{
  Mpz n1, n2, g, tmp;
  Mpz::gcdext (g, n1, n2, m1, m2);
  Mpz::mul (r, m2, n2);
  Mpz::mul (r, r, a1);
  Mpz::mul (tmp, m1, n1);
  Mpz::addmul (r, tmp, a2);
  if (g == 1UL)
  {
    Mpz::mul (tmp, m1, m2);
    Mpz::mod (r, r, tmp);
  }
  else
  {
    // TODO check n1 == n2 modulo g
    Mpz::divexact (r, r, g);
    Mpz::mul (tmp, m1, m2);
    Mpz::divexact (tmp, tmp, g);
    Mpz::mod (r, r, tmp);
  }
}

/*
 * Set r to ceil(log(|n|)^2).
 *
 * Naive implementation: use Taylor expansion of log(1-z)^2
 *      log(1-z)^2 = sum_{i=1}^{infinity}{2*Hi/(i+1) * z^(i+1)}
 *  with Hi = sum_{j=1}^{i}{1/j}
 * Argument reduction from |n| to a value z in ]0, 2[ by taking square root.
 *
 * FIXME: what precision is needed ? Do we only need 2*log2(nbits) of
 * precision (+epsilon) ?
 */
inline
void Mpz::ceil_abslog_square (Mpz &r, const Mpz &n)
{
  const size_t nbits = n.nbits();
  size_t m;

  mpf_t nf, acc, z, pow, H, t, tmp;
  mpf_init2 (nf, nbits);
  mpf_init2 (acc, nbits);
  mpf_init2 (z, nbits);
  mpf_init2 (pow, nbits);
  mpf_init2 (H, nbits);
  mpf_init2 (t, nbits);
  mpf_init2 (tmp, nbits);

  mpf_set_z (nf, n.mpz_);
  mpf_abs (nf, nf);

  /* compute nf and m such that log(|n|)^2 = 4^m*log(nf)^2 with 0<nf<2 */
  m = 0;
  for (m = 0; mpf_cmp_ui (nf, 2) >= 0; m += 1)
  {
    mpf_sqrt (nf, nf);
  }

  mpf_ui_sub (z, 1, nf); /* -1 < z < 1 */
  mpf_mul (pow, z, z);
  mpf_set_ui (acc, 0);
  mpf_set_ui (H, 1);

  for (size_t i = 1; i < 1000000; i++) /* 1e6 harcoded max iter */
  {
    mpf_set_ui (tmp, i+1);
    mpf_ui_div (tmp, 1, tmp);
    mpf_mul (t, pow, H);
    mpf_mul_2exp (t, t, 1);
    mpf_mul (t, t, tmp);
    mpf_add (acc, acc, t);

    mpf_mul (pow, pow, z);
    mpf_add (H, H, tmp);

    mpf_div (tmp, t, acc);
    mpf_abs (tmp, tmp);
    if (mpf_cmp_d (tmp, 1e-9) <= 0) /* TODO: remove hardcoded constant */
      break;
  }

  mpf_mul_2exp (acc, acc, 2*m); /* mul by 2^(2*m) = 4^m */
  mpf_ceil (acc, acc);
  mpz_set_f (r.mpz_, acc);

  mpf_clears (nf, acc, z, pow, H, t, tmp, NULL);
}

/*
 * Input: a, b and target_nlimb
 *
 * output: a and b of length <= target_nlimb and matrix U = (uij)
 *  (0 <= i,j < 2)
 * such that:
 *  - U*(a, b) = (a_input, b_input)
 *  - det U = 1
 *
 * Note: input a and b are overwritten with the output values
 *
 * Assumption: target_nlimb must be >= 1.
 *
 * Scratch variable(s): t0, t1
 */
inline
void Mpz::partial_euclid (Mpz &u00, Mpz &u01, Mpz &u10, Mpz &u11, Mpz &a,
                          Mpz &b, mp_size_t target_nlimb, Mpz &t0, Mpz &t1)
{
  int swapped = 0;

  /* ensure that ABSIZ(b) >= ABSIZ (a) */
  if (ABSIZ (b.mpz_) < ABSIZ (a.mpz_))
  {
    swap (a, b);
    swapped = 1;
  }

  /* signs are handled before swap is undone (if necessary) */
  int a_is_neg = SIZ (a.mpz_) < 0;
  int b_is_neg = SIZ (b.mpz_) < 0;

  if (SIZ (a.mpz_) == 0) /* if a == 0 => identity matrix */
  {
    u00 = 1UL; u01 = 0UL;
    u10 = 0UL; u11 = 1UL;
  }
  else
  {
    mp_size_t n = ABSIZ (b.mpz_); /* Fact: 0 < ABSIZ(a) <= ABSIZ(b) = n */

    /* get the pointer for a and b (and reallocate a if necessary) */
    mp_ptr const ap = mpz_limbs_modify (a.mpz_, n);
    mp_ptr const bp = mpz_limbs_modify (b.mpz_, n);
    mpn_zero (ap+ABSIZ(a.mpz_), n-ABSIZ(a.mpz_));

    /* realloc u10, u11 if necessary, and set u10 to 0 and u11 to 1 */
    mp_ptr const u10p = mpz_limbs_write (u10.mpz_, n+1);
    mp_ptr const u11p = mpz_limbs_write (u11.mpz_, n+1);
    mpn_zero (u10p, n+1);
    mpn_zero (u11p, n+1);
    u11p[0] = 1;
    mp_size_t un = 1;

    /* realloc u00, u01 if necessary, and set u00 to 1 and u01 to 0 */
    mp_ptr const u00p = mpz_limbs_write (u00.mpz_, n+1);
    mp_ptr const u01p = mpz_limbs_write (u01.mpz_, n+1);
    mpn_zero (u00p, n+1);
    mpn_zero (u01p, n+1);
    u00p[0] = 1;
    mp_size_t vn = 1;

    /* get the pointer for t0 and t1 (and reallocate a if necessary) */
    mp_ptr const t0p = mpz_limbs_modify (t0.mpz_, n+1);
    mp_ptr const t1p = mpz_limbs_modify (t1.mpz_, n+1);
    mpn_zero (t0p, n+1);
    mpn_zero (t1p, n+1);

    /* Loop invariant: n == ABSIZ(b) or ABSIZ(a), and both are <= n */
    while (n > target_nlimb)
    {
      hgcd_matrix1 M;
      mp_limb_t ah, al, bh, bl;

      mpn_highest_two_limbs (&ah, &al, &bh, &bl, ap, bp, n);

      /* Try an mpn_hgcd2 step */
      if (mpn_hgcd2 (ah, al, bh, bl, &M))
      {
        /* Compute  M^-1 * (a,b), the result is written in (t1p, bp) then
         * swap to obtain it in (ap ,bp).
         * n is updated to the new max (ABSIZ(a), ABSIZ(b))
         */
        n = mpn_matrix22_mul1_inverse_vector (&M, t1p, ap, bp, n);
        mpn_copyi (ap, t1p, n);
        //MP_PTR_SWAP (ap, t1p);

        /* apply (u10,u11) * M, the result is written in (t0p, u11) then
         * swap to obtain it in (u10 ,u11).
         * un is updated to the new max (ABSIZ(u10), ABSIZ(u11))
         */
        un = mpn_hgcd_mul_matrix1_vector(&M, t0p, u10p, u11p, un);
        mpn_copyi (u10p, t0p, un);
        //MP_PTR_SWAP (u10, t0p);

        /* same for u00 and u01 */
        vn = mpn_hgcd_mul_matrix1_vector(&M, t0p, u00p, u01p, vn);
        mpn_copyi (u00p, t0p, vn);
      }
      else
      {
        /* Looking at the code of mpn_hgcd2, it returns 0 in 3 cases:
         *  1. if ah < 2 or bh < 2: as the most significant bit of ah or bh
         *      must be 1, it means that the ratio a/b is too large
         *  2. (al|ah) >= (bl|bh) and highest limb of (al|ah)-(bl|bh) is < 2
         *  3. (bl|bh) >= (al|ah) and highest limb of (bl|bh)-(al|ah) is < 2
         * For case 1, we perform a div, for 2 and 3 we perform a sub.
         */
        if (bh < 2) /* a is too large compared to b */
        {
          /* do (a, b) <- (a-b*q, b) */
          mp_size_t bn = n - 1;
          MPN_NORMALIZE (bp, bn);
          mpn_tdiv_qr (t1p, ap, 0, ap, n, bp, bn);
          mp_size_t qn = n - bn + 1;
          MPN_NORMALIZE (t1p, qn);
          n = bn;

          vn = mpn_addmul (u01p, t1p, qn, u00p, vn, t0p);
          un = mpn_addmul (u11p, t1p, qn, u10p, un, t0p);
        }
        else if (ah < 2) /* b is too large compared to a */
        {
          /* do (a, b) <- (a, b-a*q) */
          mp_size_t an = n - 1;
          MPN_NORMALIZE (ap, an);
          mpn_tdiv_qr (t1p, bp, 0, bp, n, ap, an);
          mp_size_t qn = n - an + 1;
          MPN_NORMALIZE (t1p, qn);
          n = an;

          vn = mpn_addmul (u00p, t1p, qn, u01p, vn, t0p);
          un = mpn_addmul (u10p, t1p, qn, u11p, un, t0p);
        }
        else /* a and b are too close */
        {
          int c = mpn_cmp (ap, bp, n);
          if (c < 0)
          {
            /* do (a, b) <- (a, b-a) */
            mpn_sub_n (bp, bp, ap, n);

            u00p[vn] = mpn_add_n (u00p, u00p, u01p, vn);
            vn += (u00p[vn] > 0);
            u10p[un] = mpn_add_n (u10p, u10p, u11p, un);
            un += (u10p[un] > 0);
          }
          else
          {
            /* do (a, b) <- (a-b, b) */
            mpn_sub_n (ap, ap, bp, n);

            u01p[vn] = mpn_add_n (u01p, u00p, u01p, vn);
            vn += (u01p[vn] > 0);
            u11p[un] = mpn_add_n (u11p, u10p, u11p, un);
            un += (u11p[un] > 0);
          }
        }
      }
    }
    /* Note that sign is handled before reswapping a and b, if necessary */
    mpz_limbs_finish (a.mpz_, a_is_neg ? -n : n);
    mpz_limbs_finish (b.mpz_, b_is_neg ? -n : n);
    mpz_limbs_finish (u10.mpz_, b_is_neg ^ a_is_neg ? -un : un);
    mpz_limbs_finish (u11.mpz_, un);
    mpz_limbs_finish (u00.mpz_, vn);
    mpz_limbs_finish (u01.mpz_, b_is_neg ^ a_is_neg ? -vn : vn);
  }

  if (swapped)
  {
    /* swap a and b, and the two diagonal of the matrix */
    swap (u10, u01);
    swap (u11, u00);
    swap (a, b);
  }
}

/* */
inline
void Mpz::partial_euclid (Mpz &u00, Mpz &u01, Mpz &u10, Mpz &u11, Mpz &a,
                          Mpz &b, mp_size_t target_nlimb)
{
  Mpz t0, t1;
  partial_euclid (u00, u01, u10, u11, a, b, target_nlimb, t0, t1);
}

/* */
std::ostream & operator<< (std::ostream &o, const Mpz &v)
{
  return o << v.mpz_;
}

/* */
std::istream & operator>> (std::istream &i, Mpz &v)
{
  return i >> v.mpz_;
}

/* */
inline
int Mpz::richcmp (const Mpz &a, const Mpz &b)
{
  return mpz_cmp (a.mpz_, b.mpz_);
}

/* */
inline
int Mpz::richcmp (const Mpz &a, unsigned long v)
{
  return mpz_cmp_ui (a.mpz_, v);
}

/* */
inline
int Mpz::richcmp (const Mpz &a, long v)
{
  return mpz_cmp_si (a.mpz_, v);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/* */
inline
RandGen::RandGen ()
{
  gmp_randinit_default (rand_);
}

/* */
inline
RandGen::RandGen (const RandGen &other)
{
  gmp_randinit_set (rand_, other.rand_);
}

/* */
inline
RandGen::RandGen (const Mpz &seed) : RandGen()
{
  set_seed (seed);
}

/* */
inline
RandGen::~RandGen ()
{
  gmp_randclear (rand_);
}

/* */
inline void
RandGen::set_seed (const Mpz &seed)
{
  gmp_randseed (rand_, seed.mpz_);
}

/* */
inline
Mpz RandGen::random_mpz (const Mpz &m)
{
  Mpz r;
  mpz_urandomm (r.mpz_, rand_, m.mpz_);
  return r;
}

/* */
inline
Mpz RandGen::random_mpz_2exp (mp_bitcnt_t n)
{
  Mpz r;
  mpz_urandomb (r.mpz_, rand_, n);
  return r;
}

/* */
inline
unsigned char RandGen::random_uchar ()
{
  return gmp_urandomb_ui (rand_, CHAR_BIT);
}

/* */
inline
std::vector<unsigned char> RandGen::random_bytes (size_t n)
{
  std::vector<unsigned char> r(n);
  for (size_t i = 0; i < n; i++)
    r[i] = random_uchar();
  return r;
}

/* */
inline
unsigned long RandGen::random_ui (unsigned long m)
{
  return gmp_urandomm_ui (rand_, m);
}

/* */
inline
unsigned long RandGen::random_ui_2exp (mp_bitcnt_t n)
{
  return gmp_urandomb_ui (rand_, n);
}

/* */
inline
Mpz RandGen::random_negative_discriminant (mp_bitcnt_t n)
{
  Mpz D;
  unsigned long Dmod4;
  do
  {
    D = random_mpz_2exp (n);
    D.neg();
    Dmod4 = D.mod4();
  } while (Dmod4 == 2 || Dmod4 == 3);
  return D;
}

/* */
inline
bool RandGen::random_bool ()
{
  return static_cast<bool>(random_ui_2exp (1));
}

/*
 * Return a random prime with exactly nbits bits, i.e., a prime p such that
 *      2^(nbits-1) <= p < 2^nbits
 * If nbits <= 1, it always returns 2.
 */
inline
Mpz RandGen::random_prime (size_t nbits)
{
  Mpz p;

  if (nbits <= 1)
    p = 2UL;
  else if (nbits == 2)
    p = random_bool () ? 3UL : 2UL;
  else
  {
    do {
      mpz_urandomb (p.mpz_, rand_, nbits); /* random 0 <= p < 2^nbits */
      p.setbit (nbits-1);                 /* ensure that p is >= 2^(nbits-1) */
      p.nextprime ();                     /* ensure p is prime */
    } while (p.nbits() != nbits);
  }

  return p;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/* */
inline
JSF::JSF (const Mpz &n0, const Mpz &n1)
  : std::vector<uint8_t> (std::max (n0.nbits(), n1.nbits()) + 1)
{
  int d0 = 0, d1 = 0;
  int n0j = n0.tstbit (0);
  int n0jp1 = n0.tstbit (1);
  int n0jp2 = n0.tstbit (2);
  int n1j = n1.tstbit (0);
  int n1jp1 = n1.tstbit (1);
  int n1jp2 = n1.tstbit (2);

  for (size_t j = 0; j < size(); j++)
  {
    int u0, u1;

    /* bi := (di + 2*nijp1 + nij) % 4 == 2.
     * Computed now, as di may change when we need this test.
     */
    int b0 = d0 == n0j && n0jp1 ^ d0;
    int b1 = d1 == n1j && n1jp1 ^ d1;

    if (d0 == n0j)
      u0 = 0;
    else
    {
      u0 = n0jp1 ^ ((n0jp2 ^ n0jp1) & b1) ? 3 : 1;
      d0 = u0 >> 1;
    }

    if (d1 == n1j)
      u1 = 0;
    else
    {
      u1 = n1jp1 ^ ((n1jp2 ^ n1jp1) & b0) ? 3 : 1;
      d1 = u1 >> 1;
    }

    set (j, u0, u1);

    n0j = n0jp1; n0jp1 = n0jp2; n0jp2 = n0.tstbit (j+3);
    n1j = n1jp1; n1jp1 = n1jp2; n1jp2 = n1.tstbit (j+3);
  }

  if (operator[] (size()-1) == 0)
    resize (size()-1);
}

/* */
inline
uint8_t JSF::operator[] (size_t i) const
{
  if (i < size())
    return std::vector<uint8_t>::operator[] (i);
  else
    return 0U;
}

/* */
inline
void JSF::set (size_t i, int d0, int d1)
{
  std::vector<uint8_t>::operator[] (i) = (d1 << 4) | d0;
}

#undef ALLOC
#undef PTR
#undef SIZ
#undef ABS
#undef ABSIZ
#undef MPN_NORMALIZE
#undef GMP_NUMB_HIGHBIT
#undef MPN_EXTRACT_NUMB
#undef count_leading_zeros

#undef mpn_hgcd2
#undef mpn_matrix22_mul1_inverse_vector
#undef mpn_hgcd_mul_matrix1_vector
#undef mpn_binvert
#undef mpn_binvert_itch
#undef mpn_redc_n

#endif /* BICYCL_GMP_EXTRAS_INL */
