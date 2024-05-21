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
#ifndef BICYCL_QFI_INL
#define BICYCL_QFI_INL

/* */
inline
QFI::QFI () : a_(1UL), b_(1UL), c_(1UL)
{
}

/* */
inline
QFI::QFI (const Mpz &a, const Mpz &b, const Mpz &c, bool bypass_check)
    : a_(a), b_(b), c_(c)
{
  if (!bypass_check)
  {
    if (discriminant().sgn() >= 0)
      throw std::runtime_error ("form must have negative discriminant");

    Mpz g;
    Mpz::gcd (g, a_, b_);
    Mpz::gcd (g, g, c_);
    if (g != 1UL)
      throw std::runtime_error ("form must be primitive");

    if (a_.sgn() <= 0)
      throw std::runtime_error ("the coeff a of the form must be positive");

    if (Mpz::cmpabs (a_, b_) <= 0 && a_ != b_)
      throw std::runtime_error ("form must be normal");

    if (a_ > c_ || (a_ == c_ && b_.sgn() < 0))
      throw std::runtime_error ("form must be reduced");
  }
}

/* */
inline
bool QFI::operator== (const QFI &other) const
{
  return a_ == other.a_ && b_ == other.b_ && c_ == other.c_;
}

/* */
inline
QFI::QFI (const QFICompressedRepresentation &c, const Mpz &disc)
    : QFI()
{
  if (c.g.is_zero() && c.tp.is_zero() && c.b0.is_zero() && not c.is_neg)
  {
    a_ = c.ap;
    b_ = 0UL;
    set_c_from_disc (disc);
  }
  else if (c.tp.is_zero())
  {
    a_ = c.g;
    b_ = c.g;
    set_c_from_disc (disc);
  }
  else
  {
    Mpz t, x, s, sp, bp, f, l;
    Mpz::mul (a_, c.ap, c.g);
    Mpz::mul (t, c.tp, c.g);

    Mpz::mul (x, t, t);
    Mpz::mod (x, x, a_);
    Mpz::mul (x, x, disc);
    Mpz::mod (x, x, a_);  /* x <- t^2*disc mod a */

    Mpz::sqrt (s, x);
    Mpz::divexact (sp, s, c.g);
    Mpz::mod_inverse (bp, c.tp, c.ap); /* typo in the article: t' not t */
    Mpz::mul (bp, bp, sp);
    Mpz::mod (bp, bp, c.ap);

    Mpz::abs (f, c.g);
    Mpz::lcm (l, f, c.ap);
    while (l < a_)
    {
      Mpz::add (f, f, 1);
      Mpz::lcm (l, f, c.ap);
    }

    Mpz::CRT (b_, bp, c.ap, c.b0, f);
    if (c.is_neg)
      b_.neg();
    set_c_from_disc (disc);
  }
}

/* */
inline
const Mpz & QFI::a () const
{
  return a_;
}

/* */
inline
const Mpz & QFI::b () const
{
  return b_;
}

/* */
inline
const Mpz & QFI::c () const
{
  return c_;
}

/*
 * Return the discriminant of the qfi.
 */
inline
Mpz QFI::discriminant () const
{
  Mpz d;
  Mpz::mul (d, a_, c_);     /* a*c */
  Mpz::mulby4 (d, d);       /* 4*a*c */
  Mpz::submul (d, b_, b_);  /* 4*a*c - b^2 */
  d.neg();                  /* d = b^2 - 4*a*c */
  return d;
}

/*
 * return true if the form is the neutral element of the class group
 */
inline
bool QFI::is_one () const
{
  return a_.is_one(); /* assume the form is reduced */
}

/*
 * Inverse the qfi
 * The inverse of the form (a, b, c) is (a, -b, c).
 *
 * As we assumed the input is in reduced form and we want to output a reduced
 * form, there are some corner cases:
 *  - if a == c: (a, -b, c) is not reduced (we need the second coeff to be
 *    nonnegative in this case), but (c, b, a) is an equivalent reduced form. As
 *    a == c, there is nothing to do in this case.
 *  - if a == b: (a, -b, c) is not reduced (we need the second coeff to be in
 *    ]-a, a]). But (a, -b+2*a, a-b+c) = (a, b, c) is an equivalent reduced
 *    form. Again there is nothing to do in this case.
 *  - otherwise a < c and |b| < a and (a, -b, c) is in reduced form.
 */
inline
void QFI::neg ()
{
  if (a_ != c_ && a_ != b_)
    b_.neg();
}

/* */
inline
Mpz QFI::eval (const Mpz &x, const Mpz &y) const
{
  Mpz v(0UL), t;
  Mpz::mul (t, x, x);
  Mpz::mul (v, a_, t);
  Mpz::mul (t, x, y);
  Mpz::addmul (v, b_, t);
  Mpz::mul (t, y, y);
  Mpz::addmul (v, c_, t);
  return v;
}

/* */
inline
QFICompressedRepresentation QFI::compressed_repr () const
{
  Mpz zero (0UL), one(1UL);

  if (b_.is_zero())
    return QFICompressedRepresentation (a_, zero, zero, zero, false);
  else if (a_ == b_)
    return QFICompressedRepresentation (one, a_, zero, zero, false);
  else
  {
    Mpz b, g, ap, f, l, b0;
    bool is_neg = (b_.sgn() < 0);
    Mpz::abs (b, b_);

    Mpz s(b), sp (a_), t(1UL), tp(0UL), sq, q;
    Mpz::sqrt (sq, a_);
    while (s >= sq)
    {
      Mpz::fdiv_qr (q, sp, sp, s);
      Mpz::swap (sp, s);
      Mpz::submul (tp, q, t);
      Mpz::swap (tp, t);
    }

    Mpz::gcd (g, a_, t);

    Mpz::divexact (ap, a_, g);
    Mpz::divexact (tp, t, g);

    Mpz::abs (f, g);
    Mpz::lcm (l, f, ap);
    while (l < a_)
    {
      Mpz::add (f, f, 1);
      Mpz::lcm (l, f, ap);
    }

    Mpz::mod (b0, b, f);
    return QFICompressedRepresentation (ap, g, tp, b0, is_neg);
  }
}

/*
 * Lift the qfi into an qfi of discriminant l^2 times the discriminant of the
 * form.
 *
 * Assumes l is a prime power.
 *
 * Ref: Algorithm 2 (GoToNonMaxOrder) of [HJPT1998]_.
 *
 * The l-lift of f=(a,b,c) is (a, l*b, l^2*c) if a is prime to l.
 * So, first we need to move f to a an equivalent representation with the first
 * coeff prime to l. Then we apply the above formula. And finally we reduce the
 * new form.
 */
inline
void QFI::lift (const Mpz &l)
{
  prime_to (l);
  /* a stays the same, b is set to b*l and c is set to c*l^2 */
  Mpz::mul (b_, b_, l);
  Mpz::mul (c_, c_, l);
  Mpz::mul (c_, c_, l);
  reduction();
}

/*
 * Same as QFI::lift for the case of power of l=2^k.
 */
inline
void QFI::lift_2exp (unsigned int k)
{
  prime_to_2exp ();
  /* a stays the same, b is set to b*2^k and c is set to c*2^(2*k) */
  Mpz::mulby2k (b_, b_, k);
  Mpz::mulby2k (c_, c_, 2*k);
  reduction();
}

/*
 * Move the qfi into the maximal order.
 *
 * Assumptions:
 *  - l is a odd prime power
 *  - DeltaK is a fundamental discriminant and odd
 *  - the discriminant is l^2 times DeltaK
 *
 * Ref: Algorithm 3 (GoToMaxOrder) of [HJPT1998]_.
 *
 * The image of (a,b,c) is, if a is prime to l, (a, bu + adv, ...) where
 *    - 1 = u*l + v*a (so a must be prime to l).
 *    - d = discriminant % 2
 * So, first we need to move f to a an equivalent representation with the
 * first coeff prime to l. Then we apply the above formula.
 * Then, if with_reduction is true (the default), the form is finally reduced.
 */
inline
void QFI::to_maximal_order (const Mpz &l, const Mpz &DeltaK,
                            bool with_reduction=true)
{
  Mpz tmp0, g0, g1;
  prime_to (l);
  Mpz::gcdext (tmp0, g0, g1, l, a_);
  /* As we assume that DeltaK and l are odd, d equals 1 */
  Mpz::mul (b_, b_, g0);
  Mpz::addmul (b_, a_, g1);
  set_c_from_disc (DeltaK);
  if (with_reduction)
    reduction();
}

/*
 * Same as QFI::to_maximal_order for the case of power of l=2^k.
 * Assumptions:
 *  - DeltaK is a fundamental discriminant and even
 *  - the discriminant is l^2 times DeltaK
 */
inline
void QFI::to_maximal_order_2exp (unsigned int k, const Mpz &DeltaK,
                                 bool with_reduction=true)
{
  prime_to_2exp ();
  Mpz u, v;
  Mpz::mod_inverse_2k (v, a_, k, u); /* u is used as temp variable here */
  Mpz::mul (u, v, a_);
  Mpz::sub (u, u, 1);
  u.neg();
  Mpz::divby2k (u, u, k); /* u = (1-v*a)/2^k */
  /* As we assume that DeltaK and l are even, d equals 0 */
  Mpz::mul (b_, b_, u);
  set_c_from_disc (DeltaK);

  if (with_reduction)
    reduction();
}

/*
 * Compute a representative of the form in (OK / lOK)* / (Z/lZ)*.
 *
 * Assumption:
 *  - the discriminant of the form is l^2 * DeltaK
 *  - l is a prime power
 *  - l and DeltaK are odd
 *  - the form belongs to the kernel of the projection from Cl(O) to
 *    Cl(O_DeltaK) (std::invalid_argument is thrown otherwhise)
 */
inline
Mpz QFI::kernel_representative (const Mpz &l, const Mpz &DeltaK) const
{
  Mpz tmp0, tmp1, g0, g1;
  /* Move f to the O_DeltaK
   */
  QFI ft(*this);
  ft.to_maximal_order (l, DeltaK, false); /* no reduction */

  /* Reduce ft and build gamma while doing it. Each time we apply rho to the
   * form (a,b,c) gamma is multiplied by (b+sqrt(DeltaK))/(2*a).
   * We do not care about the 2*a in the denominator as at the end we are going
   * to compute a representative of gamma^(-1) in the kernel as -g1/g0 mod l. We
   * just have to remove common factor in g0 and g1 before taking the invert.
   */
  int cmp;
  g0 = 1UL; g1 = 0UL; /* start with g = 1 + 0 sqrt(DeltaK) */
  ft.normalize (tmp0, tmp1);
  while ((cmp = Mpz::cmpabs (ft.a_, ft.c_)) > 0) /* while a larger than c */
  {
    Mpz::mul (tmp0, g1, DeltaK);    /* tmp0 <- g1*DeltaK */
    Mpz::mul (g1, g1, ft.b_);       /* g1 <- g1*b */
    Mpz::add (g1, g1, g0);          /* g1 <- g1*b + g0 */
    Mpz::mul (g0, g0, ft.b_);       /* g0 <- g0*b */
    Mpz::add (g0, g0, tmp0);        /* g0 <- g0*b + g1*DeltaK */

    ft.rho (tmp0, tmp1);
  }

  /* Throw invalid_argument if f is not (1, 1, ...), the reduced neutral form
   * for odd values of DeltaK.
   */
  if (ft.a_ != 1UL or ft.b_ != 1UL)
    throw std::invalid_argument ("the form is not in the kernel");

  Mpz::gcd (tmp1, g0, g1);
  Mpz::divexact (g0, g0, tmp1);
  Mpz::divexact (g1, g1, tmp1);

  Mpz::mod_inverse (tmp0, g0, l);
  g1.neg();
  Mpz::mul (tmp0, tmp0, g1);
  Mpz::mod (tmp0, tmp0, l);

  return tmp0;
}

/*
 * Compute a representative of the form in (OK / 2^kOK)* / (Z/2^kZ)*.
 *
 * Assumption:
 *  - the discriminant of the form is 2^(2*k) * DeltaK
 *  - DeltaK is even
 *  - the form belongs to the kernel of the projection from Cl(O) to
 *    Cl(O_DeltaK) (std::invalid_argument is thrown otherwhise)
 */
inline
Mpz QFI::kernel_representative_2exp (size_t k, const Mpz &DeltaK) const
{
  Mpz tmp0, tmp1, g0, g1;
  QFI ft(*this);
  ft.to_maximal_order_2exp (k, DeltaK, false); /* no reduction */

  /* Reduce ft and build gamma while doing it. */
  int cmp;
  g0 = 1UL; g1 = 0UL; /* start with g = 1 + 0 sqrt(DeltaK) */
  ft.normalize (tmp0, tmp1);
  while ((cmp = Mpz::cmpabs (ft.a_, ft.c_)) > 0) /* while a larger than c */
  {
    Mpz::mul (tmp0, g1, DeltaK);    /* tmp0 <- g1*DeltaK */
    Mpz::mul (g1, g1, ft.b_);       /* g1 <- g1*b */
    Mpz::add (g1, g1, g0);          /* g1 <- g1*b + g0 */
    Mpz::mul (g0, g0, ft.b_);       /* g0 <- g0*b */
    Mpz::add (g0, g0, tmp0);        /* g0 <- g0*b + g1*DeltaK */

    ft.rho (tmp0, tmp1);
  }

  /* Throw invalid_argument if f is not (1, 0, ...), the reduced neutral form
   * for even values of DeltaK.
   */
  if (ft.a_ != 1UL or ft.b_.sgn() != 0)
    throw std::invalid_argument ("the form is not in the kernel");

  size_t v = g0.val2();
  Mpz::divby2k (g0, g0, v);
  Mpz::divby2k (g1, g1, v);

  Mpz::mod_inverse_2k (tmp0, g0, k, tmp1);
  Mpz::mod2k (tmp0, tmp0, k);
  g1.neg();
  Mpz::mul (tmp0, tmp0, g1);
  Mpz::mod2k (tmp0, tmp0, k);

  return tmp0;
}


/* */
inline
std::ostream & operator<< (std::ostream &o, const QFI &f)
{
  return o << "(" << f.a_ << ", " << f.b_ << ", " << f.c_ << ")";
}

/*
 * Set the coefficient c of the qfi given its discriminant.
 * Assumes that the coefficients a and b are already set and that 4*a
 * divides b^2-disc.
 */
inline
void QFI::set_c_from_disc (const Mpz &disc)
{
  Mpz::mul (c_, b_, b_);      /* b^2 */
  Mpz::sub (c_, c_, disc);    /* b^2-disc */
  Mpz::divexact (c_, c_, a_); /* (b^2-disc)/a */
  Mpz::divby4 (c_, c_);       /* (b^2-disc)/(4*a) */
}

/*
 * Normalize the qfi.
 * The form f is called normal if -a < b <= a.
 *
 * Note: for normalization (when not called via rho), 99% of the time q is
 * -1 (and 99.9% we have |q| <= 2)
 */
inline
void QFI::normalize ()
{
  Mpz q, r;
  normalize (q, r);
}

inline
void QFI::normalize (Mpz &q, Mpz &r)
{
  Mpz::cdiv_qr (q, r, b_, a_); /* b = q*a + r    and    -a < r <= 0 */
  if (q.is_odd())
    Mpz::add (r, r, a_);
  Mpz::divby2 (q, q); /* divide q by 2 */
  /* Now we have b = (2*a)*q + r    and   -a < r <= a */
  Mpz::swap (b_, r);
  Mpz::add (r, b_, r);
  Mpz::divby2 (r, r);
  Mpz::submul (c_, q, r);
}

/*
 * Apply the rho transformation to the qfi
 * rho(f=(a,b,c)) is defined as the normalization of the form (c,-b,a)
 */
inline
void QFI::rho ()
{
  Mpz q, r;
  rho (q, r);
}

inline
void QFI::rho (Mpz &q, Mpz &r)
{
  /* The 2 versions seems to be equivalent (first even seems slightly
   * faster) */
#if 1
    Mpz::swap (a_, c_);
    b_.neg();
    normalize (q, r);
#else
  Mpz::fdiv_qr (q, r, b_, c_); /* b = q*c + r    and   0 <= r < c */
  if (q.is_odd())
    Mpz::sub (r, r, c_);
  mpzc_div_q_2exp (q, q, 1); /* divide q by 2 */
  /* Now we have b = (2*c)*q + r    and   -c <= r < c,
   * It implies     -c < -b + (2*c)*q <= c
   */
  Mpz::add (b_, b_, r);
  Mpz::divby2 (b_, b_);
  Mpz::submul (a_, q, b_);
  mpz_neg (b_, r);
  Mpz::swap (a_, c_);
#endif
}

/*
 * Reduce the qfi f.
 * The form f is reduced if the form is
 *      normal (-a < b <= a)
 *    and
 *      a < c   or    a = c and b >= 0.
 */
inline
void QFI::reduction ()
{
  Mpz t0, t1;
  reduction (t0, t1);
}

inline
void QFI::reduction (Mpz &t0, Mpz &t1)
{
  int cmp;
  normalize (t0, t1);
  /* We know a and c > 0, do not consider signs for comparisons */
  while ((cmp = Mpz::cmpabs (a_, c_)) > 0) /* while a larger than c */
    rho (t0, t1);

  if (cmp == 0 && b_.sgn() < 0) /* if a == c, we need b positive */
    b_.neg();
}

/*
 * Modify the form into an equivalent qfi such that l is coprime with its first
 * coeff.
 *
 * Assumes l is a prime power.
 *
 * Ref: Algorithm 1 (FindIdealPrimeTo) of [HJPT1998]_.
 *
 * Remark: the output form is not necessarily reduced.
 */
inline
void QFI::prime_to (const Mpz &l)
{
  Mpz g;
  Mpz::gcd (g, a_, l);
  if (g > 1UL)
  {
    Mpz::gcd (g, c_, l);
    if (g > 1UL) /* transform f into (a+b+c, -b-2a, a) */
    {
      Mpz::add (c_, c_, a_);
      Mpz::add (c_, c_, b_);
      Mpz::add (b_, b_, a_);
      Mpz::add (b_, b_, a_);
      b_.neg();
      Mpz::swap (a_, c_);
    }
    else /* c_ is coprime to l: transform f into (c, -b, a) */
    {
      Mpz::swap (a_, c_);
      b_.neg();
    }
  }
  /* else do nothing if a_ is already coprime to l */
}

/*
 * Same as QFI::prime_to for the case of power of 2.
 */
inline
void QFI::prime_to_2exp ()
{
  if (a_.is_even())
  {
    if (c_.is_even()) /* transform f into (a+b+c, -b-2a, a) */
    {
      Mpz::add (c_, c_, a_);
      Mpz::add (c_, c_, b_);
      Mpz::add (b_, b_, a_);
      Mpz::add (b_, b_, a_);
      b_.neg();
      Mpz::swap (a_, c_);
    }
    else /* c_ is odd: transform f into (c, -b, a) */
    {
      Mpz::swap (a_, c_);
      b_.neg();
    }
  }
  /* else do nothing if a_ is already odd */
}

#ifdef BICYCL_WITH_TIMINGS
static uint64_t _nucomp_ttot, _nucomp_tgcdext, _nucomp_tpartial, _nucomp_treduc,
                _nudupl_ttot, _nudupl_tgcdext, _nudupl_tpartial, _nudupl_treduc;
#endif
#ifdef BICYCL_WITH_COUNTS
static uint64_t _nucomp_ncalls, _nudupl_ncalls;
#endif

/*
 * Compute the composition of two qfi.
 *
 * Input:
 *  - f1 and f2: the two input qfi
 *  - L: bound for partial reduction
 *  - negf2: negate f2 before doing the composition
 * Ouput:
 *  - r: output qfi corresponding to the composition of f1 with f2 (or
 *        f2^(-1) if negf2 is nonzero).
 *
 * Assumption: f1 and f2 must have the same discriminant.
 *
 * Note on the algorithm:
 *  Ref: "A Note on NUCOMP", Alfred J. van der Poorten in Mathematics of
 *        computation, vol 72, number 244. 2003.
 *
 *  Let (ai, bi, ci) be the coefficients of the input form fi.
 *  m and s are defined as:
 *    s = (b1+b2)/2
 *    m = (b2-b1)/2
 *
 *  The goal of the algorithm is to find a 2x4 matrix
 *
 *      Ax Bx Cx Dx
 *      Ay By Cy Dy
 *
 *  such that the minors are:
 *          minor_AB = a1, minor_AC = a2,
 *          minor_AD = s, minor_BC = m,
 *          minor_BD = c2, minor_CD = c1
 *
 *  The returned qfi is then computed as
 *    r = (By*Cy-Ay*Dy, (Ax*Dy+Ay*Dx)-(Bx*Cy+By*Cx), Bx*Cx-Ax*Dx)
 *
 *  The 2x4 matrix can be multiplied on the left by any matrix of SL2(Z).
 *
 *  The two following vectors are orthogonal to the two rows:
 *    v1 = (m, -a2, a1, 0)
 *    v2 = (c2, -s, 0, a1)
 *  This fact allows to retrieve Cx, Cy, Dx, Dy from Ax, Ay, Bx, By.
 *    v1 orthogonal to row1 <=> m*Ax - a2*Bx + a1*Cx == 0
 *    v1 orthogonal to row2 <=> m*Ay - a2*By + a1*Cy == 0
 *    v2 orthogonal to row1 <=> c2*Ax - s*Bx + a1*Dx == 0
 *    v2 orthogonal to row2 <=> c2*Ay - s*By + a1*Dy == 0
 *
 *  A possible solution matrix (used as the initial matrix) is given by
 *
 *    G  (m*v)/H+(l*By)/H   -m*x*u+y*c1+k*a2/G  -x*(v*c1+u*c2) + k*s/G
 *    0  a1/G               a2/G                s/G
 *
 *  where F = gcd (a1, a2) = u*a1 + v*a2
 *        G = gcd (F, s) = x*F + y*s
 *        H = F/G
 *        l = y*(v*c1+u*c2) + k*H
 *        k can be any integer
 *
 *  Then a partial euclidean algorithm is computed with Bx and By to reduce
 *  the size of the coefficients of the matrix.
 *  The returned form is computed with this new matrix.
 *
 * TODO: reduce the number of scratch variables
 */
inline
void QFI::nucomp (QFI &r, const QFI &f1, const QFI &f2, const Mpz &L,
                  bool negf2, OpsAuxVars &tmp)
{
#ifdef BICYCL_WITH_COUNTS
  _nucomp_ncalls++;
#endif
#ifdef BICYCL_WITH_TIMINGS
  uint64_t _time, _time0 = get_realtime_ns();
#endif
  mp_size_t Lsize = L.nlimbs(); /* depending on log2(L) may need to do -1 */

  Mpz::add (tmp.s, f1.b_, f2.b_);
  Mpz::divby2 (tmp.s, tmp.s);
  Mpz::sub (tmp.m, f2.b_, tmp.s);

  if (negf2)
  {
    Mpz::swap (tmp.s, tmp.m);
    tmp.s.neg();
    tmp.m.neg();
  }

  /* F = gcd (a1, a2) = u*a1 + v*a2 */
#ifdef BICYCL_WITH_TIMINGS
  _time = get_realtime_ns();
#endif
  Mpz::gcdext (tmp.F, tmp.u, tmp.v, f1.a_, f2.a_);
#ifdef BICYCL_WITH_TIMINGS
  _nucomp_tgcdext += get_realtime_ns() - _time;
#endif
  if (tmp.F.is_one())
  {
    tmp.Ax = 1UL;
    Mpz::mul (tmp.Bx, tmp.m, tmp.v); /* Bx = m*v */
    tmp.By = f1.a_;
  }
  else if (tmp.s.is_divisible_by (tmp.F))
  {
    tmp.Ax = tmp.F;
    Mpz::mul (tmp.Bx, tmp.m, tmp.v); /* Bx = m*v */
    Mpz::divexact (tmp.By, f1.a_, tmp.Ax);
  }
  else
  {
    Mpz::gcdext (tmp.Ax, tmp.x, tmp.y, tmp.F, tmp.s);
    Mpz::divexact (tmp.H, tmp.F, tmp.Ax);
    Mpz::mod (tmp.t0, f1.c_, tmp.H);
    Mpz::mod (tmp.t1, f2.c_, tmp.H);
    Mpz::mul (tmp.t0, tmp.t0, tmp.v);
    Mpz::addmul (tmp.t0, tmp.t1, tmp.u);
    Mpz::mod (tmp.t0, tmp.t0, tmp.H);
    Mpz::mul (tmp.t0, tmp.t0, tmp.y);
    Mpz::mod (tmp.l, tmp.t0, tmp.H);

    Mpz::divexact (tmp.By, f1.a_, tmp.Ax);

    Mpz::mul (tmp.t0, tmp.v, tmp.m);
    Mpz::addmul (tmp.t0, tmp.l, tmp.By);
    Mpz::divexact (tmp.Bx, tmp.t0, tmp.H);
  }
  Mpz::divexact (tmp.Cy, f2.a_, tmp.Ax);
  Mpz::divexact (tmp.Dy, tmp.s, tmp.Ax);

  /* Set Bx to Bx mod By */
  Mpz::fdiv_qr (tmp.q, tmp.Bx, tmp.Bx, tmp.By);

  /* Partially reduce Bx and By with Euclidean algorithm and compute
   * transformation matrix M.
   */
  tmp.by = tmp.By;
#ifdef BICYCL_WITH_TIMINGS
  _time = get_realtime_ns();
#endif
  Mpz::partial_euclid (tmp.m00, tmp.m01, tmp.m10, tmp.m11, tmp.Bx, tmp.by,
                       Lsize, tmp.t0, tmp.t1);
#ifdef BICYCL_WITH_TIMINGS
  _nucomp_tpartial += get_realtime_ns() - _time;
#endif

  /* new Ay */
  Mpz::mul (tmp.Ay, tmp.m10, tmp.Ax);
  tmp.Ay.neg ();

  /* new Cx and Cy */
  Mpz::mul (tmp.Cx, tmp.Bx, tmp.Cy);
  Mpz::submul (tmp.Cx, tmp.m, tmp.m11);
  Mpz::divexact (tmp.Cx, tmp.Cx, tmp.By);

  if (tmp.Bx.is_zero ())
  {
    Mpz::mul (tmp.Cy, f2.a_, tmp.by);
    Mpz::submul (tmp.Cy, tmp.Ay, tmp.m);
    Mpz::divexact (tmp.Cy, tmp.Cy, f1.a_);
  }
  else
  {
    Mpz::mul (tmp.Cy, tmp.Cx, tmp.by);
    Mpz::add (tmp.Cy, tmp.Cy, tmp.m);
    Mpz::divexact (tmp.Cy, tmp.Cy, tmp.Bx);
  }

  /* new Dx and Dy */
  Mpz::mul (tmp.Dx, tmp.Bx, tmp.Dy);
  Mpz::submul (tmp.Dx, f2.c_, tmp.m11);
  Mpz::divexact (tmp.Dx, tmp.Dx, tmp.By);

  Mpz::submul (tmp.Dy, tmp.Dx, tmp.m10);
  Mpz::divexact (tmp.Dy, tmp.Dy, tmp.m11);

  /* new Ax */
  Mpz::mul (tmp.Ax, tmp.m11, tmp.Ax);

  /* r.a_ = by*Cy - Ay*Dy */
  Mpz::mul (r.a_, tmp.by, tmp.Cy);
  Mpz::submul (r.a_, tmp.Ay, tmp.Dy);

  /* r.c_ = Bx*Cx - Ax*Dx */
  Mpz::mul (r.c_, tmp.Bx, tmp.Cx);
  Mpz::submul (r.c_, tmp.Ax, tmp.Dx);

  /* r.b_ = (Ax*Dy+Ay*Dx) - (Bx*Cy + by*Cx) */
  Mpz::mul (r.b_, tmp.Ax, tmp.Dy);
  Mpz::addmul (r.b_, tmp.Ay, tmp.Dx);
  Mpz::submul (r.b_, tmp.Bx, tmp.Cy);
  Mpz::submul (r.b_, tmp.by, tmp.Cx);

#ifdef BICYCL_WITH_TIMINGS
  _time = get_realtime_ns();
#endif
  r.reduction (tmp.t0, tmp.t1);
#ifdef BICYCL_WITH_TIMINGS
  _nucomp_treduc += get_realtime_ns() - _time;
  _nucomp_ttot += get_realtime_ns() - _time0;
#endif
}

/* */
inline
void QFI::nucomp (QFI &r, const QFI &f1, const QFI &f2, const Mpz &L,
                                                        bool negf2)
{
  QFI::OpsAuxVars tmp;
  nucomp (r, f1, f2, L, negf2, tmp);
}

/*
 * Special case of nucomp where f1 = f2
 * Input:
 *  - f: input qfi
 *  - L: bound for partial reduction
 * Ouput:
 *  - r: output qfi corresponding to the composition of f with itself.
 *
 * Remarks (see QFI_nucomp_scratch description for more details):
 *  Let (a, b, c) be the coefficients of the input form f.
 *
 *  The conditions on the minors become:
 *         minor_AB = minor_AC = a,
 *         minor_AD = b, minor_BC = 0,
 *         minor_BD = minor_CD = c
 *
 *  A possible solution matrix (used as the initial matrix) is given by
 *
 *    d  v*c  v*c  -u*c
 *    0  a/d  a/d  b/d
 *
 *  where d = gcd (a, b) = u*a + v*b
 *
 *  The two following vectors are orthogonal to the two rows:
 *    v1 = (0, -1, 1, 0)
 *    v2 = (c, -b, 0, a)
 *  This fact allows to retrieve Cx, Cy, Dx, Dy from Ax, Ay, Bx, By.
 *    v1 orthogonal to row1 <=> Bx == Cx
 *    v1 orthogonal to row2 <=> By == Cy
 *    v2 orthogonal to row1 <=> c*Ax - b*Bx + a*Dx == 0
 *    v2 orthogonal to row2 <=> c*Ay - b*By + a*Dy == 0
 */
inline
void QFI::nudupl (QFI &r, const QFI &f, const Mpz &L, OpsAuxVars &tmp)
{
#ifdef BICYCL_WITH_COUNTS
  _nudupl_ncalls++;
#endif
#ifdef BICYCL_WITH_TIMINGS
  uint64_t _time, _time0 = get_realtime_ns();
#endif

  mp_size_t Lsize = L.nlimbs();

  /* Ax = d = gcd(a,b) = u*a + v*b (u and v are stored in m11 and m01) */
#ifdef BICYCL_WITH_TIMINGS
  _time = get_realtime_ns();
#endif
  Mpz::gcdext (tmp.Ax, tmp.m11, tmp.m01, f.a_, f.b_);
#ifdef BICYCL_WITH_TIMINGS
  _nudupl_tgcdext += get_realtime_ns() - _time;
#endif

  /* Set By to a/d and Dy to b/d (resp. stored in a_ and b_) */
  if (!tmp.Ax.is_one ())
  {
    Mpz::divexact (r.a_, f.a_, tmp.Ax);
    Mpz::divexact (r.b_, f.b_, tmp.Ax);
  }
  else
  {
    r.a_ = f.a_;
    r.b_ = f.b_;
  }

  /* Set Dx to -u*c */
  Mpz::mul (tmp.Dx, f.c_, tmp.m11);
  tmp.Dx.neg ();

  /* Set Bx to v*c (stored in r.c_) */
  Mpz::mul (r.c_, f.c_, tmp.m01);

  /* Compute q = floor(Bx/By) and set Bx to Bx - q*By (= Bx mod By) and Dx
   * to Dx - q*Dy (i.e., apply matrix [[1, -q], [0, 1]])
   */
  Mpz::fdiv_qr (tmp.t0, r.c_, r.c_, r.a_);
  Mpz::submul (tmp.Dx, tmp.t0, r.b_);

  /* Partially reduce Bx and By with Euclidean algorithm and compute
   * transformation matrix M.
   */
#ifdef BICYCL_WITH_TIMINGS
  _time = get_realtime_ns();
#endif
  Mpz::partial_euclid (tmp.m00, tmp.m01, tmp.m10, tmp.m11, r.c_, r.a_, Lsize,
                       tmp.t0, tmp.t1);
#ifdef BICYCL_WITH_TIMINGS
  _nudupl_tpartial += get_realtime_ns() - _time;
#endif

  /* apply M^-1 to (Ax, Ay) and (Dx, Dy) (Ay is stored in t1) */
  Mpz::mul (tmp.t1, tmp.Ax, tmp.m10);
  tmp.t1.neg (); /* Ax*(-m10) + Ay*m00 (with Ay=0) */
  Mpz::mul (tmp.Ax, tmp.Ax, tmp.m11); /* Ax*m11 - Ay*m10 (with Ay=0) */

  Mpz::mul (tmp.t0, tmp.Dx, tmp.m11);
  Mpz::submul (tmp.t0, r.b_, tmp.m01); /* Dx*m11 - Dy*m01 */
  Mpz::mul (r.b_, r.b_, tmp.m00);
  Mpz::submul (r.b_, tmp.Dx, tmp.m10); /* Dx*(-m10) + Dy*m00 */
  Mpz::swap (tmp.t0, tmp.Dx);

  /* temporarily store Cy*Bx (=By*Bx) in t0 */
  Mpz::mul (tmp.t0, r.a_, r.c_);

  /* a_ = By*Cy - Ay*Dy (with By == Cy) */
  Mpz::mul (r.a_, r.a_, r.a_);
  Mpz::submul (r.a_, tmp.t1, r.b_);

  /* c_ = Bx*Cx - Ax*Dx (with Bx == Cx) */
  Mpz::mul (r.c_, r.c_, r.c_);
  Mpz::submul (r.c_, tmp.Ax, tmp.Dx);

  /* b_ = (Ax*Dy+Ay*Dx) - (Bx*Cy + By*Cx) (with Bx == Cx and By == Cy) */
  Mpz::mul (r.b_, tmp.Ax, r.b_);
  Mpz::addmul (r.b_, tmp.t1, tmp.Dx);
  Mpz::mulby2 (tmp.t0, tmp.t0); /* mul by 2 */
  Mpz::sub (r.b_, r.b_, tmp.t0);

#ifdef BICYCL_WITH_TIMINGS
  _time = get_realtime_ns();
#endif
  r.reduction (tmp.t0, tmp.t1);
#ifdef BICYCL_WITH_TIMINGS
  _nudupl_treduc += get_realtime_ns() - _time;
  _nudupl_ttot += get_realtime_ns() - _time0;
#endif
}

/* */
inline
void QFI::nudupl (QFI &r, const QFI &f, const Mpz &L)
{
  QFI::OpsAuxVars tmp;
  nudupl (r, f, L, tmp);
}

/*
 * Compute the power of a qfi by an integer.
 *
 * Input:
 *  - f: input qfi
 *  - n: integer exponent
 *  - L: bound for partial reduction in nucomp and nudupl
 * Ouput:
 *  - r: output qfi corresponding to the composition f^n.
 *
 * Algorithm: binary exponentiation with windows of bits
 */
inline
void QFI::nupow (QFI &r, const QFI &f, const Mpz &n, const Mpz &L)
{
  if (n.is_zero ())
  {
    r = ClassGroup (f.discriminant()).one();
  }
  else /* n != 0: binary exponentiation with abs(n) and handle sign after */
  {
    Mpz t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12,
        t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24;
    QFI::OpsAuxVars tmp;

    /* implem of wNAF*: a left to right wNAF (see Brian King, ACNS 2008) */
    const mp_limb_t w = 7;
    const mp_limb_t pow2w = (1UL << w); /* 2^w */
    const mp_limb_t u = (1UL<<(w-2)); /* u = 2^(w-2) */

    /* precomputation: tab[i] = f^(2*i+1)  for 0 <= i < u = 2^(w-2) */
    QFI ff(f), tab[u];

    nudupl (ff, f, L, tmp); /* f^2 */
    tab[0] = f;
    for (mp_limb_t i = 1; i < u; i++)
      nucomp (tab[i], tab[i-1], ff, L, 0, tmp); /* tab[i] <- tab[i-1]*ff */

    int j = n.nbits() - 1;
    mp_limb_t c;

    /* first digit is done out of the main loop */
    {
      /* for the first digit we know that dj=1 and c=0 */
      mp_limb_t m = n.extract_bits ((size_t) j, w);
      c = m & 0x1; /* d_{j-w+1} */
      mp_limb_t t = m + (m & 0x1); /* + d_{j-w+1} */
      size_t val2 = mpn_scan1 (&t, 0); /* note: m cannot be zero */
      size_t tau = val2 < w ? val2 : w-1;
      t >>= tau;

      r = t == 2 ? ff : tab[t>>1];
      size_t b = ((size_t) j) < w-1 ? tau+1+j-w : tau;
      for (size_t i = 0; i < b; i++)
        nudupl (r, r, L, tmp);
      j -= w;
    }

    /* main loop */
    while (j >= 0)
    {
      mp_limb_t m = n.extract_bits ((size_t) j, w);
      mp_limb_t dj = (m >> (w-1)) & 0x1;
      mp_limb_t djmwp1 = m & 0x1;

      if (c == dj)
      {
        nudupl (r, r, L, tmp);
        j -= 1;
      }
      else
      {
        int neg = c;
        mp_limb_t t = m + djmwp1;
        t = c ? (pow2w - t) : t;
        c = djmwp1;

        size_t val2 = t > 0 ? mpn_scan1 (&t, 0) : w-1;
        size_t tau = val2 < w ? val2 : w-1;
        t >>= tau;
        for (size_t i = 0; i < w-tau; i++)
          nudupl (r, r, L, tmp);
        nucomp (r, r, t == 2 ? ff : tab[t>>1], L, neg, tmp);
        size_t b = ((size_t) j) < w-1 ? tau+1+j-w : tau;
        for (size_t i = 0; i < b; i++)
          nudupl (r, r, L, tmp);
        j -= w;
      }
    }

    if (c)
      nucomp (r, r, tab[0], L, 1, tmp);

    /* toggle the result if n is negative */
    if (n.sgn () < 0)
    {
      r.neg();
    }
  }
}

/*
 * Multiple exponentiation.
 *
 * Input:
 *  - f0: input qfi
 *  - f1: input qfi
 *  - n0: integer exponent
 *  - n1: integer exponent
 *  - L: bound for partial reduction in nucomp and nudupl
 * Ouput:
 *  - r: output qfi corresponding to the composition f0^n0*f1^n1.
 *
 * Assumes: f0 and f1 have the same discriminant
 *          n0 and n1 are positive
 */
inline
void QFI::nupow (QFI &r, const QFI &f0, const Mpz &n0,
                               const QFI &f1, const Mpz &n1, const Mpz &L)
{
  if (n0.is_zero() && n1.is_zero())
  {
    r = ClassGroup (f0.discriminant()).one();
  }
  else if (n0.is_zero())
    nupow (r, f1, n1, L);
  else if (n1.is_zero())
    nupow (r, f0, n0, L);
  else /* n0*n1 != 0: use joint sparse form */
  {
    Mpz t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12,
        t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24;
    QFI::OpsAuxVars tmp;

    /* precomputation */
    QFI tab[4];

    tab[0] = f0;
    tab[1] = f1;
    nucomp (tab[2], f0, f1, L, 0, tmp);
    nucomp (tab[3], f0, f1, L, 1, tmp);

    JSF jsf (n0, n1);

    /* init r */
    uint8_t most_significant = jsf[jsf.size()-1];
    if (most_significant == 0x01)
      r = tab[0];
    else if (most_significant == 0x10)
      r = tab[1];
    else if (most_significant == 0x11)
      r = tab[2];
    /* else case not needed (TODO why digit -1 is not possible ??) */

    /* main loop (skipping first nonzero digit) */
    for (size_t j = jsf.size()-1; j > 0; j--)
    {
      uint8_t d = jsf[j-1];

      nudupl (r, r, L, tmp);
      if (d == 0x01) /* f0 */
        nucomp (r, r, tab[0], L, 0, tmp);
      else if (d == 0x03) /* f0^-1 */
        nucomp (r, r, tab[0], L, 1, tmp);
      else if (d == 0x10) /* f1 */
        nucomp (r, r, tab[1], L, 0, tmp);
      else if (d == 0x30) /* f1^-1 */
        nucomp (r, r, tab[1], L, 1, tmp);
      else if (d == 0x11) /* f0 * f1 */
        nucomp (r, r, tab[2], L, 0, tmp);
      else if (d == 0x13) /* f0^-1 * f1 */
        nucomp (r, r, tab[3], L, 1, tmp);
      else if (d == 0x31) /* f0 * f1^-1 */
        nucomp (r, r, tab[3], L, 0, tmp);
      else if (d == 0x33) /* f0^-1 * f1^-1 */
        nucomp (r, r, tab[2], L, 1, tmp);
    }
  }
}

/*
 * Multiple exponentiation.
 *
 * Input:
 *  - f: input qfi
 *  - n: integer exponent
 *  - L: bound for partial reduction in nucomp and nudupl
 *  - d, e: two positive integers such that e < d
 *  - fe: qfi such that fe = f^(2^e)
 *  - fd: qfi such that fe = f^(2^d)
 *  - fed: qfi such that fe = f^(2^(e+d))
 * Ouput:
 *  - r: output qfi corresponding to the composition f^n.
 *
 * Assumes: f0 and f1 have the same discriminant
 *          n0 and n1 are positive
 */
void QFI::nupow (QFI &r, const QFI &f, const Mpz &n, size_t d, size_t e,
                 const QFI &fe, const QFI &fd, const QFI &fed, const Mpz &L)
{
  if (n.is_zero ())
  {
    r = ClassGroup (f.discriminant()).one();
  }
  else if (n.nbits() < e)
  {
    nupow (r, f, n, L);
  }
  else /* n != 0: exponentiation with abs(n) and handle sign after */
  {
    Mpz t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12,
        t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24;
    QFI::OpsAuxVars tmp;

    Mpz::mod2k (t00, n, d);
    Mpz::divby2k (t01, n, d);

    /* */
    JSF jsf (t00, t01);

    /* precomputations */
    /* tab[i] = f^d0 * fd^d1 * fe^d2 * fed^d3
     * where [ d0, d1, d2, d3] is (i+1) written in basis 3 with digits in
     * (-1, 0, 1).
     * Example: tab[20] = f^0 * fd^1 * fe^-1 * fed^1
     *          because 21 = 0*3^0 + 1*3^1 + -1*3^2 + 1*3^3
     */
    QFI tab[40];

    tab[0] = f;
    tab[2] = fd;
    tab[8] = fe;
    tab[26] = fed;

    for (size_t B = 1, pow3 = 3, i = 0; i < 3; i++, B+=pow3, pow3*=3)
    {
      for (size_t k = 0; k < B; k++)
      {
        nucomp (tab[pow3+k],   tab[pow3-1], tab[k], L, 0, tmp);
        nucomp (tab[pow3-k-2], tab[pow3-1], tab[k], L, 1, tmp);
      }
    }

    /* */
    r = ClassGroup (f.discriminant()).one();
    for (size_t j = jsf.size(); j > 2*e; j--)
    {
      uint8_t digh = jsf[j-1];

      nudupl (r, r, L, tmp);

      if (digh != 0)
      {
        /* address decoding */
        int idx = 0;
        idx += (digh & 0x02) ?  -9 : ((digh & 0x01) ? 9 : 0);
        idx += (digh & 0x20) ? -27 : ((digh & 0x10) ? 27 : 0);

        bool s = idx < 0 ? true : false;
        idx = idx < 0 ? -idx : idx;
        nucomp (r, r, tab[idx-1], L, s, tmp);
      }
    }

    for (size_t j = 2*e; j > e; j--)
    {
      uint8_t digh = jsf[j-1];
      uint8_t digl = jsf[j-e-1];

      nudupl (r, r, L, tmp);

      if (digh != 0 || digl != 0)
      {
        /* address decoding */
        int idx = 0;
        idx += (digl & 0x02) ?  -1 : ((digl & 0x01) ?  1 : 0);
        idx += (digl & 0x20) ?  -3 : ((digl & 0x10) ?  3 : 0);
        idx += (digh & 0x02) ?  -9 : ((digh & 0x01) ?  9 : 0);
        idx += (digh & 0x20) ? -27 : ((digh & 0x10) ? 27 : 0);

        bool s = idx < 0 ? true : false;
        idx = idx < 0 ? -idx : idx;
        nucomp (r, r, tab[idx-1], L, s, tmp);
      }
    }
  }
}

/* */
template<>
inline
void OpenSSL::HashAlgo::hash (const QFI &f)
{
  hash (f.a());
  hash (f.b());
  hash (f.c());
}

/******************************************************************************/
/* */
inline
QFICompressedRepresentation::QFICompressedRepresentation (const Mpz &ap,
                                                          const Mpz &g,
                                                          const Mpz &tp,
                                                          const Mpz &b0,
                                                          bool is_neg)
  : ap(ap), g(g), tp(tp), b0(b0), is_neg(is_neg)
{
}

/* */
inline
size_t QFICompressedRepresentation::nbits () const
{
  /* +1 for is_neg which is a boolean */
  return ap.nbits() + g.nbits() + tp.nbits() + b0.nbits() + 1;
}

/* */
inline
std::ostream & operator<< (std::ostream &o,
                           const QFICompressedRepresentation &f)
{
  return o << "(" << f.ap << ", " << f.g << ", " << f.tp << ", " << f.b0 << ", "
           << (f.is_neg ? "1" : "0") << ")";
}

/******************************************************************************/
inline
ClassGroup::ClassGroup (const Mpz &D) : disc_(D), class_number_bound_ (0UL)
{
  if (D.sgn() >= 0)
    throw std::invalid_argument ("the discriminant must be negative");

  unsigned long Dmod4 = D.mod4 ();
  if (Dmod4 != 0 && Dmod4 != 1)
    throw std::invalid_argument ("the discriminant must == 0, 1 mod 4");

  /* Compute the bound for nucomp: floor(|D|^1/4)
   * TODO: Find a ref on why it is better to use |D|^1/4 instead of
   * (|D|/4)^(1/4)
   */
  Mpz::abs (default_nucomp_bound_, disc_);
  Mpz::root4th (default_nucomp_bound_, default_nucomp_bound_);
}

/* */
inline
const Mpz & ClassGroup::discriminant () const
{
  return disc_;
}

/* Return the default bound for NUCOMP */
inline
const Mpz & ClassGroup::default_nucomp_bound () const
{
  return default_nucomp_bound_;
}

/*
 * Return the neutral element of the class group, i.e., the principal qfi
 * for the discriminant of the class group.
 * It is the form [ 1, b, c ] with
 *    - b = disc % 2
 *    - c = 1/4*(b^2-disc)
 * Note that it is already in reduced form.
 */
inline
QFI ClassGroup::one () const
{
  QFI f;

  f.a_ = 1UL;
  f.b_ = disc_.is_odd() ? 1UL : 0UL;
  f.c_ = disc_.is_odd() ? 1UL : 0UL;
  Mpz::sub (f.c_, f.c_, disc_);
  Mpz::divby4 (f.c_, f.c_);

  return f;
}

/*
 * Return one of the two prime forms (l, _, _) of discriminant disc.
 * Assume l is prime and kronecker (disc, l) == 1.
 * Always choose the prime form (l, b, _) where b is the square root of
 * disc modulo 4*l in [1..l-1].
 */
inline
QFI ClassGroup::primeform (const Mpz &l) const
{
  QFI r;
  r.a_ = l;
  Mpz::sqrt_mod_prime (r.b_, disc_, l);
  if (r.b_.is_odd () != disc_.is_odd ()) /* not same parity */
    Mpz::sub (r.b_, l, r.b_);
  r.set_c_from_disc (disc_);
  r.reduction();
  return r;
}

/* /!\ not uniform, use only for tests and benchmarks */
template <size_t ngens, unsigned long coeff_bound>
inline
QFI ClassGroup::random (RandGen &randgen) const
{
  QFI r(one()), fl;
  Mpz m(coeff_bound);

  Mpz l(2UL);

  for (size_t i = 0; i < ngens; i++)
  {
    for ( ; disc_.kronecker (l) != 1; l.nextprime ());
    fl = primeform (l);
    nupow (fl, fl, randgen.random_mpz (m));
    nucomp (r, r, fl);
  }

  return r;
}

/*
 * Return the bound on the class number. The bound B satisfies
 *      B <= class_number < 2*B
 *
 * The bound is cached
 * XXX Do we need the disc to be fundamental for the formula to work ???
 */
inline
const Mpz & ClassGroup::class_number_bound () const
{
  if (class_number_bound_.is_zero()) /* computed it if not already computed */
  {
    Mpz primebound, l, tmp;

    Mpz::ceil_abslog_square (primebound, disc_);

    size_t prec = disc_.nbits() + 50;

    mpf_t acc, t, d;
    mpf_inits (acc, t, d, NULL);
    mpf_set_prec (acc, prec);
    mpf_set_prec (t, prec);
    mpf_set_prec (d, prec);

    mpf_set_ui (acc, 1);

    for (l = 2UL; l < primebound; l.nextprime())
    {
      int k = disc_.kronecker (l);
      if (k < 0)
        Mpz::add (tmp, l, -k);
      else
        Mpz::sub (tmp, l, k);
      mpf_set_z (t, static_cast<mpz_srcptr> (l));
      mpf_set_z (d, static_cast<mpz_srcptr> (tmp));
      mpf_div (t, t, d);
      mpf_mul (acc, acc, t);
    }

    Mpz::abs (tmp, disc_);
    Mpz::sqrt (tmp, tmp); /* tmp <- floor(sqrt(|D|)) */
    Mpz::mul (tmp, tmp, 21);
    mpf_div_ui (acc, acc, 88);
    mpf_set_z (t, static_cast<mpz_srcptr> (tmp));
    mpf_mul (t, t, acc);
    mpf_ceil (t, t);

    class_number_bound_ = t;

    mpf_clears (acc, t, d, NULL);
  }
  return class_number_bound_;
}

inline
void ClassGroup::nucomp (QFI &r, const QFI &f1, const QFI &f2) const
{
  QFI::nucomp (r, f1, f2, default_nucomp_bound(), 0);
}

inline
void ClassGroup::nucompinv (QFI &r, const QFI &f1, const QFI &f2) const
{
  QFI::nucomp (r, f1, f2, default_nucomp_bound(), 1);
}

inline
void ClassGroup::nudupl (QFI &r, const QFI &f) const
{
  QFI::nudupl (r, f, default_nucomp_bound());
}

/* */
inline
void ClassGroup::nudupl (QFI &r, const QFI &f, size_t niter) const
{
  Mpz L = default_nucomp_bound();
  QFI::OpsAuxVars tmp;
  if (niter > 0)
  {
    QFI::nudupl (r, f, L, tmp);
    for (size_t i = 1; i < niter; i++)
      QFI::nudupl (r, r, L, tmp);
  }
}

inline
void ClassGroup::nupow (QFI &r, const QFI &f, const Mpz &n) const
{
  QFI::nupow (r, f, n, default_nucomp_bound());
}

inline
void ClassGroup::nupow (QFI &r, const QFI &f0, const Mpz &n0,
                                   const QFI &f1, const Mpz &n1) const
{
  QFI::nupow (r, f0, n0, f1, n1, default_nucomp_bound());
}

inline
void ClassGroup::nupow (QFI &r, const QFI &f, const Mpz &n, size_t d, size_t e,
                 const QFI &fe, const QFI &fd, const QFI &fed) const
{
  QFI::nupow (r, f, n, d, e, fe, fd, fed, default_nucomp_bound());
}

#endif /* BICYCL_QFI_INL */
