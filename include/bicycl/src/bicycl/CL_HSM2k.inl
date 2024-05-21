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
#ifndef BICYCL_CL_HSM2K_INL
#define BICYCL_CL_HSM2K_INL

/**
 * \param[in] q the prime q
 * \param[in] p the prime p or 1
 * \param[in] distance positive integer
 * \param[in] compact_variant ""
 *
 */
inline
CL_HSM2k::CL_HSM2k (const Mpz &N, size_t k, unsigned int distance,
                    bool compact_variant)
  : N_(N),
    k_(k),
    Cl_DeltaK_ (compute_DeltaK (N)),
    Cl_Delta_ (compute_Delta (Cl_DeltaK_.discriminant(), k)),
    compact_variant_ (compact_variant),
    distance_ (distance)
{
  /* Checks */
  if (N_.sgn() <= 0)
    throw std::invalid_argument ("N must be positive");
  if (k_ == 0)
    throw std::invalid_argument ("k must be positive");

  /* Compute M = 2^k */
  M_ = 1UL;
  Mpz::mulby2k (M_, M_, k);

  /* Assess if we need the large message variant, i.e., if 2^(2k) > 1 - DeltaK
   */
  Mpz t;
  Mpz::mul (t, M_, M_);       /* t <- M^2 = 2^(2*k) */
  Mpz::sub (t, t, 1UL);       /* t <- 2^(2*k) - 1 */
  Mpz::add (t, t, DeltaK());  /* t <- 2^(2*k) - 1 + DeltaK  */
  large_message_variant_ = (t.sgn() > 0);

  /* Compute the generator h
   * For the non compact variant, the generator is the square of the
   * smallest primeform of Cl_Delta raised to the power M=2^k.
   * For the compact variant, we push it into Cl_DeltaK and raise it to the
   * power M=2^k.
   */
  /* smallest primeform */
  Mpz l(2UL);
  for ( ; Delta().kronecker (l) != 1; l.nextprime ());
  h_ = Cl_Delta_.primeform (l);

  /* square it */
  Cl_Delta_.nudupl (h_, h_);

  /* raise it to power M=2^k */
  raise_to_power_M (Cl_Delta_, h_);

  if (compact_variant) /* For compact variant, we need \pi(h)^M */
  {
    h_.to_maximal_order_2exp (k_+1, DeltaK());
    raise_to_power_M (Cl_DeltaK_, h_);
  }

  /* Compute the exponent_bound as class_number_bound times 2^(distance-2).
   * If distance is < 2, the default it to use distance = 42.
   */
  exponent_bound_ = Cl_DeltaK_.class_number_bound();
  distance_ = distance_ < 2 ? 42 : distance_;
  Mpz::mulby2k (exponent_bound_, exponent_bound_, distance_-2);

  /*
   * Precomputation
   */
  d_ = (encrypt_randomness_bound().nbits () + 1)/2;
  e_ = d_/2 + 1;
  h_de_precomp_ = h_;
  for (size_t i = 0; i < d_+e_; i++)
  {
    if (i == e_)
      h_e_precomp_ = h_de_precomp_;
    if (i == d_)
      h_d_precomp_ = h_de_precomp_;
    Cl_G().nudupl (h_de_precomp_, h_de_precomp_);
  }
}

/**
 */
inline
CL_HSM2k::CL_HSM2k (const Mpz &N, size_t k, unsigned int distance)
  : CL_HSM2k (N, k, distance, false)
{
}

/**
 */
inline
CL_HSM2k::CL_HSM2k (const Mpz &N, size_t k, bool compact_variant)
  : CL_HSM2k (N, k, 0U, compact_variant)
{
}

/**
 */
inline
CL_HSM2k::CL_HSM2k (const Mpz &N, size_t k)
  : CL_HSM2k (N, k, 0U)
{
}

/**
 */
inline
CL_HSM2k::CL_HSM2k (const CL_HSM2k &C, bool compact_variant)
  : CL_HSM2k(C)
{
  if (compact_variant != C.compact_variant_)
  {
    compact_variant_ = compact_variant;
    if (compact_variant)
    {
      /* we go from non compact to compact variant, we need to compute
       * \pi(h)^M
       */
      h_.to_maximal_order_2exp (k_+1, DeltaK());
      raise_to_power_M (Cl_DeltaK_, h_);
    }
    else
    {
      /* smallest primeform */
      Mpz l(2UL);
      for ( ; Delta().kronecker (l) != 1; l.nextprime ());
      h_ = Cl_Delta_.primeform (l);

      /* square it */
      Cl_Delta_.nudupl (h_, h_);

      /* raise it to power M=2^k */
      raise_to_power_M (Cl_Delta_, h_);
    }

    /* precomputations */
    h_de_precomp_ = h_;
    for (size_t i = 0; i < d_+e_; i++)
    {
      if (i == e_)
        h_e_precomp_ = h_de_precomp_;
      if (i == d_)
        h_d_precomp_ = h_de_precomp_;
      Cl_G().nudupl (h_de_precomp_, h_de_precomp_);
    }
  }
}

/**
 * @private
 */
template <class... T>
inline
CL_HSM2k::CL_HSM2k (size_t N_nbits, size_t k, RandGen &randgen, T... args)
  : CL_HSM2k (random_N (randgen, N_nbits), k, args...)
{
}

/**
 * @private
 */
template <class... T>
inline
CL_HSM2k::CL_HSM2k (SecLevel seclevel, size_t k, RandGen &randgen, T... args)
  : CL_HSM2k (seclevel.RSA_modulus_bitsize(), k, randgen, args...)
{
}

/* */
inline
bool CL_HSM2k::large_message_variant () const
{
  return large_message_variant_;
}

/* */
inline
std::ostream & operator<< (std::ostream &o, const CL_HSM2k &C)
{
  return o << "N = " << C.N() << " # " << C.N().nbits() << " bits" << std::endl
           << "k = " << C.k() << std::endl
           << "DeltaK = -8*N # " << C.DeltaK().nbits() << " bits" << std::endl
           << "Delta = 2^(2*(k+1))*DeltaK # " << C.Delta().nbits() << " bits"
           << std::endl
           << "h = " << C.h() << std::endl
           << "compact_variant = " << C.compact_variant() << std::endl
           << "large_message_variant = " << C.large_message_variant()
           << std::endl;
}

/* */
inline
const Mpz & CL_HSM2k::N () const
{
  return N_;
}

/* */
inline
size_t CL_HSM2k::k () const
{
  return k_;
}

/* */
inline
const Mpz & CL_HSM2k::M () const
{
  return M_;
}


/* */
inline
const ClassGroup & CL_HSM2k::Cl_DeltaK () const
{
  return Cl_DeltaK_;
}

/* */
inline
const ClassGroup & CL_HSM2k::Cl_Delta () const
{
  return Cl_Delta_;
}

/* */
inline
const ClassGroup & CL_HSM2k::Cl_G () const
{
  return compact_variant_ ? Cl_DeltaK_ : Cl_Delta_;
}

/* */
inline
const Mpz & CL_HSM2k::DeltaK () const
{
  return Cl_DeltaK_.discriminant();
}

/* */
inline
const Mpz & CL_HSM2k::Delta () const
{
  return Cl_Delta_.discriminant();
}

/* */
inline
const QFI & CL_HSM2k::h () const
{
  return h_;
}

/* */
inline
bool CL_HSM2k::compact_variant () const
{
  return compact_variant_;
}

/* */
inline
const Mpz & CL_HSM2k::secretkey_bound () const
{
  return exponent_bound_;
}

/* */
inline
const Mpz & CL_HSM2k::cleartext_bound () const
{
  return M_;
}

/* */
inline
const Mpz & CL_HSM2k::encrypt_randomness_bound () const
{
  return exponent_bound_;
}

/* */
inline
unsigned int CL_HSM2k::lambda_distance () const
{
  return distance_;
}

/**
 * \param[out] r the quadratic form corresponding to #gen to the power of \p e
 * \param[in] e the exponent
 */
inline
void CL_HSM2k::power_of_h (QFI &r, const Mpz &n) const
{
  Cl_G().nupow (r, h_, n, d_, e_, h_e_precomp_, h_d_precomp_, h_de_precomp_);
}

/*
 */
inline
QFI CL_HSM2k::power_of_f (const Mpz &m) const
{
  const size_t val2 = m.val2();
  if (val2 >= k_) /* m == 0 mod 2^k */
  {
    return Cl_Delta_.one();
  }
  else /* m != 0: compute Lucas chains U_m and V_m */
  {
    Mpz m00(1UL), m01(0UL);
    Mpz m10(0UL), m11(1UL);
    Mpz minusQ, t0, t1, t2, t3;

    Mpz::sub (minusQ, DeltaK(), 1UL);

    for (size_t i = k_; i > 0; i--)
    {
      /* square */
      Mpz::mul (t0, m00, m00);
      Mpz::addmul (t0, m01, m10);
      Mpz::mod2k (t0, t0, k_);

      Mpz::mul (t1, m00, m01);
      Mpz::addmul (t1, m01, m11);
      Mpz::mod2k (t1, t1, k_);

      Mpz::mul (t2, m10, m00);
      Mpz::addmul (t2, m11, m10);
      Mpz::mod2k (t2, t2, k_);

      Mpz::mul (t3, m10, m01);
      Mpz::addmul (t3, m11, m11);
      Mpz::mod2k (t3, t3, k_);

      Mpz::swap (m00, t0);
      Mpz::swap (m01, t1);
      Mpz::swap (m10, t2);
      Mpz::swap (m11, t3);

      /* mul */
      if (m.tstbit (i-1))
      {
        Mpz::mul (m00, m00, minusQ);
        Mpz::add (m00, m00, m10);
        Mpz::add (m00, m00, m10);
        Mpz::mod2k (m00, m00, k_);

        Mpz::mul (m01, m01, minusQ);
        Mpz::add (m01, m01, m11);
        Mpz::add (m01, m01, m11);
        Mpz::mod2k (m01, m01, k_);

        Mpz::swap (m00, m10);
        Mpz::swap (m01, m11);
      }
    }

    /* Vn = 2*(m00+m01) => Vn/2 = m00+m01 */
    Mpz::add (t0, m00, m01);
    /* Un = m01 */
    Mpz::divby2k (t1, m01, val2);

    Mpz::mod_inverse_2k (t2, t1, k_-val2, t3);
    Mpz::mul (t0, t0, t2);
    Mpz::mod2k_centered (t0, t0, k_-val2);

    /* a <- 2^(2*(k-val2))    [ stored in t1 ] */
    t1 = 1UL;
    Mpz::mulby2k (t1, t1, 2*(k_-val2));
    /* b <- 2^(k-val2+1) * u  [ stored in t2 ] */
    Mpz::mulby2k (t2, t0, k_-val2+1);
    /* c <- u^2 - 2^(2*val2)*Delta_K = u^2 + 2^(2*val2+3)*N  [ stored in t3 ] */
    Mpz::mulby2k (t3, N_, 2*val2+3);
    Mpz::addmul (t3, t0, t0);

    /* No need to check the form (a,b,c) */
    return QFI (t1, t2, t3, true);
  }
}

/* Assume fm is in F */
inline
Mpz CL_HSM2k::dlog_in_F (const QFI &fm) const
{
  Mpz m(0UL);
  /* Compute log(fm) for the basis f
   * The computation is done using representative in
   * (OK / 2^(k+1) OK)* / (Z / 2^(k+1) Z)*
   * A form in f corresponds to a 1+t*sqrt(DeltaK) with 0 <= t < 2^k.
   * In particular:
   *  - f corresponds to 1+sqrt(DeltaK)
   *  - fm corresponds to (1+tm*sqrt(DeltaK)) with tm = (1/u % 2^j) * 2^(k-j)
   */
  if (!fm.is_one())
  {
    Mpz tm, aux0, aux1, aux2;

    if (large_message_variant_)
    {
      tm = fm.kernel_representative_2exp (k_+1, DeltaK());
    }
    else
    {
      /* compute tm: reading j, u from b coeff */
      size_t j = fm.b().val2() - 1;
      Mpz::divby2k (aux0, fm.b(), j+1); /* aux0 <- u */
      Mpz::mod_inverse_2k (tm, aux0, j, aux1);
      Mpz::mulby2k (tm, tm, k_-j);
    }
    Mpz::mod2k (tm, tm, k_);

    Mpz t(1UL);

    for (size_t i = 0; i < k_; i++)
    {
      /* loop invariant: val2(t) = i <=> order(f) = 2^(k-i) */
      size_t val2 = tm.sgn() ? tm.val2(): SIZE_MAX;
      if (val2 == i) /* do fm and f have the same order */
      {
        m.setbit (i);
        /* (1+tm*X) = (1+tm*X)/(1+t*X) */
        F_kerphi_div (tm, t, aux0, aux1, aux2);
      }
      F_kerphi_square (t, aux0, aux1, aux2); /* (1+t*X) = (1+t*X)^2 */
    }
  }
  /* else: m is already set to the correct value: 0 */
  return m;
}

/* */
inline
void CL_HSM2k::from_Cl_DeltaK_to_Cl_Delta (QFI &f) const
{
  f.lift_2exp (k_+1);
  int chi_m4 = (f.a().mod4() == 3 || f.c().mod4() == 3) ? -1 : 1;
  raise_to_power_M (Cl_Delta(), f);
  if (chi_m4 == -1)
  {
    Mpz t(DeltaK());
    Mpz::mulby2k (t, t, 2*(k_-1));
    Mpz::sub (t, t, 1UL);
    t.neg();
    Mpz four (4UL);
    QFI g (four, four, t);
    Cl_Delta().nucomp (f, f, g);
  }
}

/* */
inline
CL_HSM2k::SecretKey CL_HSM2k::keygen (RandGen &randgen) const
{
  return SecretKey (*this, randgen);
}

inline
CL_HSM2k::PublicKey CL_HSM2k::keygen (const SecretKey &sk) const
{
  return PublicKey (*this, sk);
}

/*
 * Encrypt the plaintext using the cryptosystems described by params, the
 * public key pk and the randomness r.
 *
 * Input:
 *  params: the parameters of the cryptosystems
 *  pk: the public key
 *  m: the plaintext to encrypt
 *  r: randomness
 */
CL_HSM2k::CipherText CL_HSM2k::encrypt (const PublicKey &pk,
                                        const ClearText &m, const Mpz &r) const
{
  return CipherText (*this, pk, m, r);
}


/*
 * Same as above but without the randomness
 */
CL_HSM2k::CipherText CL_HSM2k::encrypt (const PublicKey &pk, const ClearText &m,
                                        RandGen &randgen) const
{
  return encrypt (pk, m, randgen.random_mpz (encrypt_randomness_bound()));
}

/*
 * Decrypt the ciphertext c using the cryptosystems described by params and
 * the secret key sk
 *
 * Input:
 *  sk: the secret key
 *  c: the ciphertext
 */
CL_HSM2k::ClearText CL_HSM2k::decrypt (const SecretKey &sk, const CipherText &c)
                                                                          const
{
  return ClearText (*this, sk, c);
}

/* */
inline
CL_HSM2k::CipherText CL_HSM2k::add_ciphertexts (const PublicKey &pk,
                                                const CipherText &ca,
                                                const CipherText &cb,
                                                RandGen &randgen) const
{
  Mpz r(randgen.random_mpz (encrypt_randomness_bound()));
  return add_ciphertexts (pk, ca, cb, r);
}

/* */
inline
CL_HSM2k::CipherText CL_HSM2k::add_ciphertexts (const PublicKey &pk,
                                                const CipherText &ca,
                                                const CipherText &cb,
                                                const Mpz &r) const
{
  return CipherText (*this, pk, ca, cb, r);
}

/* */
inline
CL_HSM2k::ClearText CL_HSM2k::add_cleartexts (const ClearText &ma,
                                              const ClearText &mb) const
{
  return ClearText (*this, ma, mb);
}

/* */
inline
CL_HSM2k::CipherText CL_HSM2k::scal_ciphertexts (const PublicKey &pk,
                                                 const CipherText &c,
                                                 const Mpz &s,
                                                 RandGen &randgen) const
{
  Mpz r(randgen.random_mpz (encrypt_randomness_bound()));
  return scal_ciphertexts (pk, c, s, r);
}

/* */
inline
CL_HSM2k::CipherText CL_HSM2k::scal_ciphertexts (const PublicKey &pk,
                                                 const CipherText &c,
                                                 const Mpz &s,
                                                 const Mpz &r) const
{
  return CipherText (*this, pk, c, s, r);
}

/* */
inline
CL_HSM2k::ClearText CL_HSM2k::scal_cleartexts (const ClearText &m,
                                               const Mpz &s) const
{
  return ClearText (*this, m, s);
}

/*
 */
inline
Mpz CL_HSM2k::random_N (RandGen &randgen, size_t N_nbits)
{
  Mpz p, q, N;

  /* We generate prime p, q such that (p%8, q%8) is (3,5) or (5,3).
   * It corresponds to one of the valid case of Table 1 of CL22.
   * These cases were chosen as they do not require computing kronecker symbol.
   */

  /* Generate a random prime p satisfying the conditions */
  p = randgen.random_prime ((N_nbits+1)/2);
  while (p.mod8() != 3 && p.mod8() != 5)
  {
    p.nextprime ();
  }

  /* Generate a random prime q satisfying the conditions */
  q = randgen.random_prime (N_nbits/2);
  unsigned long q_mod_8 = p.mod8() == 3 ? 5 : 3;
  while (q.mod8() != q_mod_8)
  {
    q.nextprime ();
  }

  Mpz::mul (N, p, q);
  return N;
}

/* Compute Delta_OK = -8N
 */
inline
Mpz CL_HSM2k::compute_DeltaK (const Mpz &N)
{
  Mpz d;
  Mpz::mulby2k (d, N, 3);
  d.neg ();
  return d;
}

/* Compute Delta_O2k = 2^(2*(k+1)) * Delta_OK
 */
inline
Mpz CL_HSM2k::compute_Delta (const Mpz &DeltaK, size_t k)
{
  Mpz d;
  Mpz::mulby2k (d, DeltaK, 2*(k+1));
  return d;
}

/* */
inline
void CL_HSM2k::raise_to_power_M (const ClassGroup &Cl, QFI &f) const
{
  // TODO allocate tmp var only once
  for (size_t i = 0; i < k_; i++)
  {
    Cl.nudupl (f, f);
  }
}

/*
 * Compute (1+t*sqrt(DeltaK))^2 in (OK/2^(k+1) OK)* / (Z/q^(k+1) Z)*
 */
inline
void CL_HSM2k::F_kerphi_square (Mpz &t, Mpz &aux0, Mpz &aux1, Mpz &aux2) const
{
  /* t <- 2t / (1+t^2*DeltaK) mod 2^k */
  Mpz::mul (aux0, t, t);                      /* aux0 <- t^2 */
  Mpz::mod2k (aux0, aux0, k_);                /* aux0 <- t^2 mod 2^k */
  Mpz::mul (aux0, aux0, DeltaK());            /* aux0 <- t^2*DeltaK */
  Mpz::add (aux0, aux0, 1UL);                 /* aux0 <- 1+t^2*DeltaK */
  Mpz::mod2k (aux0, aux0, k_);                /* aux0 <- 1+t^2*DeltaK mod 2k */
  Mpz::mod_inverse_2k (aux2, aux0, k_, aux1); /* aux2 <- aux0^(-1) mod 2^k */
  Mpz::mul (t, t, aux2);                      /* t <- t * aux2 */
  Mpz::mulby2k (t, t, 1);                     /* t <- 2 * t * aux2 */
  Mpz::mod2k (t, t, k_);                      /* t <- 2 * t * aux2 mod 2^k */
}

/*
 * Compute (1+s*sqrt(DeltaK))/(1+t*sqrt(DeltaK) in
 * (OK/2^(k+1) OK)* / (Z/q^(k+1) Z)*
 */
inline
void CL_HSM2k::F_kerphi_div (Mpz &t, const Mpz &s, Mpz &aux0, Mpz &aux1,
                                                              Mpz &aux2) const
{
  /* t <- (t-s) / (1-t*s*DeltaK) mod 2^k */
  Mpz::mul (aux0, t, s);                      /* aux0 <- t*s */
  Mpz::mod2k (aux0, aux0, k_);                /* aux0 <- t*s mod 2^k */
  Mpz::mul (aux0, aux0, DeltaK());            /* aux0 <- t*s*DeltaK */
  Mpz::sub (aux0, aux0, 1UL);                 /* aux0 <- t*s*DeltaK - 1 */
  aux0.neg();                                 /* aux0 <- 1-t*s*DeltaK */
  Mpz::mod2k (aux0, aux0, k_);                /* aux0 <- 1-t*s*DeltaK mod 2^k */
  Mpz::mod_inverse_2k (aux2, aux0, k_, aux1); /* aux2 <- aux0^(-1) mod 2^k */
  Mpz::sub (t, t, s);                         /* t <- t-s */
  Mpz::mul (t, t, aux2);                      /* t <- (t-s)*aux0^(-1) */
  Mpz::mod2k (t, t, k_);                      /* t <- (t-s)*aux0^(-1) mod 2^k */
}

#endif /* BICYCL_CL_HSM2K_INL */
