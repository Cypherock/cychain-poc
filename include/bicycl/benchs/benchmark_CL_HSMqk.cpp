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
#include <string>
#include <sstream>

#include "bicycl.hpp"
#include "internals.hpp"

using std::string;

using namespace BICYCL;
using BICYCL::Bench::ms;
using BICYCL::Bench::us;

/* Needed to make protected methods public */
class BenchQFI : public QFI
{
  public:
    using QFI::nudupl;
    using QFI::nucomp;
};

/* */
template <class Cryptosystem>
void CL_HSMqk_benchs_ops (const Cryptosystem &C, RandGen &randgen,
                          const string &pre)
{
  size_t niter = 100000;

  QFI::OpsAuxVars tmp;

  QFI r, r2, f;

  f = C.Cl_G().random (randgen);
  const Mpz &bound = C.Cl_G().default_nucomp_bound();

  /* nudupl */
  auto nudupl = [ &r, &f, &bound ] ()
  {
    BenchQFI::nudupl (r, f, bound);
  };
  Bench::one_function<ms, us> (nudupl, niter, "nudupl", pre);

  /* nudupl with temporaries */
  auto nudupl_tmp = [ &r, &f, &bound, &tmp ] ()
  {
    BenchQFI::nudupl (r, f, bound, tmp);
  };
  Bench::one_function<ms, us> (nudupl_tmp, niter, "nudupl w/ tmps", pre);

  /* nucomp */
  auto nucomp = [ &r2, &f, &r, &bound ] ()
  {
    BenchQFI::nucomp (r2, f, r, bound, 0);
  };
  Bench::one_function<ms, us> (nucomp, niter, "nucomp", pre);

  /* nucomp using temporaries */
  auto nucomp_tmp = [ &r2, &f, &r, &bound, &tmp ] ()
  {
    BenchQFI::nucomp (r2, f, r, bound, 0, tmp);
  };
  Bench::one_function<ms, us> (nucomp_tmp, niter, "nucomp w/ tmps", pre);

  /* lift (for compact variant) */
  if (C.compact_variant())
  {
    auto lift = [ &r, &f, C ] ()
    {
      r = f;
      C.from_Cl_DeltaK_to_Cl_Delta (r);
    };
    Bench::one_function<ms, ms> (lift, niter/100, "lift", pre);
  }

  /* Compressed representation */
  auto compress = [ &f ] ()
  {
    f.compressed_repr ();
  };
  Bench::one_function<ms, us> (compress, niter, "compressed_repr", pre);
}

/* */
template <class Cryptosystem>
void CL_HSMqk_benchs_encrypt (const Cryptosystem &C, RandGen &randgen,
                              const string &pre)
{
  const size_t niter = 100;

  typename Cryptosystem::SecretKey sk = C.keygen (randgen);
  typename Cryptosystem::PublicKey pk = C.keygen (sk);
  typename Cryptosystem::ClearText m (C, randgen);
  Mpz r (randgen.random_mpz (C.encrypt_randomness_bound()));
  typename Cryptosystem::CipherText c (C.encrypt (pk, m, r));

  /* */
  auto keygen = [ &C, &randgen ] ()
  {
    typename Cryptosystem::SecretKey sk = C.keygen (randgen);
    C.keygen (sk);
  };
  Bench::one_function<ms, ms> (keygen, niter, "keygen", pre);

  /* */
  auto encrypt = [&C, &pk, &m, &r] ()
  {
    auto c = C.encrypt (pk, m, r);
    c.c1().compressed_repr();
    c.c2().compressed_repr();
  };
  Bench::one_function<ms, ms> (encrypt, niter, "encrypt", pre);

  /* */
  auto decrypt = [&C, &sk, &c] () {C.decrypt (sk, c);};
  Bench::one_function<ms, ms> (decrypt, niter, "decrypt", pre);
}

/* */
void CL_HSMqk_ZKAoK_benchs (const CL_HSMqk &cryptosystem, OpenSSL::HashAlgo &H,
                            RandGen &randgen, const string &pre)
{
  const size_t niter = 100;

  CL_HSMqk::SecretKey sk = cryptosystem.keygen (randgen);
  CL_HSMqk::PublicKey pk = cryptosystem.keygen (sk);

  CL_HSMqk::ClearText a (cryptosystem, randgen);
  Mpz r (randgen.random_mpz (cryptosystem.encrypt_randomness_bound()));
  CL_HSMqk::CipherText c (cryptosystem.encrypt (pk, a, r));
  CL_HSMqk_ZKAoKProof p (cryptosystem, H, pk, c, a, r, randgen);

  auto proof = [&cryptosystem, &H, &pk, &c, &a, &r, &randgen] ()
  {
    CL_HSMqk_ZKAoKProof p (cryptosystem, H, pk, c, a, r, randgen);
  };
  Bench::one_function<ms, ms> (proof, niter, "ZKAoK proof", pre);

  auto verify = [&p, &cryptosystem, &H, &pk, &c] ()
  {
    p.verify (cryptosystem, H, pk, c);
  };
  Bench::one_function<ms, ms> (verify, niter, "ZKAoK verify", pre);
}

/* */
template <class Cryptosystem>
void CL_HSMqk_benchs_all (const Cryptosystem &C, OpenSSL::HashAlgo &H,
                          RandGen &randgen, const string &pre)
{
  CL_HSMqk_benchs_ops (C, randgen, pre);
  CL_HSMqk_benchs_encrypt (C, randgen, pre);
  CL_HSMqk_ZKAoK_benchs (C, H, randgen, pre);
}

/* */
int
main (int argc, char *argv[])
{
  RandGen randgen;
  randseed_from_argv (randgen, argc, argv);

  std::cout << "security level | q_nbits | variant |    operation    | "
            << "  time   | time per op. " << std::endl;

  for (const SecLevel seclevel: SecLevel::All())
  {
    /* With a random q twice as big as the security level */
    std::stringstream desc;
    OpenSSL::HashAlgo H (seclevel);
    CL_HSMqk C (2*seclevel, 1, seclevel, randgen, false);
    desc << "  " << seclevel << " bits     |   " << C.q().nbits() << "   |";

    std::string pre = desc.str() + "  none  ";

    /* */
    auto setup = [ &randgen, &seclevel ] ()
    {
      CL_HSMqk C (2*seclevel, 1, seclevel, randgen, false);
    };
    Bench::one_function<ms, ms> (setup, 10, "setup", pre);

    /* */
    CL_HSMqk_benchs_all (C, H, randgen, pre);
    std::cout << std::endl;

    /* */
    CL_HSMqk Ccompact (C.q(), C.k(), C.p(), true);
    CL_HSMqk_benchs_all (Ccompact, H, randgen, desc.str() + " compact");
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
