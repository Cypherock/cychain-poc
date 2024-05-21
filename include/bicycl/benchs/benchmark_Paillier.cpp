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

/* */
template <class Cryptosystem>
void Paillier_benchs_encrypt (const Cryptosystem &C, RandGen &randgen,
                         const string &pre)
{
  const size_t niter = 100;

  typename Cryptosystem::SecretKey sk = C.keygen (randgen);
  typename Cryptosystem::PublicKey pk = C.keygen (sk);
  typename Cryptosystem::ClearText m (C, pk, randgen);
  Mpz r (randgen.random_mpz (C.encrypt_randomness_bound(pk)));
  typename Cryptosystem::CipherText c (C.encrypt (pk, m, r));

  /* */
  auto keygen = [ &C, &randgen ] ()
  {
    typename Cryptosystem::SecretKey sk = C.keygen (randgen);
    C.keygen (sk);
  };
  Bench::one_function<ms, ms> (keygen, niter/10, "keygen", pre);

  auto encrypt = [&C, &pk, &m, &r] () {C.encrypt (pk, m, r);};
  Bench::one_function<ms, ms> (encrypt, niter, "encrypt", pre);

  auto decrypt = [&C, &pk, &sk, &c] () {C.decrypt (pk, sk, c);};
  Bench::one_function<ms, ms> (decrypt, niter, "decrypt", pre);
}

/* */
template <class Cryptosystem>
void CamenischShoup_benchs_encrypt (const Cryptosystem &C, RandGen &randgen,
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

  auto encrypt = [&C, &pk, &m, &r] () {C.encrypt (pk, m, r);};
  Bench::one_function<ms, ms> (encrypt, niter, "encrypt", pre);

  auto decrypt = [&C, &sk, &c] () {C.decrypt (sk, c);};
  Bench::one_function<ms, ms> (decrypt, niter, "decrypt", pre);
}

/* */
int
main (int argc, char *argv[])
{
  RandGen randgen;
  randseed_from_argv (randgen, argc, argv);

  std::cout << "                | sec. | N_nbits | operation      "
            << " |   time   | time per op. " << std::endl;

  for (const SecLevel seclevel: SecLevel::All())
  {
    /* Paillier */
    std::cout << std::endl;
    {
      Paillier C (seclevel);
      std::stringstream desc;
      desc << "Paillier        |  " << seclevel << " | "
           << std::setw (7) << C.N_nbits();

      /* */
      auto setup = [ &seclevel ] ()
      {
        Paillier C (seclevel);
      };
      Bench::one_function<ms, ms> (setup, 10, "setup", desc.str());

      /* */
      Paillier_benchs_encrypt (C, randgen, desc.str());
    }

    /* Camenisch-Shoup */
    std::cout << std::endl;
    {
      CamenischShoup CS (seclevel, randgen);
      std::stringstream desc;
      desc << "Camenisch-Shoup |  " << seclevel << " | "
           << std::setw (7) << CS.N().nbits();

      /* */
      auto setup = [ &randgen, &seclevel ] ()
      {
        CamenischShoup C (seclevel, randgen);
      };
      Bench::one_function<ms, ms> (setup, 10, "setup", desc.str());

      /* */
      CamenischShoup_benchs_encrypt (CS, randgen, desc.str());
    }
  }


  return EXIT_SUCCESS;
}
