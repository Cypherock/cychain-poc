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

/* */
bool Paillier_test_encryption (const Paillier &C, RandGen &randgen,
                               size_t niter, const string & pre)
{
  int ret = 1;

  Paillier::SecretKey sk = C.keygen (randgen);
  Paillier::PublicKey pk = C.keygen (sk);

  /* Test niter random encryption + decryption */
  for (size_t i = 0; i < niter; i++)
  {
    Paillier::ClearText m (C, pk, randgen);
    Paillier::CipherText c = C.encrypt (pk, m, randgen);
    Paillier::ClearText t = C.decrypt (pk, sk, c);

    ret &= (m == t);
  }

  /* Test encryption + decryption of 0 */
  Paillier::ClearText m (C, pk, Mpz (0UL));
  Paillier::CipherText c (C.encrypt (pk, m, randgen));
  Paillier::ClearText t (C.decrypt (pk, sk, c));
  ret &= (m == t);

  Test::result_line (pre + " Paillier encrypt/decrypt", ret);
  return ret;
}

/******************************************************************************/
bool check (RandGen &randgen, size_t N_nbits, size_t niter, const string &name)
{
  bool success = true;
  Paillier C (N_nbits);
  success &= Paillier_test_encryption (C, randgen, niter, name);

  CamenischShoup CS (N_nbits, randgen);
  success &= Test::test_encryption (CS, randgen, niter,
                                                  name + " Camenisch-Shoup");

  return success;
}

bool check (RandGen &randgen, SecLevel seclevel, size_t niter)
{
  bool success = true;

  std::stringstream desc;
  desc << "security " << seclevel << " bits";

  Paillier C (seclevel);
  success &= Paillier_test_encryption (C, randgen, niter, desc.str());

  CamenischShoup CS (seclevel, randgen);
  desc << " Camenisch-Shoup";
  success &= Test::test_encryption (CS, randgen, niter, desc.str());

  return success;
}

/******************************************************************************/
int
main (int argc, char *argv[])
{
  bool success = true;

  RandGen randgen;
  randseed_from_argv (randgen, argc, argv);

  success &= check (randgen, 200, 500, "small");
  success &= check (randgen, SecLevel::_112, 25);
  success &= check (randgen, SecLevel::_128, 20);
  success &= check (randgen, SecLevel::_192, 5);
  //success &= check (randgen, SecLevel::_256, 1);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
