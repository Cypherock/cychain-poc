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
bool
test_ECNIZK (SecLevel seclevel, size_t niter)
{
  std::stringstream pre;
  pre << "security " << seclevel << " bits ECNIZK";

  bool ret = true;

  OpenSSL::ECGroup E (seclevel);
  OpenSSL::HashAlgo H (seclevel);

  for (size_t i = 0; i < niter; i++)
  {
    ECNIZKProof::SecretValue s (E.random_mod_order());
    ECNIZKProof::PublicValue Q (E, s);

    ECNIZKProof proof (E, H, s, Q);

    ret &= proof.verify (E, H, Q);
  }

  Test::result_line (pre.str(), ret);
  return ret;
}

/******************************************************************************/
bool check (SecLevel seclevel, size_t niter)
{
  bool success = true;

  std::stringstream desc;
  desc << "security " << seclevel << " bits";

  ECDSA C (seclevel);
  desc << " ECDSA";
  success &= Test::test_sign (C, niter, desc.str());

  success &= test_ECNIZK (seclevel, niter);

  return success;
}

/******************************************************************************/
int
main (int argc, char *argv[])
{
  bool success = true;

  RandGen randgen;
  randseed_from_argv (randgen, argc, argv);

  Test::OverrideOpenSSLRand::WithRandGen tmp_override (randgen);

  success &= check (SecLevel::_112, 50);
  success &= check (SecLevel::_128, 50);
  success &= check (SecLevel::_192, 50);
  success &= check (SecLevel::_256, 50);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
