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

string generate_desc (const CL_HSM2k &C, const string &name)
{
  std::stringstream desc;
  desc << name << " | " << std::setw(4) << C.k()
               << " | " << std::setw(6) << C.DeltaK().nbits()
               << " | " << (C.large_message_variant() ? " true" : "false")
               << " | " << (C.compact_variant() ? "  true " : " false ")
               << " |";
  return desc.str();
}

bool CL_HSM2k_check_setup (const CL_HSM2k &C, const string &pre)
{
  bool ret = 1;
  Mpz tmp;

  /* As p and q are not given, we cannot chech validity according to Table 1
   * of CL22.
   */

  /* Check Delta_K == -8*N */
  Mpz::mulby2k (tmp, C.N(), 3);
  tmp.neg();
  ret &= (tmp == C.DeltaK());

  Mpz::mulby2k (tmp, C.DeltaK(), 2*(C.k()+1));
  ret &= (tmp == C.Delta());

  ret &= C.h().discriminant() == C.Cl_G().discriminant();

  // TODO check generator is a square
  // TODO check precomputation values

  Test::result_line (pre + " setup", ret);
  return ret;
}

/* */
bool CL_HSM2k_check_all (const CL_HSM2k &C, RandGen &randgen, size_t niter,
                         const string &name)
{
  bool success = true;
  std::cout << std::endl;
  success &= CL_HSM2k_check_setup (C, name);
  success &= Test::test_encryption (C, randgen, niter, name);
  success &= Test::test_ciphertext_ops (C, randgen, niter, name);
  return success;
}

/* */
bool check (const CL_HSM2k &C, RandGen &randgen, size_t niter,
                                                 const string &name)
{
  bool success = true;

  success &= CL_HSM2k_check_all (C, randgen, niter, generate_desc (C, name));

  CL_HSM2k Ccompact (C, true);
  success &= CL_HSM2k_check_all (Ccompact, randgen, niter,
                                                generate_desc (Ccompact, name));

  return success;
}

/* */
bool check (RandGen &randgen, size_t N_nbits, size_t k, size_t niter,
            const string &name)
{
  CL_HSM2k C (N_nbits, k, randgen, false);
  return check (C, randgen, niter, name);
}

/* */
bool check (RandGen &randgen, SecLevel seclevel, size_t k, size_t niter)
{
  std::stringstream desc;
  desc << " " << seclevel << " ";
  CL_HSM2k C (seclevel, k, randgen, false);
  return check (C, randgen, niter, desc.str());
}

/******************************************************************************/
int
main (int argc, char *argv[])
{
  bool success = true;

  RandGen randgen;
  randseed_from_argv (randgen, argc, argv);

  std::cout << "      |      | #bits  |    variants     |" << std::endl
            << " sec. |  k   | DeltaK | large | compact |" << std::endl;

  size_t k = 64;

  success &= check (randgen, 150, k, 500, "small");
  success &= check (randgen, 150, 128, 500, "small"); /* large message */
  success &= check (randgen, SecLevel::_112, k, 25);
  success &= check (randgen, SecLevel::_112, 1100, 25); /* large message */
  success &= check (randgen, SecLevel::_128, k, 20);
  success &= check (randgen, SecLevel::_192, k, 5);
  //success &= check (randgen, SecLevel::_256, k, 1);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
