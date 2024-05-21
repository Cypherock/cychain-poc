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

string generate_desc (const CL_HSMqk &C, const string &name)
{
  std::stringstream desc;
  desc << name << " | " << std::setw(5) << C.q().nbits()
               << " | " << std::setw(2) << C.k()
               << " | " << std::setw(6) << C.DeltaK().nbits()
               << " | " << (C.large_message_variant() ? " true" : "false")
               << " | " << (C.compact_variant() ? "  true " : " false ")
               << " |";
  return desc.str();
}

bool CL_HSMqk_check_setup (const CL_HSMqk &C, const string & pre)
{
  bool ret = true;
  Mpz tmp, one(1UL);

  ret &= C.q().is_prime ();
  ret &= C.p().is_prime ();

  ret &= mpz_fdiv_ui (static_cast<mpz_srcptr> (C.DeltaK()), 4) == 1;

  ret &= C.q().kronecker (C.p()) == -1;

  Mpz::mul (tmp, C.p(), C.q());
  tmp.neg();
  ret &= (tmp == C.DeltaK());

  for (size_t i = 0; i < C.k(); i++)
  {
    Mpz::mul (tmp, tmp, C.q());
    Mpz::mul (tmp, tmp, C.q());
  }
  ret &= tmp == C.Delta();

  /* Check that gen has the correct discriminant */
  ret &= C.h().discriminant() == C.Cl_G().discriminant();

  /* Check that h and f are squares */
  ret &= C.genus (C.h()) == typename CL_HSMqk::Genus ({ 1, 1 });
  ret &= C.genus (C.power_of_f(one)) == typename CL_HSMqk::Genus ({ 1, 1 });

  // TODO check precomputation values

  Test::result_line (pre + " setup", ret);
  return ret;
}

/* */
bool
test_CL_HSMqk_ZKAoK (const CL_HSMqk &C, OpenSSL::HashAlgo &H, RandGen &randgen,
                     size_t niter, const string & pre)
{
  bool ret = true;

  CL_HSMqk::SecretKey sk = C.keygen (randgen);
  CL_HSMqk::PublicKey pk = C.keygen (sk);

  for (size_t i = 0; i < niter; i++)
  {
    CL_HSMqk::ClearText a (C, randgen);
    Mpz r (randgen.random_mpz (C.encrypt_randomness_bound()));
    CL_HSMqk::CipherText c (C.encrypt (pk, a, r));
    CL_HSMqk_ZKAoKProof proof (C, H, pk, c, a, r, randgen);
    ret &= proof.verify (C, H, pk, c);
  }

  Test::result_line (pre + " ZKAoK", ret);
  return ret;
}

/* */
bool CL_HSMqk_check_all (const CL_HSMqk &C, OpenSSL::HashAlgo &H,
                         RandGen &randgen, size_t niter, const string &name)
{
  bool success = true;
  std::cout << std::endl;
  success &= CL_HSMqk_check_setup (C, name);
  success &= Test::test_encryption (C, randgen, niter, name);
  success &= Test::test_ciphertext_ops (C, randgen, niter, name);
  success &= test_CL_HSMqk_ZKAoK (C, H, randgen, niter, name);
  return success;
}

/* */
bool check (const CL_HSMqk &C, OpenSSL::HashAlgo &H, RandGen &randgen,
            size_t niter, const string &name)
{
  bool success = true;
  size_t DeltaK_nbits = C.DeltaK().nbits();

  success &= CL_HSMqk_check_all (C, H, randgen, niter, generate_desc (C, name));

  CL_HSMqk Ccompact (C, true);
  success &= CL_HSMqk_check_all (Ccompact, H, randgen, niter,
                                                generate_desc (Ccompact, name));

  CL_HSMqk Clarge ((DeltaK_nbits+C.k()-1)/C.k(), C.k(), DeltaK_nbits, randgen);
  success &= CL_HSMqk_check_all (Clarge, H, randgen, niter,
                                                  generate_desc (Clarge, name));

  CL_HSMqk Clargecompact (Clarge, true);
  success &= CL_HSMqk_check_all (Clargecompact, H, randgen, niter,
                                          generate_desc (Clargecompact, name));

  return success;
}

/* */
bool check (RandGen &randgen, size_t q_nbits, size_t k, size_t Delta_OK_nbits,
            size_t niter, const string &name)
{
  CL_HSMqk C (q_nbits, k, Delta_OK_nbits, randgen, false);
  SecLevel seclevel = SecLevel::_112; /* for hashalgo, use smallest seclevel */
  OpenSSL::HashAlgo H (seclevel);
  return check (C, H, randgen, niter, name);
}

/* */
bool check (RandGen &randgen, SecLevel seclevel, size_t k, size_t niter)
{
  std::string desc = std::string (" ") + std::to_string (seclevel) + " ";
  /* With a random q as big as the security level */
  CL_HSMqk C (seclevel, k, seclevel, randgen, false);
  OpenSSL::HashAlgo H (seclevel);
  return check (C, H, randgen, niter, desc);
}

/******************************************************************************/
int
main (int argc, char *argv[])
{
  bool success = true;

  RandGen randgen;
  randseed_from_argv (randgen, argc, argv);

  std::cout << "      | #bits |    | #bits  |    variants     |"
            << std::endl
            << " sec. |   q   | k  | DeltaK | large | compact |"
            << std::endl;

  success &= check (randgen, 50, 1, 150, 500, "small");
  success &= check (randgen, 5, 15, 150, 500, "small");
  success &= check (randgen, SecLevel::_112, 1, 25);
  success &= check (randgen, SecLevel::_112, 2, 25);
  success &= check (randgen, SecLevel::_128, 1, 20);
  success &= check (randgen, SecLevel::_192, 1, 5);
  //success &= check (randgen, SecLevel::_256, 1, 1);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
