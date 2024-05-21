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
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>

#include "bicycl.hpp"
#include "internals.hpp"

using namespace BICYCL;

using std::string;
using std::ifstream;
using std::istringstream;

/******************************************************************************/
bool test_qfi_nupow ()
{
  bool ret = true;

  #include "test_qfi_nupow.data.hpp"

  QFI r, f, fe, fd, fde, ref;
  Mpz n, nl, nh;

  for (const auto &data: Data)
  {
    bool b = true;
    std::tie(f, n, ref) = data;
    ClassGroup Cl (f.discriminant());

    /* Test first version of nupow */
    Cl.nupow (r, f, n);
    b = (r == ref);
    if (!b)
    {
      std::cout << "Failure, computed result is wrong for "
                << f << "^" << n << " (first version):" << std::endl
                << " expected " << ref << std::endl
                << "      got " << r << std::endl;
    }

    /* Precomputation to test other versions */
    size_t d = (n.nbits()+1)/2;
    size_t e = d/2 + 1;
    fde = f;
    for (size_t i = 0; i < d+e; i++)
    {
      if (i == e)
        fe = fde;
      if (i == d)
        fd = fde;
      Cl.nudupl (fde, fde);
    }

    /* Test second version of nupow */
    Mpz::mod2k (nl, n, d);
    Mpz::divby2k (nh, n, d);
    Cl.nupow (r, f, nl, fd, nh);
    b = (r == ref);
    if (!b)
    {
      std::cout << "Failure, computed result is wrong for "
                << f << "^" << n << " (second version):" << std::endl
                << " expected " << ref << std::endl
                << "      got " << r << std::endl;
    }

    /* Test third version of nupow */
    Cl.nupow (r, f, n, d, e, fe, fd, fde);
    b = (r == ref);
    if (!b)
    {
      std::cout << "Failure, computed result is wrong for "
                << f << "^" << n << " (third version):" << std::endl
                << " expected " << ref << std::endl
                << "      got " << r << std::endl;
    }

    ret &= b;
  }

  Test::result_line ("test_qfi_nupow", ret);
  return ret;
}

bool
test_class_number_bound ()
{
  bool ret = true;

  #include "test_qfi_class_number.data.hpp"

  Mpz B, h;
  ClassGroup Cl (Mpz("-3"));

  for (const auto &data: Data)
  {
    std::tie(Cl, h) = data;
    B = Cl.class_number_bound ();

    bool b = (B <= h);
    Mpz::mulby2 (B, B);
    b &= (h < B);

    if (!b)
    {
      Mpz::divby2 (B, B);
      std::cout << "Failure, bound is wrong for discriminant D="
                << Cl.discriminant() << std::endl
                << " expected B <= h < 2*B with h=" << h << std::endl
                << " got B=" << B << std::endl;
    }

    ret &= b;
  }

  Test::result_line ("test_class_number_bound", ret);
  return ret;
}

/* */
bool test_compressed_repr (RandGen &randgen, size_t niter)
{
  bool ret = true;

  for (size_t i = 0; i < niter; i++)
  {
    size_t disc_nbits = 100 + i*50;
    ClassGroup Cl (randgen.random_negative_discriminant (disc_nbits));
    QFI f (Cl.random (randgen));
    QFICompressedRepresentation fc (f.compressed_repr ());
    QFI ft (fc, Cl.discriminant());
    ret &= (f == ft);
  }

  Test::result_line ("test_compressed_repr", ret);
  return ret;
}

/******************************************************************************/
int
main (int argc, char *argv[])
{
  bool success = true;

  RandGen randgen;
  randseed_from_argv (randgen, argc, argv);

  success &= test_qfi_nupow ();
  success &= test_class_number_bound ();
  success &= test_compressed_repr (randgen, 100);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
