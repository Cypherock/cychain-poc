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

/*
 * This code illustrates how to use the BICYCL library to implement a program
 * responding to the old Chia VDF competition (for more details, see
 * https://medium.com/@chia.net/chia-vdf-competition-guide-5382e1f4bd39)
 *
 * The program reads, on the command line, a (negative) discriminant D and a
 * number n and output the coefficients a and b of the quadratic form
 * corresponding to (2, 1)^n in the class group Cl(D).
 */
#include <iostream>

#include "bicycl.hpp"

/* All BICYCL classes and functions are in the BICYCL namespace */
using namespace BICYCL;

int main (int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <disc> <niter>" << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    Mpz D (argv[1]);
    Mpz N (argv[2]);
    size_t niter = (size_t) N;

    ClassGroup Cl (D); /* build the class group of discriminant D */
    QFI f = Cl.primeform (Mpz (2UL)); /* build the prime form for p = 2 */

    Cl.nudupl (f, f, niter); /* perform niter nudupl */
    std::cout << f.a() << std::endl << f.b(); /* print the coefficients */

    return EXIT_SUCCESS;
  }
}
