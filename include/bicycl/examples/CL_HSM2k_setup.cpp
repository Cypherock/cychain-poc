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
#include <iostream>
#include <chrono>
#include <stdexcept>

#include "bicycl.hpp"

using std::string;

void usage (const string &argv0)
{
  std::cout << "Usage: " << argv0
            << " -seclevel <integer> [-compact-variant] [-k <k>] [-seed <s>]"
            << std::endl << std::endl
            << "Parameters:" << std::endl
            << "  -seclevel <integer> target security level (";

  size_t i = 0;
  for (const BICYCL::SecLevel v: BICYCL::SecLevel::All())
  {
    if (i > 0 && i+1 < BICYCL::SecLevel::All().size())
      std::cout << ", ";
    else if (i+1 == BICYCL::SecLevel::All().size())
      std::cout << " or ";
    std::cout << v;
    i++;
  }
  std::cout << ")" << std::endl
            << "  -k <k>              value used to setup CL_HSM2k" << std::endl
            << "  -seed <s>           seed for random generator" << std::endl
            << "  -compact-variant    use the compact variant" << std::endl;
}

int
main (int argc, char *argv[])
{
  BICYCL::Mpz seed;
  size_t k = 0;
  char *seclevel_str = NULL;
  BICYCL::RandGen randgen;
  bool compact_variant = false; /* by default the compact variant is not used */
  string argv0 (argv[0]);

  auto T = std::chrono::system_clock::now();
  seed = static_cast<unsigned long>(T.time_since_epoch().count());

  /* First look for options */
  argc--; argv++; /* skip binary name */
  while (argc && argv[0][0] == '-')
  {
    string argument (argv[0]);
    if (argc >= 2 && argument == "-seclevel")
    {
      seclevel_str = argv[1];
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 2 && argument == "-k")
    {
      try
      {
        BICYCL::Mpz kmpz(argv[1]);
        k = static_cast<unsigned long> (kmpz);
      }
      catch (const std::runtime_error& e)
      {
        std::cerr << "Error, '" << argv[1] << "' could not be parsed: "
                  << e.what() << std::endl;
        usage (argv0);
        return EXIT_FAILURE;
      }
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 2 && argument == "-seed")
    {
      seed = string (argv[1]);
      argv += 2;
      argc -= 2;
    }
    else if (argc >= 1 && argument == "-compact-variant")
    {
      compact_variant = true;
      argv += 1;
      argc -= 1;
    }
    else
    {
      std::cerr << "Error, unknown option '" << argument << "'" << std::endl;
      usage (argv0);
      return EXIT_FAILURE;
    }
  }

  /* If unparsed options remaining => error */
  if (argc)
  {
    std::cerr << "Error, could not parsed option '" << argv[0] << "'" << std::endl;
    usage (argv0);
    return EXIT_FAILURE;
  }

  /* */
  if (seclevel_str == NULL)
  {
    std::cerr << "Error, missing -seclevel parameter" << std::endl;
    usage (argv0);
    return EXIT_FAILURE;
  }

  /* */
  if (k == 0)
  {
    std::cerr << "Error, missing -k parameter" << std::endl;
    usage (argv0);
    return EXIT_FAILURE;
  }

  /* */
  std::cout << "# Using seed = " << seed << std::endl;
  randgen.set_seed (seed);

  /* */
  BICYCL::SecLevel seclevel (seclevel_str);
  std::cout << "# security: " << seclevel << " bits" << std::endl;

  /* */
  std::cout << BICYCL::CL_HSM2k (seclevel, k, randgen, compact_variant);

  return EXIT_SUCCESS;
}
