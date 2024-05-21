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
#ifndef BICYCL_BICYCL_HPP
#define BICYCL_BICYCL_HPP

#ifndef BICYCL_GMP_PRIMALITY_TESTS_ITERATION
#define BICYCL_GMP_PRIMALITY_TESTS_ITERATION 30
#endif

#include "bicycl/gmp_extras.hpp"
#include "bicycl/openssl_wrapper.hpp"
#include "bicycl/qfi.hpp"
#include "bicycl/seclevel.hpp"
#include "bicycl/ec.hpp"
#include "bicycl/CL_HSMqk.hpp"
#include "bicycl/CL_HSM2k.hpp"
#include "bicycl/Paillier.hpp"
#include "bicycl/Joye_Libert.hpp"
#include "bicycl/threshold_ECDSA.hpp"

namespace BICYCL
{
  #include "bicycl/seclevel.inl"
} /* namespace BICYCL */

#endif /* BICYCL_BICYCL_HPP */
