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
#ifndef BICYCL_SECLEVEL_INL
#define BICYCL_SECLEVEL_INL

/* */
inline
const std::vector<SecLevel> SecLevel::All ()
{
  return { _112, _128, _192, _256 };
}

/* */
inline
SecLevel::SecLevel (unsigned int s)
{
  switch(s)
  {
    case 112 : value_ = _112; break;
    case 128 : value_ = _128; break;
    case 192 : value_ = _192; break;
    case 256 : value_ = _256; break;
    default  : throw InvalidSecLevelException() ; break;
  }
}

/* */
inline
SecLevel::SecLevel (const std::string &s)
{
  if (s == "112")       value_ = _112;
  else if (s == "128")  value_ = _128;
  else if (s == "192")  value_ = _192;
  else if (s == "256")  value_ = _256;
  else                  throw InvalidSecLevelException();
}

/* */
inline
size_t SecLevel::RSA_modulus_bitsize () const
{
  if (value_ == _112)        return 2048;
  else if (value_ == _128)   return 3072;
  else if (value_ == _192)   return 7680;
  else if (value_ == _256)   return 15360;
  else                       throw InvalidSecLevelException();
}

/* */
inline
size_t SecLevel::discriminant_bitsize () const
{
  if (value_ == _112)        return 1348;
  else if (value_ == _128)   return 1827;
  else if (value_ == _192)   return 3598;
  else if (value_ == _256)   return 5971;
  else                       throw InvalidSecLevelException();
}

/* */
inline
unsigned int SecLevel::soundness () const
{
  return static_cast<unsigned int>(value_);
}

/* */
inline
int SecLevel::elliptic_curve_openssl_nid () const
{
  if (value_ == _112)        return OpenSSL::ECGroup::P224;
  else if (value_ == _128)   return OpenSSL::ECGroup::P256;
  else if (value_ == _192)   return OpenSSL::ECGroup::P384;
  else if (value_ == _256)   return OpenSSL::ECGroup::P521;
  else                       throw InvalidSecLevelException();
}

/* */
inline
int SecLevel::sha3_openssl_nid () const
{
  if (value_ == _112)        return OpenSSL::HashAlgo::SHA3_224;
  else if (value_ == _128)   return OpenSSL::HashAlgo::SHA3_256;
  else if (value_ == _192)   return OpenSSL::HashAlgo::SHA3_384;
  else if (value_ == _256)   return OpenSSL::HashAlgo::SHA3_512;
  else                       throw InvalidSecLevelException();
}

/* */
inline
std::ostream & operator<< (std::ostream &o, SecLevel seclevel)
{
  return o << static_cast<unsigned int>(seclevel.value_);
}

/* */
inline
std::string to_string (SecLevel seclevel)
{
  return std::to_string (static_cast<unsigned int>(seclevel.value_));
}

#endif /* BICYCL_SECLEVEL_INL */
