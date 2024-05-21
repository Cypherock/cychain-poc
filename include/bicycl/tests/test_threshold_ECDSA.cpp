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
#include <algorithm>
#include <numeric>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

#include "bicycl.hpp"
#include "internals.hpp"

using std::string;

using namespace BICYCL;

/* */
class LogicalError : public std::runtime_error
{
  public:
    using runtime_error::runtime_error;
};

/* */
std::ostream & operator<< (std::ostream &o, const std::vector<unsigned int> &v)
{
  o << '[';
  for (unsigned int i: v)
    o << ' ' << i << ',';
  return o << ']';
}

/* */
bool
test_commitments (const thresholdECDSA &C, size_t niter, const string &pre)
{
  bool ret = true;

  thresholdECDSA::Commitment c;
  thresholdECDSA::CommitmentSecret r;

  for (size_t i = 0; i < niter; i++)
  {
    OpenSSL::BN k (C.ec_group().random_mod_order());
    OpenSSL::ECPoint Q (C.ec_group(), k);
    tie(c, r) = C.commit (Q);
    ret &= C.open (c, Q, r);
  }

  string line (pre + " commitments");
  Test::result_line (line, ret);
  return ret;
}

using Keygen1 = thresholdECDSA::KeygenPart1;
using Keygen2 = thresholdECDSA::KeygenPart2;
using SecretKey = thresholdECDSA::SecretKey;
using Sign1 = thresholdECDSA::SignPart1;
using Sign2 = thresholdECDSA::SignPart2;
using Sign3 = thresholdECDSA::SignPart3;
using Sign4 = thresholdECDSA::SignPart4;
using Sign5 = thresholdECDSA::SignPart5;
using Sign6 = thresholdECDSA::SignPart6;
using Sign7 = thresholdECDSA::SignPart7;
using Sign8 = thresholdECDSA::SignPart8;
using Signature = thresholdECDSA::Signature;
template <class T>
using PMap = thresholdECDSA::ParticipantsMap<T>;

/* */
bool
test_sign (const thresholdECDSA &C, size_t niter, RandGen &randgen,
           const string &pre)
{
  bool ret = true;

  const OpenSSL::ECGroup & E = C.ec_group();

  for (size_t iter = 0; iter < niter; iter++)
  {
    unsigned int n = 2 + randgen.random_ui_2exp (3); /* in [2,9] */
    unsigned int t = 1 + randgen.random_ui (n-1); /* in [1, n-1] */
    std::vector<unsigned int> S(n);
    std::iota(S.begin(), S.end(), 0); /* S contains [0..n-1] */
    while (S.size() > t+1)
      S.erase (S.begin() + randgen.random_ui (S.size()));

    std::vector<Keygen1> data1;     data1.reserve (n);
    std::vector<Keygen2> data2;     data2.reserve (n);
    std::vector<SecretKey> sk;      sk.reserve (n);

    PMap<Sign1> signdata1;          signdata1.reserve (t);
    PMap<Sign2> signdata2;          signdata2.reserve (t);
    PMap<Sign3> signdata3;          signdata3.reserve (t);
    PMap<Sign4> signdata4;          signdata4.reserve (t);
    PMap<Sign5> signdata5;          signdata5.reserve (t);
    PMap<Sign6> signdata6;          signdata6.reserve (t);
    PMap<Sign7> signdata7;          signdata7.reserve (t);
    PMap<Sign8> signdata8;          signdata8.reserve (t);
    PMap<Signature> signature;      signature.reserve (t);

    const thresholdECDSA::Message &m = C.random_message();
    try
    {
      /***** Keygen ***********************************************************/

      /* keygen part 1 */
      for (unsigned int i = 0; i < n; i++)
      {
        data1.push_back (Keygen1 (C, n, t, i));
      }

      /* fake the communication */
      std::vector<thresholdECDSA::Commitment> CoQ;
      std::vector<OpenSSL::ECPoint> Q_vec;
      std::vector<thresholdECDSA::CommitmentSecret> CoQSec;
      std::vector<std::vector<OpenSSL::ECPoint>> V;
      std::vector<std::vector<OpenSSL::BN>> Sigma (n);
      for (unsigned int i = 0; i < n; i++)
      {
        CoQ.push_back (data1[i].commitment());
        Q_vec.push_back (OpenSSL::ECPoint (E, data1[i].Q_part()));
        CoQSec.push_back (data1[i].commitment_secret());
        V.push_back (std::vector<OpenSSL::ECPoint>());

        for (unsigned k = 0; k < t; k++)
          V[i].push_back (OpenSSL::ECPoint (E, data1[i].V(k)));

        for (unsigned j = 0; j < n; j++)
          Sigma[j].push_back (data1[i].sigma (j));
      }

      /* keygen part 2 */
      for (unsigned int i = 0; i < n; i++)
      {
        data2.push_back (Keygen2 (C, data1[i], randgen, CoQ, Q_vec, CoQSec, V,
                                                                    Sigma[i]));
      }

      /* fake the communication */
      std::vector<CL_HSMqk::PublicKey> PK;
      std::vector<ECNIZKProof> ZK;
      for (unsigned int i = 0; i < n; i++)
      {
        PK.push_back (data2[i].CL_public_key());
        ZK.push_back (ECNIZKProof (E, data2[i].zk_proof()));
      }

      /* keygen part 3 */
      for (unsigned int i = 0; i < n; i++)
      {
        sk.push_back (SecretKey (C, i, data1[i], data2[i], V, ZK, PK));
      }

      /***** Check *****/
      /* Check that all public key are the same and equal to sum of Qi */
      OpenSSL::ECPoint Q (E, data1[0].Q_part());
      for (unsigned int i = 1; i < n; i++)
      {
        E.ec_add (Q, Q, data1[i].Q_part());
      }
      for (unsigned int i = 1; i < n; i++)
      {
        if (!E.ec_point_eq (sk[i].public_key(), Q))
        {
          throw ("public key Q do not match");
        }
      }

      /* Check that Xj = xi P */
      OpenSSL::ECPoint X (E);
      for (unsigned int i = 0; i < n; i++)
      {
        E.scal_mul_gen (X, sk[i].x_part()); /* xi P */
        for (unsigned int j = 0; j < n; j++)
        {
          if (i != j && !E.ec_point_eq (X, sk[j].X(i)))
          {
            throw LogicalError ("cannot verify Xj = xi P");
          }
        }
      }

      /* Check that u = sum{u_i} can be retrieved from secret xi's */
      OpenSSL::BN u (data1[0].u_part());
      for (unsigned int i = 1; i < n; i++)
      {
        E.add_mod_order (u, u, data1[i].u_part());
      }
      OpenSSL::BN ut(0UL), l;
      for (unsigned int s: S)
      {
        l = C.lagrange_at_zero (S, s);
        E.mul_mod_order (l, l, sk[s].x_part());
        E.add_mod_order (ut, ut, l);
      }
      if (u != ut)
      {
        throw LogicalError ("cannot reconstruct u from secret {xi}");
      }

      /* Check that Q = uP */
      OpenSSL::ECPoint U (E);
      E.scal_mul_gen (U, u);
      if (!E.ec_point_eq (Q, U))
      {
        throw LogicalError ("Q is not u P");
      }

      /***** Signing **********************************************************/

      /* signing part 1 */
      for (unsigned int i: S)
      {
        signdata1.emplace (i, Sign1 (C, i, S, sk[i], randgen));
      }

      /* fake the communication */
      PMap<thresholdECDSA::Commitment> Co_map;
      PMap<CL_HSMqk::CipherText> C1;
      PMap<CL_HSMqk_ZKAoKProof> ZK1;
      for (unsigned int i: S)
      {
        Co_map.insert ({i, signdata1.at(i).commitment()});
        C1.insert ({i, signdata1.at(i).ciphertext()});
        ZK1.insert ({i, signdata1.at(i).zk_encrypt_proof()});
      }

      /* signing part 2 */
      for (unsigned int i: S)
      {
        signdata2.emplace (i, Sign2 (C, signdata1.at(i), sk[i], Co_map, C1, ZK1,
                                        randgen));
      }

      /* fake the communication */
      PMap<PMap<CL_HSMqk::CipherText>> c_kg_map;
      PMap<PMap<CL_HSMqk::CipherText>> c_kw_map;
      PMap<PMap<OpenSSL::ECPoint>> B_map;
      for (unsigned int i: S)
      {
        for (unsigned int j: S)
        {
          if (i != j)
          {
            c_kg_map[i].insert ({j, signdata2.at(j).c_kg(i)});
            c_kw_map[i].insert ({j, signdata2.at(j).c_kw(i)});
            B_map[i].insert ({j, OpenSSL::ECPoint (E, signdata2.at(j).B(i))});
          }
        }
      }

      /* signing part 3 */
      for (unsigned int i: S)
      {
        signdata3.emplace (i, Sign3 (C, signdata1.at(i), signdata2.at(i), sk[i],
                                        c_kg_map.at(i), c_kw_map.at(i),
                                        B_map.at(i)));
      }

      /* fake the communication */
      PMap<OpenSSL::BN> delta_map;
      for (unsigned int i: S)
      {
        delta_map.insert ({i, signdata3.at(i).delta_part()});
      }

      /* signing part 4 */
      for (unsigned int i: S)
      {
        signdata4.emplace (i, Sign4 (C, signdata1.at(i), delta_map));
      }

      /* check */
      OpenSSL::BN k(0UL), gamma(0UL), delta;
      for (unsigned int i: S)
      {
        OpenSSL::BN ki (signdata1.at(i).k_part());
        E.add_mod_order (k, k, ki);
        E.add_mod_order (gamma, gamma, signdata1.at(i).gamma());
      }
      E.mul_mod_order (delta, k, gamma);
      for (unsigned int i: S)
      {
        if (delta != signdata4.at(i).delta())
        {
          throw LogicalError ("delta is not equal to k * gamma");
        }
      }

      /* fake the communication */
      PMap<ECNIZKProof> zk_map;
      PMap<thresholdECDSA::CommitmentSecret> CoS_map;
      PMap<OpenSSL::ECPoint> Gamma_map;
      for (unsigned int i: S)
      {
        zk_map.insert ({i, ECNIZKProof (E, signdata1.at(i).zk_gamma())});
        CoS_map.insert ({i, signdata1.at(i).commitment_secret()});
        Gamma_map.insert ({i, OpenSSL::ECPoint (E, signdata1.at(i).Gamma())});
      }

      /* signing part 5 */
      for (unsigned int i: S)
      {
        signdata5.emplace (i, Sign5 (C, signdata1.at(i), signdata2.at(i),
                                        signdata3.at(i), signdata4.at(i),
                                        m, Gamma_map, CoS_map, zk_map));
      }

      /* check */
      for (unsigned int i: S)
      {
        if (!E.ec_point_eq(signdata5.at(i).R(), signdata5.at(S[0]).R()))
        {
          throw LogicalError ("R does not match");
        }
        if (signdata5.at(i).r() != signdata5.at(S[0]).r())
        {
          throw LogicalError ("r does not match");
        }
        if (signdata5.at(i).z() != signdata5.at(S[0]).z())
        {
          throw LogicalError ("z does not match");
        }
      }

      /* fake the communication */
      PMap<thresholdECDSA::Commitment> Co2_map;
      for (unsigned int i: S)
      {
        Co2_map.insert ({i, signdata5.at(i).commitment()});
      }

      /* signing part 6 */
      for (unsigned int i: S)
      {
        signdata6.emplace (i, Sign6 (C, signdata5.at(i), Co2_map));
      }

      /* fake the communication */
      PMap<ECNIZKAoK> aok_map;
      PMap<thresholdECDSA::CommitmentSecret> C2S_map;
      PMap<OpenSSL::ECPoint> V_map;
      PMap<OpenSSL::ECPoint> A_map;
      for (unsigned int i: S)
      {
        aok_map.insert ({i, ECNIZKAoK (E, signdata6.at(i).aok())});
        C2S_map.insert ({i, signdata5.at(i).commitment_secret()});
        V_map.insert ({i, OpenSSL::ECPoint (E, signdata5.at(i).V_part())});
        A_map.insert ({i, OpenSSL::ECPoint (E, signdata5.at(i).A_part())});
      }

      /* signing part 7 */
      for (unsigned int i: S)
      {
        signdata7.emplace (i, Sign7 (C, signdata1.at(i), signdata5.at(i),
                                        signdata6.at(i), sk.at(i), V_map, A_map,
                                        C2S_map, aok_map));
      }

      /* fake the communication */
      PMap<thresholdECDSA::Commitment> Co3_map;
      PMap<thresholdECDSA::CommitmentSecret> C3S_map;
      PMap<OpenSSL::ECPoint> U_map;
      PMap<OpenSSL::ECPoint> T_map;
      for (unsigned int i: S)
      {
        Co3_map.insert ({i, signdata7.at(i).commitment()});
        C3S_map.insert ({i, signdata7.at(i).commitment_secret()});
        U_map.insert ({i, OpenSSL::ECPoint (E, signdata7.at(i).U_part())});
        T_map.insert ({i, OpenSSL::ECPoint (E, signdata7.at(i).T_part())});
      }

      /* signing part 8 */
      for (unsigned int i: S)
      {
        signdata8.emplace (i, Sign8 (C, signdata1.at(i), signdata7.at(i),
                                        Co3_map, U_map, T_map, C3S_map));
      }

      /* fake the communication */
      PMap<OpenSSL::BN> s_map;
      for (unsigned int i: S)
      {
        s_map.insert ({i, signdata5.at(i).s_part()});
      }

      /* signing part 9 */
      for (unsigned int i: S)
      {
        signature.emplace (i, Signature(C, signdata1.at(i), signdata5.at(i),
                                                            sk[i], s_map));
      }

      /***** Final check *****/
      for (unsigned int i: S)
      {
        if (signature.at(i) != signature.at(S[0]))
        {
          throw LogicalError ("signatures does not match");
        }

        if (!C.verify (signature.at(i), sk[i].public_key(), m))
        {
          throw LogicalError ("could not verify signature");
        }
      }
    }
    catch (LogicalError &e)
    {
      std::cerr << iter << ": LogicalError, " << e.what() << std::endl;
      ret = false;
    }
    catch (thresholdECDSA::ProtocolAbortError &e)
    {
      std::cerr << iter << ": ProtocolAbortError, " << e.what() << std::endl;
      ret = false;
    }
  }

  Test::result_line (pre, ret);
  return ret;
}

/******************************************************************************/
bool check (SecLevel seclevel, RandGen &randgen, size_t niter)
{
  bool success = true;

  std::stringstream desc;
  desc << "security " << seclevel << " bits";

  thresholdECDSA C (seclevel, randgen);
  desc << " threshold ECDSA";

  success &= test_commitments (C, niter, desc.str());
  success &= test_sign (C, niter, randgen, desc.str());

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

  success &= check (SecLevel::_112, randgen, 15);
  success &= check (SecLevel::_128, randgen, 10);
  success &= check (SecLevel::_192, randgen, 5);
  success &= check (SecLevel::_256, randgen, 1);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
