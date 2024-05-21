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
#include <numeric>
#include <string>
#include <sstream>

#include "bicycl.hpp"
#include "internals.hpp"

using std::string;

using namespace BICYCL;
using BICYCL::Bench::ms;
using BICYCL::Bench::us;

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
void threshold_ECDSA_benchs_sign (const thresholdECDSA &C,
                                 RandGen &randgen, const string &pre)
{
  const OpenSSL::ECGroup & E = C.ec_group();

  for (unsigned int n = 2; n < 10; n += 2)
  {
    for (unsigned int t = 1; t < n; t++)
    {
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

      const size_t niter = 100;
      unsigned int p = S[0];

      std::stringstream desc;
      desc << pre << "| " << n << " | " << t;

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


      /* benchs keygen */
      auto keygen = [ &C, &p, &n, &t, &data1, &data2, &randgen, &CoQ, &Q_vec,
                                              &CoQSec, &V, &Sigma, &ZK, &PK ] ()
      {
        Keygen1 (C, n, t, p);
        Keygen2 (C, data1[p], randgen, CoQ, Q_vec, CoQSec, V, Sigma[p]);
        SecretKey (C, p, data1[p], data2[p], V, ZK, PK);
      };
      Bench::one_function<ms, ms> (keygen, niter, "keygen", desc.str());


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

      /* benchs sign */
      auto sign = [ &C, &p, &S, &sk, &randgen, &signdata1, &signdata2,
                    &signdata3, &signdata4, &signdata5, &signdata6, &signdata7,
                    &Co_map, &C1, &ZK1, &c_kg_map, &c_kw_map, &B_map,
                    &delta_map, &m, &Gamma_map, &CoS_map, &zk_map, Co2_map,
                    &V_map, &A_map, &C2S_map, &aok_map, &Co3_map, &U_map,
                    &T_map, &C3S_map, &s_map ] ()
      {
        Sign1 (C, p, S, sk[p], randgen);
        Sign2 (C, signdata1.at(p), sk[p], Co_map, C1, ZK1, randgen);
        Sign3 (C, signdata1.at(p), signdata2.at(p), sk[p], c_kg_map.at(p),
                                                  c_kw_map.at(p), B_map.at(p));
        Sign4 (C, signdata1.at(p), delta_map);
        Sign5 (C, signdata1.at(p), signdata2.at(p), signdata3.at(p),
                  signdata4.at(p), m, Gamma_map, CoS_map, zk_map);
        Sign6 (C, signdata5.at(p), Co2_map);
        Sign7 (C, signdata1.at(p), signdata5.at(p), signdata6.at(p),
                  sk.at(p), V_map, A_map, C2S_map, aok_map);
        Sign8 (C, signdata1.at(p), signdata7.at(p), Co3_map, U_map, T_map,
                                                                      C3S_map);
        Signature (C, signdata1.at(p), signdata5.at(p), sk[p], s_map);
      };
      Bench::one_function<ms, ms> (sign, niter, "sign", desc.str());
    }
  }
}

/* */
void threshold_ECDSA_benchs_all (const thresholdECDSA &C,
                                 RandGen &randgen, const string &pre)
{
  threshold_ECDSA_benchs_sign (C, randgen, pre);
}

/* */
int
main (int argc, char *argv[])
{
  RandGen randgen;
  randseed_from_argv (randgen, argc, argv);

  std::cout << "security level | n | t |     operation   |   time   "
            << "| time per op. " << std::endl;

  for (const SecLevel seclevel: SecLevel::All())
  {
    /* With a random q twice as big as the security level */
    std::stringstream desc;
    OpenSSL::HashAlgo H (seclevel);
    thresholdECDSA C (seclevel, randgen);
    desc << "  " << seclevel << " bits     ";
    threshold_ECDSA_benchs_all (C, randgen, desc.str());
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
