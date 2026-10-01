/*
Copyright (C) 2017-2025 The Kira Developers (see AUTHORS file)

This file is part of Kira.

Kira is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Kira is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Kira.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "kira/integral.h"
#include "kira/kira.h"
#include "kira/tools.h"
#include "pyred/integrals.h"

using namespace pyred;

using namespace std;
using namespace GiNaC;

istream& operator>>(istream& in, SEEDIntegral& per) {
  in >> per.id;
  for (unsigned i = 0; i < per.l_Indices; i++) {
    in >> per.indices[i];
  }
  in >> per.characteristics[SECTOR];
  in >> per.characteristics[TOPOLOGY];
  in >> per.characteristics[DENCOUNT];
  in >> per.characteristics[DOTS];
  in >> per.characteristics[NUM];
  in >> per.flag2;

  return in;
}

ostream& operator<<(ostream& out, const SEEDIntegral& per) {
  out << per.id << ",";
  for (unsigned i = 0; i < per.l_Indices; i++)
    out << per.indices[i] << ",";
  out << per.characteristics[SECTOR] << ",";
  out << per.characteristics[TOPOLOGY] << ",";
  out << per.characteristics[DENCOUNT] << ",";
  out << per.characteristics[DOTS] << ",";
  out << per.characteristics[NUM] << ",";
  out << per.flag2;

  return out;
}

IBPIntegral::IBPIntegral(pyred::Weight weights, unsigned l_Indices_) {
  length = 1;
  l_Indices = l_Indices_;
  auto iglback = pyred::Integral(weights);

  auto property = iglback.properties(weights);
  characteristics[SECTOR] = property.sector;
  characteristics[TOPOLOGY] = property.topology;
  characteristics[DENCOUNT] = property.lines;
  characteristics[DOTS] = property.dots;
  characteristics[NUM] = property.sps;

  for (unsigned i = 0; i < l_Indices; i++)
    indices[i] = iglback.m_powers[i];
}

void IBPIntegral::copy(const SEEDIntegral& integral, unsigned l_Indices_) {
  l_Indices = l_Indices_;
  id = integral.id;
  characteristics[TOPOLOGY] = integral.characteristics[TOPOLOGY];
  characteristics[SECTOR] = integral.characteristics[SECTOR];
  //   flag2 = integral.flag2;
}

void TESTIntegral::copy(const SEEDIntegral& integral, unsigned l_Indices_) {
  l_Indices = l_Indices_;
  id = integral.id;
  characteristics[TOPOLOGY] = integral.characteristics[TOPOLOGY];
  characteristics[SECTOR] = integral.characteristics[SECTOR];
  flag2 = integral.flag2;
}

int IBPIntegral::compare_indices(IBPIntegral& B) {
  unsigned count = 0;

  if (characteristics[TOPOLOGY] != B.characteristics[TOPOLOGY])
    return static_cast<int>(0);

  for (unsigned k = 0; k < l_Indices; k++) {
    if (indices[k] == B.indices[k]) {
      count++;
    }
  }

  return (static_cast<int>(count == l_Indices));
}

void TESTIntegral::generate_characteristics(int topology) {
  characteristics[SECTOR] = 0;
  characteristics[DENCOUNT] = 0;
  characteristics[DOTS] = 0;
  characteristics[NUM] = 0;
  characteristics[TOPOLOGY] = topology;

  for (unsigned kk = 0; kk < l_Indices; kk++) {
    int num = indices[kk];

    if (0 == num) continue;

    if (num > 0) {
      characteristics[SECTOR] += (1 << kk);
      characteristics[DENCOUNT]++;
      characteristics[DOTS] += num;
    }
    else
      characteristics[NUM] -= num;
  }
}

void IBPIntegral::generate_characteristics(int topology) {
  characteristics[SECTOR] = 0;
  characteristics[DENCOUNT] = 0;
  characteristics[DOTS] = 0;
  characteristics[NUM] = 0;
  characteristics[TOPOLOGY] = topology;

  for (unsigned kk = 0; kk < l_Indices; kk++) {
    int num = indices[kk];

    if (0 == num) continue;

    if (num > 0) {
      characteristics[SECTOR] += (1 << kk);
      characteristics[DENCOUNT]++;
      characteristics[DOTS] += num;
    }
    else
      characteristics[NUM] -= num;
  }
}

void SEEDIntegral::generate_characteristics(int topology) {
  characteristics[SECTOR] = 0;
  characteristics[DENCOUNT] = 0;
  characteristics[DOTS] = 0;
  characteristics[NUM] = 0;
  characteristics[TOPOLOGY] = topology;

  for (unsigned kk = 0; kk < l_Indices; kk++) {
    int num = indices[kk];

    if (0 == num) continue;

    if (num > 0) {
      characteristics[SECTOR] += (1 << kk);
      characteristics[DENCOUNT]++;
      characteristics[DOTS] += num;
    }
    else
      characteristics[NUM] -= num;
  }
}

istream& operator>>(istream& in, IBPIntegral& per) {
  in >> per.id;
  for (unsigned i = 0; i < per.l_Indices; i++) {
    in >> per.indices[i];
  }
  in >> per.characteristics[SECTOR];
  in >> per.characteristics[TOPOLOGY];
  in >> per.characteristics[DENCOUNT];
  in >> per.characteristics[DOTS];
  in >> per.characteristics[NUM];
  in >> per.flag2;

  return in;
}

ostream& operator<<(ostream& out, const IBPIntegral& per) {
  out << per.id << ",";
  for (unsigned i = 0; i < per.l_Indices; i++)
    out << per.indices[i] << ",";
  out << per.characteristics[SECTOR] << ",";
  out << per.characteristics[TOPOLOGY] << ",";
  out << per.characteristics[DENCOUNT] << ",";
  out << per.characteristics[DOTS] << ",";
  out << per.characteristics[NUM] << ",";
  out << per.flag2;

  return out;
}

void BaseIntegral::copy(const TESTIntegral& integral) {
  id = integral.id;
  length = integral.length;
  flag2 = integral.flag2;
  characteristics[SECTOR] = integral.characteristics[SECTOR];
  characteristics[TOPOLOGY] = integral.characteristics[TOPOLOGY];
  coefficientString = integral.coefficientString;
}

void BaseIntegral::copy(const IBPIntegral& integral) {
  id = integral.id;
  length = integral.length;
  flag2 = integral.flag2;
  characteristics[SECTOR] = integral.characteristics[SECTOR];
  characteristics[TOPOLOGY] = integral.characteristics[TOPOLOGY];
  //   string
  //   tokenString(integral.coefficientString.data(),integral.coefficientString.length());
  //   for( unsigned ii = 0; ii < integral.coefficientString2.size(); ii++)
  coefficientString = integral.coefficientString;
}

// void BaseIntegral::integral_buffer(char buffer[]) {
//
//   snprintf(buffer, sizeof(buffer), "%llu*(%s)\n", id, coefficientString.c_str());
//
// };

istream& operator>>(istream& in, BaseIntegral& per) {
  in >> per.length;
  in >> per.coefficientString;
  in >> per.id;
  in >> per.characteristics[SECTOR];
  in >> per.characteristics[TOPOLOGY];
  in >> per.flag2;
  return in;
}

ostream& operator<<(ostream& out, const BaseIntegral& per) {
  out << per.length << " ";
  out << per.coefficientString << " ";
  out << per.id << " ";
  out << per.characteristics[SECTOR] << " ";
  out << per.characteristics[TOPOLOGY] << " ";
  out << per.flag2 << endl;

  return out;
}

bool sort_rules4(BaseEquation l[], BaseEquation r[]) {
  if (l[0].equation[0].id < r[0].equation[0].id)
    return true;
  else if (l[0].equation[0].id > r[0].equation[0].id)
    return false;

  unsigned length = 0;
  if (l[0].l_Equation < r[0].l_Equation)
    length = l[0].l_Equation;
  else
    length = r[0].l_Equation;

  if (length > 1) {
    if (l[0].equation[1].id < r[0].equation[1].id)
      return true;
    else if (l[0].equation[1].id > r[0].equation[1].id)
      return false;
  }

  if (l[0].l_Equation < r[0].l_Equation)
    return true;
  else if (l[0].l_Equation > r[0].l_Equation)
    return false;

  for (unsigned it = 0; it < l[0].l_Equation; it++) {
    if (l[0].equation[it].id < r[0].equation[it].id)
      return true;
    else if (l[0].equation[it].id > r[0].equation[it].id)
      return false;
  }

  return false;
}

BaseEquation::BaseEquation() {}

BaseEquation::~BaseEquation() { delete[] equation; }

void BaseEquation::delete_IBP() { delete[] equationIBP; }

BaseEquation::BaseEquation(unsigned length_Equation_, unsigned length_Indices_)
    : l_Equation(length_Equation_), l_Indices(length_Indices_) {
  equation = new BaseIntegral[l_Equation];
}

BaseEquation::BaseEquation(const BaseEquation& ibp, pyred::Weight weights,
                           INTEGMAP& integralMap, const int topology,
                           Fermat*& fermat, int id, MASTERSMAP& mastersMap) {
  //  plant seeds
  l_Equation = ibp.l_Equation;
  l_Indices = ibp.l_Indices;

  VECTORequationTEST.resize(l_Equation);

  for (unsigned i = 0; i < l_Equation; i++) {
    VECTORequationTEST[i].l_Indices = l_Indices;
    VECTORequationTEST[i].characteristics[TOPOLOGY] = topology;
  }

  auto iglback = pyred::Integral(weights);

  //   Do not swap seed and VECTORequationTEST
  for (unsigned j = 0; j < l_Equation; j++) {
    stringstream options;
    char token[2048];

    vector<int> pows;
    pows.reserve(l_Indices);

    for (unsigned i = 0; i < l_Indices; i++) {
      VECTORequationTEST[j].indices[i] = ibp.equationIBP[j].indices[i];
      VECTORequationTEST[j].indices[i] += iglback.m_powers[i];

      pows.push_back(VECTORequationTEST[j].indices[i]);

      if (i > 0) {
        snprintf(token, sizeof(token), "+");
        options << token;
      }
      snprintf(token, sizeof(token), "(%s)*(%i)", ibp.equationIBP[j].coefficientIBP[i].data(),
              iglback.m_powers[i]);
      options << token;
    }

    snprintf(token, sizeof(token), "+(%s)*d",
              ibp.equationIBP[j].coefficientIBP[l_Indices].data());
    options << token;

    fermat[id].fermat_collect(const_cast<char*>(options.str().c_str()));

    fermat[id].fermat_calc();

    VECTORequationTEST[j].coefficientString.assign(fermat[id].g_baseout);

    auto igl = pyred::Integral(topology, std::move(pows));
    auto ID = igl.to_weight();

    VECTORequationTEST[j].id = ID;

    auto mapContent = integralMap.find(ID);

    if (mapContent != integralMap.end()) {
      VECTORequationTEST[j].characteristics[TOPOLOGY] =
          get<TOPOLOGY>(mapContent->second);
      VECTORequationTEST[j].characteristics[SECTOR] =
          get<SECTOR>(mapContent->second);
      VECTORequationTEST[j].flag2 = get<2>(mapContent->second);
    }
    else {
      auto property = igl.properties(ID);
      VECTORequationTEST[j].characteristics[SECTOR] = property.sector;
      VECTORequationTEST[j].characteristics[TOPOLOGY] = property.topology;
      VECTORequationTEST[j].flag2 = -1;
    }
  }

  eliminate_zeros(mastersMap);

  equation = new BaseIntegral[l_Equation];

  for (unsigned itt = 0; itt < l_Equation; itt++) {
    equation[itt].copy(VECTORequationTEST[itt]);
  }
}

void BaseEquation::collect_integrals(IBPVG& ibp) {
  for (unsigned it = 0; it < l_Equation; it++) {
    for (unsigned itt = it + 1; itt < l_Equation; itt++) {
      if ((*ibp[it]).id == (*ibp[itt]).id) {
        (*ibp[it]).coefficient += (*ibp[itt]).coefficient;
        (*ibp[it]).coefficient = (*ibp[it]).coefficient.expand();
        swap(*ibp[itt], *ibp[l_Equation - 1]);
        l_Equation--;
        itt--;
      }
    }

    (*ibp[it]).coefficientString = something_string((*ibp[it]).coefficient);

    if ((*ibp[it]).coefficientString == "0") {
      swap(*ibp[it], *ibp[l_Equation - 1]);
      l_Equation--;
      it--;
    }
  }

  for (unsigned it = 0; it < l_Equation; it++)
    (*ibp[it]).length = l_Equation;
}

void BaseEquation::eliminate_zeros(IBPVG& ibp, MASTERSMAP& mastersMap) {
  for (unsigned itt = 0; itt < l_Equation; itt++) {
    if ((*ibp[itt]).id != pyred::Weight(0u)) {
      continue;
    }
    swap(*ibp[itt], *ibp[l_Equation - 1]);
    l_Equation--;
    itt--;
  }

  for (unsigned itt = 0; itt < l_Equation; itt++) {
    (*ibp[itt]).length = l_Equation;

    auto mapContent = mastersMap.find((*ibp[itt]).id);

    if (mapContent != mastersMap.end()) {
      (*ibp[itt]).id = mapContent->second;
    }
  }
}

BaseEquation::BaseEquation(IBPVG& ibp, unsigned length_Indices_,
                           INTEGMAP& integralMap, MASTERSMAP& mastersMap) {
  l_Equation = static_cast<unsigned>(ibp.size());
  l_Indices = static_cast<unsigned>(length_Indices_);

  for (ItIBPVG i = ibp.begin(); i != ibp.end(); i++) {
    (*i)->l_Indices = l_Indices;
    (*i)->length = l_Equation;

    vector<int> pows;
    pows.reserve(l_Indices);

    for (unsigned j = 0; j < l_Indices; j++) {
      pows.push_back((*i)->indices[j]);
    }
    auto igl =
        pyred::Integral((*i)->characteristics[TOPOLOGY], std::move(pows));
    auto ID = igl.to_weight();
    (*i)->id = ID;

    auto mapContent = integralMap.find(ID);

    if (mapContent != integralMap.end()) {
      (*i)->flag2 = get<2>(mapContent->second);
    }
    else {
      (*i)->flag2 = -1;
    }
  }

  collect_integrals(ibp);

  eliminate_zeros(ibp, mastersMap);

  equation = new BaseIntegral[l_Equation];

  for (unsigned itt = 0; itt < l_Equation; itt++) {
    equation[itt].copy(*ibp[itt]);
  }
}

BaseEquation::BaseEquation(IBPVG& ibp, unsigned length_Indices_) {
  l_Equation = static_cast<unsigned>(ibp.size());
  l_Indices = length_Indices_;

  if (l_Equation != 0) {
    equationIBP = new IBPIntegral[l_Equation];
    equationIBP[0].l_Indices = l_Indices;
    equationIBP[0].length = l_Equation;

    int jj = 0;

    for (ItIBPVG i = ibp.begin(); i != ibp.end(); i++) {
      for (unsigned k = 0; k < l_Indices; k++) {
        equationIBP[jj].indices[k] = (*i)->indices[k];
      }

      equationIBP[jj].characteristics[TOPOLOGY] =
          (*i)->characteristics[TOPOLOGY];
      equationIBP[jj].l_Indices = l_Indices;
      equationIBP[jj].length = l_Equation;
      equationIBP[jj].coefficient = (*i)->coefficient;
      equationIBP[jj].coefficientString =
          something_string(equationIBP[jj].coefficient);

      jj++;
      delete *i;
    }
    equation = new BaseIntegral[1];

    collect_integrals();
    dismember_coef();
  }
}

void BaseEquation::write_file(ofstream& out, string& topology) {

  for (unsigned i = 0; i < l_Equation; i++) {

    out << topology <<"[";
    for (unsigned j = 0; j < l_Indices; j++) {
      out << equationIBP[i].indices[j];
      if (j < l_Indices - 1) {
        out << ",";
      }
    }
    out << "]*(" << equationIBP[i].coefficientString << ")" << endl;
  }
  out << endl;
}

void BaseEquation::collect_integrals() {
  for (unsigned it = 0; it < l_Equation; it++) {
    for (unsigned itt = it + 1; itt < l_Equation; itt++) {
      if (equationIBP[it].compare_indices(equationIBP[itt])) {
        equationIBP[it].coefficient += equationIBP[itt].coefficient;
        equationIBP[it].coefficient = equationIBP[it].coefficient.expand();
        swap(equationIBP[itt], equationIBP[l_Equation - 1]);
        l_Equation--;
        itt--;
      }
    }

    equationIBP[it].coefficientString =
        something_string(equationIBP[it].coefficient);

    if (equationIBP[it].coefficientString == "0") {
      swap(equationIBP[it], equationIBP[l_Equation - 1]);
      l_Equation--;
      it--;
    }
  }

  for (unsigned it = 0; it < l_Equation; it++)
    equationIBP[it].length = l_Equation;
}

void BaseEquation::eliminate_zeros(MASTERSMAP& mastersMap) {
  for (unsigned itt = 0; itt < l_Equation; itt++) {
    if (VECTORequationTEST[itt].coefficientString != "0") {
      if (VECTORequationTEST[itt].id != pyred::Weight(0u)) {
        continue;
      }
    }
    swap(VECTORequationTEST[itt], VECTORequationTEST[l_Equation - 1]);
    l_Equation--;
    itt--;
  }

  for (unsigned itt = 0; itt < l_Equation; itt++) {
    VECTORequationTEST[itt].length = l_Equation;

    auto mapContent = mastersMap.find(VECTORequationTEST[itt].id);

    if (mapContent != mastersMap.end()) {
      VECTORequationTEST[itt].id = mapContent->second;
    }
  }
}

int BaseEquation::sort() {
  if (l_Equation > 0) {
    for (unsigned i = 0; i < l_Equation - 1; i++) {
      for (unsigned itt = i; itt < l_Equation; itt++) {
        if ((equation[itt].id) > (equation[i].id)) {
          swap(equation[i], equation[itt]);
        }
      }
    }
    for (unsigned itt = 0; itt < l_Equation; itt++)
      equation[itt].length = l_Equation;
  }
  return l_Equation;
}

void BaseEquation::dismember_coef() {
  possymbol* indices = new possymbol[l_Indices];

  generate_symbols(indices, "a", l_Indices);

  for (unsigned i = 0; i < l_Equation; i++) {
    for (unsigned k = 0; k < l_Indices; k++) {
      ex token = diff(equationIBP[i].coefficient, indices[k]);
      equationIBP[i].coefficientIBP[k] = something_string(token);
    }
    ex token = diff(equationIBP[i].coefficient, get_symbol("d"));
    equationIBP[i].coefficientIBP[l_Indices] = something_string(token);
  }
  delete[] indices;
}
