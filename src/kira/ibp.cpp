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

#include "kira/kira.h"
#include "kira/tools.h"

using namespace std;
using namespace GiNaC;

static Loginfo& logger = Loginfo::instance();

void Kira::reduce_scal(IBPVG& ibp) {
  /* Search and replace reducible scalar products with irreducible ones:*/
  ex coef1, number;
  IBPIntegral* newIntegral;
  IBPVG ibp_reduced;
  possymbol v("v");
  for (ItIBPVG i = ibp.begin(), ends = ibp.end(); i != ends; ++i) {
    (*i)->coefficient =
        subs((*i)->coefficient.expand(), integralfamily.scal2Props[0],
             subs_options::algebraic);
    fs<lst>((*i)->coefficient, kinematic);
    (*i)->coefficient = (*i)->coefficient.expand();
    for (size_t k = 0; k < integralfamily.props.nops(); ++k) {
      if (0 != (coef1 = diff(subs((*i)->coefficient.expand(),
                                  integralfamily.propSymb[k] == v,
                                  subs_options::algebraic),
                             v))) {
        newIntegral = new IBPIntegral;
        for (size_t j = 0; j < integralfamily.props.nops(); j++) {
          newIntegral->indices[j] = (*i)->indices[j];
        }
        --(newIntegral->indices[k]);
        newIntegral->coefficient = coef1.expand();
        fs<lst>(newIntegral->coefficient, kinematic);
        newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
        ibp_reduced.push_back(newIntegral);
        (*i)->coefficient = (*i)->coefficient.subs(
            integralfamily.propSymb[k] == 0, subs_options::algebraic);
      }
    }
    if ((*i)->coefficient != 0) {
      (*i)->coefficient = (*i)->coefficient.expand();
      fs<lst>((*i)->coefficient, kinematic);
      (*i)->characteristics[TOPOLOGY] = integralfamily.topology;
      ibp_reduced.push_back((*i));
    }
    else
      delete *i;
  }
  ibp.clear();
  ibp = ibp_reduced;
}

bool sortSize(BaseEquation*& l, BaseEquation*& r) {
  if (l->l_Equation < r->l_Equation) return true;
  return false;
}

void Kira::create_IBP_helper(vector<possymbol>& var) {
  IBPVG ibp;
  IBPIntegral* newIntegral;
  ex coef1;
  /*These are the exponents of the propagators in the integrand.*/
  possymbol* indices = new possymbol[integralfamily.jule];
  generate_symbols(indices, "a", integralfamily.jule);
  possymbol* symbolicInd = new possymbol[integralfamily.jule];
  generate_symbols(symbolicInd, "b", integralfamily.jule);

  for (size_t iii = 0; iii < var.size(); iii++) {
    for (size_t ii = 0; ii < integralfamily.loopVar.size(); ii++) {
      ibp.clear();

      for (size_t i = 0; i < integralfamily.props.nops(); i++) {

        auto symbIt = find(integralfamily.symbolicIBP.begin(),
                           integralfamily.symbolicIBP.end(), i);

        if (symbIt != integralfamily.symbolicIBP.end())

          coef1 = (diff(integralfamily.props[i], integralfamily.loopVar[ii]) *
                   (-(indices[i] + symbolicInd[i])) * var[iii]);
        else

          coef1 = (diff(integralfamily.props[i], integralfamily.loopVar[ii]) *
                   (-(indices[i])) * var[iii]);
        fs<lst>(coef1, kinematic);
        if (coef1 != 0) {
          newIntegral = new IBPIntegral;
          for (size_t j = 0; j < integralfamily.props.nops(); j++) {
            newIntegral->indices[j] = 0;
          }
          newIntegral->coefficient = coef1;
          newIntegral->indices[i] = 1;
          newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
          ibp.push_back(newIntegral);
        }
      }
      if (var[iii] == integralfamily.loopVar[ii]) {
        newIntegral = new IBPIntegral;
        newIntegral->coefficient = dimension;
        for (size_t j = 0; j < integralfamily.props.nops(); j++) {
          newIntegral->indices[j] = 0;
        }
        newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
        ibp.push_back(newIntegral);
      }
      reduce_scal(ibp);

      BaseEquation* equationIBP;
      equationIBP = new BaseEquation(ibp, integralfamily.jule);

      integralfamily.identitiesIBP.push_back(equationIBP);
    }
  }

  delete[] indices;
  delete [] symbolicInd;
}

void Kira::create_LEE_vectors(vector<possymbol>& var) {
  IBPVG ibp;
  IBPIntegral* newIntegral;
  ex coef1;
  /*These are the exponents of the propagators in the integrand.*/
  possymbol* indices = new possymbol[integralfamily.jule];
  generate_symbols(indices, "a", integralfamily.jule);

  for (size_t ii = 0; ii < integralfamily.loopVar.size(); ii++) {
    ibp.clear();

    int next = ii + 1;

    if (ii == integralfamily.loopVar.size() - 1) next = 0;

    for (size_t i = 0; i < integralfamily.props.nops(); i++) {
      coef1 = (diff(integralfamily.props[i], integralfamily.loopVar[ii]) *
               (-indices[i]) * var[next]);
      fs<lst>(coef1, kinematic);
      if (coef1 != 0) {
        newIntegral = new IBPIntegral;
        for (size_t j = 0; j < integralfamily.props.nops(); j++) {
          newIntegral->indices[j] = 0;
        }
        newIntegral->coefficient = coef1;
        newIntegral->indices[i] = 1;
        newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
        ibp.push_back(newIntegral);
      }
    }
    if (var[next] == integralfamily.loopVar[ii]) {
      newIntegral = new IBPIntegral;
      newIntegral->coefficient = dimension;
      for (size_t j = 0; j < integralfamily.props.nops(); j++) {
        newIntegral->indices[j] = 0;
      }
      newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
      ibp.push_back(newIntegral);
    }
    reduce_scal(ibp);

    BaseEquation* equationIBP;
    equationIBP = new BaseEquation(ibp, integralfamily.jule);

    integralfamily.identitiesLEE.push_back(equationIBP);
  }

  ibp.clear();

  for (size_t ii = 0; ii < integralfamily.loopVar.size(); ii++) {
    int next = ii;

    for (size_t i = 0; i < integralfamily.props.nops(); i++) {
      coef1 = (diff(integralfamily.props[i], integralfamily.loopVar[ii]) *
               (-indices[i]) * var[next]);
      fs<lst>(coef1, kinematic);
      if (coef1 != 0) {
        newIntegral = new IBPIntegral;
        for (size_t j = 0; j < integralfamily.props.nops(); j++) {
          newIntegral->indices[j] = 0;
        }
        newIntegral->coefficient = coef1;
        newIntegral->indices[i] = 1;
        newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
        ibp.push_back(newIntegral);
      }
    }
    if (var[next] == integralfamily.loopVar[ii]) {
      newIntegral = new IBPIntegral;
      newIntegral->coefficient = dimension;
      for (size_t j = 0; j < integralfamily.props.nops(); j++) {
        newIntegral->indices[j] = 0;
      }
      newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
      ibp.push_back(newIntegral);
    }
  }
  reduce_scal(ibp);

  BaseEquation* equationIBP;
  equationIBP = new BaseEquation(ibp, integralfamily.jule);

  integralfamily.identitiesLEE.push_back(equationIBP);

  delete[] indices;
}

void Kira::create_LEE_vectors2(vector<possymbol>& var) {
  IBPVG ibp;
  IBPIntegral* newIntegral;
  ex coef1;
  /*These are the exponents of the propagators in the integrand.*/
  possymbol* indices = new possymbol[integralfamily.jule];
  generate_symbols(indices, "a", integralfamily.jule);

  for (size_t iii = 0; iii < var.size(); iii++) {
    //     for(size_t ii = 0; ii<integralfamily.loopVar.size(); ii++){

    ibp.clear();
    for (size_t i = 0; i < integralfamily.props.nops(); i++) {
      coef1 = (diff(integralfamily.props[i], integralfamily.loopVar[0]) *
               (-indices[i]) * var[iii]);
      fs<lst>(coef1, kinematic);
      if (coef1 != 0) {
        newIntegral = new IBPIntegral;
        for (size_t j = 0; j < integralfamily.props.nops(); j++) {
          newIntegral->indices[j] = 0;
        }
        newIntegral->coefficient = coef1;
        newIntegral->indices[i] = 1;
        newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
        ibp.push_back(newIntegral);
      }
    }
    if (var[iii] == integralfamily.loopVar[0]) {
      newIntegral = new IBPIntegral;
      newIntegral->coefficient = dimension;
      for (size_t j = 0; j < integralfamily.props.nops(); j++) {
        newIntegral->indices[j] = 0;
      }
      newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
      ibp.push_back(newIntegral);
    }
    reduce_scal(ibp);

    BaseEquation* equationIBP;
    equationIBP = new BaseEquation(ibp, integralfamily.jule);

    integralfamily.identitiesLEE.push_back(equationIBP);
    //     }
  }

  delete[] indices;
}

void Kira::create_IBP() {
  /*Generate IBP's*/
  logger << "\n*****Generate IBP identities*******************************\n";
  Clock clock;

  create_IBP_helper(integralfamily.allVar);
  create_IBP_helper(externalVar);

  std::sort(integralfamily.identitiesIBP.begin(),
            integralfamily.identitiesIBP.end(), sortSize);

  ofstream myf(
      (outputDir + "/sectormappings/" + integralfamily.name + "/IBP").c_str());

  for (auto i = integralfamily.identitiesIBP.begin();
       i != integralfamily.identitiesIBP.end(); i++) {
    (*i)[0].write_file(myf,integralfamily.name);
  };

  logger << "There are " << integralfamily.identitiesIBP.size() << " "
         << "IBP identities\n";
  logger << "( " << clock.eval_time() << " s )\n";
}

void Kira::create_LEE() {
  /*Generate LEE's*/
  logger << "\n*****Generate LEE identities*******************************\n";
  Clock clock;

  create_LEE_vectors(integralfamily.allVar);
  create_LEE_vectors2(externalVar);

  std::sort(integralfamily.identitiesLEE.begin(),
            integralfamily.identitiesLEE.end(), sortSize);

  ofstream myf(
      (outputDir + "/sectormappings/" + integralfamily.name + "/LEE").c_str());

  for (auto i = integralfamily.identitiesLEE.begin();
       i != integralfamily.identitiesLEE.end(); i++) {
    (*i)[0].write_file(myf,integralfamily.name);
  };

  logger << "There are " << integralfamily.identitiesLEE.size() << " "
         << "LEE vectors identities\n";
  logger << "( " << clock.eval_time() << " s )\n";
}

void Kira::create_LI() {
  /*Generate  Lorentz identities*/
  Clock clock;
  logger << "\n*****Generate LI identities********************************\n";
  IBPVG li;
  IBPIntegral* newIntegral;
  ex coef1, coef2;
  possymbol* indices = new possymbol[integralfamily.jule];
  generate_symbols(indices, "a", integralfamily.jule);
  possymbol* symbolicInd = new possymbol[integralfamily.jule];
  generate_symbols(symbolicInd, "b", integralfamily.jule);

  for (size_t iii = 0; iii < externalVar.size(); iii++) {
    for (size_t ii = iii + 1; ii < externalVar.size(); ii++) {
      li.clear();
      for (size_t i4 = 0; i4 < externalVar.size(); i4++) {
        for (size_t i = 0; i < integralfamily.props.nops(); i++) {
          auto symbIt = find(integralfamily.symbolicIBP.begin(),
                             integralfamily.symbolicIBP.end(), i);

          if (symbIt != integralfamily.symbolicIBP.end())

            coef1 = 2 * (expand(diff(integralfamily.props[i], externalVar[i4]) *
                                (-(indices[i] + symbolicInd[i]))) *
                         externalVar[iii]);
          else

            coef1 = 2 * (expand(diff(integralfamily.props[i], externalVar[i4]) *
                                (-(indices[i]))) *
                         externalVar[iii]);

          coef2 = externalVar[i4] * externalVar[ii];
          fs<lst>(coef1, kinematic);
          fs<lst>(coef2, kinematic);
          if (coef1 * coef2 != 0) {
            newIntegral = new IBPIntegral;
            for (size_t j = 0; j < integralfamily.props.nops(); j++) {
              newIntegral->indices[j] = 0;
            }
            newIntegral->coefficient = coef1 * coef2;
            newIntegral->indices[i] = 1;
            newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
            li.push_back(newIntegral);
          }

          symbIt = find(integralfamily.symbolicIBP.begin(),
                        integralfamily.symbolicIBP.end(), i);

          if (symbIt != integralfamily.symbolicIBP.end())

            coef1 =
                -2 * (expand(diff(integralfamily.props[i], externalVar[i4]) *
                             (-(indices[i] + symbolicInd[i]))) *
                      externalVar[ii]);

          else

            coef1 =
                -2 * (expand(diff(integralfamily.props[i], externalVar[i4]) *
                             (-(indices[i]))) *
                      externalVar[ii]);

          coef2 = externalVar[i4] * externalVar[iii];
          fs<lst>(coef1, kinematic);
          fs<lst>(coef2, kinematic);
          if (coef1 * coef2 != 0) {
            newIntegral = new IBPIntegral;
            for (size_t j = 0; j < integralfamily.props.nops(); j++) {
              newIntegral->indices[j] = 0;
            }
            newIntegral->coefficient = coef1 * coef2;
            newIntegral->indices[i] = 1;
            newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
            li.push_back(newIntegral);
          }
        }
      }
      reduce_scal(li);

      BaseEquation* equationLI;
      equationLI = new BaseEquation(li, integralfamily.jule);
      integralfamily.identitiesLI.push_back(equationLI);
    }
  }

  delete[] indices;
  delete [] symbolicInd;

  std::sort(integralfamily.identitiesLI.begin(),
            integralfamily.identitiesLI.end(), sortSize);

  ofstream myf(
      (outputDir + "/sectormappings/" + integralfamily.name + "/LI").c_str());

  for (ItVE i = integralfamily.identitiesLI.begin();
       i != integralfamily.identitiesLI.end(); i++) {
    (*i)[0].write_file(myf,integralfamily.name);
  };

  logger << "There are " << integralfamily.identitiesLI.size() << " "
         << "LI identities\n";
  logger << "( " << clock.eval_time() << " s )\n";
}
