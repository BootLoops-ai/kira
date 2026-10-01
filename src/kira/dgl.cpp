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

#include <cstddef>

#include "kira/kira.h"
#include "kira/ReadYamlFiles.h"
#include "kira/tools.h"
#include "pyred/coeff_int.h"
#include "pyred/defs.h"
#include "pyred/integrals.h"
#include "pyred/interface.h"
#include "pyred/parser.h"

// using namespace pyred;
using namespace std;
using namespace YAML;
using namespace GiNaC;

static Loginfo& logger = Loginfo::instance();


void Kira::dgl_vs_scalar_products(){

  cout << "construct DGLs with respect to the scalar products of external momenta"  << endl;
  // construct DGL with respect to scalar products of external momenta

  vector<possymbol> intvar = externalVar;

  lst mZS;
  for (size_t i = 0; i < externalVar.size(); i++) {
    for (size_t j = 0; j < externalVar.size(); j++) {
      ex token = externalVar[i] * externalVar[j];
      fs<lst>(token, kinematic);
      mZS.append(token);
    }
  }

  cout << mZS << endl;

  matrix mZS2(externalVar.size(), externalVar.size(), mZS);

  matrix mZS2inverse = mZS2.inverse();
  ex mZS2det = mZS2.determinant();

        cout << mZS2inverse << endl;
        cout << mZS2inverse(0,0) << endl;

  /*These are the exponents of the propagators in the integrand.*/
  possymbol* indices = new possymbol[integralfamily.jule];
  generate_symbols(indices, "a", integralfamily.jule);
  IBPVG ibp;
  ex coef1;

  {
    IBPIntegral* newIntegral;

    for (size_t iii = 0; iii < intvar.size(); iii++) {
      ibp.clear();
      for (size_t l = 0; l < intvar.size(); l++) {
        for (size_t ip = 0; ip < integralfamily.props.nops(); ip++) {
          coef1 = (diff(integralfamily.props[ip], intvar[iii]) *
                    (-indices[ip]) * intvar[l]);
          fs<lst>(coef1, kinematic);

          // 	    cout << coef1 << endl;

          if (coef1 != 0) {
            newIntegral = new IBPIntegral;
            for (size_t j = 0; j < integralfamily.props.nops(); j++) {
              newIntegral->indices[j] = 0;
            }
            newIntegral->coefficient = coef1 * mZS2inverse(l, iii) / 2;

            newIntegral->indices[ip] = 1;
            newIntegral->characteristics[TOPOLOGY] =
                integralfamily.topology;
            ibp.push_back(newIntegral);
          }
        }
      }
      reduce_scal(ibp);
      BaseEquation* equationIBP;
      equationIBP = new BaseEquation(ibp, integralfamily.jule);

      integralfamily.identitiesDGLmom.push_back(equationIBP);

      cout << "1) number of terms: " << equationIBP[0].l_Equation << endl;

    }
  }

  {
    IBPIntegral* newIntegral;

    for (size_t iii = 0; iii < intvar.size(); iii++) {
      for (size_t k = iii + 1; k < intvar.size(); k++) {
        ibp.clear();

        for (size_t l = 0; l < intvar.size(); l++) {
          for (size_t ip = 0; ip < integralfamily.props.nops(); ip++) {
            coef1 = (diff(integralfamily.props[ip], intvar[iii]) *
                      (-indices[ip]) * intvar[l]);
            fs<lst>(coef1, kinematic);

            // 		cout << coef1 << endl;

            if (coef1 != 0) {
              newIntegral = new IBPIntegral;
              for (size_t j = 0; j < integralfamily.props.nops(); j++) {
                newIntegral->indices[j] = 0;
              }
              newIntegral->coefficient = coef1 * mZS2inverse(l, k);
              newIntegral->indices[ip] = 1;
              newIntegral->characteristics[TOPOLOGY] =
                  integralfamily.topology;
              ibp.push_back(newIntegral);
            }
          }
        }
        reduce_scal(ibp);
        BaseEquation* equationIBP;
        equationIBP = new BaseEquation(ibp, integralfamily.jule);

        integralfamily.identitiesDGLmom.push_back(equationIBP);
        cout << "2) number of terms: " << equationIBP[0].l_Equation << endl;
      }
    }
  }

  delete[] indices;

  ofstream myf(
      (outputDir + "/sectormappings/" + integralfamily.name + "/DGLmomenta")
          .c_str());

  for (auto i = integralfamily.identitiesDGLmom.begin();
        i != integralfamily.identitiesDGLmom.end(); i++) {
    (*i)[0].write_file(myf,integralfamily.name);
  }

  cout << "done" << endl << endl;
}


void Kira::generate_dgl() {

  for (unsigned itC = 0; itC < collectReductions.size(); itC++) {

    integralfamily = topology[collectReductions[itC]];

    int numX = 0;

    lst invariants4DGL = invariants4sym;
    lst invariantsListDGL = invariantsList;
    vector<possymbol> symbolInvariantsDGL = symbolInvariants;

    cout << "get variables for DGL" << endl;
    {
      ex coef1;
      for (size_t iii = 0; iii < invar.size(); iii++) {
        cout << invar[iii] << endl;

        for (size_t ip = 0; ip < integralfamily.props.nops(); ip++) {
          if (diff(integralfamily.props[ip].expand(), invar[iii]) == 0)
            continue;

          for (int num = 1; num < 4; num++) {
            coef1 = coeff(integralfamily.props[ip].expand(), invar[iii], num);

            if (coef1 != 0) {
              int skip = 0;
              for (auto itXX : invariants4sym) {
                if (pow(invar[iii], num) == itXX.lhs()) {
                  skip = 1;
                  break;
                }
              }
              if (skip) break;

              cout << pow(invar[iii], num) << endl;

              ex newRule = pow(invar[iii], num) ==
                           get_symbol("xKiraDGL" + to_string(numX));

              invariantsListDGL.append(
                  get_symbol("xKiraDGL" + to_string(numX)));
              symbolInvariantsDGL.push_back(
                  get_symbol("xKiraDGL" + to_string(numX)));
              numX++;
              invariants4DGL.append(newRule);

              break;
            }
          }
        }
      }
    }
    cout << "invariants4DGL: " << invariants4DGL << endl;
    cout << "invariantsListDGL: " << invariantsListDGL << endl;

    dgl_vs_scalar_products();

    cout << "Take derivative of xKira... with respect to scalar products of "
            "external momenta"
         << endl;
    std::vector<std::vector<GiNaC::ex> > coefDGL;

    // Derivatives with respect to p1^2, p1*p2, ... in xKira..., ...
    {
      bS = new possymbol[kinematic.nops()];
      generate_symbols(bS, "bS", static_cast<int>(kinematic.nops()));

      lst kinematic4DGL;

      lst unkV;

      for (auto itV : invariantsList) {

        unkV.append(itV);
      }

      for (unsigned itV = 0; itV < kinematic.nops(); itV++) {
        unkV.append(bS[itV]);
      }

      for (size_t itE = 0; itE < kinematicR.nops(); itE++)
        kinematic4DGL.append(kinematicR[itE].subs(invariants4sym));

      ex solutionDGL = lsolve(kinematic4DGL, unkV);

      lst kinInv2D;
      for (size_t itX = 0; itX < invariantsList.nops(); itX++) {
        cout << solutionDGL[itX] << endl;
        kinInv2D.append(solutionDGL[itX]);
      }

      for (auto itX : kinInv2D) {
        vector<ex> coefDGLtoken;
        for (unsigned itV = 0; itV < kinematic.nops(); itV++) {
          ex diffToken = diff(itX.rhs(), bS[itV]);

          if (diffToken != 0) {
            cout << pow(diffToken, -1) << endl;
            coefDGLtoken.push_back(pow(diffToken, -1));
          }
          else {
            coefDGLtoken.push_back(diffToken);
            cout << diffToken << endl;
          }
        }
        coefDGL.push_back(coefDGLtoken);
        cout << endl;
      }

      delete[] bS;
    }

    cout << "done" << endl << endl;

    cout << "construct DGL with respect to all variables (no chain rule)"
         << endl;
    // construct DGLs partly with respect to all variables
    {
      IBPIntegral* newIntegral;
      possymbol* indices = new possymbol[integralfamily.jule];
      generate_symbols(indices, "a", integralfamily.jule);
      IBPVG ibp;
      ex coef1;

      for (size_t iii = 0; iii < symbolInvariantsDGL.size(); iii++) {
        ex coefChain = 0;
        for (auto itY : invar)
          coefChain += diff(invariants4DGL[iii].lhs(), itY);
        cout << coefChain << endl;

        ibp.clear();

        for (size_t ip = 0; ip < integralfamily.props.nops(); ip++) {
          coef1 = (diff(integralfamily.props[ip].expand().subs(
                            invariants4DGL, subs_options::algebraic),
                        symbolInvariantsDGL[iii]) *
                   (-indices[ip])) *
                  coefChain;

          fs<lst>(coef1, kinematic);

          cout << "coeff: " << coef1 << endl;

          if (coef1 != 0) {
            newIntegral = new IBPIntegral;
            for (size_t j = 0; j < integralfamily.props.nops(); j++) {
              newIntegral->indices[j] = 0;
            }
            newIntegral->coefficient = coef1;
            newIntegral->indices[ip] = 1;
            newIntegral->characteristics[TOPOLOGY] = integralfamily.topology;
            ibp.push_back(newIntegral);
          }
        }
        reduce_scal(ibp);
        BaseEquation* equationIBP;
        equationIBP = new BaseEquation(ibp, integralfamily.jule);

        integralfamily.identitiesDGLmasses.push_back(equationIBP);
        cout << "number of terms: " << equationIBP[0].l_Equation << endl;
      }

      delete[] indices;
    }
    cout << "done masses: " << integralfamily.identitiesDGLmasses.size() << endl
         << endl;

    cout << "construct DGL with respect to the xKira0 (only chain rule)"
         << endl;

    // construct DGL with respect to the invariants
    {
      IBPVG ibp;
      int countXDGL = 0;
      for (auto cDGL : coefDGL) {
        ex coefChain = 0;
        for (auto itY : invar)
          coefChain += diff(invariants4sym[countXDGL].lhs(), itY);
        cout << coefChain << endl;

        ibp.clear();
        for (size_t itCoef = 0; itCoef < cDGL.size(); itCoef++) {
          if (cDGL[itCoef] == 0) continue;

          auto iEq = integralfamily.identitiesDGLmom[itCoef];

          IBPIntegral* newIntegral;
          cout << itCoef << " : " << (iEq)[0].l_Equation << endl;
          ;

          for (unsigned i = 0; i < (iEq)[0].l_Equation; i++) {
            newIntegral = new IBPIntegral;

            for (unsigned j = 0; j < (iEq)[0].l_Indices; j++) {
              newIntegral->indices[j] = (iEq)[0].equationIBP[i].indices[j];
            }

            // 	    newIntegral->coefficientString =
            // "("+(iEq)[0].equationIBP[i].coefficientString+")*(" +
            // something_string(cDGL[itCoef])+")";

            newIntegral->coefficient =
                (iEq)[0].equationIBP[i].coefficient * cDGL[itCoef] * coefChain;

            newIntegral->characteristics[TOPOLOGY] =
                (iEq)[0].equationIBP[i].characteristics[TOPOLOGY];
            ibp.push_back(newIntegral);
          }
        }
        countXDGL++;
        reduce_scal(ibp);

        BaseEquation* equationIBP;
        equationIBP = new BaseEquation(ibp, integralfamily.jule);

        integralfamily.identitiesDGLxKira.push_back(equationIBP);
      }
      ofstream myf(
          (outputDir + "/sectormappings/" + integralfamily.name + "/DGLxKira")
              .c_str());

      for (auto i = integralfamily.identitiesDGLxKira.begin();
           i != integralfamily.identitiesDGLxKira.end(); i++) {
        (*i)[0].write_file(myf,integralfamily.name);
      }
    }

    cout << "done: " << integralfamily.identitiesDGLxKira.size() << endl
         << endl;

    // put mass derivatives and kinematics together if possible.

    cout << "Generate the final DGL for the integrals" << endl;

    {
      IBPVG ibp;
      for (size_t itDGL = 0; itDGL < integralfamily.identitiesDGLxKira.size();
           itDGL++) {

        cout << itDGL << endl;

        cout << invariants4sym.nops() << endl;

        string massDGL;
        for (auto itY : invar) {
          cout << invariants4sym[itDGL].lhs() << endl;
          cout << diff(invariants4sym[itDGL].lhs(), itY) << endl;
          if (diff(invariants4sym[itDGL].lhs(), itY) != 0) {
            massDGL = something_string(itY);
            break;
          }
        }
        cout << "mass: " << massDGL << endl;

        if (integralfamily.identitiesDGLmasses[itDGL][0].l_Equation == 0 &&
            integralfamily.identitiesDGLxKira[itDGL][0].l_Equation == 0) {
          cout << "none: " << itDGL << endl;
          continue;
        }
        else if (integralfamily.identitiesDGLxKira[itDGL][0].l_Equation == 0) {
          integralfamily.identitiesDGL.push_back(
              make_pair(massDGL, integralfamily.identitiesDGLmasses[itDGL]));

          cout << "masses alone: " << itDGL << endl;

          continue;
        }
        else if (integralfamily.identitiesDGLmasses[itDGL][0].l_Equation == 0) {
          integralfamily.identitiesDGL.push_back(
              make_pair(massDGL, integralfamily.identitiesDGLxKira[itDGL]));

          cout << "chain rule alone: " << itDGL << endl;

          continue;
        }
        cout << "combine chain rule and masses: " << itDGL << endl;

        ibp.clear();
        {
          auto iEq = integralfamily.identitiesDGLmasses[itDGL];

          IBPIntegral* newIntegral;
          cout << itDGL << " : " << (iEq)[0].l_Equation << endl;
          ;

          for (unsigned i = 0; i < (iEq)[0].l_Equation; i++) {
            newIntegral = new IBPIntegral;

            for (unsigned j = 0; j < (iEq)[0].l_Indices; j++) {
              newIntegral->indices[j] = (iEq)[0].equationIBP[i].indices[j];
            }

            newIntegral->coefficient = (iEq)[0].equationIBP[i].coefficient;

            newIntegral->characteristics[TOPOLOGY] =
                (iEq)[0].equationIBP[i].characteristics[TOPOLOGY];
            ibp.push_back(newIntegral);
          }
        }
        {
          auto iEq = integralfamily.identitiesDGLxKira[itDGL];

          IBPIntegral* newIntegral;
          cout << itDGL << " : " << (iEq)[0].l_Equation << endl;
          ;

          for (unsigned i = 0; i < (iEq)[0].l_Equation; i++) {
            newIntegral = new IBPIntegral;

            for (unsigned j = 0; j < (iEq)[0].l_Indices; j++) {
              newIntegral->indices[j] = (iEq)[0].equationIBP[i].indices[j];
            }

            newIntegral->coefficient = (iEq)[0].equationIBP[i].coefficient;

            newIntegral->characteristics[TOPOLOGY] =
                (iEq)[0].equationIBP[i].characteristics[TOPOLOGY];
            ibp.push_back(newIntegral);
          }
        }

        reduce_scal(ibp);

        BaseEquation* equationIBP;
        equationIBP = new BaseEquation(ibp, integralfamily.jule);

        integralfamily.identitiesDGL.push_back(make_pair(massDGL, equationIBP));
      }
      cout << "hmm" << endl;

      for (size_t itDGL = integralfamily.identitiesDGLxKira.size();
           itDGL < integralfamily.identitiesDGLmasses.size(); itDGL++) {
        string massDGL;
        for (auto itY : invar) {
          if (diff(invariants4DGL[itDGL].lhs(), itY) != 0) {
            massDGL = something_string(itY);
            break;
          }
        }
        cout << massDGL << endl;

        if (integralfamily.identitiesDGLmasses[itDGL][0].l_Equation != 0) {
          integralfamily.identitiesDGL.push_back(
              make_pair(massDGL, integralfamily.identitiesDGLmasses[itDGL]));

          cout << "masses alone: " << itDGL << endl;
        }
      }

      ofstream myf(
          (outputDir + "/sectormappings/" + integralfamily.name + "/DGL")
              .c_str());

      for (auto i = integralfamily.identitiesDGL.begin();
           i != integralfamily.identitiesDGL.end(); i++) {
        (*i).second[0].write_file(myf,integralfamily.name);
      }
    }

    topology[collectReductions[itC]] = integralfamily;
  }
}

std::vector<std::tuple<std::string, int, std::vector<int>, std::string> >
Kira::read_seeds_dgl(std::string& itFile) {

  ifstream mastersInput;
  mastersInput.open(itFile);

  std::vector<std::tuple<std::string, int, std::vector<int>, std::string> >
      seedsDGL;

  string line;

  logger.set_level(2);

  while (pyred::getline_normalized(mastersInput, line)) {
    size_t found;

    string massStr;
    string topo;
    string indicesStr;
    string coefStr;

    if ((found = line.find_first_of(" ")) != string::npos) {
      massStr = line.substr(0, found);
    }

    line = line.substr(found + 1, -1);

    // get string indices
    if (line.size() < 3) continue;

    line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

    // convert string indices to integer array
    int indicesCounter;
    vector<int> inputSeed;

    if ((found = line.find_first_of("#")) != string::npos) {
      logger << "skeep this line: " << line << "\n";
      continue;
    }

    if ((found = line.find_first_of("[")) == string::npos) {
      logger << "skeep this line: " << line << "\n";
      continue;
    }
    istringstream(line.substr(0, found)) >> topo;

    line.erase(line.begin(), line.begin() + found + 1);

    if ((found = line.find_first_of("]")) == string::npos) {
      logger << "skeep this line: " << line << "\n";
      continue;
    }

    indicesStr = line.substr(0, found);

    coefStr = line.substr(found + 1, -1);

    if (coefStr[0] == '*') {
      coefStr = line.substr(found + 2, -1);
    }

    if ((found = indicesStr.find_first_of(",")) != string::npos) {
      bool has_only_digits = (indicesStr.substr(0, found).find_first_not_of(
                                  "0123456789-+") == string::npos);
      if (!has_only_digits) {
        istringstream(indicesStr.substr(0, found)) >> topo;
        indicesStr = indicesStr.substr(found + 1);
      }
    }

    while ((found = indicesStr.find_first_of(",")) != string::npos) {
      istringstream(indicesStr.substr(0, found)) >> indicesCounter;
      inputSeed.push_back(indicesCounter);
      indicesStr = indicesStr.substr(found + 1);
    }

    istringstream(indicesStr) >> indicesCounter;
    inputSeed.push_back(indicesCounter);

    cout << massStr << endl;
    cout << coefStr << endl;
    cout << topo << " ";
    cout << topology[topo].topology << endl;

    for (auto itY : inputSeed) {
      cout << itY << " ";
    }
    cout << endl;
    cout << endl;

    seedsDGL.push_back(
        make_tuple(massStr, topology[topo].topology, inputSeed, coefStr));
  }

  logger.set_level(1);

  return seedsDGL;
}

void Kira::insert_seeds2DGL(
    std::vector<std::tuple<std::string, int, std::vector<int>, std::string> >&
        seedsDGL) {

  possymbol* indices = new possymbol[integralfamily.jule];
  generate_symbols(indices, "a", integralfamily.jule);

  ofstream myf(
          (outputDir + "/results/" + integralfamily.name + "/dgls")
              .c_str());

  for (auto itG : seedsDGL) {

    integralfamily = topology[collectReductions[get<1>(itG)]];
    vector<int> propsindex = get<2>(itG);
    parser symbolReader(GiNaCSymbols);

    ex coefDGLExtra = symbolReader(get<3>(itG));

    ex coefDGLExtraDD = diff(coefDGLExtra, get_symbol(get<0>(itG)));

    if (coefDGLExtra == 0) continue;

    cout << get<0>(itG) << endl;

    cout << "coefficient to take derivative: " << coefDGLExtra << endl;

    cout << "derivative taken in the coefficient: " << coefDGLExtraDD << endl;

    lst aha;

    for (size_t i = 0; i < propsindex.size(); i++) {
      aha.append(indices[i] == propsindex[i]);
    }

    cout << "indices to substitute: " << aha << endl;
    cout << "number of possible DGL: " << integralfamily.identitiesDGL.size()
         << endl;
    int pickDGL = -1;

    for (size_t itDGL = 0; itDGL < integralfamily.identitiesDGL.size();
         itDGL++) {
      if (get<0>(itG) == integralfamily.identitiesDGL[itDGL].first) {
        // 	cout << get<0>(itG) << " : " << itDGL << endl;
        pickDGL = itDGL;
        break;
      }
    }

    if (pickDGL == -1) continue;

    BaseEquation* token = integralfamily.identitiesDGL[pickDGL].second;

    int equationLength = 0;

    for (unsigned i = 0; i < token[0].l_Equation; i++) {

      ex coefDGL = token[0].equationIBP[i].coefficient;
      fs<lst>(coefDGL, aha);

      coefDGL = coefDGL * coefDGLExtra;

      if (coefDGL == 0) continue;

      cout << collectReductions[get<1>(itG)];
      cout << "[";
      for (unsigned j = 0; j < token[0].l_Indices; j++) {
        cout << (propsindex[j] + token[0].equationIBP[i].indices[j]);
        if (j < token[0].l_Indices - 1) {
          cout << ",";
        }
      }
      cout << "]*(" << coefDGL.normal() << ")" << endl;


      myf << collectReductions[get<1>(itG)];
      myf << "[";
      for (unsigned j = 0; j < token[0].l_Indices; j++) {
        myf << (propsindex[j] + token[0].equationIBP[i].indices[j]);
        if (j < token[0].l_Indices - 1) {
          myf << ",";
        }
      }
      myf << "]*(" << coefDGL.normal() << ")" << endl;

      equationLength++;
    }

    if (coefDGLExtra != 0) {

      cout << collectReductions[get<1>(itG)];
      cout << "[";
      for (unsigned j = 0; j < token[0].l_Indices; j++) {
        cout << (propsindex[j]);
        if (j < token[0].l_Indices - 1) {
          cout << ",";
        }
      }
      cout << "]*(" << coefDGLExtra.normal() << ")" << endl;


      myf << collectReductions[get<1>(itG)];
      myf << "[";
      for (unsigned j = 0; j < token[0].l_Indices; j++) {
        myf << (propsindex[j]);
        if (j < token[0].l_Indices - 1) {
          myf << ",";
        }
      }
      myf << "]*(" << coefDGLExtra.normal() << ")" << endl;
    }



    cout << "Laenge: " << equationLength << endl;
    cout << endl;
    myf << endl;
  }

  delete[] indices;
}
