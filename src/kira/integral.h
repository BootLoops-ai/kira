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

#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <tuple>
#include <unordered_map>

#include "ginac/ginac.h"
#include "pyred/defs.h"

#define SEEDSIZE 128 // Kira can handle up to 128 Propagators
#define SECTOR 0
#define TOPOLOGY 1
#define DENCOUNT 2
#define DOTS 3
#define NUM 4
#define NKEY 5
#define KEYSIZE (SEEDSIZE + 6)

class Fermat;

class SEEDIntegral {
public:
  SEEDIntegral(){};
  SEEDIntegral(unsigned l_Indices_) : l_Indices(l_Indices_){};
  ~SEEDIntegral(){};
  void generate_characteristics(int topology);
  int characteristics[NKEY];
  int indices[SEEDSIZE];
  unsigned length;
  pyred::Weight id;
  int flag2;
  unsigned l_Indices;
  int pass;
  friend std::ostream& operator<<(std::ostream& out, const SEEDIntegral& per);
  friend std::istream& operator>>(std::istream& out, SEEDIntegral& per);
};

class TESTIntegral : public SEEDIntegral {
public:
  TESTIntegral(){};
  void copy(const SEEDIntegral& integral, unsigned l_Indices_);
  void generate_characteristics(int topology);
  std::string coefficientString;
  friend std::ostream& operator<<(std::ostream& out, const TESTIntegral& per);
  friend std::istream& operator>>(std::istream& out, TESTIntegral& per);
};

class IBPIntegral : public SEEDIntegral {
public:
  IBPIntegral(){};
  IBPIntegral(unsigned l_Indices_) { l_Indices = l_Indices_; };
  IBPIntegral(pyred::Weight weights, unsigned l_Indices_);
  ~IBPIntegral(){};
  void copy(const SEEDIntegral& integral, unsigned l_Indices_);
  int compare_indices(IBPIntegral& B);
  void generate_characteristics(int topology);
  GiNaC::ex coefficient;
  std::string coefficientString;
  std::string coefficientIBP[SEEDSIZE + 1];
  friend std::ostream& operator<<(std::ostream& out, const IBPIntegral& per);
  friend std::istream& operator>>(std::istream& out, IBPIntegral& per);
};

typedef std::vector<IBPIntegral*> IBPVG;
typedef IBPVG::iterator ItIBPVG;

typedef std::unordered_map<pyred::Weight, std::tuple<int, int, int> > INTEGMAP;
typedef INTEGMAP::const_iterator INTEGMAPI;

typedef std::unordered_map<pyred::Weight, pyred::Weight> MASTERSMAP;
typedef MASTERSMAP::const_iterator MASTERSMAPI;

class BaseIntegral {
public:
  BaseIntegral(){};
  ~BaseIntegral(){};
  //   void integral_buffer(char buffer[]);
  void copy(const TESTIntegral& integral);
  void copy(const IBPIntegral& integral);
  int characteristics[2];
  std::string coefficientString;
  std::vector<char> coefficientString2;
  unsigned length;
  pyred::Weight id;
  int flag2;
  friend std::ostream& operator<<(std::ostream& out, const BaseIntegral& per);
  friend std::istream& operator>>(std::istream& out, BaseIntegral& per);
};

typedef std::vector<BaseIntegral*> VG;
typedef VG::iterator ItVG;

typedef std::vector<int*> Vintx;
typedef Vintx::iterator ItVintx;

typedef std::vector<int> Vint;
typedef Vint::iterator ItVint;

typedef std::vector<std::vector<int> > VVint;
typedef VVint::iterator ItVVint;
typedef VVint::reverse_iterator RItVVint;

class BaseEquation;

class IBPEquation;

typedef std::vector<BaseEquation*> VE;
typedef VE::iterator ItVE;

typedef std::vector<IBPEquation*> IBPVE;
typedef IBPVE::iterator ItIBPVE;

class BaseEquation {
public:
  BaseEquation();
  BaseEquation(unsigned l_Equation_, unsigned l_Indices_);
  BaseEquation(IBPVG& ibp, unsigned l_Indices_);
  BaseEquation(IBPVG& ibp, unsigned l_Indices_, INTEGMAP& integralMap,
               MASTERSMAP&);
  BaseEquation(const BaseEquation& ibp, pyred::Weight weights,
               INTEGMAP& integralMap, const int topology, Fermat*& fermat,
               int id, MASTERSMAP&);
  ~BaseEquation();
  void delete_IBP();
  void write_file(std::ofstream& soad,std::string& topology);
  void eliminate_zeros(MASTERSMAP& mastersMap);
  void eliminate_zeros(IBPVG& ibp, MASTERSMAP& mastersMap);
  int sort();
  void collect_integrals();
  void collect_integrals(IBPVG& ibp);
  void dismember_coef();
  void plant_seed(BaseEquation*& ibp, int*& seed);
  BaseIntegral* equation;
  IBPIntegral* equationIBP;
  std::vector<IBPIntegral> VECTORequationIBP;
  std::vector<SEEDIntegral> VECTORequationSEED;
  std::vector<TESTIntegral> VECTORequationTEST;
  unsigned l_Equation;
  unsigned l_Indices;
};

bool sort_rules3(BaseIntegral l[], BaseIntegral r[]);
bool sort_rules4(BaseEquation l[], BaseEquation r[]);

class symmetries {
public:
  int ing[SEEDSIZE];
  int sector;
  int externalSymmetry;
  GiNaC::lst subst;
  GiNaC::lst ruleS;
  int det;
  int nOfProps;
  int topology;
  int symDOTS;
};

class shifts {
public:
  GiNaC::lst ruleMom;
  GiNaC::lst ruleMand;
};

#endif // INTEGRAL_H
