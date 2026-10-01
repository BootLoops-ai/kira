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

#ifndef READYAMLFILES_H_
#define READYAMLFILES_H_

#include <cstddef>
#include <algorithm>
#include <tuple>

#include "ginac/ginac.h"
#include "pyred/config.h"
#include "pyred/integrals.h"
#include "yaml-cpp/yaml.h"
#include "kira/exceptions.h"
#include "kira/integral.h"
#include "kira/tools.h"

#define SIZESYM 64
// #include "trivial_sym.h"

typedef std::vector<symmetries> SYM;
typedef SYM::iterator ItSYM;

typedef std::vector<BaseEquation*> VE;
typedef VE::iterator ItVE;

struct Jobs {
  std::vector<std::tuple<std::vector<std::string> /*topologies*/,
                         std::vector<std::string> /*sectors*/, int /*rmax*/,
                         int /*smax*/, int /*dmax*/> >
      reductSpec;

  std::vector<std::tuple<std::vector<std::string> /*topologies*/,
                         std::vector<std::string> /*sectors*/, int /*rmax*/,
                         int /*smax*/, int /*dmax*/> >
      selectSpec;

  std::vector<std::pair<std::string, std::string> > sector2Reduce;
  std::vector<std::pair<std::string, std::string> > den;
  std::vector<std::pair<std::string, std::string> > num;
  std::vector<std::pair<std::string, std::string> > mandatoryFile;
  std::vector<std::string> mandatoryFileVector;
#ifdef KIRAFIREFLY
  std::string ff_recon;
  std::string factor_scan;
#endif
  std::vector<std::vector<std::string> > mandatoryRec;
  std::vector<std::vector<std::string> > optionalRec;
  std::string masters;
  std::string extraRelations;
  std::string specialname;
  std::string runBacksubstitution;
  std::string runTriangular;
  std::string pyredDatabase;
  std::string runInitiate;
  int level;
  int config = 1;
  int integralOrdering;
  bool check_masters_ = false;
  bool autoMode = false;
  std::string runSymmetries;
  std::string dataFile;
  std::string outputDir;
  std::string numericalPoints;
  std::size_t pointsPerFile = 1000;
  std::string writeNumericalSystem;
  std::string conditional;
  std::string LIflag;
  std::string symbols2numsFlag;
  std::string algebraicReconstruction;
  unsigned resumeRun;
  std::multimap<std::string, std::vector<int> > selectMastersReduction;
  std::string iterativeReduction;
  std::vector<std::string> trimmedReduction;
  std::vector<std::string> fileDenominators;
  std::vector<std::string> fileAmplitude;
  int weightMode;
  std::vector<std::string> filePrefactors;
  std::string inputSystem;
  std::tuple<std::vector<std::string>, std::size_t, int> inputSystemTuple;
  std::vector<std::tuple<std::vector<std::string>,std::vector<std::tuple<uint32_t,int>>,int,int>> truncateSP;
                         /*tuple<topologies,tspexceptions vector<sector,s>,sp,lvl>*/
  int trcMore, trcMoreD;

};

struct Kira2File {
  std::vector<std::pair<std::string, std::string> > target;
  std::vector<std::string> targetVector;
  std::vector<std::vector<std::string> > mandatoryRec;
  std::string inputDir;
  std::string reconstructMass;
  std::vector<std::tuple<std::vector<std::string> /*topologies*/,
                         std::vector<std::string> /*sectors*/, int /*rmax*/,
                         int /*smax*/, int /*dmax*/> >
      selectSpec;
};

struct Merge {
  std::vector<std::string> files2merge;
  std::string outputDir;
};

struct Dgl {
  std::vector<std::string> filesDGL;
};

struct Auxiliary {
  std::vector<std::pair<std::string,std::string>> userNumerics;
  std::tuple<std::string,std::vector<uint32_t> > permutation;
  std::string name;
  std::tuple<std::vector<std::tuple<uint32_t,int>>,int,int> truncateSP;
                         /*tuple<tspexceptions vector<sector,s>,sp,lvl>*/
  //int truncateSP;
  std::vector<uint32_t> truncateSectors, truncateSectors2;
  int trcMore, trcMoreD, splvl;
  std::vector<std::pair<std::string, std::string>> randvars;
  std::string setPrimeNNumber;
};

struct Integral_F {
  std::string name;
  std::vector<std::string> loop;
  std::vector<GiNaC::possymbol> loopVar;
  std::vector<GiNaC::possymbol> allVar;
  GiNaC::lst loop2loop, loop2loop2, loopVarList, loopVarList2;
  GiNaC::lst props, propsMomFlowA;
  std::vector<GiNaC::lst> scal2Props, propsMomFlowB;
  uint32_t jule;
  GiNaC::possymbol* propSymb;
  GiNaC::lst scalprod;
  GiNaC::ex G, Gsym, U, F; // U + F Polynomial
  GiNaC::lst FPolynom;
  std::vector<uint32_t>*mask, *allowSector;
  std::set<int>* skipSector;
  std::vector<int> zeroSector;

  VE identitiesIBP, identitiesLI, identitiesDGLmom, identitiesLEE,
      identitiesDGLmasses, identitiesDGLxKira;

  std::vector<std::pair<std::string, BaseEquation*> > identitiesDGL;

  SYM *symVec, *relVec, *symVecReverse;
  std::vector<int> sector2Reduce, bound;
  int biggestBound;
  int biggestSector2Reduce;
  std::vector<unsigned> lowestSectors;
  int topology;
  std::vector<std::pair<std::string, std::string> > propagator;
  int* invarID;
  std::vector<uint32_t> topLevelSectors;
  std::vector<int> userZeroSectors;
  std::vector<int> cutProps;
  unsigned maskCut;
  std::vector<size_t> permute;
  std::string permuteOption;
  std::vector<std::tuple<std::vector<uint32_t>, int, int, int> > reductSpec;
  std::string magic_relations;
  std::unordered_map<
      int, std::vector<std::tuple<size_t, std::vector<std::vector<int> >, int, int> > >
      symmetries;
  std::unordered_map<
      int, std::vector<std::tuple<size_t, std::vector<std::vector<int> >, int, int> > >
      symmetriesCrossed;
  int propsMomentaFlowMask;
  std::vector<int> symbolicIBP;
  std::vector<std::pair<uint32_t, uint32_t> > sectorweight;
  std::vector<std::tuple<std::vector<std::tuple<uint32_t,int>>,int,int>> truncateSP;
                         /*tuple<tspexceptions vector<sector,s>,sp,lvl>*/
  //std::vector<uint32_t> truncateSectors, truncateSectors2;
  //int truncateSP, splvl;
};

struct Kinematics {
  std::vector<std::string> im;
  std::vector<std::string> om;
  std::string rpl;
  std::pair<std::string, std::string> mc;
  std::vector<std::pair<std::string, int> > ki;
  std::vector<GiNaC::possymbol> invariantsPlaceholder;
  GiNaC::lst invariantsReplacement, invariantsReplacementRev;
  std::vector<std::pair<std::pair<std::string, std::string>, std::string> > sr;
};

namespace YAML {
template <>
struct convert<Jobs> {
  static bool decode(const Node& node, Jobs& rhs) {
    for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
      std::string token = (*it).first.as<std::string>();
      if (token != "algebraic_reconstruction" &&
          token != "input_system" &&
          token != "reduce" &&
          token != "sector_selection" &&
          token != "identities" &&
          token != "select_integrals" &&
          token != "select_masters" &&
          token != "preferred_masters" &&
          token != "extra_relations" &&
          token != "check_masters" &&
          token != "run_triangular" &&
          token != "run_back_substitution" &&
          token != "integral_ordering" &&
          token != "run_initiate" &&
          token != "run_symmetries" &&
          token != "alt_dir" &&
          token != "data_file" &&
          token != "write_numerical_system" &&
          token != "select_masters_reduction" &&
          token != "create_denominator_database" &&
          token != "amplitude_translate" &&
          token != "weight_notation" &&
          token != "conditional" &&
          token != "set_zero_sectors" &&
          token != "lorenz_invariance" &&
          token != "symbols_2_nums" &&
          token != "iterative_reduction" &&
          token != "generate_input" &&
          token != "truncate_more" &&
          token != "truncate_more_dots" &&
          token != "truncate_sp" &&
          token != "reduce_automatic" &&
          token != "numerical_points"

#ifndef KIRAFIREFLY
          && token != "pyred_database") {
#else
          && token != "insert_prefactors"
          && token != "pyred_database"
          && token != "run_firefly"
          && token != "factor_scan") {
#endif
        throw ExceptionConfigFiles("job file: Option '" + token + "' unkown");
      }
    }

    if (node["run_triangular"])
      rhs.runTriangular = node["run_triangular"].as<std::string>();
    if (node["run_back_substitution"])
      rhs.runBacksubstitution = node["run_back_substitution"].as<std::string>();
    if (node["pyred_database"])
      rhs.pyredDatabase = node["pyred_database"].as<std::string>();
    if (node["run_initiate"])
      rhs.runInitiate = node["run_initiate"].as<std::string>();
    if (node["run_symmetries"])
      rhs.runSymmetries = node["run_symmetries"].as<std::string>();
    if (node["data_file"]) rhs.dataFile = node["data_file"].as<std::string>();
    if (node["alt_dir"]) rhs.outputDir = node["alt_dir"].as<std::string>();
    if (node["numerical_points"]) {
      if (node["numerical_points"].size() == 0) {
        rhs.numericalPoints = node["numerical_points"].as<std::string>();
      }
      else {
        if (node["numerical_points"]["file"]) {
          rhs.numericalPoints = node["numerical_points"]["file"].as<std::string>();
        }
        else {
          throw ExceptionConfigFiles("job file: Option 'numerical_points' expects an input file");
        }

        if (node["numerical_points"]["points_per_output_file"]) {
          rhs.pointsPerFile = node["numerical_points"]["points_per_output_file"].as<std::size_t>();
        }
      }
    }
    if (node["write_numerical_system"])
      rhs.writeNumericalSystem =
          node["write_numerical_system"].as<std::string>();
    if (node["conditional"])
      rhs.conditional = node["conditional"].as<std::string>();

    if (node["lorenz_invariance"]) {
      rhs.LIflag = node["lorenz_invariance"].as<std::string>();
    }
    rhs.symbols2numsFlag = "true";
    if (node["symbols_2_nums"]) {
      rhs.symbols2numsFlag = node["symbols_2_nums"].as<std::string>();
    }

    if (node["algebraic_reconstruction"]) {
      rhs.algebraicReconstruction =
          node["algebraic_reconstruction"].as<std::string>();
    }

    if (node["reduce_automatic"])
      rhs.autoMode = node["reduce_automatic"].as<bool>();

    if (node["reduce"]) {
      const Node& nodeibp = node["reduce"];

      for (unsigned i = 0; i < nodeibp.size(); i++) {
        std::vector<std::string> topologies;
        std::vector<std::string> sectors;
        int rmax = -1;
        int smax = -1;
        int dmax = std::numeric_limits<int>::max();

        if (nodeibp[i]["topologies"]) {
          for (unsigned itT = 0; itT != nodeibp[i]["topologies"].size();
               itT++) {
            topologies.push_back(
                nodeibp[i]["topologies"][itT].as<std::string>());
          }
        }

        if (nodeibp[i]["sectors"]) {
          for (unsigned itT = 0; itT != nodeibp[i]["sectors"].size(); itT++){

            sectors.push_back(nodeibp[i]["sectors"][itT].as<std::string>());
          }
        }

        if (nodeibp[i]["r"] && nodeibp[i]["r"].size() == 0)
          rmax = nodeibp[i]["r"].as<int>();

        if (nodeibp[i]["s"] && nodeibp[i]["s"].size() == 0)
          smax = nodeibp[i]["s"].as<int>();

        if (nodeibp[i]["d"] && nodeibp[i]["d"].size() == 0)
          dmax = nodeibp[i]["d"].as<int>();

        if (rmax < 0) {
          throw ExceptionConfigFiles("job file: rmax in 'r: rmax' should be positive");
        }
        if (smax < 0) {
          throw ExceptionConfigFiles("job file: smax in 's: smax' should be positive");
        }
        if (dmax < 0) {
          throw ExceptionConfigFiles("job file: dmax in 'd: dmax' should be positive");
        }

        if (dmax == std::numeric_limits<int>::max()) dmax = -1;

        rhs.reductSpec.push_back(
            std::make_tuple(topologies, sectors, rmax, smax, dmax));

        if (nodeibp[i]["r"] && nodeibp[i]["r"].size() == 2)
          rhs.den.push_back(make_pair(nodeibp[i]["r"][0].as<std::string>(),
                                      nodeibp[i]["r"][1].as<std::string>()));
        if (nodeibp[i]["s"] && nodeibp[i]["s"].size() == 2)
          rhs.num.push_back(make_pair(nodeibp[i]["s"][0].as<std::string>(),
                                      nodeibp[i]["s"][1].as<std::string>()));
      }
    }

    rhs.level = std::numeric_limits<int>::max();
    if (node["generate_input"]) {
      if(node["generate_input"]["level"]){
        rhs.level = node["generate_input"]["level"].as<int>();
      }
      else if (node["generate_input"].size()==0 && node["generate_input"].as<std::string>()=="true"){
        rhs.level = 0;
      }
      else if (node["generate_input"].size()==0 && node["generate_input"].as<std::string>()=="false"){
        rhs.level = std::numeric_limits<int>::max();
      }
      else{
        throw ExceptionConfigFiles("job file: 'generate_input' only accepts 'true', 'false' or {level: <int>}");
      }
      if(rhs.level<0){
        throw ExceptionConfigFiles("job file: 'generate_input' only accepts positive {level: <int>}");
      }
      rhs.dataFile = "true";
    }


    rhs.config = 1;
    if (node["input_system"]) {

      if(node["input_system"].size()==0){
        rhs.inputSystem = node["input_system"].as<std::string>();
      }
      else if(node["input_system"]["files"]){

        std::vector<std::string> files;
        std::size_t numberOfEqs = std::numeric_limits<std::size_t>::max();
        int otf = 0;

        if(node["input_system"]["files"].size() == 0){

          files.push_back(node["input_system"]["files"].as<std::string>());
        }
        else if(node["input_system"]["files"].size() > 0){

          for (auto it: node["input_system"]["files"]) {

            files.push_back(it.as<std::string>());
          }
        }
        if(node["input_system"]["size"] && node["input_system"]["size"].size() == 0){
          numberOfEqs = node["input_system"]["size"].as<std::size_t>();
        }
        else if (node["input_system"]["size"] && node["input_system"]["size"].size() > 0){
          throw ExceptionConfigFiles("job file: 'input_system' only accepts the total number of equations\nIf multiple files or a directory containing multiple\nfiles is used, please enter the cumulative number of equations");
        }

        if(node["input_system"]["config"] && node["input_system"]["config"].size()==0){
          if(node["input_system"]["config"].as<std::string>() == "false"){
            rhs.config = 0;
          }
          else if(node["input_system"]["config"].as<std::string>() == "true"){
            rhs.config = 1;
          }
          else{
            throw ExceptionConfigFiles("job file: 'input_system': 'config' only accepts 'true' and 'false'");
          }
        }


        if(node["input_system"]["otf"] && node["input_system"]["otf"].size() == 0){
          if(node["input_system"]["otf"].as<std::string>() == "true"){
            otf = 1;
          }
          else if(node["input_system"]["otf"].as<std::string>() == "false"){
            otf = 0;
          }
          else{
            throw ExceptionConfigFiles("job file: 'input_system': 'otf' only accepts 'true' and 'false'");
          }
        }
        rhs.inputSystemTuple = make_tuple(files,numberOfEqs,otf);
      }
      else{
        throw ExceptionConfigFiles("job file: 'input_system': Unkown option");
      }
    }

    if (node["sector_selection"]) {
      for (YAML::const_iterator it = node["sector_selection"].begin();
           it != node["sector_selection"].end(); ++it) {
        std::string token = (*it).first.as<std::string>();
        if (token != "select_recursively") {
          throw ExceptionConfigFiles("job file: 'select_recursively': Option '" + token + "' unkown");
        }
      }

      if (node["sector_selection"]["select_recursively"]) {
        const Node& nodeb = node["sector_selection"]["select_recursively"];

        for (unsigned i = 0; i < nodeb.size(); i++) {
          rhs.sector2Reduce.push_back(make_pair(nodeb[i][0].as<std::string>(),
                                                nodeb[i][1].as<std::string>()));
        }
      }
    }

    if (node["identities"]) {
      for (YAML::const_iterator it = node["identities"].begin();
           it != node["identities"].end(); ++it) {
        std::string token = (*it).first.as<std::string>();
        if (token != "ibp") {
          throw ExceptionConfigFiles("job file: 'identities': Option '" + token + "' unkown");
        }
      }

      if (node["identities"]["ibp"]) {
        const Node& nodeibp = node["identities"]["ibp"];

        for (unsigned i = 0; i < nodeibp.size(); i++) {
          if (nodeibp[i]["r"] && nodeibp[i]["r"].size() == 2)
            rhs.den.push_back(make_pair(nodeibp[i]["r"][0].as<std::string>(),
                                        nodeibp[i]["r"][1].as<std::string>()));
          if (nodeibp[i]["s"] && nodeibp[i]["s"].size() == 2)
            rhs.num.push_back(make_pair(nodeibp[i]["s"][0].as<std::string>(),
                                        nodeibp[i]["s"][1].as<std::string>()));
        }
      }
    }

    if (node["select_masters"]) {
      const Node& nodeibp = node["select_masters"];

      rhs.masters = nodeibp.as<std::string>();
      if(!file_exists(rhs.masters.c_str())){
        throw ExceptionConfigFiles("job file: 'select_masters': File '" + rhs.masters + "' not found");
      };
    }

    if (node["preferred_masters"]) {
      const Node& nodeibp = node["preferred_masters"];
      rhs.masters = nodeibp.as<std::string>();
      if(!file_exists(rhs.masters.c_str())){
        throw ExceptionConfigFiles("job file: 'select_masters': File '" + rhs.masters + "' not found");
      };
    }

    if (node["extra_relations"]) {
      const Node& nodeibp = node["extra_relations"];
      rhs.extraRelations = nodeibp.as<std::string>();
      if(!file_exists(rhs.extraRelations.c_str())){
        throw ExceptionConfigFiles("job file: 'extra_relations': File '" + rhs.extraRelations + "' not found");
      };
    }

    if (node["check_masters"]) {
      if (!(node["select_masters"]) && !(node["preferred_masters"])) {
        throw ExceptionConfigFiles("job file: 'check_masters' only works together with 'select_masters' or 'preferred_masters'");
      }

      const Node& nodeibp = node["check_masters"];
      rhs.check_masters_ = nodeibp.as<bool>();
    }

    if (node["select_integrals"]) {
      for (YAML::const_iterator it = node["select_integrals"].begin();
           it != node["select_integrals"].end(); ++it) {
        std::string token = (*it).first.as<std::string>();
        if (token != "select_mandatory_list" &&
            token != "select_optional_list" &&
            token != "select_mandatory_recursively" &&
            token != "select_optional_recursively" &&
            token != "select_masters" &&
            token != "preferred_masters" &&
            token != "extra_relations" &&
            token != "check_masters") {
          throw ExceptionConfigFiles("job file: 'select_integrals': Option '" + token + "' unkown");
        }
      }

      if (node["select_integrals"]["select_masters"]) {
        const Node& nodeibp = node["select_integrals"]["select_masters"];

        rhs.masters = nodeibp.as<std::string>();

        if(!file_exists(rhs.masters.c_str())){
          throw ExceptionConfigFiles("job file: 'select_integrals': 'select_masters': File '" + rhs.masters + "' not found");
        };
      }

      if (node["select_integrals"]["extra_relations"]) {
        const Node& nodeibp = node["select_integrals"]["extra_relations"];

        rhs.extraRelations = nodeibp.as<std::string>();

        if(!file_exists(rhs.extraRelations.c_str())){
          throw ExceptionConfigFiles("job file: 'select_integrals': 'extra_relations': File '" + rhs.extraRelations + "' not found");
        };
      }

      if (node["select_integrals"]["preferred_masters"]) {
        const Node& nodeibp = node["select_integrals"]["preferred_masters"];

        rhs.masters = nodeibp.as<std::string>();

        if(!file_exists(rhs.masters.c_str())){
          throw ExceptionConfigFiles("job file: 'select_integrals': 'select_masters': File '" + rhs.masters + "' not found");
        };
      }

      if (node["select_integrals"]["select_mandatory_list"]) {
        const Node& nodeselect = node["select_integrals"]["select_mandatory_list"];

        for (unsigned i = 0; i < nodeselect.size(); i++) {
          if (nodeselect[i].size() == 1) {
            if (!file_exists(nodeselect[i][0].as<std::string>().c_str())) {
              throw ExceptionConfigFiles("job file: 'select_integrals': 'select_mandatory_list': File '" + nodeselect[i][0].as<std::string>() + "' not found");
            }
            rhs.mandatoryFileVector.push_back(
              nodeselect[i][0].as<std::string>());
          }
          else if (nodeselect[i].size() == 2) {
            if (!file_exists(nodeselect[i][1].as<std::string>().c_str())) {
              throw ExceptionConfigFiles("job file: 'select_integrals': 'select_mandatory_list': File '" + nodeselect[i][1].as<std::string>() + "' not found");
            }
            rhs.mandatoryFile.push_back(
            make_pair(nodeselect[i][0].as<std::string>(),
                      nodeselect[i][1].as<std::string>()));
          }
        }

        if(rhs.level != std::numeric_limits<int>::max()){
          throw ExceptionConfigFiles("job file: 'select_mandatory_list' not allowed for 'generate_input'\nUse 'select_mandatory_recursively' instead");
        }
      }

      if (node["select_integrals"]["select_mandatory_recursively"]) {
        const Node& nodeselect =
            node["select_integrals"]["select_mandatory_recursively"];

        for (unsigned i = 0; i < nodeselect.size(); i++) {
          if (nodeselect[i]["r"] && nodeselect[i]["s"]) {
            std::vector<std::string> topologies;
            std::vector<std::string> sectors;
            int rmax = -1;
            int smax = -1;
            int dmax = std::numeric_limits<int>::max();

            if (nodeselect[i]["topologies"]) {
              for (unsigned itT = 0; itT != nodeselect[i]["topologies"].size();
                   itT++) {
                topologies.push_back(
                    nodeselect[i]["topologies"][itT].as<std::string>());
              }
            }

            if (nodeselect[i]["sectors"]) {

              for (unsigned itT = 0; itT != nodeselect[i]["sectors"].size(); itT++){

                sectors.push_back(nodeselect[i]["sectors"][itT].as<std::string>());
              }
            }

            if (nodeselect[i]["r"] && nodeselect[i]["r"].size() == 0)
              rmax = nodeselect[i]["r"].as<int>();

            if (nodeselect[i]["s"] && nodeselect[i]["s"].size() == 0)
              smax = nodeselect[i]["s"].as<int>();

            if (nodeselect[i]["d"] && nodeselect[i]["d"].size() == 0)
              dmax = nodeselect[i]["d"].as<int>();

            if (rmax < 0) {
              throw ExceptionConfigFiles("job file: rmax in 'r: rmax' should be positive");
            }
            if (smax < 0) {
              throw ExceptionConfigFiles("job file: smax in 's: smax' should be positive");
            }
            if (dmax < 0) {
              throw ExceptionConfigFiles("job file: dmax in 'd: dmax' should be positive");
            }

            if (dmax == std::numeric_limits<int>::max())
              dmax = -1;

            rhs.selectSpec.push_back(
                std::make_tuple(topologies, sectors, rmax, smax, dmax));
          }
          else if (nodeselect[i].size() == 4) {
            std::vector<std::string> info;
            info.push_back(nodeselect[i][0].as<std::string>());
            info.push_back(nodeselect[i][1].as<std::string>());
            info.push_back(nodeselect[i][2].as<std::string>());
            info.push_back(nodeselect[i][3].as<std::string>());
            rhs.mandatoryRec.push_back(info);
          }
        }
      }

      if (node["select_integrals"]["select_optional_recursively"]) {
        const Node& nodeselect =
            node["select_integrals"]["select_optional_recursively"];
        for (unsigned i = 0; i < nodeselect.size(); i++) {
          std::vector<std::string> info;
          info.push_back(nodeselect[i][0].as<std::string>());
          info.push_back(nodeselect[i][1].as<std::string>());
          info.push_back(nodeselect[i][2].as<std::string>());
          info.push_back(nodeselect[i][3].as<std::string>());

          rhs.optionalRec.push_back(info);
        }
      }
    }

    if (node["create_denominator_database"]) {
      const Node& nodeTop = node["create_denominator_database"];

      for (unsigned iM = 0; iM < nodeTop.size(); iM++) {
        if (nodeTop[iM].size() == 0) {
          rhs.fileDenominators.push_back(nodeTop[iM].as<std::string>());
          continue;
        }
      }
    }

    if (node["amplitude_translate"]) {

      if(!(
        ((node["reduce"] && node["reduce"].size()>0)||(node["input_system"]))
          )
        ){
        throw ExceptionConfigFiles("job file: 'amplitude_translate' requires either 'reduce' or 'input_system' to be configured");
      }
      const Node& nodeTop = node["amplitude_translate"];

      for (unsigned iM = 0; iM < nodeTop.size(); iM++) {
        if (!file_exists(nodeTop[iM].as<std::string>().c_str())) {
          throw ExceptionConfigFiles("job file: 'amplitude_translate': File '" + nodeTop[iM].as<std::string>() + "' not found");
        }
        rhs.fileAmplitude.push_back(nodeTop[iM].as<std::string>());
      }
    }
    rhs.weightMode = 1;
    if(node["weight_notation"]){
      if(!(
        ((node["reduce"] && node["reduce"].size()>0)||(node["input_system"])))){
        throw ExceptionConfigFiles("job file: 'weight_notation' requires either 'reduce' or 'input_system' to be configured");
      }
      if(node["weight_notation"].as<std::string>()=="true"){
        rhs.weightMode = 1;
      }
      else if(node["weight_notation"].as<std::string>()=="false"){
        rhs.weightMode = 0;
      }
      else{
        throw ExceptionConfigFiles("job file: 'weight_notation': '" + node["weight_notation"].as<std::string>() + "' unkown");
      }

    }


    if (node["truncate_sp"]) {
      const Node& nodeTop = node["truncate_sp"];
      std::cout << nodeTop.size() << std::endl;
      for (unsigned i = 0; i < nodeTop.size(); i++) {
        std::cout << "l: " << std::endl;
        if (nodeTop[i]["l"]) {
          std::vector<std::string> topologies;
          std::vector<std::tuple<uint32_t,int>> tspexceptions;
          int splvl = 0;
          int truncateSP = -100;
          if (nodeTop[i]["topologies"]) {
            for (unsigned itT = 0; itT != nodeTop[i]["topologies"].size(); itT++) {
              topologies.push_back(
                  nodeTop[i]["topologies"][itT].as<std::string>());
            }
          }
          if (nodeTop[i]["tspexceptions"]) {
            auto iterE =nodeTop[i]["tspexceptions"];
            for (unsigned itT = 0; itT != iterE.size(); itT++){
              tspexceptions.push_back(std::make_tuple(iterE[itT][0].as<uint32_t>(),iterE[itT][1].as<int>()));
            }
          }
          // if (nodeTop[i]["sectors"]) {
          //   for (unsigned itT = 0; itT != nodeTop[i]["sectors"].size(); itT++){
          //     sectors.push_back(nodeTop[i]["sectors"][itT].as<std::string>());
          //   }
          // }
          // if (nodeTop[i]["sectors2"]) {
          //   for (unsigned itT = 0; itT != nodeTop[i]["sectors2"].size(); itT++){
          //     sectors2.push_back(nodeTop[i]["sectors2"][itT].as<std::string>());
          //   }
          // }
          if(nodeTop[i]["l"]){
            truncateSP = nodeTop[i]["l"].as<int>();
          }
          if(nodeTop[i]["lvl"]&&nodeTop[i]["lvl"].as<std::string>()=="true")
            splvl = 1;
          else if(nodeTop[i]["lvl"]&&nodeTop[i]["lvl"].as<std::string>()=="false")
            splvl = 0;
          else if(nodeTop[i]["lvl"]){
            throw ExceptionInternal("job file: not known option lvl: " + nodeTop[i]["lvl"].as<std::string>() + "\n");
          }
          rhs.truncateSP.push_back(std::make_tuple(topologies,tspexceptions,truncateSP,splvl));
        }
      }
    }

    rhs.trcMore = -1;
    if(node["truncate_more"]){
      rhs.trcMore = node["truncate_more"].as<int>();
    }
    rhs.trcMoreD = -1;
    if(node["truncate_more_dots"]){
      rhs.trcMoreD = node["truncate_more_dots"].as<int>();
    }



    if (node["insert_prefactors"]) {

      if(!(node["run_firefly"])){
        throw ExceptionConfigFiles("job file: 'insert_prefactors' only works with 'run_firefly'");
      }
      const Node& nodeTop = node["insert_prefactors"];

      for (unsigned iM = 0; iM < nodeTop.size(); iM++) {
        if (!file_exists(nodeTop[iM].as<std::string>().c_str())) {
          throw ExceptionConfigFiles("job file: 'insert_prefactors': File '" + nodeTop[iM].as<std::string>() + "' not found");
        }
        rhs.filePrefactors.push_back(nodeTop[iM].as<std::string>());
      }
    }

    if(node["iterative_reduction"]){

      rhs.iterativeReduction = node["iterative_reduction"].as<std::string>();

      if(rhs.iterativeReduction != "masterwise" && rhs.iterativeReduction != "sectorwise"){
        throw ExceptionConfigFiles("job file: 'iterative_reduction' only accepts 'masterwise' or 'sectorwise'");
      }

    }

    if (node["select_masters_reduction"]) {

      const Node& nodeTop = node["select_masters_reduction"];

      for (unsigned iM = 0; iM < nodeTop.size(); iM++) {
        if (nodeTop[iM].size() == 0) {
          if (!file_exists(nodeTop[iM].as<std::string>().c_str())) {
            throw ExceptionConfigFiles("job file: 'select_masters_reduction': File '" + nodeTop[iM].as<std::string>() + "' not found");
          }
          rhs.trimmedReduction.push_back(nodeTop[iM].as<std::string>());
          continue;
        }

        if (nodeTop[iM].size() != 2) {
          std::cout << "Wrong options in:\n";
          std::cout << "select_masters_reduction\n";
          continue;
        }
        std::vector<int> token;
        for (unsigned jM = 0; jM < nodeTop[iM][1].size(); jM++) {
          std::cout << nodeTop[iM][1][jM].as<int>() << std::endl;
          token.push_back(nodeTop[iM][1][jM].as<int>());
        }
        rhs.selectMastersReduction.insert(
            std::make_pair(nodeTop[iM][0].as<std::string>(), token));
      }
    }

    if (node["set_zero_sectors"]) {

      const Node& nodeTop = node["set_zero_sectors"];

      if(nodeTop["files"]){

        for (auto it: nodeTop["files"]) {
          if (!file_exists(it.as<std::string>().c_str())) {
            throw ExceptionConfigFiles("job file: 'set_zero_sectors': 'files': File '" + it.as<std::string>() + "' not found");
          }
          rhs.trimmedReduction.push_back(it.as<std::string>());
        }
      }
    }

    rhs.integralOrdering = 9;
    if (node["integral_ordering"]) {
      if (node["integral_ordering"].as<int>() < 1 ||
          node["integral_ordering"].as<int>() > 8) {
        throw ExceptionConfigFiles("job file: 'integral_ordering': Unkown ordering");
      }
      rhs.integralOrdering = node["integral_ordering"].as<int>();
    }

#ifdef KIRAFIREFLY
    if (node["run_firefly"])
      rhs.ff_recon = node["run_firefly"].as<std::string>();
    if (node["factor_scan"])
      rhs.factor_scan = node["factor_scan"].as<std::string>();
#endif

    return true;
  }
};

template <>
struct convert<Kira2File> {
  static bool decode(const Node& node, Kira2File& rhs) {
    for (YAML::const_iterator it = node.begin(); it != node.end(); it++) {
      std::string token = (*it).first.as<std::string>();
      if (token != "target" && token != "reconstruct_mass" &&
          token != "alt_dir") {
        throw ExceptionConfigFiles("job file: Option '" + token + "' unkown");
      }
    }

    if (node["target"]) {
      const Node& nodet = node["target"];

      for (unsigned i = 0; i < nodet.size(); i++) {

        if (nodet[i]["r"] && nodet[i]["s"]) {

          std::vector<std::string> topologies;
          std::vector<std::string> sectors;
          int rmax = -1;
          int smax = -1;
          int dmax = std::numeric_limits<int>::max();

          if (nodet[i]["topologies"]) {
            for (unsigned itT = 0; itT != nodet[i]["topologies"].size();
                 itT++) {
              topologies.push_back(
                  nodet[i]["topologies"][itT].as<std::string>());
            }
          }

          if (nodet[i]["sectors"]) {
            for (unsigned itT = 0; itT != nodet[i]["sectors"].size(); itT++){

              sectors.push_back(nodet[i]["sectors"][itT].as<std::string>());

            }
          }

          if (nodet[i]["r"] && nodet[i]["r"].size() == 0)
            rmax = nodet[i]["r"].as<int>();

          if (nodet[i]["s"] && nodet[i]["s"].size() == 0)
            smax = nodet[i]["s"].as<int>();

          if (nodet[i]["d"] && nodet[i]["d"].size() == 0)
            dmax = nodet[i]["d"].as<int>();

          if (rmax < 0) {
            throw ExceptionConfigFiles("job file: rmax in 'r: rmax' should be positive");
          }
          if (smax < 0) {
            throw ExceptionConfigFiles("job file: smax in 's: smax' should be positive");
          }
          if (dmax < 0) {
            throw ExceptionConfigFiles("job file: dmax in 'd: dmax' should be positive");
          }

          if (dmax == std::numeric_limits<int>::max()) dmax = -1;

          rhs.selectSpec.push_back(
              std::make_tuple(topologies, sectors, rmax, smax, dmax));
        }
        else if (nodet[i].size() == 2) {
          rhs.target.push_back(make_pair(nodet[i][0].as<std::string>(),
                                         nodet[i][1].as<std::string>()));
        }
        else if (nodet[i].size() == 4) {
          std::vector<std::string> info;
          info.push_back(nodet[i][0].as<std::string>());
          info.push_back(nodet[i][1].as<std::string>());
          info.push_back(nodet[i][2].as<std::string>());
          info.push_back(nodet[i][3].as<std::string>());
          rhs.mandatoryRec.push_back(info);
        }
        else if (nodet[i].size() == 1) {
          rhs.targetVector.push_back(nodet[i][0].as<std::string>());
        }
      }
    }

    if (node["alt_dir"]) rhs.inputDir = node["alt_dir"].as<std::string>();
    if (node["reconstruct_mass"])
      rhs.reconstructMass = node["reconstruct_mass"].as<std::string>();

    return true;
  }
};

template <>
struct convert<Merge> {
  static bool decode(const Node& node, Merge& rhs) {
    for (YAML::const_iterator it = node.begin(); it != node.end(); it++) {
      std::string token = (*it).first.as<std::string>();
      if (token != "files2merge" && token != "alt_dir") {
        throw ExceptionConfigFiles("job file: Option '" + token + "' unkown");
      }
    }

    if (node["files2merge"]) {
      const Node& nodet = node["files2merge"];

      for (unsigned i = 0; i < nodet.size(); i++) {
        if (!file_exists(nodet[i].as<std::string>().c_str())) {
          throw ExceptionConfigFiles("job file: 'files2merge': File '" + nodet[i].as<std::string>() + "' not found");
        }
        rhs.files2merge.push_back(nodet[i].as<std::string>());
      }
    }
    if (node["alt_dir"]) rhs.outputDir = node["alt_dir"].as<std::string>();

    return true;
  }
};

template <>
struct convert<Dgl> {
  static bool decode(const Node& node, Dgl& rhs) {
    for (YAML::const_iterator it = node.begin(); it != node.end(); it++) {
      std::string token = (*it).first.as<std::string>();
      if (token != "derive_dgl") {
        throw ExceptionConfigFiles("job file: Option '" + token + "' unkown");
      }
    }

    if (node["derive_dgl"]) {
      std::cout << "fff" << std::endl;
      const Node& nodet = node["derive_dgl"];

      for (unsigned i = 0; i < nodet.size(); i++) {
        if (!file_exists(nodet[i].as<std::string>().c_str())) {
          throw ExceptionConfigFiles("job file: 'derive_dgl': File '" + nodet[i].as<std::string>() + "' not found");
        }
        rhs.filesDGL.push_back(nodet[i].as<std::string>());
      }
    }

    return true;
  }
};

template <>
struct convert<Kinematics> {
  static bool decode(const Node& node, Kinematics& rhs) {
    for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
      std::string token = (*it).first.as<std::string>();

      if (token != "incoming_momenta" && token != "outgoing_momenta" &&
          token != "momentum_conservation" && token != "kinematic_invariants" &&
          token != "scalarproduct_rules"
          && token != "symbol_to_replace_by_one") {
        throw ExceptionConfigFiles("config/kinematics.yaml: Option '" + token + "' unkown");
      }
    }

    if (node["incoming_momenta"]) {
      const Node& nodei = node["incoming_momenta"];
      for (unsigned i = 0; i < nodei.size(); i++) {
        rhs.im.push_back(nodei[i].as<std::string>());
      }
    }
    if (node["outgoing_momenta"]) {
      const Node& nodeo = node["outgoing_momenta"];
      for (unsigned i = 0; i < nodeo.size(); i++) {
        rhs.om.push_back(nodeo[i].as<std::string>());
      }
    }

    if (node["momentum_conservation"]) {
      const Node& nodem = node["momentum_conservation"];
      if (nodem.size())
        rhs.mc =
            make_pair(nodem[0].as<std::string>(), nodem[1].as<std::string>());
    }

    if (node["kinematic_invariants"]) {
      const Node& nodek = node["kinematic_invariants"];
      for (unsigned i = 0; i < nodek.size(); i++) {
        rhs.ki.push_back(std::pair<std::string, int>(
            nodek[i][0].as<std::string>(), nodek[i][1].as<int>()));

        std::string invariantsToken = rhs.ki[i].first + "2place";

        rhs.invariantsPlaceholder.push_back(get_symbol(invariantsToken));

        rhs.invariantsReplacement.append(
            get_symbol(rhs.ki[i].first) ==
            get_symbol(invariantsToken)); // inv == inv'

        rhs.invariantsReplacementRev.append(
            get_symbol(invariantsToken) ==
            get_symbol(rhs.ki[i].first)); // inv' == inv
      }
    }
    if (node["scalarproduct_rules"]) {
      const Node& nodes = node["scalarproduct_rules"];
      for (unsigned i = 0; i < nodes.size(); i++) {
        std::pair<std::string, std::string> pure(
            nodes[i][0][0].as<std::string>(), nodes[i][0][1].as<std::string>());
        rhs.sr.push_back(
            std::pair<std::pair<std::string, std::string>, std::string>(
                pure, nodes[i][1].as<std::string>()));
      }
    }

    if (node["symbol_to_replace_by_one"]) {
      const Node noderpl = node["symbol_to_replace_by_one"];
      rhs.rpl = noderpl.as<std::string>();
    }

    return true;
  }
};

template <>
struct convert<Integral_F> {
  static bool decode(const Node& node, Integral_F& rhs) {
    for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
      std::string token = (*it).first.as<std::string>();
      if (token != "propagators" && token != "loop_momenta" &&
          token != "name" && token != "top_level_sectors" &&
          token != "magic_relations" && token != "cut_propagators" &&
          token != "permutation" && token != "permutation_option" &&
          token != "symbolic_ibp" &&
          token != "zero_sectors") {
        throw ExceptionConfigFiles("integralfamilies.yaml: Option '" + token + "' unkown");
      }
    }

    if (!node["propagators"]) {
      throw ExceptionConfigFiles("integralfamilies.yaml: Option 'propagators' missing");
    }
    if (!node["loop_momenta"]) {
      throw ExceptionConfigFiles("integralfamilies.yaml: Option 'loop_momenta' missing");
    }
    if (!node["name"]) {
      throw ExceptionConfigFiles("integralfamilies.yaml: Option 'name' missing");
    }

    if(node["name"].as<std::string>()=="Tuserweight"){
      throw ExceptionConfigFiles("integralfamilies.yaml: Topology name 'Tuserweight' is reserved for weight-bit notation");
    }

    rhs.name = node["name"].as<std::string>();


    if (node["magic_relations"])
      rhs.magic_relations = node["magic_relations"].as<std::string>();
    else
      rhs.magic_relations = "false";

    const Node& nodet = node["loop_momenta"];
    for (unsigned i = 0; i < nodet.size(); i++) {
      rhs.loop.push_back(nodet[i].as<std::string>());
      rhs.allVar.push_back(get_symbol(nodet[i].as<std::string>()));
      rhs.loopVar.push_back(get_symbol(nodet[i].as<std::string>()));
      rhs.loopVarList.append(get_symbol(nodet[i].as<std::string>()));

      std::string loop2 = rhs.loop[i] + "anything";
      rhs.loopVarList2.append(get_symbol(loop2));

      rhs.loop2loop2.append(rhs.loopVar[i] == get_symbol(loop2)); // k == k'
      rhs.loop2loop.append(get_symbol(loop2) == rhs.loopVar[i]);  // k' == k
    }

    rhs.permuteOption = "1";
    if (node["permutation"]) {
      rhs.permuteOption = "";
      const Node& nodeTop = node["permutation"];

      for (unsigned i = 0; i < nodeTop.size(); i++) {
        rhs.permute.push_back(nodeTop[i].as<int>());
      }
    }

    if (node["permutation_option"]) {
      const Node& nodeTop = node["permutation_option"];
      std::cout << "With topology: " << rhs.name << " use ";
      rhs.permuteOption = nodeTop.as<std::string>();
      if(rhs.permuteOption == "1"){
        std::cout << "permutation criteria option 1.\n";
        std::cout << "Prefer in ascending order:\n";
        std::cout << "1. massless propagators\n";
        std::cout << "2. shortest propagators\n";
      }
      else if(rhs.permuteOption == "2"){
        std::cout << "permutation criteria option 2.\n";
        std::cout << "Prefer in ascending order:\n";
        std::cout << "1. shortest propagators\n";
        std::cout << "2. massless propagators\n";
      }
      else if(rhs.permuteOption == "3"){
        std::cout << "permutation criteria option 3.\n";
        std::cout << "Prefer in ascending order:\n";
        std::cout << "1. longest propagators\n";
        std::cout << "2. massive propagators\n";
      }
      else if(rhs.permuteOption == "4"){
        std::cout << "permutation criteria option 4.\n";
        std::cout << "Prefer in ascending order:\n";
        std::cout << "1. massive propagators\n";
        std::cout << "2. shortest propagators\n";
      }
      else{
        throw ExceptionConfigFiles("integralfamilies.yaml: 'permutation_option: " + rhs.permuteOption + "' unkown");
      }
    }

    if(rhs.permuteOption.size()>0 && rhs.permute.size()>0){
      throw ExceptionConfigFiles("integralfamilies.yaml: 'permutation_option' conflicts with 'permutation'");
    }

    const Node& nodep = node["propagators"];

    for (unsigned i = 0; i < nodep.size(); i++) {
      unsigned entry = i;
//       if (rhs.permute.size() == nodep.size()) {
//         entry = rhs.permute[i] - 1;
//       }

      if (nodep[entry]["bilinear"]) {
        const Node& bilinear = nodep[entry]["bilinear"];
        rhs.propagator.push_back(std::pair<std::string, std::string>(
            "(" + bilinear[0][0].as<std::string>() + ")*(" +
                bilinear[0][1].as<std::string>() + ")",
            bilinear[1].as<std::string>()));
      }
      else {
        rhs.propagator.push_back(std::pair<std::string, std::string>(
            nodep[entry][0].as<std::string>(),
            nodep[entry][1].as<std::string>()));
      }
    }

    rhs.jule = static_cast<int>(nodep.size());

    int usedOption = 0;

    if (node["top_level_sectors"]) {
      const Node& nodeTop = node["top_level_sectors"];

      std::vector<uint32_t> tmpTopLevelSectors;
      for (size_t i = 0; i < nodeTop.size(); i++) {

        int sectorsTmp = pyred::Integral::parse_sector(nodeTop[i].as<std::string>(),rhs.jule);
        tmpTopLevelSectors.push_back(sectorsTmp);
      }

      std::sort(tmpTopLevelSectors.begin(), tmpTopLevelSectors.end(), std::greater<uint32_t>());

      std::vector<uint32_t> tmpLines;
      for (size_t i = 0; i < tmpTopLevelSectors.size(); i++) {

        uint32_t testSector = tmpTopLevelSectors[i];
        uint32_t num_ones = pyred::count_set_bits(testSector);

        int skip = 0;
        for(size_t ing = 0; ing < rhs.topLevelSectors.size(); ing++){
          if (num_ones < tmpLines[ing]){
            if (tmpTopLevelSectors[i] & rhs.topLevelSectors[ing]){
              skip = 1;
              continue;
            }
          }
          else{
            if(rhs.topLevelSectors[ing] == tmpTopLevelSectors[i]){
              skip = 1;
              continue;
            }
          }
        }

        if(skip==1)
          continue;

        tmpLines.push_back(num_ones);
        rhs.topLevelSectors.push_back(tmpTopLevelSectors[i]);
        rhs.bound.push_back(num_ones);
      }
      usedOption = 1;
    }

    if (node["zero_sectors"]) {
      const Node& nodeTop = node["zero_sectors"];

      for (unsigned i = 0; i < nodeTop.size(); i++) {

        int sectorsTmp = pyred::Integral::parse_sector(nodeTop[i].as<std::string>(),rhs.jule);

        rhs.userZeroSectors.push_back(sectorsTmp);
      }
    }

    rhs.maskCut = 0;

    if (node["cut_propagators"]) {
      const Node& nodeTop = node["cut_propagators"];

      for (unsigned i = 0; i < nodeTop.size(); i++) {
        rhs.cutProps.push_back(nodeTop[i].as<int>());
        rhs.maskCut += (1 << (nodeTop[i].as<int>() - 1));
      }
    }

    if (usedOption == 0) {
      rhs.topLevelSectors.push_back((1 << rhs.jule) - 1);
      rhs.bound.push_back(rhs.jule);
    }

    if (rhs.bound.size() > 0)
      rhs.biggestBound = rhs.bound[0];
    else
      rhs.biggestBound = 0;

    for (unsigned itB = 0; itB < rhs.bound.size(); itB++) {
      if (rhs.bound[itB] > rhs.biggestBound) rhs.biggestBound = rhs.bound[itB];
    }

    if (node["symbolic_ibp"]) {
      const Node& nodeS = node["symbolic_ibp"];
      for (unsigned i = 0; i < nodeS.size(); i++) {
        rhs.symbolicIBP.push_back((nodeS[i].as<int>()) - 1);
      }
    }
    return true;
  }
};

template <>
struct convert<IBPIntegral> {
  static bool decode(const Node& node, IBPIntegral& rhs) {
    const Node& noderef = node;

    rhs.id = pyred::Weight(noderef[0].as<std::string>());
    rhs.flag2 = noderef[1].as<int>();
    rhs.characteristics[SECTOR] = noderef[2].as<int>();
    rhs.characteristics[DOTS] = noderef[3].as<int>();
    rhs.characteristics[NUM] = noderef[4].as<int>();
    rhs.length = noderef[5].as<int>();
    for (unsigned i = 6; i < noderef.size() - 1; i++) {
      rhs.indices[i - 6] = noderef[i].as<int>();
    }
    rhs.coefficientString = noderef[noderef.size() - 1].as<std::string>();
    return true;
  }
};

template <>
struct convert<Auxiliary> {
  static bool decode(const Node& node, Auxiliary& rhs) {


    for (YAML::const_iterator it = node.begin(); it != node.end(); it++) {
      std::string token = (*it).first.as<std::string>();
      if (token != "name" && token != "permutation" &&
          token != "truncate_more" && token != "truncate_sp" &&
          token != "truncate_more_dots" && token != "user_numerics"
          && token != "user_prime" ) {
        throw ExceptionConfigFiles("job file: Option '" + token + "' unkown");
      }
    }

    if (!node["permutation"]) {
      throw ExceptionResume("sectormappings/auxiliary.yaml: 'permutation' missing");
    }
    if (!node["name"]) {
      throw ExceptionResume("sectormappings/auxiliary.yaml: 'name' missing");
    }

    if (node["user_numerics"]) {
      const Node& nodeTop = node["user_numerics"];
      for (unsigned i = 0; i < nodeTop.size(); i++) {
        rhs.userNumerics.push_back(make_pair(nodeTop[i][0].as<std::string>(),nodeTop[i][1].as<std::string>()));
      }
    }
    std::vector<uint32_t> permutation;
    if (node["permutation"]) {
      const Node& nodeTop = node["permutation"];
      for (unsigned i = 0; i < nodeTop.size(); i++) {
        permutation.push_back(nodeTop[i].as<uint32_t>());
      }
    }
    rhs.permutation=std::make_tuple(node["name"].as<std::string>(),permutation);

    int truncateSP = -100;
    int splvl = 0;
    std::vector<std::tuple<uint32_t,int>> tspexceptions;
    if (node["truncate_sp"]) {
      const Node& nodeTop = node["truncate_sp"];
      if (nodeTop["tspexceptions"]) {
        auto iterE =nodeTop["tspexceptions"];
        for (unsigned itT = 0; itT != iterE.size(); itT++){
          tspexceptions.push_back(std::make_tuple(iterE[itT][0].as<uint32_t>(),iterE[itT][1].as<int>()));
        }
      }
      if(nodeTop["l"]){
        truncateSP = nodeTop["l"].as<int>();
      }
      if(nodeTop["lvl"]&&nodeTop["lvl"].as<std::string>()=="true")
        splvl = 1;
      else if(nodeTop["lvl"]&&nodeTop["lvl"].as<std::string>()=="false")
        splvl = 0;
      else if(nodeTop["lvl"]){
        throw ExceptionInternal("job file: not known option lvl: " + nodeTop["lvl"].as<std::string>() + "\n");
      }
      rhs.truncateSP = std::make_tuple(tspexceptions,truncateSP,splvl);
      // if(nodeTop["sectors"]){
      //   for (auto it: nodeTop["sectors"]) {
      //     rhs.truncateSectors.push_back(it.as<unsigned int>());
      //   }
      // }
      // if(nodeTop["sectors2"]){
      //   for (auto it: nodeTop["sectors2"]) {
      //     rhs.truncateSectors2.push_back(it.as<unsigned int>());
      //   }
      // }
      // if(nodeTop["l"]){
      //   rhs.truncateSP = nodeTop["l"].as<int>();
      // }
      // if(nodeTop["lvl"]&&nodeTop["lvl"].as<std::string>()=="true")
      //   rhs.splvl = 1;
      // else if(nodeTop["lvl"]&&nodeTop["lvl"].as<std::string>()=="false")
      //   rhs.splvl = 0;
      // else if(nodeTop["lvl"]){
      //   throw ExceptionInternal("sectormappings/auxiliary.yaml: Not known option lvl: " + nodeTop["lvl"].as<std::string>() + "\n");
      // }
    }

    rhs.trcMore = -1;
    if(node["truncate_more"]){
      rhs.trcMore = node["truncate_more"].as<int>();
    }
    rhs.trcMoreD = -1;
    if(node["truncate_more_dots"]){
      rhs.trcMoreD = node["truncate_more_dots"].as<int>();
    }
    if(node["user_numerics"]){
      auto tmp_node = node["user_numerics"];
	    for(auto ix: tmp_node){
        rhs.randvars.push_back({ix[0].as<std::string>(), ix[1].as<std::string>()});
      }
    }
    if(node["user_prime"]){
      rhs.setPrimeNNumber = node["user_prime"].as<std::string>();
    }
    return true;
  }
};

} // namespace YAML

#endif
