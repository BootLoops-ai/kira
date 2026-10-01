/*
Copyright (C) 2017-2025 The Kira Developers (see AUTHORS file)
Modified 2026 for BootLoops by Matthew D. Schwartz; distributed under the GNU GPL v3 or later, see COPYING. The changes are described in PATCHES.md.

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

#ifndef KIRA_H_
#define KIRA_H_
#include <fstream>
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <regex>

#include "pyred/config.h"
#include "pyred/defs.h"
#include "pyred/interface.h"
#include "kira/ReadYamlFiles.h"
#include "kira/connect2kira.h"
#include "kira/dataBase.h"
#include "kira/integral.h"
#include "kira/threadpool.h"
#include "gzstream/gzstream.h"

#ifdef KIRAFIREFLY
#include <firefly/Reconstructor.hpp>

#include "kira/black_box.h"
#include "kira/black_box_numerical.h"
#endif

#define RED 20
#define SIZESYM 64

typedef std::vector<int> Vint;
typedef Vint::iterator ItVint;
typedef Vint::reverse_iterator RItVint;

typedef std::vector<std::vector<int>> VVint;
typedef VVint::iterator ItVVint;
typedef VVint::reverse_iterator RItVVint;

typedef std::vector<int*> Vintx;
typedef Vintx::iterator ItVintx;

typedef std::vector<BaseIntegral*> VG;
typedef VG::iterator ItVG;

typedef std::vector<BaseEquation*> VE;
typedef VE::iterator ItVE;

typedef std::vector<symmetries> SYM;
typedef SYM::iterator ItSYM;

typedef std::vector<std::vector<BaseIntegral*>> VVG;
typedef VVG::iterator ItVVG;

typedef std::unordered_map<pyred::Weight, std::tuple<int, int, int>> INTEGMAP;
typedef INTEGMAP::const_iterator INTEGMAPI;

typedef std::unordered_map<pyred::Weight, pyred::Weight> MASTERSMAP;
typedef MASTERSMAP::const_iterator MASTERSMAPI;

typedef std::map<pyred::Weight, std::tuple<BaseIntegral***, int>> FWRED;
typedef FWRED::iterator FWREDI;

typedef std::unordered_map<pyred::Weight, int> RDY2F;
typedef RDY2F::iterator RDY2FI;

typedef std::map<int, BaseIntegral*> map_id_to_key;
typedef map_id_to_key::iterator iterator_id_to_key;

typedef std::list<int> Lint;
typedef std::list<int>::iterator ItLint;

class ConvertResult;
class LeeCriterium;
class info;
class treatEq;
class Algebra;

struct tupelCharacteristics {
  int sector;
  int topology;
  int flag2;
};

class Kira {
public:
#ifndef KIRAFIREFLY
  Kira(std::string, int coreNumber /*, std::string algebra_*/,
       std::string pyred_config, int integralOrdering,
       std::vector<std::string> setVariable, std::vector<std::string> trim,
       std::string setSector, std::string setPrimeNNumber,
       std::string database_format, std::string pathNum, std::string auxName);
#else
  Kira(std::string, int coreNumber /*, std::string algebra_*/,
       std::string pyred_config, int integralOrdering,
       std::vector<std::string> setVariable, std::vector<std::string> trim,
       std::string setSector, std::string setPrimeNNumber,
       uint32_t bunch_size, int silent_flag_, std::string database_format,
       std::string pathNum, std::string auxName);
#endif

  friend std::ostream& operator<<(std::ostream& out, const BaseIntegral& per);
  friend std::istream& operator>>(std::istream& out, BaseIntegral& per);
  friend class ConvertResult;
  friend class BaseEquation;
  //   friend class SeedsKira;
  //   friend void initiate_seeds(Kira &kira, std::uint32_t nOfTopology,
  //   std::string InputMasters );
  friend void run_relations(Kira& kira_, std::uint32_t nOfTopology_);
  friend class Algebra;
  friend class Pak;
  friend class treatEq;
  // read base data and initiate these base data
  void read_config();
  void read_kinematics(int flag_user_defined_system);
  void init_kinematics();
  void read_integralfamilies(int flag_user_defined_system);
  void init_integralfamilies();
  void destroy_integralfamilies();
  void collect_reductions(Jobs& jobs);
  void collect_reductions_helper(
      std::string& itTopo,
      std::tuple<std::vector<std::string> /*topologies*/,
                 std::vector<std::string> /*sectors*/, int /*rmax*/, int /*smax*/,
                 int /*dmax*/>& itSpec);
  void get_topology_relations();
  void execute_jobs();

  static std::string treatcoeff2(const std::string& str);

  // create IBP and LI
  void create_IBP();
  void create_IBP_helper(std::vector<GiNaC::possymbol>& var);
  void create_LI();
  void reduce_scal(IBPVG& ibpp);

  void create_LEE();
  void create_LEE_vectors(std::vector<GiNaC::possymbol>& var);
  void create_LEE_vectors2(std::vector<GiNaC::possymbol>& var);

  // get zero sectors
  int find_zero_sectors();
  void update_permutations();

  // symmetry relations
  void search_symmetry_relations();
  int skip_symmetry(std::string topoName, int j, int it, int itt);
  int skip_symmetry_topology(unsigned op);
  int symmetry_finder(int it, int itt, std::string topoName, int j,
                      std::vector<int> array2, std::vector<int> array3, int flag, unsigned klop);
  int symmetry_finder_reverse(int it, int itt, std::string topoName, int j,
                              int array2[], int array3[], int flag);
  int prepare_symmetry();
  int symmetry_relations(unsigned op);
  int test(GiNaC::lst& matr,std::string& topoName);
  int testProps(GiNaC::lst& matr,std::string& topoName);
  void write_symmetries(const std::string otputName, SYM symVec[]);
  int read_symmetries(const std::string inputName, SYM symVec[]);
  void print_relations(const std::string name, SYM*& symVec, int itC);
  void print_relations_reverse(const std::string name, SYM*& symVec, int itC);

  void update_auxiliary_file(Jobs& jobs);

  std::vector<uint32_t> get_loop_details(Integral_F& topoX, int i);

  // setup a system of equations
  void read_integrals(
      std::string fileMasters,
      std::vector<std::tuple<std::string, std::vector<int>, std::string, int> >&
          arraySeed);
  void read_equations(
    std::string fileMasters,
    std::vector<std::vector<std::tuple<std::string, std::vector<int>, std::string, int> > >& vector_equations);
  void preferred_masters(std::string fileMasters);
  void print_equation(BaseIntegral*& integral,
                      GZSTREAM_NAMESPACE::ogzstream& output,
                      std::ofstream& outputX,
                      std::tuple<int, int /*,unsigned,int*/>& printControl,
                      std::tuple<std::string, std::string>& dest);
  void print_equationSW(BaseIntegral*& integral,
                        GZSTREAM_NAMESPACE::ogzstream& output,
                        std::ofstream& outputX,
                        std::tuple<int, int, unsigned, int>& printControl,
                        std::tuple<std::string, std::string>& dest,
                        int& rememberID);
  void print_equation2(BaseIntegral*& integral,
                       GZSTREAM_NAMESPACE::ogzstream& output);
  void write_init_system(std::vector<BaseEquation*>& setUpEq,
                         const std::string& dest, int secNUM,
                         const std::string& add);
  void write_pyred(
      std::vector<std::vector<std::pair<pyred::Weight, std::string>>>& eqs, /*
std::vector<std::vector<int > > & arraySeed,*/
      std::vector<std::size_t>& independent_eqnums);
  std::vector<std::vector<int>> read_id2int();

  void write_triangular(
      /*unsigned sectorNumber, int topology,*/ std::unordered_map<
          pyred::Weight, std::tuple<BaseIntegral***, int>>& forwardRed,
      std::vector<pyred::Weight>& reduce_please, int flagy);
  void write_triangularSW(
      unsigned sectorNumber, int topology,
      std::unordered_map<pyred::Weight, std::tuple<BaseIntegral***, int>>&
          forwardRed,
      std::vector<pyred::Weight>& reduce_please, int flagy);

  void skip_integral(GZSTREAM_NAMESPACE::igzstream& input, int eqLength);

  void trace_integral(GZSTREAM_NAMESPACE::igzstream& input, int eqLength);
  void read_integral(BaseIntegral*& integral,
                     GZSTREAM_NAMESPACE::igzstream& input, int eqlength);
  template <typename T>
  void read_integral2(
      T& input, int Eqlength,
      std::unordered_map<pyred::Weight, std::tuple<BaseIntegral***, int>>&
          forwardRed,
      std::unordered_map<pyred::Weight, int>& rdy2F);
  //   void
  //   load_pyred(std::vector<std::vector<std::pair<pyred::Weight,std::string> >
  //   > & eqs ); void load_pyred_on_the_fly(std::vector<pyred::Weight>&
  //   mandatory, std::vector<pyred::Weight>& optional);

#ifdef KIRAFIREFLY
  std::unordered_map<pyred::Weight,std::map<pyred::Weight,std::string> > prefactorEquations;
#endif
  void set_masters2zero();
  void set_masters2zero_sectorwise(std::vector<std::tuple<pyred::Weight, std::string, uint32_t> >& masterVector);
  void set_masters2zero_masterwise(std::vector<std::tuple<pyred::Weight, std::string, uint32_t> >& masterVector);
  std::vector<std::vector<pyred::Weight>> mastersSetZero;

  int load_triangular(
      unsigned k, int sectorNumber,
      std::unordered_map<pyred::Weight, std::tuple<BaseIntegral***, int>>&
          forwardRed,
      std::vector<pyred::Weight>& reduce_please,
      std::unordered_map<pyred::Weight, int>& rdy2F);
  int load_triangularVER(
      unsigned k, int sectorNumber,
      std::unordered_map<pyred::Weight, std::tuple<BaseIntegral***, int>>&
          forwardRed,
      std::vector<pyred::Weight>& reduce_please,
      std::unordered_map<pyred::Weight, int>& rdy2F);

  void check_reduced_integrals(uint32_t sectorNumber, uint32_t k,
                               DataBase databaseEQ[]);
  void load_all_integrals(uint32_t sectorNumber, uint32_t k,
                          DataBase databaseEQ[],
                          std::vector<std::pair<uint32_t, uint32_t>>& setofS,
                          std::uint32_t& reducedReihen,
                          std::uint32_t& skipReihen);
  int load_back_substitution(std::string& kiraSaveName);

  void initiate_system_of_equations(int nOfTopology);
  int setup_ibp_li_eq(std::vector<BaseEquation*>& setUpEq, VE& identities,
                      std::vector<pyred::Weight>& weights, std::string a,
                      int topology);
  int setup_sym_eq(std::vector<BaseEquation*>& setUpEq, std::string a,
                   int topologyNUM, int secNUM);
  int look_up_symmetry(std::vector<BaseEquation*>& setUpEq, IBPIntegral*& r,
                       SYM*& sym, int nOfTopology);

  // select equations
  void select_spec_helper(int itC, std::vector<pyred::SeedSpec>& initiateMAN,
                          std::tuple<std::vector<std::string> /*topologies*/,
                                     std::vector<std::string> /*sectors*/, int /*rmax*/,
                                     int /*smax*/, int /*dmax*/>& itSpec,
                          int& countCHOICE);
  void select_equations(std::vector<pyred::Weight>& mandatory,
                        std::vector<pyred::Weight>& optional, int itC,
                        Jobs& jobs, int flag_user_defined_system);
  void select_initial_integrals(std::vector<pyred::SeedSpec>& initiateSOE);
  void select_initial_integrals(std::vector<pyred::SeedSpec>& initiateSOE, int itLevel);
  void complement_initial_integrals(
      std::vector<pyred::SeedSpec>& complementSOE);

  // use pyred
  void extract_reduction_specs(std::vector<pyred::Weight>& mandatory);
  void write_seeds_to_disk(std::vector<pyred::Weight>& mandatory);
  int generate_SOE(std::vector<pyred::Weight>& selected_integrals, int generateFlag);
  void generate_input_system(std::vector<pyred::Weight>& selected_integrals, int recursivLevel, int weightFlag);
  void write_zero_equations(pyred::System& sys);
  int run_pyred(std::vector<pyred::Weight>& mandatory,
                 std::vector<pyred::Weight>& optional, pyred::System& sys, int generateFlag);
  int run_pyred_otf(std::vector<pyred::Weight>& mandatory,
                 std::vector<pyred::Weight>& optional, pyred::System& sys, std::vector<std::string> & files, int generateFlag);
  //   void get_reduced_eqs(pyred::System & sys);
  void record_masters(std::vector<pyred::Weight>& masters, std::vector<pyred::Weight>& mandatory);
  void CheckMasters(std::vector<pyred::Weight>& masters);

  // run reduction
  std::vector<std::string> collectReductions;
  std::vector<std::string> topologyNames;
  void initiate_fermat(int kira2math, int onlyOnePool, int noFermat);
  void destroy_fermat(int onlyOnePool, int noFermat);
  void complete_triangular(std::vector<pyred::Weight>& mandatory);
  void complete_triangularSW(std::vector<pyred::Weight>& mandatory);
  int insert_equation(BaseIntegral*& EqPtr, BaseIntegral*& insertEqPtr);
  //   void run_reduction( std::map<pyred::Weight, std::tuple<BaseIntegral ***,
  //   int> >& forwardRed,
  // 		      std::vector <pyred::Weight>& reduce_please,
  // 		      std::unordered_map<pyred::Weight,int>& rdy2F,
  // std::set<pyred::Weight>& setMandatory);

  void run_reduction(
      std::unordered_map<pyred::Weight, std::tuple<BaseIntegral***, int>>&
          forwardRed,
      std::vector<pyred::Weight>& reduce_please,
      std::unordered_map<pyred::Weight, int>& rdy2F,
      std::set<pyred::Weight>& setMandatory, std::tuple<int, int>& printControl);

  void parallel_pool(
      std::unordered_map<pyred::Weight, std::tuple<BaseIntegral***, int>>&
          forwardRed,
      std::vector<pyred::Weight>& reduce_please_back_up, int& red_1,
      std::unordered_map<pyred::Weight, int>& rdy2F,
      std::set<pyred::Weight>& setMandatory, pyred::Weight it,
      std::uint32_t& information);

  std::pair<int, int> get_power(std::string& stringofinterest);

  void init_numbers(std::vector<int>& numbersINIT,
                    std::vector<std::vector<std::string>>& numbersDenStrVec);

  int get_power_level(std::string& numeratorPart);

  std::pair<std::string, std::string> normsample(std::string stringofinterest,
                                                 int& valueForD, int& degree,
                                                 int& powerLevel1,
                                                 int& powerLevel2);

  std::vector<int> generate_sample_points(
      std::pair<std::string, std::string>& normsampleVar, int degree);

  std::vector<int> check_sample_points(
      std::pair<std::string, std::string>& normsampleVar,
      std::vector<int>& points);

  std::tuple<std::vector<std::string>, std::vector<std::string>,
             std::vector<int>>
  sample(std::vector<int>& numbersINIT, std::string& stringofinterest,
         std::pair<std::string, std::string>& normsampleVar, pyred::Weight& tmp1,
         pyred::Weight& tmp2);

  std::tuple<std::vector<std::string>, std::vector<std::string>,
             std::vector<int>>
  sample(std::vector<int>& numbersINIT, std::list<std::string>& ptrToString,
         std::pair<std::string, std::string>& normsampleVar, uint64_t& tmp1,
         uint64_t& tmp2);

  std::vector<std::string> normalize(std::vector<int>& numbersINIT,
                                     int valueForD,
                                     std::vector<std::string>& numerator,
                                     std::string& numNormalization);

  std::pair<std::string, int> reconstruct_function(
      std::vector<int>& numbersINIT,
      std::vector<std::vector<std::string>>& numbersDenStrVec,
      std::vector<std::string>& numerator, int& powerLevel);

  std::string reconstruct_final(std::string& reconstrNum,
                                std::string& reconstrDen);
  std::string reconstruct_final(std::list<std::string>& ptrToString);
  void sort_equation(BaseIntegral*& EqPtr);

  void check_fermat_firefly_used(int flag);
  int complete_reduction();

  int get_number_of_masterintegrals(unsigned it,
                                    std::vector<pyred::Weight>& lriMI,
                                    int** IDarray);
  unsigned addup_common_masters(
      unsigned it, unsigned j, int** IDarray, BaseIntegral newEqPtr[],
      std::vector<std::list<std::string>>& listPtrToString);
  unsigned addup_final_masters(
      unsigned it, unsigned j, int** IDarray, BaseIntegral newEqPtr[],
      std::vector<std::list<std::string>>& listPtrToString);
  void back_subs(std::string& kiraSaveName);
  void write_result(BaseIntegral*& EqPtr, int i_Integrals,
                    std::ofstream& output);
  void clean_back_subs();

  // output mathematica
  void kira_mathematica(std::string, int massRC);

#ifdef KIRAFIREFLY
  void run_firefly(/*const*/ std::vector<pyred::Weight>& mandatory_vec,
      const uint32_t mode_, const std::string topology_name,
      const int factor_scan);
  template <typename BlackBoxType>
  void load_ff_system(BlackBoxType* bb_, const std::vector<std::string>& variables);
  void prepare_factors(const std::vector<pyred::Weight>& mandatory_vec, const std::unordered_set<pyred::Weight>& masters_set_to_zero);
  void write_to_database(std::unordered_set<pyred::Weight>& masters);
  void run_numerical(const std::vector<pyred::Weight>& mandatory_vec, const std::string& topology_name, const std::string& input_file_name);
  void numerical_compute(std::queue<std::pair<std::vector<firefly::FFInt>, pyred::SystemOfEqs<firefly::FFInt>>>& results, const std::tuple<std::string, std::vector<std::string>, std::vector<std::vector<std::string>>>& set_of_points, std::size_t pos);
  template <std::size_t N>
  void numerical_compute_bunch(std::queue<std::pair<std::vector<firefly::FFInt>, pyred::SystemOfEqs<firefly::FFInt>>>& results, const std::tuple<std::string, std::vector<std::string>, std::vector<std::vector<std::string>>>& set_of_points, std::size_t pos);
  void numerical_write_results(std::queue<std::pair<std::vector<firefly::FFInt>, pyred::SystemOfEqs<firefly::FFInt>>>& results, const std::vector<pyred::Weight>& mandatory_vec, const std::string& topology_name, const std::string& input_file_name, const std::size_t points);
  void numerical_write_equation(const std::vector<pyred::Equation<firefly::FFInt>>::iterator it, std::ofstream& filestream);
#endif

  void write_amplitude_file(pyred::Weight generateID, std::vector<std::tuple<std::string, std::vector<int>, std::string, int> >&, std::string& amplitudeFile, int weightFlag);

public:
  static int s_coreNumber; // parallel-treatcoeff: pool size for per-thread fermat[tid]
private:
  std::string kiraMode;
  static Fermat* fermat;
  ThreadPool *pool, *pool2, *pool3, *pool4, *pool5, *poolBS;
  pyred::Weight core_mask;
  std::stack<int> idFermats, idCombine;
  int coreNumber /*, algebra*/, integralOrdering, sectorOrdering;
  std::string pyred_config;
  std::string outputDir;
  std::vector<std::tuple<std::string /*prime*/, std::vector<std::string> /*symbols*/, std::vector<std::vector<std::string>> /*list of sample points*/>> numericalPointsAux;
  std::string inputDir;
  std::vector<std::tuple<std::string, std::string>> setVariable;
  std::vector<std::tuple<std::string, uint32_t>> systemTrim;
  std::string setSector;
  std::string setPrimeNNumber;
  std::string databaseFormat;

  int dataFile, pyredDatabase, magicRelations, conditionalSystem, conditionalMemory,
  getOnlyMasters;
  unsigned heuristic, reconstFlag, termNumber, algebraicReconstruction;
  bool LIflag, symbols2numsFlag;
  std::size_t pointsPerFile = 1000;
  std::string pathNum, auxName;

  int calculate_coefficient(std::list<std::string>& ptrToString);
  int calculate_coefficient_term(std::list<std::string>& ptrToString,
                                 int idComb);
  int calculate_coefficient_term2(std::list<std::string> ptrToString,
                                  std::string& interpVar, int numbers,
                                  std::string& result, int idComb);
  std::string jobName, specialName;

  // maps
  unsigned numberOfEq, numberOfEqXXX;
  BaseIntegral** systemEq;
  std::unordered_map<pyred::Weight, pyred::Weight> mastersMap, mastersReMap;
  //   std::vector<BaseEquation*> setUpEq;
  //   std::vector<BaseEquation> setUpEq2;
  //   BaseEquation** setUpEq;

  // DGL
  void generate_dgl();
  void dgl_vs_scalar_products();
  std::vector<std::tuple<std::string, int, std::vector<int>, std::string>>
  read_seeds_dgl(std::string& itFile);
  void insert_seeds2DGL(
      std::vector<std::tuple<std::string, int, std::vector<int>, std::string>>&
          seedsDGL);

  // zero sectors/symmetries/shifts
  std::pair<int, GiNaC::ex> test_quadratic(GiNaC::ex& start);
  void trim_the_system();

  GiNaC::lst suby;
  std::vector<std::tuple<GiNaC::lst, GiNaC::lst, int, GiNaC::ex,
                         std::vector<std::string>>>
      externalTransf;
  GiNaC::lst invariants4sym, invariants4symRev;
  GiNaC::lst invariantsReplacement, invariantsReplacementRev, invariantsList;
  std::vector<GiNaC::possymbol> invariantsPlaceholder, symbolInvariants;
  int controlSymmetries;
  std::vector<std::string> invarMap;

  // initiate kira, base data
  Integral_F integralfamily;
  std::map<std::string, Integral_F> topology;
  std::string fermatPath;
  int numFlag, biggestBound, realSector;
  std::vector<std::tuple<std::string, std::vector<int>, std::string, int>>
      preferredMasterSectors;
  bool check_masters_ = false;
  bool autoMode = false;
  std::vector<int> num, numMin, den, denTPlus;
  std::vector<GiNaC::possymbol> externalVar, allVar, invar;

  std::vector<std::string> reductVar;
  std::string interpVar;
  std::vector<int> invarDim;
  std::vector<std::string> invarStr, invarSol;
  int massSet2OneDim;
  GiNaC::ex dimension, massSet2One;
  GiNaC::lst kinematic, kinematicR, kinematic2, kinematicOld, specialKinematics,
      momentConservation, kinematicShift, kinematicShiftB, kinematicReverse,
      mass2One, kinematicShiftR, unknownsExt;
  std::vector<GiNaC::possymbol> bSsymbols;

  GiNaC::ex mom_uno;
  GiNaC::possymbol* bS;
  GiNaC::symtab GiNaCSymbols;

#ifdef KIRAFIREFLY
  uint32_t max_bunch_size = 1;
  int silent_flag = 0;
  uint32_t reconstruction_mode =
      0; // 0: full reconstruction, 1: only back substitution
  bool multiply_factors = false;
  std::vector<std::string> symbols;

  BlackBoxKira* bb;
  firefly::Reconstructor<BlackBoxKira>* reconst;

  BlackBoxNumerical* bbNum;

  DataBase* database_reconst;
  bool first = true;
  std::unordered_map<pyred::Weight, std::pair<int, int>>
      integral_data; // sector, flag2
#endif

  int onlyBacksubstitution, onlyTriangular, onlyPyred, onlyInitiate, startJob;

  unsigned reihen, totalReihen;
  BaseIntegral*** allEq;
  std::vector<pyred::Weight> masterVectorSkip;
  std::multimap<std::string, std::vector<int>> selectMastersReduction;
  std::string iterativeReduction;
  std::vector<std::string> trimmedReduction;
  std::vector<std::tuple<uint32_t, std::string>> trimmedSectors;

  unsigned *length, *rdy2P;
  std::unordered_map<pyred::Weight, unsigned> occurrence;
  std::vector<unsigned> reduct2StartHere;
  pyred::Weight* last_reduce;
  std::unordered_map<pyred::Weight, unsigned> reverseLastReduce;
  unsigned eqnum;
};

class info {
public:
  int dots;
  int nums;
  int sector;
  int topology;
};

class kiraOutput {
public:
  std::string str[14];
  kiraOutput(int choice) {
    switch (choice) {
      case 0:
        str[0] = "FORM";
        str[1] = ".inc";
        str[2] = "";
        str[3] = "id ";
        str[4] = " + ";
        str[5] = "(";
        str[6] = ";\n\n";
        str[7] = ") = \n";
        str[8] = ")*(";
        str[9] = ")\n";
        str[10] = "0\n";
        str[11] = ";";
        break;
      case 1:
        str[0] = "Mathematica";
        str[1] = ".m";
        str[2] = "{\n";
        str[3] = "";
        str[4] = " + ";
        str[5] = "[";
        str[6] = ",\n";
        str[7] = "] -> \n";
        str[8] = "]*(";
        str[9] = ")\n";
        str[10] = "0\n";
        str[11] = "}";
        break;
      case 2:
        str[0] = "Kira";
        str[1] = ".kira";
        str[2] = "";
        str[3] = "";
        str[4] = "";
        str[5] = "[";
        str[6] = "\n";
        str[7] = "]*(-1)\n";
        str[8] = "]*(";
        str[9] = ")\n";
        str[10] = "";
        str[11] = "";
        str[12] = "*(-1)\n";
        str[13] = "*(";
        break;
      case 3:
        str[0] = "FORM";
        str[1] = ".h";
        str[2] = "";
        str[3] = "Fill ";
        str[4] = " + ";
        str[5] = "(";
        str[6] = ";\n\n";
        str[7] = ") = \n";
        str[8] = ")*(";
        str[9] = ")\n";
        str[10] = "0\n";
        str[11] = ";";
        break;
    }
  };
};

class Algebra {
public:
  Algebra(){};
  Algebra(Kira* kira_, std::string& interpVar_,
          std::vector<std::string>& reductVar_,
          std::vector<std::pair<std::string, int>>& numericVar_, int level);

  Kira* kira;
  std::string interpVar;
  std::vector<std::string> reductVar;
  std::vector<std::pair<std::string, int>> numericVar;
  int depthLevel;
  int poolLevel;

  std::string reconstruct_final(std::string& stringofinterest1,
                                std::string& stringofinterest2);

  std::pair<std::string, std::string> normsample(std::string stringofinterest,
                                                 int& valueForD, int& degree,
                                                 int& powerLevel1,
                                                 int& powerLevel2);

  std::tuple<std::vector<std::string>, std::vector<std::string>,
             std::vector<int>>
  sample(std::vector<int>& numbersINIT, std::string& stringofinterest,
         std::pair<std::string, std::string>& normsampleVar, uint64_t& tmp1,
         uint64_t& tmp2);

  std::vector<std::string> normalize(std::vector<int>& numbersINIT,
                                     int valueForD,
                                     std::vector<std::string>& numerator,
                                     std::string& numNormalization);

  std::pair<std::string, int> reconstruct_function(
      std::vector<int>& numbersINIT,
      std::vector<std::vector<std::string>>& numbersDenStrVec,
      std::vector<std::string>& numerator, int& powerLevel);

  std::vector<int> check_sample_points(
      std::pair<std::string, std::string>& normsampleVar,
      std::vector<int>& points);

  std::vector<int> generate_sample_points(
      std::pair<std::string, std::string>& normsampleVar, int degree);

  int get_power_level(std::string& numeratorPart);

  void init_numbers(std::vector<int>& numbersINIT,
                    std::vector<std::vector<std::string>>& numbersDenStrVec);

  std::pair<int, int> get_power(std::string& stringofinterest);
};

class ConvertResult {
public:
  ConvertResult(){};
  ConvertResult(Kira& kira, std::string topologyName, int topologyNumber,
                std::string& inputName, const std::string& inputDir,
                std::vector<pyred::Weight>& idOfSeed);
  ConvertResult(Kira& kira, std::string topologyName, int topologyNumber,
                const std::string& inputDir,
                std::vector<pyred::Integral>& listOfIntegrals,
                std::vector<pyred::Weight>& idOfSeed, std::string& inputName);
  ~ConvertResult();
  void look_up_seeds();
  void get_integral_id(std::vector<pyred::Weight>& idOfSeed);
  void get_integral_id_pyred(std::vector<pyred::Weight>& idOfSeed);
  int skip_integral(std::vector<int>& seed);
  void reconstruct_init(Kira& kira);
  void reconstruct_mass(Kira& kira, std::vector<DBintegral>& integralV);
  void prepare_FORM(std::vector<DBintegral>& integralV);
  void get_firefly_flag();

  int output(Kira& kira, int massReconstruction,
             std::vector<pyred::Weight>& idOfSeed, int choice);
  DataBase *database;
  int fireflyFlag;
  std::string topologyName;
  int topologyNumber;
  std::string inputDir, inputName;
  std::vector<std::regex> massReconstructRegex;
  std::vector<std::string> replaceWord;
  std::regex word_regex{"(.*)\\*prefactor\\[(.*)\\]"}, word_regex2{"(.*)\\/\\((.*)"}, word_regex3{"(.*)\\/(?![(])(.*)"};

  std::ofstream Output;
  std::ifstream seedsInput;
  std::ifstream trivialSectorInput;
  std::vector<int> trivSV;
  std::vector<std::vector<int>> arraySeed;
  std::vector<std::tuple<std::vector<int>,uint32_t> > zeroIntegrals;
  std::vector<GiNaC::possymbol> invarT;
};

#endif // KIRA_H_
