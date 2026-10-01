/*
Copyright (C) 2017-2025 The Kira Developers (see AUTHORS file)

This file is part of pyRed.

pyRed is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyRed is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyRed.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PYRED_INTERFACE_H
#define PYRED_INTERFACE_H

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "pyred/coeff_helper.h"
#include "pyred/coeff_int.h"
#include "pyred/coeff_vec.h"
#include "pyred/defs.h"
#include "pyred/gauss.h"
#include "pyred/integrals.h"
#include "pyred/keyvaluedb.h"
#include "pyred/parallel.h"
#include "pyred/ppmacros.h"
#include "pyred/relations.h"

namespace pyred {

// content[integral] = pair(eqnum, {unreduced_integrals})
using content_type =
    std::unordered_map<Weight, std::pair<std::size_t, std::vector<Weight>>>;

// content[integral] = pair(eqnum, {<unreduced_integrals,coeff/FF>})
using content_coeff_type =
    std::unordered_map<Weight, std::pair<std::size_t, std::vector<std::tuple<Weight,std::string>>>>;

float time_diff(const std::clock_t &, const std::clock_t &);
float time_diff(const std::chrono::time_point<std::chrono::system_clock> &start,
                const std::chrono::time_point<std::chrono::system_clock> &end);

class System {
public:

  class EquationGetter {
    // A callable object to retrieve equations via System::generate_retrieve().
    // Collect equations as {(weight,string), ...}.
  private:
    // Using a (copyable) pointer here is essential, because when an instance
    // of EquationGetter is bound to a std::function object, the EquationGetter
    // is copied. I.e. if we had a std::vector<eqdata> member instead,
    // the equations would be inserted into a copy of the system which is not
    // accessible through the original EquationGetter instance.
    std::shared_ptr<std::vector<eqdata>> m_sys_ptr;
  public:
    EquationGetter() : m_sys_ptr{std::make_shared<std::vector<eqdata>>()} {}
    std::vector<eqdata> &get() { return *m_sys_ptr; }
    void operator()(eqdata &&eq) { m_sys_ptr->push_back(std::move(eq)); }
  };

  System();
  System(const std::vector<eqdata> &);
  System(std::vector<eqdata> &&);
  System(const std::string &);
  System(const std::vector<std::string> &);
  std::size_t size() const;
  static void sorteq(eqdata &);
  std::vector<std::size_t> generate_solve(const std::vector<SeedSpec> &,
                                          const std::vector<SeedSpec> &,
                                          const std::vector<SeedSpec> &);
  void reserve(const std::size_t);
  void add(const eqdata &);
  void add(eqdata &&);
  void add(const std::vector<eqdata> &);
  void add(std::vector<eqdata> &&);
  void add(const std::string &, const std::vector<std::string> &fileexts = {},
           bool unsafe = false);
  void add(const std::vector<std::string> &,
           const std::vector<std::string> &fileexts = {},
           bool unsafe = false);
#define PYRED_PP_CCS_SYSTEMADD(k)                                   \
  void add(wc_equation<PYRED_PP_COEFFCLASS(k)> &&eq, std::size_t eqnum) { \
    m_numsys.PYRED_PP_COEFFCLASSMEM(k).add(std::move(eq), eqnum);   \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMADD)
  void add_forward(const eqdata &,
                   std::size_t neq = std::numeric_limits<std::size_t>::max());
  void add_forward(const std::string &infile,
                   const std::vector<std::string> &fileexts = {},
                   bool unsafe = false);
  void add_forward(const std::vector<std::string> &infiles,
                   const std::vector<std::string> &fileexts = {},
                   bool unsafe = false);
  std::vector<std::size_t> solve();
  void backward();

  void backward_coeff();
  std::vector<std::size_t> independent();
  std::vector<Weight> unreduced() const;
  std::pair<std::vector<std::size_t>, std::vector<Weight>> select(
      const std::vector<Weight> &, const std::vector<Weight> & = {});
  std::pair<std::vector<std::size_t>, std::vector<Weight>> select_coeff(
      const std::vector<Weight> &, const std::vector<Weight> & = {});

  static std::pair<std::size_t,std::vector<uint32_t>>
  file_counteqs_setweightbits(const std::vector<std::string> &infiles,
                              const std::vector<std::string> &fileexts,
                              bool forcecount = false);
  std::vector<uint32_t>
  reserve_setweightbits(const std::vector<std::string> &infiles,
                        const std::vector<std::string> &fileexts,
                        bool forcecount = false);
  std::vector<eqdata> retrieve(std::vector<std::size_t> &&eqnums);
  static std::size_t file_retrieve(
    const std::string &infile, const std::vector<std::string> &fileexts,
    std::vector<std::size_t> &&eqnums, const std::function<void(eqdata &&)> &treateq,
    const std::function<std::string(std::string &)> &treatcoeff = nullptr,
    bool unsafe = false);
  static std::size_t file_retrieve(
    const std::vector<std::string> &infiles,
    const std::vector<std::string> &fileexts, std::vector<std::size_t> &&eqnums,
    const std::function<void(eqdata &&)> &treateq,
    const std::function<std::string(std::string &)> &treatcoeff = nullptr,
    bool unsafe = false);
  static std::size_t file_retrieve_sorted(
    const std::string &infile, const std::vector<std::string> &fileexts,
    std::vector<std::size_t> &&eqnums, const std::function<void(eqdata &&)> &treateq,
    const std::function<std::string(std::string &)> &treatcoeff = nullptr,
    bool unsafe = false);
  static std::size_t file_retrieve_sorted(
    const std::vector<std::string> &infiles,
    const std::vector<std::string> &fileexts, std::vector<std::size_t> &&eqnums,
    const std::function<void(eqdata &&)> &treateq,
    const std::function<std::string(std::string &)> &treatcoeff = nullptr,
    bool unsafe = false);
  void generate_retrieve(
      const std::vector<SeedSpec> &ibp_seedspec,
      const std::vector<SeedSpec> &ibp_seedcompl,
      const std::vector<SeedSpec> &sym_seedspec,
      std::vector<std::size_t>&& eqnums,
      const std::function<void(eqdata &&)> &treateq,
      const std::function<std::string(const std::string &)> &/*treatcoeff*/
        = nullptr);
  void generate_retrieve(
      std::vector<std::size_t> &&eqnums,
      const std::function<void(eqdata &&)> &treateq,
      const std::function<std::string(const std::string &)> &/*treatcoeff*/
        = nullptr);
  std::vector<Weight> reduction_content(Weight i) const {
    return m_content.at(i).second;
  }
  const content_type &reduction_content() const { return m_content; }
  std::vector<std::tuple<Weight,std::string>> reduction_content_coeff(Weight i) const {
    return m_content_coeff.at(i).second;
  }
  const content_coeff_type &reduction_content_coeff() const {return m_content_coeff; }
  template <typename Coeff>
  const SystemOfEqs<Coeff> &get_numsys() const;
  const std::vector<eqdata> &get_sys() const { return sys; }

private:
  std::vector<eqdata> sys;
  int coeff_cls;
  // True only if the system was generated with generate(.,.),
  // i.e. seeds is set and generate() may be called without arguments.
  bool m_is_generated{false};
  // ibp and symmetry seeds. Set when generate(ibp_seeds, sym_seeds)
  // is called; used when generate() is called without arguments.
  std::tuple<std::vector<SeedSpec>, std::vector<SeedSpec>,
             std::vector<SeedSpec>>
      m_seeds;
  // unordered_map content[integral] = pair(eqnum, {unreduced_integrals})
  bool m_content_prepared;
  // Was numsys sorted by the solver before solving?
  // This happens if the solve method is used instead of the on-the-fly solver.
  // If yes, the selected equations must be sorted in the same way.
  bool m_solver_sorted;
  content_type m_content;
  content_coeff_type m_content_coeff;
  // Data members for the on-the-fly forward solver and
  // to directly add equations with finite integer coefficients.
  // SystemOfEqs cannont be packed in a union.
  // Maybe C++ 17 std::variant would be useful here.
  std::size_t neqs = 0u;
  struct {
#define PYRED_PP_CCS_NUMSYSDECL(k) \
  SystemOfEqs<PYRED_PP_COEFFCLASS(k)> PYRED_PP_COEFFCLASSMEM(k);
    PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_NUMSYSDECL)
  } m_numsys;
  struct {
#define PYRED_PP_CCS_SOLMAPDECL(k) \
  sol_map<PYRED_PP_COEFFCLASS(k)> PYRED_PP_COEFFCLASSMEM(k);
    PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SOLMAPDECL)
  } m_sols;
  struct {
#define PYRED_PP_COEFFMAPDECL(k)                          \
  std::unordered_map<PYRED_PP_COEFFCLASS(k), std::string> \
      PYRED_PP_COEFFCLASSMEM(k);
    PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_COEFFMAPDECL)
  } m_coeff_map;
  // methods
  int systemtype(bool = false) const;
  void add_file(const std::string &,
                const std::vector<std::string> &fileexts = {},
                bool unsafe = false);
  template <typename Coeff>
  SystemOfEqs<Coeff> &get_numsys();
  template <typename Coeff>
  std::unordered_map<Coeff, std::string> &get_coeff_map();
  template <typename Coeff>
  std::vector<std::size_t> generate_solve_tmpl(
    const std::vector<SeedSpec> &,
    const std::vector<SeedSpec> &,
    const std::vector<SeedSpec> &);
  template <typename Coeff>
  void add_forward_tmpl(const eqdata &, SystemOfEqs<Coeff> &, sol_map<Coeff> &,
                        std::size_t neq);
  template <typename Coeff>
  void backward_tmpl(SystemOfEqs<Coeff> &);
  template <typename Coeff>
  void backward_coeff_tmpl(SystemOfEqs<Coeff> &);
  template <typename Coeff>
  std::vector<std::size_t> solve_tmpl();
  std::unique_ptr<keyvaluedb::KeyValueDB> &get_db();
};

#define PYRED_PP_CCS_GETNUMSYS(k)                                    \
  template <>                                                        \
  inline SystemOfEqs<PYRED_PP_COEFFCLASS(k)> &System::get_numsys() { \
    return m_numsys.PYRED_PP_COEFFCLASSMEM(k);                       \
  }
PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_GETNUMSYS)

#define PYRED_PP_CCS_GETCONSTNUMSYS(k)                               \
  template <>                                                        \
  inline const SystemOfEqs<PYRED_PP_COEFFCLASS(k)> &                 \
      System::get_numsys() const {                                   \
    return m_numsys.PYRED_PP_COEFFCLASSMEM(k);                       \
  }
PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_GETCONSTNUMSYS)

#define PYRED_PP_CCS_GETCOEFFMAP(k)                              \
  template <>                                                    \
  inline std::unordered_map<PYRED_PP_COEFFCLASS(k), std::string> \
      &System::get_coeff_map() {                                 \
    return m_coeff_map.PYRED_PP_COEFFCLASSMEM(k);                \
  }
PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_GETCOEFFMAP)

template <typename Coeff>
std::vector<std::size_t> System::generate_solve_tmpl(
    const std::vector<SeedSpec> &ibp_seedspec,
    const std::vector<SeedSpec> &ibp_seedcompl,
    const std::vector<SeedSpec> &sym_seedspec) {
  auto &numsys = get_numsys<Coeff>();
  numsys.clear();
  numsys = GeneratorHelper::generate_and_solve<Coeff>(
      ibp_seedspec, ibp_seedcompl, sym_seedspec, Config::auto_symseed(),
      Config::parallel());
  auto indep_eqnums = numsys.independent();
  auto &cmap = get_coeff_map<Coeff>();
  cmap = IntegralRelations::cache<Coeff>().get();
  return indep_eqnums;
}

template <typename Coeff>
void System::add_forward_tmpl(const eqdata &eq, SystemOfEqs<Coeff> &numsys,
                              sol_map<Coeff> &sols, std::size_t neq) {
  if (neq == std::numeric_limits<std::size_t>::max()) {
    neq = neqs;
  }
  if (!(neqs++)) {
    numsys.reserve(sys.capacity());
    sys.clear();
  }
  wc_equation<Coeff> tmpeq;
  tmpeq.reserve(eq.size());
  for (const auto &ic : eq) {
    auto cf = parse_coeff<Coeff>(ic.second);
    if (cf) {
      tmpeq.emplace_back(ic.first, cf);
    }
  }
  if (numsys.size() == numsys.capacity()) {
    throw init_error("pyred::System() grew beyond reserved capacity.");
  }
  numsys.add(std::move(tmpeq), neq);
  auto &cureq = numsys.sys.back();
  cureq.solve(sols, true, numsys.get_db());
  if (!cureq.empty()) {
    sols.insert({cureq.front().first, std::ref(cureq)});
  }
  else {
    numsys.sys.pop_back();
  }
}

template <typename Coeff>
std::vector<std::size_t> System::solve_tmpl() {
  auto &numsys = get_numsys<Coeff>();
  if (sys.size() && !numsys.size()) {
    // We are dealing with a system of (weight,string) equations.
    // Need to parse the coefficients.
    auto time_parse_begin = std::clock();
    Config::log(1) << "parse coefficients" << std::flush;
    numsys = SystemOfEqs<Coeff>(sys);
    auto time_parse_end = std::clock();
    Config::log(1) << ": " << time_diff(time_parse_begin, time_parse_end)
                   << " (" << numsys.size() << " equations)" << std::endl;
  }
  else if (!sys.size() && numsys.size()) {
    // Keep the system as it is.
  }
  else {
    throw init_error("pyred::System(): equations with string coefficients "
                     "and equations with finite integer coefficients must not "
                     "be added to the same system.");
  }
  Config::log(1) << "solve forward" << std::flush;
  auto time_fwd_begin = std::clock();

  numsys.setup_insertions_db(numsys.size());
  auto nequations = numsys.size();
//  auto maxinsertions = numsys.solve(false,true);
  auto maxinsertions = numsys.solve(true,true);
  // SystemOfEqs<Coeff>::solve() sorts numsys.
  m_solver_sorted = true;
  auto indep_eqnums = numsys.independent();
  auto neqs_indep = indep_eqnums.size();

  auto time_fwd_end = std::clock();
  Config::log(1) << ":      " << time_diff(time_fwd_begin, time_fwd_end)
                 << std::endl;
  Config::log(1) << "max insertions:     " << maxinsertions << std::endl;

  if (Config::backward()) {
    backward_tmpl<Coeff>(numsys);
    auto time_solve_end = std::clock();
    Config::log(1) << "total solve time:   "
                   << time_diff(time_fwd_begin, time_solve_end) << std::endl;
  }
  else {
    numsys.clear();
  }
  Config::log(1) << std::setw(8) << nequations
                 << " equations: " << nequations - neqs_indep << " zero + "
                 << neqs_indep << " independent" << std::endl;
  return indep_eqnums;
}

template <typename Coeff>
void System::backward_tmpl(SystemOfEqs<Coeff> &numsys) {
  Config::log(1) << "solve backward" << std::flush;
  auto time_backward_begin = std::clock();
  Equation<Coeff> ban(1);
  ban.banMap=1;
  numsys.solve();
  auto time_backward_end = std::clock();
  Config::log(1) << ": " << time_diff(time_backward_begin, time_backward_end)
                 << std::endl;
  ban.banMap=0;
  // Fill the unordered_map m_content:
  // Weight -> pair(eqnum, {unreduced_integrals})
  // keep 'insertions' in the database where they are already.
  // Delete equations as soon as they have been processed.
  if (Config::insertion_tracer() != 0) {
    Config::log(1) << "prepare content for selector" << std::flush;
    auto time_content_begin = std::clock();
    for (auto &eq : numsys.sys) {
      auto itchk = m_content.insert(
        {eq.front().first, std::make_pair(eq.eqnum, std::vector<Weight>{})}
      ).first;
      auto &cont = itchk->second.second;
      cont.reserve(eq.size() - 1);
      // Empty equations have already been dropped by independent().
      for (auto it = eq.cbegin() + 1; it != eq.cend(); ++it) {
        cont.emplace_back(it->first);
      }
      eq.clear_eq();
    }
    auto time_content_end = std::clock();
    Config::log(1) << ": " << time_diff(time_content_begin, time_content_end)
                   << std::endl;
  }
  m_content_prepared = true;
  numsys.clear();
}

template <typename Coeff>
void System::backward_coeff_tmpl(SystemOfEqs<Coeff> &numsys) {
  Config::log(1) << "solve backward and save the coefficients" << std::flush;
  auto time_backward_begin = std::clock();
  numsys.solve();
  auto time_backward_end = std::clock();
  Config::log(1) << ": " << time_diff(time_backward_begin, time_backward_end)
                 << std::endl;
  // Fill the unordered_map m_content_coeff:
  // intid -> pair(eqnum, {<unreduced_integrals,coeff/FF>})
  // keep 'insertions' in the database where they are already.
  // Delete equations as soon as they have been processed.
  if (Config::insertion_tracer() != 0) {
    Config::log(1) << "prepare content for selector" << std::flush;
    auto time_content_begin = std::clock();
    for (auto &eq : numsys.sys) {
      auto itchk = m_content_coeff
                       .insert({eq.front().first,
                                std::make_pair(eq.eqnum, std::vector<std::tuple<Weight,std::string>>{})})
                       .first;
      auto &cont = itchk->second.second;
      cont.reserve(eq.size() - 1);
      // Empty equations have already been dropped by independent().
      for (auto it = eq.cbegin() + 1; it != eq.cend(); ++it) {
        std::ostringstream s1;
        s1 << it->second;
        std::string s2 = s1.str();
        cont.emplace_back(std::make_tuple(it->first,s2));
      }
      eq.clear_eq();
    }
    auto time_content_end = std::clock();
    Config::log(1) << ": " << time_diff(time_content_begin, time_content_end)
                   << std::endl;
  }
  m_content_prepared = true;
  numsys.clear();
}

} // namespace pyred

#endif
