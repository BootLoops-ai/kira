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

#ifndef PYRED_GAUSS_H
#define PYRED_GAUSS_H

#include <algorithm> // max, min, sort, copy, min_element
#include <cassert>
#include <cstddef>
#include <functional> // reference_wrapper
#include <iostream>
#include <iterator> // back_inserter
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "pyred/coeff_helper.h"
#include "pyred/defs.h"
#include "pyred/keyvaluedb.h"
#include "pyred/parser.h" // parse_coeff

namespace pyred {

template <typename Coeff>
class Equation;

template <typename Coeff>
class EquationSolver;

template <typename Coeff>
using sol_map =
    std::unordered_map<Weight, std::reference_wrapper<Equation<Coeff>>>;

template <typename Coeff>
bool cmp_wc_pair(const wc_pair<Coeff>& a, const wc_pair<Coeff>& b) {
  /*
  operator< to compare integral-coefficient pairs by integral.
  Makes std::sort place higher integrals first.
  */
  return a.first > b.first;
}

template <typename Coeff>
bool cmp_eqn(const Equation<Coeff>& a, const Equation<Coeff>& b) {
  /*
  Compare Equations (like operator<) by
  * highest integral: lower first
  * length: shorter first
  * if same highest integral and same length:
    first if lower integrals following
  * if all integrals the same:
    lower equation number first.
  * if an equation is empty, place it last.
  Makes std::sort place lower equations first.
  */
  if (a.empty()) return false;
  if (b.empty()) return true;
  if (a.front().first != b.front().first) {
    return a.front().first < b.front().first;
  }
  if (a.size() != b.size()) return a.size() < b.size();
  for (std::size_t i = 1; i != a.size(); ++i) {
    if (a.eq[i].first != b.eq[i].first) {
      return a.eq[i].first < b.eq[i].first;
    }
  }
  return a.eqnum < b.eqnum;
}

template <typename Coeff>
bool cmp_eqdata(const wc_equation<Coeff> &a, const wc_equation<Coeff> &b) {
  /*
  Compare equation data (like operator<) by
  * highest integral: lower first
  * length: shorter first
  * if same highest integral and same length:
    first if lower integrals following
  * if an equation is empty, place it last.
  Makes std::sort place lower equations first.
  NOTE:
  Equations with equal integral content are regarded as sorted (equal).
  Since there is no equation number to make the ordering unique,
  cmp_eqdata() should only be used with std::stable_sort().
  The result will then be the same as with cmp_eqn() and ascending
  equation numbers.
  */
  if (a.empty()) return false;
  if (b.empty()) return true;
  if (a.front().first != b.front().first) {
    return a.front().first < b.front().first;
  }
  if (a.size() != b.size()) return a.size() < b.size();
  for (std::size_t i = 1; i != a.size(); ++i) {
    if (a[i].first != b[i].first) {
      return a[i].first < b[i].first;
    }
  }
  return false;
}

template <typename Coeff>
bool cmp_eqn_fwd(const Equation<Coeff>& a, const Equation<Coeff>& b) {
  /*
  Compare Equations (like operator<) by
  * highest integral: higher first
  * length: shorter first
  * if same highest integral and same length:
    first if lower integrals following
  * if all integrals the same:
    lower equation number first.
  * if an equation is empty, place it last.
  This is the ordering for Kira's forward elimination.
  */
  if (a.empty()) return false;
  if (b.empty()) return true;
  if (a.front().first != b.front().first) {
    return a.front().first > b.front().first;
  }
  if (a.size() != b.size()) return a.size() < b.size();
  for (std::size_t i = 1; i != a.size(); ++i) {
    if (a.eq[i].first != b.eq[i].first) {
      return a.eq[i].first < b.eq[i].first;
    }
  }
  return a.eqnum < b.eqnum;
}

template <typename Coeff>
struct cmp_eqn_fwd_heap {
  bool operator()(const Equation<Coeff>& a, const Equation<Coeff>& b) {
    return cmp_eqn_fwd(b, a);
  }
};

template <typename Coeff>
bool cmp_eqn_vec(const wc_equation<Coeff>& a, const wc_equation<Coeff>& b) {
  /*
  Compare wc_equation (equations as vectors).
  Same as comparing Equation objects, but without an equation number.
  I.e. the ordering is not unique.
  */
  if (a.empty()) return false;
  if (b.empty()) return true;
  if (a.front().first != b.front().first) {
    return a.front().first < b.front().first;
  }
  if (a.size() != b.size()) return a.size() < b.size();
  for (std::size_t i = 1; i != a.size(); ++i) {
    if (a[i].first != b[i].first) {
      return a[i].first < b[i].first;
    }
  }
  return false;
}

template <typename Coeff>
bool cmp_eqn_psolve(const Equation<Coeff>& a, const Equation<Coeff>& b) {
  /*
  Compare Equations by
  * highest integral: lower first
  * length: longer first
  * if same highest integral and same length:
    first if higher integrals following
  * if all integrals the same:
    lower equation number first.
  * if an equation is empty, place it last.
  Makes std::sort place lower equations first.
  */
  if (a.empty()) {
    return false;
  }
  if (b.empty()) {
    return true;
  }
  if (a.front().first < b.front().first) {
    return true;
  }
  if (a.front().first > b.front().first) {
    return false;
  }
  if (a.size() < b.size()) {
    return false;
  }
  if (a.size() > b.size()) {
    return true;
  }
  for (std::size_t i = 1; i != a.size(); ++i) {
    if (a.eq[i].first < b.eq[i].first) {
      return false;
    }
    if (a.eq[i].first > b.eq[i].first) {
      return true;
    }
  }
  return a.eqnum > b.eqnum;
}

template <typename Coeff>
bool cmp_eqn_inscount(const std::unique_ptr<EquationSolver<Coeff>>& a,
                      const std::unique_ptr<EquationSolver<Coeff>>& b) {
  if (a->fwd_insertions != b->fwd_insertions) {
    return a->fwd_insertions > b->fwd_insertions;
  }
  if (a->origlength != b->origlength) {
    return a->origlength > b->origlength;
  }
  return a->eqnum > b->eqnum;
}

template <typename Coeff>
struct cmp_eqn_heap {
  bool operator()(const std::shared_ptr<EquationSolver<Coeff>>& a,
                  const std::shared_ptr<EquationSolver<Coeff>>& b) {
    //   bool operator()(const std::unique_ptr<EquationSolver<Coeff>>& a,
    //                          const std::unique_ptr<EquationSolver<Coeff>>& b)
    //                          {
    if (a->front().first < b->front().first) {
      return true;
    }
    if (a->front().first > b->front().first) {
      return false;
    }
    if (a->fwd_insertions != b->fwd_insertions) {
      return a->fwd_insertions > b->fwd_insertions;
    }
    if (a->origlength != b->origlength) {
      return a->origlength > b->origlength;
    }
    if (a->fwd_insertions == 0u) {
      for (std::size_t i = 1; i != a->eq.size(); ++i) {
        if (a->eq[i].first < b->eq[i].first) {
          return false;
        }
        if (a->eq[i].first > b->eq[i].first) {
          return true;
        }
      }
    }
    return a->eqnum > b->eqnum;
  }
};

template <typename Coeff>
class Equation {
public:
  static int banMap;
  // If it is a solution, the first element is (integral, -1).
  wc_equation<Coeff> eq;
  // equation number (should start with 1).
  // could use eqnum = 0 if unassigned
  std::size_t eqnum;
  // all constructors set eqnum
  Equation() = delete;
  Equation(std::size_t eqnum);
  template <typename T>
  Equation(const T& eqin, std::size_t eqnum);
  Equation(wc_equation<Coeff>&& eqin, std::size_t eqnum);
  Equation(EquationSolver<Coeff>&& sol);
  void move_from(EquationSolver<Coeff>&& sol);
  std::size_t solve(const sol_map<Coeff>& sols, bool insertall,
                    std::unique_ptr<keyvaluedb::KeyValueDB>&);
  bool empty() const { return eq.empty(); }
  std::size_t size() const { return eq.size(); }
  void reserve(std::size_t new_cap) { eq.reserve(new_cap); }
  wc_pair<Coeff>& at(std::size_t pos) { return eq.at(pos); }
  const wc_pair<Coeff>& at(std::size_t pos) const { return eq.at(pos); }
  wc_pair<Coeff>& front() { return eq.front(); }
  const wc_pair<Coeff>& front() const { return eq.front(); }
  wc_pair<Coeff>& back() { return eq.back(); }
  const wc_pair<Coeff>& back() const { return eq.back(); }
  typename wc_equation<Coeff>::iterator begin() { return eq.begin(); }
  typename wc_equation<Coeff>::const_iterator cbegin() const { return eq.cbegin(); }
  typename wc_equation<Coeff>::iterator end() { return eq.end(); }
  typename wc_equation<Coeff>::const_iterator cend() const { return eq.cend(); }
  typename wc_equation<Coeff>::iterator
    erase(typename wc_equation<Coeff>::iterator pos) { return eq.erase(pos); }
  void clear_eq() { eq.clear(); }
  std::vector<Weight> get_insertions(std::unique_ptr<keyvaluedb::KeyValueDB>&);
  void set_insertions(std::unique_ptr<keyvaluedb::KeyValueDB>&,
                      const std::vector<Weight>&) const;
  void clear_insertions(std::unique_ptr<keyvaluedb::KeyValueDB>&);
  void insert_one(const Equation<Coeff>&);
};

template <typename Coeff>
int Equation<Coeff>::banMap = 0;

template <typename Coeff>
class EquationSolver {
  friend struct cmp_eqn_heap<Coeff>;

private:
  wc_equation<Coeff> eq;

public:
  std::size_t eqnum;
  // temporary result
  wc_equation<Coeff> tmpeq;
  std::size_t origlength;
  // Number of inserted equations before the equation is solved.
  std::size_t fwd_insertions{0u};
  // solutions to insert as (iterator, end, normfactor)
  std::vector<std::tuple<typename wc_equation<Coeff>::const_iterator,
                         typename wc_equation<Coeff>::const_iterator, Coeff>>
      neededsols;
  std::vector<Weight> m_inserted;
  // highest integral in the solutions which has not yet been treated
  Weight hi_insols;
  EquationSolver(Equation<Coeff>&);
  void insert(const Equation<Coeff>&);
  bool proceed();
  bool compactify();
  void normalise();
  bool empty() { return tmpeq.empty(); }
  wc_pair<Coeff>& front() { return tmpeq.front(); }
};


template <typename Coeff>
class SystemOfEqs {
public:
  std::vector<Equation<Coeff>> sys;
  SystemOfEqs() = default;
  SystemOfEqs(const std::vector<eqdata>& eqs);
  SystemOfEqs(
      const std::vector<std::vector<std::pair<Weight, std::size_t>>>& system,
      const std::vector<Coeff>& funs, bool solve_otf = false);
  void reserve(std::size_t sz);
  Equation<Coeff>& add(wc_equation<Coeff>&& eqdata, std::size_t eqnum) {
    sys.emplace_back(std::move(eqdata), eqnum);
    return sys.back();
  }
  Equation<Coeff>& add(Equation<Coeff>&& eq) {
    sys.push_back(std::move(eq));
    return sys.back();
  }
  void sort();
  void sort(typename std::vector<Equation<Coeff>>::iterator,
            typename std::vector<Equation<Coeff>>::iterator);
  void sort_fwd();
  void sort_fwd(typename std::vector<Equation<Coeff>>::iterator,
                typename std::vector<Equation<Coeff>>::iterator);
  std::size_t solve(bool dosort = true, bool insertall = true);
  std::size_t psolve_map(
      std::size_t insertion_limit = std::numeric_limits<std::size_t>::max(),
      std::size_t max_solversize = std::numeric_limits<std::size_t>::max());
  // std::size_t psolve_heap();
  Equation<Coeff>& operator[](std::size_t pos) { return sys[pos]; }
  std::size_t size() const { return sys.size(); }
  bool empty() const { return sys.empty(); }
  std::size_t capacity() const { return sys.capacity(); }
  typename std::vector<Equation<Coeff>>::iterator begin() {
    return sys.begin();
  }
  typename std::vector<Equation<Coeff>>::iterator end() { return sys.end(); }
  void clear() { sys.clear(); }
  std::size_t truncate();
  std::vector<std::size_t> independent();
  std::vector<Weight> unreduced() const;
  std::unique_ptr<keyvaluedb::KeyValueDB>& get_db() { return m_db; }
  void setup_insertions_db(std::size_t);
  // Kira's forward elimination
  bool iterate_forward();
  //   std::size_t solve_forward();
  void solve_forward();
  static std::unordered_set<Weight> post_select(
      std::vector<std::vector<std::pair<Weight, std::size_t>>>& /*sys*/,
      const std::vector<Coeff>& /*coeffs*/,
      const std::unordered_set<Weight>& /*zeromasters*/);

private:
  std::unique_ptr<keyvaluedb::KeyValueDB> m_db;
};

/************
 * Equation *
 ************/

template <typename Coeff>
Equation<Coeff>::Equation(std::size_t eqnum) : eqnum(eqnum) {}

template <typename Coeff>
template <typename T>
Equation<Coeff>::Equation(const T& eqin, std::size_t eqnum) : eqnum(eqnum) {
  reserve(eqin.size());
  std::copy(eqin.cbegin(), eqin.cend(), std::back_inserter(eq));
  std::sort(begin(), end(), cmp_wc_pair<Coeff>);
  // Check that there are no duplicate integrals.
  Weight last_igl = Weight::none();
  for (const auto &w_c: eq) {
    if (w_c.first == last_igl) {
      throw input_error("Duplicate integrals in Equation().");
    }
    last_igl = w_c.first;
  }
}

template <typename Coeff>
Equation<Coeff>::Equation(wc_equation<Coeff>&& eqin, std::size_t eqnum)
    : eq(std::move(eqin)), eqnum(eqnum) {
  // Check that the integrals are sorted and that there are no duplicates.
  Weight last_igl = Weight::none();
  for (const auto &w_c: eq) {
    if (w_c.first >= last_igl) {
      throw input_error("Integrals in Equation() are not sorted/merged.");
    }
    last_igl = w_c.first;
  }
}

template <typename Coeff>
Equation<Coeff>::Equation(EquationSolver<Coeff>&& sol)
    : eq{std::move(sol.tmpeq)}, eqnum{sol.eqnum} {}

template <typename Coeff>
void Equation<Coeff>::move_from(EquationSolver<Coeff>&& sol) {
  // Does not copy the list of inserted equations.
  // They are retained, however, when the moved equation belongs
  // to the same SystemOfEqs.
  eq = std::move(sol.tmpeq);
  eqnum = sol.eqnum;
}

template <typename Coeff>
std::size_t Equation<Coeff>::solve(const sol_map<Coeff>& sols, bool insertall,
                             std::unique_ptr<keyvaluedb::KeyValueDB>& db) {
  auto solver = EquationSolver<Coeff>(*this);
  bool untreated = true;
  bool unnormalised = true;
  while (untreated) {
    if (unnormalised || insertall) {
      auto sol = sols.find(solver.tmpeq.back().first);
      if (sol != sols.cend()) {
        solver.insert(sol->second.get());
        if (unnormalised) {
          ++solver.fwd_insertions;
        }
      }
      else if (unnormalised) {
        solver.normalise();
        unnormalised = false;
      }
    }
    untreated = solver.proceed();
  }
  eq = std::move(solver.tmpeq);
  if (eq.empty()) {
    solver.fwd_insertions = 0u;
    clear_insertions(db);
  }
  else if(banMap==0){
    auto ins = get_insertions(db);
    if (ins.empty()) {
      set_insertions(db, solver.m_inserted);
    }
    else {
      ins.reserve(ins.size() + solver.m_inserted.size());
      for (const auto &iid: solver.m_inserted)
        ins.push_back(iid);
      set_insertions(db, ins);
    }
  }
  return solver.fwd_insertions;
}

template <typename Coeff>
std::vector<Weight> Equation<Coeff>::get_insertions(
    std::unique_ptr<keyvaluedb::KeyValueDB>& db) {
  return db->get(eqnum, false);
}

template <typename Coeff>
void Equation<Coeff>::set_insertions(
    std::unique_ptr<keyvaluedb::KeyValueDB>& db,
    const std::vector<Weight>& ins) const {
  db->put(eqnum, ins);
}

template <typename Coeff>
void Equation<Coeff>::clear_insertions(
    std::unique_ptr<keyvaluedb::KeyValueDB>& db) {
  db->remove(eqnum, false);
}

template <typename Coeff>
void Equation<Coeff>::insert_one(const Equation<Coeff>& other) {
  // Insert the equation 'other' to eliminate the highest integral in *this.
  // Both equations must have the same highest integral.
  // Neither equation needs to be normalised. And neither will be the result.
  assert(!this->empty() && !other.empty() &&
         this->front().first == other.front().first);
  wc_equation<Coeff> tmpeq;
  auto tit = this->cbegin();
  auto tend = this->cend();
  auto oit = other.cbegin();
  auto oend = other.cend();
  auto fac = -tit->second / oit->second;
  ++tit;
  ++oit;
  while (true) {
    if (oit == oend) {
      // This also catches the case where tit==tend at the same time.
      for (; tit != tend; ++tit) {
        tmpeq.emplace_back(*tit);
      }
      break;
    }
    else if (tit == tend) {
      for (; oit != oend; ++oit) {
        tmpeq.emplace_back(oit->first, fac * oit->second);
      }
      break;
    }
    else {
      if (tit->first == oit->first) {
        auto cf = tit->second + fac * oit->second;
        if (cf) {
          tmpeq.emplace_back(tit->first, cf);
        }
        ++tit;
        ++oit;
      }
      else if (tit->first > oit->first) {
        tmpeq.push_back(*tit);
        ++tit;
      }
      else { // tit->first < oit->first
        tmpeq.emplace_back(oit->first, fac * oit->second);
        ++oit;
      }
    }
  }
  eq = std::move(tmpeq);
}

/******************
 * EquationSolver *
 ******************/

template <typename Coeff>
EquationSolver<Coeff>::EquationSolver(Equation<Coeff>& eqin)
    : eq{std::move(eqin.eq)},
      eqnum{eqin.eqnum},
      origlength{eq.size()} {
  // Note: invalidates eqin; eqin must not be empty
  auto snd = eq.cbegin() + 1;
  if (snd != eq.cend()) {
    neededsols.emplace_back(snd, eq.cend(), Coeff(1));
    hi_insols = snd->first;
  }
  else {
    hi_insols = Weight::none();
  }
  tmpeq.push_back(eq.front());
}

template <typename Coeff>
void EquationSolver<Coeff>::insert(const Equation<Coeff>& addsol) {
  // Add addsol to a solver.
  // addsol must be the solution of the lowest (=last) integral
  // (i.e. it has coefficient -1) in tmpeq.
  // Remember which integrals were inserted.
  m_inserted.push_back(addsol.front().first);
  auto solit = addsol.cbegin() + 1;
  if (solit != addsol.cend()) {
    neededsols.emplace_back(solit, addsol.cend(), tmpeq.back().second);
    if (hi_insols == Weight::none()) {
      hi_insols = solit->first;
    }
    else {
      hi_insols = std::max(hi_insols, solit->first);
    }
  }
  tmpeq.pop_back();
}

template <typename Coeff>
bool EquationSolver<Coeff>::proceed() {
  // Calculate the next non-vanishing coefficient.
  // Return true if a non-zero term was calculated, false otherwise
  // (false implies that all inserted solutions are depleted).
  tmpeq.emplace_back(hi_insols, Coeff(0));
  while (hi_insols != Weight::none()) {
    auto& lastcoeff = tmpeq.back().second;
    Weight next_hi_insols{Weight::none()};
    auto sol = neededsols.begin();
    while (sol < neededsols.end()) {
      // add contribution to the current integral from all neededsols
      auto& solit = std::get<0>(*sol);
      auto& solend = std::get<1>(*sol);
      if (solit->first == hi_insols) {
        // sol contains the current integral -> insert
        lastcoeff += solit->second * std::get<2>(*sol);
        ++solit;
      }
      if (solit == solend) {
        // Nothing more to insert from this solution:
        // Remove it, place the solution from the end of neededsols here
        // (note that this might be the current solution)
        // and remove the latter from the end.
        *sol = neededsols.back();
        neededsols.pop_back();
      }
      else {
        if (next_hi_insols == Weight::none()) {
          next_hi_insols = solit->first;
        }
        else {
          next_hi_insols = std::max(next_hi_insols, solit->first);
        }
        ++sol;
      }
    }
    hi_insols = next_hi_insols;
    if (!lastcoeff) {
      tmpeq.back().first = hi_insols;
    }
    else {
      break;
    }
  }
  if (tmpeq.back().first == Weight::none()) {
    tmpeq.pop_back();
    return false;
  }
  return true;
}

template <typename Coeff>
bool EquationSolver<Coeff>::compactify() {
  // Insert all solutions from 'neededsols' until they are depleted.
  // Return true if the equation is non-empty, false otherwise.
  // As long as the equation is not normalised,
  // this can be used instead of proceed().
  while (proceed())
    ;
  neededsols.clear();
  eq.clear();
  if (tmpeq.empty()) return false;
  eq.swap(tmpeq);
  tmpeq.push_back(eq.front());
  auto snd = eq.cbegin() + 1;
  if (snd != eq.cend()) {
    neededsols.emplace_back(snd, eq.cend(), Coeff(1));
    hi_insols = snd->first;
  }
  else {
    hi_insols = Weight::none();
  }
  return true;
}

template <typename Coeff>
void EquationSolver<Coeff>::normalise() {
  // Adjust the prefactors of the solutions so that the equation gets solved
  // for the integral in tmpeq. Set the coefficient of the integral to -1.
  auto inversefac = Coeff(-1) / tmpeq.back().second;
  tmpeq.back().second = Coeff(-1);
  for (auto& sol : neededsols) {
    std::get<2>(sol) *= inversefac;
  }
}

/***************
 * SystemOfEqs *
 ***************/

template <typename Coeff>
SystemOfEqs<Coeff>::SystemOfEqs(const std::vector<eqdata>& eqs) {
  std::size_t neqs{0};
  reserve(eqs.size());
  for (const auto& eq : eqs) {
    wc_equation<Coeff> tmpeq;
    tmpeq.reserve(eq.size());
    for (const auto& ic : eq) {
      auto cf = parse_coeff<Coeff>(ic.second);
      if (cf) {
        tmpeq.emplace_back(ic.first, cf);
      }
    }
    // Do not std::move(tmpeq) here, because tmpeq might not be sorted.
    sys.emplace_back(tmpeq, neqs++);
  }
}

template <typename Coeff>
void SystemOfEqs<Coeff>::reserve(std::size_t sz) {
  sys.reserve(sz);
  setup_insertions_db(sz);
}

template <typename Coeff>
void SystemOfEqs<Coeff>::sort() {
  std::sort(sys.begin(), sys.end(), cmp_eqn<Coeff>);
}

template <typename Coeff>
void SystemOfEqs<Coeff>::sort_fwd() {
  std::sort(sys.begin(), sys.end(), cmp_eqn_fwd<Coeff>);
}

template <typename Coeff>
void SystemOfEqs<Coeff>::sort(
    typename std::vector<Equation<Coeff>>::iterator first,
    typename std::vector<Equation<Coeff>>::iterator last) {
  std::sort(first, last, cmp_eqn<Coeff>);
}

template <typename Coeff>
void SystemOfEqs<Coeff>::sort_fwd(
    typename std::vector<Equation<Coeff>>::iterator first,
    typename std::vector<Equation<Coeff>>::iterator last) {
  std::sort(first, last, cmp_eqn_fwd<Coeff>);
}

template <typename Coeff>
std::size_t SystemOfEqs<Coeff>::solve(bool dosort, bool insertall) {
  if (dosort) {
    sort();
  }
  sol_map<Coeff> solmap;
  std::size_t maxinsertions{0};
  for (auto& eq : sys) {
    auto ninsertions = eq.solve(solmap, insertall, m_db);
    maxinsertions = std::max(maxinsertions, ninsertions);
    if (!eq.empty()) {
      solmap.insert({eq.front().first, std::ref(eq)});
    }
  }
  return maxinsertions;
}

template <typename Coeff>
void SystemOfEqs<Coeff>::setup_insertions_db(std::size_t sz) {
  auto filename = Config::database_file().first;
  auto overwrite = Config::database_file().second;
  if (Config::insertion_tracer() == 0) {
    m_db = std::make_unique<keyvaluedb::KeyValueDiscard>("", 0, overwrite);
  }
  else if (Config::insertion_tracer() == 1) {
    m_db = std::make_unique<keyvaluedb::KeyValueVector>("", sz, overwrite);
  }
  else if (Config::insertion_tracer() == 2) {
    m_db = std::make_unique<keyvaluedb::KeyValueSQLite>(filename + ".db", sz,
                                                        overwrite);
  }
#ifdef PYRED_KCDB
  else if (Config::insertion_tracer() == 3) {
    m_db = std::make_unique<keyvaluedb::KeyValueKC>(filename + ".kch", sz,
                                                    overwrite);
  }
#endif
}

template <typename Coeff>
std::size_t SystemOfEqs<Coeff>::psolve_map(std::size_t insertion_limit,
                                           std::size_t max_solversize) {
  using sol_uptr_vec = std::vector<std::unique_ptr<EquationSolver<Coeff>>>;
  std::size_t maxinsertions{0};
  if (sys.empty()) return maxinsertions;
  std::unordered_map<Weight, sol_uptr_vec> solvers;
  std::sort(sys.begin(), sys.end(), cmp_eqn_psolve<Coeff>);
  auto sysit = sys.rbegin();
  auto solsysit = sys.rbegin();
  while (sysit < sys.rend() || !solvers.empty()) {
    Weight hi{Weight::none()};
    if (sysit < sys.rend()) {
      hi = sysit->front().first;
    }
    Weight hi_insolvers{Weight::none()};
    // Look up highest in solvers.
    for (const auto& s : solvers) {
      if (hi_insolvers == Weight::none()) {
        hi_insolvers = s.first;
      }
      else {
        hi_insolvers = std::max(hi_insolvers, s.first);
      }
    }
    typename std::unordered_map<Weight, sol_uptr_vec>::iterator mapit;
    // typename sol_uptr_vec::iterator soleqit; // doesn't work with gcc 4.7
    typename std::vector<std::unique_ptr<EquationSolver<Coeff>>>::iterator
        soleqit;
    if ((hi >= hi_insolvers && hi != Weight::none()) ||
        hi_insolvers == Weight::none()) {
      // Solve the top equation from sys,
      // because other than those in solvers it has no insertions.
      auto sol = EquationSolver<Coeff>(*sysit);
      sol.normalise();
      while (sol.proceed())
        ;
      maxinsertions = std::max(maxinsertions, sol.fwd_insertions);
      solsysit->move_from(std::move(sol));
      // Create a new entry in solvers or return an iterator
      // to the existing (vector) entry.
      mapit = solvers.emplace(hi, sol_uptr_vec{}).first;
      while (++sysit < sys.rend() && sysit->front().first == hi) {
        // Insert following equations from sys with the same
        // highest integral into solsys.
        mapit->second.emplace_back(
            std::make_unique<EquationSolver<Coeff>>(*sysit));
      }
      // Solved equation was not in solvers (see below).
      soleqit = mapit->second.end();
    }
    else {
      // Choose the equation to solve from solvers.
      mapit = solvers.find(hi_insolvers);
      soleqit = std::max_element(mapit->second.begin(), mapit->second.end(),
                                 cmp_eqn_inscount<Coeff>);
      (*soleqit)->normalise();
      while ((*soleqit)->proceed())
        ;
      maxinsertions = std::max(maxinsertions, (*soleqit)->fwd_insertions);
      solsysit->move_from(std::move(*(*soleqit)));
    }
    auto& sol = *solsysit;
    ++solsysit;
    for (auto it = mapit->second.begin(); it != mapit->second.end(); ++it) {
      // Do not insert the solution into itself.
      if (it == soleqit) continue;
      // Insert solution and move non-vanishing equations
      // to their appropriate places in solvers.
      auto& eqptr = *it;
      eqptr->insert(sol);
      if (eqptr->neededsols.size() > max_solversize) {
        eqptr->compactify();
      }
      else {
        eqptr->proceed();
      }
      ++(eqptr->fwd_insertions);
      if (!eqptr->empty() && (eqptr->fwd_insertions <= insertion_limit)) {
        // Create solver entry if it doesn't exist.
        // TODO: avoid repeated look-ups of the same element.
        auto insit =
            solvers.emplace(eqptr->tmpeq.front().first, sol_uptr_vec{}).first;
        insit->second.emplace_back(std::move(eqptr));
      }
    }
    solvers.erase(mapit);
  }
  while (solsysit < sys.rend()) {
    solsysit->eq.clear();
    ++solsysit;
  }
  return maxinsertions;
}

// template<typename Coeff>
// std::size_t SystemOfEqs<Coeff>::psolve_heap() {
//   std::size_t maxinsertions{0};
//   std::priority_queue<std::unique_ptr<EquationSolver<Coeff>>,
//                       std::vector<std::shared_ptr<EquationSolver<Coeff>>>,
// //                       std::vector<std::unique_ptr<EquationSolver<Coeff>>>,
//                       cmp_eqn_heap<Coeff>> sysheap;
//   std::vector<Equation<Coeff>> solsys;
//   solsys.reserve(sys.size());
//   for (auto& eq: sys) {
//     sysheap.emplace(new EquationSolver<Coeff>(eq));
//   }
//   while (!sysheap.empty()) {
//     auto solptr = sysheap.top();
//     sysheap.pop();
//     solptr->normalise();
//     while (solptr->proceed());
//     maxinsertions = std::max(maxinsertions, solptr->fwd_insertions);
//     solsys.emplace_back(*solptr);
//     auto hi = solsys.back().front().first;
//     while (!sysheap.empty() && sysheap.top()->front().first == hi) {
//       auto eqptr = sysheap.top();
//       sysheap.pop();
//       eqptr->insert(solsys.back());
//       ++(eqptr->fwd_insertions);
//       if (eqptr->proceed()) {
//         sysheap.push(eqptr);
//       }
//     }
//   }
//   sys = std::move(solsys);
//   return maxinsertions;
// }

template <typename Coeff>
std::size_t SystemOfEqs<Coeff>::truncate() {
  sort();
  while (!sys.empty() && sys.back().empty()) {
    sys.pop_back();
  }
  return sys.size();
}

template <typename Coeff>
std::vector<std::size_t> SystemOfEqs<Coeff>::independent() {
  auto neqs = truncate();
  std::vector<std::size_t> eqnums;
  eqnums.reserve(neqs);
  for (const auto& eq : sys) {
    eqnums.push_back(eq.eqnum); // subtract one?
  }
  std::sort(eqnums.begin(), eqnums.end());
  return eqnums;
}

template <typename Coeff>
std::vector<Weight> SystemOfEqs<Coeff>::unreduced() const {
  std::unordered_set<Weight> unred_set;
  for (const auto &eq : sys) {
    if (eq.empty()) continue;
    for (auto it = eq.cbegin() + 1; it != eq.cend(); ++it) {
      unred_set.insert(it->first);
    }
  }
  std::vector<Weight> unred_vec;
  for (const auto &id: unred_set) {
    unred_vec.push_back(id);
  }
  std::sort(unred_vec.begin(), unred_vec.end());
  return unred_vec;
}

// Construct a SystemOfEqs from equations of the form
// {(weight, coeff_pos), ...}, where coeff_pos refers to the coefficient
// funs[coeff_pos].
// NOTE: The integrals in the equations in system must be sorted.
template <typename Coeff>
SystemOfEqs<Coeff>::SystemOfEqs(
    const std::vector<std::vector<std::pair<Weight, std::size_t>>>& system,
    const std::vector<Coeff>& funs, bool solve_otf) {
  std::size_t neqs{0};
  reserve(system.size());
  if (solve_otf) {
    sol_map<Coeff> sols;
    for (const auto& eq : system) {
      wc_equation<Coeff> tmp_eq;
      tmp_eq.reserve(eq.size());
      for (const auto& el : eq) {
        tmp_eq.emplace_back(el.first, funs[el.second]);
      }
      sys.emplace_back(std::move(tmp_eq), neqs++);
      auto &cureq = sys.back();
      cureq.solve(sols, true, get_db());
      if (!cureq.empty()) {
        sols.insert({cureq.front().first, std::ref(cureq)});
      }
      else {
        sys.pop_back();
      }
    }
  }
  else {
    for (const auto& eq : system) {
      wc_equation<Coeff> tmp_eq;
      tmp_eq.reserve(eq.size());
      for (const auto& el : eq) {
        tmp_eq.emplace_back(el.first, funs[el.second]);
      }
      sys.emplace_back(std::move(tmp_eq), neqs++);
    }
  }
}


// Solve the system to determine which integrals reduce to zero
// when zeromasters are set to zero.
// Set those integrals to zero in sys, modifying the argument sys.
// Perform the forward elimination to determine equations that reduce to zero.
// Remove those equations from sys.
template <typename Coeff>
std::unordered_set<Weight> SystemOfEqs<Coeff>::post_select(
    std::vector<std::vector<std::pair<Weight, std::size_t>>>& sys,
    const std::vector<Coeff>& coeffs,
    const std::unordered_set<Weight>& zeromasters) {
  std::unordered_set<Weight> zeroweights = zeromasters;
  // remove all zeromasters from the system
  for (auto& eq : sys) {
    auto ic_it = eq.begin();
    while (ic_it != eq.end()) {
      if (zeromasters.find(ic_it->first) != zeromasters.end()) {
        ic_it = eq.erase(ic_it);
      }
      else {
        ++ic_it;
      }
    }
  }
  {
    auto numsys = SystemOfEqs<Coeff>(sys, coeffs);
    numsys.truncate();
    numsys.solve(false);
    numsys.solve();
    for (const auto& eq : numsys.sys) {
      if (eq.size() == 1u) {
        zeroweights.insert(eq.front().first);
      }
    }
  }
  // remove all zeroweights from the system
  for (auto& eq : sys) {
    auto ic_it = eq.begin();
    while (ic_it != eq.end()) {
      if (zeroweights.find(ic_it->first) != zeroweights.end()) {
        ic_it = eq.erase(ic_it);
      }
      else {
        ++ic_it;
      }
    }
  }
  auto numsys = SystemOfEqs<Coeff>(sys, coeffs);
  numsys.truncate();
  numsys.solve(false);
  std::vector<std::vector<std::pair<Weight, std::size_t>>> newsys;
  newsys.reserve(numsys.sys.size());
  for (const auto &eq : numsys.sys) {
    if (!eq.empty()) {
      newsys.push_back(std::move(sys[eq.eqnum]));
      newsys.back().shrink_to_fit();
    }
  }
  sys = std::move(newsys);
  sys.shrink_to_fit();
  return zeroweights;
}


template <typename Coeff>
void SystemOfEqs<Coeff>::solve_forward() {
  std::priority_queue<Equation<Coeff>, std::vector<Equation<Coeff>>,
                      cmp_eqn_fwd_heap<Coeff>>
      queue;
  for (auto& eq : sys) {
    queue.push(std::move(eq));
  }
  sys.clear();
  sys.reserve(queue.size());
  auto last_igl = Weight::none();
  while (!queue.empty()) {
    auto& eq = queue.top();
    auto current_igl = eq.front().first;
    if (current_igl != last_igl) {
      // New subsystem.
      // This is the equation we insert into all others of the subsystem.
      last_igl = current_igl;
      sys.push_back(eq);
      // TODO: normalise
      queue.pop();
    }
    else {
      // Same subsystem as the last equation that was moved to sys.
      auto eq_cp = eq;
      queue.pop();
      eq_cp.insert_one(sys.back());
      queue.push(std::move(eq_cp));
    }
  }
  setup_insertions_db(size());
}

} // namespace pyred

#endif
