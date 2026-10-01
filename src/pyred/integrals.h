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

#ifndef PYRED_INTEGRALS_H
#define PYRED_INTEGRALS_H

#include <array>
#include <cstddef>
#include <exception>
#include <functional> // hash
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "pyred/defs.h"

namespace YAML {
class Node;
}

namespace pyred {

std::pair<Weight, std::string> split_intcoeff(const std::string &);

std::size_t bitwidth(std::size_t);

std::size_t count_set_bits(std::size_t);

std::size_t factorial(const std::size_t);

std::size_t binomial(int n, int k);

std::size_t dotsp_combs(const std::size_t, const std::size_t, const std::size_t,
                        const std::size_t);


class PropagatorPermutation {
public:
  PropagatorPermutation(const PropagatorPermutation &) = default;
  PropagatorPermutation(PropagatorPermutation &&) = default;
  PropagatorPermutation(std::vector<std::size_t> perm,
                        std::size_t n=std::numeric_limits<std::size_t>::max());
  PropagatorPermutation(std::size_t n=std::numeric_limits<std::size_t>::max())
    : m_permutation{}, m_n{n} {}
  PropagatorPermutation & operator=(const PropagatorPermutation &) = default;
  PropagatorPermutation & operator=(PropagatorPermutation &&) = default;
  std::size_t size() const { return m_n; }
  bool is_identity() const { return m_permutation.empty(); }
  // template<typename T>
  // std::vector<T> permute(const std::vector<T> vec, bool inverse = false);
  uint32_t permute_sector(uint32_t sector, bool inverse = false) const;
  typename std::vector<std::size_t>::iterator begin() {
    return m_permutation.begin();
  }
  typename std::vector<std::size_t>::iterator end() {
    return m_permutation.end();
  }
  typename std::vector<std::size_t>::const_iterator begin() const {
    return m_permutation.begin();
  }
  typename std::vector<std::size_t>::const_iterator end() const {
    return m_permutation.end();
  }
  typename std::vector<std::size_t>::const_iterator cbegin() const {
    return m_permutation.cbegin();
  }
  typename std::vector<std::size_t>::const_iterator cend() const {
    return m_permutation.cend();
  }
private:
  std::vector<std::size_t> m_permutation;
  std::size_t m_n;
};


/******************
 * Seed generator *
 ******************/

using compos_it_type = std::vector<std::vector<uint32_t>>::const_iterator;
using compos_rit_type =
    std::vector<std::vector<uint32_t>>::const_reverse_iterator;

inline uint32_t lines_from_sn(uint32_t sn) {
  uint32_t lines{0};
  while (sn) {
    if (sn & 1) {
      ++lines;
    }
    sn >>= 1;
  }
  return lines;
}

class SeedIterator;
class Topology;
class SeedSpecSector;
class IntegralRelations; // declared in relations.h

// Recursive seed specification for a top sector
// (i.e. all lower sector are included) and a maximal number of dots and sps
// depending on the number of lines in the respective sector.
// * mindots: the maximal number of dots on the highest sector.
// * n1 (default: #lines): the number of lines below which the number of dots
//   should be increased by one for each removed line.
// * maxdots (default: infinity): the maximal number of dots on each sector.
// Example: the highest sector has L lines and mindots=d. Then all sectors
// with n1 lines or more will get d dots. Sectors with n1-k (k>0) lines will
// get min(d+k,maxdots) dots.
class SeedSpec {
public:
  enum class Recursive { full, dotsp, no };
  std::shared_ptr<Topology> m_topo;
  uint32_t m_sector;
  uint32_t m_mindots;
  uint32_t m_sps;
  uint32_t m_n1;
  uint32_t m_maxdots;
  uint32_t m_lines;
  Recursive m_recursive;
  SeedSpec(const Topology &topo, uint32_t sector, uint32_t rmax, uint32_t sps,
           int maxdots = -1, int n1 = -1, Recursive rec = Recursive::full);
  static std::vector<SeedSpecSector> expand_sector(
      const std::vector<SeedSpec> &, const std::vector<SeedSpec> & = {},
      bool include_trivial = false);
  static std::vector<Integral> list_integrals(const std::vector<SeedSpec> &,
                                              int = -1);
  static std::vector<Weight> integral_selector(
      const std::vector<SeedSpec> &, int = -1);
};

// Seed specification for a particular sector
// with given (dots,sps) corners, skipping given (dots,sps) contributions.
class SeedSpecSector {
public:
  SeedSpecSector() = default;
  SeedSpecSector(const Topology &topo, uint32_t sector,
                 const std::vector<std::pair<uint32_t, uint32_t>> & = {});
  // methods
  void add_corner(uint32_t d, uint32_t s) { m_corners.push_back({d, s}); }
  void skip(uint32_t d, uint32_t s, bool recursive) {
    m_skip.emplace_back(d, s, recursive);
  }
  bool empty() const { return m_corners.empty(); }
  // Generate all (dot,sp) combinations within m_corners, skipping m_skip.
  std::vector<std::pair<uint32_t, uint32_t>> expand() const;
  // Cut m_corners at Config::max_symdots and Config::max_symsps.
  void cut();
  // members
  std::shared_ptr<Topology> m_topo;
  uint32_t m_sector;
  std::vector<std::pair<uint32_t, uint32_t>> m_corners;
  // skip (dots,sps,recursive)
  std::vector<std::tuple<uint32_t, uint32_t, bool>> m_skip;
};

bool operator<(const SeedSpecSector &, const SeedSpecSector &);

class Seeds {
public:
  typedef SeedIterator const_iterator;
  Seeds(const Seeds &) = default;
  Seeds(Seeds &&) = default;
  Seeds(uint32_t, uint32_t, uint32_t, uint32_t);
  SeedIterator begin() const;
  SeedIterator cbegin() const;
  SeedIterator end() const;
  SeedIterator cend() const;

private:
  std::vector<bool> sector;
  compos_rit_type ppows_begin, ppows_end;
  compos_it_type spows_begin, spows_end;
};

class SeedIterator {
public:
  typedef std::ptrdiff_t difference_type;
  typedef std::vector<pow_type> value_type;
  typedef const std::vector<pow_type> *pointer;
  typedef const std::vector<pow_type> &reference;
  typedef std::input_iterator_tag iterator_category;
  SeedIterator();
  SeedIterator(const SeedIterator &) = default;
  SeedIterator(SeedIterator &&) = default;
  SeedIterator(const std::vector<bool> &, const compos_rit_type &,
               const compos_rit_type &, const compos_it_type &,
               const compos_it_type &);
  SeedIterator &operator=(const SeedIterator &) = default;
  SeedIterator &operator++();
  SeedIterator operator++(int);
  std::vector<pow_type> operator*() const;
  std::unique_ptr<std::vector<pow_type>> operator->() const;
  friend bool operator==(const SeedIterator &, const SeedIterator &);
  friend void swap(SeedIterator &, SeedIterator &);

private:
  std::vector<bool> sector;
  compos_rit_type p_begin, p_it, p_end;
  compos_it_type s_it, s_end;
};

inline SeedIterator Seeds::begin() const {
  return SeedIterator(sector, ppows_begin, ppows_end, spows_begin, spows_end);
}

inline SeedIterator Seeds::end() const { return SeedIterator(); }

inline SeedIterator Seeds::cbegin() const { return begin(); }

inline SeedIterator Seeds::cend() const { return end(); }

bool operator==(const SeedIterator &, const SeedIterator &);
bool operator!=(const SeedIterator &, const SeedIterator &);
void swap(SeedIterator &, SeedIterator &);

/******************
 * Integral class *
 ******************/

class IntegralRelations;
class Topology;

class IntegralProperties {
public:
  uint32_t topology;
  uint32_t sector;
  uint32_t lines;
  uint32_t nums;
  uint32_t dots;
  uint32_t sps;
  uint32_t ppows_weight;
  uint32_t spows_weight;
};

class Integral {
  friend class IntegralRelations;
  friend class Topology;
  friend bool operator==(const Integral &, const Integral &);
  friend bool operator!=(const Integral &, const Integral &);
  struct preferred_basis_t {
    std::vector<Integral> igls;
    // map (topoid,sector) -> vector of (Integral,string) equations
    std::unordered_map<std::pair<uint32_t, uint32_t>, std::vector<is_equation>>
        lcs;
    std::size_t n_basis_lcs;
  };
  struct external_eqs_map_t {
    // map (topoid,sector) -> vector of (Integral,string) equations
    std::unordered_map<std::pair<uint32_t, uint32_t>, std::vector<is_equation>>
        lcs;
  };
  struct external_eqs_t {
    // The number of imported form factors.
    // Used to continuously number the imported form factors.
    std::size_t nimported_ffs;
    // eqs: Both form factors and other external equations.
    std::vector<is_equation> eqs;
  };

public:
  constexpr static uint32_t no_topoid{std::numeric_limits<uint32_t>::max()};
  enum class StringFormat { pyred, mathematica, form, indices };
  // constructors
  Integral() = default;
  Integral(const Integral &) = default;
  Integral(Integral &&) = default;
  Integral(uint32_t topoid, const std::vector<pow_type> &);
  Integral(uint32_t topoid, std::vector<pow_type> &&);
  Integral(std::istringstream &, uint32_t topoid = no_topoid,
           bool force_topoid = false);
  Integral(const std::string &s) {
    std::istringstream ss{s};
    *this = Integral(ss);
  }
  Integral(Weight);
  // members
  uint32_t m_topoid;
  std::vector<pow_type> m_powers;
  // methods
  Weight to_weight() const;
  void reserve(std::size_t);
  std::size_t size() const;
  Integral zip(const std::vector<pow_type> &) const;
  Integral &operator=(const Integral &) = default;
  bool operator<(const Integral &) const;
  std::string to_string(const StringFormat & = StringFormat::pyred) const;
  // static methods
  static std::vector<std::shared_ptr<Topology>> &setup(
      uint32_t sector_ordering = 0, uint32_t dotsp_ordering = 0,
      const std::string &configdir = "./config", const std::string &ibpdir = "",
      const std::string &basisfile = "",
      const std::string &extraRelationsFile = "", const bool use_li = true);
  // Add form factors: return the list of Integral instances representing
  // the form factors in the order of their submission.
  static std::vector<Integral> add_formfactors(const std::string &ff_file);
  static std::vector<Integral> add_formfactors(
      const std::vector<is_equation> &ffs);
  static void add_equations(const std::string &eqs_file);
  static void add_equations(std::vector<is_equation> &&eqs);
  //   static void add_equations();
  static std::vector<uint32_t> assign_weight_bits(
      const std::vector<uint32_t> & = {});
  static std::vector<uint32_t> assign_weight_bits(
    uint32_t maxdots, uint32_t maxsps,
    uint32_t maxdotcombs, uint32_t maxspcombs);
  static std::vector<uint32_t> assign_weight_bits(
      const std::vector<SeedSpec> &);
  static uint32_t np() { return s_np; }
  static void np(uint32_t newnp);
  static void use_cache(int = 3);
  static int get_cache_level() { return s_cache_level; };
  static void clear_cache(int = 3);
  static void auto_clear_cache(int = 3);
  static uint32_t sector_to_sectweight(uint32_t topo, uint32_t s) {
    return s_sector_weight.at(topo).at(s).first;
  }
  static uint32_t sector_from_sectweight(uint32_t topo, uint32_t w) {
    return s_sector_weight.at(topo).at(w).second;
  }
  static bool is_zero(Weight);
  static IntegralProperties properties(Weight);
  static Seeds seeds(uint32_t, uint32_t, uint32_t);
  static std::string ibpdir() { return s_ibpdir; }
  static bool get_use_li() { return s_use_li; }
  static uint32_t sector_ordering() { return s_sector_ordering; }
  static uint32_t dotsp_ordering() { return s_dotsp_ordering; }
  static uint32_t parse_sector(std::string, uint32_t = s_np);
  static std::string bigendian_sector(uint32_t sector, uint32_t np);
  static const preferred_basis_t &get_preferred_masters();
  static const std::vector<is_equation> &get_preferred_masters(uint32_t topoid,
                                                               uint32_t sector);
  static const external_eqs_t &get_external_equations();
  static const external_eqs_map_t &get_external_equations_map();
  // static members
  // default weight bits (see static member initialisation):
  // topology | sector | dotsp1 | dotsp2 | ppows | spows || sum
  //     rest |   s_np |      5 |      4 |    17 |    13 ||  64 (or 128)
  static std::array<uint32_t, 4> s_default_weight_bits;
  static std::array<uint32_t, 6> s_weight_bits;
  static std::array<uint32_t, 6> s_weight_proj;

private:
  // static methods
  static std::pair<uint32_t, uint32_t> weight_to_dots(uint32_t, uint32_t,
                                                      uint32_t);
  static void weight_overflow(uint32_t, uint32_t, uint32_t, uint32_t, uint32_t,
                              uint32_t);
  static void set_preferred_masters(bool /*weight overflow is fatal*/);
  // static members
  // s_np is automatically initialised in the Topology or IntegralRelations
  // constructor (whichever is called first).
  static uint32_t s_np;
  static uint32_t s_sector_ordering;
  static uint32_t s_dotsp_ordering;
  static std::string s_ibpdir;
  static uint32_t s_minlines;
  // If s_auto_topologies, Integral(stringstream) will automatically define
  // a new topology when a new name is encountered.
  // When a topology is created by other means,
  // s_auto_topologies will be set to false.
  static bool s_auto_topologies;
  // s_trivialsectors and s_sector_weight per topology in vectors
  // are duplicated here (wrt. the Topology objects) for faster access
  // than via Topology::id_to_topo(m_topoid)
  static std::vector<std::unordered_set<uint32_t>> s_trivialsectors;
  static bool s_nontrivial_zerosector;
  static std::vector<std::vector<std::pair<uint32_t, uint32_t>>>
      s_sector_weight;
  // s_cache_level: 0=(no cache), 1=(weight to integral only),
  //                2=(integral to weight only), 3=(both)
  static int s_cache_level;
  // s_auto_clear_cache_level: effectively used as argument for clear_cache()
  //                           when clear_cache(0) is called.
  static int s_auto_clear_cache_level;
  // The list of preferred master integrals
  static preferred_basis_t s_preferred_masters;

  // External equations: equations that will be added to the generated
  // systems of equations. E.g. equations defining form factors.
  static external_eqs_t s_external_equations;
  static external_eqs_map_t s_external_equations_map;
  // Map preferred master integrals to their custom weight.
  static Cache<Integral, Weight> s_preferred_masters_i2w;
  // Map custom weights of preferred master integrals to their default weight.
  static std::vector<Weight> s_weight_map;
  thread_local static Cache<Integral, Weight> s_i2w_cache;
  thread_local static Cache<Weight, Integral> s_w2i_cache;
  static bool s_use_li; // use LI identities?
};

inline void Integral::reserve(std::size_t sz) { m_powers.reserve(sz); }

inline std::size_t Integral::size() const { return m_powers.size(); }

std::ostream &operator<<(std::ostream &, const Integral &);

/******************
 * Topology class *
 ******************/

// defined in relations.cpp
void delete_IntegralRelations(IntegralRelations *);

using relations_ptr = std::unique_ptr<IntegralRelations,
                                      std::function<void(IntegralRelations *)>>;

// Factory to be called in the Topology constructor.
relations_ptr import_relations(const std::string &toponame, uint32_t np,
                               uint32_t topoid);
// For a consistency check of the topology ID when the relations are attached
// to a Topology in the Topology constructor.
uint32_t get_relations_topoid(const relations_ptr &);

class TopoConfigData {
  // Parser and data container for information
  // in kinematics.yaml and integralfamilies.yaml.
private:
  class topoconfig_error : public std::exception {
  private:
    std::string msg;
  public:
    inline topoconfig_error(const std::string &s) : msg(s) {}
    virtual inline const char *what() const noexcept { return msg.c_str(); }
  };
  void check_node(const YAML::Node &node, const std::string &name);
public:
  class IntegralFamily {
  public:
    std::string m_name;
    uint32_t m_np;
    std::vector<uint32_t> m_topsectors;
    PropagatorPermutation m_permutation;
  };
  std::vector<std::string> m_invariants;
  std::string m_settoone;
  std::vector<IntegralFamily> m_integralfamilies;
  TopoConfigData() {}
  TopoConfigData(const std::string &project_dir);
  static std::vector<Integral> import_integrals(
    const std::string &filename, uint32_t topoid = Integral::no_topoid);
  static std::pair<std::vector<Integral>, std::vector<is_equation>>
  import_basis(const std::string &filename,
               uint32_t topoid = Integral::no_topoid);
};

class Topology {
  // Topology has no copy constructor/assignment.
  // These could be defined under the following conditions:
  // All non-static public members of this class must be const.
  // Otherwise a copy of a topology may be modified
  // and therefore not agree with the one stored in s_id_to_topo.
  // NOTE: Integrals must only be created after all topologies
  //       have been created. The reason for this is that Integral::s_minlines
  //       is static and global for all topologies, i.e. it might be reduced
  //       if new topologies are defined.
  //       Exception: it's ok if minlines is guaranteed to stay the same,
  //                  e.g. when topologies are defined without trivial sector
  //                  information (like with auto topology definition).
  friend class Integral;
  friend std::shared_ptr<Topology> new_topology(
      const std::string &name, const uint32_t np,
      const std::vector<uint32_t> &topsects,
      const std::vector<uint32_t> &trivialsectors,
      relations_ptr &&relations,
      const PropagatorPermutation &permutation);
  friend std::shared_ptr<Topology> new_topology(
      const std::string &name, const uint32_t np,
      const std::vector<uint32_t> &topsects,
      const PropagatorPermutation &permutation);
  friend std::ostream &operator<<(std::ostream &, const Integral &);

public:
  Topology() = delete; // implicitly deleted anyway
  Topology(const Topology &) = delete; // implicitly deleted anyway
  // members
  const std::string m_name;
  const uint32_t m_np;
  const uint32_t m_id;
  const std::vector<uint32_t> m_topsectors;
  const std::unordered_set<uint32_t> m_trivialsectors;
  const PropagatorPermutation m_permutation;
  // methods
  Integral operator()(const std::vector<pow_type> &) const;
  Integral operator()(std::vector<pow_type> &&) const;
  std::string name() const { return m_name; }
  SeedSpec seed_spec(uint32_t sector, uint32_t rmax, uint32_t smax,
                     int dmax = -1, int n1 = -1,
                     SeedSpec::Recursive rec = SeedSpec::Recursive::full) const;
  const IntegralRelations &relations() const { return *m_relations; }
  // static methods
  static std::shared_ptr<Topology> id_to_topo(uint32_t topoid) {
    return s_id_to_topo.at(topoid);
  }
  static uint32_t toponame_to_id(const std::string &toponame,
                                 bool fatal = true) {
    auto posit = s_toponame_to_id.find(toponame);
    if (posit == s_toponame_to_id.end()) {
      if (fatal) {
        throw input_error(
            std::string("Topology::toponame_to_id(): topology with name ") +
            toponame + " has not been defined.");
      }
      else {
        return Integral::no_topoid;
      }
    }
    return posit->second;
  }
  static const std::vector<std::shared_ptr<Topology>> &get_topologies() {
    return s_id_to_topo;
  }
  static std::vector<std::pair<std::string,uint32_t>> parse_topology_ordering(
    const std::string &);
  static std::vector<uint32_t> import_trivialsectors(const std::string &);
  static std::vector<std::pair<uint32_t,uint32_t>> sector_weight_table(
    uint32_t sector_ordering,
    const std::vector<uint32_t> &topsectors,
    const PropagatorPermutation &permutation);

private:
  Topology(const std::string &name, uint32_t np,
           const std::vector<uint32_t> &topsects,
           const std::vector<uint32_t> &trivialsectors,
           relations_ptr &&relations,
           const PropagatorPermutation &permutation);
  // methods
  void assign_sector_weight();
  uint32_t sector_to_sectweight(uint32_t s) { return m_sector_weight[s].first; }
  uint32_t sector_from_sectweight(uint32_t w) {
    return m_sector_weight[w].second;
  }
  // members
  std::vector<std::pair<uint32_t, uint32_t>> m_sector_weight;
  uint32_t m_minlines;
  // Cannot use std::default_delete here, because then the full declaration
  // of IntegralRelations is needed in the Topology constructor
  // when m_relations is created.
  // The actual deleter is implemented in relations.cpp.
  relations_ptr m_relations;
  // static members
  static uint32_t s_next_id;
  // use a pointer here so that references don't get invalidated
  static std::vector<std::shared_ptr<Topology>> s_id_to_topo;
  static std::unordered_map<std::string, uint32_t> s_toponame_to_id;
  static std::string s_basislc_name;
  static std::string s_formfactor_name;
  static std::string s_userweights_name;
};

inline Integral Topology::operator()(const std::vector<pow_type> &pows) const {
  return Integral(m_id, pows);
}

inline Integral Topology::operator()(std::vector<pow_type> &&pows) const {
  return Integral(m_id, std::move(pows));
}

// Topologies must only be created by these factories
// so that they get properly registered.
std::shared_ptr<Topology> new_topology(
    const std::string &name, const uint32_t np,
    const std::vector<uint32_t> &topsects,
    const std::vector<uint32_t> &trivialsectors,
    relations_ptr &&relations = nullptr,
    const PropagatorPermutation &permutation = {});
std::shared_ptr<Topology> new_topology(const std::string &name,
                                       const uint32_t np,
                                       const std::vector<uint32_t> &topsects,
                                       const PropagatorPermutation &permutation = {});

} // namespace pyred

namespace std {
template <>
struct hash<pyred::Integral> {
  size_t operator()(const pyred::Integral &igl) const {
    size_t h = igl.m_topoid + 0x9e3779b9;
    for (const auto &pow : igl.m_powers) {
      h ^= pow + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};
} // namespace std

#endif
