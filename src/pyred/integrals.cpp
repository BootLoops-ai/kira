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

#include <algorithm> // min, max, swap
#include <array>
#include <cctype> // isdigit
#include <cstddef>
#include <fstream>
#include <numeric> // accumulate
#include <tuple>
#include <unordered_map>
#include <utility>

#include "pyred/coeff_helper.h" // set invariants
#include "pyred/integrals.h"
#include "pyred/parallel.h"
#include "pyred/parser.h" // split
#include "yaml-cpp/yaml.h"

namespace pyred {

std::pair<Weight, std::string> split_intcoeff(const std::string& s) {
  static bool first = true;
  auto spair = split(s, '*', 1);
  if (spair.size() != 2) {
    throw parser_error(std::string("ill-formed integral*coefficient: ") + s);
  }
  Weight weight;
  std::string rest;
  std::istringstream ss{spair.front()};
  // Try if the integral is given as integer weight.
  auto success = bool(ss >> weight);
  if (!success) {
    // Assume that the integral is given as TOPO[indices].
    // First reset the error flag of the stream.
    ss.clear();
    weight = Integral(ss).to_weight();
  }
  else if (first) {
    first = false;
    // If this is a user-defined weight, make sure the corresponding
    // pseudo-topology gets initialised.
    Integral{weight};
  }
  pyred::getline_normalized(ss, rest);
  if (rest.find_first_not_of(' ') != std::string::npos) {
    // The remainder must be empty or only contain whitespace.
    throw parser_error(std::string("ill-formed integral*coefficient: ") + s);
  }
  return {weight, spair.back()};
}

/*****************
 * Combinatorics *
 *****************/

std::size_t bitwidth(std::size_t n) {
  // Count the number of bits which are required to store n.
  std::size_t bits{0};
  while (n) {
    ++bits;
    n >>= 1;
  }
  return bits;
}

std::size_t count_set_bits(std::size_t n) {
  // Count the number of assigned bits of n.
  std::size_t bits{0};
  while (n) {
    bits += n & 1;
    n >>= 1;
  }
  return bits;
}

std::size_t factorial(const std::size_t n) {
  std::size_t f{1};
  for (std::size_t k = 0; k != n;) {
    f *= ++k;
  }
  return f;
}

std::size_t binomial(int n, int k) {
  // Binomial coefficient 'n choose k'.
  if (k < 0 || k > n) return 0;
  int bino = 1;
  for (int i = 0; i != std::min(k, n - k); ++i) {
    bino = (bino * (n - i)) / (i + 1);
  }
  return bino;
}

std::size_t dotsp_combs(const std::size_t np, const std::size_t lines,
                        const std::size_t dots, const std::size_t sps) {
  // The number of combinations of 'dots' dots on 'lines' lines
  // and 'sps' scalar products on 'np-lines' numerators.
  if (lines == np) {
    return binomial(lines + dots - 1, dots);
  }
  return binomial(lines + dots - 1, dots) * binomial(np - lines + sps - 1, sps);
}


PropagatorPermutation::PropagatorPermutation(
    std::vector<std::size_t> perm, std::size_t n)
: m_permutation{perm}, m_n{n}
{
  // n is optional. If given, check that n = perm.size().
  // Check that perm is a permutation of 0..n-1 or 1..n.
  // In the latter case, subtract 1 from all elements in the result.
  // The identity permutation is replaced by {}.
  // Note that length information is then only available if n has been set
  // correctly (this cannot be checked).
  if (m_n == std::numeric_limits<std::size_t>::max() && perm.size()) {
    m_n = perm.size();
  }
  if (perm.empty()) {
    // An empty permutation is treated as the identity.
    m_permutation.clear();
    return;
  }
  if (perm.size() && m_n != perm.size()) {
    throw permutation_error("PropagatorPermutation: wrong permutation size.");
  }
  // Check that perm is actually a permutation.
  std::sort(perm.begin(), perm.end());
  auto minelem = perm.front();
  for (std::size_t k = 0u; k != m_n; ++k) {
    if (perm[k] != k + minelem) {
      throw permutation_error("PropagatorPermutation: is not a permutation.");
    }
  }
  // Shift all indices by -1 if we got a permutation of 1..n.
  if (minelem == 1u) {
    for (auto &el: m_permutation) {
      el -= 1;
    }
  }
  else if (minelem != 0u) {
    throw permutation_error("PropagatorPermutation: wrong minimal element.");
  }
  // Set m_permutation={} for the identity.
  bool is_identity = true;
  for (std::size_t k = 0u; k != m_permutation.size(); ++k) {
    if (m_permutation[k] != k) {
      is_identity = false;
      break;
    }
  }
  if (is_identity) {
    m_permutation.clear();
  }
}


uint32_t PropagatorPermutation::permute_sector(uint32_t sector,
                                               bool inverse) const {
  if (m_n == std::numeric_limits<std::size_t>::max()) {
    throw permutation_error("permute_sector(): unknown permutation size.");
  }
  if (m_permutation.empty()) {
    return sector;
  }
  uint32_t psect = 0u;
  if (inverse) {
    // In: a sector refering to the permuted list of propagators.
    // Out: the corresponding sector wrt. the unpermuted propagators.
    for (std::size_t pos = 0u; pos != m_permutation.size(); ++pos) {
      if (sector & (1 << pos)) {
        psect += (1 << m_permutation[pos]);
      }
    }
  }
  else {
    // In: a sector refering to the unpermuted list of propagators.
    // Out: the corresponding sector wrt. the permuted propagators.
    for (std::size_t pos = 0u; pos != m_permutation.size(); ++pos) {
      if (sector & (1 << m_permutation[pos])) {
        psect += (1 << pos);
      }
    }
  }
  return psect;
}


/****************
 * Compositions *
 ****************/

std::vector<std::vector<uint32_t>> compositions(uint32_t n, uint32_t k) {
  // Return a list of all compositions of the integer n into k parts.
  std::vector<std::vector<uint32_t>> compos;
  std::vector<std::vector<uint32_t>> tmp;
  compos.reserve(binomial(n + k - 1, k - 1));
  tmp.reserve(compos.size());
  compos.emplace_back();
  // stop condition i<k (not i!=k) to cover the case k=0.
  for (uint32_t i = 1; i < k; ++i) {
    for (const auto& c : compos) {
      auto clen =
          std::accumulate(c.cbegin(), c.cend(), static_cast<uint32_t>(0));
      for (uint32_t j = 0; j != n + 1 - clen; ++j) {
        tmp.emplace_back(c);
        tmp.back().emplace_back(j);
      }
    }
    compos.swap(tmp);
    tmp.clear();
  }
  if (k > 0) {
    for (auto& c : compos) {
      c.emplace_back(n - std::accumulate(c.cbegin(), c.cend(), 0));
    }
  }
  return compos;
}

// Cache for get_composition_id().
// Note that pointers/references to values in an unordered_map
// remain valid after insert/erase/rehashing.
#ifndef PYRED_GLOBAL_COMPOSITION_CACHE
thread_local
#endif
    std::unordered_map<std::pair<uint32_t, uint32_t>,
                       std::unique_ptr<std::pair<
                           std::unordered_map<std::vector<uint32_t>, uint32_t>,
                           std::vector<std::vector<uint32_t>>>>>
        composition_id;

const std::pair<std::unordered_map<std::vector<uint32_t>, uint32_t>,
                std::vector<std::vector<uint32_t>>>&
get_composition_id(uint32_t dots, uint32_t lines) {
  // Return the composition weight dictionary 'compo_weights'
  // and the composition from weight list 'compo_from_weight'.
  // compo_weights[compo] = weight and compo_from_weight[weight] = compo
  // for a composition of the integer 'dots' into 'lines' parts.
#ifdef PYRED_GLOBAL_COMPOSITION_CACHE
  static std::mutex mtx;
#endif
  {
#ifdef PYRED_GLOBAL_COMPOSITION_CACHE
    std::lock_guard<std::mutex> lock(mtx);
#endif
    auto resit = composition_id.find({dots, lines});
    if (resit != composition_id.cend()) {
      return *resit->second;
    }
  }
  std::unordered_map<std::vector<uint32_t>, uint32_t> compo_weights;
  std::vector<std::vector<uint32_t>> compo_from_weight;
  uint32_t weight{0};
  auto compos = compositions(dots, lines);
  for (const auto& compo : compos) {
    compo_weights.insert({compo, weight});
    compo_from_weight.push_back(compo);
    ++weight;
  }
#ifdef PYRED_GLOBAL_COMPOSITION_CACHE
  std::lock_guard<std::mutex> lock(mtx);
#endif
  // Check if the result has been added to the cache in the meantime.
  // If yes, return the result from the cache to not invalidate
  // references to it.
  auto resit = composition_id.find({dots, lines});
  if (resit != composition_id.cend()) {
    return *resit->second;
  }
  // Otherwise add it to the cache and return.
  auto it =
      composition_id
          .emplace(
              std::make_pair(dots, lines),
              std::make_unique<
                  std::pair<std::unordered_map<std::vector<uint32_t>, uint32_t>,
                            std::vector<std::vector<uint32_t>>>>(
                  compo_weights, compo_from_weight))
          .first;
  return *it->second;
}

/******************
 * Seed generator *
 ******************/

SeedSpec::SeedSpec(const Topology& topo, uint32_t sector, uint32_t rmax,
                   uint32_t sps, int maxdots, int n1, Recursive rec)
    // Make sure m_topo refers to a persistent object, not to a temporary copy.
    : m_topo(Topology::id_to_topo(topo.m_id)),
      m_sector(sector),
      m_sps(sps),
      m_maxdots(maxdots),
      m_recursive(rec) {
  m_lines = 0;
  while (sector) {
    if (sector & 1) ++m_lines;
    sector >>= 1;
  }
  if (rmax < m_lines) {
    throw input_error("Seed specification: rmax must not be lower than "
                      "the number of lines of the sector.");
  }
  m_mindots = rmax - m_lines;
  if (n1 < 0) {
    m_n1 = m_lines;
  }
  else {
    m_n1 = std::min(static_cast<uint32_t>(n1), m_lines);
  }
  if (maxdots < 0) {
    m_maxdots = std::numeric_limits<uint32_t>::max();
  }
}

std::vector<SeedSpecSector> SeedSpec::expand_sector(
    const std::vector<SeedSpec>& seedspecs,
    const std::vector<SeedSpec>& seedcompl, bool include_trivial) {
  // TODO: include target sectors of sector mappings with appropriate corners.
  std::vector<SeedSpecSector> sectorspecs;
  // The seed specifications might belong to different topologies.
  // Treat one topology after another.
  std::unordered_set<const Topology*> topos;
  for (const auto& spec : seedspecs) {
    if (spec.m_recursive == Recursive::no) {
      throw input_error("Non-recursive seeds in dots/sps are only allowed in "
                        "seed complements.");
    }
    topos.insert(spec.m_topo.get());
  }
  for (const auto topo : topos) {
    const uint32_t np = topo->m_np;
    // For this topology: map sector_number -> {number_of_lines, spec}
    // for all sectors of the topology
    // (excluding trivial sectors unless include_trivial=true).
    std::unordered_map<uint32_t, std::pair<uint32_t, SeedSpecSector>>
        spec4sector;
    const auto& zs = topo->m_trivialsectors;
    for (uint32_t sn = 1; sn != (static_cast<uint32_t>(1) << np); ++sn) {
      if (zs.find(sn) == zs.end() || include_trivial) {
        auto sncpy = sn;
        uint32_t lines{0};
        while (sncpy) {
          if (sncpy & 1) ++lines;
          sncpy >>= 1;
        }
        spec4sector.insert({sn, {lines, SeedSpecSector(*topo, sn)}});
      }
    }
    for (const auto& spec : seedspecs) {
      // for all seed specifications for the current topology
      if (spec.m_topo->m_id != topo->m_id) continue;
      // For all sub-sectors, add the (dots, sps) corners
      // to the SeedSpecSector object.
      for (auto& sn_lines_spec : spec4sector) {
        const auto sn = sn_lines_spec.first;
        // Skip if sn is not a sub-sector of spec.m_sector.
        if ((spec.m_sector & sn) != sn) {
          continue;
        }
        // If the seed spec is non-recursive in its sector (but recursive
        // in dots/sps -- fully non-recursive is not allowed here),
        // keep only this sector
        if (spec.m_recursive == Recursive::dotsp && spec.m_sector != sn) {
          continue;
        }
        const auto lines = sn_lines_spec.second.first;
        auto& sector_spec = sn_lines_spec.second.second;
        auto dots = spec.m_mindots;
//        if (lines < spec.m_n1) {
          dots = std::min(dots + spec.m_n1 - lines, spec.m_maxdots);
//        }
//        std::cout << "dots: " << spec.m_mindots << " " << spec.m_maxdots << " " << dots << std::endl;
        sector_spec.add_corner(dots, spec.m_sps);
      }
    }
    for (const auto& spec : seedcompl) {
      // for all seed skip specifications for the current topology
      if (spec.m_topo->m_id != topo->m_id) continue;
      // For all sectors: If the sector is
      // * a subsector of a sector recursive seed complement,
      // * or equal to the seed complement sector,
      // add a skip prescription to the sector.
      for (auto& sn_lines_spec : spec4sector) {
        const auto sn = sn_lines_spec.first;
        // Skip if sn is not a sub-sector of spec.m_sector, or
        // if the complement is sector non-recursive
        // and the sectors do not agree.
        if (((spec.m_sector & sn) != sn) ||
            ((spec.m_recursive != Recursive::full) && (spec.m_sector != sn))) {
          continue;
        }
        const auto lines = sn_lines_spec.second.first;
        auto& sector_spec = sn_lines_spec.second.second;
        auto dots = spec.m_mindots;
//        if (lines < spec.m_n1) {
          dots = std::min(dots + spec.m_n1 - lines, spec.m_maxdots);
//        }
//        std::cout << "compl dots: " << dots << std::endl;
        sector_spec.skip(dots, spec.m_sps, spec.m_recursive != Recursive::no);
      }
    }
    for (const auto& sn_lines_spec : spec4sector) {
      if (!sn_lines_spec.second.second.empty()) {
        sectorspecs.push_back(sn_lines_spec.second.second);
      }
    }
  }
  std::sort(sectorspecs.begin(), sectorspecs.end());
  return sectorspecs;
}

std::vector<Integral> SeedSpec::list_integrals(
    const std::vector<SeedSpec>& seedspecs, int parallel) {
  // Generate all integrals in the given seed specification.
  // Include integrals in trivial sectors (except sector 0).
  // parallel: per default (parallel=-1) use Config::parallel().
  std::vector<Integral> selected_integrals;
  auto sector_seeds = expand_sector(seedspecs, {}, true);
  if (parallel < 0) parallel = Config::parallel();
  using seed_t = decltype(sector_seeds)::value_type;
  selected_integrals =
      Distributor<seed_t, std::vector<Integral>, std::vector<Integral>>::create(
          [](const seed_t& spec4sect) {
            std::vector<Integral> sel_igls;
            for (const auto& dotsp : spec4sect.expand()) {
              for (auto seed : Seeds(spec4sect.m_topo->m_np, spec4sect.m_sector,
                                     dotsp.first, dotsp.second)) {
                sel_igls.push_back((*spec4sect.m_topo)(std::move(seed)));
              }
            }
            return sel_igls;
          },
          sector_seeds,
          [](std::vector<Integral>&& in_igls, std::vector<Integral>& acc_igls) {
            for (auto& igl : in_igls) {
              acc_igls.push_back(igl);
            }
          },
          parallel);
  return selected_integrals;
}

std::vector<Weight> SeedSpec::integral_selector(
    const std::vector<SeedSpec>& seedspecs, int parallel) {
  // Generate the weights of all integrals in the given seed specification.
  // Integrals with Weight::none() (i.e. in trivial sectors) are discarded.
  // Disable the integral/weight cache.
  // parallel: per default (parallel=-1) use Config::parallel().
  auto cache_level_bak = Integral::get_cache_level();
  Integral::use_cache(0);
  std::vector<Weight> selected_weights;
  auto sector_seeds = expand_sector(seedspecs);
  if (parallel < 0) parallel = Config::parallel();
  using seed_t = decltype(sector_seeds)::value_type;
  selected_weights =
      Distributor<seed_t, std::vector<Weight>, std::vector<Weight>>::
          create(
              [](const seed_t& spec4sect) {
                std::vector<Weight> sel_weights;
                for (const auto& dotsp : spec4sect.expand()) {
                  for (auto seed :
                       Seeds(spec4sect.m_topo->m_np, spec4sect.m_sector,
                             dotsp.first, dotsp.second)) {
                    auto weight =
                        (*spec4sect.m_topo)(std::move(seed)).to_weight();
                    if (weight != Weight::none()) {
                      sel_weights.push_back(weight);
                    }
                  }
                }
                return sel_weights;
              },
              sector_seeds,
              [](std::vector<Weight>&& in_weights,
                 std::vector<Weight>& acc_weights) {
                for (auto& w : in_weights) {
                  acc_weights.push_back(w);
                }
              },
              parallel);
  // Restore the integral/weight cache level.
  Integral::use_cache(cache_level_bak);
  return selected_weights;
}

SeedSpecSector::SeedSpecSector(
    const Topology& topo, uint32_t sector,
    const std::vector<std::pair<uint32_t, uint32_t>>& corners)
    : m_topo(Topology::id_to_topo(topo.m_id)),
      m_sector(sector),
      m_corners(corners) {}

std::vector<std::pair<uint32_t, uint32_t>> SeedSpecSector::expand() const {
  std::vector<std::pair<uint32_t, uint32_t>> dotsps;
  uint32_t dotlimit = 0;
  uint32_t splimit = 0;
  for (const auto& corner : m_corners) {
    dotlimit = std::max(dotlimit, corner.first);
    splimit = std::max(splimit, corner.second);
  }
  // For each (dots,sps) within (dotlimit,splimit),
  // first check if it's within one of the corners,
  // then check if it should be skipped.
  for (uint32_t dots = 0; dots <= dotlimit; ++dots) {
    for (uint32_t sps = 0; sps <= splimit; ++sps) {
      bool accept = false;
      for (const auto& corner : m_corners) {
        if (dots <= corner.first && sps <= corner.second) {
          accept = true;
          break;
        }
      }
      if (accept) {
        for (const auto& sk : m_skip) {
          if (std::get<2>(sk)) {
            // recursive skip
            if (dots <= std::get<0>(sk) && sps <= std::get<1>(sk)) {
              accept = false;
              break;
            }
          }
          else {
            if (dots == std::get<0>(sk) && sps == std::get<1>(sk)) {
              accept = false;
              break;
            }
          }
        }
      }
      if (accept) {
        dotsps.push_back({dots, sps});
      }
    }
  }
  return dotsps;
}

void SeedSpecSector::cut() {
  if (Config::symlimits().first < 0 || Config::symlimits().second < 0) {
    // do not symmetrise
    m_corners.clear();
    return;
  }
  std::vector<std::pair<uint32_t, uint32_t>> cutcorners;
  cutcorners.reserve(m_corners.size());
  for (auto& corner : m_corners) {
    cutcorners.push_back(
        {std::min(corner.first,
                  static_cast<uint32_t>(Config::symlimits().first)),
         std::min(corner.second,
                  static_cast<uint32_t>(Config::symlimits().second))});
  }
  m_corners = std::move(cutcorners);
}

bool operator<(const SeedSpecSector& a, const SeedSpecSector& b) {
  // TODO: Use sector ordering for comparison.
  //       Note that the sector ordering is already taken into account in
  //       GeneratorHelper::generate_ts_bunches().
  if (a.m_topo->m_id < b.m_topo->m_id) {
    return true;
  }
  else if (a.m_topo->m_id > b.m_topo->m_id) {
    return false;
  }
  return (a.m_sector < b.m_sector);
}

Seeds::Seeds(uint32_t np, uint32_t sect, uint32_t dots, uint32_t sps)
    : sector(np, false) {
  uint32_t lines{0};
  uint32_t nums;
  for (auto pr = sector.begin(); pr != sector.end(); ++pr) {
    if (sect & 1) {
      ++lines;
      *pr = true;
    }
    sect >>= 1;
  }
  nums = np - lines;
  if (nums || !sps) {
    const auto& ppows = get_composition_id(dots, lines).second;
    const auto& spows = get_composition_id(sps, nums).second;
    ppows_begin = ppows.crbegin();
    ppows_end = ppows.crend();
    spows_begin = spows.cbegin();
    spows_end = spows.cend();
  }
  else {
    // nums=0 && sps>0: there are no such seeds --> set *this to end().
    ppows_begin = compos_rit_type{};
    ppows_end = compos_rit_type{};
    spows_begin = compos_it_type{};
    spows_end = compos_it_type{};
  }
}

SeedIterator::SeedIterator()
    : sector{}, p_begin{}, p_it{}, p_end{}, s_it{}, s_end{} {}

SeedIterator::SeedIterator(const std::vector<bool>& sector,
                           const compos_rit_type& p_begin,
                           const compos_rit_type& p_end,
                           const compos_it_type& s_begin,
                           const compos_it_type& s_end)
    : sector{sector},
      p_begin{p_begin},
      p_it{p_begin},
      p_end{p_end},
      s_it{s_begin},
      s_end{s_end} {}

SeedIterator& SeedIterator::operator++() {
  if (p_it != p_end) {
    ++p_it;
  }
  if (p_it == p_end) {
    p_it = p_begin;
    ++s_it;
  }
  if (s_it == s_end) {
    // The iterator is depleted.
    p_begin = compos_rit_type{};
    p_it = compos_rit_type{};
    p_end = compos_rit_type{};
    s_it = compos_it_type{};
    s_end = compos_it_type{};
  }
  return *this;
}

SeedIterator SeedIterator::operator++(int) {
  auto it = *this;
  ++(*this);
  return it;
}

std::vector<pow_type> SeedIterator::operator*() const {
  std::vector<pow_type> powers;
  powers.reserve(sector.size());
  auto ppows_it = p_it->cbegin();
  auto spows_it = s_it->cbegin();
  for (bool pr : sector) {
    if (pr) {
      powers.push_back(1 + (*(ppows_it++)));
    }
    else {
      powers.push_back(-(*(spows_it++)));
    }
  }
  return powers;
}

std::unique_ptr<std::vector<pow_type>> SeedIterator::operator->() const {
  auto tmp = std::make_unique<std::vector<pow_type>>();
  *tmp = *(*this);
  return tmp;
}

bool operator==(const SeedIterator& lhs, const SeedIterator& rhs) {
  // Empty sectors compare equal to all sectors. This ensures the equality
  // of an initially empty seed iterator to the default constructed iterator.
  return (
      (lhs.sector.empty() || rhs.sector.empty() || lhs.sector == rhs.sector) &&
      lhs.p_begin == rhs.p_begin && lhs.p_it == rhs.p_it &&
      lhs.p_end == rhs.p_end && lhs.s_it == rhs.s_it && lhs.s_end == lhs.s_end);
}

bool operator!=(const SeedIterator& lhs, const SeedIterator& rhs) {
  return !(lhs == rhs);
}

void swap(SeedIterator& lhs, SeedIterator& rhs) {
  std::swap(lhs.sector, rhs.sector);
  std::swap(lhs.p_begin, rhs.p_begin);
  std::swap(lhs.p_it, rhs.p_it);
  std::swap(lhs.p_end, rhs.p_end);
  std::swap(lhs.s_it, rhs.s_it);
  std::swap(lhs.s_end, lhs.s_end);
}

/******************
 * Integral class *
 ******************/

// static members
uint32_t Integral::s_np{0};
uint32_t Integral::s_sector_ordering{3};
uint32_t Integral::s_dotsp_ordering{1};
std::string Integral::s_ibpdir{};
bool Integral::s_use_li{true};
std::array<uint32_t, 4> Integral::s_default_weight_bits{{5, 4, 17, 13}};
std::array<uint32_t, 6> Integral::s_weight_bits{{0, 0, 5, 4, 17, 13}};
std::array<uint32_t, 6> Integral::s_weight_proj{{0, 0, 0, 0, 0, 0}};
uint32_t Integral::s_minlines{1u << 31};
bool Integral::s_auto_topologies{true};
int Integral::s_cache_level{0};
int Integral::s_auto_clear_cache_level{0};
Integral::preferred_basis_t Integral::s_preferred_masters{{}, {}, 0u};
Integral::external_eqs_map_t Integral::s_external_equations_map{{}};
Integral::external_eqs_t Integral::s_external_equations{0u, {}};
Cache<Integral, Weight> Integral::s_preferred_masters_i2w{};
std::vector<Weight> Integral::s_weight_map{};
thread_local Cache<Integral, Weight> Integral::s_i2w_cache{};
thread_local Cache<Weight, Integral> Integral::s_w2i_cache{};
std::vector<std::unordered_set<uint32_t>> Integral::s_trivialsectors{};
bool Integral::s_nontrivial_zerosector{false};
std::vector<std::vector<std::pair<uint32_t, uint32_t>>>
    Integral::s_sector_weight{};

Integral::Integral(uint32_t topoid, const std::vector<pow_type>& pows)
    : m_topoid{topoid}, m_powers{pows} {
  if (size() != s_np) {
    throw init_error("Integral(): wrong number of indices");
  }
}

Integral::Integral(uint32_t topoid, std::vector<pow_type>&& pows)
    : m_topoid{topoid}, m_powers{std::move(pows)} {
  if (size() != s_np) {
    throw init_error("Integral(): wrong number of indices");
  }
}

Integral::Integral(std::istringstream& ss, uint32_t topoid, bool force_topoid)
  : m_topoid{topoid} {
  // Input format: "TOPO[n1,n2,..]*"
  // topoid is used as fallback, if TOPO is not part of the string.
  std::string topostr;
  pyred::getline_normalized(ss >> std::ws, topostr, '['); // discard leading whitespace
  bool define_topology = false;
  int igltype{0};
  if (!topostr.empty() && topoid != no_topoid && (!force_topoid)) {
    // If topostr is not empty, ignore topoid, unless force_topoid=true.
    topoid = no_topoid;
  }
  else if (topostr.empty()) {
    Config::log(1) << "WARNING: deprecated integral notation " << ss.str()
                   << " (missing topology name)" << std::endl;
  }
  if (topoid == no_topoid) {
    // Interpret topostr to set topoid.
    if (topostr == Topology::s_basislc_name) {
      igltype = 1;
      m_topoid = Topology::s_next_id;
    }
    else if (topostr == Topology::s_formfactor_name) {
      igltype = 2;
      m_topoid = Topology::s_next_id;
    }
    else if (!topostr.empty() && isdigit(topostr[0])) {
      // Assume that the string represents a (user-defined) weight.
      Weight weight;
      std::istringstream ss{topostr};
      if (!bool(ss >> weight)) {
        throw input_error(
            "Integral(std::istringstream &, -1): invalid weight as argument.");
      }
      *this = Integral(weight);
      return;
    }
    else {
      m_topoid = Topology::toponame_to_id(topostr, false /*non-fatal*/);
      if (m_topoid == no_topoid) {
        if (s_auto_topologies) {
          define_topology = true;
        }
        else {
          throw input_error(
              std::string(
                  "Integral(std::istringstream &, -1): topology with name \"") +
              topostr +
              "\" has not been defined and auto topology definition "
              "has been disabled (usually because topologies have already "
              "been created  by other means).");
        }
      }
    }
  }
  char sep{','};
  pow_type pow;
  reserve(s_np);
  while (sep == ',') {
    if (!bool(ss >> pow) || !bool(ss >> sep)) {
      throw parser_error("Integral parser error (possibly int32_t overflow).");
    }
    m_powers.push_back(pow);
  }
  if (igltype > 0) {
    for (uint32_t k = 1; k < s_np; ++k) {
      m_powers.push_back(0);
    }
    if (m_powers[0] <= 0) {
      throw input_error("Invalid BASISLC or FORMFACTOR string representation.");
    }
    if (igltype == 2) {
      // The integral represents a form factor.
      // m_powers[0] > 0 is guaranteed here.
      if (static_cast<uint32_t>(m_powers[0]) > s_weight_proj[4]) {
        std::ostringstream ss;
        ss << "Integral(): form factor number too large (got " << m_powers[0]
           << ", allowed: " << s_weight_proj[4]
           << " in the current weight representation.";
        throw input_error(ss.str());
      }
      m_powers[0] = -m_powers[0];
    }
    // Otherwise the integral represents a basis linear combination
    // (nothing to do).
  }
  else if (define_topology) {
    // auto topology definition
    m_topoid = new_topology(topostr, /*np*/m_powers.size(),
                            {/*topsectors*/}, {/*trivialsectors*/},
                            /*relations*/nullptr,
                            /*permutation*/{m_powers.size()})->m_id;
    // Reset s_auto_topologies to true,
    // because the Topology constructor set it to false.
    s_auto_topologies = true;
  }
  if (sep != ']') {
    throw init_error("Integral parser: ill-formatted integral");
  }
  if (size() != s_np) {
    throw init_error("Integral parser: wrong number of indices");
  }
}

Integral::Integral(Weight w) {
  // Convert an integral weight to an integral.
  // NOTE: calling this function with a weight that was not generated by
  // Integral::to_weight() with the same topology definitions and weight bits
  // is undefined behaviour. With an invalid weight it could happen
  // that the weight is mapped to the all-zero integral and then cached.
  // Getting the weight of the all-zero integral will then return the invalid
  // weight instead of Weight::none().
  // Or the program may just crash.
  if (!s_np) {
    // Use this as special case of user-defined weights.
    if (s_auto_topologies) {
      // No topologies defined yet.
      m_topoid = new_topology(Topology::s_userweights_name, /*np*/0u,
                              {/*topsectors*/}, {/*trivialsectors*/},
                              /*relations*/nullptr,
                              /*permutation*/{0u})->m_id;
      // new_topology() sets s_auto_topologies=false.
      // If user-defined weights are used, this is the only topology
      // that must ever be created.
    }
    else if (Topology::s_next_id != 1u) {
      throw input_error(
          "Integral(): Invalid topology definitions for user-defined weights.");
    }
    m_topoid = 0u;
    // A static_assert() in defs.h guarantees
    // that sizeof(Weight) is a multiple of sizeof(pow_type).
    constexpr static int pows_for_weight = sizeof(Weight)/sizeof(pow_type);
    m_powers.resize(pows_for_weight);
    *reinterpret_cast<Weight*>(m_powers.data()) = w;
    return;
  }
  s_i2w_cache.set_if_empty(s_preferred_masters_i2w);
  if (w == Weight::none()) {
    throw input_error(
        "Integral(): cannot construct an integral from Weight::none().");
  }
  if (s_cache_level & 1) {
    auto found = s_w2i_cache.lookup(w);
    if (found.second) {
      *this = found.first;
      return;
    }
  }
  auto props = properties(w);
  m_topoid = props.topology;
  m_powers.reserve(s_np);
  if (m_topoid == Topology::s_next_id) {
    for (uint32_t pos = 0; pos != s_np; ++pos) {
      m_powers.push_back(0);
    }
    if (props.sps == 0) {
      // The integral represents a basis linear combination.
      // props.dots=weight=1,2,...
      m_powers[0] = props.dots;
    }
    else {
      // The integral represents a form factor.
      // props.sps=1,2,..., weight=Weight::none()-props.sps
      m_powers[0] = -props.sps;
    }
  }
  else if (s_np == 1) {
    // Special case np=1
    m_powers.push_back(props.sector ? props.dots + 1 : -props.sps);
  }
  else {
    // * dots to the left have lower weight
    // auto ppows_it = get_composition_id(props.dots,props.lines).second.cend()[
    //                   -1-static_cast<int>(props.ppows_weight)].cbegin();
    // * dots to the right have lower weight
    auto ppows_it = get_composition_id(props.dots, props.lines)
                        .second[props.ppows_weight]
                        .cbegin();
    // * sps to the right have lower weight
    // auto spows_it = get_composition_id(props.sps,props.nums).second[
    //                   props.spows_weight].cbegin();
    // * sps to the left have lower weight
    auto spows_it =
        get_composition_id(props.sps, props.nums)
            .second.cend()[-1 - static_cast<int>(props.spows_weight)]
            .cbegin();
    const auto &perm = Topology::id_to_topo(m_topoid)->m_permutation;
    if (perm.is_identity()) {
      for (uint32_t dummy = 0; dummy != s_np; ++dummy) {
        if (props.sector & 1) {
          m_powers.push_back(1 + (*(ppows_it++)));
        }
        else {
          m_powers.push_back(-(*(spows_it++)));
        }
        props.sector >>= 1;
      }
    }
    else {
      m_powers.resize(s_np);
      for (const auto pos : perm) {
        if (props.sector & (1 << pos)) {
          m_powers[pos] = 1 + (*(ppows_it++));
        }
        else {
          m_powers[pos] = -(*(spows_it++));
        }
      }
    }
  }
  // Insert the element into the cache(s)
  // if it hasn't been cached in the meantime.
  if (s_cache_level & 2) {
    s_i2w_cache.insert(*this, w);
  }
  if (s_cache_level & 1) {
    s_w2i_cache.insert(w, *this);
  }
}

bool Integral::is_zero(Weight w) {
  return w == Weight::none() ? true : false;
}

IntegralProperties Integral::properties(Weight w) {
  // Extract topology, sector, number of lines, numerators,
  // dots, scalar products, dot and sp composition weight.
  // Return them as an IntegralProperties instance.
  // * Changes necessary to give dot composition has higher weight
  //   than sp composition:
  //   spow_w = w.project_shift(s_weight_proj[5], s_weight_bits[5]);
  //   ppow_w = w.project_shift(s_weight_proj[4], s_weight_bits[4]);
  // * Here: sp composition has higher weight than dot composition.
  if (!s_np) {
    // Special case: user-defined weights.
    return IntegralProperties{0u, 1u, 0u, 0u, 0u, 0u, 0u, 0u};
  }
  if (w <= Weight(s_preferred_masters.igls.size()
                  + s_preferred_masters.n_basis_lcs)) {
    if (w <= Weight(s_weight_map.size())) {
      // If the weight has been mapped, i.e. it belongs to a preferred
      // master integral, replace the custom weight by the default weight
      // (because the custom weight does not contain integral information).
      w = s_weight_map[w.project<uint32_t>() - 1];
    }
    else {
      // This is the weight of a basis linear combination.
      return IntegralProperties{
        Topology::s_next_id, 1u, 1u, s_np - 1u,
        w.project<uint32_t>()
        - static_cast<uint32_t>(s_preferred_masters.igls.size()),
        0u, 0u, 0u};
    }
  }
  else if (w >= Weight::amplitude(s_weight_proj[4])) {
    // s_weight_proj[4] is used as the maximal possible form factor number.
    // This is the weight of a form factor.
    return IntegralProperties{
      Topology::s_next_id, 1u, 1u, s_np - 1u, 0u,
      s_weight_proj[4] ^ w.project(s_weight_proj[4]), 0u, 0u};
  }
  auto ppow_w = w.project_shift(s_weight_proj[4], s_weight_bits[4]);
  auto spow_w = w.project_shift(s_weight_proj[5], s_weight_bits[5]);
  auto dotsp_weight_2 = w.project_shift(s_weight_proj[3], s_weight_bits[3]);
  auto dotsp_weight_1 = w.project_shift(s_weight_proj[2], s_weight_bits[2]);
  auto wforsect = w.project_shift(s_weight_proj[1], s_weight_bits[1]);
  auto topoid = w.project(s_weight_proj[0]);
  if (s_nontrivial_zerosector) {
    // Special case: remove offset for preferred masters
    // (see Integral::to_weight() for details).
    --topoid;
  }
  uint32_t sect = s_sector_weight.at(topoid)[wforsect].second;
  uint32_t sectcpy{sect};
  uint32_t lines{0u};
  for (uint32_t dummy = 0u; dummy != s_np; ++dummy) {
    if (sectcpy & 1u) ++lines;
    sectcpy >>= 1;
  }
  uint32_t nums = s_np - lines;
  uint32_t dots;
  uint32_t sps;
  if (s_np == 1) {
    // Special case np=1.
    dots = sect ? dotsp_weight_1 : 0;
    sps = sect ? 0 : dotsp_weight_1;
  }
  else {
    std::tie(dots, sps) = weight_to_dots(dotsp_weight_1, dotsp_weight_2, lines);
  }
  return IntegralProperties{topoid, sect, lines,  nums,
                            dots,   sps,  ppow_w, spow_w};
}

uint32_t Integral::parse_sector(std::string sstr, uint32_t np) {
  // Convert a string to a sector number.
  // Allowed formats:
  // * The sector number as a string,
  // * The sector number in big-endian binary notation, e.g. "b11100" = 7.
  if (!np) {
    throw init_error("parse_sector(): require number of lines.");
  }
  if (np >= 32) {
    throw init_error(
      "parse_sector(): More than 31 propagators are not supported.");
  }
  while (sstr.back() == ' ')
    sstr.pop_back();
  uint32_t sect{0u};
  bool success = true;
  try {
    // Is it a sector number?
    sect = string_to_int(sstr);
  }
  catch (const input_error&) {
    success = false;
  }
  if (success) {
    if (sect >= (1u << np)) {
      throw init_error("parse_sector(): sector number " + sstr +
                       " is too large");
    }
    return sect;
  }
  // Assume that we are dealing with the binary notation.
  success = true;
  bool bits = true;
  uint32_t sector_np{0};
  for (auto chit = sstr.rbegin(); chit != sstr.rend(); ++chit) {
    if (bits) {
      sect <<= 1;
      ++sector_np;
      if (*chit == '1') {
        ++sect;
      }
      else if (*chit == 'b') {
        bits = false;
        sect >>= 1;
        --sector_np;
      }
      else if (*chit != '0') {
        success = false;
        break;
      }
    }
    else {
      if (*chit != ' ') {
        success = false;
        break;
      }
    }
  }
  if (!success) {
    throw init_error("parse_sector(): invalid sector: '" + sstr + "'");
  }
  if (sector_np != np) {
    throw init_error("parse_sector(): '" + sstr +
                     "' has the wrong number of lines.");
  }
  return sect;
}

std::string Integral::bigendian_sector(uint32_t sector, uint32_t np) {
  // Convert a sector number to a string representing the big-endian
  // binary notation, prefixed by 'b', e.g. 7 = "b11100".
  std::string sstr;
  sstr.reserve(np + 1);
  sstr.push_back('b');
  for (uint32_t k = 0; k != np; ++k) {
    if (sector & 1) {
      sstr.push_back('1');
    }
    else {
      sstr.push_back('0');
    }
    sector >>= 1;
  }
  return sstr;
}

std::pair<uint32_t, uint32_t> Integral::weight_to_dots(
    uint32_t dotsp_1, uint32_t dotsp_2, uint32_t /*lines*/) {
  uint32_t dots{0};
  uint32_t sps{0};
  switch (s_dotsp_ordering) {
    case 1:
      dots = dotsp_1;
      sps = dotsp_2;
      break;
    case 2:
      sps = dotsp_1;
      dots = dotsp_2;
      break;
    case 3:
      dots = dotsp_2;
      sps = dotsp_1 - dotsp_2;
      break;
    case 4:
      dots = dotsp_1 - dotsp_2;
      sps = dotsp_2;
      break;
  }
  return {dots, sps};
}

void TopoConfigData::check_node(const YAML::Node& node,
                                const std::string& name) {
  if (!node) throw topoconfig_error(name + " node missing");
}

std::vector<Integral> TopoConfigData::import_integrals(
    const std::string& filename, uint32_t topoid) {
  // Import a list of integrals from a file.
  auto imp = import_basis(filename, topoid);
  if (!imp.second.empty()) {
    throw init_error(std::string("import_integrals(): invalid file format ") +
                     filename);
  }
  return imp.first;
}

std::pair<std::vector<Integral>, std::vector<is_equation>>
TopoConfigData::import_basis(const std::string& filename, uint32_t topoid) {
  // Import a list of integrals and linear combinations of integrals
  // from a file.
  // Linear combinations must be separated by empty (or comment-only) lines.
  // Integrals are assumed to belong to a linear combination, if they are
  // followed by a corfficient.
  // Integrals without a coefficient are assumed to belong to a list of
  // integrals and may not appear in linear combinations.
  // Comments starting with '#', whitespace and empty lines are ignored.
  // Leading '-' is removed (YAML list notation).
  // NOTE: Linear combinations are not sorted by integrals, because
  //       import_basis() is used during setup, i.e. possibly before
  //       integral weights can reliably be calculated.
  std::vector<Integral> basis_igls;
  std::vector<is_equation> basis_lcs;
  std::ifstream instream{filename};
  if (!instream.good()) {
    throw init_error(std::string("Failed reading integral file ") + filename);
  }
  std::string line;
  is_equation tmp_lc;
  std::unordered_set<Integral> igls_in_lc;
  // Elements of a linear combination are only allowed at the beginning of the
  // file, after an empty line, or within an already started linear combination.
  bool allow_lc{true};
  while (pyred::getline_normalized(instream, line)) {
    bool leading_minus{false};
    auto line_cpy{line};
    auto line_comment = split(line, '#', 1, false);
    if (line_comment.size() > 1) {
      line = line_comment.front();
    }
    if (line.find_first_not_of(' ') == std::string::npos) {
      // The line is empty or contains only spaces.
      // If a linear combination was read in, it is complete.
      if (!tmp_lc.empty()) {
        basis_lcs.push_back(tmp_lc);
        tmp_lc.clear();
        igls_in_lc.clear();
      }
      allow_lc = true;
      continue;
    }
    auto pos = line.find_first_not_of(' ');
    if (pos != std::string::npos && line[pos] == '-') {
      Config::log(1) << "WARNING: deprecated integral notation " << line_cpy
                     << " (leading minus) in file \"" << filename << "\"."
                     << std::endl;
      leading_minus = true;
      line[pos] = ' ';
      pos = line.find_first_not_of(' ');
    }
    // Support for the format [T,...]
    if (pos != std::string::npos && line[pos] == '[') {
      line[pos] = ' ';
      auto pos2 = line.find_first_not_of(' ');
      if (line[pos2] == '-' || std::isdigit(line[pos2])) {
        // The format is not [T,...], but [...] without topology.
        line[pos] = '[';
      }
      else {
        pos = line.find(',');
        if (pos == std::string::npos) {
          throw parser_error(std::string("Ill-formed integral string \"") +
                             line_cpy + "\" in file \"" + filename + "\".");
        }
        line[pos] = '[';
      }
    }
    std::string rest;
    std::istringstream ss{line};
    auto igl = Integral(ss, topoid);
    pyred::getline_normalized(ss, rest);
    pos = rest.find_first_not_of(' ');
    if (pos == std::string::npos) {
      // The remainder is empty or only contains whitespace:
      // This is a preferred basis integral.
      // No check for Weight::none() here. The caller should do that.
      if (!tmp_lc.empty()) {
        // We are in the middle of a linear combination!
        throw parser_error(std::string("Integral without coefficient in a "
                                       "linear combination in file \"") +
                           filename + "\".");
      }
      basis_igls.push_back(igl);
      allow_lc = false;
      continue;
    }
    else if (rest[pos] == '*') {
      // The '*' separator between integral and coefficient is optional.
      rest[pos] = ' ';
    }
    if (!allow_lc) {
      throw parser_error(std::string("Invalid starting point for a linear "
                                     "combination in file \"") +
                         filename + "\".");
    }
    // The integral is part of a linear combination
    if (leading_minus) {
      throw parser_error(std::string("Leading minus sign is not allowed in "
                                     "basis linear combinations, file \"") +
                         filename + "\".");
    }
    // Treat duplicate integrals by summing their coefficients.
    if (igls_in_lc.find(igl) == igls_in_lc.end()) {
      igls_in_lc.insert(igl);
      tmp_lc.push_back({igl, rest});
    }
    else {
      for (auto &i_s: tmp_lc) {
        if (i_s.first == igl) {
          i_s.second = i_s.second + "+" + rest;
          break;
        }
      }
    }
  }
  if (!tmp_lc.empty()) {
    // Add the last found linear combination if there is one.
    basis_lcs.push_back(tmp_lc);
  }
  return {std::move(basis_igls), std::move(basis_lcs)};
}

TopoConfigData::TopoConfigData(const std::string& config_dir) {
  // import kinematics.yaml
  std::string kinematics_file = config_dir + "/kinematics.yaml";
  YAML::Node kinematics;
  try {
    kinematics = YAML::LoadFile(kinematics_file);
  }
  catch (const YAML::BadFile& e) {
    throw init_error("Error reading file \"" + kinematics_file + "\"");
  }
  auto kinematics_node = kinematics["kinematics"];
  check_node(kinematics_node, "kinematics");
  auto invariants_node = kinematics_node["kinematic_invariants"];
  check_node(invariants_node, "invariants");
  for (const auto& inv_dim : invariants_node) {
    m_invariants.push_back(inv_dim[0].as<std::string>());
  }
  auto settoone_node = kinematics_node["symbol_to_replace_by_one"];
  if (settoone_node) {
    m_settoone = settoone_node.as<std::string>();
  }
  // import integralfamilies.yaml
  std::string ifam_file = config_dir + "/integralfamilies.yaml";
  YAML::Node integralfamilies;
  try {
    integralfamilies = YAML::LoadFile(ifam_file);
  }
  catch (const YAML::BadFile& e) {
    throw init_error("Error reading file \"" + ifam_file + "\"");
  }
  // import auxiliary.yaml which is automaticaly generated
  std::string auxiliary_file = Integral::ibpdir() + "/auxiliary.yaml";
  YAML::Node auxiliary;
  std::ifstream ifile(auxiliary_file);
  if(static_cast<bool>(ifile)){
    try {
      auxiliary = YAML::LoadFile(auxiliary_file);
    }
    catch (const YAML::BadFile& e) {
      throw init_error("Error reading file \"" + auxiliary_file + "\"");
    }
  }
  auto integralfamilies_node = integralfamilies["integralfamilies"];
  check_node(integralfamilies_node, "integralfamilies");
  for (const auto& topodef : integralfamilies_node) {
    auto name_node = topodef["name"];
    check_node(name_node, "topology name");
    auto topsectors_node = topodef["top_level_sectors"];
    auto propagators_node = topodef["propagators"];
    check_node(propagators_node, "topology propagators");
    uint32_t np = propagators_node.size();
    std::vector<uint32_t> topsectors;
    if (topsectors_node) {
      auto topsectors_str = topsectors_node.as<std::vector<std::string>>();
      for (const auto &sectstr: topsectors_str) {
        topsectors.push_back(Integral::parse_sector(sectstr, np));
      }
    }
    else {
      // If top_level_sectors is not specified, set the top sector to 2^np-1.
      // Note that this sector is usually unphysical.
      topsectors = {(static_cast<uint32_t>(1) << np) - 1};
    }
    auto permutation_node = topodef["permutation"];
    auto permutation = PropagatorPermutation({},np);
    if (permutation_node) {
      permutation = PropagatorPermutation(
        permutation_node.as<std::vector<std::size_t>>(), np);
    }
    if(auxiliary["auxiliary"]){
      auto auxiliary_node = auxiliary["auxiliary"];
      for (const auto& auxdef : auxiliary_node) {
        auto name_tmp_node = auxdef["name"];
        check_node(name_tmp_node, "topology name");
        if(name_tmp_node.as<std::string>() == name_node.as<std::string>()){
          auto permutation_tmp_node = auxdef["permutation"];
          if (permutation_tmp_node) {
            permutation = PropagatorPermutation(
            permutation_tmp_node.as<std::vector<std::size_t>>(), np);
          }
        }
      }
    }
    m_integralfamilies.push_back({name_node.as<std::string>(), np, topsectors,
                                  permutation});
  }
}

std::vector<std::shared_ptr<Topology>>& Integral::setup(
    uint32_t sector_ordering, uint32_t dotsp_ordering,
    const std::string& configdir, const std::string& ibpdir,
    const std::string& basisfile, const std::string& extraRelationsFile, 
    const bool use_li) {
  // Setting configdir will trigger topology import from Kira files.
  // If ibpdir="" (default), ibpdir=configdir.
  // sector_ordering:
  //   1: only by ascending sector number
  //   2: first by number of lines, then by ascending sector number
  //   3 (default): first by number of lines,
  //                then by ascending binary complement of the sector number.
  // dotsp_ordering:
  //   1 (default): dots have higher weight than sps.
  //                I.e. integrals with fewer dots are always simpler
  //                independent of the number of sps.
  //   2: sps have higher weight than dots.
  //   3: first order by dots+sps and if equal, dots have higher weight.
  //   4: first order by dots+sps and if equal, sps have higher weight.
  static std::mutex setup_mtx;
  std::lock_guard<std::mutex> lock(setup_mtx);
  if (Topology::s_next_id) {
    throw init_error("Integral::setup() can only be called"
                     "before creating Topology objects.");
  }
  if (sector_ordering > 3) {
    throw init_error("Integral::setup(): sector ordering must be one of"
                     " {1,2,3} or 0 to keep the default");
  }
  if (dotsp_ordering > 4) {
    throw init_error("Integral::setup(): dot-sp ordering must be one of"
                     " {1,2,3,4} or 0 to keep the default");
  }
  if (sector_ordering) s_sector_ordering = sector_ordering;
  if (dotsp_ordering) s_dotsp_ordering = dotsp_ordering;
  s_ibpdir = ibpdir;
  if (ibpdir.empty()) s_ibpdir = configdir + "/../sectormappings";
  s_use_li = use_li;
  if (!configdir.empty()) {
    // Import integral families.
    auto topoconf = TopoConfigData(configdir);
    CoeffHelper::settoone(topoconf.m_settoone);
    CoeffHelper::add_invariant("d");
    for (const auto& inv : topoconf.m_invariants) {
      CoeffHelper::add_invariant(inv);
    }
    for (const auto& ifam : topoconf.m_integralfamilies) {
      new_topology(ifam.m_name, ifam.m_np, ifam.m_topsectors,
                   ifam.m_permutation);
    }
  }
  Config::finish();
  if (!extraRelationsFile.empty()) {
    auto extra_relations = TopoConfigData::import_basis(extraRelationsFile);
    // Discard integrals in trivial sectors.
    std::vector<pow_type> lc_pseudoigl(np(), 0);
    for (const auto& lc : extra_relations.second) {
      auto maxweight = Weight(0u);
      is_equation lc_cp;
      for (const auto& i_s : lc) {
        auto w = i_s.first.to_weight();
        if (w != Weight::none()) {
          lc_cp.push_back(i_s);
          maxweight = std::max(maxweight, w);
        }
      }
      if(maxweight==Weight(0u))
        maxweight = Weight::none();
      if(!is_zero(maxweight)){
        auto props = properties(maxweight);
        auto ins = s_external_equations_map.lcs.insert({{props.topology, props.sector}, {}});
        ins.first->second.push_back(std::move(lc_cp));
        //++s_external_equations.nimported_ffs;
        //s_external_equations.eqs.push_back(std::move(lc_cp));
      }
    }
  }
  if (!basisfile.empty()) {
    auto pref_basis = TopoConfigData::import_basis(basisfile);
    s_preferred_masters.igls = std::move(pref_basis.first);
    // Group the basis linear combinations by sector in an unordered_map.
    // A linear combination belongs to the topology/sector of its highest
    // integral. Discard integrals in trivial sectors.
    int basis_lc{0};
    std::vector<pow_type> lc_pseudoigl(np(), 0);
    for (const auto& lc : pref_basis.second) {
      lc_pseudoigl[0] = ++basis_lc;
      auto maxweight = Weight(0u);
      is_equation lc_cp;
      for (const auto& i_s : lc) {
        auto w = i_s.first.to_weight();
        if (w != Weight::none()) {
          lc_cp.push_back(i_s);
          maxweight = std::max(maxweight, w);
        }
      }
      // Basis linear combinations are represented by BASISLC[n,0,0,...]
      // with n=1,2,...
      lc_cp.push_back({Integral(Topology::s_next_id, lc_pseudoigl), "-1"});
      auto props = properties(maxweight);
      auto ins =
          s_preferred_masters.lcs.insert({{props.topology, props.sector}, {}});
      ins.first->second.push_back(std::move(lc_cp));
    }
    s_preferred_masters.n_basis_lcs = basis_lc;
    // Ignore preferred masters which cause a weight overflow.
    // If they lie in the seed, weight bits must be set anyway and the
    // preferred masters will be set again.
    // If they lie outside the seed, they can be ignored.
    set_preferred_masters(false);
  }
  return Topology::s_id_to_topo;
}

void Integral::set_preferred_masters(bool fatal) {
  uint64_t mapped_weight{0u};
  std::unordered_map<Integral, Weight> preferred_masters_i2w{};
  std::vector<Weight> weight_map{};
  // Clear caches unconditionally. This is required if the weight bits changed.
  s_preferred_masters_i2w.clear();
  clear_cache(3);
  for (const auto& igl : s_preferred_masters.igls) {
    Weight weight;
    try {
      weight = igl.to_weight();
    }
    catch (const weight_error& e) {
      if (fatal) {
        throw weight_error(e.what());
      }
      // Otherwise ignore and proceed with the next integral.
      continue;
    }
    if (weight != Weight::none()) {
      // Do not populate the actual maps yet, because then the caches
      // would initialise themselves with incomplete maps
      // as soon as Integral::to_weight() is called.
      preferred_masters_i2w.insert({igl, Weight(++mapped_weight)});
      weight_map.push_back(weight);
    }
  }
  // Clear s_i2w_cache and s_w2w_cache which may have been populated
  // by Integral::to_weight(), so that s_i2w_cache will be initialised
  // with the complete s_preferred_masters_i2w.
  clear_cache(s_cache_level);
  // Assign the fully populated maps.
  s_preferred_masters_i2w.get_unsafe() = preferred_masters_i2w;
  s_weight_map = weight_map;
}

const Integral::preferred_basis_t& Integral::get_preferred_masters() {
  return s_preferred_masters;
}

const std::vector<is_equation>& Integral::get_preferred_masters(
    uint32_t topoid, uint32_t sector) {
  static std::vector<is_equation> empty{};
  auto it = s_preferred_masters.lcs.find({topoid, sector});
  if (it == s_preferred_masters.lcs.end()) {
    return empty;
  }
  else {
    return it->second;
  }
}

const Integral::external_eqs_t& Integral::get_external_equations() {
  return s_external_equations;
}

const Integral::external_eqs_map_t& Integral::get_external_equations_map() {
  return s_external_equations_map;
}

std::vector<uint32_t> Integral::assign_weight_bits(
    const std::vector<uint32_t>& weight_bits) {
  if (s_np == 1u) {
    s_weight_bits[2] = 32;
    s_weight_bits[3] = 0;
    s_weight_bits[4] = 0;
    s_weight_bits[5] = 0;
  }
  else if (!weight_bits.empty()) {
    if (weight_bits.size() != s_weight_bits.size() - 2) {
      std::ostringstream ss;
      ss << "Integral::assign_weight_bits(): weight_bits must have "
         << s_weight_bits.size() - 2 << " elements.";
      throw init_error(ss.str());
    }
    for (std::size_t n = 2; n != 6; ++n) {
      s_weight_bits[n] = weight_bits[n - 2];
    }
  }
  uint32_t nontopo_weight = s_np + s_weight_bits[2] + s_weight_bits[3] +
                            s_weight_bits[4] + s_weight_bits[5];
  if (nontopo_weight > Weight::width) {
    std::ostringstream ss;
    ss << "Integral weight representation exceeds " << Weight::width << " bits.";
    if (Weight::width == 64u) {
      ss << "\nIf you really -- and I mean *REALLY* -- need more, build Kira "
            "with the option weight_width=128.";
    }
    throw init_error(ss.str());
  }
  if (s_weight_bits[2] > 8*sizeof(s_weight_proj[2])) {
    throw input_error(std::string("dotsp1 weight exceeds projector width ") +
      "(requires " + std::to_string(s_weight_bits[2]) + " bits ).");
  }
  if (s_weight_bits[3] > 8*sizeof(s_weight_proj[3])) {
    throw input_error(std::string("dotsp2 weight exceeds projector width ") +
      "(requires " + std::to_string(s_weight_bits[3]) + " bits ).");
  }
  if (s_weight_bits[4] > 8*sizeof(s_weight_proj[4])) {
    throw input_error(std::string("ppows weight exceeds projector width ") +
      "(requires " + std::to_string(s_weight_bits[4]) + " bits ).");
  }
  if (s_weight_bits[5] > 8*sizeof(s_weight_proj[5])) {
    throw input_error(std::string("spows weight exceeds projector width ") +
      "(requires " + std::to_string(s_weight_bits[5]) + " bits ).");
  }
  // topology weight: the rest, but at most as much as fits into the projector.
  s_weight_bits[0] = Weight::width - nontopo_weight;
  if (s_weight_bits[0] > 8*sizeof(s_weight_proj[0])) {
    s_weight_bits[0] = 8*sizeof(s_weight_proj[0]);
  }
  s_weight_bits[1] = s_np; // sector weight
  for (std::size_t n = 0; n != 6; ++n) {
    if (s_weight_bits[n] == 8*sizeof(s_weight_proj[n])) {
      // Needs special treatment because 1u<<32 is undefined behaviour.
      s_weight_proj[n] = std::numeric_limits<uint32_t>::max();
    }
    else {
      s_weight_proj[n] = (static_cast<uint32_t>(1) << s_weight_bits[n]) - 1;
    }
  }
  // Set the preferred masters with the updated weight bits.
  set_preferred_masters(true);
  return {s_weight_bits[2], s_weight_bits[3], s_weight_bits[4],
          s_weight_bits[5]};
}

std::vector<uint32_t> Integral::assign_weight_bits(
    uint32_t maxdots, uint32_t maxsps,
    uint32_t maxdotcombs, uint32_t maxspcombs) {
  uint32_t dotsp_1{0};
  uint32_t dotsp_2{0};
  switch (s_dotsp_ordering) {
    case 1:
      dotsp_1 = maxdots + 1;
      dotsp_2 = maxsps + 1;
      break;
    case 2:
      dotsp_1 = maxsps + 1;
      dotsp_2 = maxdots + 1;
      break;
    case 3:
      dotsp_1 = maxdots + maxsps + 2;
      dotsp_2 = maxdots + 1;
      break;
    case 4:
      dotsp_1 = maxdots + maxsps + 2;
      dotsp_2 = maxsps + 1;
      break;
  }
  s_weight_bits[2] = bitwidth(dotsp_1);
  s_weight_bits[3] = bitwidth(dotsp_2);
  s_weight_bits[4] = bitwidth(maxdotcombs);
  s_weight_bits[5] = bitwidth(maxspcombs);
  if (Topology::s_next_id) {
    auto topo_weight_bits = bitwidth(Topology::s_next_id - 1);
    if (s_weight_bits[2] <= s_default_weight_bits[0] &&
        s_weight_bits[3] <= s_default_weight_bits[1] &&
        s_weight_bits[4] <= s_default_weight_bits[2] &&
        s_weight_bits[5] <= s_default_weight_bits[3] &&
        (topo_weight_bits + s_np + s_default_weight_bits[0] +
         s_default_weight_bits[1] + s_default_weight_bits[2] +
         s_default_weight_bits[3] <= Weight::width)) {
      // If the default weight bit representation fits (including bits for
      // the sector and the topology), use the default weight bits.
      s_weight_bits[2] = s_default_weight_bits[0];
      s_weight_bits[3] = s_default_weight_bits[1];
      s_weight_bits[4] = s_default_weight_bits[2];
      s_weight_bits[5] = s_default_weight_bits[3];
    }
    else if (topo_weight_bits + s_np + s_weight_bits[2] + s_weight_bits[3] +
              s_weight_bits[4] + s_weight_bits[5] > Weight::width) {
      std::ostringstream ss;
      ss << "Integral weight representation exceeds " << Weight::width << " bits.";
      if (Weight::width == 64u) {
        ss << "\nIf you really -- and I mean *REALLY* -- need more, build Kira "
              "with the option weight_width=128.";
      }
      throw init_error(ss.str());
    }
  }
  return assign_weight_bits();
}

std::vector<uint32_t> Integral::assign_weight_bits(
    const std::vector<SeedSpec>& seedspec) {
  if (!seedspec.empty()) {
    // Calculate the needed weight bits from the integral ordering
    // and the seed specification.
    auto maxdots4lines = std::vector<uint32_t>(s_np, 0u);
    uint32_t maxdots{0};
    uint32_t maxsps{0};
    for (const auto& spec : seedspec) {
      uint32_t dots{spec.m_mindots};
      for (uint32_t lines = spec.m_lines; lines >= s_minlines; --lines) {
        if (!lines) break;
        if (lines < spec.m_n1 && dots < spec.m_maxdots) {
          ++dots;
        }
        maxdots4lines[lines - 1] = std::max(maxdots4lines[lines - 1], dots);
      }
      maxsps = std::max(maxsps, spec.m_sps);
    }
    uint32_t maxdotcombs{0};
    // TODO: Take increments from imported IBP relations.
    //       Need to make sure that IBPs are already imported.
    uint32_t ibp_dot_increment = 1u;
    uint32_t ibp_sp_increment = 1u;
    for (uint32_t lines = 1; lines <= s_np; ++lines) {
      auto dotsonlines = maxdots4lines[lines - 1];
      maxdotcombs = std::max(
        static_cast<std::size_t>(maxdotcombs),
        binomial(lines + dotsonlines - 1 + ibp_dot_increment,
                 dotsonlines + ibp_dot_increment));
      maxdots = std::max(maxdots, dotsonlines);
    }
    auto maxspcombs = binomial(
      s_np - s_minlines + maxsps - 1 + ibp_sp_increment,
      maxsps + ibp_sp_increment);
    assign_weight_bits(maxdots, maxsps, maxdotcombs, maxspcombs);
  }
  else {
    assign_weight_bits();
  }
  return {s_weight_bits[2], s_weight_bits[3], s_weight_bits[4],
          s_weight_bits[5]};
}

void Integral::np(uint32_t newnp) {
  if (s_np != newnp) {
    if (s_np) {
      throw init_error(
          "All topologies must have the same number of propagators.");
    }
    else {
      if (newnp >= 32u) {
        throw init_error(
          "np(): More than 31 propagators are not supported.");
      }
      s_np = newnp;
      assign_weight_bits();
    }
  }
}

Seeds Integral::seeds(uint32_t sect, uint32_t dots, uint32_t sps) {
  return Seeds(s_np, sect, dots, sps);
}

Weight Integral::to_weight() const {
  // Convert an integral to an integer weight.
  if (!s_np) {
    // Special case: user-defined weights.
    return *reinterpret_cast<const Weight*>(m_powers.data());
  }
  s_i2w_cache.set_if_empty(s_preferred_masters_i2w);
  if (s_cache_level & 2 || !s_preferred_masters_i2w.empty()) {
    auto found = s_i2w_cache.lookup(*this);
    if (found.second) {
      return found.first;
    }
  }
  if (m_topoid == Topology::s_next_id) {
    for (std::size_t k = 1; k < s_np; ++k) {
      if (m_powers[k] != 0) {
        throw weight_error(
            "Invalid basis LC or form factor Integral instance.");
      }
    }
    if (m_powers[0] >= 1 &&
        m_powers[0] <= static_cast<pow_type>(s_preferred_masters.n_basis_lcs)) {
      // The integral represents a basis linear combination.
      return Weight(s_preferred_masters.igls.size() + m_powers[0]);
    }
    else if (m_powers[0] <= -1 &&
             static_cast<uint32_t>(-m_powers[0]) <= s_weight_proj[4]) {
      // s_weight_proj[4] is used as the maximal possible form factor number.
      // The integral represents a form factor.
      return Weight::amplitude(-m_powers[0]);
    }
    else {
      throw weight_error(
          "Invalid basis LC or form factor Integral instance (overflow).");
    }
  }
  uint32_t lines{0};
  uint32_t dots{0};
  uint32_t nums{0};
  uint32_t sps{0};
  uint32_t dotsp_weight_1{0};
  uint32_t dotsp_weight_2{0};
  std::vector<uint32_t> ppows;
  std::vector<uint32_t> spows;
  uint32_t sector{0};
  uint32_t linepos{0};
  const auto &perm = Topology::id_to_topo(m_topoid)->m_permutation;
  if (perm.is_identity()) {
    for (const auto p : m_powers) {
      if (p > 0) {
        sector += (1 << linepos);
        ++lines;
        dots += (p - 1);
        ppows.push_back(p - 1);
      }
      else {
        sps -= p;
        spows.push_back(-p);
      }
      ++linepos;
    }
  }
  else {
    for (const auto pos : perm) {
      auto p = m_powers[pos];
      if (p > 0) {
        sector += (1 << pos);
        ++lines;
        dots += (p - 1);
        ppows.push_back(p - 1);
      }
      else {
        sps -= p;
        spows.push_back(-p);
      }
    }
  }
  if (s_trivialsectors.at(m_topoid).find(sector) !=
      s_trivialsectors.at(m_topoid).end()) {
    return Weight::none();
  }
  auto topo_w = m_topoid;
  if (s_nontrivial_zerosector) {
    // Effectively increase topoid by 1 as an offset, so that all lower weights
    // can be used for preferred masters. This is only needed when sector 0 is
    // regarded as non-trivial (otherwise the sector weight serves as offset).
    ++topo_w;
  }
  auto weight = Weight(topo_w);
  weight.shift_set(s_weight_bits[1], s_sector_weight[m_topoid][sector].first);
  uint32_t ppows_weight;
  uint32_t spows_weight;
  if (s_np == 1) {
    // Special case np=1
    dotsp_weight_1 = sector ? dots : sps;
    dotsp_weight_2 = 0u;
    ppows_weight = 0u;
    spows_weight = 0u;
  }
  else {
    nums = s_np - lines;
    // dot and sp composition weights
    // * dots to the left have lower weight
    // const auto& cid = get_composition_id(dots,lines);
    // ppows_weight = cid.second.size() - cid.first.at(ppows) - 1;
    // * dots to the right have lower weight
    ppows_weight = get_composition_id(dots, lines).first.at(ppows);
    // * sps to the right have lower weight
    // spows_weight = get_composition_id(sps,nums).first.at(spows);
    // * sps to the left have lower weight
    const auto& cid = get_composition_id(sps, nums);
    spows_weight = cid.second.size() - cid.first.at(spows) - 1;
    // (dots, sps) ordering
    switch (s_dotsp_ordering) {
      case 1:
        dotsp_weight_1 = dots;
        dotsp_weight_2 = sps;
        break;
      case 2:
        dotsp_weight_1 = sps;
        dotsp_weight_2 = dots;
        break;
      case 3:
        dotsp_weight_1 = dots + sps;
        dotsp_weight_2 = dots;
        break;
      case 4:
        dotsp_weight_1 = dots + sps;
        dotsp_weight_2 = sps;
        break;
    }
  }
  weight.shift_set(s_weight_bits[2], dotsp_weight_1);
  weight.shift_set(s_weight_bits[3], dotsp_weight_2);
  // Changes necessary to give dot composition has higher weight
  // than sp composition:
  //   weight.shift_set(s_weight_bits[4], ppows_weight);
  //   weight.shift_set(s_weight_bits[5], spows_weight);
  // Here: sp composition has higher weight than dot composition.
  weight.shift_set(s_weight_bits[5], spows_weight);
  weight.shift_set(s_weight_bits[4], ppows_weight);
  if (s_weight_proj[0] < topo_w ||
      s_weight_proj[2] < dotsp_weight_1 ||
      s_weight_proj[3] < dotsp_weight_2 ||
      s_weight_proj[4] < ppows_weight ||
      s_weight_proj[5] < spows_weight) {
    weight_overflow(topo_w, dotsp_weight_1, dotsp_weight_2,
                    ppows_weight, spows_weight, lines);
  }
  // Insert the element into the cache(s)
  // if it hasn't been cached in the meantime.
  if (s_cache_level & 2) {
    s_i2w_cache.insert(*this, weight);
  }
  if (s_cache_level & 1 && weight != Weight::none()) {
    // Do not cache Weight::none();
    // constructing Integral(Weight::none()) is not allowed.
    s_w2i_cache.insert(weight, *this);
  }
  return weight;
}

void Integral::weight_overflow(uint32_t topo, uint32_t dotsp_1,
                               uint32_t dotsp_2, uint32_t ppow_w,
                               uint32_t spow_w, uint32_t lines) {
  // Error handling in case a sub-weight does not fit
  // in the Integral weight representation, because too few bits
  // are reserved for the particular sub-weight.
  std::ostringstream ss;
  ss << "Integral weight overflow in";
  uint32_t dots;
  uint32_t sps;
  std::tie(dots, sps) = weight_to_dots(dotsp_1, dotsp_2, lines);
  if (topo > s_weight_proj[0]) {
    ss << " topology weight (" << topo << " does not fit in "
       << s_weight_bits[0] << " bits);";
  }
  if (dotsp_1 > s_weight_proj[2]) {
    switch (s_dotsp_ordering) {
      case 1:
        ss << " dots (" << dots << " do not fit in " << s_weight_bits[2]
           << " bits);";
        break;
      case 2:
        ss << " sps (" << sps << " do not fit in " << s_weight_bits[2]
           << " bits);";
        break;
      case 3:
      case 4:
        ss << " dots+sps (" << dots << "+" << sps << " does not fit in "
           << s_weight_bits[2] << " bits);";
        break;
    }
  }
  if (dotsp_2 > s_weight_proj[3]) {
    switch (s_dotsp_ordering) {
      case 1:
      case 4:
        ss << " sps (" << sps << " do not fit in " << s_weight_bits[3]
           << " bits);";
        break;
      case 2:
      case 3:
        ss << " dots (" << dots << " do not fit in " << s_weight_bits[3]
           << " bits);";
        break;
    }
  }
  if (ppow_w > s_weight_proj[4]) {
    ss << " dot distribution (" << ppow_w << " does not fit in "
       << s_weight_bits[4] << " bits);";
  }
  if (spow_w > s_weight_proj[5]) {
    ss << " sp distribution (" << spow_w << " does not fit in "
       << s_weight_bits[5] << " bits);";
  }
  ss << "\nuse Integral::setup() to reserve more space where needed.";
  throw weight_error(ss.str());
}

Integral Integral::zip(const std::vector<pow_type>& merge) const {
  // Return a new Integral with 'm_powers' and 'merge' added componentwise.
  std::vector<pow_type> pows;
  pows.reserve(s_np);
  auto p1 = m_powers.cbegin();
  for (auto p2 : merge) {
    pows.push_back(*p1 + p2);
    ++p1;
  }
  return Integral(m_topoid, std::move(pows));
}

void Integral::use_cache(int level) { s_cache_level = level; }

void Integral::clear_cache(int level) {
  if (!level) {
    level = s_auto_clear_cache_level;
  }
  if (level & 1) {
    s_w2i_cache.clear();
  }
  if (level & 2) {
    s_i2w_cache.clear();
    s_i2w_cache.set_if_empty(s_preferred_masters_i2w);
  }
}

void Integral::auto_clear_cache(int level) { s_auto_clear_cache_level = level; }

bool Integral::operator<(const Integral& rhs) const {
  return (to_weight() < rhs.to_weight());
}

std::string Integral::to_string(const StringFormat& fmt) const {
  std::ostringstream ss;
  if (fmt == StringFormat::pyred && !Integral::np()) {
    // Special case: user-defined weight: print the weight.
    ss << this->to_weight();
    return ss.str();
  }
  if (fmt != StringFormat::indices) {
    if (m_topoid == Topology::s_next_id) {
      if (m_powers[0] > 0) {
        // The integral represents a basis linear combination.
        ss << Topology::s_basislc_name;
      }
      else { // actually m_powers[0] < 0
        // The integral represents a form factor.
        ss << Topology::s_formfactor_name;
      }
    }
    else {
      ss << Topology::id_to_topo(m_topoid)->m_name;
    }
    if (fmt == StringFormat::form) {
      ss << "(";
    }
    else {
      ss << "[";
    }
  }
  if (!Integral::np()) {
    ss << this->to_weight();
  }
  else if (m_topoid == Topology::s_next_id) {
    if (m_powers[0] > 0) {
      ss << m_powers[0];
    }
    else {
      ss << -m_powers[0];
    }
  }
  else {
    for (auto it = m_powers.cbegin(); it != m_powers.cend(); ++it) {
      if (it != m_powers.cbegin()) ss << ",";
      ss << *it;
    }
  }
  if (fmt != StringFormat::indices) {
    if (fmt == StringFormat::form) {
      ss << ")";
    }
    else {
      ss << "]";
    }
  }
  return ss.str();
}

bool operator==(const Integral& lhs, const Integral& rhs) {
  return (lhs.m_topoid == rhs.m_topoid && lhs.m_powers == rhs.m_powers);
}

bool operator!=(const Integral& lhs, const Integral& rhs) {
  return (lhs.m_topoid != rhs.m_topoid || lhs.m_powers != rhs.m_powers);
}

std::ostream& operator<<(std::ostream& out, const Integral& igl) {
  out << igl.to_string();
  return out;
}

std::vector<Integral> Integral::add_formfactors(const std::string& ff_file) {
  auto imp = TopoConfigData::import_basis(ff_file);
  if (!imp.first.empty()) {
    throw input_error(std::string("Form factor file \"") + ff_file +
                      "\": no individual integrals without factor allowed.");
  }
  return add_formfactors(imp.second);
}

std::vector<Integral> Integral::add_formfactors(
    const std::vector<is_equation>& ffs) {
  std::vector<Integral> ff_placeholders;
  std::vector<pow_type> pseudoigl(np(), 0);
  for (const auto& ff : ffs) {
    is_equation ff_eq;
    ff_eq.reserve(ff.size() + 1);
    pseudoigl[0] = -static_cast<pow_type>(++s_external_equations.nimported_ffs);
    ff_placeholders.push_back(Integral(Topology::s_next_id, pseudoigl));
    ff_eq.push_back({ff_placeholders.back(), "-1"});
    for (const auto& i_s : ff) {
      ff_eq.push_back(i_s);
    }
    s_external_equations.eqs.push_back(std::move(ff_eq));
  }
  return ff_placeholders;
}

void Integral::add_equations(const std::string& eqs_file) {
  auto imp = TopoConfigData::import_basis(eqs_file);
  if (!imp.first.empty()) {
    throw input_error(std::string("Equations file \"") + eqs_file +
                      "\": no individual integrals without factor allowed.");
  }
  add_equations(std::move(imp.second));
}

void Integral::add_equations(std::vector<is_equation>&& eqs) {
  for (auto& eq : eqs) {
    s_external_equations.eqs.push_back(std::move(eq));
  }
}

/******************
 * Topology class *
 ******************/

// static members
uint32_t Topology::s_next_id{0};
std::vector<std::shared_ptr<Topology>> Topology::s_id_to_topo{};
std::unordered_map<std::string, uint32_t> Topology::s_toponame_to_id{};
std::string Topology::s_basislc_name{"BASISLC"};
std::string Topology::s_formfactor_name{"FORMFACTOR"};
std::string Topology::s_userweights_name{"T"};

Topology::Topology(const std::string& name, const uint32_t np,
                   const std::vector<uint32_t>& topsects,
                   const std::vector<uint32_t>& trivialsectors,
                   relations_ptr&& relations,
                   const PropagatorPermutation &permutation)
    : m_name{name},
      m_np{np},
      m_id{s_next_id},
      m_topsectors{topsects},
      m_trivialsectors{trivialsectors.cbegin(), trivialsectors.cend()},
      m_permutation{permutation},
      m_relations{std::move(relations)}
{
  Integral::s_auto_topologies = false;
  // Set np if not yet set yet, check if it is the same if already set.
  // Assign weight bits.
  Integral::np(m_np);
  if (m_permutation.size() != m_np) {
    std::ostringstream ss;
    ss << "Topology constructor: np=" << m_np
       << " is inconsistent with the size " << m_permutation.size()
       << " of the permutation";
    throw init_error(ss.str());
  }
  if (m_relations && get_relations_topoid(m_relations) != m_id) {
    std::ostringstream ss;
    ss << "Topology constructor: cannot assign IntegralRelations "
       << "with topology ID " << get_relations_topoid(m_relations)
       << " to the Topology with ID " << m_id << ".";
    throw init_error(ss.str());
  }
  assign_sector_weight();
  Integral::s_trivialsectors.push_back(m_trivialsectors);
  if (m_trivialsectors.find(0u) == m_trivialsectors.end()) {
    Integral::s_nontrivial_zerosector = true;
  }
  Integral::s_sector_weight.push_back(m_sector_weight);
  s_toponame_to_id.insert({m_name, m_id});
  ++s_next_id;
}

std::vector<uint32_t> Topology::import_trivialsectors(const std::string &tsfile)
{
  // Trivial sector format: "sn1,sn2,..."
  std::vector<uint32_t> trivialsectors;
  std::ifstream instream{tsfile};
  if (!instream.good()) {
    std::ostringstream ss;
    ss << "Failed reading trivial sector file " << tsfile;
    throw init_error(ss.str());
  }
  std::string line;
  pyred::getline_normalized(instream, line);
  std::istringstream ss{line};
  uint32_t sn;
  char sep;
  do {
    ss >> sn;
    trivialsectors.push_back(sn);
  } while (static_cast<bool>(ss >> sep));
  return trivialsectors;
}

std::shared_ptr<Topology> new_topology(
    const std::string& name,
    const uint32_t np,
    const std::vector<uint32_t>& topsects,
    const std::vector<uint32_t>& trivialsectors,
    relations_ptr&& relations,
    const PropagatorPermutation &permutation) {
  // make_shared needs a public constructor.
  auto perm = permutation;
  if (perm.size() == std::numeric_limits<std::size_t>::max()) {
    // Case of a default-constructed permutation.
    perm = PropagatorPermutation(np);
  }
  Topology::s_id_to_topo.push_back(std::shared_ptr<Topology>{
      new Topology{name, np, topsects, trivialsectors, std::move(relations),
      perm}});
  return Topology::s_id_to_topo.back();
}

std::shared_ptr<Topology> new_topology(
    const std::string& name,
    const uint32_t np,
    const std::vector<uint32_t>& topsects,
    const PropagatorPermutation &permutation) {
  // Take trivial sectors and integral relations from Kira files.
  std::string tsfile = Integral::ibpdir() + "/" + name + "/trivialsector";
  auto trivialsectors = Topology::import_trivialsectors(tsfile);
  // Set np if not yet set, check if it is the same if already set.
  // Assign weight bits.
  Integral::np(np);
  relations_ptr relations = import_relations(name, np, Topology::s_next_id);
  return new_topology(name, np, topsects, trivialsectors, std::move(relations),
                      permutation);
}


std::vector<std::pair<uint32_t,uint32_t>> Topology::sector_weight_table(
    uint32_t sector_ordering,
    const std::vector<uint32_t> &topsectors,
    const PropagatorPermutation &permutation) {
  // sector weight table
  //   sector_weight[sector].first = weight
  //   sector_weight[weight].second = sector
  auto np = permutation.size();
  if (np == std::numeric_limits<std::size_t>::max()) {
    throw init_error("Topology::sector_weight_table(): "
      "default constructed propagator permutation argument is not allowed.");
  }
  std::vector<std::array<uint32_t, 3>> sect_w;
  for (uint32_t s = 0; s != (1u << np); ++s) {
    // Note that s is the sector refering to the permuted propagators,
    // if a permutation is set for the topology.
    uint32_t scpy{s};
    uint32_t lines{0};
    uint32_t lines_w{0};
    for (uint32_t l = 0; l != np; ++l) {
      // lines_w for sector_ordering==3
      if (scpy & 1) {
        ++lines;
      }
      else {
        ++lines_w;
      }
      scpy >>= 1;
      lines_w <<= 1;
    }
    // If a sector is not a subsector of top_sectors add 1<<np to its sector
    // weight so that it's higher than those of all subsectors of top_sectors.
    bool outside_topsects{true};
    for (const auto ts : topsectors) {
      const auto permuted_ts = permutation.permute_sector(ts);
      if (!(s & (~permuted_ts))) outside_topsects = false;
    }
    if (outside_topsects) {
      lines_w += 1u << np;
    }
    if (sector_ordering == 2) {
      // lines weight = sector number, when sectors are ordered
      // first by number of lines, then by lines weight.
      lines_w = s + (outside_topsects ? 1u << np : 0u);
    }
    if (sector_ordering == 1) {
      // sectors ordered by sector number only
      sect_w.push_back({{s, 0, s}});
    }
    else {
      // sectors first ordered by number of lines
      sect_w.push_back({{s, lines, lines_w}});
    }
  }
  // Sort by (number of lines, weight for fixed lines) --> global weight.
  std::sort(sect_w.begin(), sect_w.end(),
            [](const std::array<uint32_t, 3>& slw1,
               const std::array<uint32_t, 3>& slw2) {
              return slw1[1] == slw2[1] ? slw1[2] < slw2[2] : slw1[1] < slw2[1];
            });
  // sect_w is now sorted by weight.
  std::vector<std::pair<uint32_t,uint32_t>> sector_weight;
  sector_weight.resize(sect_w.size());
  for (uint32_t sw = 0; sw != sect_w.size(); ++sw) {
    auto sector = permutation.permute_sector(sect_w[sw][0], /*inverse*/true);
    sector_weight[sw].second = sector;
    sector_weight[sector].first = sw;
  }
//   for (uint32_t sw = 0; sw != sect_w.size(); ++sw) {
//     std::cout << "sector of weight " << sw << ": "
//               << sector_weight[sw].second << " "
//               << Integral::bigendian_sector(sector_weight[sw].second, np)
//               << " | weight of sector " << sw << " "
//               << Integral::bigendian_sector(sw, np) << ": "
//               << sector_weight[sw].first << std::endl;
//   }
  return sector_weight;
}


void Topology::assign_sector_weight() {
  // assign m_sector_weight and m_minlines
  m_sector_weight = sector_weight_table(
    Integral::s_sector_ordering, m_topsectors, m_permutation);
  m_minlines = m_np;
  for (uint32_t s = 0; s != (1u << m_np); ++s) {
    // m_minlines: minimal number of lines in any non-trivial sector.
    auto lines = count_set_bits(s);
    if (lines < m_minlines &&
        m_trivialsectors.find(s) == m_trivialsectors.end()) {
      m_minlines = lines;
    }
  }
  Integral::s_minlines = std::min(Integral::s_minlines, m_minlines);
}


SeedSpec Topology::seed_spec(uint32_t sector, uint32_t rmax, uint32_t sps,
                             int maxdots, int n1,
                             SeedSpec::Recursive rec) const {
  return SeedSpec(*this, sector, rmax, sps, maxdots, n1, rec);
}


std::vector<std::pair<std::string,uint32_t>> Topology::parse_topology_ordering(
    const std::string &filename) {
  std::ifstream instream{filename};
  if (!instream.good()) {
    throw input_error(
      std::string("Failed reading topology_ordering file ") + filename);
  }
  std::string line;
  std::vector<std::pair<std::string,uint32_t>> vec_topo_np;
  while (pyred::getline_normalized(instream, line)) {
    auto line_cpy{line};
    auto line_comment = split(line, '#', 1, false);
    if (line_comment.size() > 1) {
      line = line_comment.front();
    }
    auto split_line = split(line, ' '); // effectively trims all items
    if (split_line.size() != 2u) {
      throw input_error(
        std::string("Invalid topology_ordering file, line ") + line_cpy);
    }
    vec_topo_np.push_back({split_line[0], parse_int<uint32_t>(split_line[1])});
  }
  return vec_topo_np;
}

} // namespace pyred
