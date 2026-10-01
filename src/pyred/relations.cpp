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

#include <algorithm>
#include <cassert>
#include <fstream>
#include <limits>
#include <regex>
#include <unordered_set>

#include "pyred/coeff_int.h" // need std::hash<Coeff>
#include "pyred/coeff_vec.h" // for all Coeff classes
#include "pyred/parser.h"
#include "pyred/ppmacros.h"
#include "relations.h"

namespace pyred {

namespace {

std::string coeff_patt(const std::string& cf) {
  std::string out;
  out.reserve(cf.size());
  auto modify_token = [&out](const std::string& tk) {
    std::istringstream ss{tk};
    char a;
    ss >> a;
    if (a == 'a') {
      int i;
      auto found_int = bool(ss >> i);
      std::string rest;
      pyred::getline_normalized(ss, rest);
      if (found_int && rest.empty()) {
        // allow whitespace:
        // rest.find_first_not_of(' ') != std::string::npos
        out += "$";
        out += static_cast<char>(i + 48);
        return;
      }
    }
    out += tk;
  };
  std::regex re{"\\ba(\\d+)\\b"};
  std::for_each(std::sregex_token_iterator(cf.cbegin(), cf.cend(), re, {-1, 0}),
                std::sregex_token_iterator(), modify_token);
  return out;
}

std::pair<Integral, std::string> import_eqcf(const std::string& s,
                                             uint32_t topoid) {
  std::istringstream ss{s};
  std::string coeff;
  auto igl = Integral(ss, topoid, true); // force topoid to ignore "+ INT"
  ss >> std::ws;
  pyred::getline_normalized(ss, coeff);
  if (!coeff.empty() && coeff[0] == '*') {
    coeff.erase(coeff.begin());
  }
  return {std::move(igl), coeff_patt(coeff)};
}

std::string sow_coeff(const std::string& cf_tmpl,
                      const std::vector<std::string>& seed) {
  std::string coeff;
  coeff.reserve(cf_tmpl.size());
  for (auto ch = cf_tmpl.cbegin(); ch != cf_tmpl.cend(); ++ch) {
    if (*ch == '$') {
      coeff += seed[*(++ch) - 48];
    }
    else {
      coeff.push_back(*ch);
    }
  }
  return coeff;
}

} // namespace

int GeneratorHelper::countJO = 0;
int GeneratorHelper::nofbunches = 0;

std::size_t neqs_from_seed(
    const std::pair<SeedSpecSector, SeedSpecSector>& ibp_sym_spec4sect) {
  const auto& ibp_spec4sect = ibp_sym_spec4sect.first;
  const auto& sym_spec4sect = ibp_sym_spec4sect.second;
  std::size_t neqs{0};
  // Number of IBP equations.
  for (const auto& dotsp : ibp_spec4sect.expand()) {
    neqs += ibp_spec4sect.m_topo->relations().n_ibps() *
            dotsp_combs(Integral::np(), lines_from_sn(ibp_spec4sect.m_sector),
                        dotsp.first, dotsp.second);
  }
  // Number of symmetry equations.
  for (const auto& dotsp : sym_spec4sect.expand()) {
    neqs += sym_spec4sect.m_topo->relations().n_syms(sym_spec4sect.m_sector,
                                                     !dotsp.second) *
            dotsp_combs(Integral::np(), lines_from_sn(sym_spec4sect.m_sector),
                        dotsp.first, dotsp.second);
  }
  return neqs;
}

/***********************
 * MappingPrescription *
 ***********************/

std::size_t MappingPrescription::hash() const {
  std::vector<std::size_t> hashvec;
  hashvec.push_back(std::hash<bool>{}(m_dotonly));
  hashvec.push_back(std::hash<uint32_t>{}(m_target_topoid));
  hashvec.push_back(std::hash<uint32_t>{}(m_target_sector));
  hashvec.push_back(std::hash<std::vector<std::pair<uint32_t,
                    std::vector<std::pair<uint32_t,std::string>>>>>{}(m_smap));
  return std::hash<std::vector<std::size_t>>{}(hashvec);
}

bool operator==(const MappingPrescription &mp1, const MappingPrescription &mp2)
{
  if (mp1.m_dotonly != mp2.m_dotonly ||
      mp1.m_target_topoid != mp2.m_target_topoid ||
      mp1.m_target_sector != mp2.m_target_sector)
    return false;
  if (mp1.m_smap.size() != mp2.m_smap.size()) return false;
  for (std::size_t n = 0u; n != mp1.m_smap.size(); ++n) {
    if (mp1.m_smap[n].first != mp2.m_smap[n].first) return false;
    if (mp1.m_smap[n].second.size() != mp2.m_smap[n].second.size())
      return false;
    for (std::size_t m = 0u; m != mp1.m_smap[n].second.size(); ++m) {
      if (mp1.m_smap[n].second[m].first != mp2.m_smap[n].second[m].first)
        return false;
      if (parse_coeff<Coeff_int>(mp1.m_smap[n].second[m].second) !=
          parse_coeff<Coeff_int>(mp2.m_smap[n].second[m].second))
        return false;
    }
  }
  return true;
}

bool operator==(const SectorMappingPrescription &smp1,
                const SectorMappingPrescription &smp2) {
  std::unordered_map<uint32_t, std::unordered_set<MappingPrescription>>
    smap1, smap2;
  for (const auto &from_sector_mps: smp1.m_smap) {
    smap1.insert({from_sector_mps.first,
      std::unordered_set<MappingPrescription>(from_sector_mps.second.begin(),
                                              from_sector_mps.second.end())});
  }
  for (const auto &from_sector_mps: smp2.m_smap) {
    smap2.insert({from_sector_mps.first,
      std::unordered_set<MappingPrescription>(from_sector_mps.second.begin(),
                                              from_sector_mps.second.end())});
  }
  return (smap1 == smap2);
}

/*********************
 * IntegralRelations *
 *********************/

IntegralRelations::IntegralRelations(uint32_t np, uint32_t topoid)
    : m_topoid{topoid}, m_np{np} {
}

IntegralRelations::IntegralRelations(const std::string& toponame, uint32_t np,
                                     uint32_t topoid)
    : m_topoid{topoid}, m_np{np} {
  if (Integral::ibpdir().empty()) {
    throw init_error("IntegralRelations constructor: ibpdir was not set.");
  }
  import_sym(Integral::ibpdir() + "/" + toponame + "/symmetries");
  import_sym(Integral::ibpdir() + "/" + toponame + "/relations");
  import_ibp(Integral::ibpdir() + "/" + toponame + "/IBP");
  if (Integral::get_use_li()) {
    import_ibp(Integral::ibpdir() + "/" + toponame + "/LI");
  }
}

uint32_t IntegralRelations::n_syms(uint32_t sect, const bool dotonly) const {
  // Return the number of symmetry relations of sector 'sect'.
  // If 'dotonly', all symetries are counted,
  // otherwise dot-only symmetries are skipped
  // (i.e. 'dotonly' means the regarded seed has no scalar products).
  uint32_t nsyms{0};
  auto mpresc_it = symmetries().find(sect);
  if (mpresc_it != symmetries().end()) {
    if (dotonly) {
      nsyms = static_cast<uint32_t>(mpresc_it->second.size());
    }
    else {
      for (const auto& mp : mpresc_it->second) {
        if (!mp.m_dotonly) {
          ++nsyms;
        }
      }
    }
  }
  return nsyms;
}

std::size_t IntegralRelations::cache_size() {
#define PYRED_PP_RELATIONSCACHESIZE(k) +cache<PYRED_PP_COEFFCLASS(k)>().size()
  return (
      0 PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_RELATIONSCACHESIZE));
}

void IntegralRelations::import_ibp(const std::string& fn) {
  std::ifstream instream{fn};
  if (!instream.good()) {
    std::ostringstream ss;
    ss << "Failed reading IBP file " << fn;
    throw init_error(ss.str());
  }
  is_equation eqn;
  std::string line;
  while (pyred::getline_normalized(instream, line)) {
    if (line.empty() || line.front() == '0') {
      if (!eqn.empty()) {
        m_ibps.push_back(std::move(eqn));
        eqn.clear();
      }
    }
    else {
      eqn.push_back(import_eqcf(line, m_topoid));
    }
  }
  if (!eqn.empty()) {
    m_ibps.push_back(std::move(eqn));
  }
}

SectorMappingPrescription IntegralRelations::import_sym(const std::string &fn,
                                                        uint32_t from_topoid,
                                                        uint32_t np) {
  // Parse sector mapping prescriptions.
  // Notation:
  //   0 # dot-only symmetry? (step 1)
  //   ​2 9 # map from_topoid 2, from_sector 9 (step 2)
  //   0 9 # to to_topoid 0, to_sector 9 (step 3)
  //   1  0  0  0  0  0  0  0 # D_0; (step 4 from here on)
  //   0  1 -1  1 -1  0  1  0 # D_1; mapping of the Propagator D_n
  //   0  0  1  0  0  0  0  0 # D_2; of from_topoid to to_topoid
  //   0  0  0  1  0  0  0  0 # ...
  //   0  0  0  1 -1  1  0 -s
  //   1 -1  1  0  0  0  0 -s
  //   0  0  0  0  0  0  1  0 # D_{np-1}
  //
  // Propagators of the target topology: D'_j, j=1,...,n-1.
  // Introduce D'_{n} = 1
  // D_k mapping "f_{k,0} f_{k,1} ... f_{k,n-1} f_{k,n}"
  // means D_k -> \sum_{j=0}^{n} f_{k,j}*D'_j
  SectorMappingPrescription symmetries;
  std::ifstream instream{fn};
  if (!instream.good()) {
    std::ostringstream ss;
    ss << "Failed reading symmetry file " << fn;
    throw init_error(ss.str());
  }
  int step{1};
  uint32_t from_den{0};
  uint32_t linenumber{0};
  int dotonly{0};
  uint32_t from_sector{0};
  uint32_t to_topoid{0};
  uint32_t to_sector{0};
  std::vector<
      std::pair<uint32_t, std::vector<std::pair<uint32_t, std::string>>>>
      smap;
  std::string line;
  bool finished = false;
  while (!finished) {
    if (!pyred::getline_normalized(instream, line)) {
      // One last pass to add the last smap to symmetries.
      line = "";
      finished = true;
    }
    ++linenumber;
    auto tokens = split(line, ' ');
    if (tokens.empty()) {
      if (!smap.empty()) {
        // insert smap into symmetries
        if (step != 4) {
          std::ostringstream ss;
          ss << fn << ":" << linenumber << " Invalid symmetry file format.";
          throw parser_error(ss.str());
        }
        if (smap.size() != np) {
          std::ostringstream ss;
          ss << fn << ":" << linenumber << " Wrong number of mapped "
             << "propagators. Expected " << np << " for consistency, "
             << "found " << smap.size() << ".";
          throw parser_error(ss.str());
        }
        // Sort smap such that propataor mappings with fewer terms come earlier.
        std::sort(
            smap.begin(), smap.end(),
            [](decltype(smap)::value_type& a, decltype(smap)::value_type& b) {
              return a.second.size() < b.second.size();
            });
        auto it = symmetries.m_smap.insert({from_sector, {}}).first;
        it->second.push_back(MappingPrescription{
            static_cast<bool>(dotonly), to_topoid, to_sector, std::move(smap)});
        smap.clear();
        step = 1;
        from_den = 0;
      }
    }
    else if (tokens.front().front() != '#') {
      // is not a comment
      if (step == 1) {
        if (tokens.size() != 1) {
          std::ostringstream ss;
          ss << fn << ":" << linenumber << " Invalid symmetry file format."
             << " Expected one number, found " << tokens.size() << ".";
          throw parser_error(ss.str());
        }
        try {
          dotonly = parse_int<int>(tokens[0]);
        }
        catch (const parser_error& err) {
          std::ostringstream ss;
          ss << fn << ":" << linenumber << " " << err.message();
          throw parser_error(ss.str());
        }
        ++step;
      }
      else if (step == 2 || step == 3) {
        if (tokens.size() != 2) {
          std::ostringstream ss;
          ss << fn << ":" << linenumber << " Invalid symmetry file format."
             << " Expected two numbers, found " << tokens.size() << ".";
          throw parser_error(ss.str());
        }
        uint32_t topoid;
        uint32_t sector;
        try {
          topoid = parse_int<uint32_t>(tokens[0]);
          sector = parse_int<uint32_t>(tokens[1]);
        }
        catch (const parser_error& err) {
          std::ostringstream ss;
          ss << fn << ":" << linenumber << " " << err.message();
          throw parser_error(ss.str());
        }
        if (step == 2) {
          // topoid means from_topoid
          if (topoid != from_topoid) {
            std::ostringstream ss;
            ss << fn << ":" << linenumber << " Mapping prescription does "
               << "not match the topology for which it is imported. Is "
               << topoid << ", expected " << from_topoid << ".";
            throw init_error(ss.str());
          }
          from_sector = sector;
        }
        else {
          // TODO: check if to_topoid (and from_sector and to_sector) exists.
          to_topoid = topoid;
          to_sector = sector;
        }
        ++step;
      }
      else {
        // step 4
        if (tokens.size() != np + 1) {
          std::ostringstream ss;
          ss << fn << ":" << linenumber << " Invalid symmetry file format."
             << " Expected " << np + 1 << " coefficients for consistency, "
             << "found " << tokens.size() << ".";
          throw parser_error(ss.str());
        }
        std::vector<std::pair<uint32_t, std::string>> pmap;
        for (uint32_t k = 0; k != np + 1; ++k) {
          if (tokens[k] != "0") {
            pmap.emplace_back(k, tokens[k]);
          }
        }
        smap.emplace_back(from_den, std::move(pmap));
        ++from_den;
      }
    }
  }
  return symmetries;
}

void IntegralRelations::import_sym(const std::string& fn) {
  auto symmetries = import_sym(fn, m_topoid, m_np);
  for (const auto &sym: symmetries.m_smap) {
    auto it = m_symmetries.m_smap.insert({sym.first, {}}).first;
    for (const auto &mp: sym.second) {
      it->second.push_back(mp);
    }
  }
}

std::vector<eqdata> IntegralRelations::sow_ibp_streq(
    const std::vector<pow_type>& seed) const {
  std::vector<eqdata> numibps;
  std::vector<std::string> s_seed;
  s_seed.reserve(seed.size());
  for (auto p : seed) {
    s_seed.push_back(std::to_string(p));
  }
  for (const auto& ibp_tmpl : m_ibps) {
    eqdata crop;
    crop.reserve(ibp_tmpl.size());
    for (const auto& ic : ibp_tmpl) {
      auto weight = ic.first.zip(seed).to_weight();
      auto coeff = sow_coeff(ic.second, s_seed);
      if (weight != Weight::none() && coeff != "0") {
        crop.emplace_back(weight, coeff);
      }
    }
    numibps.push_back(std::move(crop));
  }
  return numibps;
}

std::vector<eqdata> IntegralRelations::sow_sym_streq(
    const std::vector<pow_type>& seed) const {
  // Symmetrise an integral with propagator powers 'seed'.
  if (seed.size() != m_np) {
    throw input_error(
        "IntegralRelations::sow_sym_streq(): wrong number of seed powers.");
  }
  // Calculate the sector number and the number of scalar products
  // from the seed.
  uint32_t sector{0};
  uint32_t sps{0};
  {
    uint32_t linepos{0};
    for (const auto p : seed) {
      if (p > 0) {
        sector |= (1 << linepos);
      }
      else {
        sps -= p;
      }
      ++linepos;
    }
  }
  {
    // Check if the seed integral belongs to a trivial sector.
    const auto& zs = Integral::s_trivialsectors.at(m_topoid);
    if (zs.count(sector)) {
      return {};
    }
  }
  const auto mapping_presc_it = symmetries().find(sector);
  if (mapping_presc_it == symmetries().end()) {
    // No mapping for this sector is known: return zero equations.
    // std::cout << "sow_sym_streq(): no mapping for sector "
    //           << sector << "." << std::endl;
    return {};
  }
  const auto& mapping_presc_vec = mapping_presc_it->second;
  std::vector<eqdata> symmeqs;
  symmeqs.reserve(mapping_presc_vec.size());
  for (const auto& mapping_presc : mapping_presc_vec) {
    // Loop over all mapping prescriptions for 'sector' to
    // generate one equation for each mapping.
    if (mapping_presc.m_dotonly && sps) {
      continue;
    }
    const auto mapped_topoid = mapping_presc.m_target_topoid;
    const auto mapped_sector = mapping_presc.m_target_sector;
    const auto& smap = mapping_presc.m_smap;
    const auto& zero_sects = Integral::s_trivialsectors.at(mapped_topoid);
    auto mapped_perm = std::vector<pow_type>(m_np, 0);
    auto fac_perm = std::string{"1"};
    bool permute{true};
    // temporary equation: map {<{powers}, sector>: coefficient, ...}
    std::unordered_map<std::pair<std::vector<pow_type>, uint32_t>, std::string>
        symmeq, next_mapped;
    for (const auto& pmap : smap) {
      const auto oldpos = pmap.first;
      const auto& pshifts = pmap.second;
      if (pshifts.size() == 1) {
        // These are always grouped at the beginning of smap.
        mapped_perm[pshifts.front().first] = seed[oldpos];
        if (pshifts.front().second != "1" && seed[oldpos]) {
          fac_perm = fac_perm + "*(" + pshifts.front().second + ")";
          if (seed[oldpos] != 1) {
            fac_perm = fac_perm + "^" + std::to_string(-seed[oldpos]);
          }
        }
        continue;
      }
      if (permute) {
        // From now on all pmap have pshifts.size()>1,
        // i.e. the permutation part is already done.
        // Convert the permuted integral to an expression.
        // symmeq is the right hand side of the equation seed==symmeq.
        symmeq.emplace(
            std::make_pair(std::move(mapped_perm), uint32_t{mapped_sector}),
            std::move(fac_perm));
        permute = false;
      }
      if (seed[oldpos] > 0) {
        std::ostringstream ss;
        ss << "IntegralRelations::sow_sym_streq(): invalid mapping of ("
           << m_topoid << "," << sector << ") -> (" << mapped_topoid << ","
           << mapped_sector << "): encountered a linear combination of "
           << "propagators in the denominator.";
        throw init_error(ss.str());
      }
      // Map scalar products which are not trivially permuted.
      for (pow_type ppow = 0; ppow != seed[oldpos]; --ppow) {
        next_mapped.clear();
        for (const auto& pshift : pshifts) {
          const auto shiftpos = pshift.first;
          const auto& shiftfac = pshift.second;
          for (const auto& fi_fac : symmeq) {
            const auto& fi = fi_fac.first;
            const auto& fac = fi_fac.second;
            auto mm_sn = fi.second;
            auto mm_powers = fi.first;
            if (shiftpos != m_np) {
              // shiftpos==m_np is a pseudo-propagator D_{m_np}==1.
              if (mm_powers[shiftpos] == 1) {
                // Pinching the line 'shiftpos' lowers the sector.
                mm_sn &= ~(1 << shiftpos);
                if (zero_sects.count(mm_sn)) {
                  // The lowered sector is trivial.
                  continue;
                }
              }
              --mm_powers[shiftpos];
            }
            auto inserted = next_mapped.emplace(
                std::make_pair(std::move(mm_powers), std::move(mm_sn)),
                std::string{});
            inserted.first->second += "+(" + fac + ")*(" + shiftfac + ")";
          }
        }
        symmeq = std::move(next_mapped);
      }
    }
    if (permute) {
      // If the mapping is a pure permutation,
      // still need to convert the permuted integral to an expression.
      symmeq.emplace(
          std::make_pair(std::move(mapped_perm), uint32_t{mapped_sector}),
          std::move(fac_perm));
    }
    eqdata eq;
    // Add the left-hand side integral.
    if (mapped_topoid == m_topoid) {
      // Symmetry equation:
      // all integrals belong to the same topology 'mapped_topoid'.
      auto inserted =
          symmeq.emplace(std::make_pair(seed, sector), std::string{});
      inserted.first->second += "-1";
      eq.reserve(symmeq.size());
    }
    else {
      // The left-hand side integral is the only one which belongs
      // to the topology 'm_topoid'.
      eq.reserve(symmeq.size() + 1);
      eq.emplace_back(Integral(m_topoid, seed).to_weight(), std::string("-1"));
    }
    // Convert symmeq into a weight-string equation.
    for (auto& powsec_fac : symmeq) {
      auto w = Integral(mapped_topoid, powsec_fac.first.first).to_weight();
      assert(w != Weight::none());
      eq.emplace_back(std::move(w), std::move(powsec_fac.second));
    }
    symmeqs.push_back(std::move(eq));
    // Do not sort here, but in IntegralRelations::parse_streq().
  }
  return symmeqs;
}

void delete_IntegralRelations(IntegralRelations* rel) { delete rel; }

relations_ptr import_relations(const std::string& toponame, uint32_t np,
                               uint32_t topoid) {
  // make_unique cannot be used with a custom deleter.
  return {new IntegralRelations(toponame, np, topoid),
          delete_IntegralRelations};
}

uint32_t get_relations_topoid(const relations_ptr& rel) {
  return rel->topoid();
}

std::vector<std::vector<std::pair<uint32_t, uint32_t>>>
GeneratorHelper::generate_ts_bunches(
    const GeneratorHelper::specmap_t& specmap) {
  // Create bunches of (topo,sector).
  // Vector of {(t,s1),(t,s2),...} where all s_i have the same number of lines.
  // First sorted by topo, then by number of lines, then by sector ordering.
  //
  // Strategy: each bunch contains all sectors with a specific number of lines
  // from a specific topology.
  // If the generator with "global sorting", i.e. sorting equations across
  // several seed sectors, is used, this could be improved to create smaller
  // bunches (to save memory).
  //
  // Special case sector_ordering=1: no grouping by lines.
  // Return one bunch with all sectors in order.
  if (Integral::sector_ordering() == 1) {
    std::vector<std::pair<uint32_t, uint32_t>> ts_vect;
    ts_vect.reserve(specmap.size());
    for (const auto& iss : specmap) {
      ts_vect.push_back(iss.first);
    }
    // Sort by (topoid,sector_weight)
    std::sort(ts_vect.begin(), ts_vect.end(),
              [](const std::pair<uint32_t, uint32_t>& a,
                 const std::pair<uint32_t, uint32_t>& b) {
                if (a.first == b.first) { // topo equal
                  // compare sector weight
                  return (Integral::sector_to_sectweight(a.first, a.second) <
                          Integral::sector_to_sectweight(b.first, b.second));
                }
                else { // topo different
                  return a.first < b.first;
                }
              });
    return {std::move(ts_vect)};
  }
  // All other sector orderings, i.e. those where sectors are first ordered
  // by number of lines.
  // Group sectors by (topoid,#lines) in a map.
  std::unordered_map<std::pair<uint32_t, uint32_t>, std::vector<uint32_t>>
      tl_sects;
  for (const auto& iss : specmap) {
    auto lines = lines_from_sn(iss.first.second);
    auto it = tl_sects.insert({{iss.first.first, lines}, {}}).first;
    it->second.push_back(iss.first.second);
  }
  // Copy the map to a pair(key,value) vector
  // and sort the sector lists.
  std::vector<std::pair<std::pair<uint32_t, uint32_t>, std::vector<uint32_t>>>
      tls_vect;
  tls_vect.reserve(tl_sects.size());
  for (auto& tl_s : tl_sects) {
    tls_vect.push_back(std::move(tl_s));
    auto& svect = tls_vect.back().second;
    auto topo = tl_s.first.first;
    // Sort sectors by sector weight
    std::sort(svect.begin(), svect.end(),
              [&topo](const uint32_t& a, const uint32_t& b) {
                return (Integral::sector_to_sectweight(topo, a) <
                        Integral::sector_to_sectweight(topo, b));
              });
  }
  // Sort by (topoid,lines)
  std::sort(
      tls_vect.begin(), tls_vect.end(),
      [](const std::pair<std::pair<uint32_t, uint32_t>, std::vector<uint32_t>>&
             a,
         const std::pair<std::pair<uint32_t, uint32_t>, std::vector<uint32_t>>&
             b) {
        if (a.first.first == b.first.first) {     // topo equal
          return a.first.second < b.first.second; // compare lines
        }
        else { // topo different
          return a.first.first < b.first.first;
        }
      });
  // Pair (topoid,sector) and assign the bunches.
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> bunches;
  bunches.reserve(tls_vect.size());
  for (const auto& tl_s : tls_vect) {
    std::vector<std::pair<uint32_t, uint32_t>> ts_vect;
    ts_vect.reserve(tl_s.second.size());
    for (const auto& s : tl_s.second) {
      ts_vect.push_back({tl_s.first.first, s});
    }
    bunches.push_back({std::move(ts_vect), {}});
  }
  return bunches;
}

std::pair<GeneratorHelper::gen_sol_bunches_t, std::size_t>
GeneratorHelper::GeneratorParallelMaster_init(const specmap_t& specmap) {
  // Used to generate the 'to solve' part of
  // GeneratorParallelMaster::m_bunches
  // and the maximal size of the system.
  GeneratorHelper::gen_sol_bunches_t bunches;
  {
    auto ts_bunches = generate_ts_bunches(specmap);
    bunches.reserve(ts_bunches.size());
    for (auto& bunch : ts_bunches) {
      bunches.emplace_back(std::move(bunch),
                           std::vector<std::pair<uint32_t, uint32_t>>{});
    }
  }
  // Three "lookahead" strategies.
  if (Config::lookahead() <= 0) {
    for (auto& b : bunches) {
      // Schedule the sectors to solve right after their bunch
      // has been generated.
      b.second = b.first;
    }
  }
  else if (Config::lookahead() == 1) {
    // Schedule the sectors to solve after the next bunch
    // has been generated.
    for (auto bit = bunches.begin(); bit != bunches.end() - 1; ++bit) {
      (bit + 1)->second = bit->first;
    }
    // Add the sectors from the last bunch to its own 'tosolve'.
    for (const auto& ts : bunches.back().first) {
      bunches.back().second.push_back(ts);
    }
  }
  else {
    // lookahead == 2: solve after everything has been generated.
    std::vector<std::pair<uint32_t, uint32_t>> ts_vect;
    for (const auto& b : bunches) {
      for (const auto& ts : b.first) {
        ts_vect.push_back(ts);
      }
    }
    bunches.back().second = std::move(ts_vect);
    // Sort sectors by sector number (not by number of lines).
    //     std::sort(ts_vect.begin(), ts_vect.end(),
    //               [](const std::pair<uint32_t,uint32_t> &a,
    //                 const std::pair<uint32_t,uint32_t> &b) {
    //                 return a.second < b.second;
    //               });
    //     bunches.clear();
    //     bunches.push_back({ts_vect,ts_vect});
  }
  // Done creating bunches.
  // Determine the maximal system size
  // (i.e. the size in case no equation vanishes).
  std::size_t maxsyssz = Integral::get_preferred_masters().n_basis_lcs +
                         Integral::get_external_equations().eqs.size();
  for (const auto& b : bunches) {
    for (const auto& ts : b.first) {
      auto neqs4sect = neqs_from_seed(specmap.at(ts));
      maxsyssz += neqs4sect;
    }
  }
  return {std::move(bunches), std::move(maxsyssz)};
}

GeneratorHelper::specmap_t GeneratorHelper::seedspecs_to_specmap(
    const std::vector<SeedSpec>& ibp_seedspec,
    const std::vector<SeedSpec>& ibp_seedcompl,
    const std::vector<SeedSpec>& sym_seedspec, const bool sym_autoseed) {
  // sym_autoseed
  //   If sym_autoseed=true and sym_seedspec is empty:
  //   automatically select the symmetry seed by cutting the ibp seed
  //   at Config::s_symlimits. If sym_seedspec is not empty,
  //   always use it as it is.
  auto ibp_spec4sects = SeedSpec::expand_sector(ibp_seedspec, ibp_seedcompl);
  auto sym_spec4sects = SeedSpec::expand_sector(sym_seedspec);
  // Note that there is no overlap in 'ibp_spec4sects',
  // i.e. there is only one SeedSpecSector for any (topoid,sector).
  // Also, they are ordered by (topoid,sector).
  if (sym_autoseed && sym_seedspec.empty()) {
    // Auto symmetry seed selection: use the ibp seeds without complements
    // and cut at Config::s_symlimits.
    sym_spec4sects = SeedSpec::expand_sector(ibp_seedspec);
    for (auto& ss : sym_spec4sects) {
      ss.cut();
    }
  }
  // Group seed specifications for ibps and symmetries into
  // pairs -- one pair per (topoid,sector). Store them in a map
  // with (topoid,sector) as keys.
  specmap_t ibp_sym_spec_map;
  auto ibp_spec_it = ibp_spec4sects.cbegin();
  auto sym_spec_it = sym_spec4sects.cbegin();
  while (ibp_spec_it != ibp_spec4sects.cend() ||
         sym_spec_it != sym_spec4sects.cend()) {
    if (ibp_spec_it == ibp_spec4sects.cend() ||
        (sym_spec_it != sym_spec4sects.cend() && *sym_spec_it < *ibp_spec_it)) {
      // ibp_spec empty
      ibp_sym_spec_map.insert(
          {{sym_spec_it->m_topo->m_id, sym_spec_it->m_sector},
           {SeedSpecSector(*sym_spec_it->m_topo, sym_spec_it->m_sector),
            *sym_spec_it}});
      ++sym_spec_it;
    }
    else if (sym_spec_it == sym_spec4sects.cend() ||
             (ibp_spec_it != ibp_spec4sects.cend() &&
              *ibp_spec_it < *sym_spec_it)) {
      // sym_spec empty
      ibp_sym_spec_map.insert(
          {{ibp_spec_it->m_topo->m_id, ibp_spec_it->m_sector},
           {*ibp_spec_it,
            SeedSpecSector(*ibp_spec_it->m_topo, ibp_spec_it->m_sector)}});
      ++ibp_spec_it;
    }
    else {
      // ibp_spec and sym_spec have the same (topoid,sector).
      ibp_sym_spec_map.insert(
          {{ibp_spec_it->m_topo->m_id, ibp_spec_it->m_sector},
           {*ibp_spec_it, *sym_spec_it}});
      ++ibp_spec_it;
      ++sym_spec_it;
    }
  }
  return ibp_sym_spec_map;
}

} // namespace pyred
