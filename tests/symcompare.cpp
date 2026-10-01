
// Usage: symcompare --mode=1|2|3 /path/to/project_dir /path/to/reference_dir job_id
// Reference results are taken from results_<job_id>

// Compare
// mode >= 1:
//   sectormappings/TOPONAME/trivialsector
//   sectormappings/TOPONAME/symmetries
//   sectormappings/TOPONAME/relations
//   TODO: compare IBP and LI.
// mode >= 2:
//   results/TOPONAME/masters
// mode >= 3:
//   results/TOPONAME/masters.final

#include <string>
#include <iostream>
#include <vector>
#include <unordered_set>

#include "pyred/integrals.h"
#include "pyred/relations.h"


int print_check(bool test, const std::string &name) {
  int error = 0;
  std::cout << " " << name;
  if (test) {
    std::cout << " ok";
  }
  else {
    ++error;
    std::cout << " MISMATCH";
  }
  return error;
}


int main(int argc, char *argv[]) {
  int mode = 3; // 1: sym; 2: masters; 3: masters.final
  std::string chkcfg;
  std::string refcfg;
  std::string job_id;
  int posarg = 0;
  for (int narg = 1; narg < argc; ++narg) {
    auto arg = std::string(argv[narg]);
    if (arg == "--mode=1") { mode = 1; }
    else if (arg == "--mode=2") { mode = 2; }
    else if (arg == "--mode=3") { mode = 3; }
    else {
      ++posarg;
      if (posarg == 1) { chkcfg = arg; }
      else if (posarg == 2) { refcfg = arg; }
      else if (posarg == 3) { job_id = arg; }
    }
  }
  if (posarg != 3) {
    std::cout << ("Exactly two positional arguments are required: "
      "the Kira config dirs to compare.\n"
      "The only other allowed argument is --mode=1|2|3.") << std::endl;
    return 1;
  }
  std::cout << "Compare Kira topology data" << std::endl;
  std::cout << "  " << chkcfg << std::endl;
  std::cout << "  " << refcfg << std::endl;

  auto topoorder_chk = pyred::Topology::parse_topology_ordering(
    chkcfg + "/sectormappings/topology_ordering");
  auto topoorder_ref = pyred::Topology::parse_topology_ordering(
    refcfg + "/sectormappings/topology_ordering");
  if (topoorder_chk != topoorder_ref) {
    std::cout << "topology_ordering doesn't match." << std::endl;
    return 1;
  }
//   for (const auto &toponame_np: topoorder_ref) {
//     pyred::new_topology(toponame_np.first,
//                         toponame_np.second,
//                         {(static_cast<uint32_t>(1) << toponame_np.second) - 1});
//   }

  int error = 0;
  for (uint32_t topoid = 0u; topoid != topoorder_chk.size(); ++topoid) {
    auto toponame = topoorder_chk[topoid].first;
    auto np = topoorder_chk[topoid].second;

    auto ts_chk = pyred::Topology::import_trivialsectors(
      chkcfg + "/sectormappings/" + toponame + "/trivialsector");
    auto ts_ref = pyred::Topology::import_trivialsectors(
      refcfg + "/sectormappings/" + toponame + "/trivialsector");

    // symmetries
    auto sym_chk = pyred::IntegralRelations::import_sym(
      chkcfg + "/sectormappings/" + toponame + "/symmetries", topoid, np);
    auto sym_ref = pyred::IntegralRelations::import_sym(
      refcfg + "/sectormappings/" + toponame + "/symmetries", topoid, np);

    // mappings
    auto map_chk = pyred::IntegralRelations::import_sym(
      chkcfg + "/sectormappings/" + toponame + "/relations", topoid, np);
    auto map_ref = pyred::IntegralRelations::import_sym(
      refcfg + "/sectormappings/" + toponame + "/relations", topoid, np);

    // masters
    bool masters_success = false;
    if (mode >= 2) {
      auto chkmasters = pyred::TopoConfigData::import_integrals(
        chkcfg + "/results/" + toponame + "/masters");
      auto refmasters = pyred::TopoConfigData::import_integrals(
        refcfg + "/results_" + job_id + "/" + toponame + "/masters");
      std::unordered_set<pyred::Integral> chkmasters_set(
        chkmasters.begin(), chkmasters.end());
      std::unordered_set<pyred::Integral> refmasters_set(
        refmasters.begin(), refmasters.end());
      masters_success = (chkmasters_set == refmasters_set);
    }

    // masters.final
    bool mastersfinal_success = false;
    if (mode >= 3) {
      auto chkmasters = pyred::TopoConfigData::import_integrals(
        chkcfg + "/results/" + toponame + "/masters.final");
      auto refmasters = pyred::TopoConfigData::import_integrals(
        refcfg + "/results_" + job_id + "/" + toponame + "/masters.final");
      std::unordered_set<pyred::Integral> chkmasters_set(
        chkmasters.begin(), chkmasters.end());
      std::unordered_set<pyred::Integral> refmasters_set(
        refmasters.begin(), refmasters.end());
      mastersfinal_success = (chkmasters_set == refmasters_set);
    }

    std::cout << "Topology " << toponame << ":";
    error += print_check(ts_chk == ts_ref, "trivialsectors");
    error += print_check(sym_chk == sym_ref, "sym");
    error += print_check(map_chk == map_ref, "map");
    if (mode >= 2) error += print_check(masters_success, "masters");
    if (mode >= 3) error += print_check(mastersfinal_success, "masters.final");
    std::cout << std::endl;
  }
  return error ? 1 : 0;
}
