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

#include <unistd.h> // getopt
#include <sys/stat.h>
#include <dirent.h>

#include <fstream>
#include <stack>
#include <stdexcept>

#include "pyred/interface.h"
#include "pyred/parser.h"
#include "gzstream/gzstream.h"


namespace pyred {

float time_diff(const std::clock_t& start, const std::clock_t& end) {
  return (float(end - start) / CLOCKS_PER_SEC);
}

float time_diff(const std::chrono::time_point<std::chrono::system_clock>& start,
                const std::chrono::time_point<std::chrono::system_clock>& end) {
  return static_cast<std::chrono::duration<float>>(end - start).count();
}

/***************************
 * Logger and LockedLogger *
 ***************************/

Logger::Logger(int cout_verb, int file_verb)
    : m_verbosity{std::make_pair(cout_verb, file_verb)}, m_logfilestream{} {}

std::pair<int, int> Logger::verbosity(int cout_verb, int file_verb) {
  // The lowest verbosity is -1 (suppress all output).
  // file_verb=-2 (default) sets file_verb=cout_verb.
  if (cout_verb < -1) {
    throw input_error("Verbosity level <-1 is not allowed.");
  }
  m_verbosity.first = cout_verb;
  m_verbosity.second = cout_verb;
  if (file_verb >= -1) m_verbosity.second = file_verb;
  return m_verbosity;
}

void Logger::attach_logfile(const std::string& fn, const FileOpenMode mode) {
  if (mode == FileOpenMode::create && keyvaluedb::KeyValueDB::file_exists(fn)) {
    throw input_error("Log file " + fn + " already exists.");
  }
  auto omode = std::ios::trunc;
  if (mode == FileOpenMode::append) {
    omode = std::ios::app;
  }
  m_logfilestream = std::ofstream(fn, omode);
  if (!m_logfilestream.good()) {
    throw input_error("Log file " + fn + " is not writable.");
  }
}

void Logger::detach_logfile() {
  if (m_logfilestream.is_open()) {
    m_logfilestream.close();
  }
}

LockedLogger Logger::operator()(int cout_lev, int file_lev) {
  return LockedLogger(*this, cout_lev, file_lev);
}

LockedLogger::LockedLogger(Logger& lg, int cout_lev, int file_lev)
    : m_verbosity{lg.m_verbosity},
      m_level{std::make_pair(cout_lev, file_lev)},
      m_logfilestream{lg.m_logfilestream},
      m_lck{lg.m_mtx} {
  if (file_lev < 0) m_level.second = cout_lev;
  if (m_level.second < 0) {
    throw input_error("Logger: log level <0 is not allowed.");
  }
  std::cout << std::flush;
}

LockedLogger& LockedLogger::operator<<(std::ostream& (*omanip)(std::ostream&)) {
  if (m_level.first <= m_verbosity.first) std::cout << (*omanip);
  if (m_level.second <= m_verbosity.second && m_logfilestream.is_open()) {
    m_logfilestream << (*omanip);
  }
  return *this;
}

/**********
 * Config *
 **********/

Logger Config::s_log_val{1};
int Config::s_coeff_cls_val{1};
bool Config::s_backward_val{true};
int Config::s_parallel{1};
bool Config::s_auto_symseed{true};
std::pair<int, int> Config::s_symlimits{std::numeric_limits<int>::max(), 2};
int Config::s_lookahead{1};
int Config::s_insertion_tracer{1};
std::pair<std::string, bool> Config::s_database_file{"./insertions", false};
bool Config::s_finished{false};
bool Config::s_johanntrick{false};

void Config::coeff_cls(int ccls) {
  if (ccls < 1 || ccls > PYRED_PP_NCOEFFCLASSES) {
    std::ostringstream ss;
    ss << "Invalid coefficient class number (must be 1..."
       << PYRED_PP_NCOEFFCLASSES << "): " << ccls;
    throw init_error(ss.str());
  }
  s_coeff_cls_val = ccls;
}

void Config::backward(bool bkwd) { s_backward_val = bkwd; }

std::pair<int, int> Config::verbosity(int verb, int file_verb) {
  return s_log_val.verbosity(verb, file_verb);
}

void Config::parallel(int n) {
  if (n < 0) n = 1;
  s_parallel = n;
}

void Config::auto_symseed(const bool ass) { s_auto_symseed = ass; }

void Config::symlimits(const int maxdots, const int maxsps) {
  s_symlimits = {maxdots, maxsps};
}

void Config::lookahead(int lah) {
  if (lah < -1 || lah > 2) {
    throw init_error("Config::lookahead may only be initialised to -1,0,1,2.");
  }
  s_lookahead = lah;
}

void Config::insertion_tracer(int mode) {
  if (mode < 0 || mode > 4) {
    throw init_error("Config::insertion_tracer may only be initialised to "
                     "0 (off), 1 (in memory, default), 2 (SQLite), "
                     "3 (Kyoto Cabinet), "
                     "4 (Kyoto Cabinet if available, otherwise SQLite).");
  }
#ifdef PYRED_KCDB
  if (mode == 4) mode = 3;
#else
  if (mode == 4) mode = 2;
  if (mode == 3) {
    throw init_error("To use Config::insertion_tracer=3, pyRed must be built "
                     "with Kyoto Cabinet database support.");
  }
#endif
  s_insertion_tracer = mode;
}

std::pair<std::string, std::string> Config::parse(const std::string& conf) {
  auto tmp = split(conf, ' ');
  tmp.insert(tmp.begin(), ""); // getopt() ignores the first element
  return parse(tmp);
}

std::pair<std::string, std::string> Config::parse(
    const std::vector<std::string>& conf) {
  if (s_finished) {
    throw init_error(
        "Config::parse() error: finalised config cannot be changed anymore.");
  }
  // convert conf to char**
  int argc = conf.size();
  std::vector<char*> argv;
  argv.reserve(argc);
  for (const auto& opt : conf) {
    argv.push_back(const_cast<char*>(opt.c_str()));
  }
  int opt;
  std::string infile{};
  std::string outfile{};
  optind = 1;
  while ((opt = getopt(argc, argv.data(), "1:c:i:o:s:b:t:v:")) != -1) {
    /*
     * -1:string symbol to set to one
     * -c:int use coefficient array of this size
     * -i:string declare symbol (comma separated or multiple -i options)
     * -o:string output file name
     * -s:int integer seed for random number generator
     * -b:int 0/1 to turn off/on backward insertion (default=1)
     */
    int seed;
    switch (opt) {
      case '1':
        CoeffHelper::settoone(std::string(optarg));
        break;
      case 'c':
        try {
          CoeffHelper::coeff_n(string_to_int(std::string(optarg)));
        }
        catch (std::invalid_argument& e) {
          std::ostringstream ss;
          ss << "Option c argument must be integer, but is " << optarg
             << std::endl;
          throw init_error(ss.str());
        }
        coeff_cls(2);
        break;
      case 'i':
        for (const auto& var : split(std::string(optarg), ',')) {
          CoeffHelper::add_invariant(var);
        }
        break;
      case 'o':
        outfile = optarg;
        break;
      case 's':
        try {
          seed = string_to_int(std::string(optarg));
        }
        catch (std::invalid_argument& e) {
          std::ostringstream ss;
          ss << "Seed argument must be integer, but is " << optarg << std::endl;
          throw init_error(ss.str());
        }
        CoeffHelper::random_seed(seed);
        break;
      case 'b':
        if (std::string(optarg) == "0") {
          backward(false);
        }
        else if (std::string(optarg) == "1") {
          backward(true);
        }
        else {
          std::ostringstream ss;
          ss << "Option -b must either be 0 or 1, but is " << optarg
             << std::endl;
          throw init_error(ss.str());
        }
        break;
      case 'v':
        try {
          verbosity(string_to_int(std::string(optarg)));
        }
        catch (std::invalid_argument& e) {
          std::ostringstream ss;
          ss << "Verbosity argument must be integer, but is " << optarg
             << std::endl;
          throw init_error(ss.str());
        }
        break;
      default:
        std::ostringstream ss;
        ss << "Invalid argument " << opt << std::endl;
        throw init_error(ss.str());
    }
  }
  if (optind == argc - 1) {
    infile = argv[optind];
  }
  else if (optind < argc - 1) {
    std::ostringstream ss;
    ss << "Too many arguments" << std::endl;
    throw init_error(ss.str());
  }
  return {infile, outfile};
}

void Config::finish() {
#define PYRED_PP_CCS_CONFIGFINISH(k)                    \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls() == k + 1) { \
    PYRED_PP_COEFFCLASS(k)::init();                     \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_CONFIGFINISH)
  s_finished = true;
}

/**************
 * FileSystem *
 **************/

// FileType::regular if fname is the name of a regular file,
// FileType::directory for a directory, FileType::other for other files.
// Throw input_error if the file doesn't exist.
FileSystem::FileType FileSystem::get_filetype(const std::string &fname) {
  struct stat st;
  auto ftype = FileType::other;
  if (stat(fname.c_str(), &st) == 0) {
    if (st.st_mode & S_IFREG) {
      ftype = FileType::regular;
    }
    else if (st.st_mode & S_IFDIR) {
      ftype = FileType::directory;
    }
  }
  else {
    throw input_error(std::string(
      "ERROR in get_filetype(): could not read file \"") + fname + "\".");
  }
  return ftype;
}

// If file_or_dir is a regular file, return a vector with it as the only
// element. If it is a directory, return a sorted vector of all names of
// regular files in it, including the directory name as prefix,
// excluding files starting with '.'.
// If fileexts is non-empty given, select only files from the directory that end
// with one of the fileexts. Ignored if file_or_dir is a regular file.
std::vector<std::string> FileSystem::get_filenames(
    const std::string &file_or_dir,
    const std::vector<std::string> &fileexts) {
  auto ftype = get_filetype(file_or_dir);
  if (ftype == FileType::regular) {
    return {file_or_dir};
  }
  else if (ftype == FileType::other) {
    throw input_error(std::string("ERROR in get_filenames(): \"") +
      file_or_dir + "\" is neither a regular file, nor a directory.");
  }
  // file_or_dir is a directory
  std::vector<std::string> filenames;
  DIR *dir;
  struct dirent *entry;
  if ((dir = opendir(file_or_dir.c_str())) != NULL) {
    while ((entry = readdir(dir)) != NULL) {
      auto fullname = file_or_dir + "/" + (entry->d_name);
      if (entry->d_name[0] != '.' &&
          get_filetype(fullname) == FileType::regular) {
        bool selectfile = fileexts.empty();
        for (const auto &fileext: fileexts) {
          if (fullname.size() >= fileext.size() &&
              fullname.compare(fullname.size() - fileext.size(),
                               fileext.size(), fileext) == 0) {
            selectfile = true;
          }
        }
        if (selectfile) {
          filenames.push_back(fullname);
        }
      }
    }
    closedir(dir);
  }
  else {
    throw input_error(std::string(
      "ERROR in get_filenames(): could not open directory \"") +
      file_or_dir + "\".");
  }
  std::sort(filenames.begin(), filenames.end());
  return filenames;
}

/**********
 * System *
 **********/

System::System()
    : sys{},
      coeff_cls{Config::coeff_cls()},
      m_content_prepared{false},
      m_solver_sorted{false} {}

System::System(const std::vector<eqdata>& eqs)
    : sys{eqs},
      coeff_cls{Config::coeff_cls()},
      m_content_prepared{false},
      m_solver_sorted{false} {}

System::System(std::vector<eqdata>&& eqs)
    : sys{std::move(eqs)},
      coeff_cls{Config::coeff_cls()},
      m_content_prepared{false},
      m_solver_sorted{false} {}

System::System(const std::string& infile)
    : sys{},
      coeff_cls{Config::coeff_cls()},
      m_content_prepared{false},
      m_solver_sorted{false}
{
  add(infile);
}

System::System(const std::vector<std::string>& infiles)
    : sys{},
      coeff_cls{Config::coeff_cls()},
      m_content_prepared{false},
      m_solver_sorted{false}
{
  add(infiles);
}

std::size_t System::size() const {
#define PYRED_PP_CCS_SYSTEMSIZE(k)                       \
  PYRED_PP_IF_TRUE(k, else) if (systemtype() == k + 1) { \
    return m_numsys.PYRED_PP_COEFFCLASSMEM(k).size();    \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMSIZE)
  else if (systemtype() == 0) {
    return sys.size();
  }
  else { // systemtype() == -1: empty system
    return 0;
  }
}

int System::systemtype(bool fatal_empty) const {
#define PYRED_PP_CCS_NUMSYSEMPTY(k) &&m_numsys.PYRED_PP_COEFFCLASSMEM(k).empty()
  if (!sys.empty() // &&
       PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_NUMSYSEMPTY)) {
    // String coefficients.
    return 0;
  }
#define PYRED_PP_CCS_NUMSYSNONEMPTY(k) \
  +(m_numsys.PYRED_PP_COEFFCLASSMEM(k).empty() ? 0 : 1)
  else if (sys.empty() && (1 == 0 // +
                           PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES,
                                            PYRED_PP_CCS_NUMSYSNONEMPTY))) {
    // Exactly one of the numsys is non-empty.
    return coeff_cls;
  }
  else if (sys.empty() // &&
           PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_NUMSYSEMPTY)) {
    if (fatal_empty) {
      throw init_error("System is empty.");
    }
    else {
      // Empty system.
      return -1;
    }
  }
  else {
    throw init_error("System mixes different coefficient types.");
  }
}

void System::sorteq(eqdata &eq) {
  // Sort equation by integral weights.
  // Remove duplicate integrals, adding their coefficients.
  std::sort(eq.begin(), eq.end(), cmp_wc_pair<std::string>);
  if (eq.size() < 2) return;
  std::size_t pos1 = 0u;
  std::size_t pos2 = 1u;
  while (pos2 != eq.size()) {
    if (eq[pos1].first == eq[pos2].first) {
      eq[pos1].second = eq[pos1].second + "+" + eq[pos2].second;
      eq.erase(eq.begin() + pos2);
      continue;
    }
    ++pos1;
    ++pos2;
  }
}

std::vector<std::size_t> System::generate_solve(
    const std::vector<SeedSpec>& ibp_seedspec,
    const std::vector<SeedSpec>& ibp_seedcompl,
    const std::vector<SeedSpec>& sym_seedspec) {
  if (Integral::sector_ordering() == 1 &&
      (Config::lookahead() == 0 || Config::lookahead() == 1)) {
    // In sector ordering 1, sectors are not ordered by number of lines.
    // Therefore the on-the-fly generator/solver which sorts equations
    // across seed sectors (i.e. lookahead != -1) can only be used with
    // lookahead=2, i.e. the entire system is generated before it is solved.
    throw init_error(
        "System::generate_solve() with sectors ordering 1 "
        "(i.e. by sector number) can only be used with lookahead=-1 "
        "(no sorting across seed sectors) or lookahead=2 "
        "(generate the entire system before sorting).");
  }
  m_is_generated = true;
  m_seeds = std::make_tuple(ibp_seedspec, ibp_seedcompl, sym_seedspec);
  std::vector<std::size_t> indep_eqnums;
#define PYRED_PP_CCS_SYSTEMGENERATESOLVE(k)                     \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) {           \
    indep_eqnums = generate_solve_tmpl<PYRED_PP_COEFFCLASS(k)>( \
        ibp_seedspec, ibp_seedcompl, sym_seedspec);             \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMGENERATESOLVE)
  if (Config::backward()) {
    backward();
  }
  return indep_eqnums;
}

void System::generate_retrieve(
    const std::vector<SeedSpec>& ibp_seedspec,
    const std::vector<SeedSpec>& ibp_seedcompl,
    const std::vector<SeedSpec>& sym_seedspec,
    std::vector<std::size_t>&& eqnums,
    const std::function<void(eqdata&&)>& treateq,
    const std::function<std::string(const std::string&)>& treatcoeff) {
  // Generate a system of equation for the given seed, insert coefficients
  // as strings, optionally simplified by treatcoeff,
  // and apply treateq to all equations (i.e. no selection).
#define PYRED_PP_CCS_SYSTEMGENERATERETRIEVE1(k)                            \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) {                      \
    GeneratorHelper::generate_and_retrieve<PYRED_PP_COEFFCLASS(k)>(        \
        ibp_seedspec, ibp_seedcompl, sym_seedspec, Config::auto_symseed(), \
        std::move(eqnums), treatcoeff, treateq, Config::parallel());       \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMGENERATERETRIEVE1)
}

void System::generate_retrieve(
    std::vector<std::size_t>&& eqnums, const std::function<void(eqdata&&)>& treateq,
    const std::function<std::string(const std::string&)>& treatcoeff) {
  // Re-generate the system, insert coefficients as strings,
  // optionally simplified by treatcoeff,
  // and apply treateq to the (selected) equations with numbers in eqnums.
  if (!m_is_generated) {
    throw init_error("System::generate_retrieve(): seeds have not been set for "
                     "re-generation.");
  }
#define PYRED_PP_CCS_SYSTEMGENERATERETRIEVE2(k)                           \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) {                     \
    GeneratorHelper::generate_and_retrieve<PYRED_PP_COEFFCLASS(k)>(       \
        std::get<0>(m_seeds), std::get<1>(m_seeds), std::get<2>(m_seeds), \
        Config::auto_symseed(), std::move(eqnums), treatcoeff, treateq,   \
        Config::parallel());                                              \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMGENERATERETRIEVE2)
}

void System::reserve(const std::size_t sz) { sys.reserve(sz); }

void System::add(const eqdata& eq) { sys.push_back(eq); }

void System::add(eqdata&& eq) { sys.push_back(std::move(eq)); }

void System::add(const std::vector<eqdata>& eqs) {
  reserve(sys.size() + eqs.size());
  for (const auto& eq : eqs) {
    sys.push_back(eq);
  }
}

void System::add(std::vector<eqdata>&& eqs) {
  reserve(sys.size() + eqs.size());
  for (auto& eq : eqs) {
    sys.push_back(std::move(eq));
  }
}

void System::add_forward(const eqdata& eq, std::size_t neq) {
#define PYRED_PP_CCS_SYSTEMADDFORWARD(k)              \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) { \
    add_forward_tmpl<PYRED_PP_COEFFCLASS(k)>(         \
        eq, get_numsys<PYRED_PP_COEFFCLASS(k)>(),     \
        m_sols.PYRED_PP_COEFFCLASSMEM(k), neq);       \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMADDFORWARD)
}

void System::add_forward(const std::string &infile,
                         const std::vector<std::string> &fileexts,
                         bool unsafe) {
  add_forward(std::vector<std::string>{infile}, fileexts, unsafe);
}

void System::add_forward(const std::vector<std::string> &infiles,
                         const std::vector<std::string> &fileexts,
                         bool unsafe) {
  std::vector<std::string> infilesext;
  for (const std::string &infile: infiles) {
    for (const auto &fname: FileSystem::get_filenames(infile, fileexts)) {
      infilesext.push_back(fname);
    }
  }
  if (!sys.capacity()) {
    reserve(System::file_retrieve(infilesext, fileexts, {}, [](eqdata &&){}));
  }
  std::string line;
  eqdata tmpeq;
  for (const auto &fname : infilesext) {
    // This also works with uncompressed files.
    GZSTREAM_NAMESPACE::igzstream instream{fname.c_str()};
    if (!instream.good()) {
      std::ostringstream ss;
      ss << "System::add_forward(): failed reading file \"" << fname << "\""
         << std::endl;
      throw input_error(ss.str());
    }
    while (pyred::getline_normalized(instream, line)) {
      if (line.find_first_not_of(' ') == std::string::npos) {
        if (!tmpeq.empty()) {
          if (!unsafe) sorteq(tmpeq);
          add_forward(tmpeq);
          tmpeq.clear();
        }
      }
      else {
        auto wc = split_intcoeff(line);
        if (wc.first != Weight::none()) {
          tmpeq.push_back(std::move(wc));
        }
      }
    }
    if (!tmpeq.empty()) {
      if (!unsafe) sorteq(tmpeq);
      add_forward(tmpeq);
      tmpeq.clear();
    }
    instream.close();
  }
}

void System::add(const std::string& infile,
                 const std::vector<std::string> &fileexts,
                 bool unsafe) {
  Config::log(1) << "import equations" << std::flush;
  auto starttime = std::clock();
  for (const auto &fname: FileSystem::get_filenames(infile, fileexts)) {
    add_file(fname, {}, unsafe);
  }
  auto endtime = std::clock();
  Config::log(1) << ":   " << time_diff(starttime, endtime) << std::endl;
}

void System::add(const std::vector<std::string>& infiles,
                 const std::vector<std::string> &fileexts,
                 bool unsafe) {
  Config::log(1) << "import equations" << std::flush;
  auto starttime = std::clock();
  for (const std::string &infile: infiles) {
    for (const auto &fname: FileSystem::get_filenames(infile, fileexts)) {
      add_file(fname, {}, unsafe);
    }
  }
  auto endtime = std::clock();
  Config::log(1) << ":   " << time_diff(starttime, endtime) << std::endl;
}

void System::add_file(const std::string& infile,
                      const std::vector<std::string> &fileexts,
                      bool unsafe) {
  std::string line;
  eqdata tmpeq;
  if (infile.empty()) {
    // read equations from stdin
    while (pyred::getline_normalized(std::cin, line)) {
      if (line.find_first_not_of(' ') == std::string::npos && !tmpeq.empty()) {
        if (!unsafe) sorteq(tmpeq);
        sys.push_back(std::move(tmpeq));
        tmpeq.clear();
      }
      else {
        auto wc = split_intcoeff(line);
        if (wc.first != Weight::none()) {
          tmpeq.push_back(std::move(wc));
        }
      }
    }
    if (!tmpeq.empty()) {
      if (!unsafe) sorteq(tmpeq);
      sys.push_back(std::move(tmpeq));
    }
  }
  else {
    for (const auto &fname : FileSystem::get_filenames(infile, fileexts)) {
      // This also works with uncompressed files.
      GZSTREAM_NAMESPACE::igzstream instream{fname.c_str()};
      if (!instream.good()) {
        std::ostringstream ss;
        ss << "failed reading file \"" << fname << "\"" << std::endl;
        throw input_error(ss.str());
      }
      while (pyred::getline_normalized(instream, line)) {
        if (line.find_first_not_of(' ') == std::string::npos) {
          if (!tmpeq.empty()) {
            if (!unsafe) sorteq(tmpeq);
            sys.push_back(std::move(tmpeq));
            tmpeq.clear();
          }
        }
        else {
          auto wc = split_intcoeff(line);
          if (wc.first != Weight::none()) {
            tmpeq.push_back(std::move(wc));
          }
        }
      }
      if (!tmpeq.empty()) {
        if (!unsafe) sorteq(tmpeq);
        sys.push_back(std::move(tmpeq));
        tmpeq.clear();
      }
      instream.close();
    }
  }
}

std::vector<std::size_t> System::solve() {
  std::vector<std::size_t> indep_eqnums;
#define PYRED_PP_CCS_SYSTEMSOLVE(k)                      \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) {    \
    indep_eqnums = solve_tmpl<PYRED_PP_COEFFCLASS(k)>(); \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMSOLVE)
  return indep_eqnums;
}

void System::backward() {
#define PYRED_PP_CCS_SYSTEMBACKWARD(k)                \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) { \
    m_sols.PYRED_PP_COEFFCLASSMEM(k).clear();         \
    backward_tmpl<PYRED_PP_COEFFCLASS(k)>(            \
        get_numsys<PYRED_PP_COEFFCLASS(k)>());        \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMBACKWARD)
}


void System::backward_coeff() {
#define PYRED_PP_CCS_SYSTEMBACKWARD_COEFF(k)                \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) { \
    m_sols.PYRED_PP_COEFFCLASSMEM(k).clear();         \
    backward_coeff_tmpl<PYRED_PP_COEFFCLASS(k)>(            \
        get_numsys<PYRED_PP_COEFFCLASS(k)>());        \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMBACKWARD_COEFF)
}

std::vector<std::size_t> System::independent() {
  std::vector<std::size_t> indep_eqnums;
  std::size_t nequations{0};
#define PYRED_PP_CCS_SYSTEMINDEPENDENT(k)                \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) {    \
    auto& numsys = get_numsys<PYRED_PP_COEFFCLASS(k)>(); \
    nequations = numsys.size();                          \
    indep_eqnums = numsys.independent();                 \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMINDEPENDENT)
  auto neqs_indep = indep_eqnums.size();
  Config::log(1) << std::setw(8) << nequations
                 << " equations: " << nequations - neqs_indep << " zero + "
                 << neqs_indep << " independent" << std::endl;
  return indep_eqnums;
}

std::vector<Weight> System::unreduced() const {
  if (!m_content_prepared) {
    throw init_error("System::unreduced(): Content map has not been prepared. "
                     "The system must be solved with backward insertion.");
  }
  std::unordered_set<Weight> unred_set;
  std::vector<Weight> unred_vec;
  for (const auto &eq_cont: m_content) {
    for (const auto &unred_int: eq_cont.second.second) {
      unred_set.insert(unred_int);
    }
  }
  for (const auto &unred_int: unred_set) {
    unred_vec.push_back(unred_int);
  }
  std::sort(unred_vec.begin(), unred_vec.end());
  return unred_vec;
}

std::unique_ptr<keyvaluedb::KeyValueDB>& System::get_db() {
#define PYRED_PP_CCS_SYSTEMGETDB(k)                       \
  PYRED_PP_IF_TRUE(k, else) if (coeff_cls == k + 1) {     \
    return get_numsys<PYRED_PP_COEFFCLASS(k)>().get_db(); \
  }
  PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_CCS_SYSTEMGETDB)
  else {
    // Avoid compiler warning
    throw init_error("Oops! This should never happen.");
  }
}

std::pair<std::vector<std::size_t>, std::vector<Weight>> System::select(
    const std::vector<Weight>& mandatory, const std::vector<Weight>& if_reduced) {
  /*
   * 'mandatory':
   * - Vector of integrals for which a reduction formula is required.
   * 'if_reduced':
   * - Vector of integrals to reduce, if there is a reduction formula known.
   * return: pair(neededeqs, unreduced)
   * 'neededeqs'
   * - The numbers of all equations which are needed to reduce integrals
   *   from 'mandatory' and 'if_reduced'.
   * 'unreduced':
   * - All integrals which appear in reduction formulas of mandatory integrals.
   * - All mandatory integrals for which there is no reduction formula known.
   */
  if (!m_content_prepared) {
    throw init_error("Content map has not been prepared for selection. "
                     "The system must be solved with backward insertion.");
  }
  // Stack of untreated integrals for which the dependencies must be looked up.
  std::stack<Weight, std::vector<Weight>> intstack;
  std::unordered_set<Weight> neededset;
  std::unordered_set<Weight> mandatoryset;
  // Set of integrals which are regarded as master integrals.
  std::unordered_set<Weight> unreducedset;
  for (const auto &iid: mandatory) {
    if (iid != Weight::none()) neededset.insert(iid);
  }
  mandatoryset = neededset;
  for (const auto &iid: if_reduced) {
    if (iid != Weight::none()) neededset.insert(iid);
  }
  for (const auto &iid: neededset)
    intstack.push(iid);
  while (!intstack.empty()) {
    auto next = intstack.top();
    intstack.pop();
    const auto cntit = m_content.find(next);
    if (cntit != m_content.cend()) {
      // A reduction formula is known.
      // Depends on equations with the numbers in the respective 'insertions'.
      // Get them from the insertions database.
      auto deps = get_db()->get(cntit->second.first /*eqnum*/, false);
      for (const auto& dep : deps) {
        // For each integral on which 'next' depends:
        // add it to the set of needed integrals.
        // If it was not already in the set, push it onto the stack.
        auto isnew = neededset.insert(dep).second;
        if (isnew) intstack.push(dep);
      }
      if (mandatoryset.find(next) != mandatoryset.cend()) {
        // If 'next' is in 'mandatory', treat integral content masters.
        for (const auto& iid : cntit->second.second /*unreduced*/) {
          unreducedset.insert(iid);
        }
      }
    }
    else {
      // No reduction formula is known.
      // Dependencies always correspond to (inserted) equations,
      // i.e. this branch is only reached for unreduced integrals
      // in 'mandatory' and 'if_reduced'.
      if (mandatoryset.find(next) != mandatoryset.cend()) {
        // 'next' is in 'mandatory', but there is no reduction formula.
        // -> Regard it as a master integral.
        unreducedset.insert(next);
      }
      else {
        // 'next' is in 'if_reduced', but there is no reduction formula.
        // -> Remove it from neededset. It won't be added again,
        //    because it was never inserted anywhere.
        // Do not regard it as a master integral.
        neededset.erase(next);
      }
    }
  }
  // For each integral id in 'neededset':
  // push the corresponding equation number on 'neededeqs',
  // if a reduction formula exists.
  std::vector<std::size_t> neededeqs;
  neededeqs.reserve(neededset.size());
  for (const auto &iid: neededset) {
    const auto cntit = m_content.find(iid);
    if (cntit != m_content.cend()) {
      neededeqs.push_back(cntit->second.first /*eqnum*/);
    }
  }
  std::sort(neededeqs.begin(), neededeqs.end());
  // Convert 'unreducedset' to a vector.
  std::vector<Weight> unreduced;
  for (const auto &iid: unreducedset)
    unreduced.push_back(iid);
  std::sort(unreduced.begin(), unreduced.end());
  return {neededeqs, unreduced};
}

std::pair<std::vector<size_t>, std::vector<Weight>> System::select_coeff(
    const std::vector<Weight>& mandatory, const std::vector<Weight>& if_reduced) {
  /*
   * 'mandatory':
   * - Vector of integrals for which a reduction formula is required.
   * 'if_reduced':
   * - Vector of integrals to reduce, if there is a reduction formula known.
   * return: pair(neededeqs, unreduced)
   * 'neededeqs'
   * - The numbers of all equations which are needed to reduce integrals
   *   from 'mandatory' and 'if_reduced'.
   * 'unreduced':
   * - All integrals which appear in reduction formulas of mandatory integrals.
   * - All mandatory integrals for which there is no reduction formula known.
   */
  if (!m_content_prepared) {
    throw init_error("Content map has not been prepared for selection. "
                     "The system must be solved with backward insertion.");
  }
  // Stack of untreated integrals for which the dependencies must be looked up.
  std::stack<Weight, std::vector<Weight>> intstack;
  std::unordered_set<Weight> neededset;
  std::unordered_set<Weight> mandatoryset;
  // Set of integrals which are regarded as master integrals.
  std::unordered_set<Weight> unreducedset;
  for (const auto iid : mandatory) {
    if (iid != Weight::none()) neededset.insert(iid);
  }
  mandatoryset = neededset;
  for (const auto iid : if_reduced) {
    if (iid != Weight::none()) neededset.insert(iid);
  }
  for (const auto iid : neededset)
    intstack.push(iid);
  while (!intstack.empty()) {
    auto next = intstack.top();
    intstack.pop();
    const auto cntit = m_content_coeff.find(next);
    if (cntit != m_content_coeff.cend()) {
      // A reduction formula is known.
      // Depends on equations with the numbers in the respective 'insertions'.
      // Get them from the insertions database.
      auto deps = get_db()->get(cntit->second.first /*eqnum*/, false);
      for (const auto& dep : deps) {
        // For each integral on which 'next' depends:
        // add it to the set of needed integrals.
        // If it was not already in the set, push it onto the stack.
        auto isnew = neededset.insert(dep).second;
        if (isnew) intstack.push(dep);
      }
      if (mandatoryset.find(next) != mandatoryset.cend()) {
        // If 'next' is in 'mandatory', treat integral content masters.
        for (const auto& iid : cntit->second.second /*unreduced*/) {
          unreducedset.insert(std::get<0>(iid));
        }
      }
    }
    else {
      // No reduction formula is known.
      // Dependencies always correspond to (inserted) equations,
      // i.e. this branch is only reached for unreduced integrals
      // in 'mandatory' and 'if_reduced'.
      if (mandatoryset.find(next) != mandatoryset.cend()) {
        // 'next' is in 'mandatory', but there is no reduction formula.
        // -> Regard it as a master integral.
        unreducedset.insert(next);
      }
      else {
        // 'next' is in 'if_reduced', but there is no reduction formula.
        // -> Remove it from neededset. It won't be added again,
        //    because it was never inserted anywhere.
        // Do not regard it as a master integral.
        neededset.erase(next);
      }
    }
  }
  // For each integral id in 'neededset':
  // push the corresponding equation number on 'neededeqs',
  // if a reduction formula exists.
  std::vector<size_t> neededeqs;
  neededeqs.reserve(neededset.size());
  for (const auto iid : neededset) {
    const auto cntit = m_content_coeff.find(iid);
    if (cntit != m_content_coeff.cend()) {
      neededeqs.push_back(cntit->second.first /*eqnum*/);
    }
  }
  std::sort(neededeqs.begin(), neededeqs.end());
  // Convert 'unreducedset' to a vector.
  std::vector<Weight> unreduced;
  for (const auto iid : unreducedset)
    unreduced.push_back(iid);
  std::sort(unreduced.begin(), unreduced.end());
  return {neededeqs, unreduced};
}

std::vector<eqdata> System::retrieve(std::vector<size_t>&& eqnums) {
  if (sys.empty()) {
    throw init_error("System::retrieve(): empty system.");
  }
  if (eqnums.empty()) {
    return {};
  }
  std::sort(eqnums.begin(), eqnums.end());
  std::vector<eqdata> result;
  std::size_t neq = 0u;
  auto selit = eqnums.begin();
  auto selend = eqnums.end();
  for (auto& eq : sys) {
    if (neq == *selit) {
      result.push_back(std::move(eq));
      if (++selit == selend) break;
    }
    eq.clear();
    ++neq;
  }
  sys.clear();
  if (m_solver_sorted) {
    std::stable_sort(result.begin(), result.end(), cmp_eqdata<std::string>);
  }
  return result;
}

std::size_t System::file_retrieve(
    const std::string &infile,
    const std::vector<std::string> &fileexts,
    std::vector<std::size_t> &&eqnums,
    const std::function<void(eqdata &&)> &treateq,
    const std::function<std::string(std::string &)> &treatcoeff,
    bool unsafe) {
  return file_retrieve(std::vector<std::string>{infile}, fileexts,
                       std::move(eqnums), treateq, treatcoeff, unsafe);
}

std::size_t System::file_retrieve(
    const std::vector<std::string> &infiles,
    const std::vector<std::string> &fileexts,
    std::vector<std::size_t> &&eqnums,
    const std::function<void(eqdata &&)> &treateq,
    const std::function<std::string(std::string &)> &treatcoeff,
    bool unsafe) {

  class TreatEqCoeff {
    std::function<void(eqdata &&)> m_treateq;
    std::function<std::string(std::string &)> m_treatcoeff;
#define PYRED_PP_FILE_RETRIEVE_CACHE(k) \
    std::unordered_map<PYRED_PP_COEFFCLASS(k),std::string> m_coeff_cache_##k;
    PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_FILE_RETRIEVE_CACHE)
  public:
    TreatEqCoeff(const std::function<void(eqdata &&)> &treateq,
                 const std::function<std::string(std::string &)> &treatcoeff)
    : m_treateq{treateq}
    , m_treatcoeff{treatcoeff}
    {}
    void operator()(eqdata &&eq) {
      if (m_treatcoeff) {
#define PYRED_PP_FILE_RETRIEVE_TREAT(k)                                      \
        PYRED_PP_IF_TRUE(k, else) if (Config::coeff_cls() == k + 1) {        \
          for (auto &ifg_strcf: eq) {                                        \
            auto cf = parse_coeff<PYRED_PP_COEFFCLASS(k)>(ifg_strcf.second); \
            auto it_new = m_coeff_cache_##k.insert({cf, ifg_strcf.second});  \
            if (it_new.second) {                                             \
              ifg_strcf.second = m_treatcoeff(ifg_strcf.second);             \
              it_new.first->second = ifg_strcf.second;                       \
              ifg_strcf.second = ifg_strcf.second;                           \
            }                                                                \
            else {                                                           \
              ifg_strcf.second = it_new.first->second;                       \
            }                                                                \
          }                                                                  \
        }
        PYRED_PP_REPEAT0(PYRED_PP_NCOEFFCLASSES, PYRED_PP_FILE_RETRIEVE_TREAT)
      }
      m_treateq(std::move(eq));
    }
  };
  auto treateq_obj = TreatEqCoeff{treateq, treatcoeff};

  // Read equation files any apply treateq to equations with numbers in eqnums.
  // If eqnums is empty, count equations and return the number.
  std::vector<std::string> infilesext;
  for (const std::string &infile: infiles) {
    for (const auto &fname: FileSystem::get_filenames(infile, fileexts)) {
      infilesext.push_back(fname);
    }
  }
  std::sort(eqnums.begin(), eqnums.end());
  std::string line;
  std::size_t neq{0u};
  eqdata tmpeq;
  bool skip_this_eq = false;
  const bool count_only = eqnums.empty();
  const bool select_all = (eqnums.size() == 1u &&
                           eqnums[0] == std::numeric_limits<std::size_t>::max());
  auto eqnum_it = eqnums.begin();
  for (const auto &fname : infilesext) {
    if (!count_only && eqnum_it == eqnums.end()) {
      break;
    }
    // This also works with uncompressed files.
    GZSTREAM_NAMESPACE::igzstream instream{fname.c_str()};
    if (!instream.good()) {
      std::ostringstream ss;
      ss << "System::file_retrieve(): failed reading file \"" << fname << "\""
         << std::endl;
      throw input_error(ss.str());
    }
    while (pyred::getline_normalized(instream, line)) {
      if (line.find_first_not_of(' ') == std::string::npos) {
        // line is empty or contains only whitespace
        if (!tmpeq.empty()) {
          // Note that if count_only, we never get here, because tmpeq is empty.
          if (neq == *eqnum_it || select_all) {
            if (!unsafe) sorteq(tmpeq);
            treateq_obj(std::move(tmpeq));
            if (!select_all) ++eqnum_it;
            if (eqnum_it == eqnums.end()) {
              tmpeq.clear();
              break;
            }
          }
          tmpeq.clear();
          skip_this_eq = false;
          ++neq;
        }
        else if (skip_this_eq) {
          skip_this_eq = false;
          ++neq;
        }
      }
      else {
        if (skip_this_eq) continue;
        auto wc = split_intcoeff(line);
        if (wc.first != Weight::none()) {
          if ((!count_only && neq == *eqnum_it) || select_all) {
            tmpeq.push_back(std::move(wc));
          }
          else {
            skip_this_eq = true;
          }
        }
      }
    }
    if (!tmpeq.empty()) {
      if (neq == *eqnum_it || select_all) {
        if (!unsafe) sorteq(tmpeq);
        treateq_obj(std::move(tmpeq));
        if (!select_all) ++eqnum_it;
      }
      tmpeq.clear();
      skip_this_eq = false;
      ++neq;
    }
    else if (skip_this_eq) {
      skip_this_eq = false;
      ++neq;
    }
    instream.close();
  }
  if (!select_all && eqnum_it != eqnums.end()) {
    throw input_error("System::file_retrieve(): equation number does not "
                      "exist in equation files.");
  }
  return neq;
}

std::size_t System::file_retrieve_sorted(
    const std::string &infile,
    const std::vector<std::string> &fileexts,
    std::vector<std::size_t> &&eqnums,
    const std::function<void(eqdata &&)> &treateq,
    const std::function<std::string(std::string &)> &treatcoeff,
    bool unsafe)
{
  return file_retrieve_sorted(std::vector<std::string>{infile}, fileexts,
                              std::move(eqnums), treateq, treatcoeff, unsafe);
}

std::size_t System::file_retrieve_sorted(
    const std::vector<std::string> &infiles,
    const std::vector<std::string> &fileexts,
    std::vector<std::size_t> &&eqnums,
    const std::function<void(eqdata &&)> &treateq,
    const std::function<std::string(std::string &)> &treatcoeff,
    bool unsafe)
{
  auto eqgetter = System::EquationGetter();
  auto neq = System::file_retrieve(infiles, fileexts, std::move(eqnums),
                                   treateq, treatcoeff, unsafe);
  auto eqs = std::move(eqgetter.get());
  std::stable_sort(eqs.begin(), eqs.end(), cmp_eqdata<std::string>);
  for (auto &eq: eqs) {
    treateq(std::move(eq));
  }
  return neq;
}

// static
std::pair<std::size_t,std::vector<uint32_t>>
System::file_counteqs_setweightbits(const std::vector<std::string> &infiles,
                                    const std::vector<std::string> &fileexts,
                                    bool forcecount) {
  // forcecount: count equations in case of user-defined weights,
  //             i.e. if no weight bits are needed.
  // Otherwise neq=std::numeric_limits<std::size_t>::max() is returned.
  // Parse all files to count equations and determine weight bits.
  std::vector<std::string> infilesext;
  for (const std::string &infile: infiles) {
    for (const auto &fname: FileSystem::get_filenames(infile, fileexts)) {
      infilesext.push_back(fname);
    }
  }
  std::string line;
  std::size_t neq = 0u;
  std::vector<uint32_t> weight_bits;
  auto maxdots4lines = std::vector<uint32_t>(Integral::np(), 0u);
  auto maxsps4nums = std::vector<uint32_t>(Integral::np(), 0u);
  // For user-defined systems without config,
  // np must first be determined from an integral in the file.
  bool np_unset = !Integral::np();
  bool within_equation = false;
  bool userweights = false;
  for (const auto &fname : infilesext) {
    // This also works with uncompressed files.
    GZSTREAM_NAMESPACE::igzstream instream{fname.c_str()};
    if (!instream.good()) {
      std::ostringstream ss;
      ss << "System::file_neqs_weights(): failed reading file \""
         << fname << "\"" << std::endl;
      throw input_error(ss.str());
    }
    while (pyred::getline_normalized(instream, line)) {
      if (line.find_first_not_of(' ') == std::string::npos) {
        // line is empty or contains only whitespace
        if (within_equation) {
          ++neq;
          within_equation = false;
        }
      }
      else {
        within_equation = true;
        if (userweights) {
          continue;
        }
        std::istringstream ss{line};
        auto igl = Integral(ss);
        if (!Integral::np()) {
          // Special case of user-defined weights: no need for weight bits.
          // --> Skip parsing integrals from here and only count equations.
          userweights = true;
          weight_bits = {0u,0u,0u,0u};
          if (!forcecount) {
            neq = std::numeric_limits<std::size_t>::max();
            return {neq, weight_bits};
          }
          continue;
        }
        if (np_unset) {
          maxdots4lines.resize(Integral::np(), 0u);
          maxsps4nums.resize(Integral::np(), 0u);
          np_unset = false;
        }
        uint32_t lines = 0u;
        uint32_t dots = 0u;
        uint32_t nums = 0u;
        uint32_t sps = 0u;
        for (const auto pow: igl.m_powers) {
          if (pow > 0) {
            ++lines;
            dots += static_cast<uint32_t>(pow - 1);
          }
          else {
            ++nums;
            sps += static_cast<uint32_t>(-pow);
          }
        }
        if (dots) {
          maxdots4lines[lines-1] = std::max(maxdots4lines[lines-1], dots);
        }
        if (sps) {
          maxsps4nums[nums-1] = std::max(maxsps4nums[nums-1], sps);
        }
      }
    }
    if (within_equation) {
      ++neq;
    }
    instream.close();
  }
  uint32_t maxdots = 0u;
  uint32_t maxsps = 0u;
  uint32_t maxdotcombs = 1u;
  uint32_t maxspcombs = 1u;
  for (uint32_t lines = 0u; lines <= Integral::np(); ++lines) {
    if (lines) {
      auto dotsonlines = maxdots4lines[lines-1];
      maxdotcombs = std::max(static_cast<std::size_t>(maxdotcombs),
                             binomial(lines + dotsonlines - 1, dotsonlines));
      maxdots = std::max(maxdots, dotsonlines);
    }
    auto nums = Integral::np() - lines;
    if (nums) {
      auto spsonnums = maxsps4nums[nums-1];
      maxspcombs = std::max(static_cast<std::size_t>(maxspcombs),
                            binomial(nums + spsonnums - 1, spsonnums));
      maxsps = std::max(maxsps, spsonnums);
    }
  }
  if (!userweights) {
    weight_bits = Integral::assign_weight_bits(
      maxdots, maxsps, maxdotcombs, maxspcombs);
  }
  return {neq, weight_bits};
}

std::vector<uint32_t>
System::reserve_setweightbits(const std::vector<std::string> &infiles,
                              const std::vector<std::string> &fileexts,
                              bool forcecount) {
  auto neqs_weightbits = file_counteqs_setweightbits(
    infiles, fileexts, forcecount);
  if (neqs_weightbits.first != std::numeric_limits<std::size_t>::max()) {
    sys.reserve(neqs_weightbits.first);
  }
  return neqs_weightbits.second;
}

} // namespace pyred
