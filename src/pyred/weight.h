
#ifndef PYRED_WEIGHT_H
#define PYRED_WEIGHT_H

#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <iostream>
#include <sstream>
#include <exception>


namespace pyred {

using std::uint64_t;

class weight_error : public std::exception {
private:
  std::string msg;
public:
  inline weight_error(const std::string& s) : msg(s) {}
  virtual inline const char* what() const noexcept { return msg.c_str(); }
};

class Weight {
  friend class Integral;
  friend bool operator==(const Weight &, const Weight &);
  friend bool operator!=(const Weight &, const Weight &);
  friend bool operator<(const Weight &, const Weight &);
  friend bool operator>(const Weight &, const Weight &);
  friend bool operator<=(const Weight &, const Weight &);
  friend bool operator>=(const Weight &, const Weight &);
  friend std::ostream & operator<<(std::ostream &, const Weight &);
  friend std::istream & operator>>(std::istream &, Weight &);

private:
#ifdef WEIGHT128
  using weight_t = unsigned __int128;
#else
  using weight_t = uint64_t;
#endif
  weight_t m_weight;

  Weight & shift_set(uint32_t shift, const weight_t &set) {
    (m_weight <<= shift) |= set;
    return *this;
  }

  // NOTE: This template is not activated for 128-bit projectors,
  //       because std::is_unsigned<unsigned __int128>::value is false.
  template <typename UInt, typename = std::enable_if_t<std::is_unsigned<UInt>::value>>
  UInt project(const UInt &proj = std::numeric_limits<UInt>::max()) const {
    return static_cast<UInt>(m_weight & proj);
  }
  template <typename UInt, typename = std::enable_if_t<std::is_unsigned<UInt>::value>>
  UInt project_shift(const UInt &proj, uint32_t shift) {
    auto res = static_cast<UInt>(m_weight & proj);
    m_weight >>= shift;
    return res;
  }

public:

#ifdef WEIGHT128
  // std::numeric_limits may not be specialised for unsigned __int128.
  constexpr static unsigned int width{128};
  constexpr explicit Weight()
  : m_weight{(static_cast<weight_t>(std::numeric_limits<uint64_t>::max()) << 64) | std::numeric_limits<uint64_t>::max()}
  {}
#else
  constexpr static unsigned int width{std::numeric_limits<weight_t>::digits};
  constexpr explicit Weight(): m_weight{std::numeric_limits<weight_t>::max()} {}
#endif

  constexpr Weight(const Weight &) = default;
  constexpr explicit Weight(const weight_t &w): m_weight{w} {}
  explicit Weight(const std::string &s) {
    std::istringstream ss{s};
    auto success = bool(ss >> *this);
    ss >> std::ws;
    if (!ss.eof() || !success) {
      throw weight_error(
        std::string("Weight(const std::string &): invalid input '") + s + "'.");
    }
  }

  constexpr static Weight none() { return Weight(); }
  constexpr static Weight amplitude(const weight_t &w) { return Weight(~w); }
  constexpr Weight & operator=(const Weight &) = default;

  std::string to_string() const {
#ifdef WEIGHT128
    auto w = m_weight;
    char buf[39];
    char *pos = std::end(buf);
    do {
      *--pos = '0' + (w % 10);
      w /= 10;
    } while (w);
    return std::string(pos, std::end(buf));
#else
    return std::to_string(m_weight);
#endif
  }

  std::size_t hash() const {
#ifdef WEIGHT128
    const auto lo = std::hash<uint64_t>{}(static_cast<uint64_t>(
      m_weight & std::numeric_limits<uint64_t>::max()));
    const auto hi = std::hash<uint64_t>{}(static_cast<uint64_t>(
      m_weight >> 64));
    std::size_t h{0u};
    h = pyred::hash_combine(h, hi);
    h = pyred::hash_combine(h, lo);
    return h;
#else
    return std::hash<uint64_t>{}(m_weight);
#endif
  }
};

inline bool operator==(const Weight &a, const Weight &b) {
  return (a.m_weight == b.m_weight);
}
inline bool operator!=(const Weight &a, const Weight &b) {
  return (a.m_weight != b.m_weight);
}
inline bool operator<(const Weight &a, const Weight &b) {
  return (a.m_weight < b.m_weight);
}
inline bool operator>(const Weight &a, const Weight &b) {
  return (a.m_weight > b.m_weight);
}
inline bool operator<=(const Weight &a, const Weight &b) {
  return (a.m_weight <= b.m_weight);
}
inline bool operator>=(const Weight &a, const Weight &b) {
  return (a.m_weight >= b.m_weight);
}

inline std::ostream & operator<<(std::ostream &out, const Weight &w) {
#ifdef WEIGHT128
  out << w.to_string();
#else
  out << w.m_weight;
#endif
  return out;
}

inline std::istream & operator>>(std::istream &in, Weight &w) {
#ifdef WEIGHT128
  w.m_weight = 0u;
  char ch;
  bool found = false;
  in >> std::ws;
  while (in.get(ch)) {
    if (ch < '0' || ch > '9') {
      in.unget();
      break;
    }
    w.m_weight = 10 * w.m_weight + (ch - '0');
    found = true;
  }
  // clear failbit if eof
  if (in.eof() && !in.bad()) in.clear(std::ios_base::eofbit);
  if (!found) in.setstate(std::ios_base::failbit);
#else
  in >> w.m_weight;
#endif
  return in;
}

} // namespace pyred


namespace std {
  template <>
  struct hash<pyred::Weight> {
    size_t operator()(const pyred::Weight &w) const {
      return w.hash();
    }
  };
} // namespace std

#endif
