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

#include <map>
#include <utility>

#include "pyred/parser.h"

namespace pyred {

std::vector<std::string> split(const std::string& s, char delim, int n,
                               bool discardempty) {
  // Split the string s at delimiters delim.
  // If 'discardempty' (default=true), discard empty items
  // arising from multiple adjacent delimitors.
  // Apply at most n splits if n>0 (not counting discarded items).
  std::vector<std::string> elems;
  std::istringstream ss{s};
  std::string item;
  int k = 0;
  while (pyred::getline_normalized(ss, item, delim)) {
    if (discardempty && item.empty()) continue;
    elems.push_back(item);
    if (n > 0 && ++k >= n) {
      if (pyred::getline_normalized(ss, item)) {
        if (!item.empty()) {
          elems.push_back(item);
        }
      }
      break;
    }
  }
  return elems;
}

// (char)operator -> {(int)priority, (bool)is_left_associate}
const std::map<char, std::pair<int, bool>> operator_props = {
    {'+', {0, true}}, {'-', {0, true}}, {'*', {1, true}},
    {'/', {1, true}}, {'~', {2, true}}, // unary '-'
    {'!', {4, true}}, // high priority unary '-' to be used after '^'
    {'^', {3, false}}};

std::vector<char> get_operators(const std::map<char, std::pair<int, bool>>& m) {
  // Extract operators from operator_props map.
  std::vector<char> keys;
  keys.reserve(m.size());
  for (const auto& p : m) {
    keys.push_back(p.first);
  }
  return keys;
}
const std::vector<char> operators = get_operators(operator_props);

bool is_left_assoc(const char op) {
  return operator_props.find(op)->second.second;
}

int op_priority_diff(const char op1, const char op2) {
  return (operator_props.find(op1)->second.first -
          operator_props.find(op2)->second.first);
}

bool is_operator(const char ch) {
  for (const auto op : operators) {
    if (ch == op) {
      return true;
    }
  }
  return false;
}

std::vector<std::string> tokenize(const std::string& expr) {
  // Split string expr into tokens at boundaries between
  // numbers/symbols and (char) operators.
  // Ignore spaces unless they separate two arguments, which causes an error.
  // Check for invalid structures.
  // Note that the high priority unary minus '!' is not set here:
  // In parse_coeff(), '~' after '^' is replaced by '!'.
  std::vector<std::string> tokens;
  std::string arg = "";
  bool last_was_arg = false;
  for (const char ch : expr) {
    if (is_operator(ch) || ch == '(' || ch == ')' || ch == ' ') {
      if (!arg.empty()) {
        if (last_was_arg) {
          throw parser_error("two successive arguments without operator");
        }
        else if (!tokens.empty() && tokens.back().front() == ')') {
          throw parser_error("argument following right parenthesis");
        }
        tokens.push_back(arg);
        arg.clear();
        last_was_arg = true;
      }
      if (ch != ' ') {
        if (tokens.empty()) {
          // leading operator: ignore if '+', unary if '-', replace by '~',
          // error if other operator or right parenthesis
          if (ch == '-') {
            tokens.push_back(std::string("~"));
          }
          else if (ch == '*' || ch == '/' || ch == '^' || ch == ')') {
            throw parser_error("invalid leading operator or parenthesis");
          }
          else if (ch != '+') {
            tokens.push_back(std::string(1, ch));
          }
        }
        else if (!last_was_arg) {
          // if the previous token was also an operator/parenthesis
          auto pt = tokens.back().front();
          if (ch == '+' &&
              (pt == '*' || pt == '/' || pt == '^' || pt == '(' || pt == '~')) {
            // ignore unary '+': *+ /+ ^+ (+ ~+
          }
          else if (ch == '-' &&
                   (pt == '*' || pt == '/' || pt == '^' || pt == '(')) {
            // unary '-', replace by '~': *- /- ^- (-
            tokens.push_back(std::string("~"));
          }
          else if ((ch == '-' && pt == '~') ||
                   ((ch == '+' || ch == '-') && (pt == '+' || pt == '-'))) {
            if (ch == '-' && pt == '~') {
              // unary '-', following '~' i.e. another unary '-':
              // neutralise the '~', i.e. unary '--' == '~-' --> ''
              tokens.pop_back();
            }
            else if (pt == ch) {
              // '++', '--' --> '+'
              tokens.back() = std::string("+");
            }
            else {
              // '-+', '+-' --> '-'
              tokens.back() = std::string("-");
            }
          }
          else if ((ch == '(' && pt != ')') || (pt == ')' && ch != '(')) {
            // remaining allowed operator pairs:
            // +( -( *( /( ^( ((
            // )+ )- )* )/ )^ ))
            tokens.push_back(std::string(1, ch));
          }
          else {
            std::ostringstream ss;
            ss << "invalid combination of adjacent operators: \"" << pt << ch
               << "\"";
            throw parser_error(ss.str());
          }
        }
        else {
          if (ch == '(') {
            std::ostringstream ss;
            ss << "left parenthesis following argument \"" << tokens.back()
               << "\"";
            throw parser_error(ss.str());
          }
          tokens.push_back(std::string(1, ch));
        }
        last_was_arg = false;
      }
    }
    else {
      arg.push_back(ch);
    }
  }
  if (!arg.empty()) {
    if (last_was_arg) {
      throw parser_error("two successive arguments without operator");
    }
    tokens.push_back(arg);
  }
  if (tokens.empty()) {
    throw init_error("Empty expression");
  }
  if (is_operator(tokens.back().front())) {
    throw parser_error("Operator which needs an argument to the right");
  }
  return tokens;
}

} // namespace pyred
