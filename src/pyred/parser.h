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

#ifndef PYRED_PARSER_H
#define PYRED_PARSER_H

#include <mutex>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pyred/defs.h"

namespace pyred {

template <typename Coeff>
using arg_stack_t = std::stack<Coeff, std::vector<Coeff>>;

bool is_left_assoc(const char);
int op_priority_diff(const char, const char);
bool is_operator(const char);
std::vector<std::string> split(const std::string& s, char delim, int n = 0,
                               bool discardempty = true);
std::vector<std::string> tokenize(const std::string&);

template <typename T>
T parse_int(const std::string& s) {
  T result;
  std::istringstream ss{s};
  auto success = static_cast<bool>(ss >> result);
  if (!(success && ss.rdbuf()->in_avail() == 0)) {
    throw parser_error(std::string("String to integer parsing failed for \"") +
                       s + "\".");
  }
  return result;
}

template <typename Coeff>
void stack_eval(arg_stack_t<Coeff>& arg_stack, const char op) {
  if (arg_stack.empty()) {
    throw parser_error("operator is lacking arguments");
  }
  Coeff rarg = arg_stack.top();
  arg_stack.pop();
  switch (op) {
    case '+':
      arg_stack.top() += rarg;
      break;
    case '-':
      arg_stack.top() -= rarg;
      break;
    case '~': // unary '-'
    case '!': // high priority unary '-'
      arg_stack.emplace(-rarg);
      break;
    case '*':
      arg_stack.top() *= rarg;
      break;
    case '/':
      arg_stack.top() /= rarg;
      break;
    case '^':
      arg_stack.top() = std::move(pow(arg_stack.top(), rarg));
      break;
    default:
      throw parser_error(std::string("invalid operator ") + std::string(1, op));
  }
  return;
}

template <typename Coeff>
Coeff parse_coeff(const std::string& expr = "") {
  // Shunting-yard algorithm with on-the-fly evaluation.
  // Parenthesis are not supported.
  // parse_coeff() or parse_coeff("") clears the cache.
  thread_local static Cache<std::string, Coeff> s_cache;
  if (expr.empty()) {
    s_cache.clear();
  }
  auto coeff_found = s_cache.lookup(expr);
  if (coeff_found.second) {
    return coeff_found.first;
  }
  arg_stack_t<Coeff> arg_stack;
  std::stack<char, std::vector<char>> op_stack;
  std::vector<std::string> tokens;
  try {
    tokens = tokenize(expr);
  }
  catch (const parser_error& e) {
    throw parser_error(e.message() + std::string(" in expression \"") + expr +
                       "\"");
  }
  for (auto token_it = tokens.begin(); token_it != tokens.end(); ++token_it) {
    auto token = *token_it;
    auto curop = token.front();
    if (is_operator(curop)) {
      if (curop == '^' && (token_it + 1)->front() == '~') {
        // '^' followed by unary minus: replace low priority unary minus '~'
        // by high priority unary minus '!' to be avaluated before '^'.
        // Note that there is always a token after '^'.
        *(token_it + 1) = "!";
      }
      if (!op_stack.empty()) {
        auto prevop = op_stack.top();
        // While the token on top of op_stack is an operator prevop and
        // * curop is left-associative
        //   and its precedence is the same as that of prevop
        // * or curop has lower precedence than prevop:
        while (is_operator(prevop) && ((is_left_assoc(curop) &&
                                        op_priority_diff(curop, prevop) == 0) ||
                                       (op_priority_diff(curop, prevop) < 0))) {
          stack_eval<Coeff>(arg_stack, prevop);
          op_stack.pop();
          if (!op_stack.empty()) {
            prevop = op_stack.top();
          }
          else {
            break;
          }
        }
      }
      // push token onto the operator op_stack.
      op_stack.push(curop);
    }
    else if (curop == '(') {
      op_stack.push(curop);
    }
    else if (curop == ')') {
      // push operators from op_stack onto arg_stack (i.e. evaluate)
      // until a left parenthesis is encountered.
      while (op_stack.top() != '(') {
        stack_eval<Coeff>(arg_stack, op_stack.top());
        op_stack.pop();
        if (op_stack.empty()) {
          // op_stack ran empty without encountering
          // a matching left parenthesis
          throw parser_error(std::string("mismatched parenthesis in ") + expr);
        }
      }
      // pop the left parenthesis from op_stack
      if (!op_stack.empty()) op_stack.pop();
    }
    else {
      // token is an argument, i.e. a number or a symbol
      arg_stack.push(token);
    }
  }
  while (!op_stack.empty()) {
    // apply operators from op_stack to arguments
    if (op_stack.top() == '(' || op_stack.top() == ')') {
      throw parser_error(std::string("mismatched parenthesis in ") + expr);
    }
    stack_eval<Coeff>(arg_stack, op_stack.top());
    op_stack.pop();
  }
  // fill the cache
  s_cache.insert(expr, arg_stack.top());
  return arg_stack.top();
}

} // namespace pyred

#endif
