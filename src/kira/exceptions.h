/*
Copyright (C) 2017-2025 The Kira Developers (see AUTHORS file)

This file is part of Kira.

Kira is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Kira is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Kira.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <exception>
#include <string>

class ExceptionCommandLine : public std::exception {
  std::string full_msg;
public:
  ExceptionCommandLine(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionCommandLine: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

class ExceptionConfigFiles : public std::exception {
  std::string full_msg;
public:
  ExceptionConfigFiles(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionConfigFiles: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

class ExceptionDisk : public std::exception {
  std::string full_msg;
public:
  ExceptionDisk(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionDisk: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

class ExceptionFermat : public std::exception {
  std::string full_msg;
public:
  ExceptionFermat(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionFermat: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

class ExceptionInternal : public std::exception {
  std::string full_msg;
public:
  ExceptionInternal(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionInternal: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

class ExceptionPipe : public std::exception {
  std::string full_msg;
public:
  ExceptionPipe(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionPipe: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

class ExceptionResume : public std::exception {
  std::string full_msg;
public:
  ExceptionResume(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionResume: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

class ExceptionRuntime : public std::exception {
  std::string full_msg;
public:
  ExceptionRuntime(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionRuntime: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

// not a real exception, just quit the program and return 0
class ExceptionQuit : public std::exception {
  std::string full_msg;
public:
  ExceptionQuit(const std::string& msg_ = "Unspecified error") {
    full_msg = "Kira::ExceptionQuit: " + msg_;
  };
  const char* what() const throw () {
    return full_msg.c_str();
  }
};

#endif
