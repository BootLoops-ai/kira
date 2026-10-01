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

#ifndef CONNECT2KIRA_H
#define CONNECT2KIRA_H
#include <signal.h>

#include <string>
#include <vector>

#include "ginac/ginac.h"
#include "kira/integral.h"
#define DELTA_OUT 134217728
// 268435456
// 536870912
// 1073741824
// 1024
#define FROMFERMATBUFSIZE 2048

class Connect2Kira {
public:
  Connect2Kira() : g_to(0), g_from(0){};
  void pipe_kira();
  void close_pipe();
  void setup(std::vector<std::string> &args);

  ~Connect2Kira() {
    for (size_t i = 0; i < argvLength - 1; i++) {
      delete[] argv[i];
    }
    delete[] argv;
  }

protected:
  unsigned int argvLength;
  char **argv;
  int fdin[2], fdout[2];
  size_t length;
  pid_t g_childpid;
  int status;
  FILE *g_to, *g_from;
};

class Perl2Kira : public Connect2Kira {
public:
  Perl2Kira(const std::vector<GiNaC::possymbol> &invar, GiNaC::ex massONE,
            const std::vector<int> &invarDIM, int DIM);
  Perl2Kira(int flagFireFly, int prefactorFlag);
  ~Perl2Kira(){};
  void read_pipe(std::string &output);
  void put_pipe(std::string &input);
};


/* This is originally from gateToFermat.
 * Note that gateToFermat used to be GPLv2 only, but the author agreed to
 * license it under GPLv2 or later so that it could be included in Kira.
 *
 * This is the wrapper to the program FERMAT
 * http://www.bway.net/~lewis
 * This file is a part of the program gateToFermat.
 * Copyright (C) Mikhail Tentyukov <tentukov@physik.uni-bielefeld.de>
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */
class Fermat : public Connect2Kira {
public:
  Fermat() {
    g_fullout = NULL;
    g_stopout = NULL;
    g_baseout = NULL;
  };
  ~Fermat();
  void start_fermat(std::string &fermatPath, char *pvars);
  void set_variable(char *pvars);
  void unset_variable(char *pvars);
  void set_numeric(char *pvars, int numeric);
  void set_numeric2(char *pvars_, std::string &numeric);
  void unset_numeric(char *pvars);
  void fermat_collect(char *buffer);
  int fermat_calc(int optional = 0);
  void close_calc(int mustCleanup);
  /*The output buffer (baseout)*/
  char *g_baseout;

private:
  char *pvars;
  /*The inline function places one char to the output buffer with possible
  expansion of the buffer:*/
  void add_to_out(char ch);
  /*reads the stream 'from' up to the line 'terminator' (only 'thesize' first
    characters are compared):*/
  void read_up(char *terminator, int thesize);
  /*Starts Fermat and makes some initializations:*/
  char g_fbuf[FROMFERMATBUFSIZE];
  /*stopout points to the end of allocated space, fullout points to
  the end of space filled by an actual content:*/
  char *g_fullout, *g_stopout;
  /*The result from Fermat's cancel command:*/
  std::string cancel_ret;
  bool closed = true;
};

#endif // CONNECT2KIRA_H
