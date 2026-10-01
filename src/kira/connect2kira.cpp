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

#ifdef __linux__
#include <signal.h>
#include <sys/prctl.h>
#endif

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>

#include <fstream>
#include <mutex>
#include <sstream>
#include <string>

#include "kira/connect2kira.h"
#include "kira/exceptions.h"
#include "kira/tools.h"

using namespace std;
using namespace GiNaC;

#define READEND 0
#define WRITEND 1

static Loginfo &logger = Loginfo::instance();

std::mutex mutexFermat;

void Connect2Kira::pipe_kira() {
  if (pipe(fdin) == -1) {
    throw ExceptionPipe("Pipe to fdin failed");
  }

  if (pipe(fdout) == -1) {
    throw ExceptionPipe("Pipe to fdout failed");
  }

  switch (g_childpid = fork()) {
    case -1:
      throw ExceptionPipe("Could not fork, decrease 'DELTA_OUT' in 'connect2kira.h' or run Kira with less cores");

    case 0: // child

#ifdef __linux__
      prctl(PR_SET_PDEATHSIG, SIGKILL);
#endif
      // close pipes not required by the child:
      close(fdin[READEND]);
      close(fdout[WRITEND]);
      // connect stdin and stdout with the pipe:
      dup2(fdout[READEND], 0);
      dup2(fdin[WRITEND], 1);
      close(2);
      open("/dev/null", O_WRONLY);

      execvp(argv[0], argv);
      throw ExceptionPipe(std::string(argv[0]) + " returned unexpectedly, must be an execvp error");

    default: // parent
      // close pipes not required by the parent:
      close(fdin[WRITEND]);
      close(fdout[READEND]);
      g_to = fdopen(fdout[WRITEND], "w");
      g_from = fdopen(fdin[READEND], "r");
      if (g_to == NULL || g_from == NULL) {
        if (g_to != NULL) fclose(g_to);
        if (g_from != NULL) fclose(g_from);
        kill(g_childpid, SIGKILL);
        waitpid(g_childpid, &status, 0);
        throw ExceptionFermat("Could not open pipe");
      }
  }
}

void Connect2Kira::setup(vector<string> &argList) {
  argvLength = static_cast<unsigned>(argList.size()) + 1;

  argv = new char *[argvLength];

  for (size_t j = 0; j < argvLength - 1; j++) {
    argv[j] = new char[512];
    strcpy(argv[j], argList[j].c_str());
  }

  argv[argvLength - 1] = NULL;
}

void Connect2Kira::close_pipe() { fclose(g_to); }

Perl2Kira::Perl2Kira(int flagFireFly, int prefactorFlag) {
  vector<string> argList;

  string perlPath = "perl";

  stringstream options;

  char token[2048];

  options << "-ne ";

  if(prefactorFlag == 0){
    sprintf(token, "s/(.*)\\*prefactor(.*)/$1/g; s/(.*)\\/\\((.*)/num($1)*den($2/g;");
    options << token;

    if(flagFireFly == 0){
      sprintf(token, "s/(.*)\\/(?![(])(.*)/num($1)*den($2)/g;");
      options << token;
    }
  }
  else if(prefactorFlag == 1){
    sprintf(token, "s/(.*)\\*prefactor\\[(.*)\\]/$2/g;");
    options << token;
  }

  options << " print;\n";

  argList.push_back(perlPath);
  argList.push_back(options.str());

  setup(argList);

}

Perl2Kira::Perl2Kira(const vector<possymbol> &invarT, ex massOne,
                     const vector<int> &invarDim, int Dim) {
  vector<string> argList;

  string perlPath = "perl";

  stringstream options;

  char token[2048];

  options << "-ne ";


  for (size_t i = 0; i < invarT.size(); i++) {
    string strA, strB;
    strA = something_string(invarT[i]);
    strB = something_string(massOne);
    if (strA == strB) continue;

    sprintf(token, "s/(?<![a-zA-Z0-9])%s(?![a-zA-Z0-9])/(%s\\/%s^(%i\\/%i))/g;",
            strA.c_str(), strA.c_str(), strB.c_str(), invarDim[i], Dim);
    options << token;
  }

  options << " print;\n";

  argList.push_back(perlPath);
  argList.push_back(options.str());

  setup(argList);
}

void Perl2Kira::put_pipe(string &input) {
  input += "\n";
  fputs(input.c_str(), g_to);
  putc('\n', g_to);
  input.pop_back();
}

void Perl2Kira::read_pipe(string &receiveOutput) {
  char readBuffer[90];
  int status;
  receiveOutput.clear();

  while (1) {
    unsigned bytesRead = static_cast<unsigned>(
        fread(readBuffer, 1, sizeof(readBuffer) - 1, g_from));
    if (bytesRead <= 0) break;

    readBuffer[bytesRead] = '\0';
    receiveOutput += readBuffer;
  }
  fclose(g_from);
  wait(&status);
  receiveOutput.erase(receiveOutput.end() - 2, receiveOutput.end());
}


Fermat::~Fermat() { close_calc(1); }

void Fermat::start_fermat(string &fermatPath, char *pvars_) {
  vector<string> argList;
  argList.push_back(fermatPath);
  setup(argList);

  /*Allocate output buffer:*/
  g_baseout = new char[DELTA_OUT];
  if (g_baseout != NULL) {
    g_fullout = g_baseout;
    g_stopout = g_baseout + DELTA_OUT;
  }
  else {
    throw ExceptionFermat("Memory allocation error");
  }

  pvars = pvars_;

  pipe_kira();

  if (g_childpid == 0)
    throw ExceptionFermat("Fermat died unexpectedly");
  /*Fermat is running*/

  /*Switch off floating point representation:*/
  fputs("&d\n0\n", g_to);
  fflush(g_to);
  read_up((char *)"0", 1);

  if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
    throw ExceptionFermat("Failed to switch off floating point representation");
  }

  /*Switch off poly_divide:*/
  fputs("&(_t=0)\n", g_to);
  fflush(g_to);
  read_up((char *)"0", 1);
  if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
    throw ExceptionFermat("Failed to switch off 'poly_divide'");
  }
  /*Set prompt as '\n':*/
  fputs("&(M=' ')\n", g_to);
  fflush(g_to);
  read_up((char *)"> Prompt", 8);
  read_up((char *)"Elapsed", 7);
  if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
    throw ExceptionFermat("Failed to set prompt to \'\\n\'");
  }

  /*Switch off timing:*/
  fputs("&(t=0)\n", g_to);
  fflush(g_to);
  read_up((char *)"0", 1);
  if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
    throw ExceptionFermat("Failed to switch off timing");
  }

  /*Switch on "ugly" printing: no spaces in long integers;
    do not omit '*' for multiplication:*/
  fputs("&U\n", g_to);
  fflush(g_to);
  read_up((char *)"0", 1);
  if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
    throw ExceptionFermat("Failed to switch off ugly printing");
  }

  /*Switch off suppression of output g_to terminal of long polys.:*/
  fputs("&(_s=0)\n", g_to);
  fflush(g_to);
  read_up((char *)"0", 1);
  if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
    throw ExceptionFermat("Failed to switch off suppression of \noutput g_to terminal of long polys");
  }

  /*Change monomial multiply system (Thanks to Fabian Lange and Robert H. Lewis
   * for pointing out this feature):*/
  fputs("&(_o=12)\n", g_to);
  fflush(g_to);
  read_up((char *)"0", 1);
  if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
    throw ExceptionFermat("Failed to change monomial multiply system");
  }

  /*Check whether the cancel command outputs 0 or 1:*/
  fputs("a:=1\n", g_to);
  fflush(g_to);
  read_up((char *)"1", 1);
  fermat_collect((char *)"@a");
  fermat_calc();
  cancel_ret = g_baseout;

  /*Set polymomial variables:*/
  /*pvars looks like "a\nb\nc\n\n":*/
  char *ch;
  while (*pvars > '\n') {
    /*Copy the variable up g_to '\n' into g_fbuf:*/
    /* Modification by A.Smirnov: for compatibility with Mathematica*/
    for (*(ch = g_fbuf) = '\0';
         ((*ch = *pvars) > '\n' && (*ch = *pvars) != ','); ch++, pvars++)
      if (ch - g_fbuf >= FROMFERMATBUFSIZE) {
        throw ExceptionFermat("'FROMFERMATBUFSIZE' to small");
      }
    *ch = '\0';
    if (*pvars != '\0') pvars++;
    /*Set g_fbuf as a polymomial variable:*/

    if (isupper(*g_fbuf)) {
      throw ExceptionFermat("Fermat does not like capital letters in variable names, please rename '" + std::string(g_fbuf) + "'");
    }

    fprintf(g_to, "&(J=%s)\n", g_fbuf);

    fflush(g_to);
    read_up((char *)"0", 1);

    if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
      throw ExceptionFermat("Failed to set variable '" + std::string(g_fbuf) + "'");
    }
  } /*while( *pvars > '\n')*/

  closed = false;
}

void Fermat::unset_variable(char *pvars_) {
  char *var = pvars_;

  /*Set polymomial variables:*/
  /*var looks like "a\nb\nc\n\n":*/
  char *ch;
  while (*var > '\n') {
    /*Copy the variable up g_to '\n' into g_fbuf:*/
    /* Modification by A.Smirnov: for compatibility with Mathematica*/
    for (*(ch = g_fbuf) = '\0'; ((*ch = *var) > '\n' && (*ch = *var) != ',');
         ch++, var++)
      if (ch - g_fbuf >= FROMFERMATBUFSIZE) {
        throw ExceptionFermat("'FROMFERMATBUFSIZE' to small");
      }
    *ch = '\0';
    if (*var != '\0') var++;
    /*Set g_fbuf as a polymomial variable:*/

    if (isupper(*g_fbuf)) {
      throw ExceptionFermat("Fermat does not like capital letters in variable names, please rename '" + std::string(g_fbuf) + "'");
    }

    fprintf(g_to, "&(J=-%s)\n", g_fbuf);
    //     fermat_calc();
    //     fermat_calc();
    fflush(g_to);
    read_up((char *)"0", 1);

    if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
      throw ExceptionFermat("Failed to unset variable '" + std::string(g_fbuf) + "'");
    }
  } /*while( *var > '\n')*/
}

void Fermat::set_variable(char *pvars_) {
  char *var = pvars_;

  /*Set polymomial variables:*/
  /*var looks like "a\nb\nc\n\n":*/
  char *ch;
  while (*var > '\n') {
    /*Copy the variable up g_to '\n' into g_fbuf:*/
    /* Modification by A.Smirnov: for compatibility with Mathematica*/
    for (*(ch = g_fbuf) = '\0'; ((*ch = *var) > '\n' && (*ch = *var) != ',');
         ch++, var++)
      if (ch - g_fbuf >= FROMFERMATBUFSIZE) {
        throw ExceptionFermat("'FROMFERMATBUFSIZE' to small");
      }
    *ch = '\0';
    if (*var != '\0') var++;
    /*Set g_fbuf as a polymomial variable:*/

    if (isupper(*g_fbuf)) {
      throw ExceptionFermat("Fermat does not like capital letters in variable names, please rename '" + std::string(g_fbuf) + "'");
    }

    fprintf(g_to, "&(J=%s)\n", g_fbuf);

    fflush(g_to);
    read_up((char *)"0", 1);

    if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
      throw ExceptionFermat("Failed to set variable '" + std::string(g_fbuf) + "'");
    }
  } /*while( *var > '\n')*/
}

void Fermat::set_numeric(char *pvars_, int numeric) {
  char *var = pvars_;

  /*Set polymomial variables:*/
  /*var looks like "a\nb\nc\n\n":*/
  char *ch;
  while (*var > '\n') {
    /*Copy the variable up g_to '\n' into g_fbuf:*/
    /* Modification by A.Smirnov: for compatibility with Mathematica*/
    for (*(ch = g_fbuf) = '\0'; ((*ch = *var) > '\n' && (*ch = *var) != ',');
         ch++, var++)
      if (ch - g_fbuf >= FROMFERMATBUFSIZE) {
        throw ExceptionFermat("'FROMFERMATBUFSIZE' to small");
      }
    *ch = '\0';
    if (*var != '\0') var++;
    /*Set g_fbuf as a polymomial variable:*/

    if (isupper(*g_fbuf)) {
      throw ExceptionFermat("Fermat does not like capital letters in variable names, please rename '" + std::string(g_fbuf) + "'");
    }

    fprintf(g_to, "%s:=%i\n", g_fbuf, numeric);

    fflush(g_to);
    read_up((char *)const_cast<char *>(int_string(numeric).c_str()),
            int_string(numeric).size());

    if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
      throw ExceptionFermat("Failed to set numerics for variable '" + std::string(g_fbuf) + "'");
    }
  } /*while( *var > '\n')*/
}

void Fermat::set_numeric2(char *pvars_, std::string &numeric) {
  char *var = pvars_;

  /*Set polymomial variables:*/
  /*var looks like "a\nb\nc\n\n":*/
  char *ch;
  while (*var > '\n') {
    /*Copy the variable up g_to '\n' into g_fbuf:*/
    /* Modification by A.Smirnov: for compatibility with Mathematica*/
    for (*(ch = g_fbuf) = '\0'; ((*ch = *var) > '\n' && (*ch = *var) != ',');
         ch++, var++)
      if (ch - g_fbuf >= FROMFERMATBUFSIZE) {
        throw ExceptionFermat("'FROMFERMATBUFSIZE' to small");
      }
    *ch = '\0';
    if (*var != '\0') var++;
    /*Set g_fbuf as a polymomial variable:*/

    if (isupper(*g_fbuf)) {
      throw ExceptionFermat("Fermat does not like capital letters in variable names, please rename '" + std::string(g_fbuf) + "'");
    }

    fprintf(g_to, "%s:=%s\n", g_fbuf, numeric.c_str());

    fflush(g_to);
    read_up((char *)const_cast<char *>(int_string(numeric).c_str()),
            int_string(numeric).size());

    if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
      throw ExceptionFermat("Failed to set numerics for variable '" + std::string(g_fbuf) + "'");
    }
  } /*while( *var > '\n')*/
}

void Fermat::unset_numeric(char *pvars_) {
  char *var = pvars_;

  /*Set polymomial variables:*/
  /*var looks like "a\nb\nc\n\n":*/
  char *ch;
  while (*var > '\n') {
    /*Copy the variable up g_to '\n' into g_fbuf:*/
    /* Modification by A.Smirnov: for compatibility with Mathematica*/
    for (*(ch = g_fbuf) = '\0'; ((*ch = *var) > '\n' && (*ch = *var) != ',');
         ch++, var++)
      if (ch - g_fbuf >= FROMFERMATBUFSIZE) {
        throw ExceptionFermat("'FROMFERMATBUFSIZE' to small");
      }
    *ch = '\0';
    if (*var != '\0') var++;
    /*Set g_fbuf as a polymomial variable:*/

    if (isupper(*g_fbuf)) {
      throw ExceptionFermat("Fermat does not like capital letters in variable names, please rename '" + std::string(g_fbuf) + "'");
    }

    fprintf(g_to, "@(%s)\n", g_fbuf);

    fflush(g_to);
    read_up((char *)cancel_ret.c_str(), cancel_ret.length());

    if (fgets(g_fbuf, FROMFERMATBUFSIZE, g_from) == NULL) {
      throw ExceptionFermat("Failed to unset numerics for variable '" + std::string(g_fbuf) + "'");
    }
  } /*while( *var > '\n')*/
}

/*This function places one char to the output buffer with possible
  expansion of the buffer:*/
void Fermat::add_to_out(char ch) {
  if (g_fullout >= g_stopout) {
    unsigned l = g_stopout - g_baseout + DELTA_OUT;
    char *ptr = new char[l];
    memcpy(ptr, g_baseout, g_stopout - g_baseout);

    if (ptr == NULL) {
      throw ExceptionFermat("Memory allocation failed");
    }
    g_fullout = ptr + (g_fullout - g_baseout);
    g_stopout = ptr + l;
    delete[] g_baseout;
    g_baseout = ptr;
  } /*if (g_fullout >= g_stopout)*/
  *g_fullout++ = ch;

} /*addtoout*/

/*reads the stream 'from' up to the line 'terminator' (only 'thesize' first
  characters are compared):*/
void Fermat::read_up(char *terminator, int thesize) {
  char *c;
  for (;;) {
    do {
      for (c = fgets(g_fbuf, FROMFERMATBUFSIZE, g_from); *c <= ' '; c++)
        if (*c == '\0') break;
    } while (*c <= ' ');
    if (strncmp(terminator, c, thesize) == 0) return;
  }
} /*read_up*/

void Fermat::fermat_collect(char *buf1) { fputs(buf1, g_to); }

int Fermat::fermat_calc(int optional) {
  char *c;

  putc('\n', g_to); /*stroke the line*/
  fflush(g_to);     /*Now Fermat starts to work*/

  do {
    c = fgets(g_fbuf, FROMFERMATBUFSIZE, g_from);
    for (; *c <= ' '; c++) {
      if (*c == '\0') break;
    }
  } while ((*c <= ' '));

  int success = 1;
  stringstream error1;
  stringstream error2;
  if (*c == '*') {
    std::unique_lock<std::mutex> lckCritical(mutexFermat, std::defer_lock);
    // critical section mutexFermat
    lckCritical.lock();

    error1 << "Interface to Fermat:  Fermat error: \n";
    do {
      /*ignore '`' and spaces:*/
      for (; *c != '\n'; c++)
        switch (*c) {
          default:
            error2 << *c;
        } /*for(;*c!='\n';c++)switch(*c)*/
    } while (*(c = fgets(g_fbuf, FROMFERMATBUFSIZE, g_from)) != '\n');

    if (error2.str().substr(0, 46) ==
        "*** Inappropriate symbol: / can't divide by 0.") {
      success = 2;
    }
    else {
      success = 0;
    }
    lckCritical.unlock();
  }

  if (success == 1) {
    /*initialize the output buffer:*/
    g_fullout = g_baseout;


    if (c != NULL) {
      int skip = 0;
      do {
        skip = 0;

        /*ignore '`' and spaces:*/
        for (; *c != '\n'; c++){
          if (*c == '\0') {
//             if(strlen(g_fbuf)==FROMFERMATBUFSIZE-1)
              skip=1;
            break;
          }
          switch (*c) {
            case ' ':
            case '`':
              continue;
            default:
              add_to_out(*c);
          } /*for(;*c!='\n';c++)switch(*c)*/
        }
      } while (*(c = fgets(g_fbuf, FROMFERMATBUFSIZE, g_from)) != '\n' || skip==1);
      /*the empty line is the Fermat prompt*/
    }

    add_to_out('\0'); /*Complete the line*/

    return 1;
  }
  if (success == 0) {
    throw ExceptionFermat(error1.str() + " " + error2.str());
  }
  if (success == 2 && optional == 1) {
    //     c = fgets(g_fbuf,FROMFERMATBUFSIZE,g_from);
    //     do{
    //       for(;*c!='\n';c++);
    //     } while( *(c=fgets(g_fbuf,FROMFERMATBUFSIZE,g_from))!='\n' );
    c = fgets(g_fbuf, FROMFERMATBUFSIZE, g_from);
    do {
      /*ignore '`' and spaces:*/
      for (; *c != '\n'; c++)
        switch (*c) {
          case ' ':
          case '`':
            continue;
          default:
            error2 << *c;
        } /*for(;*c!='\n';c++)switch(*c)*/
    } while (*(c = fgets(g_fbuf, FROMFERMATBUFSIZE, g_from)) != '\n');

    //     logger << error1.str()<< error2.str() << "\n";
    fputs("123456789101112131415", g_to);
    putc('\n', g_to);
    fflush(g_to);
    read_up((char *)"123456789101112131415", 21);

    return 0;
  }
  else {
    throw ExceptionFermat(error1.str() + " " + error2.str());
  }
}

/*mustCleanup == 0 -- no allocated memory free:*/
void Fermat::close_calc(int mustCleanup) {
  // check if already closed before
  if (!closed) {
    closed = true;
  }
  else {
    return;
  }

  /*Stop Fermat:*/
  fputs("&q\n", g_to);
  fflush(g_to);
  fclose(g_to);
  fclose(g_from);

  // give Fermat some time to properly close itself
  std::this_thread::sleep_for(1s);

  // check if Fermat closed itself properly
  pid_t tmp = waitpid(g_childpid, &status, WNOHANG);

  if (tmp == g_childpid) {
    if (!WIFEXITED(status) && !WIFSIGNALED(status)) {
      // Fermat still runs, force termination
      kill(g_childpid, SIGKILL);
      waitpid(g_childpid, &status, 0);
    }
  }
  else if (tmp == 0) {
    // Fermat still runs, force termination
    kill(g_childpid, SIGKILL);
    waitpid(g_childpid, &status, 0);
  }

  if (mustCleanup) {
    delete[] g_baseout;
    g_baseout = NULL, g_fullout = NULL, g_stopout = NULL;
  }

} /*closeCalc*/
