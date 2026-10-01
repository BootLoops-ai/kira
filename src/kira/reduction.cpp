/*
Copyright (C) 2017-2025 The Kira Developers (see AUTHORS file)
Modified 2026 for BootLoops by Matthew D. Schwartz; distributed under the GNU GPL v3 or later, see COPYING. The changes are described in PATCHES.md.

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

#include <pthread.h>

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <sstream>
#include <tuple>
#include <iomanip>

#include "kira/connect2kira.h"
#include "kira/createSeeds.h"
#include "kira/dataBase.h"
#include "kira/kira.h"
#include "kira/tools.h"

#include "pyred/defs.h"
#include "pyred/integrals.h"

using namespace std;
using namespace GiNaC;

static Loginfo& logger = Loginfo::instance();

static pthread_mutex_t mPreserveFermat;
static pthread_cond_t cPreserveFermat;

static pthread_mutex_t mPreserveCombine;
static pthread_cond_t cPreserveCombine;

static pthread_mutex_t mPreserveResults;
static pthread_cond_t cPreserveResults;

std::mutex mutexCritical;

string build_string(list<string> ptrToString) {
  string str = "";
  for (auto it = ptrToString.begin(); it != ptrToString.end(); ++it) {
    if (it != ptrToString.begin()) {
      //       fermat[id].fermat_collect(const_cast<char*>("+"));
      str += "+";
    }
    //     fermat[id].fermat_collect(const_cast<char*>((*it).c_str()));
    str += (*it);
  }
  return str;
}

void Algebra::init_numbers(
    std::vector<int>& numbersINIT,
    std::vector<std::vector<std::string> >& numbersDenStrVec) {
  size_t sizeNumbers = numbersINIT.size();

  for (size_t itS = 1; itS < sizeNumbers; itS++) {
    vector<uint64_t> numbersDen;
    numbersDen.push_back(numbersINIT[itS] - numbersINIT[itS - 1]);

    vector<string> numbersDenStr;

    pthread_mutex_lock(&mPreserveFermat);
    while ((*kira).idFermats.empty())
      pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
    int id = (*kira).idFermats.top();
    (*kira).idFermats.pop();
    pthread_mutex_unlock(&mPreserveFermat);

    (*kira).fermat[id].fermat_collect(
        const_cast<char*>(int_string(numbersINIT[itS]).c_str()));
    (*kira).fermat[id].fermat_collect(
        const_cast<char*>(("-" + int_string(numbersINIT[itS - 1])).c_str()));
    (*kira).fermat[id].fermat_calc();

    string testTmp = (*kira).fermat[id].g_baseout;
    numbersDenStr.push_back(testTmp);

    pthread_mutex_lock(&mPreserveFermat);
    (*kira).idFermats.push(id);
    pthread_cond_signal(&cPreserveFermat);
    pthread_mutex_unlock(&mPreserveFermat);

    for (size_t it = 1; it < itS; it++) {
      pthread_mutex_lock(&mPreserveFermat);
      while ((*kira).idFermats.empty())
        pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
      int id = (*kira).idFermats.top();
      (*kira).idFermats.pop();
      pthread_mutex_unlock(&mPreserveFermat);

      numbersDen.push_back(numbersDen.back() *
                           (numbersINIT[itS] - numbersINIT[itS - 1 - it]));

      (*kira).fermat[id].fermat_collect(
          const_cast<char*>(numbersDenStr.back().c_str()));
      (*kira).fermat[id].fermat_collect(
          const_cast<char*>(("*(" + int_string(numbersINIT[itS]) + "-" +
                             int_string(numbersINIT[itS - 1 - it]) + ")")
                                .c_str()));

      (*kira).fermat[id].fermat_calc();

      string testTmp = (*kira).fermat[id].g_baseout;
      numbersDenStr.push_back(testTmp);

      pthread_mutex_lock(&mPreserveFermat);
      (*kira).idFermats.push(id);
      pthread_cond_signal(&cPreserveFermat);
      pthread_mutex_unlock(&mPreserveFermat);
    }
    numbersDenStrVec.push_back(numbersDenStr);
  }
}

uint64_t countterms(string& token) {
  char checkPlus = '+';
  char checkMinus = '-';
  char exclude = '(';
  uint64_t count = 0;

  for (size_t i = 1; i < token.size(); i++) {
    if ((token[i] == checkPlus || token[i] == checkMinus) &&
        token[i - 1] != exclude) {
      ++count;
    }
  }
  return count;
}

int Algebra::get_power_level(std::string& numeratorPart) {
  int value1 = 0;
  int skip = 0;
  int tmpValue1 = 0;
  {
    string foundStr1;
    size_t posExp1;

    vector<string> testString = {"*" + interpVar + "^", "+" + interpVar + "^",
                                 "-" + interpVar + "^", "(" + interpVar + "^"};

    for (auto it : testString) {
      size_t posVar = numeratorPart.find(it);

      if (posVar != std::string::npos) {

        foundStr1 = numeratorPart.substr(posVar + 2 + (interpVar.size()));
        posExp1 = foundStr1.find_first_not_of("0123456789");

        if (posExp1 != std::string::npos) {
          tmpValue1 = stoi(foundStr1.substr(0, posExp1));
          skip = 1;
          if (tmpValue1 > value1) value1 = tmpValue1;
        }
      }
    }
  }

  if (!skip) {
    string foundStr1;
    size_t posExp1;

    vector<string> testString = {"*" + interpVar, "+" + interpVar,
                                 "-" + interpVar, "(" + interpVar};

    for (auto it : testString) {
      size_t posVar = numeratorPart.find(it);

      if (posVar != std::string::npos) {
        foundStr1 = numeratorPart.substr(posVar + 1 + (interpVar.size()));

        posExp1 = foundStr1.find_first_of(")+-");

        if (posExp1 == 0) {
          value1 = 1;
          break;
        }
      }
    }
  }

  return value1;
}

std::pair<int, int> Algebra::get_power(std::string& stringofinterest) {
  string numeratorPart;
  string denominatorPart;

  std::size_t pos;

  pos = stringofinterest.find("/");
  if (pos != std::string::npos) {
    numeratorPart = "(" + stringofinterest.substr(0, pos) + ")";
    denominatorPart = "(" + stringofinterest.substr(pos + 1) + ")";
  }
  else {
    numeratorPart = "(" + stringofinterest + ")";
    denominatorPart = "(1)";
  }

  int value1 = get_power_level(numeratorPart);
  int value2 = get_power_level(denominatorPart);

  return make_pair(value1, value2);
}

bool sortReconst(const pair<string, int>& i, const pair<string, int>& j) {
  return (i.second < j.second);
}

bool myfunction(const string& i, const string& j) {
  return (i.length() < j.length());
}

std::pair<std::string, std::string> Algebra::normsample(
    std::string stringofinterest, int& valueForD, int& degree, int& powerLevel1,
    int& powerLevel2) {

  vector<pair<string, int> > normalization;

  uint32_t countLoop = 0;
  uint32_t loopSIZE = 9;
  std::mutex m;
  std::condition_variable cond_var;
  bool processed = false;

  for (uint32_t itR = 0; itR < loopSIZE; itR++) {
    (*kira).pool5[poolLevel].enqueue([this, itR, &valueForD, &stringofinterest,
                                      &normalization, &countLoop, &processed,
                                      &loopSIZE, &cond_var, &m]() {
      char var[100];
      int sampleVariable = valueForD + itR;

      pthread_mutex_lock(&mPreserveFermat);
      while ((*kira).idFermats.empty())
        pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
      int id = (*kira).idFermats.top();
      (*kira).idFermats.pop();
      pthread_mutex_unlock(&mPreserveFermat);

      var[0] = '\0';
      snprintf(var, sizeof(var), "%s\n", "d");
      (*kira).fermat[id].unset_variable(var);
      (*kira).fermat[id].set_numeric(var, sampleVariable);

      for (auto itX = numericVar.rbegin(); itX != numericVar.rend(); itX++) {
        var[0] = '\0';
        snprintf(var, sizeof(var), "%s\n", itX->first.c_str());
        (*kira).fermat[id].unset_variable(var);
        (*kira).fermat[id].set_numeric(var, itX->second);
      }

      for (size_t itX = 0; itX < reductVar.size(); itX++) {
        var[0] = '\0';
        snprintf(var, sizeof(var), "%s\n", reductVar[itX].c_str());
        (*kira).fermat[id].unset_variable(var);
        (*kira).fermat[id].set_numeric(var, 3 * itX + 5 * sampleVariable);
      }

      (*kira).fermat[id].fermat_collect(
          const_cast<char*>(stringofinterest.c_str()));
      if ((*kira).fermat[id].fermat_calc(1)) {
        {
          std::lock_guard<std::mutex> lock(m);

          string tmp = (*kira).fermat[id].g_baseout;
          normalization.push_back(make_pair(tmp, sampleVariable));
        }
      }

      var[0] = '\0';
      snprintf(var, sizeof(var), "%s\n", "d");
      (*kira).fermat[id].unset_numeric(var);
      (*kira).fermat[id].set_variable(var);

      for (size_t itX = 0; itX < reductVar.size(); itX++) {
        var[0] = '\0';
        snprintf(var, sizeof(var), "%s\n", reductVar[itX].c_str());
        (*kira).fermat[id].unset_numeric(var);
        (*kira).fermat[id].set_variable(var);
      }

      for (auto itX = numericVar.rbegin(); itX != numericVar.rend(); itX++) {
        var[0] = '\0';
        snprintf(var, sizeof(var), "%s\n", itX->first.c_str());
        (*kira).fermat[id].unset_numeric(var);
        (*kira).fermat[id].set_variable(var);
      }

      var[0] = '\0';
      snprintf(var, sizeof(var), "%s\n", interpVar.c_str());
      (*kira).fermat[id].unset_variable(var);
      (*kira).fermat[id].set_variable(var);

      pthread_mutex_lock(&mPreserveFermat);
      (*kira).idFermats.push(id);
      pthread_cond_signal(&cPreserveFermat);
      pthread_mutex_unlock(&mPreserveFermat);
      {
        std::lock_guard<std::mutex> lock(m);
        countLoop++;
        if (countLoop == loopSIZE) {
          processed = true;
          cond_var.notify_one();
        }
      }
    });
  }

  {
    std::unique_lock<std::mutex> lock(m);
    cond_var.wait(lock, [&processed] { return processed; });
  }

  sort(normalization.begin(), normalization.end(), sortReconst);

  std::size_t pos;

  vector<string> numNormalization;
  vector<string> denNormalization;

  vector<pair<uint64_t, uint64_t> > nOfTerms;

  for (auto itX : normalization) {
    pos = itX.first.find("/");
    if (pos != std::string::npos) {
      numNormalization.push_back(itX.first.substr(0, pos));
      denNormalization.push_back(itX.first.substr(pos + 1));
    }
    else {
      numNormalization.push_back(itX.first);
      denNormalization.push_back("1");
    }
    nOfTerms.push_back(make_pair(countterms(numNormalization.back()),
                                 countterms(denNormalization.back())));
  }

  uint64_t tmp1 = nOfTerms.front().first;
  uint64_t tmp2 = nOfTerms.front().second;
  int makeChoice = 0;
  for (size_t itX = 1; itX < nOfTerms.size(); itX++) {
    if (nOfTerms[itX].first != tmp1 || nOfTerms[itX].second != tmp2) {
      if (nOfTerms[itX].first >= tmp1 && nOfTerms[itX].second >= tmp2) {
        tmp1 = nOfTerms[itX].first;
        tmp2 = nOfTerms[itX].second;
        makeChoice = itX;
      }
    }
  }

  valueForD = normalization[makeChoice].second;

  powerLevel1 = get_power_level(numNormalization[makeChoice]);

  powerLevel2 = get_power_level(denNormalization[makeChoice]);

  degree = powerLevel1 > powerLevel2 ? powerLevel1 : powerLevel2;
  degree++; // t^0
  degree++; // additional to control things
  if (degree < 10) degree = 10;

  return make_pair(numNormalization[makeChoice], denNormalization[makeChoice]);
}

std::tuple<std::vector<std::string>, std::vector<std::string>,
           std::vector<int> >
Algebra::sample(std::vector<int>& numbersINIT, std::string& stringofinterest,
                std::pair<std::string, std::string>& normsampleVar,
                uint64_t& tmp1, uint64_t& tmp2) {
  vector<string> numerator;
  vector<string> denominator;

  char var[100];

  var[0] = '\0';
  snprintf(var, sizeof(var), "%s\n", interpVar.c_str());

  vector<pair<string, int> > samples;
  vector<int> dummyNumbers;
  vector<int> erase;

  uint32_t countLoop = 0;
  uint32_t loopSIZE = numbersINIT.size();
  std::mutex m;
  std::condition_variable cond_var;
  bool processed = false;

  for (size_t i = 0; i < numbersINIT.size(); i++) {
    (*kira).pool5[poolLevel].enqueue([this, i, &numbersINIT, &var,
                                      &stringofinterest, &dummyNumbers,
                                      &samples, &erase, &countLoop, &processed,
                                      &loopSIZE, &cond_var, &m]() {

      if (reductVar.size() > 0 && depthLevel > 0) {

        auto parseNumericVar = numericVar;
        parseNumericVar.push_back(make_pair(interpVar, numbersINIT[i]));

        Algebra reconst(kira, interpVar, reductVar, parseNumericVar,
                        depthLevel);
        string plup = "0";
        auto testTmp = reconst.reconstruct_final(stringofinterest, plup);

        {
          std::lock_guard<std::mutex> lock(m);
          samples.push_back(make_pair(testTmp, i));
        }
      }
      else {
        pthread_mutex_lock(&mPreserveFermat);
        while ((*kira).idFermats.empty())
          pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
        int id = (*kira).idFermats.top();
        (*kira).idFermats.pop();
        pthread_mutex_unlock(&mPreserveFermat);

        char varN[100];

        for (auto itX = numericVar.rbegin(); itX != numericVar.rend(); itX++) {
          varN[0] = '\0';
          snprintf(varN, sizeof(varN), "%s\n", itX->first.c_str());
          (*kira).fermat[id].unset_variable(varN);
          (*kira).fermat[id].set_numeric(varN, itX->second);
        }

        (*kira).fermat[id].unset_variable(var);
        (*kira).fermat[id].set_numeric(var, numbersINIT[i]);

        (*kira).fermat[id].fermat_collect(
            const_cast<char*>((stringofinterest).c_str()));
        if (!(*kira).fermat[id].fermat_calc(1)) {
          {
            std::lock_guard<std::mutex> lock(m);
            if (dummyNumbers.size() == 0)
              dummyNumbers.push_back(numbersINIT.back() + 1);
            else
              dummyNumbers.push_back(dummyNumbers.back() + 1);

            erase.push_back(i);
          }
        }
        else {
          {
            std::lock_guard<std::mutex> lock(m);
            string testTmp = (*kira).fermat[id].g_baseout;
            samples.push_back(make_pair(testTmp, i));
          }
        }
        (*kira).fermat[id].unset_numeric(var);
        (*kira).fermat[id].set_variable(var);

        for (auto itX = numericVar.rbegin(); itX != numericVar.rend(); itX++) {
          varN[0] = '\0';
          snprintf(varN, sizeof(varN), "%s\n", itX->first.c_str());
          (*kira).fermat[id].unset_numeric(varN);
          (*kira).fermat[id].set_variable(varN);
        }

        pthread_mutex_lock(&mPreserveFermat);
        (*kira).idFermats.push(id);
        pthread_cond_signal(&cPreserveFermat);
        pthread_mutex_unlock(&mPreserveFermat);
      }

      {
        std::lock_guard<std::mutex> lock(m);
        countLoop++;
        if (countLoop == loopSIZE) {
          processed = true;
          cond_var.notify_one();
        }
      }
    });
  }

  {
    std::unique_lock<std::mutex> lock(m);
    cond_var.wait(lock, [&processed] { return processed; });
  }

  sort(dummyNumbers.begin(), dummyNumbers.end());
  sort(samples.begin(), samples.end(), sortReconst);

  for (auto itX : erase)
    numbersINIT.erase(numbersINIT.begin() + itX);

  vector<pair<uint64_t, uint64_t> > nOfTerms;

  std::size_t pos;

  for (auto itsamples : samples) {
    pos = itsamples.first.find(")/(");
    if (pos != std::string::npos) {
      numerator.push_back(itsamples.first.substr(0, pos + 1));
      denominator.push_back(itsamples.first.substr(pos + 2));
    }
    else {
      numerator.push_back(itsamples.first);
      denominator.push_back("1");
    }
    nOfTerms.push_back(make_pair(countterms(numerator.back()),
                                 countterms(denominator.back())));
  }

  if (tmp1 == 0 && tmp2 == 0) {
    tmp1 = nOfTerms.front().first;
    tmp2 = nOfTerms.front().second;
  }

  if (tmp1 < nOfTerms.front().first || tmp1 < nOfTerms.front().first) {
    throw ExceptionInternal("Algebra::sample(): Check tmp1 and tmp2");
  }

  for (size_t itX = 1; itX < nOfTerms.size(); itX++) {
    if (nOfTerms[itX].first != tmp1 || nOfTerms[itX].second != tmp2) {
      if (nOfTerms[itX].first > tmp1 || nOfTerms[itX].second > tmp2) {
        tmp1 = nOfTerms[itX].first;
        tmp2 = nOfTerms[itX].second;
      }
    }
  }

  vector<int> collectBadPoints;
  for (size_t itX = 0; itX < nOfTerms.size(); itX++) {
    if (tmp1 != nOfTerms[itX].first || tmp2 != nOfTerms[itX].second) {
      collectBadPoints.push_back(itX);
    }
  }

  int shift = 0;

  for (size_t itX = 0; itX < collectBadPoints.size(); itX++) {
    if (dummyNumbers.size() == 0)
      dummyNumbers.push_back(numbersINIT.back() + 1);
    else
      dummyNumbers.push_back(dummyNumbers.back() + 1);

    numbersINIT.erase(numbersINIT.begin() + collectBadPoints[itX] - shift);
    numerator.erase(numerator.begin() + collectBadPoints[itX] - shift);
    denominator.erase(denominator.begin() + collectBadPoints[itX] - shift);
    shift++;
  }

  dummyNumbers = check_sample_points(normsampleVar, dummyNumbers);

  return make_tuple(numerator, denominator, dummyNumbers);
}

std::vector<std::string> Algebra::normalize(std::vector<int>& numbersINIT,
                                            int valueForD,
                                            std::vector<std::string>& numerator,
                                            std::string& numNormalization) {
  char var[100];

  vector<string> prefactor;

  pthread_mutex_lock(&mPreserveFermat);
  while ((*kira).idFermats.empty())
    pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
  int id = (*kira).idFermats.top();
  (*kira).idFermats.pop();
  pthread_mutex_unlock(&mPreserveFermat);

  var[0] = '\0';
  snprintf(var, sizeof(var), "%s\n", "d");
  (*kira).fermat[id].unset_variable(var);
  (*kira).fermat[id].set_numeric(var, valueForD);

  for (size_t itX = 0; itX < reductVar.size(); itX++) {
    var[0] = '\0';
    snprintf(var,sizeof(var), "%s\n", reductVar[itX].c_str());
    (*kira).fermat[id].unset_variable(var);
    (*kira).fermat[id].set_numeric(var, 3 * itX + 5 * valueForD);
  }

  for (auto itX = numericVar.rbegin(); itX != numericVar.rend(); itX++) {
    var[0] = '\0';
    snprintf(var, sizeof(var), "%s\n", itX->first.c_str());
    (*kira).fermat[id].unset_variable(var);
    (*kira).fermat[id].set_numeric(var, itX->second);
  }

  var[0] = '\0';
  snprintf(var, sizeof(var), "%s\n", interpVar.c_str());
  (*kira).fermat[id].unset_variable(var);

  for (size_t i = 0; i < numbersINIT.size(); i++) {
    (*kira).fermat[id].set_numeric(var, numbersINIT[i]);
    (*kira).fermat[id].fermat_collect(const_cast<char*>(
        ("(" + numerator[i] + ")/(" + numNormalization + ")").c_str()));
    (*kira).fermat[id].fermat_calc();
    string testTmp = (*kira).fermat[id].g_baseout;

    prefactor.push_back(testTmp);
  }

  var[0] = '\0';
  snprintf(var, sizeof(var), "%s\n", "d");
  (*kira).fermat[id].unset_numeric(var);
  (*kira).fermat[id].set_variable(var);

  for (size_t itX = 0; itX < reductVar.size(); itX++) {
    var[0] = '\0';
    snprintf(var, sizeof(var), "%s\n", reductVar[itX].c_str());
    (*kira).fermat[id].unset_numeric(var);
    (*kira).fermat[id].set_variable(var);
  }

  var[0] = '\0';
  snprintf(var, sizeof(var), "%s\n", interpVar.c_str());
  (*kira).fermat[id].unset_numeric(var);
  (*kira).fermat[id].set_variable(var);

  for (auto itX = numericVar.rbegin(); itX != numericVar.rend(); itX++) {
    var[0] = '\0';
    snprintf(var, sizeof(var), "%s\n", itX->first.c_str());
    (*kira).fermat[id].unset_numeric(var);
    (*kira).fermat[id].set_variable(var);
  }

  pthread_mutex_lock(&mPreserveFermat);
  (*kira).idFermats.push(id);
  pthread_cond_signal(&cPreserveFermat);
  pthread_mutex_unlock(&mPreserveFermat);

  for (size_t it = 0; it < prefactor.size(); it++) {
    numerator[it] = "(" + numerator[it] + ")/(" + prefactor[it] + ")";
  }
  return numerator;
}

std::pair<std::string, int> Algebra::reconstruct_function(
    std::vector<int>& numbersINIT,
    std::vector<std::vector<std::string> >& numbersDenStrVec,
    std::vector<std::string>& numerator, int& powerLevel) {
  string reconstructed;

  vector<string> interpolation;

  interpolation.push_back("(" + numerator.front() + ")");

  for (int it = 1; it < /*numerator.size()*/ powerLevel + 2; it++) {
    pthread_mutex_lock(&mPreserveFermat);
    while ((*kira).idFermats.empty())
      pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
    int id = (*kira).idFermats.top();
    (*kira).idFermats.pop();
    pthread_mutex_unlock(&mPreserveFermat);

    string tokenStr = "";

    for (size_t itN = 0; itN < numbersDenStrVec[it - 1].size(); itN++) {
      if (itN == 0) {
        tokenStr +=
            ("(" + numerator[it] + "-" + interpolation[0] + ")/" +
             numbersDenStrVec[it - 1][numbersDenStrVec[it - 1].size() - 1]);

        (*kira).fermat[id].fermat_collect(const_cast<char*>(
            ("(" + numerator[it] + "-" + interpolation[0] + ")/" +
             numbersDenStrVec[it - 1][numbersDenStrVec[it - 1].size() - 1])
                .c_str()));
      }
      else {
        tokenStr += ("-(" + interpolation[itN] + ")/" +
                     numbersDenStrVec[it - 1][numbersDenStrVec[it - 1].size() -
                                              1 - itN]);

        (*kira).fermat[id].fermat_collect(const_cast<char*>(
            ("-(" + interpolation[itN] + ")/" +
             numbersDenStrVec[it - 1]
                             [numbersDenStrVec[it - 1].size() - 1 - itN])
                .c_str()));
      }
    }
    (*kira).fermat[id].fermat_calc();

    string testTmp = (*kira).fermat[id].g_baseout;

    pthread_mutex_lock(&mPreserveFermat);
    (*kira).idFermats.push(id);
    pthread_cond_signal(&cPreserveFermat);
    pthread_mutex_unlock(&mPreserveFermat);

    interpolation.push_back(testTmp);
  }

  if (interpolation.back() != "0") return make_pair("0", 0);
  pthread_mutex_lock(&mPreserveFermat);
  while ((*kira).idFermats.empty())
    pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
  int id = (*kira).idFermats.top();
  (*kira).idFermats.pop();
  pthread_mutex_unlock(&mPreserveFermat);

  for (size_t it = 0; it < interpolation.size(); it++) {
    if (it != 0) {
      (*kira).fermat[id].fermat_collect(const_cast<char*>("+"));
    }
    (*kira).fermat[id].fermat_collect(
        const_cast<char*>(("(" + interpolation[it] + ")").c_str()));

    //     string coef="";
    //     coef+="("+interpolation[it]+")";

    for (size_t itS = 0; itS < it; itS++) {
      (*kira).fermat[id].fermat_collect(const_cast<char*>(
          ("*(" + interpVar + "-" + int_string(numbersINIT[itS]) + ")")
              .c_str()));

      //       coef+="*("+interpVar+"-"+int_string(numbersINIT[itS])+")";
    }
    //     ptrToString.push_front(coef);
  }

  (*kira).fermat[id].fermat_calc();
  reconstructed = (*kira).fermat[id].g_baseout;

  pthread_mutex_lock(&mPreserveFermat);
  (*kira).idFermats.push(id);
  pthread_cond_signal(&cPreserveFermat);
  pthread_mutex_unlock(&mPreserveFermat);

  return make_pair(reconstructed, 1);
}

std::vector<int> Algebra::check_sample_points(
    std::pair<std::string, std::string>& normsampleVar,
    std::vector<int>& points) {
  char var[100];
  vector<int> numbersINIT;

  for (size_t i = 0; i < points.size(); i++) {
    pthread_mutex_lock(&mPreserveFermat);
    while ((*kira).idFermats.empty())
      pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
    int id = (*kira).idFermats.top();
    (*kira).idFermats.pop();
    pthread_mutex_unlock(&mPreserveFermat);

    var[0] = '\0';
    snprintf(var, sizeof(var), "%s\n", interpVar.c_str());
    (*kira).fermat[id].unset_variable(var);

    (*kira).fermat[id].set_numeric(var, points[i]);
    (*kira).fermat[id].fermat_collect(
        const_cast<char*>(normsampleVar.first.c_str()));

    (*kira).fermat[id].fermat_calc();

    string numTmp = (*kira).fermat[id].g_baseout;

    (*kira).fermat[id].fermat_collect(
        const_cast<char*>(normsampleVar.second.c_str()));
    (*kira).fermat[id].fermat_calc();

    string denTmp = (*kira).fermat[id].g_baseout;

    if (numTmp != "0" && denTmp != "0") {
      numbersINIT.push_back(points[i]);
    }
    else {
      points.push_back(points.back() + 1);
    }
    (*kira).fermat[id].unset_numeric(var);
    (*kira).fermat[id].set_variable(var);

    pthread_mutex_lock(&mPreserveFermat);
    (*kira).idFermats.push(id);
    pthread_cond_signal(&cPreserveFermat);
    pthread_mutex_unlock(&mPreserveFermat);
  }

  return numbersINIT;
}

std::vector<int> Algebra::generate_sample_points(
    std::pair<std::string, std::string>& normsampleVar, int degree) {
  char var[100];
  vector<int> numbersINIT;

  for (int i = 1; i < degree + 1; i++) {
    pthread_mutex_lock(&mPreserveFermat);
    while ((*kira).idFermats.empty())
      pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
    int id = (*kira).idFermats.top();
    (*kira).idFermats.pop();
    pthread_mutex_unlock(&mPreserveFermat);

    var[0] = '\0';
    snprintf(var, sizeof(var), "%s\n", interpVar.c_str());
    (*kira).fermat[id].unset_variable(var);

    (*kira).fermat[id].set_numeric(var, i * 7);
    (*kira).fermat[id].fermat_collect(
        const_cast<char*>(normsampleVar.first.c_str()));

    (*kira).fermat[id].fermat_calc();

    string numTmp = (*kira).fermat[id].g_baseout;

    (*kira).fermat[id].fermat_collect(
        const_cast<char*>(normsampleVar.second.c_str()));
    (*kira).fermat[id].fermat_calc();

    string denTmp = (*kira).fermat[id].g_baseout;

    if (numTmp != "0" && denTmp != "0") {
      numbersINIT.push_back(i * 7);
    }
    else
      ++degree;

    (*kira).fermat[id].unset_numeric(var);
    (*kira).fermat[id].set_variable(var);

    pthread_mutex_lock(&mPreserveFermat);
    (*kira).idFermats.push(id);
    pthread_cond_signal(&cPreserveFermat);
    pthread_mutex_unlock(&mPreserveFermat);
  }

  return numbersINIT;
}

std::string Algebra::reconstruct_final(std::string& stringofinterest1,
                                       std::string& /*stringofinterest2*/) {

  string stringofinterest = stringofinterest1 /*+ "+" + stringofinterest2*/;

  int valueForD = 7;
  int degree;
  int powerLevel1;
  int powerLevel2;

  pair<string, int> reconstrNum, reconstrDen;
  int countRec = 0;

  while (1) {
    Clock clockNorm;

    auto normsampleVar = normsample(stringofinterest, valueForD, degree,
                                    powerLevel1, powerLevel2);


    if (normsampleVar.first == "0") {
      return "0";
    }

    vector<int> numbersINIT2 = generate_sample_points(normsampleVar, degree);

    degree = powerLevel1 > powerLevel2 ? powerLevel1 : powerLevel2;
    degree++; // t^0
    degree++;

    vector<int> numbersINIT;
    for (int itX = 0; itX < degree; itX++) {
      numbersINIT.push_back(numbersINIT2[itX]);
    }

    if (degree == 2) {
      pthread_mutex_lock(&mPreserveFermat);
      while ((*kira).idFermats.empty())
        pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
      int id = (*kira).idFermats.top();
      (*kira).idFermats.pop();
      pthread_mutex_unlock(&mPreserveFermat);

      (*kira).fermat[id].fermat_collect(
          const_cast<char*>(stringofinterest.c_str()));
      (*kira).fermat[id].fermat_calc();
      string str = (*kira).fermat[id].g_baseout;

      pthread_mutex_lock(&mPreserveFermat);
      (*kira).idFermats.push(id);
      pthread_cond_signal(&cPreserveFermat);
      pthread_mutex_unlock(&mPreserveFermat);

      return str;
    }

    uint64_t tmp1 = 0;
    uint64_t tmp2 = 0;
    Clock clock;
    auto sampleVar =
        sample(numbersINIT, stringofinterest, normsampleVar, tmp1, tmp2);

    vector<int> samplePoints = get<2>(sampleVar);

    pair<vector<string>, vector<string> > sampleCoef =
        make_pair(get<0>(sampleVar), get<1>(sampleVar));

    //int count = 0;
    while (samplePoints.size() > 0) {
      sampleVar =
          sample(samplePoints, stringofinterest, normsampleVar, tmp1, tmp2);

      for (auto itX : samplePoints)
        numbersINIT.push_back(itX);

      for (auto itX : get<0>(sampleVar))
        sampleCoef.first.push_back(itX);

      for (auto itX : get<1>(sampleVar))
        sampleCoef.second.push_back(itX);

      samplePoints = get<2>(sampleVar);
      //count++;
    }

    std::vector<std::string> numerator, denominator;
    {
      uint32_t countLoop = 0;
      uint32_t loopSIZE = 2;
      std::mutex m;
      std::condition_variable cond_var;
      bool processed = false;

      (*kira).pool5[poolLevel].enqueue(
          [this, valueForD, &numerator, &numbersINIT, &sampleCoef,
           &normsampleVar, &countLoop, &processed, &loopSIZE, &cond_var, &m]() {
            numerator = normalize(numbersINIT, valueForD, sampleCoef.first,
                                  normsampleVar.first);
            {
              std::lock_guard<std::mutex> lock(m);
              countLoop++;
              if (countLoop == loopSIZE) {
                processed = true;
                cond_var.notify_one();
              }
            }
          });

      (*kira).pool5[poolLevel].enqueue(
          [this, valueForD, &denominator, &numbersINIT, &sampleCoef,
           &normsampleVar, &countLoop, &processed, &loopSIZE, &cond_var, &m]() {
            denominator = normalize(numbersINIT, valueForD, sampleCoef.second,
                                    normsampleVar.second);
            {
              std::lock_guard<std::mutex> lock(m);
              countLoop++;
              if (countLoop == loopSIZE) {
                processed = true;
                cond_var.notify_one();
              }
            }
          });

      {
        std::unique_lock<std::mutex> lock(m);
        cond_var.wait(lock, [&processed] { return processed; });
      }
    }

    std::vector<std::vector<std::string> > numbersDenStrVec;
    init_numbers(numbersINIT, numbersDenStrVec);

    {
      uint32_t countLoop = 0;
      uint32_t loopSIZE = 2;
      std::mutex m;
      std::condition_variable cond_var;
      bool processed = false;

      (*kira).pool5[poolLevel].enqueue(
          [this, &reconstrNum, &numerator, &numbersINIT, &numbersDenStrVec,
           &powerLevel1, &countLoop, &processed, &loopSIZE, &cond_var, &m]() {
            reconstrNum = reconstruct_function(numbersINIT, numbersDenStrVec,
                                               numerator, powerLevel1);
            {
              std::lock_guard<std::mutex> lock(m);
              countLoop++;
              if (countLoop == loopSIZE) {
                processed = true;
                cond_var.notify_one();
              }
            }
          });

      (*kira).pool5[poolLevel].enqueue(
          [this, &reconstrDen, &denominator, &numbersINIT, &numbersDenStrVec,
           &powerLevel2, &countLoop, &processed, &loopSIZE, &cond_var, &m]() {
            reconstrDen = reconstruct_function(numbersINIT, numbersDenStrVec,
                                               denominator, powerLevel2);
            {
              std::lock_guard<std::mutex> lock(m);
              countLoop++;
              if (countLoop == loopSIZE) {
                processed = true;
                cond_var.notify_one();
              }
            }
          });

      {
        std::unique_lock<std::mutex> lock(m);
        cond_var.wait(lock, [&processed] { return processed; });
      }
    }

    if (reconstrNum.second == 1 && reconstrDen.second == 1) break;
    //     logger << "run failed, repeat\n";
    valueForD++;
    countRec++;
    if (countRec > 9) {
      throw ExceptionInternal("Algebra::reconstruct_final(): Failed");
    }
  }

  pthread_mutex_lock(&mPreserveFermat);
  while ((*kira).idFermats.empty())
    pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
  int id = (*kira).idFermats.top();
  (*kira).idFermats.pop();
  pthread_mutex_unlock(&mPreserveFermat);

  (*kira).fermat[id].fermat_collect(const_cast<char*>(
      ("(" + reconstrNum.first + ")/(" + reconstrDen.first + ")").c_str()));
  (*kira).fermat[id].fermat_calc();
  string resonstructedResult = (*kira).fermat[id].g_baseout;

  pthread_mutex_lock(&mPreserveFermat);
  (*kira).idFermats.push(id);
  pthread_cond_signal(&cPreserveFermat);
  pthread_mutex_unlock(&mPreserveFermat);


  return resonstructedResult;
}

Algebra::Algebra(Kira* kira_, std::string& /*interpVar_*/,
                 std::vector<std::string>& reductVar_,
                 std::vector<std::pair<std::string, int> >& numericVar_,
                 int level_) {
  kira = kira_;

  if(level_ == 0){
    depthLevel = 0;
  }
  else{
    depthLevel = level_-1;
  }

  if (reductVar_.size() == 0) {
    throw ExceptionInternal("Algebra::reconstruct_final(): reductVar is zero");
  }


  interpVar = reductVar_[reductVar_.size() - 1];

  for (size_t it = 0; it < reductVar_.size() - 1; it++) {
    reductVar.push_back(reductVar_[it]);
  }

  poolLevel = numericVar_.size();

  for (auto it : numericVar_) {
    numericVar.push_back(it);
  }
}

void Kira::initiate_fermat(int kira2Math = 0, int onlyOnePool = 0,
                           int noFermat = 0) {
  pool = new ThreadPool;
  pool->initialize(coreNumber);

  if (!noFermat) {
    char var[100], token[100];

    var[0] = '\0';
    snprintf(var, sizeof(var), "d\n");

    //int countHeuristic = 0;

    for (auto itInvar = invar.begin(); itInvar != invar.end(); ++itInvar) {
      if (something_string(*itInvar) == something_string(massSet2One) &&
          !kira2Math)
        continue;

      if((*itInvar) == get_symbol("d")){
        continue;
      }

      reductVar.push_back(something_string(*itInvar));

      snprintf(token, sizeof(token), "%s\n", something_string(*itInvar).c_str());
      strcat(var, token);
      //countHeuristic++;
    }

    reconstFlag = 0;
    termNumber = 0;

    if (reductVar.size() > 0) {
      interpVar = reductVar.back();
    }

    if (reductVar.size() > 0 && reductVar.size() < 4) {
      heuristic = 55000;
      reconstFlag = 150000;
      termNumber = 150;
    }
    else if (reductVar.size() == 4) {
      heuristic = 10000;
      reconstFlag = 700000;
      termNumber = 150;
    }
    else {
      heuristic = 1000;
      reconstFlag = 700000;
      termNumber = 150;
    }

    if (algebraicReconstruction == 0) {
      reconstFlag = 0;
      logger << "Algebraic reconstruction is switched off.\n";
    }
    else {
      logger << "Algebraic reconstruction is switched on. If the longest\n";
      logger << "coefficient is longer than " << reconstFlag
             << " characters or the sum\n";
      logger << "contains more than " << termNumber
             << " terms, algebraic reconstruction\n";
      logger << "will be used.\n";
    }

    fermat = new Fermat[coreNumber];
    Kira::s_coreNumber = coreNumber; // parallel-treatcoeff: expose pool size to static treatcoeff2
    for (int i = 0; i < coreNumber; i++) {
      fermat[i].start_fermat(fermatPath, var);
      idFermats.push(i);
    }
    /*Do variable replacements*/
    if(symbols2numsFlag){
      for (auto it : setVariable) {
       var[0] = '\0';
        snprintf(var, sizeof(var), "%s", get<0>(it).c_str());
        logger << "set: " << var << " = " << get<1>(it) << "\n";
        // parallel-treatcoeff: apply --set_value to EVERY pool instance, not
        // just [0], since treatcoeff2 now dispatches to fermat[tid].
        for (int i = 0; i < coreNumber; i++) {
          fermat[i].unset_variable(var);
          fermat[i].set_numeric2(var, get<1>(it));
        }
      }
    }
    pthread_mutex_init(&mPreserveFermat, NULL);
    pthread_cond_init(&cPreserveFermat, NULL);

    pthread_mutex_init(&mPreserveCombine, NULL);
    pthread_cond_init(&cPreserveCombine, NULL);

    pthread_mutex_init(&mPreserveResults, NULL);
    pthread_cond_init(&cPreserveResults, NULL);
  }

  if (!onlyOnePool) {
    pool2 = new ThreadPool;
    pool2->initialize(coreNumber);

    pool3 = new ThreadPool[coreNumber];
    for (int i = 0; i < coreNumber; i++) {
      pool3[i].initialize(coreNumber);
      idCombine.push(i);
    }

    pool5 = new ThreadPool[2];
    for (int i = 0; i < 2; i++)
      pool5[i].initialize(coreNumber);

    if (algebraicReconstruction != 0) {
      logger << "If Kira hangs or gives a segmenatation fault,"
             << " an experimental\n feature in Kira is not working.\n";
      if (setVariable.size() == 0) {
        if (reductVar.size() > 0) {
          string stringofinterest1 = "0";
          string stringofinterest2 = "0";
          vector<pair<string, int> > numericVar;
          Algebra reconst(this, interpVar, reductVar, numericVar, 3);
          auto finalResult =
              reconst.reconstruct_final(stringofinterest1, stringofinterest2);
        }
      }
    }
  }
}

void Kira::destroy_fermat(int onlyOnePool = 0, int noFermat = 0) {
  if (!noFermat) {
    int count = 0;
    for (int id = 0; id < coreNumber; ++id) {
      pool->enqueue([id, &count]() {
        fermat[id].close_calc(1);

        pthread_mutex_lock(&mPreserveFermat);
        ++count;
        pthread_cond_signal(&cPreserveFermat);
        pthread_mutex_unlock(&mPreserveFermat);
      });
    }

    pthread_mutex_lock(&mPreserveFermat);
    while (count != coreNumber) {
      pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
    }
    pthread_mutex_unlock(&mPreserveFermat);

    delete[] fermat;

    pthread_mutex_destroy(&mPreserveFermat);
    pthread_cond_destroy(&cPreserveFermat);

    pthread_mutex_destroy(&mPreserveCombine);
    pthread_cond_destroy(&cPreserveCombine);

    pthread_mutex_destroy(&mPreserveResults);
    pthread_cond_destroy(&cPreserveResults);
  }

  delete pool;

  if (!onlyOnePool) {
    delete pool2;
    delete[] pool3;
    delete[] pool5;
  }
  interpVar.clear();
  reductVar.clear();
  idCombine = stack<int>();
  idFermats = stack<int>();
}

void fermat2flushTimes(string& a, string& b, string& c, Fermat& fermat) {
  fermat.fermat_collect(const_cast<char*>("-("));
  fermat.fermat_collect(const_cast<char*>(a.c_str()));
  fermat.fermat_collect(const_cast<char*>(")/("));
  fermat.fermat_collect(const_cast<char*>(b.c_str()));
  fermat.fermat_collect(const_cast<char*>(")*("));
  fermat.fermat_collect(const_cast<char*>(c.c_str()));
  fermat.fermat_collect(const_cast<char*>(")"));
}

void fermat2flushPlus(string& d, Fermat& fermat) {
  fermat.fermat_collect(const_cast<char*>("+("));
  fermat.fermat_collect(const_cast<char*>(d.c_str()));
  fermat.fermat_collect(const_cast<char*>(")"));
}

int Kira::calculate_coefficient_term(list<string>& ptrToString, int idComb) {
  pthread_mutex_t mPreservePool;
  pthread_cond_t cPreservePool;

  pthread_mutex_init(&mPreservePool, NULL);
  pthread_cond_init(&cPreservePool, NULL);

  int numberOfCoefficients = ptrToString.size();
  string mth;
  int id = -1;

  if (numberOfCoefficients == 1) {
    pthread_mutex_lock(&mPreserveFermat);
    while (idFermats.empty())
      pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
    id = idFermats.top();
    idFermats.pop();
    pthread_mutex_unlock(&mPreserveFermat);

    if (id >= 0) {
      fermat[id].fermat_collect(
          const_cast<char*>((ptrToString.front()).c_str()));
      ptrToString.pop_front();
      fermat[id].fermat_calc();
      mth.assign(fermat[id].g_baseout);
      ptrToString.push_front(mth);
      pthread_mutex_lock(&mPreserveFermat);
      idFermats.push(id);
      pthread_cond_signal(&cPreserveFermat);
      pthread_mutex_unlock(&mPreserveFermat);
      return 1;
    }
  }

  int active_threads = 0;

  //   coreNumber = 11;

  unsigned totalCoeff = ((numberOfCoefficients / 3) < coreNumber)
                            ? (numberOfCoefficients / 3)
                            : coreNumber;

  if (totalCoeff == 0) totalCoeff = 1;

  unsigned countcoeff = 0;

  std::mutex m;
  std::condition_variable cond_var;
  bool notified = false;

  for (unsigned itNum = 0; itNum < totalCoeff; itNum++) {
    pool3[idComb].enqueue([this, &countcoeff, totalCoeff, &cond_var, &notified,
                           &m, &numberOfCoefficients, &active_threads,
                           &ptrToString]() {
      int id = -1;
      int leave = 0;
      while (numberOfCoefficients > 1 && !leave) {
        std::unique_lock<std::mutex> lckCritical(mutexCritical,
                                                 std::defer_lock);
        lckCritical.lock();
        numberOfCoefficients = (ptrToString.size());
        id = -1;
        if (numberOfCoefficients > 3 ||
            (active_threads == 0 && numberOfCoefficients > 2)) {
          active_threads++;
          pthread_mutex_lock(&mPreserveFermat);
          while (idFermats.empty())
            pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
          id = idFermats.top();
          idFermats.pop();
          pthread_mutex_unlock(&mPreserveFermat);

          if (id >= 0) {
            fermat[id].fermat_collect(
                const_cast<char*>((ptrToString.front()).c_str()));
            ptrToString.pop_front();
            fermat[id].fermat_collect(const_cast<char*>("+"));
            fermat[id].fermat_collect(
                const_cast<char*>((ptrToString.front()).c_str()));
            ptrToString.pop_front();
          }
        }
        else {
          leave = 1;
        }

        lckCritical.unlock();

        if (id >= 0) {
          fermat[id].fermat_calc();
          string mth;
          mth.assign(fermat[id].g_baseout);

          pthread_mutex_lock(&mPreserveFermat);
          idFermats.push(id);
          pthread_cond_signal(&cPreserveFermat);
          pthread_mutex_unlock(&mPreserveFermat);

          lckCritical.lock();
          ptrToString.push_front(mth);
          ptrToString.sort(myfunction);
          numberOfCoefficients = (ptrToString.size());
          active_threads--;
          lckCritical.unlock();
        }
        id = -1;
      }
      {
        std::lock_guard<std::mutex> lock(m);
        countcoeff++;
        if (countcoeff == totalCoeff) {
          notified = true;
          cond_var.notify_one();
        }
      }
    });
  }

  {
    std::unique_lock<std::mutex> lock(m);
    cond_var.wait(lock, [&notified] { return notified; });
  }

  pthread_mutex_destroy(&mPreservePool);
  pthread_cond_destroy(&cPreservePool);

  pthread_mutex_lock(&mPreserveFermat);
  while (idFermats.empty())
    pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
  id = idFermats.top();
  idFermats.pop();
  pthread_mutex_unlock(&mPreserveFermat);

  //   string abc, def;

  if (id >= 0) {
    fermat[id].fermat_collect(const_cast<char*>((ptrToString.front()).c_str()));
    //     abc = ptrToString.front();
    ptrToString.pop_front();
    fermat[id].fermat_collect(const_cast<char*>("+"));
    fermat[id].fermat_collect(const_cast<char*>((ptrToString.front()).c_str()));
    //     def = ptrToString.front();
    ptrToString.pop_front();
  }

  if (id >= 0) {
    fermat[id].fermat_calc();
    string mth;
    mth.assign(fermat[id].g_baseout);

    pthread_mutex_lock(&mPreserveFermat);
    idFermats.push(id);
    pthread_cond_signal(&cPreserveFermat);
    pthread_mutex_unlock(&mPreserveFermat);

    ptrToString.push_front(mth);
    ptrToString.sort(myfunction);
    numberOfCoefficients = (ptrToString.size());
    active_threads--;
  }

  return 1;
}

int Kira::calculate_coefficient_term2(list<string> ptrToString,
                                      string& interpVar, int numbers,
                                      string& result, int idComb) {
  pthread_mutex_t mPreservePool;
  pthread_cond_t cPreservePool;

  pthread_mutex_init(&mPreservePool, NULL);
  pthread_cond_init(&cPreservePool, NULL);

  int numberOfCoefficients = ptrToString.size();
  string mth;
  int id = -1;

  char var[100];

  var[0] = '\0';
  snprintf(var, sizeof(var), "%s\n", interpVar.c_str());

  if (numberOfCoefficients == 1) {
    pthread_mutex_lock(&mPreserveFermat);
    while (idFermats.empty())
      pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
    id = idFermats.top();
    idFermats.pop();
    pthread_mutex_unlock(&mPreserveFermat);

    fermat[id].unset_variable(var);
    fermat[id].set_numeric(var, numbers);

    if (id >= 0) {
      fermat[id].fermat_collect(
          const_cast<char*>((ptrToString.front()).c_str()));
      ptrToString.pop_front();

      if (!fermat[id].fermat_calc(1)) {
        fermat[id].unset_numeric(var);
        fermat[id].set_variable(var);

        pthread_mutex_lock(&mPreserveFermat);
        idFermats.push(id);
        pthread_cond_signal(&cPreserveFermat);
        pthread_mutex_unlock(&mPreserveFermat);
        return 0;
      }

      mth.assign(fermat[id].g_baseout);
      ptrToString.push_front(mth);

      fermat[id].unset_numeric(var);
      fermat[id].set_variable(var);

      pthread_mutex_lock(&mPreserveFermat);
      idFermats.push(id);
      pthread_cond_signal(&cPreserveFermat);
      pthread_mutex_unlock(&mPreserveFermat);
    }
    result = ptrToString.front();
    return 1;
  }

  int active_threads = 0;

  //   coreNumber = 11;

  unsigned totalCoeff = ((numberOfCoefficients / 3) < coreNumber)
                            ? (numberOfCoefficients / 3)
                            : coreNumber;

  if (totalCoeff == 0) totalCoeff = 1;

  unsigned countcoeff = 0;

  std::mutex m;
  std::condition_variable cond_var;
  bool notified = false;
  int broken = 0;

  for (unsigned itNum = 0; itNum < totalCoeff; itNum++) {
    pool3[idComb].enqueue([this, &broken, &var, &numbers, &countcoeff,
                           totalCoeff, &cond_var, &notified, &m,
                           &numberOfCoefficients, &active_threads,
                           &ptrToString]() {
      int id = -1;
      int leave = 0;
      while (numberOfCoefficients > 1 && !leave) {
        std::unique_lock<std::mutex> lckCritical(mutexCritical,
                                                 std::defer_lock);
        lckCritical.lock();
        numberOfCoefficients = (ptrToString.size());
        id = -1;
        if (numberOfCoefficients > 3 ||
            (active_threads == 0 && numberOfCoefficients > 2)) {
          active_threads++;
          pthread_mutex_lock(&mPreserveFermat);
          while (idFermats.empty())
            pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
          id = idFermats.top();
          idFermats.pop();
          pthread_mutex_unlock(&mPreserveFermat);

          fermat[id].unset_variable(var);
          fermat[id].set_numeric(var, numbers);
          fermat[id].fermat_collect(
              const_cast<char*>((ptrToString.front()).c_str()));
          ptrToString.pop_front();
          fermat[id].fermat_collect(const_cast<char*>("+"));
          fermat[id].fermat_collect(
              const_cast<char*>((ptrToString.front()).c_str()));
          ptrToString.pop_front();
        }
        else {
          leave = 1;
        }

        lckCritical.unlock();

        if (id >= 0) {
          if (!fermat[id].fermat_calc(1)) {
            leave = 1;
            broken = 1;
          }

          string mth;
          mth.assign(fermat[id].g_baseout);

          fermat[id].unset_numeric(var);
          fermat[id].set_variable(var);

          pthread_mutex_lock(&mPreserveFermat);
          idFermats.push(id);
          pthread_cond_signal(&cPreserveFermat);
          pthread_mutex_unlock(&mPreserveFermat);

          lckCritical.lock();
          ptrToString.push_front(mth);
          ptrToString.sort(myfunction);
          numberOfCoefficients = (ptrToString.size());
          active_threads--;
          lckCritical.unlock();
        }
        id = -1;
      }
      {
        std::lock_guard<std::mutex> lock(m);
        countcoeff++;
        if (countcoeff == totalCoeff) {
          notified = true;
          cond_var.notify_one();
        }
      }
    });
  }

  {
    std::unique_lock<std::mutex> lock(m);
    cond_var.wait(lock, [&notified] { return notified; });
  }

  pthread_mutex_destroy(&mPreservePool);
  pthread_cond_destroy(&cPreservePool);

  if (broken) {
    return 0;
  }

  pthread_mutex_lock(&mPreserveFermat);
  while (idFermats.empty())
    pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
  id = idFermats.top();
  idFermats.pop();
  pthread_mutex_unlock(&mPreserveFermat);

  fermat[id].unset_variable(var);
  fermat[id].set_numeric(var, numbers);

  if (id >= 0) {
    fermat[id].fermat_collect(const_cast<char*>((ptrToString.front()).c_str()));
    //     abc = ptrToString.front();
    ptrToString.pop_front();
    fermat[id].fermat_collect(const_cast<char*>("+"));
    fermat[id].fermat_collect(const_cast<char*>((ptrToString.front()).c_str()));
    //     def = ptrToString.front();
    ptrToString.pop_front();

    if (!fermat[id].fermat_calc(1)) {
      fermat[id].unset_numeric(var);
      fermat[id].set_variable(var);

      pthread_mutex_lock(&mPreserveFermat);
      idFermats.push(id);
      pthread_cond_signal(&cPreserveFermat);
      pthread_mutex_unlock(&mPreserveFermat);
      return 0;
    }

    string mth;
    mth.assign(fermat[id].g_baseout);

    fermat[id].unset_numeric(var);
    fermat[id].set_variable(var);

    pthread_mutex_lock(&mPreserveFermat);
    idFermats.push(id);
    pthread_cond_signal(&cPreserveFermat);
    pthread_mutex_unlock(&mPreserveFermat);

    ptrToString.push_front(mth);
    ptrToString.sort(myfunction);
    numberOfCoefficients = (ptrToString.size());
    active_threads--;
  }

  //   coreNumber = 1;

  result = ptrToString.front();
  return 1;
}

int Kira::calculate_coefficient(list<string>& ptrToString) {
  int id = -1;

  pthread_mutex_lock(&mPreserveFermat);
  while (idFermats.empty())
    pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
  id = idFermats.top();
  idFermats.pop();
  pthread_mutex_unlock(&mPreserveFermat);
  if (id >= 0) {
    for (auto it = ptrToString.begin(); it != ptrToString.end(); ++it) {
      if (it != ptrToString.begin()) {
        fermat[id].fermat_collect(const_cast<char*>("+"));
      }
      fermat[id].fermat_collect(const_cast<char*>((*it).c_str()));
    }
  }
  return id;
}

void feedFermat(Fermat*& fermat, std::stack<int>& idFermats, string& str,
                std::stack<std::tuple<std::string, int> >& ai, int inID) {
  //   if(pyred::parse_coeff<pyred::Coeff_int>(str).hash()!=0){
  //     pthread_mutex_lock(&mPreserveFermat);
  //       ai.push(make_tuple(str,inID));
  //     pthread_mutex_unlock(&mPreserveFermat);
  //   }
  //   else{
  //     pthread_mutex_lock(&mPreserveFermat);
  //       ai.push(make_tuple("0",inID));
  //     pthread_mutex_unlock(&mPreserveFermat);
  //   }

  int id = -1;

  pthread_mutex_lock(&mPreserveFermat);
  while (idFermats.empty())
    pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
  id = idFermats.top();
  idFermats.pop();
  pthread_mutex_unlock(&mPreserveFermat);

  if (id < 0) {
    throw ExceptionInternal("feedFermat(): Could not get Fermat instance");
  }

  fermat[id].fermat_collect(const_cast<char*>((str.c_str())));
  fermat[id].fermat_calc();

  pthread_mutex_lock(&mPreserveFermat);
  if (inID > -1) {
    ai.push(make_tuple(fermat[id].g_baseout, inID));
  }
  idFermats.push(id);
  pthread_cond_signal(&cPreserveFermat);
  pthread_mutex_unlock(&mPreserveFermat);
}

int Kira::insert_equation(BaseIntegral*& EqPtr, BaseIntegral*& insertEqPtr) {
  unsigned srange = EqPtr[0].length;
  unsigned rrange = insertEqPtr[0].length;

  vector<std::tuple<pyred::Weight, uint32_t, uint32_t, int, std::string, int> >
      collectTerms;

  uint32_t counter = 0;
  unsigned it = 1;
  for (unsigned itt = 1; itt < srange; itt++) {
    while (it < rrange && EqPtr[itt].id < insertEqPtr[it].id) {
      collectTerms.push_back(make_tuple(
          insertEqPtr[it].id, insertEqPtr[it].characteristics[SECTOR],
          insertEqPtr[it].characteristics[TOPOLOGY], (insertEqPtr[it].flag2),
          ("-(" + insertEqPtr[it].coefficientString + ")/(" +
           insertEqPtr[0].coefficientString + ")*(" +
           EqPtr[0].coefficientString + ")"),
          1));
      it++;
      counter++;
    }

    if (it < rrange && EqPtr[itt].id == insertEqPtr[it].id) {
      collectTerms.push_back(make_tuple(
          insertEqPtr[it].id, insertEqPtr[it].characteristics[SECTOR],
          insertEqPtr[it].characteristics[TOPOLOGY], insertEqPtr[it].flag2,
          ("-(" + insertEqPtr[it].coefficientString + ")/(" +
           insertEqPtr[0].coefficientString + ")*(" +
           EqPtr[0].coefficientString + ")+" + "(" +
           EqPtr[itt].coefficientString + ")"),
          1));

      it++;
      counter++;
      continue;
    }

    collectTerms.push_back(
        make_tuple(EqPtr[itt].id, EqPtr[itt].characteristics[SECTOR],
                   EqPtr[itt].characteristics[TOPOLOGY], (EqPtr[itt].flag2),
                   EqPtr[itt].coefficientString, 0));
  }

  while (it < rrange) {
    collectTerms.push_back(make_tuple(
        insertEqPtr[it].id, insertEqPtr[it].characteristics[SECTOR],
        insertEqPtr[it].characteristics[TOPOLOGY], (insertEqPtr[it].flag2),
        ("-(" + insertEqPtr[it].coefficientString + ")/(" +
         insertEqPtr[0].coefficientString + ")*(" + EqPtr[0].coefficientString +
         ")"),
        1));
    counter++;
    it++;
  }

  srange = collectTerms.size();

  //   vector<std::tuple <uint64_t,uint32_t,uint32_t,int,std::string,int> >
  //   collectTerms2;
  //
  //   it = 1;
  //   for (unsigned itt = 1; itt < EqPtr[0].length; itt++){
  //
  //     while(it < rrange && EqPtr[itt].id < insertEqPtr[it].id){
  //
  //       collectTerms2.push_back(
  // 	make_tuple(insertEqPtr[it].id,
  // 		   insertEqPtr[it].characteristics[SECTOR],
  // insertEqPtr[it].characteristics[TOPOLOGY], (insertEqPtr[it].flag2), ("-(" +
  // insertEqPtr[it].coefficientString + ")/(" +
  // insertEqPtr[0].coefficientString + ")*(" + EqPtr[0].coefficientString +
  // ")"),1));
  //       it++;
  // //       counter++;
  //     }
  //
  //     if( it < rrange && EqPtr[itt].id == insertEqPtr[it].id){
  //
  //       collectTerms2.push_back(
  // 	make_tuple(insertEqPtr[it].id,
  // 		   insertEqPtr[it].characteristics[SECTOR],
  // insertEqPtr[it].characteristics[TOPOLOGY], 		   insertEqPtr[it].flag2,
  // 		   ("-(" + insertEqPtr[it].coefficientString +
  // ")/("+insertEqPtr[0].coefficientString+")*("+EqPtr[0].coefficientString+")+"
  // + "(" + EqPtr[itt].coefficientString + ")"),1));
  //
  //       it++;
  // //       counter++;
  //       continue;
  //     }
  //
  //     collectTerms2.push_back(
  //       make_tuple(EqPtr[itt].id,
  // 		   EqPtr[itt].characteristics[SECTOR],
  // EqPtr[itt].characteristics[TOPOLOGY], (EqPtr[itt].flag2),
  // EqPtr[itt].coefficientString,0));
  //   }
  //
  //   while(it<rrange){
  //
  //     collectTerms2.push_back(
  // 	make_tuple(insertEqPtr[it].id,
  // 	       insertEqPtr[it].characteristics[SECTOR],
  // insertEqPtr[it].characteristics[TOPOLOGY], 	       (insertEqPtr[it].flag2),
  // 	       ("-(" + insertEqPtr[it].coefficientString + ")/(" +
  // insertEqPtr[0].coefficientString + ")*(" + EqPtr[0].coefficientString +
  // ")"),1));
  // //     counter++;
  //     it++;
  //   }

  {
    uint32_t counterW = 0;
    std::mutex m;
    std::condition_variable cond_var;
    bool processed = false;
    stack<tuple<string, int> > ai;

    if (counter != 0) {
      for (size_t it = 0; it < collectTerms.size(); it++) {
        if (get<5>(collectTerms[it]) == 1) {
          pool2->enqueue([this, &collectTerms, it, &ai, &counterW, &processed,
                          &counter, &cond_var, &m]() {
            feedFermat(fermat, idFermats, get<4>(collectTerms[it]), ai, it);
            {
              std::lock_guard<std::mutex> lock(m);
              counterW++;
              if (counterW == counter) {
                processed = true;
                cond_var.notify_one();
              }
            }
          });
        }
      }
      {
        std::unique_lock<std::mutex> lock(m);
        cond_var.wait(lock, [&processed] { return processed; });
      }
    }

    delete[] EqPtr;

    if (srange != 0) {
      EqPtr = new BaseIntegral[srange];
      unsigned sizeEq = 0;

      while (!ai.empty()) {
        get<4>(collectTerms[get<1>(ai.top())]) = get<0>(ai.top());
        ai.pop();
      }

      for (size_t it = 0; it < collectTerms.size(); it++) {
        if ((get<4>(collectTerms[it]) != "0")) {
          EqPtr[sizeEq].coefficientString = get<4>(collectTerms[it]);
          EqPtr[sizeEq].id = get<0>(collectTerms[it]);
          EqPtr[sizeEq].characteristics[SECTOR] = get<1>(collectTerms[it]);
          EqPtr[sizeEq].characteristics[TOPOLOGY] = get<2>(collectTerms[it]);
          EqPtr[sizeEq].flag2 = get<3>(collectTerms[it]);
          sizeEq++;
        }
      }

      srange = sizeEq;

      for (unsigned it = 0; it < srange; it++) {
        EqPtr[it].length = srange;
      }
    }
  }

  if (srange != 0)
    return 1;
  else
    return 0;
}

void Kira::sort_equation(BaseIntegral*& EqPtr) {

  unsigned srange = EqPtr[0].length;

  for (unsigned i = 0; i < srange - 1; i++) {
    for (unsigned ii = i + 1; ii < srange; ii++) {
      if (EqPtr[ii].id > EqPtr[i].id) {
        swap(EqPtr[i], EqPtr[ii]);
      }
    }
  }

  EqPtr[0].length = srange;

#ifdef KIRA_DEBUG
  if (srange == 0) {
    logger.set_level(2);
    logger << "This equation is zero: \n";
    logger.set_level(1);
  }
  else if (srange < 2) {
    logger.set_level(2);
    logger << "This integral is zero: \n";
    for (int i = 0; i < srange; i++) {
      logger << "Int[" << EqPtr[i].id << "," << EqPtr[i].flag2 << "]"
             << "*";
      logger << EqPtr[i].coefficientString << "\n";
    }
    logger << "--------------------------------"
           << "\n";
    logger.set_level(1);
  }
#endif
}

void Kira::parallel_pool(
    std::unordered_map<pyred::Weight, std::tuple<BaseIntegral***, int> >&
        forwardRed,
    std::vector<pyred::Weight>& reduce_please_back_up, int& red_1,
    std::unordered_map<pyred::Weight, int>& rdy2F,
    std::set<pyred::Weight>& /*setMandatory*/, pyred::Weight it,
    uint32_t& information) {
  BaseIntegral*** fwEQ;
  unsigned lengthEQ;

  pthread_mutex_lock(&mPreserveResults);

  while (((rdy2F[it] >> 0) & 1))
    pthread_cond_wait(&cPreserveResults, &mPreserveResults);

  rdy2F[it] |= 1 << 0;
  auto mapContent = forwardRed.find(it);
  std::tie(fwEQ, lengthEQ) = mapContent->second;

  get<1>(mapContent->second) = 1;
  pthread_mutex_unlock(&mPreserveResults);

  for (unsigned k = 1; k < lengthEQ; k++) {
    int ende54 = 0;
    int grenze = (fwEQ[0][0][0].length < fwEQ[0][k][0].length)
                     ? fwEQ[0][0][0].length
                     : fwEQ[0][k][0].length;

    if (ende54 == 0 && fwEQ[0][0][0].length > fwEQ[0][k][0].length) {
      swap(fwEQ[0][0], fwEQ[0][k]);
      ende54 = 1;
    }
    if (ende54 == 0 && fwEQ[0][0][0].length < fwEQ[0][k][0].length) {
      ende54 = 1;
    }

    if (ende54 != 0) continue;

    for (int i = 1; i < grenze; i++) {
      if ((fwEQ[0][0][i].id > fwEQ[0][k][i].id)) {
        swap(fwEQ[0][0], fwEQ[0][k]);
        ende54 = 1;
        break;
      }
      if ((fwEQ[0][0][i].id < fwEQ[0][k][i].id)) {
        ende54 = 1;
        break;
      }
    }
  }

  for (uint32_t jtt = 1; jtt < lengthEQ; jtt++) {
    insert_equation(fwEQ[0][jtt], fwEQ[0][0]);
  }

  for (uint32_t jtt = 1; jtt < lengthEQ; jtt++) {
    pyred::Weight ID = fwEQ[0][jtt][0].id;
    pthread_mutex_lock(&mPreserveResults);

    auto rdy2FContent = rdy2F.find(ID);

    if (!(rdy2FContent != rdy2F.end())) {
      rdy2F.insert(pair<pyred::Weight, int>(ID, 0));
    }

    while (((rdy2F[ID] >> 0) & 1))
      pthread_cond_wait(&cPreserveResults, &mPreserveResults);
    rdy2F[ID] |= 1 << 0;

    auto mapContent = forwardRed.find(ID);
    if (mapContent != forwardRed.end()) {
      if (get<1>(mapContent->second) >= RED) {
        BaseIntegral*** put2map = new BaseIntegral**[1];
        put2map[0] = new BaseIntegral*[get<1>(mapContent->second) + 1];

        for (int kIt = 0; kIt < get<1>(mapContent->second); kIt++) {
          put2map[0][kIt] = get<0>(mapContent->second)[0][kIt];
        }

        delete[] get<0>(mapContent->second)[0];

        get<0>(mapContent->second)[0] = put2map[0];
      }

      get<0>(mapContent->second)[0][get<1>(mapContent->second)] = fwEQ[0][jtt];
      get<1>(mapContent->second)++;

      if (get<1>(mapContent->second) == 2) {
        reduce_please_back_up.push_back(get<0>(mapContent->second)[0][0]->id);
      }

      red_1++;
    }
    else {
      BaseIntegral*** put2map = new BaseIntegral**[1];
      put2map[0] = new BaseIntegral*[RED];
      put2map[0][0] = fwEQ[0][jtt];
      forwardRed.insert(pair<pyred::Weight, tuple<BaseIntegral***, int> >(
          ID, make_tuple(put2map, 1)));
    }
    pthread_mutex_unlock(&mPreserveResults);

    std::unique_lock<std::mutex> lckCritical(mutexCritical, std::defer_lock);
    // critical section
    lckCritical.lock();
    information++;
    lckCritical.unlock();

    pthread_mutex_lock(&mPreserveResults);
    rdy2F[ID] &= ~(1 << 0);
    pthread_cond_broadcast(&cPreserveResults);
    pthread_mutex_unlock(&mPreserveResults);
  }
  pthread_mutex_lock(&mPreserveResults);
  rdy2F[(it)] &= ~(1 << 0);
  pthread_cond_broadcast(&cPreserveResults);
  pthread_mutex_unlock(&mPreserveResults);
}

void Kira::run_reduction(
    std::unordered_map<pyred::Weight, std::tuple<BaseIntegral***, int> >&
        forwardRed,
    std::vector<pyred::Weight>& reduce_please,
    std::unordered_map<pyred::Weight, int>& rdy2F,
    std::set<pyred::Weight>& setMandatory, tuple<int, int>& printControl) {
  ofstream myfile;
  double finalt = 0;

  uint32_t information = 0;
  int red_1 = -1;
  uint32_t whileloop = 0;

  logger << "\nIterate the system until it is in triangular form: \n";

  std::vector<pyred::Weight> reduce_please_back_up;

  int counter = 0;

  do {
    Clock clock;

    red_1 = 0;

    logger << "iteration: " << whileloop << "; ";
    logger << "items left: " << reduce_please.size() << "; ";

    uint32_t countLoop = 0;
    uint32_t reduceSIZE = reduce_please.size();
    std::mutex m;
    std::condition_variable cond_var;
    bool processed = false;

    if (reduceSIZE != 0) {
      for (auto it = reduce_please.begin(); it < reduce_please.end(); it++) {
        pool->enqueue([this, it, &forwardRed, &reduce_please_back_up, &red_1,
                       &rdy2F, &setMandatory, &countLoop, &information,
                       &processed, &reduceSIZE, &cond_var, &m]() {
          parallel_pool(forwardRed, reduce_please_back_up, red_1, rdy2F,
                        setMandatory, *it, information);

          {
            std::lock_guard<std::mutex> lock(m);
            countLoop++;
            if (countLoop == reduceSIZE) {
              processed = true;
              cond_var.notify_one();
            }
          }
        });
      }

      {
        std::unique_lock<std::mutex> lock(m);
        cond_var.wait(lock, [&processed] { return processed; });
      }
    }


   

    reduce_please.clear();
    reduce_please = reduce_please_back_up;
    reduce_please_back_up.clear();
    
/*    
    {
    reduce_please.clear();
    for (auto toIt = forwardRed.begin(); toIt != forwardRed.end(); toIt++) {
      if (get<1>(toIt->second) > 0) {
        reduce_please.push_back(get<0>(toIt->second)[0][0]->id);
      }
    }
    auto endZero = unique(reduce_please.begin(), reduce_please.end());
    reduce_please.resize(distance(reduce_please.begin(), endZero));
    sort(reduce_please.begin(), reduce_please.end());
    
    for (auto it = reduce_please.begin(); it != reduce_please.end(); it++) {
      BaseIntegral*** fwEQ;
      unsigned lengthEq;
      auto mapContent = forwardRed.find(*it);
      tie(fwEQ, lengthEq) = mapContent->second;
      if(lengthEq==1)
        continue;
      size_t bigCoeff=0;
      for (uint32_t jtt = 0; jtt < lengthEq; jtt++) {
          for (uint32_t j = 0; j < fwEQ[0][jtt][0].length; j++) {
            pyred::Integral(fwEQ[0][jtt][j].id).to_string();
            fwEQ[0][jtt][j].coefficientString;
            if(bigCoeff<fwEQ[0][jtt][j].coefficientString.size())
              bigCoeff=fwEQ[0][jtt][j].coefficientString.size();
          }
        }
        if(bigCoeff<1000)
          reduce_please_back_up.push_back(*it);
      }
      reduce_please.clear();
      reduce_please = reduce_please_back_up;
      reduce_please_back_up.clear();
    }
    
    
//    if(whileloop==18){
    if(red_1==0){
    reduce_please.clear();
    for (auto toIt = forwardRed.begin(); toIt != forwardRed.end(); toIt++) {
      if (get<1>(toIt->second) > 0) {
        reduce_please.push_back(get<0>(toIt->second)[0][0]->id);
      }
    }
    auto endZero = unique(reduce_please.begin(), reduce_please.end());
    reduce_please.resize(distance(reduce_please.begin(), endZero));
    sort(reduce_please.begin(), reduce_please.end());

    
    
    size_t countEq = 0;
    ofstream myfile("johann.fw/sys.johann");
    for (auto it = reduce_please.begin(); it != reduce_please.end(); it++) {
//    for (auto toIt = forwardRed.rbegin(); toIt != forwardRed.rend(); toIt++) {
      BaseIntegral*** fwEQ;
      unsigned lengthEq;
      auto mapContent = forwardRed.find(*it);
      tie(fwEQ, lengthEq) = mapContent->second;
//      tie(fwEQ, lengthEq) = toIt->second;
  
      countEq+=lengthEq;
      
      for (uint32_t jtt = 0; jtt < lengthEq; jtt++) {
        for (uint32_t j = 0; j < fwEQ[0][jtt][0].length; j++) {
//          cout << "bug1: "<< countEq << " " << jtt << " " << j << endl;
//          if(pyred::Integral::is_zero(fwEQ[0][jtt][j].id))
//            continue;
//          cout << "bug2" << endl;
//          cout << fwEQ[0][jtt][j].id << endl;
          myfile << pyred::Integral(fwEQ[0][jtt][j].id).to_string() <<"*(";
//          cout << endl << "bug3" << endl;
          myfile << fwEQ[0][jtt][j].coefficientString <<")" <<endl;
//          cout << "bug4" << endl;
        }
        myfile << endl;
      }
    }
    myfile.close();
    cout << "number of equations after one iteration: " << countEq << endl;
    
    exit(-1);
    }
    */

    if (red_1 == 0) {
      for (auto toIt = forwardRed.begin(); toIt != forwardRed.end();
           toIt++) { // collect all keys to iterate
        if (get<1>(toIt->second) > 0) {
          reduce_please.push_back(get<0>(toIt->second)[0][0]->id);
        }
      }
    }
    finalt += clock.eval_time();

    whileloop++;

    logger << "time: ( " << finalt << " s )\n";

  } while (red_1);

  logger << "redundant: " << counter << "; insertions " << information << "\n";
  // Triangular form is complete, norm the coefficients and write them to the
  // hard disk.
  Clock clock;

  std::atomic<uint32_t> countbad(0);
  std::atomic<uint32_t> countgood(0);
  logger << "All coefficients will be normalized\n";

  int flagSetMandatory = 1;
  if (setMandatory.size()) flagSetMandatory = 0;

  sort(reduce_please.begin(), reduce_please.end());

  for (auto it = reduce_please.rbegin(); it != reduce_please.rend(); it++) {
    BaseIntegral*** fwEQ;
    unsigned lengthEq;
    auto mapContent = forwardRed.find(*it);
    tie(fwEQ, lengthEq) = mapContent->second;

    auto setContent = setMandatory.find(fwEQ[0][0][0].id);

    if ((lengthEq /*&& fwEQ[0][0][0].flag2 == 0*/) &&
        (setContent != setMandatory.end() || flagSetMandatory)) {
      countgood++;

      for (unsigned j = 1; j < fwEQ[0][0][0].length; j++) {
        setMandatory.insert(fwEQ[0][0][j].id);
      }
    }
    else if (lengthEq && (!((get<0>(printControl) != -1) &&
                            (fwEQ[0][0][0].characteristics[SECTOR] !=
                                 get<0>(printControl) ||
                             fwEQ[0][0][0].characteristics[TOPOLOGY] !=
                                 get<1>(printControl))))) {
      countbad++;

      for (unsigned j = 0; j < fwEQ[0][0][0].length; j++) {
        fwEQ[0][0][j].coefficientString.clear();
      }
      delete[] fwEQ[0][0];
      fwEQ[0][0] = 0;
      delete[] fwEQ[0];
      fwEQ[0] = 0;
      delete[] fwEQ;
      get<1>(mapContent->second) = 0;
      forwardRed.erase(mapContent);
    }
    else {
      //       logger << fwEQ[0][0][0].characteristics[SECTOR] << " : " <<
      //       fwEQ[0][0][0].characteristics[TOPOLOGY] << "\n";
    }
  }

  uint32_t triangSysSize = forwardRed.size();
  uint32_t countTriang = 0;
  std::mutex m;
  std::condition_variable cond_var;
  bool processed = false;
  if (triangSysSize != 0) {
    for (auto it = forwardRed.begin(); it != forwardRed.end(); it++) {
      pool->enqueue([this, it, &triangSysSize, &countTriang, &processed,
                     &cond_var, &m]() {
        BaseIntegral*** fwEQ;
        unsigned lengthEq;
        tie(fwEQ, lengthEq) = it->second;

        int idFermat = -1;

        pthread_mutex_lock(&mPreserveFermat);
        while (idFermats.empty())
          pthread_cond_wait(&cPreserveFermat, &mPreserveFermat);
        idFermat = idFermats.top();
        idFermats.pop();
        pthread_mutex_unlock(&mPreserveFermat);

        for (unsigned j = 1; j < fwEQ[0][0][0].length; j++) {
          fermat[idFermat].fermat_collect(const_cast<char*>("-("));
          fermat[idFermat].fermat_collect(
              const_cast<char*>(fwEQ[0][0][j].coefficientString.c_str()));
          fermat[idFermat].fermat_collect(const_cast<char*>(")/("));
          fermat[idFermat].fermat_collect(
              const_cast<char*>(fwEQ[0][0][0].coefficientString.c_str()));
          fermat[idFermat].fermat_collect(const_cast<char*>(")"));
          fermat[idFermat].fermat_calc();
          fwEQ[0][0][j].coefficientString.assign(fermat[idFermat].g_baseout);
        }
        fwEQ[0][0][0].coefficientString = "1";
        pthread_mutex_lock(&mPreserveFermat);
        idFermats.push(idFermat);
        pthread_cond_signal(&cPreserveFermat);
        pthread_mutex_unlock(&mPreserveFermat);
        {
          std::lock_guard<std::mutex> lock(m);
          countTriang++;
          if (countTriang == triangSysSize) {
            processed = true;
            cond_var.notify_one();
          }
        }
      });
    }

    {
      std::unique_lock<std::mutex> lock(m);
      cond_var.wait(lock, [&processed] { return processed; });
    }
  }

  finalt += clock.eval_time();
  logger << "\nEquations out of bounds rmax and smax: " << countbad << "\n";
  logger << "Equations in the bounds rmax and smax: " << countgood << "\n";
  logger << "( " << finalt << " s )"
         << "\n";
}

void Kira::write_result(BaseIntegral*& EqPtr, int i_Integrals,
                        std::ofstream& output) {
  output << "- Eq: " << std::endl;

  for (int i = 0; i < i_Integrals; i++) {
    //     auto mapContent = mastersReMap.find(EqPtr[i].id);
    //
    //     if(mapContent != mastersReMap.end()){
    //
    //       output << "  - ["<<mapContent->second<<",";
    //     }
    //     else{

    //       output << "  - ["<<EqPtr[i].id<<",";
    //     }

    output << "  - [" << EqPtr[i].id << ",";
    output << EqPtr[i].flag2 << ",";
    output << EqPtr[i].characteristics[SECTOR] << ",";
    output << EqPtr[i].characteristics[TOPOLOGY] << ",";
    output << EqPtr[0].length << ",";
    output << "\"";
    output << EqPtr[i].coefficientString;
    output << "\"]" << std::endl;
  }
}

int Kira::get_number_of_masterintegrals(unsigned it,
                                        vector<pyred::Weight>& lriMI,
                                        int** IDarray) {
  vector<int> subsInt;
  int lritID = 0;
  pthread_mutex_lock(&mPreserveResults);

  for (unsigned ktt = 1; ktt < allEq[it][0][0].length; ktt++) {
    pyred::Weight ID = allEq[it][0][ktt].id;

    auto itID = reverseLastReduce.find(ID);

    if (itID != reverseLastReduce.end()) {
      unsigned eqnum = itID->second;

      while (((rdy2P[eqnum] >> 0) & 1)) {
        pthread_cond_wait(&cPreserveResults, &mPreserveResults);
      }
    }
  }
  pthread_mutex_unlock(&mPreserveResults);

  for (unsigned ktt = 1; ktt < allEq[it][0][0].length; ktt++) {
    pyred::Weight ID = allEq[it][0][ktt].id;
    auto itID = reverseLastReduce.find(ID);
    if (itID != reverseLastReduce.end()) {
      unsigned eqnum = itID->second;
      lritID++; // count non master integrals
      IDarray[ktt] = new int[allEq[eqnum][0][0].length];
      IDarray[ktt][0] = allEq[eqnum][0][0].length;
      for (int i = 1; i < IDarray[ktt][0]; i++) {
        IDarray[ktt][i] = i;
        lriMI.push_back(allEq[eqnum][0][i].id);
      }
    }
    else {
      IDarray[ktt] = new int[1];
      IDarray[ktt][0] = 0;
      subsInt.push_back(ktt); // book keeping of masterintegrals
      lriMI.push_back(ID);
    }
  }

  lritID = allEq[it][0][0].length -
           lritID; // Number of masterintegrals in the first equation
  IDarray[0] = new int[lritID];
  IDarray[0][0] = lritID;
  for (int i = 1; i < IDarray[0][0]; i++) {
    IDarray[0][i] = subsInt[i - 1];
  }
  sort(lriMI.begin(), lriMI.end());
  lriMI.resize(distance(lriMI.begin(), unique(lriMI.begin(), lriMI.end())));

  return 1;
}

unsigned Kira::addup_common_masters(unsigned it, unsigned j, int** IDarray,
                                    BaseIntegral newEqPtr[],
                                    vector<list<string> >& listPtrToString) {
  // Add up all master integrals, which exist in [it] first.
  if (IDarray[0][0] > 1) {
    for (int jt = IDarray[0][0] - 1; jt > 0; jt--) {
      int JD = IDarray[0][jt];

      list<string> ptrToString;
      string tokenStr = "(";
      tokenStr += allEq[it][0][JD].coefficientString + ")";
      ptrToString.push_back(tokenStr);

      swap(IDarray[0][jt], IDarray[0][(IDarray[0][0] - 1)]);
      IDarray[0][0]--;

      // Check the id of the masterintegral[JD] in eq[it] is the same as in
      // eq[ID2]. If yes then add up this masterintegral which exists in all
      // eq[ID2]
      for (unsigned it2 = 1; it2 < allEq[it][0][0].length; it2++) {
        if (IDarray[it2][0]) {
          pyred::Weight ID2 = allEq[it][0][it2].id;
          auto itID = reverseLastReduce.find(ID2);
          unsigned eqnum = itID->second;
          // eq which will be substituted: (master integrals)
          for (int jt2 = IDarray[it2][0] - 1; jt2 > 0; jt2--) {
            int JD2 = IDarray[it2][jt2];
            if (allEq[eqnum][0][JD2].id == allEq[it][0][JD].id) {
              tokenStr = "(";
              tokenStr += allEq[eqnum][0][JD2].coefficientString + ")*(" +
                          allEq[it][0][it2].coefficientString + ")";
              ptrToString.push_back(tokenStr);
              swap(IDarray[it2][jt2], IDarray[it2][(IDarray[it2][0] - 1)]);
              IDarray[it2][0]--;
            }
          }
        }
      }
      listPtrToString.push_back(ptrToString);
      newEqPtr[j] = allEq[it][0][JD];
      j++;
    }
  }
  return j;
}

unsigned Kira::addup_final_masters(unsigned it, unsigned j, int** IDarray,
                                   BaseIntegral newEqPtr[],
                                   vector<list<string> >& listPtrToString) {
  // Add up coefficients of the integrals with the same id: JD
  for (unsigned itK = 1; itK < allEq[it][0][0].length; itK++) {
    if (IDarray[itK][0] > 1) {
      pyred::Weight ID = allEq[it][0][itK].id;
      auto itID = reverseLastReduce.find(ID);
      unsigned eqnum = itID->second;

      for (int jt = IDarray[itK][0] - 1; jt > 0; jt--) {
        int JD = IDarray[itK][jt];

        list<string> ptrToString;
        string tokenStr = "(";
        tokenStr += allEq[eqnum][0][JD].coefficientString + ")*(" +
                    allEq[it][0][itK].coefficientString + ")";
        ptrToString.push_back(tokenStr);

        swap(IDarray[itK][jt], IDarray[itK][(IDarray[itK][0] - 1)]);
        IDarray[itK][0]--;
        for (unsigned it2 = itK + 1; it2 < allEq[it][0][0].length; it2++) {
          if (IDarray[it2][0]) {
            pyred::Weight ID2 = allEq[it][0][it2].id;
            auto itID2 = reverseLastReduce.find(ID2);
            unsigned eqnum2 = itID2->second;
            for (int jt2 = IDarray[it2][0] - 1; jt2 > 0; jt2--) {
              int JD2 = IDarray[it2][jt2];
              if (allEq[eqnum2][0][JD2].id == allEq[eqnum][0][JD].id) {
                tokenStr = "(";
                tokenStr += allEq[eqnum2][0][JD2].coefficientString + ")*(" +
                            allEq[it][0][it2].coefficientString + ")";
                ptrToString.push_back(tokenStr);

                swap(IDarray[it2][jt2], IDarray[it2][(IDarray[it2][0] - 1)]);
                IDarray[it2][0]--;
              }
            }
          }
        }
        listPtrToString.push_back(ptrToString);
        newEqPtr[j] = allEq[eqnum][0][JD];
        j++;
      }
    }
  }
  return j;
}

double timeTotal = 0;

void Kira::back_subs(std::string& kiraSaveName) {
  Clock clock;
  logger << "\nBack substitution: \n";
  ofstream myfile;

  //   int iterationCount=0;
  int reducedIntegralCount = 0;

  if (dataFile == 1) {
    string nameResult = outputDir + "/results/" + integralfamily.name + "/" + kiraSaveName;
    myfile.open(nameResult.c_str());
  }
  std::set<pyred::Weight> masterIntegral;

  DataBase database(outputDir + "/results/" + kiraSaveName + ".db");

  double lastTime = 300;

  vector<pyred::Weight> allResults;
  vector<pyred::Weight> finalResults;

  unsigned countReduct = 0;

  std::mutex m;
  std::condition_variable cond_var;
  bool processed = false;

  unsigned biggestSize =
      (reduct2StartHere.size() < static_cast<unsigned>(coreNumber))
          ? (reduct2StartHere.size())
          : coreNumber;

  biggestSize = reduct2StartHere.size();

  poolBS = new ThreadPool;
  poolBS->initialize(biggestSize);

  pool4 = new ThreadPool[biggestSize];
  for (size_t i = 0; i < biggestSize; i++)
    pool4[i].initialize(1);

  for (unsigned itReduct = 0; itReduct < biggestSize; itReduct++) {
    unsigned itStart = /*itReduct*/ reduct2StartHere[itReduct];
    unsigned itEnd;
    if (itReduct < reduct2StartHere.size() - 1)
      itEnd = reduct2StartHere[itReduct + 1];
    else
      itEnd = reihen;

    poolBS->enqueue([this, biggestSize, itReduct, &cond_var, &processed, &m,
                     &countReduct, itStart, itEnd, &allResults, &finalResults,
                     &masterIntegral, &lastTime, &reducedIntegralCount, &myfile,
                     &clock, &database]() {
      std::mutex mPool4;
      std::condition_variable cond_varPool4;
      bool processedPool4 = false;
      unsigned countPool4 = 0;
      for (unsigned it = itStart; it < /*reihen*/ itEnd;
           it++ /*=biggestSize*/) {
        pool4[itReduct /*%coreNumber*/].enqueue(
            [this, &cond_varPool4, &processedPool4, &mPool4, &countPool4,
             itStart, itEnd, &allResults, &finalResults, &masterIntegral,
             &lastTime, &reducedIntegralCount, &myfile, &clock, &database, it]()

            {
              int pass2R = 0;

              pthread_mutex_lock(&mPreserveResults);
              if (((rdy2P[it] >> 1) & 1)) {
                pass2R = 1;
                rdy2P[it] &= ~(1 << 1);
              }
              pthread_mutex_unlock(&mPreserveResults);

              if (pass2R) { // continue is not allowed not in a loop

                pthread_mutex_lock(&mPreserveResults);
                if (((rdy2P[it] >> 2) & 1)) {
                  pass2R = 0;
                }
                pthread_mutex_unlock(&mPreserveResults);

                int** IDarray;
                BaseIntegral* newEqPtr;
                unsigned j = 0;

                if (pass2R == 1) {
                  IDarray = new int*[allEq[it][0][0].length];

                  vector<pyred::Weight> lriMI;
                  get_number_of_masterintegrals(it, lriMI, IDarray);
                  newEqPtr = new BaseIntegral[lriMI.size() + 1];

                  newEqPtr[j] = allEq[it][0][0];
                  j++;

                  vector<list<string> > listPtrToString;

                  j = addup_common_masters(it, j, IDarray, newEqPtr,
                                           listPtrToString);
                  j = addup_final_masters(it, j, IDarray, newEqPtr,
                                          listPtrToString);

                  uint32_t numberOfCoefficients = (listPtrToString.size());

                  if (numberOfCoefficients != 0) {
                    uint32_t countNumberofCoeff = 0;

                    std::mutex m2;
                    std::condition_variable cond_var2;
                    bool ready = false;

                    for (int jIT = (numberOfCoefficients - 1); jIT > -1;
                         jIT--) {
                      pool2->enqueue([this, &numberOfCoefficients, &cond_var2,
                                      &ready, &m2, &newEqPtr,
                                      &countNumberofCoeff, &listPtrToString,
                                      jIT]() {
                        listPtrToString[jIT].sort(myfunction);

                        if (listPtrToString[jIT].back().size() >
                            heuristic ) {

                          if (reconstFlag &&
                              ((listPtrToString[jIT].back().size() >
                                    reconstFlag &&
                                listPtrToString[jIT].size() > termNumber) ||
                               (listPtrToString[jIT].back().size() >
                                    reconstFlag * 4 &&
                                listPtrToString[jIT].size() > termNumber / 2) ||
                               (listPtrToString[jIT].back().size() >
                                reconstFlag * 15))) {

                            string abc = build_string(listPtrToString[jIT]);
                            string def = "0";



                            Clock clock333;
                            abc = build_string(listPtrToString[jIT]);
                            def = "0";
                            vector<pair<string, int> > numericVar22;
                            Algebra reconst22(this, interpVar, reductVar,
                                              numericVar22, 2);

                            newEqPtr[jIT + 1].coefficientString =
                                reconst22.reconstruct_final(abc, def);

                          }
                          else {
                            int idComb = -1;
                            pthread_mutex_lock(&mPreserveCombine);
                            while (idCombine.empty())
                              pthread_cond_wait(&cPreserveCombine,
                                                &mPreserveCombine);
                            idComb = idCombine.top();
                            idCombine.pop();
                            pthread_mutex_unlock(&mPreserveCombine);

                            calculate_coefficient_term(listPtrToString[jIT],
                                                       idComb);

                            pthread_mutex_lock(&mPreserveCombine);
                            idCombine.push(idComb);
                            pthread_cond_signal(&cPreserveCombine);
                            pthread_mutex_unlock(&mPreserveCombine);

                            newEqPtr[jIT + 1].coefficientString =
                                listPtrToString[jIT].front();
                          }
                        }
                        else {
                          int id = -1;

                          id = calculate_coefficient(listPtrToString[jIT]);

                          if (id < 0) {
                            throw ExceptionInternal("Kira::back_subs(): Could not get Fermat instance");
                          }

                          fermat[id].fermat_calc();

                          newEqPtr[jIT + 1].coefficientString.assign(
                              fermat[id].g_baseout);

                          pthread_mutex_lock(&mPreserveFermat);
                          idFermats.push(id);
                          pthread_cond_signal(&cPreserveFermat);
                          pthread_mutex_unlock(&mPreserveFermat);
                        }

                        {
                          std::lock_guard<std::mutex> lock(m2);

                          countNumberofCoeff++;

                          if (countNumberofCoeff == numberOfCoefficients) {
                            ready = true;

                            cond_var2.notify_one();
                          }
                        }
                      });
                    }

                    {
                      std::unique_lock<std::mutex> lock(m2);
                      cond_var2.wait(lock, [&ready] { return ready; });
                    }

                    for (unsigned jIT = 1; jIT < j; jIT++) {
                      if (newEqPtr[jIT].coefficientString == "0") {
                        swap(newEqPtr[jIT], newEqPtr[j - 1]);
                        j--;
                        jIT--;
                      }
                    }
                  }
                }
                pthread_mutex_lock(&mPreserveResults);

                for (unsigned ktt = 0; ktt < allEq[it][0][0].length; ktt++) {
                  pyred::Weight ID = allEq[it][0][ktt].id;

                  auto ocContent = occurrence.find(ID);

                  if (ocContent != occurrence.end()) (ocContent->second)--;

                  auto itID = reverseLastReduce.find(ID);

                  if ((ocContent->second) == 0 &&
                      itID != reverseLastReduce.end()) {
                    finalResults.push_back(ID);
                    reducedIntegralCount++;
                  }
                }

                if (pass2R == 1) {
                  for (unsigned ktt = 0; ktt < allEq[it][0][0].length; ktt++)
                    delete[] IDarray[ktt];

                  delete[] IDarray;
                  newEqPtr[0].length = j;
                  swap(newEqPtr, allEq[it][0]);
                  delete[] newEqPtr;
                  newEqPtr = 0;
                  sort_equation(allEq[it][0]);
                  rdy2P[it] &= ~(1 << 0);
                }

                allResults.push_back(allEq[it][0][0].id);

                double nowTime = clock.eval_time();

                if (lastTime < (nowTime - 300)) {
                  lastTime = nowTime;

                  if (allResults.size() > 0) {
                    database.begin_transaction();
                    database.prepare_backsubstitution();

                    for (size_t itfR = 0, end1 = allResults.size();
                          itfR < end1; itfR++) {
                      pyred::Weight ID = allResults[itfR];
                      auto itID = reverseLastReduce.find(ID);
                      if (itID == reverseLastReduce.end()) continue;

                      unsigned eqnum = itID->second;

                      if (length[eqnum] == 0) continue;

                      for (unsigned i = 0; i < allEq[eqnum][0][0].length;
                            i++) {
                        if (i > 0) {
                          masterIntegral.insert(allEq[eqnum][0][i].id);
                        }

                        database.bind_equation(allEq[eqnum][0][i],
                                                allEq[eqnum][0][0].length, ID,
                                                mastersReMap);
                      }
                      if (dataFile == 1)
                        write_result(allEq[eqnum][0],
                                      allEq[eqnum][0][0].length, myfile);
                    }

                    database.commit_transaction();

                    for (size_t itfR = 0, end1 = finalResults.size();
                          itfR < end1; itfR++) {
                      pyred::Weight ID = finalResults[itfR];
                      auto itID = reverseLastReduce.find(ID);
                      if (itID == reverseLastReduce.end()) continue;

                      unsigned eqnum = itID->second;
                      delete[] allEq[eqnum][0];
                      allEq[eqnum][0] = 0;
                      delete[] allEq[eqnum];
                      allEq[eqnum] = 0;
                      length[eqnum] = 0;
                    }

                    finalResults.clear();

                    allResults.clear();

                    logger
                        << "\n\nEquations reduced: " << reducedIntegralCount
                        << " / " << reihen << "\n";
                    logger << "Kira committed results ";
                    logger << "( " << clock.eval_time() << " s )\n";
                  }
                }

                pthread_cond_broadcast(&cPreserveResults);
                pthread_mutex_unlock(&mPreserveResults);

                std::unique_lock<std::mutex> lckCritical(mutexCritical,
                                                         std::defer_lock);
                // critical section
                lckCritical.lock();
                load_bar(reducedIntegralCount /*++iterationCount*/, reihen, 1,
                         1000);
                lckCritical.unlock();
              }
              {
                std::lock_guard<std::mutex> lockPool4(mPool4);
                countPool4++;

                if ((countPool4 /**biggestSize*/ + itStart) /*>=*/ ==
                    itEnd /*reihen*/) {
                  processedPool4 = true;
                  cond_varPool4.notify_one();
                }
              }
            });
      }
      {
        std::unique_lock<std::mutex> lockPool4(mPool4);
        cond_varPool4.wait(lockPool4,
                           [&processedPool4] { return processedPool4; });
      }

      pthread_mutex_lock(&mPreserveResults);

      if (allResults.size() > 0) {
        database.begin_transaction();
        database.prepare_backsubstitution();

        for (size_t itfR = 0, end1 = allResults.size(); itfR < end1; itfR++) {
          pyred::Weight ID = allResults[itfR];
          auto itID = reverseLastReduce.find(ID);
          if (itID == reverseLastReduce.end()) continue;

          unsigned eqnum = itID->second;

          if (length[eqnum] == 0) continue;

          for (unsigned i = 0; i < allEq[eqnum][0][0].length; i++) {
            if (i > 0) {
              masterIntegral.insert(allEq[eqnum][0][i].id);
            }

            database.bind_equation(allEq[eqnum][0][i],
                                   allEq[eqnum][0][0].length, ID, mastersReMap);
          }
          if (dataFile == 1)
            write_result(allEq[eqnum][0], allEq[eqnum][0][0].length, myfile);
        }

        database.commit_transaction();

        for (size_t itfR = 0, end1 = finalResults.size(); itfR < end1; itfR++) {
          pyred::Weight ID = finalResults[itfR];
          auto itID = reverseLastReduce.find(ID);
          if (itID == reverseLastReduce.end()) continue;

          unsigned eqnum = itID->second;
          delete[] allEq[eqnum][0];
          allEq[eqnum][0] = 0;
          delete[] allEq[eqnum];
          allEq[eqnum] = 0;
          length[eqnum] = 0;
        }

        finalResults.clear();

        allResults.clear();

        logger << "\n\nEquations reduced: " << reducedIntegralCount << " / "
               << reihen << "\n";
        logger << "Kira committed results ";
        logger << "( " << clock.eval_time() << " s )\n";
      }
      pthread_cond_broadcast(&cPreserveResults);
      pthread_mutex_unlock(&mPreserveResults);
      {
        std::lock_guard<std::mutex> lock(m);
        countReduct++;
        if (countReduct == biggestSize) {
          processed = true;
          cond_var.notify_one();
        }
      }
    });
  }
  {
    std::unique_lock<std::mutex> lock(m);
    cond_var.wait(lock, [&processed] { return processed; });
  }

  delete poolBS;
  delete[] pool4;

  if (dataFile == 1) myfile.close();

  logger << "\n" << reducedIntegralCount << " equations were reduced\n";

  vector<pyred::Weight> masterVector;
  string masterFileName =
      outputDir + "/results/" + integralfamily.name + "/masters";
  for (const auto &igl: pyred::TopoConfigData::import_integrals(masterFileName))
  {
    masterVector.push_back(igl.to_weight());
  }

  string pathMI =
      outputDir + "/results/" + integralfamily.name + "/masters.final";
  myfile.open(pathMI.c_str());

  logger << "\n*****Master integrals***********************************\n";

  uint32_t countMaster = 0;
  uint32_t countEmptyMaster = 0;

  std::vector<pyred::Weight> iid_vec;
  for (auto itg : masterIntegral) {
    auto itM = find(masterVector.begin(), masterVector.end(), itg);
    if (itM != masterVector.end()) {
      iid_vec.push_back(itg);
      countMaster++;
    }
  }
  for (auto itg : masterVector) {
    auto itM = find(masterIntegral.begin(), masterIntegral.end(), itg);
    if (itM == masterIntegral.end()) {
      iid_vec.push_back(itg);
      countMaster++;
      countEmptyMaster++;
    }
  }
  auto masters_str = masters_list_str(iid_vec);
  logger << masters_str;
  myfile << masters_str;

  if (countEmptyMaster != 0) {
    logger << "\n*****Master integrals which do not appear on the right**\n";
    logger << "*****hand side of the selected integrals****************\n";

    iid_vec.clear();
    for (auto itg : masterVector) {
      auto itM = find(masterIntegral.begin(), masterIntegral.end(), itg);

      if (itM == masterIntegral.end()) {
        iid_vec.push_back(itg);
      }
    }
    logger << masters_list_str(iid_vec);
  }

  logger << "\nTotal number of master integrals: " << countMaster << "\n";

  uint32_t countNonMaster = 0;

  for (auto itg : masterIntegral) {
    auto itM = find(masterVector.begin(), masterVector.end(), itg);

    if (itM == masterVector.end()) {
      countNonMaster++;
    }
  }

  if (countNonMaster != 0) {
    logger << "\n*****Not selected integrals******\n";
    iid_vec.clear();
    for (auto itg : masterIntegral) {
      auto itM = find(masterVector.begin(), masterVector.end(), itg);

      if (itM == masterVector.end()) {
        iid_vec.push_back(itg);
      }
    }
    logger << masters_list_str(iid_vec);
    logger << "\nNumber of not selected integrals: "
           << (masterIntegral.size() - countMaster) << "\n";
  }

  logger << "( " << clock.eval_time() << " s )\n";
}
