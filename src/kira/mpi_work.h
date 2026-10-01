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

#ifndef KIRA_MPI_WORK_H_
#define KIRA_MPI_WORK_H_

#include "pyred/defs.h"
#include <mpi.h>


class MPIHelper {
  // A singleton to initialise MPI and hold the MPI data type for pyred::Weight.
private:
  int m_provided;
  MPI_Datatype m_weight_type;
  MPIHelper() {
    MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &m_provided);
    // A static_assert() in pyRed guarantees
    // that sizeof(Weight) is a multiple of sizeof(uint64_t).
    int blocklength = sizeof(pyred::Weight)/sizeof(uint64_t);
    int count = 1;
    int blocklengths[1] = {blocklength}; // treat unsigned __int128 as uint64_t[2]
    MPI_Aint displacements[1] = {0};
    MPI_Datatype types[1] = {MPI_UINT64_T};
    MPI_Type_create_struct(count, blocklengths, displacements, types, &m_weight_type);
    MPI_Type_commit(&m_weight_type);
  }
  static MPIHelper & get() {
    static MPIHelper instance;
    return instance;
  }
public:
  MPIHelper(const MPIHelper &) = delete;
  MPIHelper & operator=(const MPIHelper &) = delete;
  static int init() {
    return get().m_provided;
  }
  static MPI_Datatype & weight_type() {
    return get().m_weight_type;
  }
  ~MPIHelper() {
    MPI_Type_free(&weight_type());
    MPI_Finalize();
  }
};

void mpi_work(const uint32_t threads, const uint32_t bunch_size);

#endif // KIRA_MPI_WORK_H_
