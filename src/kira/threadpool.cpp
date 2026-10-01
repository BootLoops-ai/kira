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

#include "kira/threadpool.h"

// the constructor just launches some amount of workers
void ThreadPool::initialize(const uint32_t threads) {
  for (uint32_t i = 0; i < threads; ++i) {
    workers.emplace_back(std::thread(Worker(*this)));
  }
}

// the destructor joins all threads
ThreadPool::~ThreadPool() {
  // stop all threads
  {
    std::lock_guard<std::mutex> lock(queue_mutex);
    stop = true;
  }

  condition.notify_all();

  // join them
  for (size_t i = 0; i < workers.size(); ++i) {
    workers[i].join();
  }
}

void Worker::operator()() {
  std::function<void()> task;
  while (true) {
    { // acquire lock
      std::unique_lock<std::mutex> lock(pool.queue_mutex);

      // look for a work item
      while (!pool.stop && pool.tasks.empty()) {
        // if there are none wait for notification
        pool.condition.wait(lock);
      }

      if (pool.stop) {
        // exit if the pool is stopped
        return;
      }

      // get the task from the queue
      task = pool.tasks.front();
      pool.tasks.pop_front();
    } // release lock

    // execute the task
    task();
  }
}
