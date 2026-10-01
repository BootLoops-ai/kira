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

#ifndef KIRA_THREADPOOL_H_
#define KIRA_THREADPOOL_H_

#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

class Kira;
class ThreadPool;

// our worker thread objects
class Worker {
public:
  Worker(ThreadPool &s) : pool(s) {}
  void operator()();

private:
  friend class Kira;
  ThreadPool &pool;
};

// the actual thread pool
class ThreadPool {
public:
  ThreadPool() : stop(false){};
  void initialize(const uint32_t);
  template <class F>
  void enqueue(const F &f);
  ~ThreadPool();

private:
  friend class Worker;
  friend class Kira;

  // need to keep track of threads so we can join them
  std::vector<std::thread> workers;

  // the task queue
  std::deque<std::function<void()>> tasks;

  // synchronization
  std::mutex queue_mutex;
  std::condition_variable condition;
  bool stop;
};

// add new work item to the pool
template <class F>
void ThreadPool::enqueue(const F &f) {
  { // acquire lock
    std::lock_guard<std::mutex> lock(queue_mutex);

    // add the task
    tasks.emplace_back(std::function<void()>(f));
  } // release lock

  // wake up one thread
  condition.notify_one();
}

#endif // KIRA_THREADPOOL_H_
