#include "ThreadPool.h"
#include <stdexcept>

ThreadPool::ThreadPool(size_t threads) : stop(false), active_tasks(0) {
    for(size_t i = 0; i < threads; ++i) {
        workers.emplace_back([this] {
            while(true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock, [this]{ return this->stop || !this->tasks.empty(); });
                    if(this->stop && this->tasks.empty())
                        return;
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }

                ++active_tasks;
                task();
                --active_tasks;
                completion_condition.notify_all();
            }
        });
    }
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();
}

void ThreadPool::wait() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    completion_condition.wait(lock, [this]() { return tasks.empty() && active_tasks == 0; });
}
