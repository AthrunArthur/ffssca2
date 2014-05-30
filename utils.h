#ifndef UTILS_H
#define UTILS_H
#include <functional>
#include <chrono>



int count_elapse_microsecond(const std::function<void ()> & f);

#endif // UTILS_H
