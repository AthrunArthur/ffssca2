#include "defs.h"
#include <chrono>


void prefix_sums(LONG_T *input, LONG_T* result, LONG_T* p, LONG_T n) {

    LONG_T i;

#ifdef _OPENMP
    LONG_T j, r, start, end, add_value;
    int tid, nthreads;

    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    r = n/nthreads;

    result[0] = 0;

    #pragma omp for
    for (i=1; i<n+1; i++)
        result[i] = input[i-1];

    start =  tid*r + 1;
    end   = (tid+1)*r;

    if (tid == nthreads-1)
        end = n+1;

    for (j=start; j<end; j++)
        result[j] = input[j-1] + result[j-1];

    p[tid] = result[end-1];

    #pragma omp barrier

    if (tid == 0) {
        for (j=1; j<nthreads; j++)
            p[j] += p[j-1];
    }

    #pragma omp barrier

    if (tid>0) {
        add_value=p[tid-1];
        for (j=start-1; j<end; j++)
            result[j] += add_value;
    }

#else

    result[0] = 0;
    for (i=1; i<n+1; i++) {
        result[i] = result[i-1] + input[i-1];
    }

#endif

}
int count_elapse_microsecond(const std::function<void ()> & f)
{
    using namespace std::chrono;

    time_point<system_clock> start, end;

    start = system_clock::now();
    f();
    end = system_clock::now();
    return duration_cast<microseconds>(end-start).count();
}
