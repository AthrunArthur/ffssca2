
#include "defs.h"
#include "sprng.h"
#include "utils.h"
#include <chrono>
#include <ff.h>
using ff::scope_guard;


boost::shared_ptr<DOUBLE_T[]>  seq_get_bc(graph g)
{
    int seed = 2387;
    int tid = 0;
    int nthreads = 1;
    int * sprng_stream;

    auto bc_m = boost::shared_ptr<DOUBLE_T[]>(new DOUBLE_T[g.n]);
    auto bc = bc_m.get();


    scope_guard __f1([&sprng_stream, tid, nthreads, seed](){
        sprng_stream = init_sprng(0, tid, nthreads, seed, SPRNG_DEFAULT);
    }, [&sprng_stream](){
        free(sprng_stream);
    });

    using namespace std::chrono;
    time_point<system_clock> time_start, time_end;
    time_start = system_clock::now();

    int numV = 1<<g.k4approx;
    LONG_T start, end;

    auto srcs_m = boost::shared_ptr<LONG_T[]>(new LONG_T[g.n]);
    auto sig_m = boost::shared_ptr<double[]>(new double[g.n]);
    auto d_m = boost::shared_ptr<double[]>(new double[g.n]);
    auto del_m = boost::shared_ptr<double[]>(new double[g.n]);
    auto S_m = boost::shared_ptr<VERT_T []>(new VERT_T[g.n]);

    auto srcs = srcs_m.get();
    auto sig = sig_m.get();
    auto d = d_m.get();
    auto del = del_m.get();
    auto S = S_m.get();
    auto g_numEdges = g.numEdges.get();
    auto g_weight = g.weight.get();
    auto g_endV = g.endV.get();

    typedef std::vector<VERT_T> array_list;
    typedef std::shared_ptr<array_list> array_list_ptr;
    typedef std::vector<array_list_ptr> dyn_list;
    dyn_list P;
    for(int i = 0; i < g.n; ++i){
        P.push_back(array_list_ptr(new array_list()));
        d[i] = -1;
        del[i] = 0;
        bc[i] = 0;
        srcs[i] = i;
    }
    //permute srcs
    for(int i = 0; i < g.n; ++i)
    {
        int j = g.n * sprng(sprng_stream);
        if (i != j)
        {
            std::swap(srcs[i], srcs[j]);
        }
    }


////////////////
    int num_traversals = 0;
    for (int p=0; p<g.n; p++) {
        int i = srcs[p];
        if (g_numEdges[i+1] - g_numEdges[i] == 0) {
            continue;
        } else {
            num_traversals++;
        }

        if (num_traversals == numV + 1) {
            break;
        }

        start = 0;
        S[0] = i;
        end = 1;
        sig[i] = 1;
        d[i] = 0;
        while (end - start > 0)
        {
            int v = S[start];

            //std::cout<<"edges start from "<< v<<std::endl;
            for(int j = g_numEdges[v]; j<g_numEdges[v+1]; j++)
            {
                if ((g_weight[j] & 7) == 0)
                    continue;
                int w = g_endV[j];
                //std::cout<<"\t"<<w<<std::endl;
                if(v!= w){
                    if(d[w] < 0){
                        S[end++] = w;
                        d[w] = d[v] + 1;
                        sig[w] = 0;
                    }
                    if(d[w] == d[v] + 1)
                    {
                        sig[w] = sig[w] + sig[v];
                        P[w]->push_back(v);
                    }
                }
            }
            start ++;
        }
        for (LONG_T j=end-1; j>0; j--)
        {
            VERT_T w = S[j];
            array_list & pw = *P[w];
            for(auto v : pw)
            {
                del[v] = del[v] + sig[v] * (1 + del[w])/sig[w];
            }
            bc[w] += del[w];
        }
        for (LONG_T j=end-1; j>=0; j--)
        {
            VERT_T v = S[j];
            d[v] = -1;
            del[v] = 0;
            P[v]->clear();
        }
        /*
        time_point<system_clock> tiny_time4 = system_clock::now();
        std::cout<<"time\t"<<duration_cast<microseconds>(tiny_time2-tiny_time1).count()
                   <<"\t"<<duration_cast<microseconds>(tiny_time3-tiny_time2).count()
                     <<"\t"<<duration_cast<microseconds>(tiny_time4-tiny_time3).count()<<std::endl;
                     */
    }
    time_end = system_clock::now();
    ////////////////////////
    auto elapsed_time = duration_cast<microseconds>(time_end-time_start).count();
    std::cout<<"seq_get_bc elapsed_time : "<<elapsed_time<<std::endl;
    return bc_m;
}
