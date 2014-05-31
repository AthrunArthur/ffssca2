
#include "defs.h"
#include "sprng.h"
#include "utils.h"
#include <chrono>
#include "ff.h"

using namespace ff;



boost::shared_ptr<DOUBLE_T[]>  ff_get_bc(graph g)
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

    ff::para<> f;
    f([](){});

    using namespace std::chrono;
    time_point<system_clock> time_start, time_end;
    time_start = system_clock::now();

    int numV = 1<<g.k4approx;
    typedef boost::shared_ptr<ff::accumulator<DOUBLE_T> > fake_bc_ptr;
    std::vector<fake_bc_ptr> bcs;
    auto srcs_m = boost::shared_ptr<LONG_T[]>(new LONG_T[g.n]);
    auto srcs = srcs_m.get();

    paragroup pinit;
    ff::mutex m;
    pinit.for_each(0, g.n, [&srcs, &bc, &bcs, &m, &g, &sprng_stream](int i){
        srcs[i] = g.n * sprng(sprng_stream);
        bc[i] = 0;
        m.lock();
        bcs.push_back(fake_bc_ptr(new ff::accumulator<DOUBLE_T>(0, [](const DOUBLE_T & x, const DOUBLE_T & y){return x + y;})));
        m.unlock();
    });



    std::vector<VERT_T> to_traverse_verts;
    int num_traversals = 0;
    auto g_numEdges = g.numEdges.get();
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
        to_traverse_verts.push_back(i);
    }
    ff_wait(all(pinit));
    ff::paragroup pg;
    pg.for_each(to_traverse_verts.begin(), to_traverse_verts.end(), [&g, &bcs](VERT_T i){
        thread_local static auto sig_m = boost::shared_ptr<double[]>(new double[g.n]);
        thread_local static auto d_m = boost::shared_ptr<double[]>(new double[g.n]);
        thread_local static auto del_m = boost::shared_ptr<double[]>(new double[g.n]);
        thread_local static auto S_m = boost::shared_ptr<VERT_T []>(new VERT_T[g.n]);
        typedef std::vector<VERT_T> array_list;
        typedef std::shared_ptr<array_list> array_list_ptr;
        typedef std::vector<array_list_ptr> dyn_list;
        thread_local static dyn_list P;

        thread_local static auto sig = sig_m.get();
        thread_local static auto d = d_m.get();
        thread_local static auto del = del_m.get();
        thread_local static auto S = S_m.get();
        thread_local static auto g_numEdges = g.numEdges.get();
        thread_local static auto g_weight = g.weight.get();
        thread_local static auto g_endV = g.endV.get();

        for(int i = 0; i < g.n; ++i){
            P.push_back(array_list_ptr(new array_list()));
            d[i] = -1;
            del[i] = 0;
        }

        LONG_T start = 0, end = 1;
        S[0] = i;
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
            bcs[w]->increase(del[w]);
        }
        for (LONG_T j=end-1; j>=0; j--)
        {
            VERT_T v = S[j];
            d[v] = -1;
            del[v] = 0;
            P[v]->clear();
        }
    });

    ff::ff_wait(all(pg));

    ff::paragroup preduceresults;
    preduceresults.for_each(0, g.n, [&bc, &bcs](int i){
        bc[i] = bcs[i]->get();
    });
    ff::ff_wait(all(preduceresults));
    time_end = system_clock::now();
    ////////////////////////
    auto elapsed_time = duration_cast<microseconds>(time_end-time_start).count();
    std::cout<<"ff_get_bc elapsed_time : "<<elapsed_time<<std::endl;
    return bc_m;
}

