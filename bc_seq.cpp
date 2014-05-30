
#include "defs.h"
#include "sprng.h"
#include "utils.h"

boost::shared_ptr<DOUBLE_T[]>  seq_get_bc(graph g)
{
    int seed = 2387;
    int tid = 0;
    int nthreads = 1;
    int * sprng_stream;
    const int MAX_NUM_PHASES = 50;

    auto bc = boost::shared_ptr<DOUBLE_T[]>(new DOUBLE_T[g.n]);


    scope_guard __f1([&sprng_stream, tid, nthreads, seed](){
        sprng_stream = init_sprng(0, tid, nthreads, seed, SPRNG_DEFAULT);
    }, [&sprng_stream](){
        free(sprng_stream);
    });

    int numV = 1<<g.k4approx;
    LONG_T start, end;
    auto srcs = boost::shared_ptr<LONG_T[]>(new LONG_T[g.n]);
    auto sig = boost::shared_ptr<double[]>(new double[g.n]);
    auto d = boost::shared_ptr<double[]>(new double[g.n]);
    auto del = boost::shared_ptr<double[]>(new double[g.n]);
    auto S = boost::shared_ptr<VERT_T []>(new VERT_T[g.n]);

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

    std::cout<<"done permute!"<<std::endl;
    int num_traversals = 0;
    for (int p=0; p<g.n; p++) {
        int i = srcs[p];
        if (g.numEdges[i+1] - g.numEdges[i] == 0) {
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

            for(int j = g.numEdges[v]; j<g.numEdges[v+1]; j++)
            {
                if (g.weight[j] & 7 == 0)
                    continue;
                int w = g.endV[j];
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
    }

    return bc;

}



double betweennessCentrality(graph* G, DOUBLE_T* BC) {

    VERT_T *S;         /* stack of vertices in the order of non-decreasing
                          distance from s. Also used to implicitly
                          represent the BFS queue */
    plist* P;          /* predecessors of a vertex v on shortest paths from s */
    DOUBLE_T* sig;     /* No. of shortest paths */
    LONG_T* d;         /* Length of the shortest path between every pair */
    DOUBLE_T* del;     /* dependency of vertices */
    LONG_T *in_degree, *numEdges, *pSums;
    LONG_T *pListMem;
    LONG_T* Srcs;
    LONG_T *start, *end;
    LONG_T MAX_NUM_PHASES;
    LONG_T *psCount;
    int seed = 2387;


    VERT_T *myS, *myS_t;
    LONG_T myS_size;
    LONG_T i, j, k, p, count, myCount;
    LONG_T v, w, vert;
    LONG_T numV, num_traversals, n, m, phase_num;
    LONG_T tid, nthreads;
    int* stream;
    tid = 0;
    nthreads = 1;

    /* numV: no. of vertices to run BFS from = 2^K4approx */
    numV = 1<<G->k4approx;;
    n = G->n;
    m = G->m;

    /* Permute vertices */
    if (tid == 0) {
        Srcs = (LONG_T *) malloc(n*sizeof(LONG_T));
    }


    /* Initialize RNG stream */
    stream = init_sprng(0, tid, nthreads, seed, SPRNG_DEFAULT);

    for (i=0; i<n; i++) {
        Srcs[i] = i;
    }

    for (i=0; i<n; i++) {
        j = n*sprng(stream);
        if (i != j) {
                    k = Srcs[i];
                    Srcs[i] = Srcs[j];
                    Srcs[j] = k;
        }
    }


    /* Start timing code from here */
    if (tid == 0) {
        //elapsed_time = get_seconds();
        MAX_NUM_PHASES = 50;
    }

    /* Initialize predecessor lists */

    /* The size of the predecessor list of each vertex is bounded by
       its in-degree. So we first compute the in-degree of every
       vertex */

    if (tid == 0) {
        P   = (plist  *) calloc(n, sizeof(plist));
        in_degree = (LONG_T *) calloc(n+1, sizeof(LONG_T));
        numEdges = (LONG_T *) malloc((n+1)*sizeof(LONG_T));
        pSums = (LONG_T *) malloc(nthreads*sizeof(LONG_T));
    }

    for (i=0; i<m; i++) {
        v = G->endV[i];
        in_degree[v]++;
    }

    prefix_sums(in_degree, numEdges, pSums, n);

    if (tid == 0) {
        pListMem = (LONG_T *) malloc(m*sizeof(LONG_T));
    }

    for (i=0; i<n; i++) {
        P[i].list = pListMem + numEdges[i];
        P[i].degree = in_degree[i];
        P[i].count = 0;
    }

    /* Allocate shared memory */
    if (tid == 0) {
        free(in_degree);
        free(numEdges);
        free(pSums);

        S   = (VERT_T *) malloc(n*sizeof(VERT_T));
        sig = (DOUBLE_T *) malloc(n*sizeof(DOUBLE_T));
        d   = (LONG_T *) malloc(n*sizeof(LONG_T));
        del = (DOUBLE_T *) calloc(n, sizeof(DOUBLE_T));

        start = (LONG_T *) malloc(MAX_NUM_PHASES*sizeof(LONG_T));
        end = (LONG_T *) malloc(MAX_NUM_PHASES*sizeof(LONG_T));
        psCount = (LONG_T *) malloc((nthreads+1)*sizeof(LONG_T));
    }

    /* local memory for each thread */
    myS_size = (2*n)/nthreads;
    myS = (LONG_T *) malloc(myS_size*sizeof(LONG_T));
    num_traversals = 0;
    myCount = 0;

    for (i=0; i<n; i++) {
        d[i] = -1;
    }

    for (p=0; p<n; p++) {

        i = Srcs[p];
        if (G->numEdges[i+1] - G->numEdges[i] == 0) {
            continue;
        } else {
            num_traversals++;
        }

        if (num_traversals == numV + 1) {
          printf("traversal num : %d\n", num_traversals);
            break;
        }

        if (tid == 0) {
            sig[i] = 1;
            d[i] = 0;
            S[0] = i;
            start[0] = 0;
            end[0] = 1;
        }

        count = 1;
        phase_num = 0;


        while (end[phase_num] - start[phase_num] > 0) {

            myCount = 0;
            for (vert = start[phase_num]; vert < end[phase_num]; vert++) {
                v = S[vert];
                for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {

                        w = G->endV[j];
                        if (v != w) {

                                /* w found for the first time? */
                                if (d[w] == -1) {
                                    if (myS_size == myCount) {
                                        /* Resize myS */
                                        myS_t = (LONG_T *)
                                            malloc(2*myS_size*sizeof(VERT_T));
                                        memcpy(myS_t, myS, myS_size*sizeof(VERT_T));
                                        free(myS);
                                        myS = myS_t;
                                        myS_size = 2*myS_size;
                                    }
                                    myS[myCount++] = w;
                                    d[w] = d[v] + 1;
                                    sig[w] = sig[v];
                                    P[w].list[P[w].count++] = v;
                                } else if (d[w] == d[v] + 1) {
                                    sig[w] += sig[v];
                                    P[w].list[P[w].count++] = v;
                                }

                        }
                }
             }
            /* Merge all local stacks for next iteration */
            phase_num++;

            psCount[tid+1] = myCount;


            if (tid == 0) {
                start[phase_num] = end[phase_num-1];
                psCount[0] = start[phase_num];
                for(k=1; k<=nthreads; k++) {
                    psCount[k] = psCount[k-1] + psCount[k];
                }
                end[phase_num] = psCount[nthreads];
            }


            for (k = psCount[tid]; k < psCount[tid+1]; k++) {
                S[k] = myS[k-psCount[tid]];
            }

            count = end[phase_num];
        }

        phase_num--;


        while (phase_num > 0) {
            for (j=start[phase_num]; j<end[phase_num]; j++) {
                w = S[j];
                for (k = 0; k<P[w].count; k++) {
                    v = P[w].list[k];
                    del[v] = del[v] + sig[v]*(1+del[w])/sig[w];
                }
                BC[w] += del[w];
            }

            phase_num--;

        }


        for (j=0; j<count; j++) {
            w = S[j];
            d[w] = -1;
            del[w] = 0;
            P[w].count = 0;
        }



    }

    free(myS);

    if (tid == 0) {
        free(S);
        free(pListMem);
        free(P);
        free(sig);
        free(d);
        free(del);
        free(start);
        free(end);
        free(psCount);
        //elapsed_time = get_seconds() - elapsed_time;
        free(Srcs);
    }

    free_sprng(stream);
    return 0;//elapsed_time;
}
