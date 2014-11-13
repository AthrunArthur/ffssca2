#include "defs.h"
#ifdef CACHE_EVAL
#include <papi.h>
#include <assert.h>
#endif


int main(int argc, char *argv[])
{
    std::string s = "/home/sherry/FF_VS/SSCA2v2.2/dumpbc.data";
    //std::string s = "/Users/xuepengfan/git/SSCA2v2.2/dumpbc.data";
    graph g = read_data(s);
    std::cout<<"G.n : "<<g.n<<std::endl;
    std::cout<<"G.m : "<<g.m<<std::endl;
    std::cout<<"k4approx: "<<g.k4approx<<std::endl;

#ifdef CACHE_EVAL
    /*Add papi to trace cache miss*/
    int EventSet,retVal;
    long long startRecords[1], endRecords[1];
    //long long startRecords[2], endRecords[2];
    retVal = PAPI_library_init(PAPI_VER_CURRENT);
    assert(retVal == PAPI_VER_CURRENT);
    EventSet = PAPI_NULL;
    retVal = PAPI_create_eventset(&EventSet);
    assert(retVal == PAPI_OK);
    //L1 TCM & TCA
    retVal = PAPI_add_event(EventSet, PAPI_L1_TCM);
    assert(retVal == PAPI_OK);
    //retVal = PAPI_add_event(EventSet, PAPI_L1_TCA);
    //assert(retVal == PAPI_OK);
    
    retVal = PAPI_start(EventSet);
    assert(retVal == PAPI_OK);
    retVal = PAPI_read(EventSet, startRecords);
    assert(retVal == PAPI_OK);
    /*Add papi to trace cache miss*/
#endif

    auto bc = ff_get_bc(g);
    //auto bc = seq_get_bc(g);
    //bc = ff_get_bc(g);

#ifdef CACHE_EVAL
    /*Stop papi trace*/
    retVal = PAPI_stop(EventSet, endRecords);
    assert(retVal == PAPI_OK);
    retVal = PAPI_cleanup_eventset(EventSet);
    assert(retVal == PAPI_OK);
    retVal = PAPI_destroy_eventset(&EventSet);
    assert(retVal == PAPI_OK);
    PAPI_shutdown(); 
    //L1 result
    std::cout << "L1 total cache miss = " << endRecords[0] - startRecords[0] << std::endl;
    //std::cout << "L1 total cache access = " << endRecords[0] - startRecords[0] << std::endl;
    //std::cout << "L1 total cache access = " << endRecords[1] - startRecords[1] << std::endl;
    /*Stop papi trace*/
#endif 
    //auto bc = boost::shared_ptr<DOUBLE_T[]>(new DOUBLE_T[g.n]);
    //betweennessCentrality(&g, bc.get());
    double sum1 = 0;
    double sum_origin = 0;
    for(int i = 0; i < g.n; i++)
    {
        sum1 += bc[i];
        sum_origin += g.BC[i];
    }

    std::cout<<"BC FF: "<<sum1<<std::endl;
    std::cout<<"BC Origin: "<<sum_origin<<std::endl;
    return 0;
}
