
#include <iostream>
#include <cstdint>
#include <vector>
#include <queue>
#include <stack>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

typedef int VERT_T;
typedef int LONG_T;
typedef VERT_T WEIGHT_T;
typedef double DOUBLE_T;


struct graph
{
    LONG_T n;
    LONG_T m; /* AA: number of edges */

    /* Directed edges out of vertex vi (say, k edges -- v1, v2, ... vk)
     * are stored in the contiguous block endV[numEdges[i] .. numEdges[i+1]]
     * So, numEdges[i+1] - numEdges[i] = k in this case */
    boost::shared_array<VERT_T> endV;
    boost::shared_array<LONG_T> numEdges;
    boost::shared_array<WEIGHT_T> weight;

    //VERT_T* endV; /* num : m*/
    //LONG_T* numEdges; /* AA : num- n*/
    //WEIGHT_T* weight; /* AA: num: m*/


    ////////////These are for storage conv
    boost::shared_array<DOUBLE_T> BC;
    int k4approx;
};

/* Predecessor list data structure for Kernel 4 */
struct plist
{
    VERT_T* list;
    VERT_T count;
    LONG_T degree;
};


graph read_data(const std::string & fp);
//double betweennessCentrality(graph* G, DOUBLE_T* BC);

//return BC
boost::shared_ptr<DOUBLE_T[]>  seq_get_bc(graph g);
boost::shared_ptr<DOUBLE_T[]>  ff_get_bc(graph g);

void prefix_sums(LONG_T*, LONG_T*, LONG_T*, LONG_T);


