#include "defs.h"


int main(int argc, char *argv[])
{
    std::string s = "/Users/xuepengfan/git/SSCA2v2.2/dumpbc.data";
    graph g = read_data(s);
    std::cout<<"G.n : "<<g.n<<std::endl;
    std::cout<<"G.m : "<<g.m<<std::endl;
    std::cout<<"k4approx: "<<g.k4approx<<std::endl;


    //auto bc = seq_get_bc(g);
    auto bc = seq_get_bc(g);
    bc = ff_get_bc(g);

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
