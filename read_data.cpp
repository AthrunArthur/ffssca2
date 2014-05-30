#include "defs.h"

/*
 *   FILE *fp = fopen("dumpbc.data", "wb");
  fwrite(&(G->n), sizeof(G->n),1, fp);
  fwrite(&(G->m), sizeof(G->m),1, fp);
  fwrite(G->endV, sizeof(VERT_T), m, fp);
  fwrite(G->numEdges, sizeof(LONG_T), n, fp);
  fwrite(G->weight, sizeof(WEIGHT_T), m, fp);
  fwrite(BC, sizeof(DOUBLE_T), n, fp);
  fwrite(&K4approx, sizeof(K4approx),1, fp);
*/
#include <stdio.h>

graph read_data(const std::string & fpath)
{
    graph g;
    FILE * fp = fopen(fpath.c_str(), "rb");
    fread(&(g.n), sizeof(g.n), 1, fp);
    fread(&(g.m), sizeof(g.m), 1, fp);
    g.endV = boost::shared_array<VERT_T>(new VERT_T[g.m]);
    fread(g.endV.get(), sizeof(VERT_T), g.m, fp);
    g.numEdges = boost::shared_array<LONG_T>(new LONG_T[g.n]);
    fread(g.numEdges.get(), sizeof(LONG_T), g.n, fp);

    g.weight = boost::shared_array<WEIGHT_T>(new WEIGHT_T[g.m]);
    fread(g.weight.get(), sizeof(WEIGHT_T), g.m, fp);

    g.BC = boost::shared_array<DOUBLE_T>(new DOUBLE_T[g.n]);
    fread(g.BC.get(), sizeof(DOUBLE_T), g.n, fp);

    fread(&g.k4approx, sizeof(g.k4approx), 1, fp);
    fclose(fp);
    return g;
}
