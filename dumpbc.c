#include "defs.h"
#include <stdio.h>
/*THis file is for SSCA2 */

/* To use this function to dump BC result, you should do the following things:
 * 1. add the declaration in defs.h, like this:
 *      double dumpbc(graph *, DOUBLE_T *);
 * 2. append dumpbc.o to OBJS in Makefile
 * 3. call dumpbc in SSCA2.c: main after betweennessCentrality(G, BC);
 */
double dumpbc(graph *G, DOUBLE_T * BC)
{

  printf("G->n :%d (size: %d)\n", G->n, sizeof(G->n));
  printf("G->m :%d (size: %d)\n", G->m, sizeof(G->m));
  printf("k4approx :%d (size: %d)\n", K4approx, sizeof(K4approx));
  LONG_T n = G->n;
  LONG_T m = G->m;
  FILE *fp = fopen("dumpbc.data", "wb");
  fwrite(&(G->n), sizeof(G->n),1, fp);
  fwrite(&(G->m), sizeof(G->m),1, fp);
  fwrite(G->endV, sizeof(VERT_T), m, fp);
  fwrite(G->numEdges, sizeof(LONG_T), n, fp);
  fwrite(G->weight, sizeof(WEIGHT_T), m, fp);
  fwrite(BC, sizeof(DOUBLE_T), n, fp);
  fwrite(&K4approx, sizeof(K4approx),1, fp);

  fclose(fp);
  return 1.0;
}
