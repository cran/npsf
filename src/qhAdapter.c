//#include "libqhull.h"
#include "qhull_a.h"

int convhull(double *points, int n, int d, int *vertexIndices){
	char *options;
	if (d <= 3){
		options = "qhull Qt QbB Pp";
	}else{
		options = "qhull Qt Qx Qs QbB Pp";
	}
	//FILE *errfile = fopen("err.txt", "w");
	//FILE *tmpstdout = fopen("out.txt", "w");
	int exitcode = qh_new_qhull (d ,n , points, 0, options, NULL, 
		NULL);
	//int exitcode = qh_new_qhull (d ,n , points, 0, options, tmpstdout, 
	//	errfile);
	//fclose(tmpstdout);
	//fclose(errfile);
	if (!exitcode){
		//facetT *facet;
		vertexT *vertex, **vertexp;
		//unsigned int numFacets = qh num_facets;
		//std::vector<int> vertices(numFacets*d);
		int counter = 1;
		FORALLvertices{
			vertexIndices[counter++] = qh_pointid(vertex->point);
		}
		*vertexIndices = counter - 1;
/*		int i=0;
		FORALLfacets {
			int j=0;
			FOREACHvertex_ (facet->vertices) {
				if (j < d){
					vertices[i+n * j++] = 1 + qh_pointid(vertex->point);
				}
			}
			i++;
		}
		std::sort(vertices.begin(), vertices.end());
		int counter = 2;
		int curVertex = vertices[0];vertexIndices[1] = curVertex;
		for (int i = 0; i < numFacets * d; i++){
			if (vertices[i] > curVertex){
				curVertex = vertices[i];
				vertexIndices[counter++] = curVertex;
			}
		}
		*vertexIndices = counter - 1;*/
	}
	qh_freeqhull(qh_ALL);
	return exitcode;
}
