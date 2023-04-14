#include "FEM.h"
int main()
{
	Mesh* grid = new Mesh();
	grid->MakeMesh();
	VectorFEM* fem = new VectorFEM(grid);
	fem->SolveElliptic();
	return 0;
}