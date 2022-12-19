#include "FEM.h"
int main()
{
	Mesh* grid = new Mesh();
	grid->MakeMesh();
	FEM *fem = new FEM(grid);
	fem->SolveElliptic();
	return 0;
}