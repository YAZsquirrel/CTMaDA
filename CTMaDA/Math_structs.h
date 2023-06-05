#pragma once

#include "Gridmaker.h"
#include <iostream>
#include <fstream>
#include <vector>

typedef double real;

namespace maths {
   enum MatrixFormat { Dense = 1, SparseRowColumn, SparseProfile, SparseRow };

   struct Matrix {
      std::vector<real> l, u, di, gg;
      std::vector<int> ig, jg;
      std::vector<std::vector<real>> dense;
      size_t dim = 0;
      MatrixFormat format = MatrixFormat::Dense;

      Matrix(MatrixFormat _format = Dense) : format(_format)
      {}
      ~Matrix()
      {
         l.clear(), u.clear(), di.clear(), gg.clear(), ig.clear(), jg.clear();
         for (int i = 0; i < dense.size(); i++)
            dense[i].clear();
         dense.clear();
      }
   };


   Matrix* MakeSparseFormat(int localsize, size_t size, Mesh* mesh);
   Matrix* MakeSparseFormat_withEdges(int localsize, size_t size, Mesh* mesh);

   void copy(std::vector<real>& to, std::vector<real>& from);
   void copy(std::vector<int>& to, std::vector<int>& from);
   real scalar(std::vector<real>& v, std::vector<real>& u);
   void AddElement(Matrix* M, int i, int j, real elem);
   void MatxVec(std::vector<real>& v, Matrix* A, std::vector<real>& b);
   void SolveSLAE_LOS(Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void WriteMatrix(Matrix* M);

   void MatSymmetrisation(Matrix* M, std::vector<real>& b, int i);
   Matrix* MakeHolessky(Matrix* A);
   void SolveForL(std::vector<real>& q, std::vector<real>& b, Matrix* M);
   void SolveForU(std::vector<real>& q, std::vector<real>& b, Matrix* M);

   real SLAEResidualOutput(std::vector<real>& q, maths::Matrix* M, std::vector<real>& b);

}