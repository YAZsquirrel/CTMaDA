#pragma once
//#define DEBUG
#define DEBUG2
#include <iomanip>
#include <cmath>
#include <iostream>
#include <functional>
#include <list>
#include "Math_structs.h"

using namespace maths;

class FEM
{
private:
   size_t num_of_knots, num_of_FE, un;

   //real localM2d[4][4];
   real localM[4][4]; // 8*8
   real localG[4][4];
   real localA[4][4];
   real reversed_J[2][2];

   real J2D[2][2];
   real Jgrad_i[2];
   real gradi[2];
   real Jgrad_j[2];
   real gradj[2];
   inline real det_J();
   real prime_by_var(int what, int varOnFE, int knot_num[4], real ksi, real etta);
   inline int mu(int index);
   inline int v(int index);
   real W(int index, real alpha);
   real d_phi(int index, int what, real ksi, real etta);
   inline real phi(int index, real ksi, real etta);
   void calc_grad(int ij, int index, real ksi, real etta);

#ifdef DEBUG2
   real f(knot& k, element2D& e);
   real bound1func(knot& k, int n_mat);
   real bound2func(knot& k, int n_mat);
   real bound3func(knot& k, int n_mat);
   real bound3funcbeta(knot& k, int n_mat);
#endif // DEBUG

   void check_test();
   void AddFirstBounds();
   void AddSecondBounds();
   void AddToA(element2D& hexa);
   void CreateSLAE();
   void CreateM(element2D& hexa);
   void CreateG(element2D& hexa);
   void Createb(element2D& hexa);

   Matrix* A;
   std::vector<real>b;
   std::vector<real> q;
   Mesh* mesh;

   real Integrate(const std::function<real(real, real, int, int, int[4])> f, int i, int j, int knot_num[4]);

   //real Integrate2D(const std::function<real(real, real, int, int, int[4])> f, int i, int j, int knot_num[4]);

   std::function<real(real, real, int, int, int[4])> Gij;
   std::function<real(real, real, int, int, int[4])> Mij;

public:
   size_t GetKnotsNum() { return num_of_knots; }
   size_t GetHexasNum() { return num_of_FE; }

   std::vector<real>& GetKnots() { return q; };
   FEM(Mesh* _mesh);
   void SolveElliptic();
   void Output(std::ofstream& out);
};
