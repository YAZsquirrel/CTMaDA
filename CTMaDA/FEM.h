#pragma once
//#define DEBUG
#define DEBUG2
#define DEBUG_VECTOR_MKE

#include <iomanip>
#include <cmath>
#include <iostream>
#include <functional>
#include <list>
#include "Math_structs.h"

using namespace maths;
using integr_f = std::function<real(real, real, int, int, int[4])>;

class FEM
{
private:
   size_t num_of_knots, num_of_FEs, num_of_edges, un;

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

   //real Integrate2D(const std::function<real(real, real, int, int, int[4])> calc_f, int i, int j, int knot_num[4]);

   std::function<real(real, real, int, int, int[4])> Gij;
   std::function<real(real, real, int, int, int[4])> Mij;

public:
   size_t GetKnotsNum() { return num_of_knots; }
   size_t GetHexasNum() { return num_of_FEs; }
   size_t GetEdgesNum() { return num_of_edges; }

   std::vector<real>& GetKnots() { return q; };
   FEM(Mesh* _mesh);
   void SolveElliptic();
   void Output(std::ofstream& out);
};

#ifdef DEBUG_VECTOR_MKE


class VectorFEM
{
   size_t num_of_knots, num_of_FEs, num_of_edges, un;

   //real localM2d[4][4];
   real localC[4][4]{}; // 8*8
   real localG[4][4]{};
   real localA[4][4]{};
   real localb[4]{};
   real rJ[2][2]{};

   real J2D[2][2]{};

   real constexpr det_J();
   void calc_J(int edge_num[4], real ksi, real etta);
   void calc_vphi(int index, real ksi, real etta, real vphi[2]);
   void get_global_xy(real ksi, real etta, real xy[2], element2D& elem);
   real phi(int index, real ksi, real etta, int knot_num[4]);
   void get_func_by_local_coords(element2D& elem, real ksi, real etta, real sumq[2], real xy[2]);

#ifdef DEBUG2
   void calc_f(edge& k, int n_test, real F[2]);
   real bound1func(edge& k, int n_mat);
   real bound2func(bound& k, int n_mat);
#endif // DEBUG

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

   real Integrate(const integr_f f, int i, int j, int knot_num[4]);

   //real Integrate2D(const std::function<real(real, real, int, int, int[4])> calc_f, int i, int j, int knot_num[4]);

   integr_f dij;
   integr_f Mij;

public:
   size_t GetKnotsNum() { return num_of_knots; }
   size_t GetHexasNum() { return num_of_FEs; }
   size_t GetEdgesNum() { return num_of_edges; }

   std::vector<real>& GetKnots() { return q; };
   VectorFEM(Mesh* _mesh);
   void SolveElliptic();
   void Output(int point_per_FE_sqred);
   void CheckOnErrors();



};

#endif // DEBUG_VECTOR_MKE