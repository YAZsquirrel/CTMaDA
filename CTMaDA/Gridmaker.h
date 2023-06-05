#pragma once

#include <vector>
typedef double real;


struct knot
{
   //unsigned int knot_num;
   knot(real _x, real _y, real _z) : x(_x), y(_y), z(_z) {}
   knot(real _x, real _y) : x(_x), y(_y), z(0) {}
   knot() : x(0.0), y(0.0), z(0.0) {}
   real x, y, z;
   std::vector<int> edges_num;
   knot& operator=(const knot& k) {
      if (this == &k) return *this;
      x = k.x;
      y = k.y;
      z = k.z;
      return *this;
   }

   bool operator==(const knot* k) {

      bool equal = (abs(k->x - x) < 1e-10 &&
         abs(k->y - y) < 1e-10 &&
         abs(k->z - z) < 1e-10) || *this == k;
      return equal;
   }

   bool operator==(const knot& k) {
      bool equal = (abs(k.x - x) < 1e-10 &&
         abs(k.y - y) < 1e-10 &&
         abs(k.z - z) < 1e-10);
      return equal;
   }

   bool operator!=(const knot& k) {
      bool equal = (abs(k.x - x) < 1e-10 &&
         abs(k.y - y) < 1e-10 &&
         abs(k.z - z) < 1e-10);
      return !equal;
   }

   bool operator!=(const knot* k) {
      bool equal = (abs(k->x - x) < 1e-10 &&
         abs(k->y - y) < 1e-10 &&
         abs(k->z - z) < 1e-10) || *this == k;

      return !equal;
   }
};

struct edge {
   int knots_num[2]{};
   std::vector<int> elems_num;
   knot knots[2];

   edge& operator=(const edge& elem) {
      if (this == &elem) return *this;

      for (int i = 0; i < 2; i++)
         knots_num[i] = elem.knots_num[i];
      return *this;
   }

   
};

struct bound : edge{
   knot n;
   int edge_num;
   int n_mat = -1;
   int n_test = -1;
   real value1 = 0;      // ug, th, ub
   real value2 = 0;      // beta
};

struct element2D {
   real lam = 0., gam = 0.;
   int n_test = 0;
   const static int local_knots_num = 4;
   int n_mat = -1;
   int knots_num[4]{};
   knot knots[4];
   int edge_nums[4]{};

   bool containsKnot(int n)
   {
      bool found = false;
      for (int i = 0; i < local_knots_num && !found; i++)
         found = knots_num[i] == n;
      return found;
   }
   bool containsKnot(knot& k)
   {
      bool found = false;
      for (int i = 0; i < local_knots_num && !found; i++)
         found = knots[i] == k;
      return found;
   }

   element2D& operator=(const element2D& elem) {
      if (this == &elem) return *this;
      lam = elem.lam;
      gam = elem.gam;
      n_test = elem.n_test;
      n_mat = elem.n_mat;
      for (int i = 0; i < local_knots_num; i++)
         knots_num[i] = elem.knots_num[i];
      for (int i = 0; i < 4; i++)
         edge_nums[i] = elem.edge_nums[i];
      return *this;
   }

   element2D(real gamma, real lambda, int _n_mat_test, int knots_nums[4])
      : lam(lambda), gam(gamma), n_test(_n_mat_test)
   {
      for (int i = 0; i < local_knots_num; i++)
         knots_num[i] = knots_nums[i];
   }
   element2D() : lam(0), gam(0), n_test(0) {}
};

class Mesh
{
public:
   std::vector<element2D> elems;
   std::vector<knot> knots;
   std::vector<bound> bounds1;
   std::vector<bound> bounds2;
   std::vector<bound> bounds3;
   std::vector<edge> edges;
   void MakeMesh();
   knot& Cross(knot& k1, knot& k2);
   real length(knot& k1, knot& k2);
   real length(edge& e);
   bool onSegment(knot& p, knot& k1, knot& k2);
   bool isInTriangle(knot& k1, knot& k2, knot& k3, knot& p);
   int GetEdgeNumByKnots(knot& k1, knot& k2);

private:
   void SetBoundConds();
   void SetElemParameters();
   void RemoveNullKnots();
   void FindEdges();

   void output();
   struct material
   {
      real gam, lam;
      int n_test;
      material( real _gam, real _l, int _n_test)
         :  gam(_gam), lam(_l), n_test(_n_test) {}
   };
   struct area
   {
      int n_mat;
      int X1, X2, Y1, Y2;
      area(int nmat, int x1, int x2, int y1, int y2)
         : n_mat(nmat), X1(x1), X2(x2), Y1(y1), Y2(y2) {}
   };
   struct boundEdge
   {
      int n_bound;
      int p[2];
      int q;
      int axis;
      real v1, v2;
   };

   std::vector<material> mats;
   std::vector<area> areas;
   std::vector<boundEdge> bes;
   std::vector<int> IX, IY;
   std::vector<knot> XY_full;
   std::vector<knot> XY;
   int nX, nY, x_size, y_size;

   //void Set???();

};

