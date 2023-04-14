#include "FEM.h"
#include <fstream>
#include <set>
#define MAX_INNER_POINTS_IN_FE_TO_OUTPUT 5

VectorFEM::VectorFEM(Mesh* _mesh)
{
   mesh = _mesh;
   num_of_knots = mesh->knots.size();
   num_of_FEs = mesh->elems.size();
   num_of_edges = mesh->edges.size();

   this->A = MakeSparseFormat_withEdges(4, num_of_edges, mesh);
   q.resize(num_of_edges, 0.);
   b.resize(num_of_edges, 0.);

   Mij = [this](real ksi, real etta, int i, int j, int knot_num[4])
   {
      calc_J(knot_num, ksi, etta);
      // J^-1
      real det = det_J();
      {
         rvrs_J[0][0] = J2D[1][1] / det;
         rvrs_J[1][1] = J2D[0][0] / det;
         rvrs_J[1][0] = -J2D[1][0] / det;
         rvrs_J[0][1] = -J2D[0][1] / det;
      }
      real pi[2] = {0, 0}, pj[2] = {0, 0};
      calc_vphi(i, ksi, etta, pi);
      calc_vphi(j, ksi, etta, pj);

      real rJpi[2] = { rvrs_J[0][0] * pi[0] + rvrs_J[0][1] * pi[1],
                      rvrs_J[1][0] * pi[0] + rvrs_J[1][1] * pi[1]}, 

           rJpj[2] = { rvrs_J[0][0] * pj[0] + rvrs_J[0][1] * pj[1],
                      rvrs_J[1][0] * pj[0] + rvrs_J[1][1] * pj[1]};
      real res = rJpi[0] * rJpj[0] + rJpi[1] * rJpj[1];

      return res* abs(det);
   };

   dij = [this](real ksi, real etta, int i, int j, int knot_num[4])
   {
      calc_J(knot_num, ksi, etta);
      return 1. / abs(det_J());
   };

}

void VectorFEM::SolveElliptic()
{
   CreateSLAE();
   SolveSLAE_LOS(A, q, b);

   Output(MAX_INNER_POINTS_IN_FE_TO_OUTPUT);
}

void VectorFEM::Output(int point_per_FE_sqred)
{

   std::ofstream out("Q.txt");
   out << std::setprecision(7);

   for (int i = 0; i < num_of_edges; i++)
   {
      out << std::defaultfloat;
      out << i;
      out.width(30);
      out << std::scientific;
      out.width(15);
      out << q[i];
      out.width(15);
      out << "\n";
   }
   out.close();

   out.open("elements.txt");

   for (int i = 0; i < num_of_FEs; i++)
   {
      element2D& elem = mesh->elems[i];
      out << std::fixed;
      for (int i = 0; i < element2D::local_knots_num; i++)
         out << elem.knots_num[i] << " ";
      out << "\n";
   }
   out.close();


   out.open("Result_vectors.txt");
   out << std::setprecision(7);

   std::set<int> visited_knots;
   std::set<int> visited_edges;

   for (int i = 0; i < num_of_FEs; i++)
   {
      // FEs
      element2D& elem = mesh->elems[i];

      for (int k = 1; k < point_per_FE_sqred - 1; k++)
         for (int n = 1; n < point_per_FE_sqred - 1; n++)
         {
            real ksi = 2. * k / (point_per_FE_sqred - 1.) - 1., etta = 2. * n / (point_per_FE_sqred - 1.) - 1.;
            real sumq[2] = {0, 0};
            real xy[2] = { 0, 0 };

            get_func_by_local_coords(elem, ksi, etta, sumq, xy);

            out << xy[0] << ' ' << xy[1] << ' ' << sumq[0] << ' ' << sumq[1] << '\n';
         }

      // edges
      for (int e = 0; e < 4; e++)
      {
         if (visited_edges.contains(elem.edge_nums[e])) continue;

         visited_edges.insert(elem.edge_nums[e]);
         for (int k = 1; k < point_per_FE_sqred - 1; k++)
         {
            real ksi, etta;
            real xy[2] = { 0, 0 };
            real sumq[2] = { 0, 0 };
            switch (e)
            {
               case 0: ksi = -1; etta = 2. * k / (point_per_FE_sqred - 1.) - 1.; break;
               case 1: ksi = 1; etta =  2. * k / (point_per_FE_sqred - 1.) - 1.; break;
               case 2: etta = -1; ksi = 2. * k / (point_per_FE_sqred - 1.) - 1.; break;
               case 3: etta = 1; ksi =  2. * k / (point_per_FE_sqred - 1.) - 1.; break;
            }
            get_func_by_local_coords(elem, ksi, etta, sumq, xy);

            out << xy[0] << ' ' << xy[1] << ' ' << sumq[0] << ' ' << sumq[1] << '\n';
         }      
      }

      //knots
      for (int k = 0; k < 4; k++)
      {
         if (visited_knots.contains(elem.knots_num[k])) continue;

         visited_knots.insert(elem.knots_num[k]);

         real ksi = k % 2 ? 1 : -1, etta = k / 2 ? 1 : -1;
         real xy[2] = { 0, 0 };
         real sumq[2] = { 0, 0 };
         
         get_func_by_local_coords(elem, ksi, etta, sumq, xy);

         out << xy[0] << ' ' << xy[1] << ' ' << sumq[0] << ' ' << sumq[1] << '\n';
      }
   }

   out.close();
}

void VectorFEM::get_func_by_local_coords(element2D& elem, real ksi, real etta, real sumq[2], real xy[2])
{
   calc_J(elem.knots_num, ksi, etta);
   {
      real det = det_J();
      rvrs_J[0][0] = J2D[1][1] / det;
      rvrs_J[1][1] = J2D[0][0] / det;
      rvrs_J[1][0] = -J2D[1][0] / det;
      rvrs_J[0][1] = -J2D[0][1] / det;
   }

   for (int p = 0; p < 4; p++)
   {
      real vphi[2] = { 0, 0 };
      calc_vphi(p, ksi, etta, vphi);
      sumq[0] += q[elem.edge_nums[p]] * vphi[0];//(rvrs_J[0][0] * vphi[0] + rvrs_J[0][1] * vphi[1]);
      sumq[1] += q[elem.edge_nums[p]] * vphi[1];//(rvrs_J[1][0] * vphi[0] + rvrs_J[1][1] * vphi[1]);
   }

   get_global_xy(ksi, etta, xy, elem);
}

void VectorFEM::AddFirstBounds()
{
   for (auto& bound : mesh->bounds1)
   {
      int e_num = bound.edge_num;
      
      A->di[e_num] = 1.;
      for (int j = A->ig[e_num]; j < A->ig[e_num + 1]; j++)
         A->l[j] = 0.;

      for (int ii = 0; ii < A->dim; ii++)                // ���� �� ��������
         for (int j = A->ig[ii]; j < A->ig[ii + 1]; j++)   // ���� ��������� � �������
            if (A->jg[j] == e_num)          // � ������ ������ �������?
               A->u[j] = 0.;
      
      b[e_num] = bound1func(bound, (int)round(bound.value1)/*bound.n_test*/);

      MatSymmetrisation(A, b, e_num);
   }
}

void VectorFEM::AddSecondBounds()
{
   for (auto& bound : mesh->bounds2)
   {
      real l = mesh->length(bound);
      real xl = bound.knots[1].x - bound.knots[0].x, 
           yl = bound.knots[1].y - bound.knots[0].y;
      int be_num = bound.edge_num;
      b[be_num] += bound2func(bound, (int)round(bound.value1)) * 2. / l * (-bound.n.x * xl + bound.n.y * yl);

   }
}

void VectorFEM::CreateSLAE()
{
   element2D hexa;
   for (int i = 0; i < num_of_FEs; i++)
   {
      hexa = mesh->elems[i];
      CreateG(hexa);
      CreateM(hexa);
      AddToA(hexa);
      Createb(hexa);
   }

   AddSecondBounds();
   AddFirstBounds();
}

void VectorFEM::AddToA(element2D& hexa)
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         AddElement(A, hexa.edge_nums[i], hexa.edge_nums[j], hexa.lam * localG[i][j] + hexa.gam * localC[i][j]);
}

void VectorFEM::CreateM(element2D& hexa)
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localC[i][j] = Integrate(Mij, i, j, hexa.knots_num);
}

void VectorFEM::CreateG(element2D& hexa)
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localG[i][j] = Integrate(dij, i, j, hexa.knots_num) / 4. * (i == j || i + j == 3 ? 1. : -1.);
}

void VectorFEM::Createb(element2D& hexa)
{
   real localb[4]{};
   real f_[4]{};
   real t[2]{};
   real F[2]{};

   for (int i = 0; i < 4; i++)
   {
      edge e = mesh->edges[hexa.edge_nums[i]];
      knot k1 = mesh->knots[e.knots_num[0]], 
           k2 = mesh->knots[e.knots_num[1]];
      real l = mesh->length(e);

      t[0] = (k2.x - k1.x) / l;
      t[1] = (k2.y - k1.y) / l;
      calc_f(e, hexa.n_test, F);
      f_[i] = F[0] * t[0] + F[1] * t[1];
   }

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localb[i] += localC[i][j] * f_[j];

   for (int i = 0; i < 4; i++)
      b[hexa.edge_nums[i]] += localb[i];
}

void VectorFEM::calc_vphi(int index, real ksi, real etta, real vphi[2])
{
   switch (index)
   {
      case 0: vphi[1] = (1. - ksi) / 2.; break;
      case 1: vphi[1] = (1. + ksi) / 2.; break;
      case 2: vphi[0] = (1. - etta) / 2.; break;
      case 3: vphi[0] = (1. + etta) / 2.; break;
   }
}

void VectorFEM::get_global_xy(real ksi, real etta, real xy[2], element2D& elem)
{
   
   real x[4] = {mesh->knots[elem.knots_num[0]].x, mesh->knots[elem.knots_num[1]].x, mesh->knots[elem.knots_num[2]].x, mesh->knots[elem.knots_num[3]].x },
        y[4] = {mesh->knots[elem.knots_num[0]].y, mesh->knots[elem.knots_num[1]].y, mesh->knots[elem.knots_num[2]].y, mesh->knots[elem.knots_num[3]].y },
        p[4] = {phi(0, ksi, etta, elem.knots_num), phi(1, ksi, etta, elem.knots_num), phi(2, ksi, etta, elem.knots_num), phi(3, ksi, etta, elem.knots_num) };
   for (int i = 0; i < 4; i++)
   {
      xy[0] += x[i] * p[i];
      xy[1] += y[i] * p[i];
   }
}

real VectorFEM::phi(int index, real ksi, real etta, int knot_num[4])
{
   int sk = index % 2 ? 1 : -1, se = index / 2 ? 1 : -1;
   return (1. + se * etta) * (1. + sk * ksi) / 4.;
}

real VectorFEM::det_J()
{
   return J2D[1][1] * J2D[0][0] - J2D[1][0] * J2D[0][1];
}

void VectorFEM::calc_J(int knots_num[4], real ksi, real etta)
{
   real m1 = (1. - ksi) / 4., 
        m2 = (1. - etta) / 4.,
        m3 = (1. + ksi) / 4.,
        m4 = (1. + etta) / 4.;
   real v0[2] = {mesh->knots[knots_num[0]].x, mesh->knots[knots_num[0]].y};
   real v1[2] = {mesh->knots[knots_num[1]].x, mesh->knots[knots_num[1]].y};
   real v2[2] = {mesh->knots[knots_num[2]].x, mesh->knots[knots_num[2]].y};
   real v3[2] = {mesh->knots[knots_num[3]].x, mesh->knots[knots_num[3]].y};

   J2D[0][0] = m2 * (v1[0] - v0[0]) + m4 * (v3[0] - v2[0]); J2D[1][0] = m2 * (v1[1] - v0[1]) + m4 * (v3[1] - v2[1]);
   J2D[0][1] = m1 * (v2[0] - v0[0]) + m3 * (v3[0] - v1[0]); J2D[1][1] = m1 * (v2[1] - v0[1]) + m3 * (v3[1] - v1[1]);
}

real VectorFEM::Integrate(const std::function<real(real, real, int, int, int[4])> f, int i, int j, int knot_nums[4])
{
   const int nKnot = 3;//5; // Knots num

   const real xj[nKnot]
      = { -.7745966692414833, 0., .7745966692414833 }; // sqrt(0.6)
   //= { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// Scales
   //           0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

   const real qj[nKnot]
      = { .55555555555555555, .8888888888888888, .55555555555555555 };
   //= { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// Weights
   //               (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };

   real result = 0.;
   for (int ix = 0; ix < nKnot; ix++)
      for (int iy = 0; iy < nKnot; iy++)
         result += qj[ix] * qj[iy] * (f(xj[ix], xj[iy], i, j, knot_nums));
   return result;
}

