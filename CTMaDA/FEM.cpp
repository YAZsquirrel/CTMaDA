#include "FEM.h"
#include <fstream>

FEM::FEM(Mesh* _mesh)
{
    mesh = _mesh;
    num_of_knots = mesh->knots.size();
    num_of_FE = mesh->elems.size();
    this->A = MakeSparseFormat(4, num_of_knots, mesh);
    q.resize(num_of_knots, 0.);
    b.resize(num_of_knots, 0.);

    Mij = [this](real ksi, real etta, int i, int j, int knot_num[4])
    {
       for (int ip = 0; ip < 2; ip++)
          for (int jp = 0; jp < 2; jp++)
             J2D[ip][jp] = prime_by_var(ip, jp, knot_num, ksi, etta);

       return phi(i, ksi, etta) * phi(j, ksi, etta) * abs(det_J());
    };

    Gij = [this](real ksi, real etta, int i, int j, int knot_num[4])
    {
       for (int ip = 0; ip < 2; ip++)
          Jgrad_i[ip] = Jgrad_j[ip] = 0.;

       for (int ip = 0; ip < 2; ip++)                                       // | dx/d(ksi)   | dy/d(ksi) 
          for (int jp = 0; jp < 2; jp++)                                    // | dx/d(etta)  | dy/d(etta)
             J2D[ip][jp] = prime_by_var(ip, jp, knot_num, ksi, etta);  

       // J^-1
       {
          real det = det_J();
          reversed_J[0][0] = J2D[1][1] / det;
          reversed_J[1][1] = J2D[0][0] / det;
          reversed_J[1][0] = -J2D[1][0] / det;
          reversed_J[0][1] = -J2D[0][1] / det;
       }

       // grad(phi(ksi, etta, theta))
       calc_grad(1, i, ksi, etta);
       calc_grad(2, j, ksi, etta);

       // J^-1 * grad(phi)
       for (int ip = 0; ip < 2; ip++)
          for (int jp = 0; jp < 2; jp++)
          {
             Jgrad_i[ip] += reversed_J[ip][jp] * gradi[jp];
             Jgrad_j[ip] += reversed_J[ip][jp] * gradj[jp];
          }

       // Jgrad_i^T * Jgrad_j
       real res = 0;
       for (int ip = 0; ip < 2; ip++)
          res += Jgrad_i[ip] * Jgrad_j[ip];
       return res * abs(det_J());
    };

 }

 void FEM::SolveElliptic()
 {
   CreateSLAE();
   SolveSLAE_LOS(A, q, b);
   std::ofstream out("Result.txt");

#ifdef DEBUG2
   check_test();
#endif // DEBUG
   Output(out);
   out.close();
}

void FEM::Output(std::ofstream& out)
{
   //out.scientific;
   //out.precision(15);
   //std::cout.scientific;
   //std::cout.precision(15);
   out.setf(std::ios::right);
   out << "| x" << std::fixed;
   out.width(15);
   out << "| y";
   out.width(15);
   out << "| q";
   out.width(15);
   out << "|\n";
   //std::cout << title;
   out << std::setprecision(7);

   for (int i = 0; i < num_of_knots; i++)
   {
      out << std::defaultfloat;
      out << mesh->knots[i].x;
      out.width(15);
      out << mesh->knots[i].y;
      out.width(30);
      out << std::scientific;
      out.width(15);
      out<< q[i];
      out.width(15);
      out << "\n";
   }
}

void FEM::check_test()
{
   real sum = 0.;
   real sumq = 0.;
   for (int i = 0; i < q.size(); i++)
   {  
      sum += abs(q[i] - bound1func(mesh->knots[i], mesh->bounds1[0].n_mat));
      sumq += abs(bound1func(mesh->knots[i], mesh->bounds1[0].n_mat));
   }

   std::cout << "Abs error: " << sum / sumq << "\nRelative error: " << sum;
}

void FEM::AddFirstBounds()
{
   for (auto& cond : mesh->bounds1)
   {
      for (int i = 0; i < 2; i++)
      {
         A->di[cond.e.knots_num[i]] = 1.;
         for (int j = A->ig[cond.e.knots_num[i]]; j < A->ig[cond.e.knots_num[i] + 1]; j++)
            A->l[j] = 0.;
         for (int ii = 0; ii < A->dim; ii++)                // идем по столбцам
            for (int j = A->ig[ii]; j < A->ig[ii + 1]; j++)   // идем элементам в столбце
               if (A->jg[j] == cond.e.knots_num[i])          // в нужной строке элемент?
                  A->u[j] = 0.;
#ifdef DEBUG2
         b[cond.e.knots_num[i]] = bound1func(mesh->knots[cond.e.knots_num[i]], cond.n_test);
#elif
         b[cond.e.knots_num[i]] = cond.value1;
#endif // DEBUG2


         MatSymmetrisation(A, b, cond.e.knots_num[i]);
      }
   }
}

void FEM::AddSecondBounds()
{
   for (auto& bound : mesh->bounds2)
   {  
      real h = mesh->length(bound.e.knots[0], bound.e.knots[1]);
      real M[2][2] = {{2,1},
                      {1,2}};

#ifdef DEBUG2

      for (int i = 0; i < 2; i++)
         b[bound.e.knots_num[i]] += bound2func(bound.e.knots[i], (int)round(bound.value1)) * (M[i][0] + M[i][1]) * h / 6;
#elif
      for (int i = 0; i < 4; i++)
         b[bound.e.knots_num[i]] += bound.value1 * (localM[i][0] + localM[i][1]) / h;

#endif // DEBUG2

   }
}

void FEM::CreateSLAE()
{
   element2D hexa;
   for (int i = 0; i < num_of_FE; i++)
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

void FEM::AddToA(element2D& hexa)
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         AddElement(A, hexa.knots_num[i], hexa.knots_num[j], localG[i][j] + localM[i][j]);
}

void FEM::CreateM(element2D& hexa)
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localM[i][j] = hexa.gam * Integrate(Mij, i, j, hexa.knots_num);
}
void FEM::CreateG(element2D& hexa)
{
  for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++)
        localG[i][j] = hexa.lam * Integrate(Gij, i, j, hexa.knots_num);
}

void FEM::Createb(element2D& hexa)
{
   real localb[4]{};
   real f_[4]{};
   for (int i = 0; i < 4; i++)
      f_[i] = f(mesh->knots[hexa.knots_num[i]], hexa);

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localb[i] += localM[i][j] * f_[j];

   for (int i = 0; i < 4; i++)
      b[hexa.knots_num[i]] += localb[i];
}

int FEM::mu(int index)
{
   return index % 2;
}

int FEM::v(int index)
{
   return index / 2;
}

real FEM::W(int index, real alpha)
{
   if (!index) return 1. - alpha;
   return alpha;
}

real FEM::d_phi(int index, int var, real ksi, real etta)
{
   real d_phi = 0.;
   switch (var)
   {
   case 0:    // ksi
   {
      d_phi = W(v(index), etta);
      if (!mu(index)) d_phi *= -1;
      break;
   }
   case 1:     // etha
   {
      d_phi = W(mu(index), ksi);
      if (!v(index)) d_phi *= -1;
      break;
   }
   }

   return d_phi;
}

real FEM::prime_by_var(int varOnSquare, int varOnFE, int knot_num[4], real ksi, real etta)
{
   real var = 0.;
   for (int i = 0; i < 4; i++)
   {
      switch (varOnFE)
      {
      case 0: var += mesh->knots[knot_num[i]].x * d_phi(i, varOnSquare, ksi, etta); break;
      case 1: var += mesh->knots[knot_num[i]].y * d_phi(i, varOnSquare, ksi, etta); break;
      }
   }
   return var;
}

real FEM::phi(int index, real ksi, real etta)
{
   return  W(mu(index), ksi) * W(v(index), etta);
}

real FEM::det_J()
{
   return J2D[1][1] * J2D[0][0] - J2D[1][0] * J2D[0][1];
}

real FEM::Integrate(const std::function<real(real, real, int, int, int[4])> f, int i, int j, int knot_nums[4])
{
   const int nKnot = 3;//5; // Knots num

   const real xj[nKnot]
      = { .7745966692414833, 0., -.7745966692414833 }; // sqrt(0.6)
   //= { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// Scales
   //           0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

   const real qj[nKnot]
      = { .55555555555555555, .8888888888888888, .55555555555555555 };
   //= { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// Weights
   //               (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };

   real result = 0.;
   for (int ix = 0; ix < nKnot; ix++)
      for (int iy = 0; iy < nKnot; iy++)
          result += qj[ix] * qj[iy]  * (f(.5 + xj[ix] * .5, .5 + xj[iy] * .5, i, j, knot_nums));
   return result / 4.;
}

void FEM::calc_grad(int ij, int index, real ksi, real etta)
{
   switch (ij)
   {
   case 1:
      for (int i = 0; i < 2; i++)
         gradi[i] = d_phi(index, i, ksi, etta);
      break;
   case 2:
      for (int i = 0; i < 2; i++)
         gradj[i] = d_phi(index, i, ksi, etta);
      break;
   }
}

