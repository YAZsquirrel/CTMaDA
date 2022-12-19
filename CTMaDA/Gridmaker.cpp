#include "Gridmaker.h"
#include <fstream>

knot& intersection(knot p1, knot p2, knot q1, knot q2)
{
   knot k;
   real d = (p1.x - p2.x) * (q1.y - q2.y) - (q1.x - q2.x) * (p1.y - p2.y);
   k.x = ((p1.x * p2.y - p1.y * p2.x ) * (q1.x - q2.x) - (q1.x * q2.y - q1.y * q2.x) * (p1.x - p2.x)) / d;
   k.y = ((p1.x * p2.y - p1.y * p2.x ) * (q1.y - q2.y) - (q1.x * q2.y - q1.y * q2.x) * (p1.y - p2.y)) / d;

   return k;
}

void Mesh::MakeMesh()
{
   std::ifstream fmat("Materials.txt");
   int k;
   fmat >> k;
   mats.reserve(k);

   for (int i = 0; i < k; i++)
   {
      real g, l;
      int n_test;
      fmat >> g >> l >> n_test;

      material mat(g, l, n_test);
      mats.push_back(mat);
   }
   fmat.close();

   std::ifstream fgrid("GridDescription.txt");
   fgrid >> nX >> nY;

   std::vector<int> Xs, Ys;
   std::vector<real> Xr, Yr;
   XY.reserve(nX * nY);
   Xs.reserve(nX - 1);
   Ys.reserve(nY - 1);
   Xr.reserve(nX - 1);
   Yr.reserve(nY - 1);

   for (size_t i = 0; i < nX * nY; i++)
   {
      real x, y;
      fgrid >> x >> y;
      knot k(x, y);
      XY.push_back(k);
   }

   for (size_t i = 0; i < nX - 1; i++)
   {
      real xs;
      fgrid >> xs;
      Xs.push_back(xs);
   }
   for (size_t i = 0; i < nY - 1; i++)
   {
      real ys;
      fgrid >> ys;
      Ys.push_back(ys);
   }
   for (size_t i = 0; i < nX - 1; i++)
   {
      real xr;
      fgrid >> xr;
      Xr.push_back(xr);
   }
   for (size_t i = 0; i < nY - 1; i++)
   {
      real yr;
      fgrid >> yr;
      Yr.push_back(yr);
   }

   int nBound;
   fgrid >> nBound;
   for (size_t i = 0; i < nBound; i++)
   {
      int bound_n;
      real value1, value2 = 0.;
      int p1, p2,  q;
      int axis;

      fgrid >> bound_n;
      if (bound_n == 3)
         fgrid >> value1 >> value2 >> p1 >> p2 >> q >> axis;
      else
         fgrid >> value1 >> p1 >> p2 >>  q >> axis;
      
      boundEdge b;

      b.p[0] = p1;
      b.p[1] = p2;
      b.q = q;

      b.v1 = value1;
      b.v2 = value2;
      b.n_bound = bound_n;
      b.axis = axis;
      bes.push_back(b);
   }


   int nAreas;
   fgrid >> nAreas;

   for (size_t i = 0; i < nAreas; i++)
   {
      int mat_n;
      int x1, x2, y1, y2;
      fgrid >> mat_n >> x1 >> x2 >> y1 >> y2;
      area a(mat_n, x1, x2, y1, y2);
      areas.push_back(a);
   }
   fgrid.close();

   //код для создания сетки
   elems.reserve((nX - 1) * (nY - 1));

   IX.reserve(nX);
   IY.reserve(nY);
   IX.push_back(0);
   IY.push_back(0);
   x_size = nX; y_size = nY;
   for (size_t i = 0; i < Xs.size(); i++)
   {
      x_size += Xs[i] - 1;
      IX.push_back(IX[i] + Xs[i]);
   }
   for (size_t i = 0; i < Ys.size(); i++)
   {
      y_size += Ys[i] - 1;
      IY.push_back(IY[i] + Ys[i]);
   }
   XY_full.resize(x_size * y_size, knot(INFINITY, INFINITY));

   for (size_t iyg = 0; iyg < nY; iyg++)
   {
      //X_g.push_back(knot(XY[iyg * nX].x, XY[iyg * nX].y)); 
      XY_full[IY[iyg] * x_size] = knot(XY[iyg * nX].x, XY[iyg * nX].y);

      for (size_t ixg = 1; ixg < nX; ixg++)
      {
         real x[2] = { XY[iyg * nX + (ixg - 1)].x, 
                       XY[iyg * nX + ixg].x };
         real y[2] = { XY[iyg * nX + (ixg - 1)].y,
                       XY[iyg * nX + ixg].y };
         
         bool isRegular = !(abs(Xr[ixg - 1] - 1.) > 1e-10);
         real xl = x[0],
              xr = x[1],
              xh = xr - xl;
         real xb1 = isRegular ? xh / Xs[ixg - 1] :
            xh * (1. - Xr[ixg - 1]) / (1. - pow(Xr[ixg - 1], Xs[ixg - 1]));
         real x1 = xl;
         //Z_full.push_back(z);
         for (size_t n = 0; n < Xs[ixg - 1]; n++)
         {
            real scale = isRegular ? 1. : pow(Xr[ixg - 1], n);
            x1 += xb1 * scale;
            real tx = (xr - x1) / (xr - xl);
            //     ty = (yr1 - y1) / (yr1 - yl1);
            knot xy = knot(x1, y[0] + (y[1] - y[0]) * tx);
            //X_g.push_back(xy);
            XY_full[IY[iyg] * x_size + IX[ixg - 1] + 1 + n] = xy;
         }
      }
   }

   for (size_t ixg = 0; ixg < nX; ixg++)
   {
      for (size_t iyg = 1; iyg < nY; iyg++)
      {
         real x[2] = { XY[(iyg - 1) * nX + ixg].x,
                       XY[iyg * nX + ixg].x };
         real y[2] = { XY[(iyg - 1) * nX + ixg].y,
                       XY[iyg * nX + ixg].y };

         bool isRegular = !(abs(Yr[iyg - 1] - 1.) > 1e-10);
         real yl = y[0],
            yr = y[1],
            yh = yr - yl;
         real yb1 = isRegular ? yh / Ys[iyg - 1] :
            yh * (1. - Yr[iyg - 1]) / (1. - pow(Yr[iyg - 1], Ys[iyg - 1]));
         real y1 = yl;
         //Z_full.push_back(z);
         for (size_t n = 0; n < Ys[iyg - 1]; n++)
         {
            real scale = isRegular ? 1. : pow(Yr[iyg - 1], n);
            y1 += yb1 * scale;
            real ty = (yr - y1) / (yr - yl);
            //     ty = (yr1 - y1) / (yr1 - yl1);
            knot xy = knot(x[0] + (x[1] - x[0]) * ty, y1);
            //Y_g.push_back(xy);

            XY_full[(IY[iyg - 1] + 1 + n) * x_size + IX[ixg]] = xy;
         }
      }
   }

   knot inf = knot(INFINITY, INFINITY);
   for (size_t iy = 1; iy < nY; iy++)
      for (size_t iys = 0; iys < Ys[iy - 1]; iys++)
         for (size_t ix = 1; ix < nX; ix++)
            for (size_t ixs = 0; ixs < Xs[ix - 1]; ixs++)
            {
               size_t i = IX[ix - 1] + ixs, j = IY[iy - 1] + iys;
               if (XY_full[i + j * x_size].x == INFINITY || XY_full[i + j * x_size].y == INFINITY)
               {
                  knot kl = XY_full[j * x_size + IX[ix - 1]],
                       kr = XY_full[j * x_size + IX[ix]],
                       kb = XY_full[IY[iy - 1] * x_size + i],
                       ku = XY_full[IY[iy] * x_size + i] ;

                  knot xy = intersection(kl, kr, kb, ku);
                  XY_full[i + j * x_size] = xy;
               }
            }

   knots.reserve(XY_full.size());

   for (size_t ixy = 0; ixy < XY_full.size(); ixy++)
   {
      knot k = knot(XY_full[ixy].x, XY_full[ixy].y);
      knots.push_back(k);
   }

   for (size_t j = 0; j < y_size - 1; j++)
      for (size_t i = 0; i < x_size - 1; i++)
      {
         int n_mat = -1;
         int ks[4] { i      +  j      * x_size, 
                    (i + 1) +  j      * x_size,  
                     i      + (j + 1) * x_size, 
                    (i + 1) + (j + 1) * x_size};
         element2D e;
         for (size_t i = 0; i < 4; i++)
            e.knots_num[i] = ks[i];
         elems.push_back(e);
      }

   for (auto& el : elems)
      for (size_t i = 0; i < 4; i++)
         el.knots[i] = knots[el.knots_num[i]];
   
   SetElemParameters();
   SetBoundConds();
   RemoveNullKnots();
   XY.clear();
   XY_full.clear();

   output();
}

knot& Mesh::Cross(knot& k1, knot& k2)
{
   knot k = knot(k1.y * k2.z - k1.z * k2.y, k1.x * k2.z - k1.z * k2.x, k1.x * k2.y - k1.y * k2.x);
   return k;
}

inline real Mesh::length(knot& k1, knot& k2)
{
   return sqrt(pow(k1.x - k2.x, 2) + pow(k1.y - k2.y, 2));
}

inline bool Mesh::onSegment(knot& p, knot& k1, knot& k2)
{
   return abs( length(k1, k2) - (length(p, k1) + length(p, k2))) < 1e-10;
}

bool Mesh::isInTriangle(knot& k1, knot& k2, knot& k3, knot& p)
{
   knot k12 = knot(k2.x - k1.x, k2.y - k1.y), 
        k23 = knot(k3.x - k2.x, k3.y - k2.y), 
        k31 = knot(k1.x - k3.x, k1.y - k3.y), 
        pk1 = knot(p.x - k1.x, p.y - k1.y), 
        pk2 = knot(p.x - k2.x, p.y - k2.y), 
        pk3 = knot(p.x - k3.x, p.y - k3.y);
   real pk12 = Cross(pk1, k12).z, pk23 = Cross(pk2, k23).z, pk31 = Cross(pk3, k31).z;

   return (pk12 < 1e-12 && pk23 < 1e-12 && pk31 < 1e-12) || (pk12 >= 1e-12 && pk23 >= 1e-12 && pk31 >= 1e-12);
}

void Mesh::output()
{
   std::ofstream of1("knots.txt");
   of1 << knots.size() << '\n';
   for (size_t i = 0; i < knots.size(); i++)
      of1 << knots[i].x << " " << knots[i].y << " " << knots[i].z << '\n';
   of1.close();

   std::ofstream of2("elements.txt");
   of2 << elems.size() << '\n';
   for (size_t i = 0; i < elems.size(); i++)
   {
      of2 << elems[i].n_mat << ' ';
      for (size_t j = 0; j < 4; j++)
         of2 << elems[i].knots_num[j] << " ";
      of2 << '\n';
   }
   of2.close();

   std::ofstream of3("bounds1.txt");
   of3 << bounds1.size() << '\n';
   for (auto& b : bounds1)
      of3 << b.knots_num[0] << " " << b.knots_num[1] << " " << b.value1 << '\n';
   of3.close();

   std::ofstream of4("bounds2.txt");
   of4 << bounds2.size() << '\n';
   for (auto& b : bounds2)
      of4 << b.knots_num[0] << " " << b.knots_num[1] << " " << b.value1 << '\n';
   of4.close();

   std::ofstream of5("bounds3.txt");
   of5 << bounds3.size() << '\n';
   for (auto& b : bounds3)
      of5 << b.knots_num[0] << " " << b.knots_num[1] << " " << b.value1 << " " << b.value2 << '\n';
   of5.close();

}

void Mesh::SetBoundConds()
{
   for (auto& e : elems)
   {
      for (auto& be : bes)
      {
         int ks[4] = { e.knots_num[0], e.knots_num[1],
                       e.knots_num[2], e.knots_num[3]};

         bool onEdge = false;
         int index = 0;

         knot k_g[2];
         int kn[2];

         for (size_t i = be.p[0]; i < be.p[1] && !onEdge; i++)
         {
            int p1 = be.axis == 1 ? IX[i]     : IY[i], 
                p2 = be.axis == 1 ? IX[i + 1] : IY[i + 1],
                q = be.axis == 1 ? IX[be.q] : IY[be.q];
                
            if (be.axis == 1) // y
            {
               k_g[0] = XY_full[p1 + q * x_size];
               k_g[1] = XY_full[p2 + q * x_size];
               // v
               if (onSegment(knots[ks[0]], k_g[0], k_g[1]) && onSegment(knots[ks[1]], k_g[0], k_g[1]))
               {
                  kn[0] = ks[0];
                  kn[1] = ks[1];
                  onEdge = true;
               }
               // ^
               if (onSegment(knots[ks[2]], k_g[0], k_g[1]) && onSegment(knots[ks[3]], k_g[0], k_g[1]))
               {
                  kn[0] = ks[2];
                  kn[1] = ks[3];
                  onEdge = true;
               }
            }
            else
            {
               k_g[0] = XY_full[p1 * x_size + q];
               k_g[1] = XY_full[p2 * x_size + q];
               // >
               if (onSegment(knots[ks[1]], k_g[0], k_g[1]) && onSegment(knots[ks[3]], k_g[0], k_g[1]))
               {
                  kn[0] = ks[1];
                  kn[1] = ks[3];
                  onEdge = true;
               }
               // <
               if (onSegment(knots[ks[0]], k_g[0], k_g[1]) && onSegment(knots[ks[2]], k_g[0], k_g[1]))
               {
                  kn[0] = ks[0];
                  kn[1] = ks[2];
                  onEdge = true;
               }
            }

         }
         

         if (onEdge)
         {
            bound b;
            b.knots_num[0] = kn[0];
            b.knots_num[1] = kn[1];
            b.knots[0] = knots[kn[0]];
            b.knots[1] = knots[kn[1]];

            b.value1 = be.v1;
            b.value2 = be.v2;
            b.n_mat = e.n_mat;
            b.n_test = e.n_test;
            switch (be.n_bound)
            {
            case 1:
               bounds1.push_back(b);
               break;
            case 2:
               bounds2.push_back(b);
               break;
            case 3:
               bounds3.push_back(b);
               break;
            }
         }
      }
   }

}

void Mesh::SetElemParameters()
{
   for (auto& e : elems)
   {
      int n_mat = -1;
      knot center = knot((e.knots[0].x + e.knots[1].x + e.knots[2].x + e.knots[3].x) / 4.,
                         (e.knots[0].y + e.knots[1].y + e.knots[2].y + e.knots[3].y) / 4.);

      for (auto& a : areas)
      {
         knot k_g[4];
         int kn[4];

         bool isIn = false;

         for (size_t i = a.X1; i < a.X2 && !isIn; i++)
         {
            for (size_t j = a.Y1; j < a.Y2 && !isIn; j++)
            {
               int p1 = IX[i], p2 = IX[i + 1], 
                   t1 = IY[j], t2 = IY[j + 1];
               k_g[0] = knots[p1 + t1 * x_size];
               k_g[1] = knots[p2 + t1 * x_size];
               k_g[2] = knots[p1 + t2 * x_size];
               k_g[3] = knots[p2 + t2 * x_size];

               if (isInTriangle(k_g[0], k_g[1], k_g[2], center) ||
                   isInTriangle(k_g[1], k_g[2], k_g[3], center))
               {
                  n_mat = a.n_mat; isIn = true;
               }
               else n_mat = -1;
            }
         }
         if (isIn) break;
      }

      if (n_mat > -1)
      {
         e.n_mat = n_mat;
         e.gam = mats[n_mat].gam;
         e.lam = mats[n_mat].lam;
         e.n_test = mats[n_mat].n_test;

      }
   }

   std::vector<int> toRemove;
   for (size_t i = 0; i < elems.size(); i++)
      if (elems[i].n_mat == -1)
         toRemove.push_back(i);

   for (size_t j = 0; j < toRemove.size(); j++)
   {
      elems.erase(elems.begin() + toRemove[j]);
      for (size_t jr = j; jr < toRemove.size(); jr++)
         toRemove[jr]--;
   }

}

void Mesh::RemoveNullKnots()
{
   //if there's knots that are orphans
   std::vector<int> toRemove;
   toRemove.reserve(knots.size() / 2);

   for (auto& k : knots)
   {
      bool isOrphan = true;
      for (auto& el : elems)
         if (el.containsKnot(k))
         {
            isOrphan = false;
            break;
         }

      for (int i = 0; i < knots.size() && isOrphan; i++)
         if (&knots[i] == &k)
         {
            toRemove.push_back(i);
            break;
         }
   }

   //remove 'em all
   for (size_t j = 0; j < toRemove.size(); j++)
   {
      knots.erase(knots.begin() + toRemove[j]);
      for (size_t jr = j; jr < toRemove.size(); jr++)
         toRemove[jr]--;
   }

   //renumerate
   for (auto& el : elems)
      for (size_t k = 0; k < 4; k++)
         if (el.knots_num[k] >= knots.size())
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == el.knots[k])
               {
                  el.knots_num[k] = i;
                  break;
               }
         }
         else if (knots[el.knots_num[k]] != el.knots[k])
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == el.knots[k])
               {
                  el.knots_num[k] = i;
                  break;
               }
         }

   for (auto& b : bounds1)
      for (size_t k = 0; k < 2; k++)
         if (b.knots_num[k] >= knots.size())
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == b.knots[k])
               {
                  b.knots_num[k] = i;
                  break;
               }
         }
         else if (knots[b.knots_num[k]] != b.knots[k])
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == b.knots[k])
               {
                  b.knots_num[k] = i;
                  break;
               }
         }

   for (auto& b : bounds2)
      for (size_t k = 0; k < 2; k++)
         if (b.knots_num[k] >= knots.size())
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == b.knots[k])
               {
                  b.knots_num[k] = i;
                  break;
               }
         }
         else if (knots[b.knots_num[k]] != b.knots[k])
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == b.knots[k])
               {
                  b.knots_num[k] = i;
                  break;
               }
         }

   for (auto& b : bounds3)
      for (size_t k = 0; k < 2; k++)
         if (b.knots_num[k] >= knots.size())
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == b.knots[k])
               {
                  b.knots_num[k] = i;
                  break;
               }
         }
         else if (knots[b.knots_num[k]] != b.knots[k])
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == b.knots[k])
               {
                  b.knots_num[k] = i;
                  break;
               }
         }
}
