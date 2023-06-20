#include "FEM.h"
#include <math.h>

#ifdef DEBUG2
#pragma region scalar


//namespace FEM
//{ }
real FEM::f(knot& k, element2D& e)
{
   switch (e.n_test)
   {
   case 0: return k.y * k.x;


   default:
      return 0;
   }
}

real FEM::bound1func(knot& k, int n_test)
{
   switch (n_test)
   {
   case 0: return k.y * k.x;

   default:
      return 0;
   }
}

real FEM::bound2func(knot& k, int n_test)
{
   switch (n_test)
   {

   case 10: return k.x; // dy
   case 11: return -k.x; // -dy
   case 20: return k.y; // dx
   case 21: return -k.y; // -dx

   default:
      return 0;
   }
}

real FEM::bound3func(knot& k, int n_test)
{
   switch (n_test)
   {
   case 1: return 2. + k.y * k.y;

   default:
      return 0;
   }
}

real FEM::bound3funcbeta(knot& k, int n_test)
{
   switch (n_test)
   {
   case 1: return 2. + k.y * k.y;

   default:
      return 0;
   }
}

#pragma endregion

#ifdef DEBUG_VECTOR_MKE

void VectorFEM::calc_f(edge& k, int n_test, real F[])
{
   real x = (k.knots[0].x + k.knots[1].x) / 2.,
        y = (k.knots[0].y + k.knots[1].y) / 2.;
   calc_f_in_point(x, y, n_test,F);
}

void VectorFEM::calc_f_in_point(real x, real y, int n_test, real F[2])
{
   switch (n_test)
   {
   case 10:
      F[0] = 5;
      F[1] = 5; break;

   case 20:
      F[0] = 1 + y;
      F[1] = 0; break;

   case 21:
      F[0] = 0;
      F[1] = x; break;

   case 30:
      F[0] = 1 + y * y - 2;
      F[1] = 0; break;
   }
}

real VectorFEM::bound1func(bound& b, int n_test)
{
   knot k1 = b.knots[0],
        k2 = b.knots[1];

   real x = (k1.x + k2.x) / 2.,
        y = (k1.y + k2.y) / 2.;
   real l = mesh->length(b);

   real t[2] = {-b.n.y, b.n.x};//{ (k2.x - k1.x) / l, (k2.y - k1.y) / l };
   real Q[2]{};
   switch (n_test)
   {
      case 10: Q[0] = 5;
   				Q[1] = 5; break;

      case 20: Q[0] = 1 + y;
   				Q[1] = 0; break;

      case 21: Q[0] = 0;
   				Q[1] = x; break;

   	case 30: Q[0] = 1 + y * y;
					Q[1] = 0; break;

   default:
      return 0;
   }


   return t[0] * Q[0] + t[1] * Q[1];
}

real VectorFEM::bound2func(bound& b, int n_test) // THETA = (dAy/dx - dAx/dy), THETA' = 1/mu * THETA
{
   knot k1 = b.knots[0],
      k2 = b.knots[1];

   real x = (k1.x + k2.x) / 2.,
      y = (k1.y + k2.y) / 2.;

   switch (n_test)
   {
   case 20: return -1; 
   case 21: return 1;
   case 30: return -2 * y;

   case 10: return 0; break; // <
   case 11: ; break; // v
   case 12: ; break; // >
   case 13: ; break; // ^

   default:
      return 0;
   }
}
#endif


#endif

