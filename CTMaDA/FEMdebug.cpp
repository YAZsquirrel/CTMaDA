#include "FEM.h"
#include <math.h>

#ifdef DEBUG2
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


#ifdef DEBUG_VECTOR_MKE

void VectorFEM::calc_f(edge& k, int n_test, real F[])
{
   real x = (k.knots[0].x + k.knots[1].x) / 2.,
        y = (k.knots[0].y + k.knots[1].y) / 2.;
   switch (n_test)
   {
   case 10: 
      F[0] = 0;
      F[1] = 5; break;

   case 20:
      F[0] = 1 + y;
      F[1] = x; break;
   }
}

real VectorFEM::bound1func(bound& b, int n_test)
{
   real x = (b.knots[0].x + b.knots[1].x) / 2.,
        y = (b.knots[0].y + b.knots[1].y) / 2.;

   knot k1 = b.knots[0],
        k2 = b.knots[1];
   real l = mesh->length(b);

   real t[2] = { (k2.x - k1.x) / l, (k2.y - k1.y) / l };
   real Q[2];
   switch (n_test)
   {
      case 10: Q[0] = 0; Q[1] = 5; break;

      case 20: Q[0] = 1 + y; Q[1] = 0; break;
      case 21: Q[0] = x; Q[1] = 0; break;

   default:
      return 0;
   }


   return t[0] * Q[0] + t[1] * Q[1];
}

real VectorFEM::bound2func(bound& b, int n_test) // THETA = (dAy/dx - dAx/dy), THETA' = 1/mu * THETA
{

   switch (n_test)
   {
   case 20: return 0; 

   case 10: ; break; // <
   case 11: ; break; // v
   case 12: ; break; // >
   case 13: ; break; // ^

   default:
      return 0;
   }
}
#endif


#endif

