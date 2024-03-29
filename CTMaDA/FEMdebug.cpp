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
#endif // DEBUG
