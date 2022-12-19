#include "FEM.h"
#include <math.h>

#ifdef DEBUG2
real FEM::f(knot& k, element2D& e)
{
   switch (e.n_test)
   {
   case 0: return k.x;

   case 1: return -6 * k.x + k.x * k.x * k.x;	// check M
   case 2: return k.x;	// check G + M
   case 3: return k.y;
   case 4: return 51;

      return 0;
   }
}

real FEM::bound1func(knot& k, int n_test)
{
   switch (n_test)
   {
   case 0: return k.x;  // r

   case 1: return k.x * k.x * k.x;  // r
   case 2: return k.x;	// r
   case 3: return k.y;  // z
   case 4: return 51;   // 51

   default:
      return 0;
   }
}

real FEM::bound2func(knot& k, int n_test)
{
   switch (n_test)
   {
   case 5: return k.y;

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
