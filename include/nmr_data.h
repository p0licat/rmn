/*
    Header for solving NMR-related problems using the GmpMatrix class.
*/

#include <iostream>
#include "gmp_matrix.h"

#ifndef __NMR_DATA__ // header guard
#define __NMR_DATA__ // __NMR_DATA__

typedef mpfr::mpreal mpr;

namespace NMR_DATA {

    GmpMatrix GenerateMagnetizationMatrix(const int& dim);



}



#endif // __NMR_DATA__



