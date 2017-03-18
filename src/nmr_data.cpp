
#include "../include/nmr_data.h"



GmpMatrix NMR_DATA::GenerateMagnetizationMatrix(const int& dim)
{
    GmpMatrix Tau(dim, 1);  // vector Tau
    GmpMatrix T2(dim, 1);   // vector T2

    // dTau = 0.5 * 10^(-3)
    mpr dTau = 0.5 * mpfr::pow(10.0, -3.0);

    // Tau[i] = (i+1) * dTau
    for (int line = 0; line < dim; ++line)
        Tau.setElement(line, 0, (line+1) * dTau);

    mpr T2min = Tau.getElement(0, 0);           // T2min =      Tau[0]
    mpr T2max = 2 * Tau.getElement(dim-1, 0);   // T2max =  2 * Tau[n-1]

    mpr lgT2min = mpfr::log10(T2min);   // lgT2min = log10(T2min)
    mpr lgT2max = mpfr::log10(T2max);   // lgT2max = log10(T2max)

    mpr dlgT2 = (lgT2max - lgT2min) / (dim - 1);

    for (int line = 0; line < dim; ++line)
        T2.setElement(line, 0, mpfr::pow(10, lgT2min + line * dlgT2));

    GmpMatrix M(dim, dim);

    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            M.setElement(i, j, mpfr::exp(-1 * Tau.getElement(i, 0) / T2.getElement(j, 0)));

    return M;
}

