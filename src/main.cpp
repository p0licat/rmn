//
//    C++ standard headers.
//
#include <iostream>
#include <string>
//
//    Project headers.
//
#include "../include/gmp_matrix.h"
#include "../tests/tests.h"
#include "../include/nmr_data.h"
using namespace std;

//#define NO_TESTS                  // comment to run tests
#define DATA_FILE "CadereN010.dat"  // input file
#define OUTPUT_FILE "result.dat"    // output file

// REVERSE_DEFINES
#define REV_SIZE 500

/*
Attempts to open file with the name 'filename' and
counts the number of entries in that file ( lines of separated floating point numbers ).

Example:
    int rval = getInputLength('data.txt');
    would return -1 if file does not exist or cannot be opened
    or            n where n is the number of entries in that file
Sample data.txt:
    .010621052631578948	3.8984240551981477
    .021142105263157896	2.383230633518354
    .03166315789473685	1.5564231628647518
    .042184210526315796	1.0634494659959695
    .05270526315789474	.7518145984230009
    .0632263157894737	.5460431397768444
    .07374736842105264	.4054471268565628
    .08426842105263159	.3066768350620069
    .09478947368421053	.23566438359265135
    .10531052631578948	.18359445343369396
    .11583157894736844	.14476026868399805
*/
int getInputLength(std::string filename)
{
    ifstream fin(filename.c_str());
    if ( !fin.is_open() )
        return -1;

    int counter = 0;
    std::string a, b;
    while (fin >> a >> b)
    {
        if ( ( a[0] == '.' || a[1] == '.' || a[2] == '.' || a[3] == '.' ) &&  ( b[0] == '.' || b[1] == '.' || b[2] == '.' || b[3] == '.') )
            counter++;
    }
    fin.close();
    return counter;
}

void Solve(int entries)
{
    GmpMatrix tau(entries, 1);
    GmpMatrix M(1, entries);
    std::cout << "constructed 1d matrix\n"; // TODO
    ifstream fin(DATA_FILE);

    std::string a, b;
    int ctr = 0;
    while ( fin >> a >> b )
    {
        if ( ( a[0] == '.' || a[1] == '.' || a[2] == '.' || a[3] == '.' ) &&  ( b[0] == '.' || b[1] == '.' || b[2] == '.' || b[3] == '.') )
        {
            mpfr::mpreal tau_i = a;
            mpfr::mpreal M_i = b;

            tau.setElement(ctr, 0, tau_i);
            M.setElement(0, ctr, M_i);
        }
        ctr++;
    }
    std::cout << "read data.\n"; // TODO
    mpfr::mpreal T2min = tau.getElement(0, 0);
    mpfr::mpreal T2max = 2*tau.getElement(entries-1, 0);

    mpfr::mpreal lgT2min = mpfr::log10(T2min);
    mpfr::mpreal lgT2max = mpfr::log10(T2max);

    mpfr::mpreal dlgT2 = (lgT2max - lgT2min) / (double)(entries-1);

    GmpMatrix T2(1, entries);
    for (int j = 0; j < entries; ++j)
    {
        mpfr::mpreal e = mpfr::pow("10.0", (lgT2min + j*dlgT2));
        T2.setElement(0, j, e);
    }

    GmpMatrix A(entries, entries);
    std::cout << "Constructed n*n matrix.\n"; // TODO
    for (int i = 0; i < entries; ++i)
    {
        for (int j = 0; j < entries; ++j)
        {
            mpfr::mpreal e = mpfr::exp((-1*tau.getElement(i, 0))/T2.getElement(0, j));
            A.setElement(j, i, e);
        }
    }

    GmpMatrix Acpy(entries, entries);
    Acpy.copyMatrix(A);

    A.inverse_mpreal();

    GmpMatrix res(entries, entries);

    GmpMatrix::multiplyMatrix(Acpy, A, res);
    res.printToConsole();

    std::cout << "Calculated inverse. \n"; // TODO
    GmpMatrix F(1, entries);
    if ( GmpMatrix::multiplyMatrix(M, A, F) == false )
        std::cout << "Invalid proportions.\n";

    // F is known

    ofstream fout(OUTPUT_FILE);

    for (int i = 0; i < entries; ++i)
    {
        mpfr::mpreal a = F.getElement(0, i);
        mpfr::mpreal b = T2.getElement(0, i);

        fout << b << ' ' << a << '\n';
    }

    fout.close();
    fin.close();
}

void ReverseSolve(int entries)
{
    ofstream fout("reverse_result.out");
    mpfr_set_default_prec(10000);
    int N = REV_SIZE;
    mpfr::mpreal tmin = 0.0001;
    mpfr::mpreal tmax = 0.2;
    mpfr::mpreal dtau = (tmax - tmin)/(N-1);

    mpfr::mpreal tau[REV_SIZE];
    for (int i = 0; i < N; ++i)
        tau[i] = tmin + i*dtau;


    typedef mpfr::mpreal mpr;

    mpr T2min = tau[0];
    mpr T2max = 2*tau[N-1];

    mpr T2lgmin = mpfr::log10(T2min);
    mpr T2lgmax = mpfr::log10(T2max);

    mpr dlgT2 = ( T2lgmax - T2lgmin ) / ( N - 1 );

    mpr T2[REV_SIZE];
    for (int i = 0; i < N; ++i)
        T2[i] = mpfr::pow("10.0", T2lgmin + i*dlgT2);

    mpr A1 = 1;
    mpr A2 = 2;
    mpr F[REV_SIZE];

    for (int i = 0; i < N; ++i)
        F[i] = A1 * mpfr::exp( -1*mpfr::pow((i-0.25*N), 2) / mpfr::pow((0.1*N), 2) ) +
               A2 * mpfr::exp( -1*mpfr::pow((i-0.6*N), 2) / mpfr::pow((0.1*N), 2) );

    mpr sigTau[REV_SIZE];
    for (int j = 0; j < N; ++j)
        sigTau[j] = F[j] * mpfr::exp(-1*tau[j] / T2[j] );

    mpr S[REV_SIZE];
    for (int i = 0; i < N; ++i)
    {
        mpr sig = 0;
        for (int j = 0; j < N; ++j)
            sig += F[j] * exp(-1*tau[i] / T2[j]);
        S[i] = sig;
    }

    for (int i = 0; i < N; ++i)
        fout << tau[i] << ' ' << S[i].toString("%.32Rf") << '\n';

    fout.close();
}

int main()
{
#ifndef NO_TESTS
    if ( test_func::runTests() == false )
    {
        std::cout << "Can not continue, some tests failed!\n";
        return 0;
    }
#endif
/*
    int entries = getInputLength(DATA_FILE);
    std::cout << "File contains " << entries << " valid entries. \n";

    if ( entries == -1 )
    {
        std::cout << "Aborting...\n";
        return 0;
    }
*/
    NMR_DATA::GenerateMagnetizationMatrix(20);

    return 0;
}
