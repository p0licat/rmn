/*
    Tests for the gmp_matrix header.
    Checks the functionality of the class GmpMatrix.
*/
#ifndef TESTS_H // header guard
#define TESTS_H // TESTS_H

//
// C++ standard libraries.
//
#include <iostream>     // input/output
#include <string>       // uses std::strings

//
// Project headers.
//

// uncomment any of these ( at most one ) to run different kinds of test routines
//#define EXTREME
//#define COMPREHENSIVE
#define QUICK

#include "../include/gmp_matrix.h" // uses gmp_matrix

/*
    Tests header. For testing, several input files with matrices from sizes 3*3 to 1000*1000 have been generated
    using Octave/Matlab. The corresponding inverses have also been generated.

    The testing is done in two steps:
        1 - Testing by comparing with known result.
        2 - Testing by multiplying the standard matrix by its inverse and comparing result to the identity matrix.
*/

namespace test_func {

// utility functions
// compare first n characters from two strings
bool checkMatchFirstN(const int& n, const std::string& first, const std::string& other);

// standard tests by comparison with expected results
bool test_n_by_n(int n, std::string filename, std::string inverse_filename, std::string dumpfile);

// mathematical tests with identity matrix, with expected precision greater
// than 'decimals' decimal places

bool test_n_by_n_ident(const int& n, const int& decimals, const std::string& filename, const std::string& dumpfile);

// master function to run all tests
bool runTests();

} // namespace test_func

#endif // TESTS_H
