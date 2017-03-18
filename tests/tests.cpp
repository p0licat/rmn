/*
    Implementation of tests.h
*/

#include "tests.h"

/*
    Compares the strings a and b. Depending on which is longer, the function
    will compare at least the first min(a.size(), b.size()) characters.
    If the position of the first mismatch is less than or equal to n, it will return FALSE.
    If the position of the first mismatch is strictly greater than n, it will return TRUE.

    Input:
        n - (const int&) number of characters to compare
        first - (const std::string&) first string
        other - (const std::string&) second string

*/
bool test_func::checkMatchFirstN(const int& n, const std::string& first, const std::string& other)
{
    int upto = std::min(first.size(), other.size()); // iterate shortest
    int miss = -1; // position of first mismatch

    // iterate shortest string
    for (int i = 0; i < upto; ++i)
    {
        // if strings don't match
        if ( first[i] != other[i] )
        {
            miss = i; // store position
            break; // stop loop
        }
    }

    // if whole string matched don't leave miss -1
    (miss == -1) ? miss = upto - 1:0;

    // if strings didn't match or miss was lower than n
    if ( miss == -1 || miss <= n )
        return false; // don't match with enough precision
    return true; // sufficiently match
}

/*
    Standard test function. Tests an n by n matrix by comparing the inverse calculated
    from the matrix inside filename with the already known inverse matrix from filename.

    Input:
        n                   - size of matrix
        filename            - name of the file containing matrix to inverse
        inverse_filename    - name of the file containing the inverse of the matrix inside filename
        dumpfile            - file in which to dump data for analysis
    Output:
        true if inverse is good, false if it is not

    This function also dumps the contents of the computed inverse and read inverse inside a
    dump file.
*/
bool test_func::test_n_by_n(int n, std::string filename, std::string inverse_filename, std::string dumpfile)
{
    std::cout << "Launching " << n << " by " << n << " test.\n\n";
    std::cout << "Creating A(" << n << ", " << n << ") matrix.\n";

    GmpMatrix A(n, n); // creating a GmpMatrix object of dimension n * n

    std::cout << "A(" << n << ", " << n << ") created. Printing: \n";
    A.printToConsole(); // printing A to console

    std::cout << "Creating matrix invA(" << n << ", " << n << ") for storing the inverse.\n";
    std::cout << "Creating matrix result(" << n << ", " << n << ") for storing the pre-generated inverse.\n";
    GmpMatrix invA(n, n);   // creating GmpMatrix object
    GmpMatrix result(n, n); // creating GmpMatrix object
    std::cout << "Done!\n";

    std::cout << "Reading A from file" << filename << "\n";
    A.readFromFile(filename); // reading A from file = 'filename'
    std::cout << "Done! Printing A(" << n << ", " << n << "):\n";
    A.printToConsole(); // print newly read A to console

    std::cout << "Reading expected result from mat3by3inverse.txt into 'resut(3, 3)':\n";
    result.readFromFile(inverse_filename); // read the inverse from file = 'inverse_filename'
    std::cout << "Done! Printing result(" << n << ", " << n << "): \n";
    result.printToConsole(); // print inverse to console

    std::cout << "Copying matrix A(" << n << ", " << n << ") to invA(" << n << ", " << n << ")... ";
    invA.copyMatrix(A); // copy A into invA in order to invert invA
    std::cout << "Done!\n";
    std::cout << "Calculating inverse of A(" << n << ", " << n << ")... ";
    invA.inverse_mpreal(); // inverse invA
    std::cout << "Done! \n";
    std::cout << "Printing invA(" << n << ", " << n << "): \n";
    invA.printToConsole(); // print inverse to console

    std::cout << "Comparing the computed inverse invA(" << n << ", " << n << ") with expected result result(" << n << ", " << n << ")... ";
    for (int i = 0; i < A.getN(); ++i) // for every line
        for (int j = 0; j < A.getM(); ++j ) // for every column
            // compare first N digits
            if ( checkMatchFirstN(4, result.getElement(i, j).toString(), invA.getElement(i, j).toString() ) == false ) {
                std::cout << "Error!\nResults do not match!\n";
                return false;
            }
    std::cout << "Done!\nResults match.\n";

    // contents are also dumped to dumpfile so that they can be visualized
    std::cout << "Dumping test results to " << dumpfile << " ...";

    std::ofstream t_out(dumpfile);
    t_out << "<p_graph>\n";
    t_out << "#computed\n";
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            t_out << "X[" << i << "][" << j << "]: " << invA.getElement(i, j) << '\n';
    }
    t_out << "#expected\n";
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            t_out << "X[" << i << "][" << j << "]: " << invA.getElement(i, j) << '\n';
    }
    t_out << "</p_graph>\n";

    t_out.close();
    std::cout << "Done!\n";
    return true;
}

/*
    Tests the functionality of the inverse and multiply operations from the GmpMatrix class,
    by performing a series of calculations that lead to the identity matrix. This requires a
    matrix input file, with n lines and n entries per line.

    Input:
        n - (const int&) the size of the matrix ( n * n )
        decimals - (const int&) number of decimals for which to verify precision
        filename - (const std::string&) name of the file to read the matrix from
        dumpfile - (const std::string&) debugging file in which to dump data
*/
bool test_func::test_n_by_n_ident(const int& n, const int& decimals, const std::string& filename, const std::string& dumpfile)
{
    std::cout << "Generating matrix A(" << n << ", " << n << ") ... ";
    GmpMatrix A(n, n);
    std::cout << "Done!\n";

    std::cout << "Generating matrix invA(" << n << ", " << n << ") ... ";
    GmpMatrix invA(n, n);
    std::cout << "Done!\n";

    std::cout << "Reading A from file ... TODO";
    A.readFromFile(filename);
    std::cout << "Done!\n";

    std::cout << "Copying matrix A(" << n << ", " << n << ") into invA(" << n << ", " << n << ") ... ";
    invA.copyMatrix(A);
    std::cout << "Done!\n";

    std::cout << "Calculating the inverse into invA(" << n << ", " << n << ") ... ";
    invA.inverse_mpreal();
    std::cout << "Done!\n";

    std::cout << "Creating matrix mulMat(" << n << ", " << n << ") ... ";
    GmpMatrix mulMat(n, n);
    std::cout << "Done!\n";

    std::cout << "Multiplying A(" << n << ", " << n << ") with its inverse invA(" << n << ", " << n << ") and storing the result into mulMat(" << n << ", " << n << ") ... ";
    GmpMatrix::multiplyMatrix(A, invA, mulMat);
    std::cout << "Done!\n";

    // these mpreal variables are used to check that the main diagonal of the result is 1.000...
    // and that the other elements are 0.000 within a specified precision.
    //
    // the result of the multiplication should be the identity matrix I(n)

    mpfr::mpreal one = 1.0;     // 1.0000... ( precision specified in gmp_matrix.h )
    mpfr::mpreal zero = 0.0;    // 0.0000... ( precision specified in gmp_matrix.h )

    // test result to up to 300 decimal places approximation to one or zero for identity matrix
    std::cout << "Comparing elements of mulMat(" << n << ", " << n << ") with the identity matrix I(" << n << ") \
    with precision of " << decimals << " decimals ... ";

    std::string format_tostring = "%." + std::to_string(decimals) + "Rf";
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if ( i == j )
            {
                if ( checkMatchFirstN(decimals, mpfr::abs(mulMat.getElement(i, j)).toString(format_tostring), one.toString(format_tostring)) == false )
                {
                    std::cout << "Test ERROR! Following strings didn't match: \n";
                    std::cout << "At position i:" << i << ", j:" << j << '\n';
                    std::cout << mpfr::abs(mulMat.getElement(i, j)).toString(format_tostring) << '\n';
                    std::cout << one.toString(format_tostring) << '\n';
                    std::cout << "TERMINATING!\n";
                    return false;
                }
            }
            else
                if ( checkMatchFirstN(decimals, mpfr::abs(mulMat.getElement(i, j)).toString(format_tostring), zero.toString(format_tostring)) == false )
                {
                    std::cout << "Test ERROR! Following strings didn't match: \n";
                    std::cout << "At position i:" << i << ", j:" << j << '\n';
                    std::cout << mpfr::abs(mulMat.getElement(i, j)).toString(format_tostring) << '\n';
                    std::cout << zero.toString(format_tostring) << '\n';
                    std::cout << "TERMINATING!\n";
                    return false;
                }

    std::cout << "Done! \nMatrix mulMat(" << n << ", " << n << ") and I(" << n << ") match.\n";
    return true;
}

/*
    Main function that is ran in order to test the functionality.
*/
bool test_func::runTests()
{
#ifdef QUICK
    // compare with matlab
    if ( test_n_by_n(3, "../resources/tests/mat3by3.txt", "../resources/tests/mat3by3inverse.txt", "../outputs/tests_output.dat") == false )
    {
        std::cout << "3by3 failed\n";
        return false;
    }

    if ( test_n_by_n(4, "../resources/tests/mat4by4.txt", "../resources/tests/mat4by4inverse.txt", "../outputs/tests_output.dat") == false )
    {
        std::cout << "4 by 4 failed\n";
        return false;
    }
    #ifdef COMPREHENSIVE
    // compare with matlab
    if ( test_n_by_n(100, "../resources/tests/mat100by100.txt", "../resources/tests/mat100by100inverse.txt", "../outputs/tests_output.dat") == false )
    {
        std::cout << "100 by 100 failed\n";
        return false;
    }

    // compare by identity matrix
    if ( test_n_by_n_ident(100, 300, "../resources/tests/mat100by100.txt", "../outputs/tests_output.dat") == false )
    {
        std::cout << "Identity test 100*100 failed.\n";
        return false;
    }
    #endif
    std::cout << "All tests ran!\n";

#endif
    return true;
}
