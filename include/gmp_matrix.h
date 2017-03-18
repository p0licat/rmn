#ifndef GMPMATRIX_H // header guard
#define GMPMATRIX_H // GMPMATRIX_H

//
// C standard libraries.
//
#include <stdio.h>  // standard input/output, mpfr also uses printf
//
// C++ standard libraries.
//
#include <fstream>  // file input output
#include <iostream> // standard input output
#include <vector>   // STL data structure used

//
// Additional libraries.
//
#include <mpfr.h>   // GNU multiple precision library/multiple-precision correct rounding.


//
// Wrapper for mpfr.h library.
//
#include "mpreal.h"

#define DEFAULT_PREC 50000  // default precision used by MPFR

typedef std::vector<mpfr::mpreal> mpVector; // mpVector typedef
typedef std::vector<mpVector> mpMatrix;     // mpMatrix typedef

/*
    GmpMatrix is the implementation of a matrix data structure.
    This is represented internally as a std::vector of std::vectors.
    The elements of the vectors are of mpfr::mpreal type.
*/
class GmpMatrix
{
    // Matrix is of size n * m
    // The class handles internally a
    // mpMatrix(n)(m), as seen in the typedef.
    //
    // This means n lines by m columns.
    private:
        int n;  // lines
        int m;  // columns
        mpMatrix elements; // internal data structure

    public:
        // constructors
        GmpMatrix() {}          // void constructor override
        GmpMatrix(int n, int m) // constructor with specified size
        {
            mpfr_set_default_prec(DEFAULT_PREC); // the default_prec of MPFR data type is set
                                                 // as #defined in this file by DEFAULT_PREC

            this->n = n;    // set private member n to value of param n ( lin )
            this->m = m;    // set private member m to value of param m ( col )

            this->elements.resize(n);   // resize the mpMatrix to length n ( lines )

            // the loop below creates the columns for every line
            for (int i = 0; i < n; ++i) // for every line
            {


                //TODO: WHAT THE FUCK
                this->elements[i].resize(m); // give that line size m ( column )
                mpVector e_push(m); // create mpVector of size m

                for (int j = 0; j < m; ++j) // for every element of e_push
                {
                    mpfr::mpreal e;         // create mpreal element
                    e_push.push_back(e);    // and push it into e_push
                }
            }
        }

        ~GmpMatrix() {} // void destructor

        // getters
        int getN() { return this->n; }
        int getM() { return this->m; }

        // get element at index(i, j)
        mpfr::mpreal getElement(int i, int j) {
            return this->elements[i][j];
        }

        // setters
        void setElement(int i, int j, mpfr::mpreal e);

        // utility
        void copyMatrix(GmpMatrix& B);
        void readFromFile(const std::string& filename);
        void printToConsole();
        void swapElement(int i, int j, int i2, int j2);
        static bool multiplyMatrix(GmpMatrix& src1, GmpMatrix& src2, GmpMatrix& dest);

        //arithmetic
        mpfr::mpreal inverse_mpreal();

};

#endif // GMPMATRIX_H
