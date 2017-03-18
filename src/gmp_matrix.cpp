/*
    Implementation of the GmpMatrix class, as defined in gmp_matrix.h
*/
#include "../include/gmp_matrix.h"

/*
    Read data from a file.
    Note that data must fit into the matrix, this method performs no resizes.

    Input:
        filename - (const std::string&) filename
*/
void GmpMatrix::readFromFile(const std::string& filename)
{
    // try to open 'fin' stream
    std::ifstream fin(filename.c_str());

    if ( !fin )    // if file could not be opened
        return;    // stop execution

    std::string line;               // string for storing lines, declared here to avoid redeclarations
    int i = 0;                      // counter for lines

    // while a line can be read
    while ( getline(fin, line) )
    {
        if (line[0] == '#' && line[0] != ' ')   // treat lines starting with '#' as comments
            continue;                           // continue loop

        std::string currentNum = "";    // current number that is being processed
        int j = 0;                      // current column in line i

        for (unsigned int k = 1; k < line.size(); k++) // for every character in the 'line' std::string
        {
            if ( line[k] == ' ' || line[k] == '\n' || k == line.size() - 1) // if last character
            {
                this->elements[i][j] = currentNum;  // assign the current number to the current element
                                                    // note that mpfr has a constructor which uses strings
                                                    // and that is being used here
                                                    //
                                                    // type(elements[i][j]) is mpreal
                                                    // type(currentNum) is std::string

                currentNum = ""; // empty the current number string
                j++;             // increment element counter
            }
            else
                currentNum += line[k];  // add another character to currentNum std::string
                                        // this should be an alphanum() character
        }
        i++;    // increment line number
    }

    fin.close();
}

// swap elements from matrix
void GmpMatrix::swapElement(int i, int j, int i2, int j2)
{
    std::swap(this->elements[i][j], this->elements[i2][j2]);
}

// set element in matrix
void GmpMatrix::setElement(int i, int j, mpfr::mpreal e)
{
    this->elements[i][j] = e;
}

// copy this into matrix b
void GmpMatrix::copyMatrix(GmpMatrix& B)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            this->elements[i][j] = B.getElement(i, j);
}

// print matrix to console
void GmpMatrix::printToConsole()
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
            std::cout << "X[" << i << "][" << j << "]: " << this->elements[i][j] << '\n';
            //mpfr_printf("X[%d][%d]: %.64Rf\n", i, j, this->elements[i][j].mpfr_ptr());
    }
    std::cout << '\n';
}

// multiply two matrices
bool GmpMatrix::multiplyMatrix(GmpMatrix& src1, GmpMatrix& src2, GmpMatrix& dest)
{
    // number of columns of first matrix must equal
    // the number of lines of the second matrix
    if ( src1.getM() != src2.getN() )
        return false;

    // matrix is dynamically allocated, so it has to be constructed properly from the start
    // changing the parameters is not helpful and reconstructing it is too much trouble
    // TODO: multiplication which returns GmpMatrix type?
    if ( dest.getN() != src1.getN() || dest.getM() != src2.getM() )
        return false;

    int srcCols = src2.getM();
    int srcRows = src1.getN();

    for (int i = 0; i < srcRows; ++i) // for every row of src1
    {

        for (int j = 0; j < srcCols; ++j) // for every column of src2
        {
            mpfr::mpreal res = 0;
            for (int k = 0; k < srcCols; ++k) // for every row of src2
            {
                res = res + src1.getElement(i, k) * src2.getElement(k, j);
            }
            dest.setElement(i, j, res);
        }
    }

    return true;
}

/*
    Overwrites current matrix with its inverse.
    This is done using the Gauss-Jordan method.
*/
mpfr::mpreal GmpMatrix::inverse_mpreal()
{
    GmpMatrix b(n, m); // create matrix b

    // initialize with 0
    // TODO: optimize?
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            b.setElement(i, j, 0.0);

    mpfr::mpreal det = 1.0; // initial determinant is 1.0
    mpfr::mpreal amax, t;

    int i, imax, j, jmax, k;

    std::vector<int> ipivot(n);
    std::vector<int> jpivot(n);
    std::vector<int> npivot(n);

    npivot.assign(n, 0);

    imax = 0; jmax = 0;

    for (k = 0; k < n; ++k)
    {
        amax = 0.0;

        for (int i = 0; i < n; ++i)
        {
            if (npivot[i] != 1)
                for (j = 0; j < n; ++j)
                    if ( npivot[j] == 0 )
                    {
                        mpfr::mpreal absres = mpfr::abs(this->elements[i][k]);
                        if ( amax < absres )
                        {
                            amax = absres;
                            imax = i;
                            jmax = j;
                        }
                    }
                    else if ( npivot[j] > 1 )
                    {
                        printf("Singular matrix! Not invertible.\n");
                        return 1;
                    }
        }

        ipivot[k] = imax;
        jpivot[k] = jmax;
        npivot[jmax]++;

        if ( amax == 0.0 )
        {
            printf("Singular matrix! Not invertible.\n");
            return 2;
        }

        if ( imax != jmax )
        {
            det *= -1;
            for (j = 0; j < n; ++j)
                this->swapElement(imax, j, jmax, j);
            for (j = 0; j < m; ++j)
                b.swapElement(imax, j, jmax, j);
        }

        det *= this->elements[jmax][jmax];

        t = 1 / this->elements[jmax][jmax];
        this->elements[jmax][jmax] = 1;

        for (j = 0; j < n; ++j)
            this->elements[jmax][j] *= t;

        for (j = 0; j < m; ++j)
            b.setElement(jmax, j, b.getElement(jmax, j) * t);

        for (i = 0; i < n; ++i)
        {
            if ( i != jmax )
            {
                t = this->elements[i][jmax];
                this->elements[i][jmax] = 0.0;

                mpfr::mpreal aux;

                for (j = 0; j < n; ++j)
                    this->elements[i][j] -= this->elements[jmax][j] * t;

                for (j = 0; j < m; ++j)
                    b.setElement( i, j, b.getElement(i, j) - b.getElement(jmax, j) * t);
            }
        }
    }

    for (k = n - 1; k >= 0; --k)
        if ( ipivot[k] != jpivot[k] )
            for (i = 0; i < n; ++i)
                this->swapElement(i, ipivot[k], i, jpivot[k]);

    return det;
}
