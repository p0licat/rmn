//
// Created by sanda on 1/23/2016.
//

#ifndef BIGNUMBER_BNMATRIX_H
#define BIGNUMBER_BNMATRIX_H

#define NORMAL_BN_LIST -1
#define COMPACT_BN_LIST -2
//#define PRECISSION_BN_LIST 3

#include "BigNumber.h"
#include <fstream>

using namespace std;

#define OUPUTMatrix "Matrix"

class BNMatrix {
private:
	int dim1;
	int dim2;
	BigNumber *p;

	bool BNMextendTO(BNMatrix &, const int, const int);

public:
	// Default Constructor
	BNMatrix();

	//  Overload Constructor
	BNMatrix(const int, const int);


	// Copy Constuctor
	BNMatrix(const BNMatrix&);

	// Destructor
	~BNMatrix();

	bool extendTO(const int, const int);

	int getDim1() const {
		return dim1;
	}

	void setDim1(int dim1) {
		BNMatrix::dim1 = dim1;
	}

	int getDim2() const {
		return dim2;
	}

	void setDim2(int dim2) {
		BNMatrix::dim2 = dim2;
	}

	BigNumber &getElement(const int, const int) const;

	void setElement(const int, const int, const BigNumber&);
	void setElement(const int, const int, const double);
	void setElement(const int, const int, const string);

	void copyValue_whenFits(BNMatrix &, const BNMatrix &);
	bool copyValue(BNMatrix &, const BNMatrix &);
	bool expandDim(BNMatrix &, int, int);

	BNMatrix matrixMultiply(const BNMatrix&, const BNMatrix&);


	void list();
	void compactList();
	void withPrecissionList(const int);

	bool toFile(const int);
	bool toFile(const char*, const int);

	bool fromFile(const char*);


	BigNumber getDeterminant() const;
	BNMatrix getCofactorMatrix(const int, const int) const;
	BNMatrix getInverse() const;



	//    Assignment operator
	BNMatrix& operator=(const BNMatrix&);

	// Overloaded Operators
	BNMatrix operator +(const BNMatrix&);
	BNMatrix operator *(const BNMatrix&);


};


#endif //BIGNUMBER_BNMATRIX_H
