//
// Created by sanda on 1/23/2016.
//

#include "BNMatrix.h"

BNMatrix::BNMatrix() {
	dim1 = 0;
	dim2 = 0;
	p = NULL;
}



BNMatrix::BNMatrix(const int m, const int n) {
	dim1 = m;
	dim2 = n;
	p = new BigNumber[dim1*dim2];

	for (int i = 0; i < dim1; ++i) {
		for (int j = 0; j < dim2; ++j) {
			p[i*dim2 + j] = "0.0";
		}
	}
}

BNMatrix::BNMatrix(const BNMatrix &BNM) {
	dim1 = BNM.dim1;
	dim2 = BNM.dim2;
	//    cout <<" copy constructor!... for "<<BNM.p[1]<<endl;

	p = new BigNumber[dim1*dim2];
	for (int i = 0; i < dim1; ++i) {
		for (int j = 0; j < dim2; ++j) {
			p[i*dim2 + j] = BNM.p[i*dim2 + j];
		}
	}
	//    std::copy(BNM.p,BNM.p+(dim1*dim2), p);
}

BNMatrix::~BNMatrix() {
	// clean up allocated memory
	delete[] p;
	p = NULL;
}


bool BNMatrix::BNMextendTO(BNMatrix &BNM, const int m, const int n) {
	if (m<BNM.dim1 || n<BNM.dim2) return false;

	BigNumber *np = new BigNumber[m*n];

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			np[i*n + j] = "0.0";
		}
	}

	for (int i = 0; i < BNM.dim1; ++i) {
		for (int j = 0; j < BNM.dim2; ++j) {
			np[i*n + j] = BNM.p[i*BNM.dim2 + j];
		}
	}

	delete[] BNM.p;
	BNM.p = np;


	BNM.dim1 = m;
	BNM.dim2 = n;

	return true;

}

bool BNMatrix::extendTO(const int m, const int n) {
	return BNMextendTO(*this, m, n);
}



BigNumber &BNMatrix::getElement(const int i, const int j) const {
	return p[(i - 1)*dim2 + j - 1];
}

void BNMatrix::setElement(const int i, const int j, const BigNumber &BN) {
	p[(i - 1)*dim2 + j - 1] = BN;

}

void BNMatrix::setElement(const int i, const int j, const double d) {
	p[(i - 1)*dim2 + j - 1] = d;
}

void BNMatrix::setElement(const int i, const int j, const string BNstring) {
	p[(i - 1)*dim2 + j - 1] = BNstring;
}

bool BNMatrix::toFile(const int how) {
	return toFile(OUPUTMatrix, how);
}

bool BNMatrix::toFile(const char *filename, const int how) {
	ofstream outFile;

	//    cout << filename << endl;
	outFile.open(filename, ofstream::out);
	if (outFile.is_open()) {
		//        cout << "file is opened!" << endl;
		//        outFile << "ceva\n";
		for (int i = 0; i < dim1; ++i) {
			for (int j = 0; j < dim2; ++j) {
				//                outFile << " (" <<i << "," <<j <<") " ;
				if (how == -1) { outFile << p[i * dim2 + j]; }
				else {
					streambuf *coutbuf = cout.rdbuf();
					cout.rdbuf(outFile.rdbuf());
					if (how == -2) p[i * dim2 + j].compactList();
					if (how > 0) p[i * dim2 + j].withPrecissionList(how);
					cout.rdbuf(coutbuf);
				}
				outFile << " ";
			}
			outFile << endl;
		}
		outFile.close();
		return true;
	}

	cerr << "Unable to open the file!!" << endl;
	return false;
}

bool BNMatrix::fromFile(const char *filename) {
	ifstream inFile;
	std::string output;
	int d1, d2, i = 0;

	inFile.clear();

	//    cout << "Start from: " << filename<< endl;
	inFile.open(filename, ios::in);

	if (inFile.is_open()) {
		getline(inFile, output);
		istringstream iss(output);
		iss >> d1 >> d2;
		//        cout <<"dim: ("<<d1<<","<<d2<<")"<<endl;
		BNMextendTO(*this, d1, d2);
		while (!inFile.eof()) {
			getline(inFile, output);
			istringstream iss(output);
			string word;
			while (iss >> word) {
				p[i] = word;
				//                cout << "i:["<<i <<"]"<<word<<endl;
				i++;
			}

		}
		inFile.close();
		return true;
	}

	cerr << "Unable to open the file: " << filename << endl;
	return false;
}

void BNMatrix::list() {
	for (int i = 0; i < dim1; ++i) {
		for (int j = 0; j < dim2; ++j) {
			cout << "    " << p[i * dim2 + j];
		}
		cout << endl;
	}
}

void BNMatrix::compactList() {
	for (int i = 0; i < dim1; ++i) {
		for (int j = 0; j < dim2; ++j) {
			cout << "    ";
			p[i * dim2 + j].compactList();
		}
		cout << endl;
	}
}

void BNMatrix::withPrecissionList(const int precision) {
	for (int i = 0; i < dim1; ++i) {
		for (int j = 0; j < dim2; ++j) {
			cout << "    " << "(" << i * dim2 + j << ")";//i:["<<i * dim2 + j<<"]";
			p[i * dim2 + j].withPrecissionList(precision);
		}
		cout << endl;
	}
}

void BNMatrix::copyValue_whenFits(BNMatrix &BNM, const BNMatrix &oBNM) {
	for (int i = 0; i < BNM.dim1; ++i) {
		for (int j = 0; j < BNM.dim2; ++j) {
			BNM.p[i*BNM.dim2 + j] = "0.0";
		}
	}
	for (int i = 0; i < oBNM.dim1; ++i) {
		for (int j = 0; j < oBNM.dim2; ++j) {
			BNM.p[i*oBNM.dim2 + j] = oBNM.p[i*oBNM.dim2 + j];
		}
	}
}

// internal modification of a BNMatrix
bool BNMatrix::expandDim(BNMatrix &BNM, int m, int n) {
	BigNumber *np = new BigNumber[m*n];

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			np[i*n + j] = "0.0";
		}
	}

	delete[] BNM.p;
	BNM.p = np;


	BNM.dim1 = m;
	BNM.dim2 = n;

	return true;
}

bool BNMatrix::copyValue(BNMatrix &BNM, const BNMatrix &oBNM) {
	bool expanded = true;
	if (BNM.dim1<oBNM.dim1 || BNM.dim2<oBNM.dim2) expanded = expandDim(BNM, oBNM.dim1, oBNM.dim2);
	if (expanded) {
		copyValue_whenFits(BNM, oBNM);
		return true;
	}
	else return false;
}


BNMatrix &BNMatrix::operator=(const BNMatrix &oBNM) {
	if (this != &oBNM) {
		//        dim1 = oBNM.dim1;
		//        dim2 = oBNM.dim2;

		copyValue(*this, oBNM);
	}
	return *this;
}


BNMatrix BNMatrix::operator +(const BNMatrix &oBNM) {
	if ((dim1 == oBNM.dim1) && (dim2 == oBNM.dim2)) {
		//        cout <<"they are"<< endl;
		BNMatrix rez(dim1, dim2);
		for (int i = 0; i < dim1; ++i) {
			for (int j = 0; j < dim2; ++j) {
				rez.p[i * dim2 + j] = p[i * dim2 + j] + oBNM.p[i * dim2 + j];
				//                cout <<"("<<i<<","<<j<<"):"<<p[i * dim2 + j]<<"+"<<oBNM.p[i * dim2 + j]<<"="<<rez.p[i * dim2 + j]<<endl;
			}
		}
		return rez;
	}
	BNMatrix rez(0, 0);
	return rez;
}


BNMatrix BNMatrix::operator*(const BNMatrix &oBNM) {
	if (this->dim1 != oBNM.dim2 || this->dim2 != oBNM.dim1)
	{
		BNMatrix rez;
		return rez;
	}
	else { return BNMatrix::matrixMultiply(*this, oBNM); }

	BNMatrix rez(this->dim1, oBNM.dim2);

	// parcurgem matricea rezultat
	for (int i = 0; i < this->dim1; ++i)
	{
		for (int j = 0; j < oBNM.dim2; ++j)
		{
			// fiecare element din rezultat se afla dintr-o suma a
			// elementelor celorlalte matrici
			BigNumber currentSum(0, 0); // facem suma corespunzatoare

			for (int k = 0; k < this->dim2; ++k) // parcurgem la this coloanele si la oBNM liniile
			{
				currentSum += this->p[i * this->dim2 + k] * oBNM.p[k * this->dim1 + j];
			}
			rez.p[i * this->dim1 + j] = currentSum;
		}
	}

	return rez;
}

BNMatrix BNMatrix::matrixMultiply(const BNMatrix& BNM, const BNMatrix& oBNM)
{
	// pentru inmultirea matematica a doua matrici,
	// A(n, m) si B(x, y)
	// trebuie sa avem m = x si n = y
	// matricea rezultat va avea marimea Rez(n, y)
	// pentru ca n = y inseamna ca Rez este patratica, sau Rez(n, n)

	if (BNM.dim1 != oBNM.dim2 || BNM.dim2 != oBNM.dim1)
	{
		BNMatrix rez;
		return rez;
	}

	BNMatrix rez(BNM.dim1, oBNM.dim2);

	// parcurgem matricea rezultat
	for (int i = 0; i < BNM.dim1; ++i)
	{
		for (int j = 0; j < oBNM.dim2; ++j)
		{
			// fiecare element din rezultat se afla dintr-o suma a
			// elementelor celorlalte matrici
			BigNumber currentSum(0, 0); // facem suma corespunzatoare

			for (int k = 0; k < BNM.dim2; ++k) // parcurgem la this coloanele si la oBNM liniile
			{
				currentSum += BNM.p[i * BNM.dim2 + k] * oBNM.p[k * BNM.dim1 + j];
			}
			rez.p[i * BNM.dim1 + j] = currentSum;
		}
	}

	return rez;
}



BigNumber BNMatrix::getDeterminant() const {
	if (dim1 == dim2) {
		BigNumber deter, a1, a2, temp;
		deter = 0;
		if (dim1 >= 3) {
			for (int i = 1; i <= dim1; ++i) {
				a1 = 1;
				a2 = 1;
				for (int j = 1; j <= dim1; ++j) {
					temp = getElement((i + j - 2) % dim1 + 1, j);
					a1 *= temp; //getElement((i+j-2)%dim1 + 1, j);
								//                    cout <<"a1*=element(" << (i+j-2)%dim1 + 1 << ", " << j << ")="<< a1<< endl;
					temp = getElement((i + j - 2) % dim1 + 1, dim1 - j + 1);
					a2 *= temp; //getElement((i+j-2)%dim1 + 1, dim1-j+1);
								//                    cout << endl<<"a1*=(" << (i+j-2)%dim1 + 1 << ", " << j << ")="<< a1<<"  a2*=(" << (i+j-2)%dim1 + 1 << ", " <<  dim1-j+1 << "):"<<temp<<"="<< a2;
				}
				temp = a1 - a2;
				//                cout <<"---- deter="<<deter;
				deter += temp;
				//                cout << endl<<"deter+=(a1-a2)="<<a1<<"-"<<a2<<"="<< temp <<"; deter="<<deter;
			}
		}
		else {
			if (dim1 == 2) {
				deter = getElement(1, 1)*getElement(2, 2) - getElement(1, 2)*getElement(2, 1);
			}
			else {
				deter = getElement(1, 1);
			}
		}
		return deter;
	}
	else {
		cout << "I need a square matrix!";
	}
	BigNumber nullval(0, 0); //TODO: problem
	return nullval;
}

BNMatrix BNMatrix::getCofactorMatrix(const int x, const int y) const {
	BNMatrix rez(dim1 - 1, dim2 - 1);

	int k, l;
	k = 0;
	for (int i = 0; i < dim1; i++) {
		l = 0;
		if (i != (x - 1)) {
			for (int j = 0; j < dim2; j++) {
				if (j != (y - 1)) {
					rez.p[k*rez.dim2 + l] = p[i*dim2 + j];
					l++;
				}
			}
			k++;
		}
	}
	return rez;
}

BNMatrix BNMatrix::getInverse() const {
	BNMatrix rez(dim1, dim2), C(dim1 - 1, dim2 - 1);
	BigNumber det, det1, power;

	det = getDeterminant();
	// cout << "Odet:"; det.listBN(); cout <<endl <<endl;

	for (int i = 0; i < dim1; ++i) {
		for (int j = 0; j < dim2; ++j) {
			C = getCofactorMatrix(i + 1, j + 1);
			//            C.withPrecissionList(2);
			det1 = C.getDeterminant();
			// cout << "det:"; det1.listBN(); cout <<endl <<endl;
			power = pow(-1, i + j + 2);
			rez.p[j*dim2 + i] = (power*det1) / det;
			// cout << "elem: ("<<i+1<<","<<j+1<<")"; rez.p[i*dim2+j].listBN();
			// cout << endl << endl;


		}
	}
	return rez;
}

