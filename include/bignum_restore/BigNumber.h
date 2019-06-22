//
// Created by sanda on 12/20/2015.
//
#include <iostream>
#include <exception>
#include <math.h>
#include <string.h>
#include <sstream>
#include <algorithm>

#ifndef BIGNUMBER_Len
#define BIGNUMBER_LenInt 150
#define BIGNUMBER_LenFrac 150
#endif

#ifndef BIGNUMBER_maxLen
#define BIGNUMBER_maxLenInt 150
#define BIGNUMBER_maxLenFrac 200
#endif

#ifndef BIGNUMBER_BIGNUMBER_H
#define BIGNUMBER_BIGNUMBER_H

class BigNumber {

private:
	int *pFrac;
	int *pInt;

	int lenInt, lenFrac;
	bool sign = true; // TRUE - positive (+) & FALSE - negative (-)

	bool BNisZero(const BigNumber&);
	bool BNisOne(const BigNumber&);

	long BNgetDigitSum(const BigNumber&);


	int* setIntwithDouble(double, int);
	int* setFracwithDouble(double, int);

	void init(double, const int, const int);
	void makeZero();

	bool isGreater(const BigNumber&, const BigNumber&)const;
	bool isEqual(const BigNumber&, const BigNumber&)const;

	int add(int *, int *, int, int, int, bool);
	void intoadd(BigNumber&, const BigNumber&);
	void adding(BigNumber&, const BigNumber&, BigNumber&);
	int sub(int *, int *, int, int, int, bool);
	void intosub(BigNumber&, const BigNumber&);
	void subtracting(BigNumber&, const BigNumber&, BigNumber&);
	void mul(BigNumber&, BigNumber&, BigNumber&);
	void div(BigNumber&, BigNumber&, BigNumber&);

	void shiftRIGHT(BigNumber &, int);
	void shiftLEFT(BigNumber &, int);

	void swapValue(BigNumber&);
	bool copyValue(BigNumber &, const BigNumber &);
	void copyValue_whenFits(BigNumber &, const BigNumber &);
	void add1byPower(BigNumber &, int);
	bool expandINTpart(BigNumber &, int);
	bool expandFRApart(BigNumber &, int);
	int BNgetSignificatDecimals(const BigNumber &);
	int BNgetSignificatIntPart(const BigNumber &);

public:

	//    Default Constructor
	BigNumber();

	bool isSign() const;
	void setSign(const bool);

	//    Overload Constructor
	BigNumber(const int, const int);

	BigNumber(double);
	BigNumber(double, const int, const int);

	//    Copy Constructor
	BigNumber(const BigNumber &BN);

	//  Destructor
	~BigNumber();

	friend std::ostream& operator<< (std::ostream &out, BigNumber &BN);

	int getLenInt();
	int getLenFrac();
	//    void setLen(const int, const int);

	bool isZero();
	bool isOne();

	//    Assignment operator
	BigNumber& operator=(const BigNumber&);
	BigNumber& operator=(const double);
	BigNumber& operator=(const std::string);

	//  Comparison overloaded operators
	bool operator == (const BigNumber&) const;
	bool operator > (const BigNumber&) const;
	bool operator < (const BigNumber&) const;
	bool operator >= (const BigNumber&) const;
	bool operator <= (const BigNumber&) const;
	bool operator != (const BigNumber&) const;

	//  Mathematical overloaded operators
	BigNumber operator +(const BigNumber&);
	BigNumber operator -(const BigNumber&);
	BigNumber operator >> (const int); // shift to left (multiplication by (ten * @param))
	BigNumber operator <<(const int); // shift to left (multiplication by (ten * @param))
	BigNumber operator *(BigNumber&);
	BigNumber operator /(BigNumber&);


	BigNumber operator +=(const BigNumber&);
	BigNumber operator -=(const BigNumber&);
	BigNumber operator >>=(const int);
	BigNumber operator <<=(const int);
	BigNumber operator *=(BigNumber&);
	BigNumber operator /=(BigNumber&);

	int getSignificatDecimals();
	int getSignificatIntPart();
	int getFistZerousInDecimal();

	void compactList();
	void withPrecissionList(const int);

};


#endif //BIGNUMBER_BIGNUMBER_H
