//
// Created by Sanda Dragos on 12/20/2015.
//


#include "BigNumber.h"

//default constructor
BigNumber::BigNumber() {
	pInt = NULL;
	pFrac = NULL;
	lenInt = 0;
	lenFrac = 0;
}

void BigNumber::makeZero() {
	for (int i = 0; i < lenInt; ++i) {
		pInt[i] = 0;
	}
	for (int j = 0; j < lenFrac; ++j) {
		pFrac[j] = 0;
	}
}

BigNumber::BigNumber(const int lInt, const int lFrac) {
	lenInt = lInt;
	lenFrac = lFrac;
	pInt = new int[lenInt];
	pFrac = new int[lenFrac];
	makeZero();
}


// create BigNumber out of a double
int* BigNumber::setIntwithDouble(double d, int lInt) {

	int *pI = new int[lInt];
	int NOdigits = (int)((d != 0) ? log10(d) + 1 : 1);

	if (NOdigits>lInt) {// the number is bigger than its INT representation -> throw exception0
		throw "The number is bigger than its INT representation!";
	}
	else {
		for (int i = lInt - 1; i > lInt - 1 - NOdigits; i--) {
			int num = (int)(d / pow(10, lInt - i - 1));
			pI[i] = num % 10;
		}
		for (int j = lInt - 1 - NOdigits; j >= 0; j--) { pI[j] = 0; }
	}
	return pI;

}

int* BigNumber::setFracwithDouble(double d, int lFrac) {
	int *pF = new int[lFrac];
	std::ostringstream strs;
	strs << d;
	std::string str = strs.str();
	int NOdigits = str.size() - 2;

	if (NOdigits>lFrac) {// the number is bigger than its FRAC representation -> throw exception
		throw "The number is bigger than its FRAC representation! "; //+ std::to_string(d);
	}
	else {
		for (int i = 0; i < NOdigits; ++i) {
			int num = (int)(d*pow(10, i + 1));
			pF[i] = ((int)num) % 10;
		}
		for (int j = NOdigits; j < lFrac; j++) { pF[j] = 0; }
	}
	return pF;
}

void BigNumber::init(double d, const int lInt, const int lFrac) {
	if (d<0) {
		setSign(false);
		d = -d;
	}
	else setSign(true);

	try {
		lenInt = (lInt <= BIGNUMBER_maxLenInt) ? lInt : BIGNUMBER_maxLenInt;
		lenFrac = (lFrac <= BIGNUMBER_maxLenFrac) ? lFrac : BIGNUMBER_maxLenFrac;
		pInt = setIntwithDouble((long)d, lenInt);
		pFrac = setFracwithDouble(d - (long)d, lenFrac);
	}
	catch (const char* msg) {
		std::cerr << msg << std::endl;
		exit(1);
	}
}

BigNumber::BigNumber(double d) {
	init(d, BIGNUMBER_LenInt, BIGNUMBER_LenFrac);
}

BigNumber::BigNumber(double d, const int lInt, const int lFrac) {
	init(d, lInt, lFrac);
}

BigNumber::~BigNumber() {
	// clean up allocated memory
	delete pInt;
	pInt = NULL;
	delete pFrac;
	pFrac = NULL;
}

std::ostream &operator<<(std::ostream &out, BigNumber &BN) {
	// Since operator<< is a friend of the Point class, we can access
	// Point's members directly.
	out.clear();
	if (!BN.isSign()) out << "-";
	for (int i = 0; i < BN.lenInt; ++i) { out << BN.pInt[i]; }
	out << ".";
	for (int i = 0; i < BN.lenFrac; ++i) { out << BN.pFrac[i]; }

	return out;
}
// copy constructor
// create BigNumber out of a BigNumber
BigNumber::BigNumber(const BigNumber &BN) {
	sign = BN.sign;
	lenInt = BN.lenInt;
	lenFrac = BN.lenFrac;

	pInt = new int[lenInt];
	pFrac = new int[lenFrac];

	//for (int i = 0; i < lenInt; ++i) {
	//pInt[i]=BN.pInt[i];
	//}
	std::copy(BN.pInt, BN.pInt + BN.lenInt, pInt);
	//for (int j = 0; j < lenFrac; ++j) {
	//pFrac[j]=BN.pFrac[j];
	//}
	std::copy(BN.pFrac, BN.pFrac + BN.lenFrac, pFrac);

}

// setters & getters
bool BigNumber::isSign() const {
	return sign;
}

void BigNumber::setSign(const bool s) {
	sign = s;
}

int BigNumber::getLenInt() {
	return lenInt;
}

int BigNumber::getLenFrac() {
	return lenFrac;
}

// check whether a BigNumber zero or one
bool BigNumber::BNisZero(const BigNumber &oBN) {
	for (int i = 0; i < oBN.lenInt; ++i) {
		if (oBN.pInt[i] != 0)  return false;
	}
	for (int j = 0; j < oBN.lenFrac; ++j) {
		if (oBN.pFrac[j] != 0) return false;
	}
	return true;
}

bool BigNumber::isZero() {
	return BNisZero(*this);
}

bool BigNumber::BNisOne(const BigNumber &oBN) {
	if (oBN.pInt[oBN.lenInt - 1] != 1) return false;
	else {
		for (int i = 0; i < oBN.lenInt - 1; ++i) {
			if (oBN.pInt[i] != 0) return false;
		}
		for (int j = 0; j < oBN.lenFrac; ++j) {
			if (oBN.pFrac[j] != 0) return false;
		}
	}
	return true;
}

bool BigNumber::isOne() {
	return BNisOne(*this);
}

// internal modification of a BigNumber
bool BigNumber::expandINTpart(BigNumber &BN, int with) {
	if ((BN.lenInt + with) <= BIGNUMBER_maxLenInt) {
		int *np = new int[BN.lenInt + with];
		for (int i = BN.lenInt + with - 1; i >= with; i--) np[i] = BN.pInt[i - with];
		for (int i = with - 1; i >= 0; i--) np[i] = 0;
		delete[] BN.pInt;
		BN.pInt = np;
		BN.lenInt = BN.lenInt + with;
		return true;
	}
	else {
		//        std::cerr << "Number wants to go over the limit set for the integer part!"<< std::endl;
		return false;
	}
}

bool BigNumber::expandFRApart(BigNumber &BN, int with) {
	if ((BN.lenFrac + with) <= BIGNUMBER_maxLenFrac) {
		int *np = new int[BN.lenFrac + with];
		for (int i = 0; i < BN.lenFrac; i++) np[i] = BN.pFrac[i];
		for (int i = BN.lenFrac; i < BN.lenFrac + with; i++) np[i] = 0;
		delete[] BN.pFrac; // TODO: problem, corrupted heap
		BN.pFrac = np;
		BN.lenFrac = BN.lenFrac + with;
		return true;
	}
	else {
		//        std::cerr << "Number wants to go over the limit set for the fractional part!"<< std::endl;
		return false;
		//exit(101);
	}
}

void BigNumber::copyValue_whenFits(BigNumber &BN, const BigNumber &oBN) {
	for (int i = 0; i < BN.lenInt - oBN.lenInt; ++i) {
		BN.pInt[i] = 0;
	}
	for (int i = BN.lenInt - oBN.lenInt; i < BN.lenInt; ++i) {
		BN.pInt[i] = oBN.pInt[i - BN.lenInt + oBN.lenInt];
	}
	for (int j = 0; j < oBN.lenFrac; ++j) {
		BN.pFrac[j] = oBN.pFrac[j];
	}
	for (int j = oBN.lenFrac; j < BN.lenFrac; ++j) {
		BN.pFrac[j] = 0;
	}
}

bool BigNumber::copyValue(BigNumber &BN, const BigNumber &oBN) {
	bool expanded = true;
	if (BN.lenInt<oBN.lenInt) expanded = expandINTpart(BN, oBN.lenInt - BN.lenInt);
	if (BN.lenFrac<oBN.lenFrac && expanded) expanded = expandFRApart(BN, oBN.lenFrac - BN.lenFrac);
	if (expanded) {
		copyValue_whenFits(BN, oBN);
		return true;
	}
	else return false;
}


BigNumber &BigNumber::operator=(const BigNumber &oBN) {
	if (this != &oBN) {
		sign = oBN.sign;
		if ((oBN.lenInt + lenInt > 0) && (oBN.lenFrac + lenFrac > 0)) {
			copyValue(*this, oBN);
		}
	}
	return *this;
}

BigNumber &BigNumber::operator=(const double d) {
	if ((lenInt>0) && (lenFrac>0)) init(d, lenInt, lenFrac);
	else init(d, BIGNUMBER_LenInt, BIGNUMBER_LenFrac);
	return *this;
}

BigNumber &BigNumber::operator=(const std::string num) {
	lenInt = BIGNUMBER_LenInt;
	lenFrac = BIGNUMBER_LenFrac;
	pInt = new int[lenInt];
	pFrac = new int[lenFrac];
	makeZero();

	size_t found = num.find_first_of("+-");
	int j = (found != std::string::npos) ? found + 1 : 0;
	sign = (j != 0) ? false : true;

	found = num.find_first_of('.');
	int k = (found != std::string::npos) ? found : num.length();
	//    std::cout <<"len-init: "<<k<<std::endl;
	int ii = lenInt - 1;
	for (int i = k - 1; i >= j; --i) {
		pInt[ii--] = num[i] - '0';
		//        std::cout <<"i-init: "<<ii<<" digit:"<< num[i]-'0'<<std::endl;
	}
	ii = 0;
	for (unsigned int j = k + 1; j < num.length(); ++j) {
		pFrac[ii++] = num[j] - '0';
	}

	return *this;
}

bool BigNumber::isEqual(const BigNumber &BN, const BigNumber &BN2)const {
	int beg1 = 0, beg2 = 0;
	if (BN.lenInt>BN2.lenInt) {
		for (int i = 0; i < BN.lenInt - BN2.lenInt; ++i) {
			if (BN.pInt[i] != 0) return true;
		}
		beg1 = BN.lenInt - BN2.lenInt;
	}
	else if (BN.lenInt<BN2.lenInt) {
		for (int i = 0; i < BN2.lenInt - BN.lenInt; ++i) {
			if (BN2.pInt[i] != 0) return false;
		}
		beg2 = BN2.lenInt - BN.lenInt;
	}
	for (int i = 0; i < std::min(BN.lenInt, BN2.lenInt); ++i) {
		if (BN.pInt[i + beg1] != BN2.pInt[i + beg2]) return false;
	}
	for (int j = 0; j < std::min(BN.lenFrac, BN2.lenFrac); ++j) {
		if (BN.pFrac[j] != BN2.pFrac[j]) return false;
	}
	if (BN.lenFrac>BN2.lenFrac) {
		for (int j = BN2.lenFrac; j < BN.lenFrac; ++j) {
			if (BN.pFrac[j] != 0) return false;
		}
	}
	else if (BN.lenFrac<BN2.lenFrac) {
		for (int j = BN.lenFrac; j < BN2.lenFrac; ++j) {
			if (BN2.pFrac[j] != 0) return false;
		}
	}
	return true;
}

bool BigNumber::operator==(const BigNumber &oBN) const {
	if (sign != oBN.sign)  return false;
	else return isEqual(*this, oBN);

}

bool BigNumber::isGreater(const BigNumber &BN, const BigNumber &BN2) const {
	int beg1 = 0, beg2 = 0;
	if (BN.lenInt>BN2.lenInt) {
		for (int i = 0; i < BN.lenInt - BN2.lenInt; ++i) {
			if (BN.pInt[i] != 0) return true;
		}
		beg1 = BN.lenInt - BN2.lenInt;
	}
	else if (BN.lenInt<BN2.lenInt) {
		for (int i = 0; i < BN2.lenInt - BN.lenInt; ++i) {
			if (BN2.pInt[i] != 0) return false;
		}
		beg2 = BN2.lenInt - BN.lenInt;
	}
	for (int i = 0; i < std::min(BN.lenInt, BN2.lenInt); ++i) {
		if (BN.pInt[i + beg1]>BN2.pInt[i + beg2]) return true;
		if (BN.pInt[i + beg1]<BN2.pInt[i + beg2]) return false;
	}
	for (int j = 0; j < std::min(BN.lenFrac, BN2.lenFrac); ++j) {
		if (BN.pFrac[j]>BN2.pFrac[j]) return true;
		if (BN.pFrac[j]<BN2.pFrac[j]) return false;
	}
	return false;
}

bool BigNumber::operator >(const BigNumber &oBN) const {
	if (sign && !oBN.sign) return true;
	else return isGreater(*this, oBN);
}

bool BigNumber::operator <(const BigNumber &oBN) const {
	if (!sign && oBN.sign) return true;
	else {
		if (!isGreater(*this, oBN) && !isEqual(*this, oBN)) return true;
	}
	return false;
}

bool BigNumber::operator >=(const BigNumber &oBN) const {
	return (isGreater(*this, oBN) || isEqual(*this, oBN));
}

bool BigNumber::operator <=(const BigNumber &oBN) const {
	return (!isGreater(*this, oBN));
}

bool BigNumber::operator !=(const BigNumber &oBN) const {
	return (!isEqual(*this, oBN));
}

void BigNumber::swapValue(BigNumber &oBN) {
	int aux;

	if ((lenInt == oBN.lenInt) && (lenFrac == oBN.lenFrac)) {
		for (int i = 0; i < lenInt; ++i) {
			aux = oBN.pInt[i];
			oBN.pInt[i] = pInt[i];
			pInt[i] = aux;
		}
		for (int j = 0; j < lenFrac; ++j) {
			aux = oBN.pFrac[j];
			oBN.pFrac[j] = pFrac[j];
			pFrac[j] = aux;
		}
	}
	else { std::cerr << "Not the same size!" << std::endl; }
} //??

  //    we presume that len > olen ALWAYS!!!
int BigNumber::add(int *p, int *o, int len, int olen, int transport, bool frac) {
	int sum = 0;

	if (frac) {
		transport = 0;
		for (int i = olen - 1; i >= 0; i--) {
			sum = p[i] + o[i] + transport;
			p[i] = sum % 10;
			transport = sum / 10;
		}
	}
	else {
		for (int i = olen - 1; i >= 0; i--) {
			sum = p[i - olen + len] + o[i] + transport;
			p[i - olen + len] = sum % 10;
			transport = sum / 10;
		}
		for (int j = len - olen - 1; j >= 0; j--) {
			sum = p[j] + transport;
			p[j] = sum % 10;
			transport = sum / 10;
		}
	}
	return transport;
}

void BigNumber::intoadd(BigNumber &rez, const BigNumber &oBN) {
	int transport;

	if (rez.lenInt<oBN.lenInt) expandINTpart(rez, oBN.lenInt - rez.lenInt);
	if (rez.lenFrac<oBN.lenFrac) expandFRApart(rez, oBN.lenFrac - rez.lenFrac);


	transport = add(rez.pFrac, oBN.pFrac, rez.lenFrac, oBN.lenFrac, 0, true);
	if (transport>10) {
		std::cerr << "!!! transport1 ERROR: " << transport << "; when adding" << std::endl;
		//        std::cout <<"!!! transport1 ERROR: " << transport << "; when adding"<<rez <<" and "<< std::endl;
	}
	else transport = add(rez.pInt, oBN.pInt, rez.lenInt, oBN.lenInt, transport, false);
	if (transport>0) {
		if (transport<10) { bool exp = expandINTpart(rez, 1); if (exp) rez.pInt[0] = transport; }
		else {
			std::cerr << "!!! transport2 ERROR: " << transport << "; when adding" << std::endl;
			//        std::cout <<"!!! transport2 ERROR: " << transport << "; when adding"<<rez <<" and "<< std::endl;
		}
	}
}

void BigNumber::adding(BigNumber &BN, const BigNumber &oBN, BigNumber &rez) {
	copyValue(rez, BN);
	rez.sign = BN.sign;
	intoadd(rez, oBN);
	//        std::cout <<"same sign="<<std::max(lenFrac,oBN.lenFrac)<<std::endl;

}

BigNumber BigNumber::operator +(const BigNumber &oBN) {
	BigNumber rez(std::max(lenInt, oBN.lenInt), std::max(lenFrac, oBN.lenFrac));

	if (BNisZero(*this)) rez = oBN;
	else if (BNisZero(oBN)) rez = *this;
	else
		//    std::cout <<"in add=("<<rez.lenInt<<", "<<rez.lenFrac<<")"<<std::endl;
		if (sign == oBN.sign) { // same sign
			adding(*this, oBN, rez);
		}
		else { // different sign => we substract
			subtracting(*this, oBN, rez);
		}
		return rez;
}

//    we presume that len > olen ALWAYS!!!
int BigNumber::sub(int *p, int *o, int len, int olen, int transport, bool frac) {
	int dif, step, borrow;

	if (frac) {
		borrow = 0;
		for (int i = olen - 1; i >= 0; i--) {
			dif = p[i] - o[i] - borrow;
			if (dif >= 0) {
				p[i] = dif;
				borrow = 0;
			}
			else {
				p[i] = dif + 10;
				borrow = 1;
			}
		}
	}
	else {
		borrow = transport;
		for (int i = olen - 1; i >= 0; i--) {
			dif = p[i - olen + len] - o[i] - borrow;
			if (dif >= 0) {
				p[i - olen + len] = dif;
				borrow = 0;
			}
			else {
				p[i - olen + len] = dif + 10;
				borrow = 1;
			}
		}
		for (int j = len - olen - 1; j >= 0; j--) {
			dif = p[j] - borrow;
			if (dif >= 0) {
				p[j] = dif;
				borrow = 0;
			}
			else {
				p[j] = dif + 10;
				borrow = 1;
			}
		}
	}

	//    step=(frac)? olen-1: len-olen-1;
	//    int lenDiff=(!frac && (len>olen))? len-olen :0;
	//
	////    std::cout <<"frac="<<frac<<"; step="<<step<<"; len="<<len<<"; olen="<<olen<<std::endl;
	//
	//    //    for frac is the difference between the vectors
	//    //    for NOT frac is the common part
	//    for (int i = len-1; i >step ; --i) {
	////        std::cout <<"first for i="<<i<<std::endl;
	//        borrow=(frac)?transport:o[i-lenDiff]+transport;
	//        if(p[i]>=borrow){
	//            dif=p[i]-borrow;
	//            p[i]=dif;
	//            transport=0;
	//        } else {
	//            dif=p[i]-borrow+10;
	//            p[i]=dif;
	//            transport=1;
	//        }
	//    }
	//
	//    //    for frac is the common part of the two the vectors
	//    //    for NOT frac is the difference between
	//    for (int i = step; i >=0 ; --i) {
	////        std::cout <<"second for i="<<i<<std::endl;
	//        borrow=(frac)?o[i]+transport:transport;
	//        if(p[i]>=borrow){
	//            dif=p[i]-borrow;
	//            p[i]=dif;
	//            transport=0;
	//        } else {
	//            dif=p[i]-borrow+10;
	//            p[i]=dif;
	//            transport=1;
	//        }
	//    }

	return borrow;
}

void BigNumber::intosub(BigNumber &rez, const BigNumber &oBN) {
	int transport;

	if (rez.lenInt<oBN.lenInt) expandINTpart(rez, oBN.lenInt - rez.lenInt);
	if (rez.lenFrac<oBN.lenFrac) expandFRApart(rez, oBN.lenFrac - rez.lenFrac);

	transport = sub(rez.pFrac, oBN.pFrac, rez.lenFrac, oBN.lenFrac, 0, true);
	transport = sub(rez.pInt, oBN.pInt, rez.lenInt, oBN.lenInt, transport, false);
	if (transport != 0) { bool exp = expandINTpart(rez, 1); if (exp) rez.pInt[0] = transport; }
}

void BigNumber::subtracting(BigNumber &BN, const BigNumber &oBN, BigNumber &rez) {

	if (isGreater(*this, oBN)) {
		rez.sign = BN.sign;
		copyValue(rez, *this);
		intosub(rez, oBN);

	}
	else {
		rez.sign = !BN.sign;
		copyValue(rez, oBN);

		intosub(rez, BN);
	}

}

BigNumber BigNumber::operator-(const BigNumber &oBN) {
	BigNumber rez(std::max(lenInt, oBN.lenInt), std::max(lenFrac, oBN.lenFrac));
	if (sign == oBN.sign) { // same sign
		subtracting(*this, oBN, rez);
	}
	else { // different signs => we add
		adding(*this, oBN, rez);
		rez.sign = !oBN.sign;
	}
	return rez;
}

BigNumber BigNumber::operator +=(const BigNumber &oBN) {
	BigNumber rez;
	if (BNisZero(*this)) rez = oBN;
	else if (BNisZero(oBN)) rez = *this;
	else
		if (sign == oBN.sign) {
			if (lenInt<oBN.lenInt) expandINTpart(*this, oBN.lenInt - lenInt);
			if (lenFrac<oBN.lenFrac) expandFRApart(*this, oBN.lenFrac - lenFrac);

			adding(*this, oBN, rez);

		}
		else {// different signs => we substract
			subtracting(*this, oBN, rez);
			//        sign=(*this>oBN)? oBN.sign : !oBN.sign;
		}
		copyValue(*this, rez);
		sign = rez.sign;
		return *this;

}

BigNumber BigNumber::operator-=(const BigNumber &oBN) {
	if (sign == oBN.sign) {
		sign = (*this>oBN) ? oBN.sign : !oBN.sign;

		intosub(*this, oBN);
	}
	else { // different signs => we add
		intoadd(*this, oBN);
		sign = !oBN.sign;
	}
	return *this;
}



void BigNumber::shiftRIGHT(BigNumber &BN, int n) {
	//    bool search=true;
	//    int i = BN.lenFrac-1;
	//    while((i >BN.lenFrac - n)&&search) {
	//        if(BN.pFrac[i]!=0) search=false;
	//        i--;
	//    }
	int i = BNgetSignificatDecimals(BN);
	if (i<n) expandFRApart(BN, n - i + 1);

	//    std::cout << "we can shift to right: " << i <<" m="<<m<<std::endl;

	for (int i = BN.lenFrac - 1 - n; i >= 0; i--) {
		//        std::cout << "frac[" << i+n <<"]=frac["<<i<<"]"<<std::endl;
		BN.pFrac[i + n] = BN.pFrac[i];
	}

	for (int j = 0; j <n; ++j) {
		//        std::cout << "frac[" << j <<"]=int["<<BN.lenInt-n+j<<"]"<<std::endl;
		BN.pFrac[j] = (BN.lenInt >= n - j) ? BN.pInt[BN.lenInt - n + j] : 0;
	}
	for (int k = BN.lenInt - 1; k >= 0; k--) {
		BN.pInt[k] = (k >= n) ? BN.pInt[k - n] : 0;
	}
}

// right shift of BigNumber with dim
BigNumber BigNumber::operator >> (const int dim) {
	BigNumber rez(lenInt, lenFrac);
	if (dim>0) {
		rez.sign = sign;
		copyValue(rez, *this);
		shiftRIGHT(rez, dim);
	}
	else copyValue(rez, *this);
	return rez;
}

BigNumber BigNumber::operator>>=(const int dim) {
	if (dim>0) {
		shiftRIGHT(*this, dim);
	}
	return *this;
}

void BigNumber::shiftLEFT(BigNumber &BN, int n) {
	//    bool search=true;
	//    int i = 0;
	//    while((i <BN.lenInt)&&search) {
	//        if(BN.pInt[i]!=0) search=false;
	//        i++;
	//    }
	int i = BNgetSignificatIntPart(BN);
	if (i<n) expandINTpart(BN, n - i + 1);

	for (int i = 0; i < BN.lenInt - n; ++i) {
		BN.pInt[i] = BN.pInt[i + n];
	}
	for (int j = BN.lenInt - n; j <BN.lenInt; ++j) {
		//        std::cout << "int[" << j <<"]=frac["<<j-BN.lenInt+n<<"]"<<std::endl;
		BN.pInt[j] = (j - BN.lenInt + n<BN.lenFrac) ? BN.pFrac[j - BN.lenInt + n] : 0;
	}

	for (int k = 0; k < BN.lenFrac; ++k) {
		BN.pFrac[k] = (k + n<BN.lenFrac) ? BN.pFrac[k + n] : 0;
	}
}

// left shift of BigNumber with dim
BigNumber BigNumber::operator<<(const int dim) {
	BigNumber rez(lenInt, lenFrac);
	if (dim>0) {
		rez.sign = sign;
		copyValue(rez, *this);
		shiftLEFT(rez, dim);
	}
	else copyValue(rez, *this);
	return rez;
}

BigNumber BigNumber::operator<<=(const int dim) {
	if (dim>0) {
		shiftLEFT(*this, dim);
	}
	return *this;
}

long BigNumber::BNgetDigitSum(const BigNumber &oBN) {
	long rez = 0;
	for (int i = 0; i < oBN.lenInt; ++i) { rez += oBN.pInt[i]; }
	for (int i = 0; i < oBN.lenFrac; ++i) { rez += oBN.pFrac[i]; }

	return rez;
}

void BigNumber::mul(BigNumber &BN1, BigNumber &BN2, BigNumber& rez) {
	if (rez.lenInt<BN1.lenInt) expandINTpart(rez, BN1.lenInt - rez.lenInt);
	if (rez.lenFrac<BN1.lenFrac) expandFRApart(rez, BN1.lenFrac - rez.lenFrac);

	//    std::cout << "BN1:"<< BN1 << "+ BN2:"<<BN2 <<"=rez:"<<rez << std::endl;

	//    INT
	for (int i = 0; i < BN2.lenInt; ++i) {
		//        rez<<=1;
		if (!BNisZero(rez)) { shiftLEFT(rez, 1); }
		for (int j = 1; j <= BN2.pInt[i]; ++j) {
			//            rez+=BN1;
			intoadd(rez, BN1);
			//            std::cout << "i:"<< j  <<"; rez:"<<rez<<"; BN1:"<<BN1 << std::endl;
		}
	}
	//    std::cout <<"dupaINT=rez:"<<rez << std::endl;
	BigNumber rezFrac(rez.lenInt, rez.lenFrac);
	rezFrac.sign = rez.sign;

	//    FRAC
	for (int i = BN2.lenFrac - 1; i >= 0; i--) {
		for (int j = 0; j < BN2.pFrac[i]; ++j) {
			////            rez+=BN1;
			intoadd(rezFrac, BN1);
		}
		////        rez>>=1;
		if (!BNisZero(rezFrac)) { shiftRIGHT(rezFrac, 1); }
	}
	//    std::cout <<"dupaINT=rezFrac:"<<rezFrac << std::endl;
	////            rez+=rezFrac;
	intoadd(rez, rezFrac);
	//   std::cout <<"FINAL=rez:"<<rez << std::endl;
}


BigNumber BigNumber::operator*(BigNumber &oBN) {
	BigNumber rez(std::max(lenInt, oBN.lenInt), std::max(lenFrac, oBN.lenFrac));

	if (isZero() || BNisZero(oBN)) return rez;

	if (isOne()) { copyValue(rez, oBN); return rez; }
	if (BNisOne(oBN)) { copyValue(rez, *this); return rez; }

	if (BNgetDigitSum(*this)>BNgetDigitSum(oBN)) {
		rez.sign = sign;
		mul(*this, oBN, rez);
		//        std::cout <<"mul(this, oBN, rez)"<<rez << std::endl;
	}
	else {
		rez.sign = oBN.sign;
		mul(oBN, *this, rez);
		//        std::cout <<"mul(oBN, this, rez)"<<rez << std::endl;

	}

	rez.sign = sign&&oBN.sign;
	return rez;

	//    std::cout <<"mul "<<"done"<< std::endl;}

}


BigNumber BigNumber::operator*=(BigNumber &oBN) {

	if (isZero() || BNisZero(oBN)) { makeZero(); return *this; }
	sign = (sign == oBN.sign) ? true : sign&&oBN.sign;
	if (isOne()) { swapValue(oBN); return *this; }
	if (BNisOne(oBN)) return *this;

	BigNumber rez(std::max(lenInt, oBN.lenInt), std::max(lenFrac, oBN.lenFrac));
	if (BNgetDigitSum(*this)>BNgetDigitSum(oBN)) mul(*this, oBN, rez);
	else mul(oBN, *this, rez);

	//    std::cout << "mul_rez:"<<*this<<"*"<<oBN<<"="<<rez<< std::endl;

	//    swapValue(rez);
	if (!copyValue(*this, rez)) std::cout << "COULD NOT COPPY!!!!!" << std::endl;
	return *this;
}



void BigNumber::add1byPower(BigNumber &BN, int power) {
	int j;
	bool toAdd, expand = true;

	if (power >= 0) { // in the INT part
		if (BN.lenInt - 1 - power < 0) expand = expandINTpart(BN, power - BN.lenInt + 1);
		if (expand) {
			//            std::cout << "==" << power << " len=" << BN.lenInt << std::endl;
			if (BN.pInt[BN.lenInt - 1 - power] == 9) {
				j = BN.lenInt - 1 - power;
				toAdd = true;
				while (toAdd) {
					if (j == 0) {
						expandINTpart(BN, 1);
						j = 1;
					}
					BN.pInt[j] = 0;
					j--;
					if (BN.pInt[j] != 9) {
						BN.pInt[j]++;
						toAdd = false;
					}
				}
			}
			else BN.pInt[BN.lenInt - 1 - power]++;
		}
	}
	else { // in the FRAC part
		   //            std::cout << "power=" << power << " len=" << BN.lenFrac << std::endl;
		if ((-power) > BN.lenFrac) expand = expandFRApart(BN, -power - BN.lenFrac);
		if (expand) {
			if (BN.pFrac[-power - 1] == 9) {
				j = -power - 1;
				toAdd = true;
				while (toAdd) {
					if (j == 0) {
						BN.pFrac[j] = 0;
						add1byPower(BN, 0);
						toAdd = false;
					}
					if (j >= 1) {
						BN.pInt[j] = 0;
						j--;

					}
					if (BN.pInt[j] != 9) {
						BN.pInt[j]++;
						toAdd = false;
					}
				}
			}
			else BN.pFrac[-power - 1]++;
		}
	}
}

void BigNumber::div(BigNumber &BN1, BigNumber &BN2, BigNumber& rez) {
	BigNumber aux(BN1);
	int power = 0;

	while (!BNisZero(aux)) {
		while ((isGreater(aux, BN2) || isEqual(aux, BN2)) && -power<BIGNUMBER_maxLenFrac) {

			//            aux-=BN2;
			intosub(aux, BN2);

			//            rez+=1;
			add1byPower(rez, power);
			//            std::cout <<"aux="<<aux<<"; BN2="<<BN2<<"; rez="<<rez<<"; rez1="<<rez1 <<"; "<<std::endl;
		}
		power--;
		shiftLEFT(aux, 1);
		//        std::cout <<"power="<<power<<std::endl;
	}
}

BigNumber BigNumber::operator/(BigNumber &oBN) {
	BigNumber rez(std::max(lenInt, oBN.lenInt), std::max(lenFrac, oBN.lenFrac));

	if (isZero()) return rez;
	if (BNisZero(oBN)) { std::cerr << "Divide by zero!"; return rez; }
	rez.sign = sign&&oBN.sign;
	if (BNisOne(oBN)) { rez = *this; return rez; }
	if (isEqual(*this, oBN)) { rez.pInt[rez.lenInt - 1] = 1; return rez; }
	//    std::cerr << "here!";
	div(*this, oBN, rez);

	return rez;

}

BigNumber BigNumber::operator/=(BigNumber &oBN) {
	BigNumber rez(std::max(lenInt, oBN.lenInt), std::max(lenFrac, oBN.lenFrac));

	if (isZero()) return rez;
	if (BNisZero(oBN)) { std::cerr << "Divide by zero!"; return rez; }
	rez.sign = sign&&oBN.sign;
	if (BNisOne(oBN)) { rez = *this; return rez; }
	if (isEqual(*this, oBN)) { rez.pInt[rez.lenInt - 1] = 1; return rez; }
	//    std::cerr << "here!";
	div(*this, oBN, rez);

	if (!copyValue(*this, rez)) std::cout << "COULD NOT COPPY!!!!!" << std::endl;
	return *this;

}

int BigNumber::getSignificatDecimals() {
	return BNgetSignificatDecimals(*this);
}

int BigNumber::BNgetSignificatDecimals(const BigNumber &BN) {
	bool search = true;
	int i = BN.lenFrac - 1;
	while ((i >= 0) && search) {
		if (BN.pFrac[i] != 0) search = false;
		i--;
	}

	if (search) return 0;
	else return i + 2;
}

int BigNumber::getSignificatIntPart() {
	return BNgetSignificatIntPart(*this);
}

int BigNumber::BNgetSignificatIntPart(const BigNumber &BN) {
	bool search = true;
	int i = BN.lenInt - 1;
	while ((i >0) && search) {
		if (BN.pInt[i] != 0) search = false;
		i--;
	}

	return i;
}

int BigNumber::getFistZerousInDecimal() {
	bool search = true;
	int i = 0;
	while ((i <lenFrac) && search) {
		if (pFrac[i] != 0) search = false;
		i++;
	}

	return i;
}


void BigNumber::compactList() {
	if (!sign) std::cout << "-";
	if (BNgetSignificatIntPart(*this) != 0) {
		for (int i = 0; i < lenInt; ++i) {
			if (pInt[i] != 0) std::cout << pInt[i];
		}
	}
	else std::cout << "0";

	if (BNgetSignificatDecimals(*this) != 0) {
		std::cout << ".";
		for (int j = 0; j < lenFrac; ++j) {
			if (pFrac[j] != 0) std::cout << pFrac[j];
		}
	}
	//    } else std::cout << "0";

	//    std::cout <<"("<<BNgetSignificatIntPart(*this)<<","<<BNgetSignificatDecimals(*this)<<")";
	//    std::cout <<"("<<lenInt<<","<<lenFrac<<")";
}

void BigNumber::withPrecissionList(const int precision) {
	if (!sign) std::cout << "-";
	if (BNgetSignificatIntPart(*this) != 0) {
		for (int i = 0; i < lenInt; ++i) {
			if (pInt[i] != 0) std::cout << pInt[i];
		}
	}
	else std::cout << "0";


	std::cout << ".";
	for (int j = 0; j < std::min(lenFrac, precision); ++j) {
		std::cout << pFrac[j];
	}
}
