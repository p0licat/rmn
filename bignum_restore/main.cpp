#include <iostream>

//#include "BigNumber.h"
#include "BNMatrix.h"

using namespace std;

BigNumber compute_Exp(int precision, BigNumber x) {


	BigNumber rez, n, term, one, aux;

	rez = 1.0;
	term = 1.0;
	n = 0.0;
	one = 1.0;


	while (term.getFistZerousInDecimal() - 1<precision) {
		n += one;
		//        cout <<"; i="<<n<< endl ;
		aux = x / n;
		//        cout <<"; aux="<<aux<< endl ;
		term *= aux;
		//        cout <<"; term="<<term<< endl ;
		rez += term;
		//        cout <<"; rez="<<rez<< endl ;
	}

	return rez;
}

BNMatrix BuildExpMatrix(const int dim) {
	BNMatrix rez(dim, dim);
	double Tau, T2, TauMin, TauMax, T2Min, T2Max, dTau, dT2, Log10TauMin, Log10TauMax, Log10T2Min, Log10T2Max;
	int n, nTau, nT2;
	BigNumber A, B;

	nTau = dim;								// nr de puncte in timp
	nT2 = nTau;							// nr de puncte in T2
	n = nTau;							// ordinul matricii
	TauMin = 1e-3;							// timpul minim
	TauMax = 1e-2;							// timpul maxim
	T2Min = 1e-3;							// timpul minim
	T2Max = 1e-2;							// timpul maxim
	Log10TauMin = log10(TauMin);
	Log10TauMax = log10(TauMax);
	Log10T2Min = log10(T2Min);
	Log10T2Max = log10(T2Max);

	dTau = fabs(Log10TauMax - Log10TauMin) / (nTau - 1.0);
	dT2 = fabs(Log10T2Max - Log10T2Min) / (nT2 - 1.0);

	for (int i = 1; i <= nTau; i++) {
		Tau = pow(10, Log10TauMin + dTau*(i - 1));
		for (int j = 1; j <= nT2; j++) {
			T2 = pow(10, Log10T2Min + dT2*(j - 1));
			//            cout << endl <<"("<<i<<","<<j<<")"<<Tau <<" - " << T2<<endl;
			A = (-Tau / T2);
			// A.listBN();
			B = compute_Exp(150, A);
			// cout  <<"A:"<< A.getSign() <<A.getIntPart()<<"."<<A.getFracPart();
			//            cout << " elem:(" << i <<","<< j<<"="<<B<<endl;

			// B.listBN();
			rez.setElement(i, j, B);
		}
	}


	return rez;

}

//BNMatrix BuildExpMatrix(const int dim);

int main() {
	//    BigNumber x;
	//    x=-0.5;
	//    compute_Exp(150,x).withPrecissionList(150);

	BuildExpMatrix(3).toFile("out.txt", 30);

	//    , y, z, a1, a2, deter;
	//    x=1;
	//    y=7;
	//    z=5;
	//    a1=x*y; a1*=z;
	//    x=4; z=8;
	//    a2=x*y; a2*=z;
	//    a2=58;
	//    deter=a2; a1="-82";
	//    deter+=a1;
	//    cout << "det="<<a2<<"*"<<a1<<"="<< deter<<endl;

	//    cout << "Hello, World!" << endl;
	//    BigNumber X(-2.63219,1,5);
	//    BigNumber Y(1,7); BigNumber Z;
	////    X=-9.934765734;
	//    Y=5.2;
	//    Z="02.71828182845904523536028747135266249775724709369995";
	////    cout <<"X="<< X <<"; ("<<X.getLenFrac() <<")"<< endl ;
	//    cout <<"Y="<< Y <<"; ("<<Y.getLenFrac() <<")"<< endl ;
	//    cout <<"Z="<< Z <<"; ("<<Z.getLenFrac() <<")"<< endl ;


	//    Y.compactList();
	//    Z.compactList();

	//    if(Y.isZero()) cout << "Is ZERO!" << endl;
	//    if(X.isOne()) cout << "Is ONE!" << endl;
	//    if(Y!=Z) cout << "Y!=Z" << endl;
	//    else cout << "Y==Z" << endl;
	//
	//    cout << "---------------------------"<< endl ;
	//
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;
	//    adunare
	//    Z=X+Y;
	//    cout << "X+Y=" << Z<< endl;
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;
	////    scadere
	//    Z=Y-X;
	//    cout << "Y-X=" << Z<< endl;
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;
	//
	//    X-=Y;
	//    cout <<"X-=Y;X="<< X <<"; ("<<X.getLenFrac() <<")"<< endl ;
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;
	//
	//    Z=X>>3;
	//    cout << "X>>3=" << Z<< endl;
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;
	//
	//    Z=X<<3;
	//    cout << "X<<3=" << Z<< endl;
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;
	////    cout <<"X="<< X <<"; ("<<X.getLenFrac() <<")"<< endl ;

	cout << "---------------------------" << endl;

	//    compute_E(50);
	//    X="1.0";
	//    Y="40320.0";
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;
	//
	//    Z=X*Y;
	//    cout << "X*Y=" << Z<< endl;
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;
	//
	//    Z=X/Y;
	//    cout << "X/Y=" << Z<< endl;
	//    cout <<"X="<< X <<"; Y="<<Y <<"; Z="<<Z<< endl ;

	//    cout <<"Z="<< Z <<"; ("<<Z.getLenFrac() <<")"<< endl ;


	//    ifstream file;
	//    file.open("out.txt");      //open a file
	//    char output[1000];
	//    file>>output;       //write to it
	//    cout<<output;
	//    file.close();               //close it

	//    BNMatrix X(3,3);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//    BNMatrix X(2,2);
	//
	//    X.setElement(1,1,"0");
	//    X.setElement(1,2,"2");
	////    X.setElement(1,3,"3");
	//    X.setElement(2,1,"4");
	//    X.setElement(2,2,"5");
	////    X.setElement(2,3,"6");
	////    X.setElement(3,1,"7");
	////    X.setElement(3,2,"8");
	////    X.setElement(3,3,"9");
	//
	////    X.list();
	////    X.compactList();
	//    cout<<"X:";
	//    X.withPrecissionList(2);
	////
	////    X.listBNM();
	//    X.toFile("D:\\sanda\\_COLATERAL\\BigNumberMatrix_RaduFechete\\2015-12-20_rewrite\\out.txt", 3);
	//    BNMatrix Y;
	//    Y.fromFile("D:\\sanda\\_COLATERAL\\BigNumberMatrix_RaduFechete\\2015-12-20_rewrite\\input2x2.txt");
	//
	//    cout <<"dim1:" << Y.getDim1()<<"dim2:" << Y.getDim2()<<endl;
	//    cout<<"Y:";
	//    Y.withPrecissionList(2);
	//
	//    BNMatrix Z;
	////    cout<<"Z=";
	////    Z.withPrecissionList(2);
	////    Z=X;
	////    Z=Y+X;
	//    Z=X*Y;
	////    Z=X+Y;
	//    cout<<"Z:";
	////    cout <<"dim1:" << Z.getDim1()<<"dim2:" << Z.getDim2()<<endl;
	//    Z.withPrecissionList(2);
	////    Z.toFile("D:\\sanda\\_COLATERAL\\BigNumberMatrix_RaduFechete\\2015-12-20_rewrite\\out.txt", 3);
	//
	//    cout << "determinant Y:";
	//    BigNumber det;
	//    det=0;
	//    det=Y.getDeterminant();
	//    cout << det <<endl;
	//
	//
	//    X=Y.getInverse();
	//    Z=Y*X;
	//    X.withPrecissionList(2);
	//    cout <<endl;
	//    Z.withPrecissionList(2);
	//

	cout << "Done" << endl;
	return 0;
}

//BNMatrix BuildExpMatrix(const int dim){
//
//}