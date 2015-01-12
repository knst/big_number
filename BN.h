/*
 * File:   BN.h
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:16
 */

#ifndef _BN_H
#define	_BN_H

//TODO: переделать функции со сдвигами!
//TODO: переделка функций сдвига не повлияла на скорость!

#include <string>
#include <vector>

using namespace std;

#define ull	unsigned long long
#define ll	long long

#define bt	unsigned char
#define bt2	unsigned short
#define bt2s	short
#define bt4	unsigned long
#define bsize	256
#define bmax	255

//#define bt	unsigned short
//#define bt2	unsigned int
//#define bt2s	int
//#define bt4	unsigned long long
//#define bsize	(bt2)65536
//#define bmax	(bt2)65535

/*#define bt	unsigned long
#define bt2	unsigned long long
#define bsize	(bt2)(4294967296)
#define bmax	(bt2)(4294967295)*/

#define bz	(sizeof(bt))
#define bz2	(sizeof(bt2))
#define bz8	(sizeof(bt)*8)

#define caracuba_const  50

class BN {
private:
    bt *ba;
    int bc;                                 //количество выделенных баз
    int rbc;                                //количество реально используемых баз
    void GetMemory(int value=2);
    void FreeMemory();
    int Norm();                             //пересчитать rbc

public:
    BN reduction_barrett_precomputation()const;
private:
    bt qCompute(bt,int,const BN&)const;
    bool lessorequal(const BN&,const int&)const;//считается, что bn справа дополнен shift нулями
    BN subequalshift(BN&,int)const;		//считается, что bn справа дополнен shift нулями
    BN mulMontgomery(const BN& bn, const BN& mod, bt m1) const;         //MP
    BN transformationMontgomery(const BN & mod, bt m1) const;           //MR
    BN karatsuba_add(const BN & bn, int start_1, int count_1, int start_2, int count_2) const;      //для складывания частей числа для Карацубы
    BN add_appr(const BN & bn, int mul_bt);       // x += y.mulbt(mul_bt2);
    BN karatsubaRecursive(const BN & bn, int start, int len) const;
public:
    BN();
    BN(ull basecount,const int &value);         //[value] 0:fill 0; 1:fill 1; -1:rand().
    BN(ull x);
    BN(const BN&);
    BN(const BN&, int start, int count=-1);
    BN(const string &,const int & base = 0);	// base - 0: Hex, 1: Dec
    BN   mulbt(const int&)const;		// x * base^t
    BN   divbt(const int&)const;		// x / base^t
    BN   modbt(const int&)const;		// x % base^t
    BN   mulbase(const bt&)const;
    BN   mulbaseappr(const bt&);
    BN   divbase(const bt&)const;
    BN   divbaseappr(const bt&);
    BN   modbase(const bt&)const;
    BN   modbaseappr(const bt&);
    BN   sub(const BN&)const;		// result = *this - bn. Только если *this >=bn !!!
    BN & operator = (const BN&);
    BN   operator + (const BN&)const;
    BN & operator ++();
    BN   operator - (const BN&)const;	// result = abs(*this - bn)
    BN & operator --();
    BN   operator * (const BN&)const;
    BN   fast_mul (const BN&)const;	//быстрый столбик
    BN   karatsuba (const BN&)const;
    BN   karatsuba_old (const BN&)const;
    BN   operator / (const BN&)const;
    BN   operator % (const BN&)const;
    BN   operator >>(const int&)const;
    BN   operator <<(const int&)const;
    bool operator < (const BN&)const;
    bool operator <=(const BN&)const;
    bool operator > (const BN&)const;
    bool operator >=(const BN&)const;
    bool operator ==(const BN&)const;
    bool operator !=(const BN&)const;
    bt   operator [](const int)const;
    int basecount()const;
    int bitcount()const;


    BN reduction_montgomery (const BN& mod, bt m1, BN T ) const;
    BN reduction_barrett    (const BN&,const BN&	) const;
    BN reduction_special    (const BN&                  ) const;
    BN Pow(ull)const;
    BN PowMod(ull,BN)const;
    BN PowMod(BN,BN)const;
    BN PowModBarrett(BN, BN)const;
    BN expRightToLeft(BN, BN)const;
    BN expLeftToRight(BN, BN)const;
    vector <BN> expLeftToRightK_aryPrecomputation(BN)const;
    BN expLeftToRightK_ary(BN, BN, vector <BN>)const;
    vector <BN> expLeftToRightK_aryVarPrecomputation(BN, int K)const;
    BN expLeftToRightK_aryVar(BN, BN, vector <BN>, int K)const;
    vector <BN> expLeftToRightK_aryModifPrecomputation(BN)const;
    BN expLeftToRightK_aryMod(BN, BN, vector <BN> )const;
    vector <BN> expSlidingWindowPrecomputation(BN, int)const;
    BN expSlidingWindow(BN,BN, vector <BN>, int k)const;
    BN expMontgomery(BN exponent, BN mod) const;

    //for best result:
    vector <BN> expBest_SlidePrecomp(BN mod) const;
    BN expBest_Slide(BN exponent, BN mod, vector <BN> g) const;

    BN Sqrt()const;
    BN Qrt()const;
    int countzeroright()const;
    bool bitI(unsigned int i)const;
    operator ull()const;
    bool is0()const;
    bool isEven()const;
    void Print(bool newstr=true)const;
    void PrintHex(bool newstr=true)const;
    void PrintDec(bool newstr=true)const;
    ~BN();
};

#endif	/* _BN_H */

