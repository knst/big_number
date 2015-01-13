/*
 * File:   BN.h
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:16
 */

#ifndef _BN_H
#define _BN_H

#include <string>
#include <vector>

using namespace std;

#ifndef DOUBLE_BASE
using bt = uint8_t;
using bt2 = uint16_t;
using bt2s = int16_t;
using bt4 = uint32_t;
constexpr bt2 bsize = 256;
constexpr bt2 bmax = 255;

#else // DOUBLE_BASE

using bt = uint16_t;
using bt2 = uint32_t;
using bt2s = int32_t;
using bt4 = uint64_t;
constexpr bt2 bsize = 65536;
constexpr bt2 bmax = 65535;

#endif

constexpr bt bz = sizeof(bt);
constexpr bt bz2 = sizeof(bt2);
constexpr bt bz8 = sizeof(bt) * 8;

constexpr size_t karacuba_const = 50;

class BN {
private:
    vector<bt> ba;
    // used memory (base count)
    size_t rbc;
    void InitMemory(int type = 2);
    int Norm();                             //пересчитать rbc

public:
    BN reduction_barrett_precomputation()const;

private:
    bt qCompute(bt,int,const BN&)const;
    bool lessorequal(const BN&,const int&)const;//считается, что bn справа дополнен shift нулями
    BN subequalshift(BN&,int)const;                //считается, что bn справа дополнен shift нулями
    BN mulMontgomery(const BN& bn, const BN& mod, bt m1) const;         //MP
    BN transformationMontgomery(const BN & mod, bt m1) const;           //MR
    BN karatsuba_add(const BN & bn, int start_1, int count_1, int start_2, int count_2) const;      //для складывания частей числа для Карацубы
    BN add_appr(const BN & bn, int mul_bt);       // x += y.mulbt(mul_bt2);
    BN karatsubaRecursive(const BN & bn, int start, int len) const;
public:
    BN();

    //value. 0:fill 0, 1:fill 1, -1:rand()
    BN(uint64_t basecount, int type);

    BN(uint64_t x);
    BN(const BN&);
    BN(const BN&, int start, int count=-1);
    BN(const string &, const int & base = 0);        // base - 0: Hex, 1: Dec
    BN   mulbt(const int&)const;                // x * base^t
    BN   divbt(const int&)const;                // x / base^t
    BN   modbt(const int&)const;                // x % base^t
    BN   mulbase(const bt&)const;
    BN   mulbaseappr(const bt&);
    BN   divbase(const bt&)const;
    BN   divbaseappr(const bt&);
    BN   modbase(const bt&)const;
    BN   modbaseappr(const bt&);

    // This function return (this - bn).
    // Constrains: this >= bn
    BN   sub(const BN&)const;

    BN & operator = (const BN&);
    BN   operator + (const BN&)const;
    BN & operator ++();
    BN   operator - (const BN&)const;        // result = abs(*this - bn)
    BN & operator --();
    BN   operator * (const BN&)const;
    BN   fast_mul (const BN&)const;        //быстрый столбик
    BN   karatsuba (const BN&)const;
    BN   karatsuba_old (const BN&)const;
    BN   operator / (const BN&)const;
    BN   operator % (const BN&)const;
    BN   operator >>(int shift)const;
    BN   operator <<(int shift) const;
    bool operator < (const BN&)const;
    bool operator <=(const BN&)const;
    bool operator > (const BN&)const;
    bool operator >=(const BN&)const;
    bool operator ==(const BN&)const;
    bool operator !=(const BN&)const;
    bt   operator [](const int)const;
    int basecount()const;
    int bitcount()const;


    BN reduction_montgomery(const BN& mod, bt m1, const BN& T) const;
    BN reduction_barrett(const BN& mod,const BN& mu) const;
    BN reduction_special(const BN& mod) const;
    BN Pow(uint64_t)const;
    BN PowMod(uint64_t power, const BN& mod) const;
    BN PowMod(const BN& power, const BN& mod) const;
    BN PowModBarrett(const BN& power, const BN& mod) const;
    BN expRightToLeft(const BN& power, const BN& mod) const;
    BN expLeftToRight(const BN& power, const BN& mod) const;
    vector <BN> expLeftToRightK_aryPrecomputation(const BN& mod) const;
    BN expLeftToRightK_ary(const BN& exponent, const BN& mod, const vector<BN>& g) const;
    vector <BN> expLeftToRightK_aryVarPrecomputation(const BN& mod, int K) const;
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
    bool bitI(size_t i)const;
    operator uint64_t()const;
    bool is0() const;
    bool isEven() const;
    void Print(bool newstr=true)const;
    void PrintHex(bool newstr=true)const;
    void PrintDec(bool newstr=true)const;
};

#endif /* _BN_H */

