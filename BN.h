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

#ifndef DOUBLE_BASE
using bt = uint8_t;
using bt2 = uint16_t;
using bt2s = int16_t;
using bt4 = uint32_t;
constexpr bt2 bsize = 256;
constexpr bt bmax = 255;

#else // DOUBLE_BASE

using bt = uint16_t;
using bt2 = uint32_t;
using bt2s = int32_t;
using bt4 = uint64_t;
constexpr bt2 bsize = 65536;
constexpr bt bmax = 65535;

#endif

constexpr bt bz = sizeof(bt);
constexpr bt bz8 = sizeof(bt) * 8;

class BN {
public:
    static const BN bn0();
    static const BN bn1();

public:
    BN();

    //value. 0:fill 0, 1:fill 1, -1:rand()
    BN(uint64_t basecount, int type);

    explicit BN(uint64_t x);
    BN(const BN&);
    BN(BN&& bn);
    BN(const std::vector<bt>&, size_t rbc = 0);
    BN(const BN&, size_t start, size_t count = 0);
    BN(const std::string &, const int & base = 0);        // base - 0: Hex, 1: Dec

    void swap(BN& bn);

    // This function return this * base^t
    BN   mulbt(size_t t) const;

    // This function return this / base^t
    BN   divbt(size_t t) const;

    // This function return this % base^t
    BN   modbt(size_t mod) const;

    const BN mulbase(const bt&)const;
    BN& mulbaseappr(const bt&);
    const BN divbase(const bt&)const;
    BN& divbaseappr(const bt&);
    const BN modbase(const bt&)const;
    BN& modbaseappr(const bt&);

    // This function subtraction "bn" from "this".
    // Constrains: this >= bn
    // TODO: work slow
    void subappr(const BN&);

    BN & operator = (const BN&);
    BN & operator = (BN&&);
    const BN operator + (const BN&)const;
    BN & operator ++();

    // This function return (this - bn).
    // Constrains: this >= bn
    //
    const BN operator - (const BN&)const;

    BN & operator --();
    const BN operator * (const BN&)const;
    const BN fast_mul (const BN&)const;        //быстрый столбик
    const BN karatsuba (const BN&)const;
    const BN karatsuba_old (const BN&)const;
    void divmod(const BN& bn, BN& div, BN& mod) const;
    const BN operator / (const BN&)const;
    const BN operator % (const BN&)const;
    const BN operator >>(int shift)const;
    const BN operator <<(int shift) const;
    bool operator < (const BN&)const;
    bool operator <=(const BN&)const;
    bool operator > (const BN&)const;
    bool operator >=(const BN&)const;
    bool operator ==(const BN&)const;
    bool operator !=(const BN&)const;
    bt   operator [](size_t index_base)const;
    size_t digitCount()const;
    size_t bitCount()const;


    BN reduction_barrett(const BN& mod,const BN& mu) const;
    BN reduction_special(const BN& mod) const;
    BN Pow(uint64_t)const;
    BN PowMod(uint64_t power, const BN& mod) const;
    BN PowMod(const BN& power, const BN& mod) const;
    BN PowModBarrett(const BN& power, const BN& mod) const;
    BN expRightToLeft(const BN& power, const BN& mod) const;
    BN expLeftToRight(const BN& power, const BN& mod) const;
    std::vector <BN> expLeftToRightK_aryPrecomputation(const BN& mod) const;
    BN expLeftToRightK_ary(const BN& exponent, const BN& mod, const std::vector<BN>& g) const;
    std::vector <BN> expLeftToRightK_aryVarPrecomputation(const BN& mod, int K) const;
    BN expLeftToRightK_aryVar(BN, BN, std::vector <BN>, int K)const;
    std::vector <BN> expLeftToRightK_aryModifPrecomputation(BN)const;
    BN expLeftToRightK_aryMod(BN, BN, std::vector <BN> )const;
    std::vector <BN> expSlidingWindowPrecomputation(BN, int)const;
    BN expSlidingWindow(BN,BN, std::vector <BN>, int k)const;

    //for best result:
    std::vector <BN> expBest_SlidePrecomp(BN mod) const;
    BN expBest_Slide(BN exponent, BN mod, std::vector <BN> g) const;

    BN Sqrt()const;
    BN Qrt()const;
    BN fastQrt()const;
    int countzeroright()const;
    bool bitI(size_t i)const;
    operator uint64_t()const;
    bool is0() const;
    bool isEven() const;
    void PrintHex(bool newstr=true)const;
    void PrintDec(bool newstr=true)const;

private:
    void InitMemory(int type);
    // Normalization of BN: changing rbc in according to value.
    void Norm();

    BN reduction_barrett_precomputation()const;

    BN karatsuba_add(size_t start, size_t count) const;      //для складывания частей числа для Карацубы
    // x += y.mulbt(mul_bt);
    BN karatsubaRecursive(const BN& bn, size_t start, size_t len) const;

private:
    // used memory (base count)
    size_t rbc;
    // data
    std::vector<bt> ba;

};

#endif /* _BN_H */

