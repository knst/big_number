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

#define DOUBLE_BASE

#ifndef DOUBLE_BASE
using bt = uint8_t;
using bt2 = uint16_t;
using bt2s = int16_t;
using bt4 = uint32_t;
using bt4s = int32_t;
constexpr bt2 bsize = 256;
constexpr bt bmax = 255;

#else // DOUBLE_BASE

using bt = uint16_t;
using bt2 = uint32_t;
using bt2s = int32_t;
using bt4 = uint64_t;
using bt4s = int64_t;
constexpr bt2 bsize = 65536;
constexpr bt bmax = 65535;

#endif

constexpr bt bz = sizeof(bt);
constexpr bt bz8 = sizeof(bt) * 8;

enum class RadixBN { hex, dec };

class BN {
public:
    static const BN bn0() noexcept;
    static const BN bn1() noexcept;

    // Generate random number with maximal size in byteCount.
    // Warning, if byteCount isn't divisor of sizeof(bt), than
    // count of bytes could be trunced.
    // Minimal size is one "bt".
    static BN makeRandom(size_t byteCount);
public:
    BN();

    //value. 0: fill 0
    BN(uint64_t basecount, int type);

    explicit BN(uint64_t x);
    BN(const BN&);
    BN(BN&& bn) noexcept;

    BN(std::vector<bt>&&) noexcept;
    BN(const std::vector<bt>&);

    BN(std::vector<bt>&&, size_t rbc);
    BN(const std::vector<bt>&, size_t rbc);

    BN(const std::string &, RadixBN = RadixBN::hex);

    void swap(BN& bn) noexcept;

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

    BN & operator = (const BN&);
    BN & operator = (BN&&) noexcept;

    const BN operator + (const BN& ) const;
    BN& operator += (const BN& );
    BN & operator ++();

    // This function return (this - bn).
    // Constrains: this >= bn
    // If this < bn, than result like in boost library.
    //
    const BN operator - (const BN&)const;
    BN& operator -= (const BN&);

    BN & operator --() noexcept;

    const BN operator * (const BN&)const;
    // classic multiplication O(n*n)
    const BN classicMultiplication(const BN&) const;
    // quick multiplication O(n*n)
    const BN fastMultiplication(const BN&) const;
    const BN karatsubaMultiplication(const BN&) const;
    void divmod(const BN& bn, BN& div, BN& mod) const;
    const BN operator / (const BN&)const;
    const BN operator % (const BN&)const;

    const BN operator >> (size_t shift)const;
    BN& operator >>= (size_t shift);
    const BN operator << (size_t shift) const;
    BN& operator <<= (size_t shift);

    bool operator <  (const BN&) const noexcept;
    bool operator <= (const BN&) const noexcept;
    bool operator >  (const BN&) const noexcept;
    bool operator >= (const BN&) const noexcept;
    bool operator == (const BN&) const noexcept;
    bool operator != (const BN&) const noexcept;
    bt   operator [] (size_t index_base) const noexcept;

    size_t digitCount() const noexcept;
    size_t bitCount() const noexcept;


    BN reductionBarrett(const BN& mod,const BN& mu) const;
    BN Pow(uint64_t)const;
    BN PowMod(uint64_t power, const BN& mod) const;
    BN PowMod(const BN& power, const BN& mod) const;

    BN PowModBarrett(const BN& power, const BN& mod) const;

    BN expRightToLeft(const BN& power, const BN& mod) const;

    std::vector<BN> expLeftToRightK_aryPrecomputation(const BN& mod) const;
    BN expLeftToRightK_ary(const BN& exponent, const BN& mod, const std::vector<BN>& g) const;

    std::vector<BN> expLeftToRightK_aryVarPrecomputation(const BN& mod, size_t K) const;
    BN expLeftToRightK_aryVar(const BN&, const BN&, const std::vector<BN>&, size_t K) const;

    std::vector<BN> expLeftToRightK_aryModifPrecomputation(const BN&) const;
    BN expLeftToRightK_aryMod(const BN&, const BN&, const std::vector<BN>&) const;

    std::vector <BN> expSlidingWindowPrecomputation(const BN&, size_t K) const;
    BN expSlidingWindow(const BN&, const BN&, const std::vector<BN>&, size_t k) const;

    //for best result:
    std::vector <BN> expBest_SlidePrecomp(const BN& mod) const;
    BN expBest_Slide(const BN& exponent, const BN& mod, const std::vector<BN>& g) const;

    BN Sqrt()const;
    BN Qrt()const;
    BN fastQrt()const;
    size_t countzeroright() const noexcept;
    bool bitI(size_t i) const noexcept;
    uint64_t get64() const noexcept;
    bool is0() const noexcept;
    bool isEven() const noexcept;

    const std::vector<bt> raw() const noexcept;

private:
    // data
    std::vector<bt> ba;
};

std::string to_string(BN bn);
std::string to_hexstring(const BN& bn);

#endif /* _BN_H */

